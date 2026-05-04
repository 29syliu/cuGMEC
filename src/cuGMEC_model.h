/*
 * cuGMEC
 *
 * Copyright (C) 2025 Shiyang Liu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <mpi.h>
#include <nccl.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cub/cub.cuh>
#include <cudss.h>
#include <cufft.h>
#include <cusparse_v2.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "cuGMEC_setup.h"

#define RESET "\033[0m"
#define RED "\033[31m"
#define BLUE "\033[34m"
#define CYAN "\033[36m"
#define BLACK "\033[30m"
#define WHITE "\033[37m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define MAGENTA "\033[35m"
#define BOLDRED "\033[1m\033[31m"
#define BOLDBLUE "\033[1m\033[34m"
#define BOLDCYAN "\033[1m\033[36m"
#define BOLDBLACK "\033[1m\033[30m"
#define BOLDWHITE "\033[1m\033[37m"
#define BOLDGREEN "\033[1m\033[32m"
#define BOLDYELLOW "\033[1m\033[33m"
#define BOLDMAGENTA "\033[1m\033[35m"

// clang-format off

#define CHECK_CALL(cmd, ok, label, fmt, arg)                                                                           \
    do {                                                                                                               \
        auto err = (cmd);                                                                                              \
        if (err != (ok)) {                                                                                             \
            printf("Failed: " label " error %s:%d '" fmt "'\n", __FILE__, __LINE__, (arg));                            \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define MPICHECK(cmd)      CHECK_CALL(cmd, MPI_SUCCESS,             "MPI",      "%d", err)
#define NCCLCHECK(cmd)     CHECK_CALL(cmd, ncclSuccess,             "NCCL",     "%s", ncclGetErrorString(err))
#define CUDACHECK(cmd)     CHECK_CALL(cmd, cudaSuccess,             "CUDA",     "%s", cudaGetErrorString(err))
#define CUDSSCHECK(cmd)    CHECK_CALL(cmd, CUDSS_STATUS_SUCCESS,    "CUDSS",    "%d", err)
#define CUSPARSECHECK(cmd) CHECK_CALL(cmd, CUSPARSE_STATUS_SUCCESS, "CUSPARSE", "%d", err)
#define CUFFTCHECK(cmd)    CHECK_CALL(cmd, CUFFT_SUCCESS,           "CUFFT",    "%d", err)
struct HybridModelConfig {
    struct HyperDiffPair  { bool enable; double coef; };

    struct Scale          { int devNums, hostNums, hostId, localId, gridNx, gridNy, gridNz, gridGhost, ppcNums, tubes; } scale;
    struct Normalization  { double B0, L0, VA0, RHO0, RHO1, PSITMAX, dt; } norm;
    struct HyperDiffusion { HyperDiffPair A, Phi, dNe, dTe, dPi, dPa, dPb; } hyperDiff;
    struct SpeciesConfig  { bool enable; double mass, charge, beta, vmin, vmax, vb, deltaV, lambda0, deltaLambda2; } ion, alpha, beam;
    struct Time           { int total, diag, output; } time;
    struct Filter         { int leftN, rightN; } filter;
    struct FFTSize        { int time, batch, freq; };  FFTSize nFFT, mFFT;
    struct DiagFlags      { bool amplitude, frequency, Eparallel, density, diffusivity, ZFDrive, checkNAN; } diag;
    struct OutputFlags    { bool Phi, A, dNe, dTe, dPi, dPa, dPb; } output;
};
// clang-format on

template <typename mhdReal, typename picReal>
class HybridModel {

  public:
    HybridModel(const HybridModelConfig& cfg)
        :

          devNums{cfg.scale.devNums}, hostNums{cfg.scale.hostNums}, hostId{cfg.scale.hostId},
          localId{cfg.scale.localId}, gridNx{cfg.scale.gridNx}, gridNy{cfg.scale.gridNy}, gridNz{cfg.scale.gridNz},
          gridGhost{cfg.scale.gridGhost}, ppcNums{cfg.scale.ppcNums}, tubes{cfg.scale.tubes},

          hostNy{gridNy / hostNums}, devNy{hostNy / devNums}, gridNxz{gridNx * gridNz}, gridNxPlusGhost{gridNx},
          gridNyPlusGhost{gridNy + 2 * gridGhost}, gridNzPlusGhost{gridNz + 2 * gridGhost}, cellNx{gridNxPlusGhost - 1},
          cellNy{gridNyPlusGhost - 1}, cellNz{gridNzPlusGhost - 1}, cellNxz{cellNx * cellNz},
          picHost{gridNx * gridNy * gridNz / hostNums * ppcNums}, picDev{picHost / devNums}, gridDx{1.0 / (gridNx - 1)},
          gridDy{2.0 * PI / gridNy}, gridDz{2.0 * PI / tubes / gridNz}, x0{0.0}, x1{1.0}, x0PlusGhost{x0},
          x1PlusGhost{x1}, y0{-PI}, y1{PI}, y0PlusGhost{y0 - (gridGhost - 0.5) * gridDy},
          y1PlusGhost{y1 + (gridGhost - 0.5) * gridDy}, z0{-PI / tubes}, z1{PI / tubes},
          z0PlusGhost{z0 - (gridGhost - 0.5) * gridDz}, z1PlusGhost{z1 + (gridGhost - 0.5) * gridDz},

          B0{cfg.norm.B0}, L0{cfg.norm.L0}, VA0{cfg.norm.VA0}, RHO0{cfg.norm.RHO0}, RHO1{cfg.norm.RHO1},
          PSITMAX{cfg.norm.PSITMAX}, NormQE{QE / (B0 * L0 * L0 / MU0 / VA0)}, dt{cfg.norm.dt},

          nablaPerp2A{cfg.hyperDiff.A}, nablaPerp2Phi{cfg.hyperDiff.Phi}, nablaPerp2dNe{cfg.hyperDiff.dNe},
          nablaPerp2dTe{cfg.hyperDiff.dTe}, nablaPerp2dPi{cfg.hyperDiff.dPi}, nablaPerp2dPa{cfg.hyperDiff.dPa},
          nablaPerp2dPb{cfg.hyperDiff.dPb},

          ifIon{cfg.ion.enable}, IonMass{cfg.ion.mass}, IonChar{cfg.ion.charge}, IonBeta{cfg.ion.beta},
          IonVmin{cfg.ion.vmin}, IonVmax{cfg.ion.vmax}, IonVb{cfg.ion.vb}, IonDeltaV{cfg.ion.deltaV},
          IonLambda0{cfg.ion.lambda0}, IonDeltaLambda2{cfg.ion.deltaLambda2},

          ifAlpha{cfg.alpha.enable}, AlphaMass{cfg.alpha.mass}, AlphaChar{cfg.alpha.charge}, AlphaBeta{cfg.alpha.beta},
          AlphaVmin{cfg.alpha.vmin}, AlphaVmax{cfg.alpha.vmax}, AlphaVb{cfg.alpha.vb}, AlphaDeltaV{cfg.alpha.deltaV},
          AlphaLambda0{cfg.alpha.lambda0}, AlphaDeltaLambda2{cfg.alpha.deltaLambda2},

          ifBeam{cfg.beam.enable}, BeamMass{cfg.beam.mass}, BeamChar{cfg.beam.charge}, BeamBeta{cfg.beam.beta},
          BeamVmin{cfg.beam.vmin}, BeamVmax{cfg.beam.vmax}, BeamVb{cfg.beam.vb}, BeamDeltaV{cfg.beam.deltaV},
          BeamLambda0{cfg.beam.lambda0}, BeamDeltaLambda2{cfg.beam.deltaLambda2},

          totalSteps{cfg.time.total}, diagSteps{cfg.time.diag}, outputSteps{cfg.time.output}, leftN{cfg.filter.leftN},
          rightN{cfg.filter.rightN}, nFFTTimeSize{cfg.nFFT.time}, nFFTBatchSize{cfg.nFFT.batch},
          nFFTFreqSize{cfg.nFFT.freq}, mFFTTimeSize{cfg.mFFT.time}, mFFTBatchSize{cfg.mFFT.batch},
          mFFTFreqSize{cfg.mFFT.freq},

          ifDiagAmplitude{cfg.diag.amplitude}, ifDiagFrequency{cfg.diag.frequency}, ifDiagEparallel{cfg.diag.Eparallel},
          ifDiagDensity{cfg.diag.density}, ifDiagDiffusivity{cfg.diag.diffusivity}, ifDiagZFDrive{cfg.diag.ZFDrive},
          ifCheckNAN{cfg.diag.checkNAN},

          ifOutputPhi{cfg.output.Phi}, ifOutputA{cfg.output.A}, ifOutputdNe{cfg.output.dNe},
          ifOutputdTe{cfg.output.dTe}, ifOutputdPi{cfg.output.dPi}, ifOutputdPa{cfg.output.dPa},
          ifOutputdPb{cfg.output.dPb} {

        if (hostId == 0) {

            {
                std::cout << BOLDRED << std::endl;
                std::cout << "  ______   _     _    ______   __    __   ________  ______" << std::endl;
                std::cout << " / _____| | |   | |  / _____| |  \\  /  | |  _____/ / _____|" << std::endl;
                std::cout << "| |       | |   | | | |       |   \\/   | | |      | |" << std::endl;
                std::cout << "| |       | |   | | | |       |        | | |___   | |" << std::endl;
                std::cout << "| |       | |   | | | |    _  | |\\  /| | |  ___|  | |" << std::endl;
                std::cout << "| |       | |   | | | |   | | | | \\/ | | | |      | |" << std::endl;
                std::cout << "| |_____  | |___| | | |___| | | |    | | | |_____ | |_____" << std::endl;
                std::cout << " \\______|  \\_____/   \\______| |_|    |_| |_______\\ \\______|" << std::endl;
                std::cout << RESET << std::endl;
            }

            std::cout << std::endl;
            std::cout << BOLDYELLOW << "Start: Initialize Gyrokinetic-MHD Hybrid Model." << std::endl;
            std::cout << "gridNx: " << gridNx << ", gridNy: " << gridNy << ", gridNz: " << gridNz << "." << std::endl;
            std::cout << "hostNums: " << hostNums << ", devNums: " << devNums << "." << std::endl;
            std::cout << "picHost: " << picHost << ", picDev: " << picDev << "." << std::endl;
            logDone();
        }
    }

    ~HybridModel() {}
    // clang-format off
	/*------------------------------------------------------Model Parameters------------------------------------------------------*/

	const double QE = 1.6021766208e-19;
	const double MP = 1.672621637e-27;
	const double PI = 3.1415926535897932;
	const double MU0 = 4.0 * PI * 1.0e-7;
	const double KEV = 1000.0 * QE;

	const int devNums;
	const int hostNums, hostId, localId;
	const int gridNx, gridNy, gridNz;
	const int gridGhost, ppcNums, tubes;

	const int hostNy, devNy, gridNxz;
	const int gridNxPlusGhost, gridNyPlusGhost, gridNzPlusGhost;
	const int cellNx, cellNy, cellNz, cellNxz;
	const int picHost, picDev;
	const double gridDx, gridDy, gridDz;
	const double x0, x1, x0PlusGhost, x1PlusGhost;
	const double y0, y1, y0PlusGhost, y1PlusGhost;
	const double z0, z1, z0PlusGhost, z1PlusGhost;

	const double B0, L0, VA0, RHO0, RHO1, PSITMAX, NormQE, dt;

	const HybridModelConfig::HyperDiffPair nablaPerp2A;
	const HybridModelConfig::HyperDiffPair nablaPerp2Phi;
	const HybridModelConfig::HyperDiffPair nablaPerp2dNe;
	const HybridModelConfig::HyperDiffPair nablaPerp2dTe;
	const HybridModelConfig::HyperDiffPair nablaPerp2dPi;
	const HybridModelConfig::HyperDiffPair nablaPerp2dPa;
	const HybridModelConfig::HyperDiffPair nablaPerp2dPb;

	const bool ifIon;
	const double IonMass;
	const double IonChar;
	const double IonBeta;
	const double IonVmin;
	const double IonVmax;
	const double IonVb;
	const double IonDeltaV;
	const double IonLambda0;
	const double IonDeltaLambda2;

	const bool ifAlpha;
	const double AlphaMass;
	const double AlphaChar;
	const double AlphaBeta;
	const double AlphaVmin;
	const double AlphaVmax;
	const double AlphaVb;
	const double AlphaDeltaV;
	const double AlphaLambda0;
	const double AlphaDeltaLambda2;

	const bool ifBeam;
	const double BeamMass;
	const double BeamChar;
	const double BeamBeta;
	const double BeamVmin;
	const double BeamVmax;
	const double BeamVb;
	const double BeamDeltaV;
	const double BeamLambda0;
	const double BeamDeltaLambda2;

	/*---------------------------------------------------MHD Equilibrium on CPU---------------------------------------------------*/

	double** q; double** q_px;
	double** psip; double** psip_px;

	double** Ni; double** Ni_px;
	double** Ti; double** Ti_px;
	double** Pi; double** Pi_px;
	double** Ne; double** Ne_px;
	double** Te; double** Te_px;
	double** Pe; double** Pe_px;
	double** Na; double** Na_px;
	double** Ta; double** Ta_px;
	double** Nb; double** Nb_px;
	double** Tb; double** Tb_px;

	double** B; double** B_px; double** B_py;
	double** B_px2; double** B_pxy; double** B_py2;

	double** J; double** J_px; double** J_py;
	double** Bny; double** Bny_px; double** Bny_py;
	double** Va; double** Va_px; double** Va_py;
	double** Rho; double** Rho_px; double** Rho_py;
	double** JpB; double** JpB_px; double** JpB_py;
	double** R; double** Z;

	double** gconxx; double** gconxx_px; double** gconxx_py;
	double** gconxy; double** gconxy_px; double** gconxy_py;
	double** gconxz; double** gconxz_px; double** gconxz_py;
	double** gconyy; double** gconyy_px; double** gconyy_py;
	double** gconyz; double** gconyz_px; double** gconyz_py;
	double** gconzz; double** gconzz_px; double** gconzz_py;
	double** gcovxx; double** gcovxx_px; double** gcovxx_py;
	double** gcovxy; double** gcovxy_px; double** gcovxy_py;
	double** gcovxz; double** gcovxz_px; double** gcovxz_py;
	double** gcovyy; double** gcovyy_px; double** gcovyy_py;
	double** gcovyz; double** gcovyz_px; double** gcovyz_py;
	double** gcovzz; double** gcovzz_px; double** gcovzz_py;

	double** SFAconxx; double** SFAconxx_px; double** SFAconxx_py;
	double** SFAconxy; double** SFAconxy_px; double** SFAconxy_py;
	double** SFAconxz; double** SFAconxz_px; double** SFAconxz_py;
	double** SFAconyy; double** SFAconyy_px; double** SFAconyy_py;
	double** SFAconyz; double** SFAconyz_px; double** SFAconyz_py;
	double** SFAconzz; double** SFAconzz_px; double** SFAconzz_py;
	double** SFAcovxx; double** SFAcovxx_px; double** SFAcovxx_py;
	double** SFAcovxy; double** SFAcovxy_px; double** SFAcovxy_py;
	double** SFAcovxz; double** SFAcovxz_px; double** SFAcovxz_py;
	double** SFAcovyy; double** SFAcovyy_px; double** SFAcovyy_py;
	double** SFAcovyz; double** SFAcovyz_px; double** SFAcovyz_py;
	double** SFAcovzz; double** SFAcovzz_px; double** SFAcovzz_py;

	/*--------------------------------------------------MHD Perturbation on CPU---------------------------------------------------*/

	mhdReal** h_qtheta;

	mhdReal*** h_w;
	mhdReal*** h_A;
	mhdReal*** h_dNe;
	mhdReal*** h_dTe;
	mhdReal*** h_Phi;
	mhdReal*** h_dJpB;
	mhdReal*** h_dPe;

	/*-----------------------------------------------Coefficient Compression on CPU-----------------------------------------------*/

	/*-------------------------------------------Linear-------------------------------------------*/

	mhdReal** h_A_w; mhdReal** h_A_px_w; mhdReal** h_A_py_w; mhdReal** h_A_pz_w;
	mhdReal** h_dJpB_w; mhdReal** h_dJpB_px_w; mhdReal** h_dJpB_py_w; mhdReal** h_dJpB_pz_w;
	mhdReal** h_dP_w; mhdReal** h_dP_px_w; mhdReal** h_dP_py_w; mhdReal** h_dP_pz_w;
	mhdReal** h_w_py_w; mhdReal** h_w_pz_w; mhdReal** h_w_Phi;

	double** h_Phi_w;
	double** h_Phi_px_w; double** h_Phi_pz_w;
	double** h_Phi_px2_w;  double** h_Phi_pxz_w; double** h_Phi_pz2_w;

	double** h_A_resistive;
	double** h_A_px_resistive; double** h_A_pz_resistive;
	double** h_A_px2_resistive; double** h_A_pxz_resistive; double** h_A_pz2_resistive;

	double** h_F_perp2;
	double** h_F_px_perp2; double** h_F_pz_perp2;
	double** h_F_px2_perp2; double** h_F_pxz_perp2; double** h_F_pz2_perp2;

	mhdReal** h_A_dJpB;
	mhdReal** h_A_px_dJpB;  mhdReal** h_A_pz_dJpB;
	mhdReal** h_A_px2_dJpB;  mhdReal** h_A_pxz_dJpB; mhdReal** h_A_pz2_dJpB;

	mhdReal** h_Phi_A; mhdReal** h_Phi_px_A; mhdReal** h_Phi_py_A; mhdReal** h_Phi_pz_A;
	mhdReal** h_dNe_A; mhdReal** h_dNe_px_A; mhdReal** h_dNe_py_A; mhdReal** h_dNe_pz_A;
	mhdReal** h_A_A; mhdReal** h_A_px_A; mhdReal** h_A_py_A; mhdReal** h_A_pz_A;

	mhdReal** h_Phi_dNe; mhdReal** h_Phi_px_dNe; mhdReal** h_Phi_py_dNe; mhdReal** h_Phi_pz_dNe;
	mhdReal** h_dPe_dNe; mhdReal** h_dPe_px_dNe; mhdReal** h_dPe_py_dNe; mhdReal** h_dPe_pz_dNe;
	mhdReal** h_dJpB_dNe; mhdReal** h_dJpB_px_dNe; mhdReal** h_dJpB_py_dNe; mhdReal** h_dJpB_pz_dNe;
	mhdReal** h_A_dNe; mhdReal** h_A_px_dNe; mhdReal** h_A_py_dNe; mhdReal** h_A_pz_dNe;

	mhdReal** h_Phi_dTe; mhdReal** h_Phi_px_dTe; mhdReal** h_Phi_py_dTe; mhdReal** h_Phi_pz_dTe;
	mhdReal** h_dTe_dTe; mhdReal** h_dTe_px_dTe; mhdReal** h_dTe_py_dTe; mhdReal** h_dTe_pz_dTe;
	mhdReal** h_dNe_dTe; mhdReal** h_dNe_px_dTe; mhdReal** h_dNe_py_dTe; mhdReal** h_dNe_pz_dTe;

	mhdReal** h_Ne0; mhdReal** h_Te0;
	mhdReal** h_Ne0_px; mhdReal** h_Te0_px; mhdReal** h_Pe0_px;

	mhdReal*** h_F2perp2;
	mhdReal*** h_A2dJpB;
	mhdReal*** h_Phi2w;
	mhdReal*** h_wdPAdJpB2w;
	mhdReal*** h_APhidNe2A;
	mhdReal*** h_dPePhiAdJpB2dNe;
	mhdReal*** h_PhidTedNe2dTe;

	/*-----------------------------------------Nonlinear------------------------------------------*/

	mhdReal*** h_wPhi_w; mhdReal*** h_AdJpB_w;
	mhdReal*** h_PhiA_A; mhdReal*** h_NeA_A;
	mhdReal*** h_AdJpB_dNe; mhdReal*** h_dNePhi_dNe;
	mhdReal*** h_PhiTe_dTe; mhdReal*** h_PhiTeA_dTe;

	/*-------------------------------------------------------Matrix on CPU--------------------------------------------------------*/

	std::vector<int> matrix_i;
	std::vector<int> matrix_j;
	std::vector<mhdReal> matrix_v;

	/*--------------------------------------------------MHD Perturbation on GPU---------------------------------------------------*/

	mhdReal** d_qtheta;

	mhdReal** d_w;
	mhdReal** d_A;
	mhdReal** d_dNe;
	mhdReal** d_dTe;
	mhdReal** d_Phi;
	mhdReal** d_dJpB;
	mhdReal** d_dPe;

	mhdReal** d_w_beg; mhdReal** d_w_midl; mhdReal** d_w_midr; mhdReal** d_w_end;
	mhdReal** d_A_beg; mhdReal** d_A_midl; mhdReal** d_A_midr; mhdReal** d_A_end;
	mhdReal** d_dNe_beg; mhdReal** d_dNe_midl; mhdReal** d_dNe_midr; mhdReal** d_dNe_end;
	mhdReal** d_dTe_beg; mhdReal** d_dTe_midl; mhdReal** d_dTe_midr; mhdReal** d_dTe_end;
	mhdReal** d_Phi_midl; mhdReal** d_Phi_midr;
	mhdReal** d_dJpB_midl; mhdReal** d_dJpB_midr;
	mhdReal** d_dPe_midl; mhdReal** d_dPe_midr;
	mhdReal** d_dPi_midl; mhdReal** d_dPi_midr;
	mhdReal** d_dPa_midl; mhdReal** d_dPa_midr;
	mhdReal** d_dPb_midl; mhdReal** d_dPb_midr;
	mhdReal** d_Apt_midl; mhdReal** d_Apt_midr;

	/*-----------------------------------------------Coefficient Compression on GPU-----------------------------------------------*/

	/*-------------------------------------------Linear-------------------------------------------*/

	mhdReal** d_A_w; mhdReal** d_A_px_w; mhdReal** d_A_py_w; mhdReal** d_A_pz_w;
	mhdReal** d_dJpB_w; mhdReal** d_dJpB_px_w; mhdReal** d_dJpB_py_w; mhdReal** d_dJpB_pz_w;
	mhdReal** d_dP_w; mhdReal** d_dP_px_w; mhdReal** d_dP_py_w; mhdReal** d_dP_pz_w;
	mhdReal** d_w_py_w; mhdReal** d_w_pz_w; mhdReal** d_w2Phi;

	mhdReal** d_F_perp2;
	mhdReal** d_F_px_perp2; mhdReal** d_F_pz_perp2;
	mhdReal** d_F_px2_perp2; mhdReal** d_F_pxz_perp2; mhdReal** d_F_pz2_perp2;

	mhdReal** d_A_dJpB;
	mhdReal** d_A_px_dJpB; mhdReal** d_A_pz_dJpB;
	mhdReal** d_A_px2_dJpB;  mhdReal** d_A_pxz_dJpB; mhdReal** d_A_pz2_dJpB;

	mhdReal** d_Phi_A; mhdReal** d_Phi_px_A; mhdReal** d_Phi_py_A; mhdReal** d_Phi_pz_A;
	mhdReal** d_dNe_A; mhdReal** d_dNe_px_A; mhdReal** d_dNe_py_A; mhdReal** d_dNe_pz_A;
	mhdReal** d_A_A; mhdReal** d_A_px_A; mhdReal** d_A_py_A; mhdReal** d_A_pz_A;

	mhdReal** d_Phi_dNe; mhdReal** d_Phi_px_dNe; mhdReal** d_Phi_py_dNe; mhdReal** d_Phi_pz_dNe;
	mhdReal** d_dPe_dNe; mhdReal** d_dPe_px_dNe; mhdReal** d_dPe_py_dNe; mhdReal** d_dPe_pz_dNe;
	mhdReal** d_dJpB_dNe; mhdReal** d_dJpB_px_dNe; mhdReal** d_dJpB_py_dNe; mhdReal** d_dJpB_pz_dNe;
	mhdReal** d_A_dNe; mhdReal** d_A_px_dNe; mhdReal** d_A_py_dNe; mhdReal** d_A_pz_dNe;

	mhdReal** d_Phi_dTe; mhdReal** d_Phi_px_dTe; mhdReal** d_Phi_py_dTe; mhdReal** d_Phi_pz_dTe;
	mhdReal** d_dTe_dTe; mhdReal** d_dTe_px_dTe; mhdReal** d_dTe_py_dTe; mhdReal** d_dTe_pz_dTe;
	mhdReal** d_dNe_dTe; mhdReal** d_dNe_px_dTe; mhdReal** d_dNe_py_dTe; mhdReal** d_dNe_pz_dTe;

	mhdReal** d_Ne0; mhdReal** d_Te0;
	mhdReal** d_Ne0_px; mhdReal** d_Te0_px; mhdReal** d_Pe0_px;

	mhdReal** d_F2perp2;
	mhdReal** d_A2dJpB;
	mhdReal** d_Phi2w;
	mhdReal** d_wdPAdJpB2w;
	mhdReal** d_APhidNe2A;
	mhdReal** d_dPePhiAdJpB2dNe;
	mhdReal** d_PhidTedNe2dTe;

	/*-----------------------------------------Nonlinear------------------------------------------*/

	mhdReal** d_wPhi_w; mhdReal** d_AdJpB_w;
	mhdReal** d_PhiA_A; mhdReal** d_NeA_A;
	mhdReal** d_AdJpB_dNe; mhdReal** d_dNePhi_dNe;
	mhdReal** d_PhiTe_dTe; mhdReal** d_PhiTeA_dTe;

	/*-------------------------------------------------------Matrix on GPU--------------------------------------------------------*/

	int** d_laplacianCsrR;
	int** d_laplacianCsrC;
	mhdReal** d_laplacianCsrV;

	int** d_resistiveCsrR;
	int** d_resistiveCsrC;
	mhdReal** d_resistiveCsrV;

	int** d_PhiCsrR;
	int** d_PhiCsrC;
	mhdReal** d_PhiCsrV;

	int** d_dNeCsrR;
	int** d_dNeCsrC;
	mhdReal** d_dNeCsrV;

	int** d_dTeCsrR;
	int** d_dTeCsrC;
	mhdReal** d_dTeCsrV;

	int** d_dPiCsrR;
	int** d_dPiCsrC;
	mhdReal** d_dPiCsrV;

	int** d_dPaCsrR;
	int** d_dPaCsrC;
	mhdReal** d_dPaCsrV;

	int** d_dPbCsrR;
	int** d_dPbCsrC;
	mhdReal** d_dPbCsrV;

	std::vector<std::vector<cudaStream_t>> cudaStreams;
	std::vector<std::vector<cudssHandle_t>> cudssHandles;

	std::vector<std::vector<cudssConfig_t>> laplacianConfigs;
	std::vector<std::vector<cudssData_t>> laplacianDatas;
	std::vector<std::vector<cudssMatrix_t>> laplacianAs;
	std::vector<std::vector<cudssMatrix_t>> laplacianXs;
	std::vector<std::vector<cudssMatrix_t>> laplacianBs;

	std::vector<std::vector<cudssConfig_t>> resistiveConfigs;
	std::vector<std::vector<cudssData_t>> resistiveDatas;
	std::vector<std::vector<cudssMatrix_t>> resistiveAs;
	std::vector<std::vector<cudssMatrix_t>> resistiveXs;
	std::vector<std::vector<cudssMatrix_t>> resistiveBs;

	std::vector<std::vector<cudssConfig_t>> PhiConfigs;
	std::vector<std::vector<cudssData_t>> PhiDatas;
	std::vector<std::vector<cudssMatrix_t>> PhiAs;
	std::vector<std::vector<cudssMatrix_t>> PhiXs;
	std::vector<std::vector<cudssMatrix_t>> PhiBs;

	std::vector<std::vector<cudssConfig_t>> dNeConfigs;
	std::vector<std::vector<cudssData_t>> dNeDatas;
	std::vector<std::vector<cudssMatrix_t>> dNeAs;
	std::vector<std::vector<cudssMatrix_t>> dNeXs;
	std::vector<std::vector<cudssMatrix_t>> dNeBs;

	std::vector<std::vector<cudssConfig_t>> dTeConfigs;
	std::vector<std::vector<cudssData_t>> dTeDatas;
	std::vector<std::vector<cudssMatrix_t>> dTeAs;
	std::vector<std::vector<cudssMatrix_t>> dTeXs;
	std::vector<std::vector<cudssMatrix_t>> dTeBs;

	std::vector<std::vector<cudssConfig_t>> dPiConfigs;
	std::vector<std::vector<cudssData_t>> dPiDatas;
	std::vector<std::vector<cudssMatrix_t>> dPiAs;
	std::vector<std::vector<cudssMatrix_t>> dPiXs;
	std::vector<std::vector<cudssMatrix_t>> dPiBs;

	std::vector<std::vector<cudssConfig_t>> dPaConfigs;
	std::vector<std::vector<cudssData_t>> dPaDatas;
	std::vector<std::vector<cudssMatrix_t>> dPaAs;
	std::vector<std::vector<cudssMatrix_t>> dPaXs;
	std::vector<std::vector<cudssMatrix_t>> dPaBs;

	std::vector<std::vector<cudssConfig_t>> dPbConfigs;
	std::vector<std::vector<cudssData_t>> dPbDatas;
	std::vector<std::vector<cudssMatrix_t>> dPbAs;
	std::vector<std::vector<cudssMatrix_t>> dPbXs;
	std::vector<std::vector<cudssMatrix_t>> dPbBs;

	/*------------------------------------------------------Particle in Cell------------------------------------------------------*/

	mhdReal** d_pic_Phi_py_A; mhdReal** d_pic_dNe_py_A;
	mhdReal** d_pic_A_A;  mhdReal** d_pic_A_py_A; mhdReal** d_pic_A_pz_A;
	mhdReal** d_pic_APhidNe2A; mhdReal** d_pic_PhiA_A; mhdReal** d_pic_NeA_A;

	picReal IonConst;
	picReal AlphaConst;
	picReal BeamConst;
	picReal IonEPphiLambda[6];
	picReal AlphaEPphiLambda[6];
	picReal BeamEPphiLambda[6];

	mhdReal*** h_globalPi;
	mhdReal*** h_globalPa;
	mhdReal*** h_globalPb;

	mhdReal** d_globalA; mhdReal** d_globalPhi; mhdReal** d_globalApt;
	mhdReal** d_globalPi; mhdReal** d_globalPa; mhdReal** d_globalPb;

	picReal** h_pic1d; picReal** h_pic2d; picReal** h_pic3d;
	picReal** d_pic1d; picReal** d_pic2d; picReal** d_pic3d;

	int** h_Ion_offsets;
	int** h_Ion_keys;
	picReal** h_Ion_values;

	int** h_Alpha_offsets;
	int** h_Alpha_keys;
	picReal** h_Alpha_values;

	int** h_Beam_offsets;
	int** h_Beam_keys;
	picReal** h_Beam_values;

	int** d_Ion_offsets;
	int** d_Ion_keys_in;
	int** d_Ion_keys_out;
	picReal** d_Ion_values_in;
	picReal** d_Ion_values_out;

	int** d_Alpha_offsets;
	int** d_Alpha_keys_in;
	int** d_Alpha_keys_out;
	picReal** d_Alpha_values_in;
	picReal** d_Alpha_values_out;

	int** d_Beam_offsets;
	int** d_Beam_keys_in;
	int** d_Beam_keys_out;
	picReal** d_Beam_values_in;
	picReal** d_Beam_values_out;

	picReal** h_IonPhaseSpaceMapping;
	picReal** h_AlphaPhaseSpaceMapping;
	picReal** h_BeamPhaseSpaceMapping;

	picReal** d_IonPhaseSpaceMapping;
	picReal** d_AlphaPhaseSpaceMapping;
	picReal** d_BeamPhaseSpaceMapping;

	/*-------------------------------------------------------Runtime Config-------------------------------------------------------*/

	const int totalSteps, diagSteps, outputSteps;
	const int leftN, rightN;
	const int nFFTTimeSize, nFFTBatchSize, nFFTFreqSize;
	const int mFFTTimeSize, mFFTBatchSize, mFFTFreqSize;

	const bool ifDiagAmplitude, ifDiagFrequency, ifDiagEparallel, ifDiagDensity, ifDiagDiffusivity, ifDiagZFDrive, ifCheckNAN;
	const bool ifOutputPhi, ifOutputA, ifOutputdNe, ifOutputdTe, ifOutputdPi, ifOutputdPa, ifOutputdPb;

	/*---------------------------------------------------Diagnostic on CPU/GPU----------------------------------------------------*/

	std::vector<mhdReal> h_amplitude; std::vector<mhdReal> h_frequency;
	std::vector<mhdReal> h_modeReal; std::vector<mhdReal> h_modeImag;
	std::vector<mhdReal> h_Epara;    std::vector<mhdReal> h_EparaES;
	std::vector<mhdReal> h_MaxwellDrive; std::vector<mhdReal> h_ReynoldsDrive; std::vector<mhdReal> h_dwdtTotal;
	std::vector<mhdReal> h_IonDensity;   std::vector<mhdReal> h_AlphaDensity;   std::vector<mhdReal> h_BeamDensity;
	std::vector<mhdReal> h_IonDiffusivity; std::vector<mhdReal> h_AlphaDiffusivity; std::vector<mhdReal> h_BeamDiffusivity;

	mhdReal** d_amplitude; mhdReal** d_frequency;
	mhdReal** d_modeReal;  mhdReal** d_modeImag;
	mhdReal** d_Epara;     mhdReal** d_EparaES;
	mhdReal** d_MaxwellDrive; mhdReal** d_ReynoldsDrive; mhdReal** d_dwdtTotal;
	mhdReal** d_Maxwell;      mhdReal** d_Reynolds;
	mhdReal** d_IonDensity;   mhdReal** d_AlphaDensity;   mhdReal** d_BeamDensity;
	mhdReal** d_IonDiffusivity; mhdReal** d_AlphaDiffusivity; mhdReal** d_BeamDiffusivity;

	int** d_NANFlag;

	/*--------------------------------------------------Output Totals on CPU/GPU--------------------------------------------------*/

	std::vector<mhdReal> h_totalPhi; std::vector<mhdReal> h_totalA;
	std::vector<mhdReal> h_totaldNe; std::vector<mhdReal> h_totaldTe;
	std::vector<mhdReal> h_totaldPi; std::vector<mhdReal> h_totaldPa; std::vector<mhdReal> h_totaldPb;
    std::vector<mhdReal> h_IonPhaseDeltaF; std::vector<mhdReal> h_AlphaPhaseDeltaF; std::vector<mhdReal> h_BeamPhaseDeltaF;
    std::vector<mhdReal> h_IonPitchDeltaF; std::vector<mhdReal> h_AlphaPitchDeltaF; std::vector<mhdReal> h_BeamPitchDeltaF;
    std::vector<mhdReal> h_IonPhasePower; std::vector<mhdReal> h_AlphaPhasePower; std::vector<mhdReal> h_BeamPhasePower;
    std::vector<mhdReal> h_IonPitchPower; std::vector<mhdReal> h_AlphaPitchPower; std::vector<mhdReal> h_BeamPitchPower;

	mhdReal** d_totalPhi; mhdReal** d_totalA;
	mhdReal** d_totaldNe; mhdReal** d_totaldTe;
	mhdReal** d_totaldPi; mhdReal** d_totaldPa; mhdReal** d_totaldPb;
    mhdReal** d_IonPhaseDeltaF; mhdReal** d_AlphaPhaseDeltaF; mhdReal** d_BeamPhaseDeltaF;
    mhdReal** d_IonPitchDeltaF; mhdReal** d_AlphaPitchDeltaF; mhdReal** d_BeamPitchDeltaF;
    mhdReal** d_IonPhasePower; mhdReal** d_AlphaPhasePower; mhdReal** d_BeamPhasePower;
    mhdReal** d_IonPitchPower; mhdReal** d_AlphaPitchPower; mhdReal** d_BeamPitchPower;

	/*------------------------------------------------Select Mode NM Buffer on GPU------------------------------------------------*/

	mhdReal** d_Phi_yxz; mhdReal** d_Phi_xzy;
	mhdReal** d_A_yxz;   mhdReal** d_A_xzy;
	mhdReal** d_dNe_yxz; mhdReal** d_dNe_xzy;
	mhdReal** d_dTe_yxz; mhdReal** d_dTe_xzy;
	mhdReal** d_dPi_yxz; mhdReal** d_dPi_xzy;
	mhdReal** d_dPa_yxz; mhdReal** d_dPa_xzy;
	mhdReal** d_dPb_yxz; mhdReal** d_dPb_xzy;

	/*----------------------------------------------cuFFT Plan and Frequency Buffer-----------------------------------------------*/

	std::vector<cufftHandle> d_nPlanR2Cs; std::vector<cufftHandle> d_nPlanC2Rs;
	std::vector<cufftHandle> d_mPlanR2Cs; std::vector<cufftHandle> d_mPlanC2Rs;
	cufftDoubleComplex** d_nFreqd; cufftComplex** d_nFreqf;
	cufftDoubleComplex** d_mFreqd; cufftComplex** d_mFreqf;

	/*---------------------------------------------------CUB Radix Sort Storage---------------------------------------------------*/

	void** d_Ion_storage;   std::vector<size_t> d_Ion_storage_bytes;
	void** d_Alpha_storage; std::vector<size_t> d_Alpha_storage_bytes;
	void** d_Beam_storage;  std::vector<size_t> d_Beam_storage_bytes;

	/*--------------------------------------------------------CUDA Events---------------------------------------------------------*/

	std::vector<cudaEvent_t> startEvents; std::vector<cudaEvent_t> endEvents;
	std::vector<float> elapsedTime;
    // clang-format on
    /*----------------------------------------------------------Function----------------------------------------------------------*/

    void logStart(const char* msg) const {
        if (hostId == 0)
            std::cout << BOLDYELLOW << "Start: " << msg << std::endl;
    }
    void logDone() const {
        if (hostId == 0)
            std::cout << BOLDGREEN << "Done." << RESET << std::endl << std::endl;
    }
    class Rand01 {

      public:
        Rand01() : mt_gen{std::random_device()()}, rnd_dist{0.0, 1.0} {}
        ~Rand01() {}
        double operator()() { return rnd_dist(mt_gen); }

        std::mt19937_64 mt_gen;
        std::uniform_real_distribution<double> rnd_dist;
    };
    class RandNormal {

      public:
        RandNormal() : mt_gen{std::random_device()()}, norm_dist{0.0, 1.0} {}
        ~RandNormal() {}
        double operator()() { return norm_dist(mt_gen); }

        std::mt19937_64 mt_gen;
        std::normal_distribution<double> norm_dist;
    };
    class Allocator {

      public:
        Allocator() {}
        ~Allocator() {}

        /*-----------------------Allocate and Release Arrays on Host and Device-----------------------*/

        template <typename T, typename... Ts>
        void allocateHostArrays(size_t dim1, size_t dim2, T**& hostArray, Ts**&... hostArrays) {

            hostArray = new T*[dim1]();
            hostArray[0] = new T[dim1 * dim2]();

            for (int i = 1; i < dim1; i++)
                hostArray[i] = hostArray[i - 1] + dim2;

            if constexpr (sizeof...(hostArrays) > 0)
                allocateHostArrays(dim1, dim2, hostArrays...);
        }
        template <typename T, typename... Ts>
        void allocateHostArrays(size_t dim1, size_t dim2, size_t dim3, T***& hostArray, Ts***&... hostArrays) {

            hostArray = new T**[dim1]();
            hostArray[0] = new T*[dim1 * dim2]();
            hostArray[0][0] = new T[dim1 * dim2 * dim3]();

            for (int i = 1; i < dim1; i++)
                hostArray[i] = hostArray[i - 1] + dim2;
            for (int i = 1; i < dim1 * dim2; i++)
                hostArray[0][i] = hostArray[0][i - 1] + dim3;

            if constexpr (sizeof...(hostArrays) > 0)
                allocateHostArrays(dim1, dim2, dim3, hostArrays...);
        }
        template <typename T, typename... Ts>
        void releaseHostArrays(T**& hostArray, Ts**&... hostArrays) {

            delete[] hostArray[0];
            delete[] hostArray;

            if constexpr (sizeof...(hostArrays) > 0)
                releaseHostArrays(hostArrays...);
        }
        template <typename T, typename... Ts>
        void releaseHostArrays(T***& hostArray, Ts***&... hostArrays) {

            delete[] hostArray[0][0];
            delete[] hostArray[0];
            delete[] hostArray;

            if constexpr (sizeof...(hostArrays) > 0)
                releaseHostArrays(hostArrays...);
        }

        template <typename T, typename... Ts>
        void allocateDeviceArrays(int localId, int devNums, size_t size, T**& devArray, Ts**&... devArrays) {

            devArray = new T*[devNums]();

            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(cudaMalloc(reinterpret_cast<void**>(&devArray[i]), sizeof(T) * size));
                CUDACHECK(cudaMemsetAsync(devArray[i], 0, sizeof(T) * size, 0));
            }

            if constexpr (sizeof...(devArrays) > 0)
                allocateDeviceArrays(localId, devNums, size, devArrays...);
        }
        template <typename T, typename... Ts>
        void releaseDeviceArrays(int localId, int devNums, T**& devArray, Ts**&... devArrays) {

            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(cudaFree(devArray[i]));
            }
            delete[] devArray;

            if constexpr (sizeof...(devArrays) > 0)
                releaseDeviceArrays(localId, devNums, devArrays...);
        }

        /*----------------------------Memory Copy Between Host and Device-----------------------------*/

        template <typename T, typename... Ts>
        void hostToDevice(size_t size, size_t devOffset, size_t hostOffset, T*& devArray, T*& hostArray,
                          Ts*&... devHostArrays) {

            CUDACHECK(
                cudaMemcpy(devArray + devOffset, hostArray + hostOffset, sizeof(T) * size, cudaMemcpyHostToDevice));
            if constexpr (sizeof...(devHostArrays) > 0)
                hostToDevice(size, devOffset, hostOffset, devHostArrays...);
        }
        template <typename T, typename... Ts>
        void deviceToHost(size_t size, size_t devOffset, size_t hostOffset, T*& devArray, T*& hostArray,
                          Ts*&... devHostArrays) {

            CUDACHECK(
                cudaMemcpy(hostArray + hostOffset, devArray + devOffset, sizeof(T) * size, cudaMemcpyDeviceToHost));
            if constexpr (sizeof...(devHostArrays) > 0)
                deviceToHost(size, devOffset, hostOffset, devHostArrays...);
        }

        /*--------------------------------Help Initialize Equilibrium---------------------------------*/

        template <typename T, typename... Ts>
        void binaryToHost(size_t offset, size_t dim1, size_t dim2, std::vector<double>& binary, T**& hostArray,
                          Ts**&... hostArrays) {

            for (int i = 0; i < dim1; i++) {
                for (int j = 0; j < dim2; j++) {
                    hostArray[i][j] = binary[offset + i * dim2 + j];
                }
            }

            if constexpr (sizeof...(hostArrays) > 0)
                binaryToHost(offset + dim1 * dim2, dim1, dim2, binary, hostArrays...);
        }
    };

    void setup() {

        /*-------------------------------------Output Directories-------------------------------------*/

        if (hostId == 0) {
            std::filesystem::create_directories(initialDir);
            std::filesystem::create_directories(finalDir);
        }

        /*-------------------------------------Memory Allocation--------------------------------------*/

        allocateHostMemory();
        allocateDeviceMemory();

        /*------------------------------MHD Equilibrium and Perturbation------------------------------*/

        loadMHDEquilibrium(inputDir + "/" + MHDCollocated);
        if constexpr (std::is_same_v<ifContinue, trueType>)
            loadMHDPerturbation(inputDir + "/MHDContinue_" + std::to_string(continueSteps) + ".bin");
        else
            computeMHDPerturbation<perturbLeftN, perturbRightN>(perturbRadialIndex, perturbWidth, perturbAmplitude,
                                                                initialDir);

        /*------------------------------------Input Files Snapshot------------------------------------*/

        if (hostId == 0) {
            namespace fs = std::filesystem;
            fs::copy_file(inputDir + "/" + MHDCollocated, finalDir + "/" + MHDCollocated,
                          fs::copy_options::overwrite_existing);
            if constexpr (std::is_same_v<ifStaggered, trueType>)
                fs::copy_file(inputDir + "/" + MHDStaggered, finalDir + "/" + MHDStaggered,
                              fs::copy_options::overwrite_existing);

            if constexpr (std::is_same_v<ifContinue, falseType> || continueSteps == 0) {
                fs::copy_file(inputDir + "/" + MHDCollocated, initialDir + "/" + MHDCollocated,
                              fs::copy_options::overwrite_existing);
                if constexpr (std::is_same_v<ifStaggered, trueType>)
                    fs::copy_file(inputDir + "/" + MHDStaggered, initialDir + "/" + MHDStaggered,
                                  fs::copy_options::overwrite_existing);
            }
            if constexpr (std::is_same_v<ifContinue, trueType> && continueSteps == 0)
                fs::copy_file(inputDir + "/MHDContinue_0.bin", initialDir + "/MHDContinue_0.bin",
                              fs::copy_options::overwrite_existing);
        }

        /*-----------------------------Sparse Matrices on Collocated Grid-----------------------------*/

        compressCollocatedCoefficient();
        computeSparseMatrix<Laplacian, trueType, trueType>();
        if constexpr (std::is_same_v<ifNablaPerp2A, trueType> && std::is_same_v<ifStaggered, falseType>)
            computeSparseMatrix<Resistive, trueType, trueType>();
        if constexpr (std::is_same_v<ifNablaPerp2Phi, trueType>)
            computeSparseMatrix<Perp2Phi, trueType, trueType>();
        if constexpr (std::is_same_v<ifNablaPerp2dNe, trueType>)
            computeSparseMatrix<Perp2dNe, trueType, trueType>();
        if constexpr (std::is_same_v<ifNablaPerp2dTe, trueType>)
            computeSparseMatrix<Perp2dTe, trueType, trueType>();
        if constexpr (std::is_same_v<ifNablaPerp2dPi, trueType>)
            computeSparseMatrix<Perp2dPi, trueType, trueType>();
        if constexpr (std::is_same_v<ifNablaPerp2dPa, trueType>)
            computeSparseMatrix<Perp2dPa, trueType, trueType>();
        if constexpr (std::is_same_v<ifNablaPerp2dPb, trueType>)
            computeSparseMatrix<Perp2dPb, trueType, trueType>();

        /*----------------------------------Phase-Space Diagnostics-----------------------------------*/

        if constexpr (std::is_same_v<ifContinue, falseType>) {

            if constexpr (std::is_same_v<ifOutputPhaceSpaceF0, trueType>) {
                if constexpr (std::is_same_v<::ifIon, trueType>)
                    computePhaseSpaceF0<Ion, IonType, IonSpace, IonVelocity, gridE, gridPphi, gridLambda, ppcPhase>(
                        initialDir);
                if constexpr (std::is_same_v<::ifAlpha, trueType>)
                    computePhaseSpaceF0<Alpha, AlphaType, AlphaSpace, AlphaVelocity, gridE, gridPphi, gridLambda,
                                        ppcPhase>(initialDir);
                if constexpr (std::is_same_v<::ifBeam, trueType>)
                    computePhaseSpaceF0<Beam, BeamType, BeamSpace, BeamVelocity, gridE, gridPphi, gridLambda, ppcPhase>(
                        initialDir);
            }

            if constexpr (std::is_same_v<ifOutputPitchSpaceF0, trueType>) {
                if constexpr (std::is_same_v<::ifIon, trueType>)
                    computePitchSpaceF0<Ion, IonType, IonSpace, IonVelocity, gridVpara, gridVperp, ppcPitch>(
                        initialDir);
                if constexpr (std::is_same_v<::ifAlpha, trueType>)
                    computePitchSpaceF0<Alpha, AlphaType, AlphaSpace, AlphaVelocity, gridVpara, gridVperp, ppcPitch>(
                        initialDir);
                if constexpr (std::is_same_v<::ifBeam, trueType>)
                    computePitchSpaceF0<Beam, BeamType, BeamSpace, BeamVelocity, gridVpara, gridVperp, ppcPitch>(
                        initialDir);
            }

            if constexpr (std::is_same_v<ifOutputPhaceSpaceJacobian, trueType>) {
                if constexpr (std::is_same_v<::ifIon, trueType>)
                    computePhaseSpaceJacobian<Ion, gridE, gridPphi, gridLambda, ppcPhase>(initialDir);
                if constexpr (std::is_same_v<::ifAlpha, trueType>)
                    computePhaseSpaceJacobian<Alpha, gridE, gridPphi, gridLambda, ppcPhase>(initialDir);
                if constexpr (std::is_same_v<::ifBeam, trueType>)
                    computePhaseSpaceJacobian<Beam, gridE, gridPphi, gridLambda, ppcPhase>(initialDir);
            }

            if constexpr (std::is_same_v<ifOutputPitchSpaceJacobian, trueType>) {
                if constexpr (std::is_same_v<::ifIon, trueType>)
                    computePitchSpaceJacobian<Ion, gridVpara, gridVperp, ppcPitch>(initialDir);
                if constexpr (std::is_same_v<::ifAlpha, trueType>)
                    computePitchSpaceJacobian<Alpha, gridVpara, gridVperp, ppcPitch>(initialDir);
                if constexpr (std::is_same_v<::ifBeam, trueType>)
                    computePitchSpaceJacobian<Beam, gridVpara, gridVperp, ppcPitch>(initialDir);
            }

            if constexpr (std::is_same_v<ifOutputPhaceSpaceFrequency, trueType>) {

                auto diagPhaseOrbit = [&]<picType species, typename... Guards>(picReal** d_mapping, picReal** h_mapping,
                                                                               const std::string& mappingFile) {
                    if constexpr (allTrue<Guards...>) {
                        const char* displayName = "";
                        if constexpr (species == Ion)
                            displayName = "thermal ions";
                        if constexpr (species == Alpha)
                            displayName = "alpha particles";
                        if constexpr (species == Beam)
                            displayName = "beam particles";

                        loadPhaseSpaceMapping<species>(inputDir + "/" + mappingFile);
                        if (hostId == 0)
                            std::cout << BOLDYELLOW << "Start: Compute equilibrium orbit of " << displayName << "."
                                      << RESET << std::endl;

                        for (int idx = 0; idx < 20000; idx++) {
                            for (int i = 0; i < devNums; i++) {
                                CUDACHECK(cudaSetDevice(localId * devNums + i));
                                PICDiagOrbit<1, species, picReal>
                                    <<<gridE * gridPphi * gridLambda * 2 / pptNums / PICBlockDimx, PICBlockDimx, 0,
                                       0>>>(d_pic1d[i], d_pic2d[i], d_mapping[i]);
                            }
                        }

                        if (hostId == 0) {
                            std::cout << BOLDGREEN << "Done." << RESET << std::endl;
                            std::cout << std::endl;
                        }
                        for (int i = 0; i < devNums; i++) {
                            CUDACHECK(cudaSetDevice(localId * devNums + i));
                            cudaMemcpy(h_mapping[i], d_mapping[i],
                                       sizeof(picReal) * gridE * gridPphi * gridLambda * 2 * 13,
                                       cudaMemcpyDeviceToHost);
                        }
                        computePhaseSpaceFrequency<species>(inputDir + "/" + mappingFile, initialDir);
                    }
                };

                diagPhaseOrbit.template operator()<Ion, ::ifIon>(d_IonPhaseSpaceMapping, h_IonPhaseSpaceMapping,
                                                                 IonPhaseSpaceMapping);
                diagPhaseOrbit.template operator()<Alpha, ::ifAlpha>(d_AlphaPhaseSpaceMapping, h_AlphaPhaseSpaceMapping,
                                                                     AlphaPhaseSpaceMapping);
                diagPhaseOrbit.template operator()<Beam, ::ifBeam>(d_BeamPhaseSpaceMapping, h_BeamPhaseSpaceMapping,
                                                                   BeamPhaseSpaceMapping);
            }
        }

        /*----------------------------------Particle Initialization-----------------------------------*/

        if constexpr (std::is_same_v<ifContinue, trueType>) {
            loadParticles(inputDir + "/PICContinue_" + std::to_string(continueSteps) + ".bin");
            if constexpr (continueSteps == 0) {
                if constexpr (std::is_same_v<::ifIon, trueType>)
                    computeEquilibriumPressure<Ion>(initialDir);
                if constexpr (std::is_same_v<::ifAlpha, trueType>)
                    computeEquilibriumPressure<Alpha>(initialDir);
                if constexpr (std::is_same_v<::ifBeam, trueType>)
                    computeEquilibriumPressure<Beam>(initialDir);
            }
        } else {
            if constexpr (std::is_same_v<::ifIon, trueType>) {
                loadParticles<Ion, IonType, IonSpace, IonVelocity>();
                computeEquilibriumPressure<Ion>(initialDir);
            }
            if constexpr (std::is_same_v<::ifAlpha, trueType>) {
                loadParticles<Alpha, AlphaType, AlphaSpace, AlphaVelocity>();
                computeEquilibriumPressure<Alpha>(initialDir);
            }
            if constexpr (std::is_same_v<::ifBeam, trueType>) {
                loadParticles<Beam, BeamType, BeamSpace, BeamVelocity>();
                computeEquilibriumPressure<Beam>(initialDir);
            }
        }

        if constexpr (std::is_same_v<::ifIon, trueType>)
            computePhaseSpaceRange<Ion>();
        if constexpr (std::is_same_v<::ifAlpha, trueType>)
            computePhaseSpaceRange<Alpha>();
        if constexpr (std::is_same_v<::ifBeam, trueType>)
            computePhaseSpaceRange<Beam>();

        /*----------------------------------Constant Symbols Upload-----------------------------------*/

        for (int i = 0; i < devNums; i++) {
            CUDACHECK(cudaSetDevice(localId * devNums + i));
            if constexpr (std::is_same_v<::ifIon, trueType>)
                CUDACHECK(cudaMemcpyToSymbol(::IonConst, &IonConst, sizeof(picReal)));
            if constexpr (std::is_same_v<::ifAlpha, trueType>)
                CUDACHECK(cudaMemcpyToSymbol(::AlphaConst, &AlphaConst, sizeof(picReal)));
            if constexpr (std::is_same_v<::ifBeam, trueType>)
                CUDACHECK(cudaMemcpyToSymbol(::BeamConst, &BeamConst, sizeof(picReal)));
            if constexpr (std::is_same_v<::ifIon, trueType>)
                CUDACHECK(cudaMemcpyToSymbol(::IonEPphiLambda, IonEPphiLambda, sizeof(IonEPphiLambda)));
            if constexpr (std::is_same_v<::ifAlpha, trueType>)
                CUDACHECK(cudaMemcpyToSymbol(::AlphaEPphiLambda, AlphaEPphiLambda, sizeof(AlphaEPphiLambda)));
            if constexpr (std::is_same_v<::ifBeam, trueType>)
                CUDACHECK(cudaMemcpyToSymbol(::BeamEPphiLambda, BeamEPphiLambda, sizeof(BeamEPphiLambda)));
        }

        /*-----------------------------------Host to Device Upload------------------------------------*/

        memcpyHostToDevice();

        /*------------------------------Sparse Matrix on Staggered Grid-------------------------------*/

        if constexpr (std::is_same_v<ifStaggered, trueType>) {
            loadMHDEquilibrium(inputDir + "/" + MHDStaggered);
            compressStaggeredCoefficient();
            if constexpr (std::is_same_v<ifNablaPerp2A, trueType>)
                computeSparseMatrix<Resistive, trueType, trueType>();
        }
    }

    void allocateHostMemory() {

        logStart("Allocate host memory.");

        Allocator HostAllocator;

        /*-----------------------------------MHD Equilibrium on CPU-----------------------------------*/

        HostAllocator.allocateHostArrays(
            gridNx, gridNy, q, q_px, psip, psip_px, Ni, Ni_px, Ti, Ti_px, Pi, Pi_px, Ne, Ne_px, Te, Te_px, Pe, Pe_px,
            Na, Na_px, Ta, Ta_px, Nb, Nb_px, Tb, Tb_px, B, B_px, B_py, B_px2, B_pxy, B_py2, J, J_px, J_py, Bny, Bny_px,
            Bny_py, Va, Va_px, Va_py, Rho, Rho_px, Rho_py, JpB, JpB_px, JpB_py, R, Z, gconxx, gconxx_px, gconxx_py,
            gconxy, gconxy_px, gconxy_py, gconxz, gconxz_px, gconxz_py, gconyy, gconyy_px, gconyy_py, gconyz, gconyz_px,
            gconyz_py, gconzz, gconzz_px, gconzz_py, gcovxx, gcovxx_px, gcovxx_py, gcovxy, gcovxy_px, gcovxy_py, gcovxz,
            gcovxz_px, gcovxz_py, gcovyy, gcovyy_px, gcovyy_py, gcovyz, gcovyz_px, gcovyz_py, gcovzz, gcovzz_px,
            gcovzz_py);

        HostAllocator.allocateHostArrays(
            gridNx, gridNyPlusGhost, SFAconxx, SFAconxx_px, SFAconxx_py, SFAconxy, SFAconxy_px, SFAconxy_py, SFAconxz,
            SFAconxz_px, SFAconxz_py, SFAconyy, SFAconyy_px, SFAconyy_py, SFAconyz, SFAconyz_px, SFAconyz_py, SFAconzz,
            SFAconzz_px, SFAconzz_py, SFAcovxx, SFAcovxx_px, SFAcovxx_py, SFAcovxy, SFAcovxy_px, SFAcovxy_py, SFAcovxz,
            SFAcovxz_px, SFAcovxz_py, SFAcovyy, SFAcovyy_px, SFAcovyy_py, SFAcovyz, SFAcovyz_px, SFAcovyz_py, SFAcovzz,
            SFAcovzz_px, SFAcovzz_py);

        /*----------------------------------MHD Perturbation on CPU-----------------------------------*/

        HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx, h_qtheta);
        HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx, gridNz, h_w, h_A, h_dNe, h_dTe, h_Phi, h_dJpB, h_dPe);

        /*-------------------------------Coefficient Compression on CPU-------------------------------*/

        /*---------------------------Linear---------------------------*/

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_A_w, h_A_px_w, h_A_py_w, h_A_pz_w, h_dJpB_w, h_dJpB_px_w,
                                         h_dJpB_py_w, h_dJpB_pz_w, h_dP_w, h_dP_px_w, h_dP_py_w, h_dP_pz_w, h_w_py_w,
                                         h_w_pz_w, h_w_Phi);

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_Phi_w, h_Phi_px_w, h_Phi_pz_w, h_Phi_px2_w, h_Phi_pxz_w,
                                         h_Phi_pz2_w);

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_A_resistive, h_A_px_resistive, h_A_pz_resistive,
                                         h_A_px2_resistive, h_A_pxz_resistive, h_A_pz2_resistive);

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_F_perp2, h_F_px_perp2, h_F_pz_perp2, h_F_px2_perp2,
                                         h_F_pxz_perp2, h_F_pz2_perp2);

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_A_dJpB, h_A_px_dJpB, h_A_pz_dJpB, h_A_px2_dJpB, h_A_pxz_dJpB,
                                         h_A_pz2_dJpB);

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_Phi_A, h_Phi_px_A, h_Phi_py_A, h_Phi_pz_A, h_dNe_A,
                                         h_dNe_px_A, h_dNe_py_A, h_dNe_pz_A, h_A_A, h_A_px_A, h_A_py_A, h_A_pz_A);

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_Phi_dNe, h_Phi_px_dNe, h_Phi_py_dNe, h_Phi_pz_dNe, h_dPe_dNe,
                                         h_dPe_px_dNe, h_dPe_py_dNe, h_dPe_pz_dNe, h_dJpB_dNe, h_dJpB_px_dNe,
                                         h_dJpB_py_dNe, h_dJpB_pz_dNe, h_A_dNe, h_A_px_dNe, h_A_py_dNe, h_A_pz_dNe);

        HostAllocator.allocateHostArrays(gridNy, gridNx, h_Phi_dTe, h_Phi_px_dTe, h_Phi_py_dTe, h_Phi_pz_dTe, h_dTe_dTe,
                                         h_dTe_px_dTe, h_dTe_py_dTe, h_dTe_pz_dTe, h_dNe_dTe, h_dNe_px_dTe,
                                         h_dNe_py_dTe, h_dNe_pz_dTe, h_Ne0, h_Te0, h_Ne0_px, h_Te0_px, h_Pe0_px);

        HostAllocator.allocateHostArrays(gridNy, gridNx, 5, h_F2perp2);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_A2dJpB);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 5, h_Phi2w);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 10, h_wdPAdJpB2w);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 5, h_APhidNe2A);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 11, h_dPePhiAdJpB2dNe);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_PhidTedNe2dTe);

        /*-------------------------Nonlinear--------------------------*/

        HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_wPhi_w);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 9, h_AdJpB_w);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 9, h_PhiA_A, h_NeA_A);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 9, h_AdJpB_dNe, h_dNePhi_dNe);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_PhiTe_dTe);
        HostAllocator.allocateHostArrays(gridNy, gridNx, 18, h_PhiTeA_dTe);

        /*--------------------------------------Particle in Cell--------------------------------------*/

        if (ifIon) {

            HostAllocator.allocateHostArrays(devNums, 8, h_Ion_offsets);
            HostAllocator.allocateHostArrays(devNums, (size_t)picDev * 7, h_Ion_keys);
            HostAllocator.allocateHostArrays(devNums, (size_t)picDev * 7, h_Ion_values);
        }

        if (ifAlpha) {

            HostAllocator.allocateHostArrays(devNums, 8, h_Alpha_offsets);
            HostAllocator.allocateHostArrays(devNums, (size_t)picDev * 7, h_Alpha_keys);
            HostAllocator.allocateHostArrays(devNums, (size_t)picDev * 7, h_Alpha_values);
        }

        if (ifBeam) {

            HostAllocator.allocateHostArrays(devNums, 8, h_Beam_offsets);
            HostAllocator.allocateHostArrays(devNums, (size_t)picDev * 7, h_Beam_keys);
            HostAllocator.allocateHostArrays(devNums, (size_t)picDev * 7, h_Beam_values);
        }

        HostAllocator.allocateHostArrays(cellNx, 30, h_pic1d);
        HostAllocator.allocateHostArrays(cellNy * cellNx, 72, h_pic2d);
        HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx * gridNzPlusGhost * 8, h_pic3d);

        HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx, gridNz, h_globalPi, h_globalPa, h_globalPb);

        /*-------------------------------------Diagnostic on CPU--------------------------------------*/

        const size_t diagLen = (size_t)(totalSteps / diagSteps + 1) * gridNx;
        const size_t diagLenN = diagLen * (rightN - leftN + 1);

        h_amplitude.resize(diagLenN);
        h_modeReal.resize(diagLenN);
        h_modeImag.resize(diagLenN);
        h_frequency.resize(diagLen);
        h_Epara.resize(diagLen);
        h_EparaES.resize(diagLen);
        h_MaxwellDrive.resize(diagLen);
        h_ReynoldsDrive.resize(diagLen);
        h_dwdtTotal.resize(diagLen);
        h_IonDensity.resize(diagLen);
        h_AlphaDensity.resize(diagLen);
        h_BeamDensity.resize(diagLen);
        h_IonDiffusivity.resize(diagLenN);
        h_AlphaDiffusivity.resize(diagLenN);
        h_BeamDiffusivity.resize(diagLenN);

        /*------------------------------------Output Totals on CPU------------------------------------*/

        const size_t outputLen = (size_t)(totalSteps / outputSteps + 1) * gridNy * gridNxz;

        h_totalPhi.resize(outputLen);
        h_totalA.resize(outputLen);
        h_totaldNe.resize(outputLen);
        h_totaldTe.resize(outputLen);
        h_totaldPi.resize(outputLen);
        h_totaldPa.resize(outputLen);
        h_totaldPb.resize(outputLen);

        const size_t phaseDeltaFLen = (size_t)(totalSteps / outputSteps + 1) * gridE * gridPphi * gridLambda;

        h_IonPhaseDeltaF.resize(phaseDeltaFLen);
        h_AlphaPhaseDeltaF.resize(phaseDeltaFLen);
        h_BeamPhaseDeltaF.resize(phaseDeltaFLen);

        const size_t phasePowerLen =
            (size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridE * gridPphi * gridLambda;

        h_IonPhasePower.resize(phasePowerLen);
        h_AlphaPhasePower.resize(phasePowerLen);
        h_BeamPhasePower.resize(phasePowerLen);

        const size_t pitchDeltaFLen = (size_t)(totalSteps / outputSteps + 1) * gridVpara * gridVperp;

        h_IonPitchDeltaF.resize(pitchDeltaFLen);
        h_AlphaPitchDeltaF.resize(pitchDeltaFLen);
        h_BeamPitchDeltaF.resize(pitchDeltaFLen);

        const size_t pitchPowerLen =
            (size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridVpara * gridVperp;

        h_IonPitchPower.resize(pitchPowerLen);
        h_AlphaPitchPower.resize(pitchPowerLen);
        h_BeamPitchPower.resize(pitchPowerLen);

        /*-------------------------------------CUB Storage Bytes--------------------------------------*/

        d_Ion_storage_bytes.resize(devNums);
        d_Alpha_storage_bytes.resize(devNums);
        d_Beam_storage_bytes.resize(devNums);

        /*-----------------------------cuFFT Plan Handle / Event Vectors------------------------------*/

        d_nPlanR2Cs.resize(devNums);
        d_nPlanC2Rs.resize(devNums);
        d_mPlanR2Cs.resize(devNums);
        d_mPlanC2Rs.resize(devNums);
        startEvents.resize(devNums);
        endEvents.resize(devNums);
        elapsedTime.resize(devNums);

        if (hostId == 0) {
            logDone();
        }
    }
    void releaseHostMemory() {

        logStart("Release host memory.");

        Allocator HostAllocator;

        /*-----------------------------------MHD Equilibrium on CPU-----------------------------------*/

        HostAllocator.releaseHostArrays(
            q, q_px, psip, psip_px, Ni, Ni_px, Ti, Ti_px, Pi, Pi_px, Ne, Ne_px, Te, Te_px, Pe, Pe_px, Na, Na_px, Ta,
            Ta_px, Nb, Nb_px, Tb, Tb_px, B, B_px, B_py, B_px2, B_pxy, B_py2, J, J_px, J_py, Bny, Bny_px, Bny_py, Va,
            Va_px, Va_py, Rho, Rho_px, Rho_py, JpB, JpB_px, JpB_py, R, Z, gconxx, gconxx_px, gconxx_py, gconxy,
            gconxy_px, gconxy_py, gconxz, gconxz_px, gconxz_py, gconyy, gconyy_px, gconyy_py, gconyz, gconyz_px,
            gconyz_py, gconzz, gconzz_px, gconzz_py, gcovxx, gcovxx_px, gcovxx_py, gcovxy, gcovxy_px, gcovxy_py, gcovxz,
            gcovxz_px, gcovxz_py, gcovyy, gcovyy_px, gcovyy_py, gcovyz, gcovyz_px, gcovyz_py, gcovzz, gcovzz_px,
            gcovzz_py);

        HostAllocator.releaseHostArrays(
            SFAconxx, SFAconxx_px, SFAconxx_py, SFAconxy, SFAconxy_px, SFAconxy_py, SFAconxz, SFAconxz_px, SFAconxz_py,
            SFAconyy, SFAconyy_px, SFAconyy_py, SFAconyz, SFAconyz_px, SFAconyz_py, SFAconzz, SFAconzz_px, SFAconzz_py,
            SFAcovxx, SFAcovxx_px, SFAcovxx_py, SFAcovxy, SFAcovxy_px, SFAcovxy_py, SFAcovxz, SFAcovxz_px, SFAcovxz_py,
            SFAcovyy, SFAcovyy_px, SFAcovyy_py, SFAcovyz, SFAcovyz_px, SFAcovyz_py, SFAcovzz, SFAcovzz_px, SFAcovzz_py);

        /*----------------------------------MHD Perturbation on CPU-----------------------------------*/

        HostAllocator.releaseHostArrays(h_qtheta);
        HostAllocator.releaseHostArrays(h_w, h_A, h_dNe, h_dTe, h_Phi, h_dJpB, h_dPe);

        /*-------------------------------Coefficient Compression on CPU-------------------------------*/

        /*---------------------------Linear---------------------------*/

        HostAllocator.releaseHostArrays(h_A_w, h_A_px_w, h_A_py_w, h_A_pz_w, h_dJpB_w, h_dJpB_px_w, h_dJpB_py_w,
                                        h_dJpB_pz_w, h_dP_w, h_dP_px_w, h_dP_py_w, h_dP_pz_w, h_w_py_w, h_w_pz_w,
                                        h_w_Phi);

        HostAllocator.releaseHostArrays(h_Phi_w, h_Phi_px_w, h_Phi_pz_w, h_Phi_px2_w, h_Phi_pxz_w, h_Phi_pz2_w);

        HostAllocator.releaseHostArrays(h_A_resistive, h_A_px_resistive, h_A_pz_resistive, h_A_px2_resistive,
                                        h_A_pxz_resistive, h_A_pz2_resistive);

        HostAllocator.releaseHostArrays(h_F_perp2, h_F_px_perp2, h_F_pz_perp2, h_F_px2_perp2, h_F_pxz_perp2,
                                        h_F_pz2_perp2);

        HostAllocator.releaseHostArrays(h_A_dJpB, h_A_px_dJpB, h_A_pz_dJpB, h_A_px2_dJpB, h_A_pxz_dJpB, h_A_pz2_dJpB);

        HostAllocator.releaseHostArrays(h_Phi_A, h_Phi_px_A, h_Phi_py_A, h_Phi_pz_A, h_dNe_A, h_dNe_px_A, h_dNe_py_A,
                                        h_dNe_pz_A, h_A_A, h_A_px_A, h_A_py_A, h_A_pz_A);

        HostAllocator.releaseHostArrays(h_Phi_dNe, h_Phi_px_dNe, h_Phi_py_dNe, h_Phi_pz_dNe, h_dPe_dNe, h_dPe_px_dNe,
                                        h_dPe_py_dNe, h_dPe_pz_dNe, h_dJpB_dNe, h_dJpB_px_dNe, h_dJpB_py_dNe,
                                        h_dJpB_pz_dNe, h_A_dNe, h_A_px_dNe, h_A_py_dNe, h_A_pz_dNe);

        HostAllocator.releaseHostArrays(h_Phi_dTe, h_Phi_px_dTe, h_Phi_py_dTe, h_Phi_pz_dTe, h_dTe_dTe, h_dTe_px_dTe,
                                        h_dTe_py_dTe, h_dTe_pz_dTe, h_dNe_dTe, h_dNe_px_dTe, h_dNe_py_dTe, h_dNe_pz_dTe,
                                        h_Ne0, h_Te0, h_Ne0_px, h_Te0_px, h_Pe0_px);

        HostAllocator.releaseHostArrays(h_F2perp2);
        HostAllocator.releaseHostArrays(h_A2dJpB);
        HostAllocator.releaseHostArrays(h_Phi2w);
        HostAllocator.releaseHostArrays(h_wdPAdJpB2w);
        HostAllocator.releaseHostArrays(h_APhidNe2A);
        HostAllocator.releaseHostArrays(h_dPePhiAdJpB2dNe);
        HostAllocator.releaseHostArrays(h_PhidTedNe2dTe);

        /*-------------------------Nonlinear--------------------------*/

        HostAllocator.releaseHostArrays(h_wPhi_w, h_AdJpB_w, h_PhiA_A, h_NeA_A, h_AdJpB_dNe, h_dNePhi_dNe, h_PhiTe_dTe,
                                        h_PhiTeA_dTe);

        /*--------------------------------------Particle in Cell--------------------------------------*/

        if (ifIon) {

            HostAllocator.releaseHostArrays(h_Ion_offsets);
            HostAllocator.releaseHostArrays(h_Ion_keys);
            HostAllocator.releaseHostArrays(h_Ion_values);
        }

        if (ifAlpha) {

            HostAllocator.releaseHostArrays(h_Alpha_offsets);
            HostAllocator.releaseHostArrays(h_Alpha_keys);
            HostAllocator.releaseHostArrays(h_Alpha_values);
        }

        if (ifBeam) {

            HostAllocator.releaseHostArrays(h_Beam_offsets);
            HostAllocator.releaseHostArrays(h_Beam_keys);
            HostAllocator.releaseHostArrays(h_Beam_values);
        }

        HostAllocator.releaseHostArrays(h_pic1d, h_pic2d, h_pic3d);
        HostAllocator.releaseHostArrays(h_globalPi, h_globalPa, h_globalPb);

        if (hostId == 0) {
            logDone();
        }
    }
    void allocateDeviceMemory() {

        logStart("Allocate device memory.");

        Allocator DeviceAllocator;

        /*----------------------------------MHD Perturbation on GPU-----------------------------------*/

        DeviceAllocator.allocateDeviceArrays(localId, devNums, (devNy + 2 * gridGhost) * gridNx, d_qtheta);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, gridNy * gridNxz, d_w, d_A, d_dNe, d_dTe, d_Phi, d_dJpB,
                                             d_dPe);

        DeviceAllocator.allocateDeviceArrays(
            localId, devNums, (devNy + 2 * gridGhost) * gridNxz, d_w_beg, d_w_midl, d_w_midr, d_w_end, d_A_beg,
            d_A_midl, d_A_midr, d_A_end, d_dNe_beg, d_dNe_midl, d_dNe_midr, d_dNe_end, d_dTe_beg, d_dTe_midl,
            d_dTe_midr, d_dTe_end, d_Phi_midl, d_Phi_midr, d_dJpB_midl, d_dJpB_midr, d_dPe_midl, d_dPe_midr, d_dPi_midl,
            d_dPi_midr, d_dPa_midl, d_dPa_midr, d_dPb_midl, d_dPb_midr, d_Apt_midl, d_Apt_midr);

        /*-------------------------------Coefficient Compression on GPU-------------------------------*/

        /*---------------------------Linear---------------------------*/

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx, d_A_w, d_A_px_w, d_A_py_w, d_A_pz_w,
                                             d_dJpB_w, d_dJpB_px_w, d_dJpB_py_w, d_dJpB_pz_w, d_dP_w, d_dP_px_w,
                                             d_dP_py_w, d_dP_pz_w, d_w_py_w, d_w_pz_w, d_w2Phi);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx, d_F_perp2, d_F_px_perp2, d_F_pz_perp2,
                                             d_F_px2_perp2, d_F_pxz_perp2, d_F_pz2_perp2);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx, d_A_dJpB, d_A_px_dJpB, d_A_pz_dJpB,
                                             d_A_px2_dJpB, d_A_pxz_dJpB, d_A_pz2_dJpB);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx, d_Phi_A, d_Phi_px_A, d_Phi_py_A,
                                             d_Phi_pz_A, d_dNe_A, d_dNe_px_A, d_dNe_py_A, d_dNe_pz_A, d_A_A, d_A_px_A,
                                             d_A_py_A, d_A_pz_A);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx, d_Phi_dNe, d_Phi_px_dNe, d_Phi_py_dNe,
                                             d_Phi_pz_dNe, d_dPe_dNe, d_dPe_px_dNe, d_dPe_py_dNe, d_dPe_pz_dNe,
                                             d_dJpB_dNe, d_dJpB_px_dNe, d_dJpB_py_dNe, d_dJpB_pz_dNe, d_A_dNe,
                                             d_A_px_dNe, d_A_py_dNe, d_A_pz_dNe);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx, d_Phi_dTe, d_Phi_px_dTe, d_Phi_py_dTe,
                                             d_Phi_pz_dTe, d_dTe_dTe, d_dTe_px_dTe, d_dTe_py_dTe, d_dTe_pz_dTe,
                                             d_dNe_dTe, d_dNe_px_dTe, d_dNe_py_dTe, d_dNe_pz_dTe, d_Ne0, d_Te0,
                                             d_Ne0_px, d_Te0_px, d_Pe0_px);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, d_F2perp2);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_A2dJpB);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, d_Phi2w);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 10, d_wdPAdJpB2w);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, d_APhidNe2A);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 11, d_dPePhiAdJpB2dNe);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_PhidTedNe2dTe);

        /*-------------------------Nonlinear--------------------------*/

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_wPhi_w);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, d_AdJpB_w);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, d_PhiA_A, d_NeA_A);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, d_AdJpB_dNe, d_dNePhi_dNe);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_PhiTe_dTe);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 18, d_PhiTeA_dTe);

        /*--------------------------------------Particle in Cell--------------------------------------*/

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx, d_pic_Phi_py_A, d_pic_dNe_py_A,
                                             d_pic_A_A, d_pic_A_py_A, d_pic_A_pz_A);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, d_pic_APhidNe2A);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, d_pic_PhiA_A, d_pic_NeA_A);

        if (ifIon) {

            DeviceAllocator.allocateDeviceArrays(localId, devNums, 8, d_Ion_offsets);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Ion_keys_in, d_Ion_keys_out);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Ion_values_in,
                                                 d_Ion_values_out);
        }

        if (ifAlpha) {

            DeviceAllocator.allocateDeviceArrays(localId, devNums, 8, d_Alpha_offsets);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Alpha_keys_in,
                                                 d_Alpha_keys_out);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Alpha_values_in,
                                                 d_Alpha_values_out);
        }

        if (ifBeam) {

            DeviceAllocator.allocateDeviceArrays(localId, devNums, 8, d_Beam_offsets);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Beam_keys_in, d_Beam_keys_out);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Beam_values_in,
                                                 d_Beam_values_out);
        }

        DeviceAllocator.allocateDeviceArrays(localId, devNums, cellNx * 30, d_pic1d);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, cellNy * cellNx * 72, d_pic2d);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, gridNyPlusGhost * gridNx * gridNzPlusGhost * 8, d_pic3d);

        DeviceAllocator.allocateDeviceArrays(localId, devNums, gridNyPlusGhost * gridNxz, d_globalA, d_globalPhi,
                                             d_globalApt, d_globalPi, d_globalPa, d_globalPb);

        /*-------------------------------------Diagnostic on GPU--------------------------------------*/

        if constexpr (std::is_same_v<::ifDiagAmplitude, trueType>) {
            DeviceAllocator.allocateDeviceArrays(localId, devNums, h_amplitude.size(), d_amplitude);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, h_modeReal.size(), d_modeReal, d_modeImag);
        }

        if constexpr (std::is_same_v<::ifDiagFrequency, trueType>)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, h_frequency.size(), d_frequency);

        if constexpr (std::is_same_v<::ifDiagEparallel, trueType>)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, h_Epara.size(), d_Epara, d_EparaES);

        if constexpr (std::is_same_v<::ifDiagZFDrive, trueType>) {
            DeviceAllocator.allocateDeviceArrays(localId, devNums, h_MaxwellDrive.size(), d_MaxwellDrive,
                                                 d_ReynoldsDrive, d_dwdtTotal);
        }

        DeviceAllocator.allocateDeviceArrays(localId, devNums, (devNy + 2 * gridGhost) * gridNxz, d_Maxwell,
                                             d_Reynolds);

        if constexpr (std::is_same_v<::ifDiagDensity, trueType>) {
            if (ifIon)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_IonDensity.size(), d_IonDensity);
            if (ifAlpha)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_AlphaDensity.size(), d_AlphaDensity);
            if (ifBeam)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_BeamDensity.size(), d_BeamDensity);
        }

        if constexpr (std::is_same_v<::ifDiagDiffusivity, trueType>) {
            if (ifIon)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_IonDiffusivity.size(), d_IonDiffusivity);
            if (ifAlpha)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_AlphaDiffusivity.size(), d_AlphaDiffusivity);
            if (ifBeam)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_BeamDiffusivity.size(), d_BeamDiffusivity);
        }

        DeviceAllocator.allocateDeviceArrays(localId, devNums, 1, d_NANFlag);

        /*------------------------------------Output Totals on GPU------------------------------------*/

        const size_t outputTotalSize = h_totalPhi.size();
        if (ifOutputPhi)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, outputTotalSize, d_totalPhi);
        if (ifOutputA)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, outputTotalSize, d_totalA);
        if (ifOutputdNe)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, outputTotalSize, d_totaldNe);
        if (ifOutputdTe)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, outputTotalSize, d_totaldTe);
        if (ifOutputdPi)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, outputTotalSize, d_totaldPi);
        if (ifOutputdPa)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, outputTotalSize, d_totaldPa);
        if (ifOutputdPb)
            DeviceAllocator.allocateDeviceArrays(localId, devNums, outputTotalSize, d_totaldPb);

        if constexpr (std::is_same_v<ifOutputPhaceSpaceDeltaF, trueType>) {
            if (ifIon)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_IonPhaseDeltaF.size(), d_IonPhaseDeltaF);
            if (ifAlpha)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_AlphaPhaseDeltaF.size(), d_AlphaPhaseDeltaF);
            if (ifBeam)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_BeamPhaseDeltaF.size(), d_BeamPhaseDeltaF);
        }

        if constexpr (std::is_same_v<ifOutputPhaceSpacePower, trueType>) {
            if (ifIon)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_IonPhasePower.size(), d_IonPhasePower);
            if (ifAlpha)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_AlphaPhasePower.size(), d_AlphaPhasePower);
            if (ifBeam)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_BeamPhasePower.size(), d_BeamPhasePower);
        }

        if constexpr (std::is_same_v<ifOutputPitchSpaceDeltaF, trueType>) {
            if (ifIon)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_IonPitchDeltaF.size(), d_IonPitchDeltaF);
            if (ifAlpha)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_AlphaPitchDeltaF.size(), d_AlphaPitchDeltaF);
            if (ifBeam)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_BeamPitchDeltaF.size(), d_BeamPitchDeltaF);
        }

        if constexpr (std::is_same_v<ifOutputPitchSpacePower, trueType>) {
            if (ifIon)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_IonPitchPower.size(), d_IonPitchPower);
            if (ifAlpha)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_AlphaPitchPower.size(), d_AlphaPitchPower);
            if (ifBeam)
                DeviceAllocator.allocateDeviceArrays(localId, devNums, h_BeamPitchPower.size(), d_BeamPitchPower);
        }

        /*--------------------------------Select Mode NM Buffer on GPU--------------------------------*/

        DeviceAllocator.allocateDeviceArrays(localId, devNums, gridNy * gridNxz, d_Phi_yxz, d_Phi_xzy, d_A_yxz, d_A_xzy,
                                             d_dNe_yxz, d_dNe_xzy, d_dTe_yxz, d_dTe_xzy, d_dPi_yxz, d_dPi_xzy,
                                             d_dPa_yxz, d_dPa_xzy, d_dPb_yxz, d_dPb_xzy);

        /*------------------------------cuFFT Plan and Frequency Buffer-------------------------------*/

        for (int i = 0; i < devNums; i++) {
            CUDACHECK(cudaSetDevice(localId * devNums + i));
            if constexpr (std::is_same_v<mhdReal, double>) {
                CUFFTCHECK(cufftPlan1d(&d_nPlanR2Cs[i], nFFTTimeSize, CUFFT_D2Z, nFFTBatchSize));
                CUFFTCHECK(cufftPlan1d(&d_nPlanC2Rs[i], nFFTTimeSize, CUFFT_Z2D, nFFTBatchSize));
                CUFFTCHECK(cufftPlan1d(&d_mPlanR2Cs[i], mFFTTimeSize, CUFFT_D2Z, mFFTBatchSize));
                CUFFTCHECK(cufftPlan1d(&d_mPlanC2Rs[i], mFFTTimeSize, CUFFT_Z2D, mFFTBatchSize));
            } else {
                CUFFTCHECK(cufftPlan1d(&d_nPlanR2Cs[i], nFFTTimeSize, CUFFT_R2C, nFFTBatchSize));
                CUFFTCHECK(cufftPlan1d(&d_nPlanC2Rs[i], nFFTTimeSize, CUFFT_C2R, nFFTBatchSize));
                CUFFTCHECK(cufftPlan1d(&d_mPlanR2Cs[i], mFFTTimeSize, CUFFT_R2C, mFFTBatchSize));
                CUFFTCHECK(cufftPlan1d(&d_mPlanC2Rs[i], mFFTTimeSize, CUFFT_C2R, mFFTBatchSize));
            }
            CUFFTCHECK(cufftSetStream(d_nPlanR2Cs[i], 0));
            CUFFTCHECK(cufftSetStream(d_nPlanC2Rs[i], 0));
            CUFFTCHECK(cufftSetStream(d_mPlanR2Cs[i], 0));
            CUFFTCHECK(cufftSetStream(d_mPlanC2Rs[i], 0));
        }

        if constexpr (std::is_same_v<mhdReal, double>) {
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)nFFTBatchSize * nFFTFreqSize, d_nFreqd);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)mFFTBatchSize * mFFTFreqSize, d_mFreqd);
        } else {
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)nFFTBatchSize * nFFTFreqSize, d_nFreqf);
            DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)mFFTBatchSize * mFFTFreqSize, d_mFreqf);
        }

        /*------------------------------------CUB Radix Sort Outer------------------------------------*/

        d_Ion_storage = new void*[devNums]();
        d_Alpha_storage = new void*[devNums]();
        d_Beam_storage = new void*[devNums]();

        /*----------------------------------------CUDA Events-----------------------------------------*/

        for (int i = 0; i < devNums; i++) {
            CUDACHECK(cudaSetDevice(localId * devNums + i));
            CUDACHECK(cudaEventCreate(&startEvents[i]));
            CUDACHECK(cudaEventCreate(&endEvents[i]));
        }

        if (hostId == 0) {
            size_t avail, total, used;
            CUDACHECK(cudaSetDevice(localId * devNums));
            CUDACHECK(cudaMemGetInfo(&avail, &total));
            used = total - avail;
            std::cout << BOLDYELLOW << "Device memory used: " << (double)used / 1024 / 1024 / 1024 << " GB." << RESET
                      << std::endl;
        }

        if (hostId == 0) {
            logDone();
        }
    }
    void releaseDeviceMemory() {

        logStart("Release device memory.");

        Allocator DeviceAllocator;

        /*----------------------------------MHD Perturbation on GPU-----------------------------------*/

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_qtheta);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_w, d_A, d_dNe, d_dTe, d_Phi, d_dJpB, d_dPe);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_w_beg, d_w_midl, d_w_midr, d_w_end, d_A_beg, d_A_midl,
                                            d_A_midr, d_A_end, d_dNe_beg, d_dNe_midl, d_dNe_midr, d_dNe_end, d_dTe_beg,
                                            d_dTe_midl, d_dTe_midr, d_dTe_end, d_Phi_midl, d_Phi_midr, d_dJpB_midl,
                                            d_dJpB_midr, d_dPe_midl, d_dPe_midr, d_dPi_midl, d_dPi_midr, d_dPa_midl,
                                            d_dPa_midr, d_dPb_midl, d_dPb_midr, d_Apt_midl, d_Apt_midr);

        /*-------------------------------Coefficient Compression on GPU-------------------------------*/

        /*---------------------------Linear---------------------------*/

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_A_w, d_A_px_w, d_A_py_w, d_A_pz_w, d_dJpB_w,
                                            d_dJpB_px_w, d_dJpB_py_w, d_dJpB_pz_w, d_dP_w, d_dP_px_w, d_dP_py_w,
                                            d_dP_pz_w, d_w_py_w, d_w_pz_w, d_w2Phi);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_F_perp2, d_F_px_perp2, d_F_pz_perp2, d_F_px2_perp2,
                                            d_F_pxz_perp2, d_F_pz2_perp2);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_A_dJpB, d_A_px_dJpB, d_A_pz_dJpB, d_A_px2_dJpB,
                                            d_A_pxz_dJpB, d_A_pz2_dJpB);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Phi_A, d_Phi_px_A, d_Phi_py_A, d_Phi_pz_A, d_dNe_A,
                                            d_dNe_px_A, d_dNe_py_A, d_dNe_pz_A, d_A_A, d_A_px_A, d_A_py_A, d_A_pz_A);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Phi_dNe, d_Phi_px_dNe, d_Phi_py_dNe, d_Phi_pz_dNe,
                                            d_dPe_dNe, d_dPe_px_dNe, d_dPe_py_dNe, d_dPe_pz_dNe, d_dJpB_dNe,
                                            d_dJpB_px_dNe, d_dJpB_py_dNe, d_dJpB_pz_dNe, d_A_dNe, d_A_px_dNe,
                                            d_A_py_dNe, d_A_pz_dNe);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Phi_dTe, d_Phi_px_dTe, d_Phi_py_dTe, d_Phi_pz_dTe,
                                            d_dTe_dTe, d_dTe_px_dTe, d_dTe_py_dTe, d_dTe_pz_dTe, d_dNe_dTe,
                                            d_dNe_px_dTe, d_dNe_py_dTe, d_dNe_pz_dTe, d_Ne0, d_Te0, d_Ne0_px, d_Te0_px,
                                            d_Pe0_px);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_F2perp2);
        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_A2dJpB);
        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Phi2w);
        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_wdPAdJpB2w);
        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_APhidNe2A);
        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPePhiAdJpB2dNe);
        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_PhidTedNe2dTe);

        /*-------------------------Nonlinear--------------------------*/

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_wPhi_w, d_AdJpB_w, d_PhiA_A, d_NeA_A, d_AdJpB_dNe,
                                            d_dNePhi_dNe, d_PhiTe_dTe, d_PhiTeA_dTe);

        /*-----------------------------------Inverse Matrix on GPU------------------------------------*/

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_laplacianCsrR, d_laplacianCsrC, d_laplacianCsrV);

        if (nablaPerp2A.enable)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_resistiveCsrR, d_resistiveCsrC, d_resistiveCsrV);
        if (nablaPerp2Phi.enable)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_PhiCsrR, d_PhiCsrC, d_PhiCsrV);
        if (nablaPerp2dNe.enable)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dNeCsrR, d_dNeCsrC, d_dNeCsrV);
        if (nablaPerp2dTe.enable)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dTeCsrR, d_dTeCsrC, d_dTeCsrV);
        if (nablaPerp2dPi.enable)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPiCsrR, d_dPiCsrC, d_dPiCsrV);
        if (nablaPerp2dPa.enable)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPaCsrR, d_dPaCsrC, d_dPaCsrV);
        if (nablaPerp2dPb.enable)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPbCsrR, d_dPbCsrC, d_dPbCsrV);

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));

            for (int j = 0; j < devNy; j++) {

                CUDSSCHECK(cudssConfigDestroy(laplacianConfigs[i][j]));
                CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], laplacianDatas[i][j]));
                CUDSSCHECK(cudssMatrixDestroy(laplacianAs[i][j]));
                CUDSSCHECK(cudssMatrixDestroy(laplacianXs[i][j]));
                CUDSSCHECK(cudssMatrixDestroy(laplacianBs[i][j]));

                if (nablaPerp2A.enable) {

                    CUDSSCHECK(cudssConfigDestroy(resistiveConfigs[i][j]));
                    CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], resistiveDatas[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(resistiveAs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(resistiveXs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(resistiveBs[i][j]));
                }

                if (nablaPerp2Phi.enable) {

                    CUDSSCHECK(cudssConfigDestroy(PhiConfigs[i][j]));
                    CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], PhiDatas[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(PhiAs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(PhiXs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(PhiBs[i][j]));
                }

                if (nablaPerp2dNe.enable) {

                    CUDSSCHECK(cudssConfigDestroy(dNeConfigs[i][j]));
                    CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dNeDatas[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dNeAs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dNeXs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dNeBs[i][j]));
                }

                if (nablaPerp2dTe.enable) {

                    CUDSSCHECK(cudssConfigDestroy(dTeConfigs[i][j]));
                    CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dTeDatas[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dTeAs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dTeXs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dTeBs[i][j]));
                }

                if (nablaPerp2dPi.enable) {

                    CUDSSCHECK(cudssConfigDestroy(dPiConfigs[i][j]));
                    CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dPiDatas[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPiAs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPiXs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPiBs[i][j]));
                }

                if (nablaPerp2dPa.enable) {

                    CUDSSCHECK(cudssConfigDestroy(dPaConfigs[i][j]));
                    CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dPaDatas[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPaAs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPaXs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPaBs[i][j]));
                }

                if (nablaPerp2dPb.enable) {

                    CUDSSCHECK(cudssConfigDestroy(dPbConfigs[i][j]));
                    CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dPbDatas[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPbAs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPbXs[i][j]));
                    CUDSSCHECK(cudssMatrixDestroy(dPbBs[i][j]));
                }

                CUDSSCHECK(cudssDestroy(cudssHandles[i][j]));
                CUDACHECK(cudaStreamDestroy(cudaStreams[i][j]));
            }
        }

        /*--------------------------------------Particle in Cell--------------------------------------*/

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_pic_Phi_py_A, d_pic_dNe_py_A, d_pic_A_A, d_pic_A_py_A,
                                            d_pic_A_pz_A);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_pic_APhidNe2A);

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_pic_PhiA_A, d_pic_NeA_A);

        if (ifIon) {

            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Ion_offsets);
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Ion_keys_in, d_Ion_keys_out);
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Ion_values_in, d_Ion_values_out);
        }

        if (ifAlpha) {

            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Alpha_offsets);
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Alpha_keys_in, d_Alpha_keys_out);
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Alpha_values_in, d_Alpha_values_out);
        }

        if (ifBeam) {

            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Beam_offsets);
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Beam_keys_in, d_Beam_keys_out);
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Beam_values_in, d_Beam_values_out);
        }

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_pic1d, d_pic2d, d_pic3d);
        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_globalA, d_globalPhi, d_globalApt, d_globalPi,
                                            d_globalPa, d_globalPb);

        /*-------------------------------------Diagnostic on GPU--------------------------------------*/

        if constexpr (std::is_same_v<::ifDiagAmplitude, trueType>)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_amplitude, d_modeReal, d_modeImag);

        if constexpr (std::is_same_v<::ifDiagFrequency, trueType>)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_frequency);

        if constexpr (std::is_same_v<::ifDiagEparallel, trueType>)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Epara, d_EparaES);

        if constexpr (std::is_same_v<::ifDiagZFDrive, trueType>) {
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_MaxwellDrive, d_ReynoldsDrive, d_dwdtTotal);
        }

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Maxwell, d_Reynolds);

        if constexpr (std::is_same_v<::ifDiagDensity, trueType>) {
            if (ifIon)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonDensity);
            if (ifAlpha)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaDensity);
            if (ifBeam)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamDensity);
        }

        if constexpr (std::is_same_v<::ifDiagDiffusivity, trueType>) {
            if (ifIon)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonDiffusivity);
            if (ifAlpha)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaDiffusivity);
            if (ifBeam)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamDiffusivity);
        }

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_NANFlag);

        /*------------------------------------Output Totals on GPU------------------------------------*/

        if (ifOutputPhi)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_totalPhi);
        if (ifOutputA)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_totalA);
        if (ifOutputdNe)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_totaldNe);
        if (ifOutputdTe)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_totaldTe);
        if (ifOutputdPi)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_totaldPi);
        if (ifOutputdPa)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_totaldPa);
        if (ifOutputdPb)
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_totaldPb);

        if constexpr (std::is_same_v<ifOutputPhaceSpaceDeltaF, trueType>) {
            if (ifIon)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonPhaseDeltaF);
            if (ifAlpha)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaPhaseDeltaF);
            if (ifBeam)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamPhaseDeltaF);
        }

        if constexpr (std::is_same_v<ifOutputPhaceSpacePower, trueType>) {
            if (ifIon)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonPhasePower);
            if (ifAlpha)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaPhasePower);
            if (ifBeam)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamPhasePower);
        }

        if constexpr (std::is_same_v<ifOutputPitchSpaceDeltaF, trueType>) {
            if (ifIon)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonPitchDeltaF);
            if (ifAlpha)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaPitchDeltaF);
            if (ifBeam)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamPitchDeltaF);
        }

        if constexpr (std::is_same_v<ifOutputPitchSpacePower, trueType>) {
            if (ifIon)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonPitchPower);
            if (ifAlpha)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaPitchPower);
            if (ifBeam)
                DeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamPitchPower);
        }

        /*--------------------------------Select Mode NM Buffer on GPU--------------------------------*/

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Phi_yxz, d_Phi_xzy, d_A_yxz, d_A_xzy, d_dNe_yxz,
                                            d_dNe_xzy, d_dTe_yxz, d_dTe_xzy, d_dPi_yxz, d_dPi_xzy, d_dPa_yxz, d_dPa_xzy,
                                            d_dPb_yxz, d_dPb_xzy);

        /*------------------------------cuFFT Plan and Frequency Buffer-------------------------------*/

        if constexpr (std::is_same_v<mhdReal, double>) {
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_nFreqd, d_mFreqd);
        } else {
            DeviceAllocator.releaseDeviceArrays(localId, devNums, d_nFreqf, d_mFreqf);
        }

        for (int i = 0; i < devNums; i++) {
            CUDACHECK(cudaSetDevice(localId * devNums + i));
            CUFFTCHECK(cufftDestroy(d_nPlanR2Cs[i]));
            CUFFTCHECK(cufftDestroy(d_nPlanC2Rs[i]));
            CUFFTCHECK(cufftDestroy(d_mPlanR2Cs[i]));
            CUFFTCHECK(cufftDestroy(d_mPlanC2Rs[i]));
        }

        /*---------------------------------------CUB Radix Sort---------------------------------------*/

        if (ifIon) {
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(cudaFree(d_Ion_storage[i]));
            }
        }
        if (ifAlpha) {
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(cudaFree(d_Alpha_storage[i]));
            }
        }
        if (ifBeam) {
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(cudaFree(d_Beam_storage[i]));
            }
        }
        delete[] d_Ion_storage;
        delete[] d_Alpha_storage;
        delete[] d_Beam_storage;

        /*----------------------------------------CUDA Events-----------------------------------------*/

        for (int i = 0; i < devNums; i++) {
            CUDACHECK(cudaSetDevice(localId * devNums + i));
            CUDACHECK(cudaEventDestroy(startEvents[i]));
            CUDACHECK(cudaEventDestroy(endEvents[i]));
        }

        if (hostId == 0) {
            logDone();
        }
    }
    void memcpyHostToDevice() {

        logStart("Memory copy from host to device.");

        Allocator H2DAllocator;

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));

            H2DAllocator.hostToDevice((devNy + 2 * gridGhost) * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      d_qtheta[i], h_qtheta[0]);

            /*-------------------------------------------Linear-------------------------------------------*/

            /*-------------------------Vorticity--------------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      // Perturbed Parallel Vector Potential in Vorticity
                                      d_A_w[i], h_A_w[0], d_A_px_w[i], h_A_px_w[0], d_A_py_w[i], h_A_py_w[0],
                                      d_A_pz_w[i], h_A_pz_w[0],
                                      // Perturbed Parallel Current in Vorticity
                                      d_dJpB_w[i], h_dJpB_w[0], d_dJpB_px_w[i], h_dJpB_px_w[0], d_dJpB_py_w[i],
                                      h_dJpB_py_w[0], d_dJpB_pz_w[i], h_dJpB_pz_w[0],
                                      // Perturbed Pressure in Vorticity
                                      d_dP_w[i], h_dP_w[0], d_dP_px_w[i], h_dP_px_w[0], d_dP_py_w[i], h_dP_py_w[0],
                                      d_dP_pz_w[i], h_dP_pz_w[0],
                                      // Ion Diamagnetic Drift and Finite Larmor Radius in Vorticity
                                      d_w_py_w[i], h_w_py_w[0], d_w_pz_w[i], h_w_pz_w[0], d_w2Phi[i], h_w_Phi[0]);

            /*-----------------Perturbed Parallel Current-----------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      // Perturbed Parallel Vector Potential in Perturbed Parallel Current
                                      d_A_dJpB[i], h_A_dJpB[0], d_A_px_dJpB[i], h_A_px_dJpB[0], d_A_pz_dJpB[i],
                                      h_A_pz_dJpB[0], d_A_px2_dJpB[i], h_A_px2_dJpB[0], d_A_pxz_dJpB[i],
                                      h_A_pxz_dJpB[0], d_A_pz2_dJpB[i], h_A_pz2_dJpB[0]);

            /*------------Perturbed Parallel Vector Potential-------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      // Perturbed Electric Potential in Perturbed Parallel Vector Potential
                                      d_Phi_A[i], h_Phi_A[0], d_Phi_px_A[i], h_Phi_px_A[0], d_Phi_py_A[i],
                                      h_Phi_py_A[0], d_Phi_pz_A[i], h_Phi_pz_A[0],
                                      // Perturbed Density in Perturbed Parallel Vector Potential
                                      d_dNe_A[i], h_dNe_A[0], d_dNe_px_A[i], h_dNe_px_A[0], d_dNe_py_A[i],
                                      h_dNe_py_A[0], d_dNe_pz_A[i], h_dNe_pz_A[0],
                                      // Perturbed Parallel Vector Potential in Perturbed Parallel Vector Potential
                                      d_A_A[i], h_A_A[0], d_A_px_A[i], h_A_px_A[0], d_A_py_A[i], h_A_py_A[0],
                                      d_A_pz_A[i], h_A_pz_A[0]);

            /*---------------------Perturbed Density----------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      // Perturbed Electric Potential in Perturbed Density
                                      d_Phi_dNe[i], h_Phi_dNe[0], d_Phi_px_dNe[i], h_Phi_px_dNe[0], d_Phi_py_dNe[i],
                                      h_Phi_py_dNe[0], d_Phi_pz_dNe[i], h_Phi_pz_dNe[0],
                                      // Perturbed Pressure in Perturbed Density
                                      d_dPe_dNe[i], h_dPe_dNe[0], d_dPe_px_dNe[i], h_dPe_px_dNe[0], d_dPe_py_dNe[i],
                                      h_dPe_py_dNe[0], d_dPe_pz_dNe[i], h_dPe_pz_dNe[0],
                                      // Perturbed Parallel Current in Perturbed Density
                                      d_dJpB_dNe[i], h_dJpB_dNe[0], d_dJpB_px_dNe[i], h_dJpB_px_dNe[0],
                                      d_dJpB_py_dNe[i], h_dJpB_py_dNe[0], d_dJpB_pz_dNe[i], h_dJpB_pz_dNe[0],
                                      // Perturbed Parallel Vector Potential in Perturbed Density
                                      d_A_dNe[i], h_A_dNe[0], d_A_px_dNe[i], h_A_px_dNe[0], d_A_py_dNe[i],
                                      h_A_py_dNe[0], d_A_pz_dNe[i], h_A_pz_dNe[0]);

            /*-------------------Perturbed Temperature--------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      // Perturbed Electric Potential in Perturbed Temperature
                                      d_Phi_dTe[i], h_Phi_dTe[0], d_Phi_px_dTe[i], h_Phi_px_dTe[0], d_Phi_py_dTe[i],
                                      h_Phi_py_dTe[0], d_Phi_pz_dTe[i], h_Phi_pz_dTe[0],
                                      // Perturbed Temperature in Perturbed Temperature
                                      d_dTe_dTe[i], h_dTe_dTe[0], d_dTe_px_dTe[i], h_dTe_px_dTe[0], d_dTe_py_dTe[i],
                                      h_dTe_py_dTe[0], d_dTe_pz_dTe[i], h_dTe_pz_dTe[0],
                                      // Perturbed Density in Perturbed Temperature
                                      d_dNe_dTe[i], h_dNe_dTe[0], d_dNe_px_dTe[i], h_dNe_px_dTe[0], d_dNe_py_dTe[i],
                                      h_dNe_py_dTe[0], d_dNe_pz_dTe[i], h_dNe_pz_dTe[0]);

            /*------------Equilibrium Density and Temperature-------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx, d_Ne0[i], h_Ne0[0],
                                      d_Te0[i], h_Te0[0], d_Ne0_px[i], h_Ne0_px[0], d_Te0_px[i], h_Te0_px[0],
                                      d_Pe0_px[i], h_Pe0_px[0]);

            /*--------------------------wAdNedTe--------------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_F2perp2[i],
                                      h_F2perp2[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_A2dJpB[i],
                                      h_A2dJpB[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_Phi2w[i],
                                      h_Phi2w[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 10, 0, (hostId * hostNy + i * devNy) * gridNx * 10,
                                      d_wdPAdJpB2w[i], h_wdPAdJpB2w[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_APhidNe2A[i],
                                      h_APhidNe2A[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 11, 0, (hostId * hostNy + i * devNy) * gridNx * 11,
                                      d_dPePhiAdJpB2dNe[i], h_dPePhiAdJpB2dNe[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6,
                                      d_PhidTedNe2dTe[i], h_PhidTedNe2dTe[0][0]);

            /*-----------------------------------------Nonlinear------------------------------------------*/

            /*-------------------------Vorticity--------------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_wPhi_w[i],
                                      h_wPhi_w[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_AdJpB_w[i],
                                      h_AdJpB_w[0][0]);

            /*------------Perturbed Parallel Vector Potential-------------*/

            H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_PhiA_A[i],
                                      h_PhiA_A[0][0], d_NeA_A[i], h_NeA_A[0][0]);

            /*---------------------Perturbed Density----------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_AdJpB_dNe[i],
                                      h_AdJpB_dNe[0][0], d_dNePhi_dNe[i], h_dNePhi_dNe[0][0]);

            /*-------------------Perturbed Temperature--------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_PhiTe_dTe[i],
                                      h_PhiTe_dTe[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 18, 0, (hostId * hostNy + i * devNy) * gridNx * 18,
                                      d_PhiTeA_dTe[i], h_PhiTeA_dTe[0][0]);

            /*--------------------------------------------MHD---------------------------------------------*/

            H2DAllocator.hostToDevice(
                (devNy + 2 * gridGhost) * gridNxz, 0, (hostId * hostNy + i * devNy) * gridNxz, d_Phi_midl[i],
                h_Phi[0][0], d_Phi_midr[i], h_Phi[0][0], d_dJpB_midl[i], h_dJpB[0][0], d_dJpB_midr[i], h_dJpB[0][0],
                d_dPe_midl[i], h_dPe[0][0], d_dPe_midr[i], h_dPe[0][0], d_dPi_midl[i], h_globalPi[0][0], d_dPi_midr[i],
                h_globalPi[0][0], d_dPa_midl[i], h_globalPa[0][0], d_dPa_midr[i], h_globalPa[0][0], d_dPb_midl[i],
                h_globalPb[0][0], d_dPb_midr[i], h_globalPb[0][0], d_Apt_midl[i], h_dPe[0][0], d_Apt_midr[i],
                h_dPe[0][0], d_w_beg[i], h_w[0][0], d_w_midl[i], h_w[0][0], d_w_midr[i], h_dPe[0][0], d_w_end[i],
                h_dPe[0][0], d_A_beg[i], h_A[0][0], d_A_midl[i], h_A[0][0], d_A_midr[i], h_dPe[0][0], d_A_end[i],
                h_dPe[0][0], d_dNe_beg[i], h_dNe[0][0], d_dNe_midl[i], h_dNe[0][0], d_dNe_midr[i], h_dPe[0][0],
                d_dNe_end[i], h_dPe[0][0], d_dTe_beg[i], h_dTe[0][0], d_dTe_midl[i], h_dTe[0][0], d_dTe_midr[i],
                h_dPe[0][0], d_dTe_end[i], h_dPe[0][0]);

            /*--------------------------------------------PIC---------------------------------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx, d_pic_Phi_py_A[i],
                                      h_Phi_py_A[0], d_pic_dNe_py_A[i], h_dNe_py_A[0], d_pic_A_A[i], h_A_A[0],
                                      d_pic_A_py_A[i], h_A_py_A[0], d_pic_A_pz_A[i], h_A_pz_A[0]);

            H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5,
                                      d_pic_APhidNe2A[i], h_APhidNe2A[0][0]);

            H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9,
                                      d_pic_PhiA_A[i], h_PhiA_A[0][0], d_pic_NeA_A[i], h_NeA_A[0][0]);

            if (ifIon) {

                H2DAllocator.hostToDevice(8, 0, i * 8, d_Ion_offsets[i], h_Ion_offsets[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_keys_in[i],
                                          h_Ion_keys[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_keys_out[i],
                                          h_Ion_keys[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_values_in[i],
                                          h_Ion_values[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_values_out[i],
                                          h_Ion_values[0]);
            }

            if (ifAlpha) {

                H2DAllocator.hostToDevice(8, 0, i * 8, d_Alpha_offsets[i], h_Alpha_offsets[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_keys_in[i],
                                          h_Alpha_keys[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_keys_out[i],
                                          h_Alpha_keys[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_values_in[i],
                                          h_Alpha_values[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_values_out[i],
                                          h_Alpha_values[0]);
            }

            if (ifBeam) {

                H2DAllocator.hostToDevice(8, 0, i * 8, d_Beam_offsets[i], h_Beam_offsets[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_keys_in[i],
                                          h_Beam_keys[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_keys_out[i],
                                          h_Beam_keys[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_values_in[i],
                                          h_Beam_values[0]);
                H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_values_out[i],
                                          h_Beam_values[0]);
            }

            H2DAllocator.hostToDevice(cellNx * 30, 0, 0, d_pic1d[i], h_pic1d[0]);
            H2DAllocator.hostToDevice(cellNy * cellNx * 72, 0, 0, d_pic2d[i], h_pic2d[0]);
            H2DAllocator.hostToDevice(gridNyPlusGhost * gridNx * gridNzPlusGhost * 8, 0, 0, d_pic3d[i], h_pic3d[0]);

            H2DAllocator.hostToDevice(gridNyPlusGhost * gridNxz, 0, 0, d_globalPi[i], h_dPe[0][0], d_globalPa[i],
                                      h_dPe[0][0], d_globalPb[i], h_dPe[0][0]);

            H2DAllocator.hostToDevice(gridNyPlusGhost * gridNxz, 0, 0, d_globalA[i], h_dPe[0][0], d_globalPhi[i],
                                      h_dPe[0][0], d_globalApt[i], h_dPe[0][0]);
        }

        if (hostId == 0) {
            size_t avail, total, used;
            CUDACHECK(cudaSetDevice(localId * devNums));
            CUDACHECK(cudaMemGetInfo(&avail, &total));
            used = total - avail;
            std::cout << BOLDYELLOW << "Device memory used: " << (double)used / 1024 / 1024 / 1024 << " GB." << RESET
                      << std::endl;
            logDone();
        }

        /*-----------------------------------CUB Radix Sort Storage-----------------------------------*/

        logStart("Allocate memory on device for sorting particles.");

        auto setupRadixSort = [&](bool enable, void** storage, std::vector<size_t>& storageBytes, int** keys_in,
                                  int** keys_out, picReal** values_in, picReal** values_out, int** offsets) {
            if (!enable)
                return;
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                cub::DeviceSegmentedRadixSort::SortPairs(storage[i], storageBytes[i], keys_out[i], keys_in[i],
                                                         values_out[i], values_in[i], picDev * 7, 7, offsets[i],
                                                         offsets[i] + 1);
                CUDACHECK(cudaMalloc(&storage[i], storageBytes[i]));
                if (!toBool<ifContinue> || continueSteps == 0)
                    cub::DeviceSegmentedRadixSort::SortPairs(storage[i], storageBytes[i], keys_out[i], keys_in[i],
                                                             values_out[i], values_in[i], picDev * 7, 7, offsets[i],
                                                             offsets[i] + 1);
            }
        };
        setupRadixSort(ifIon, d_Ion_storage, d_Ion_storage_bytes, d_Ion_keys_in, d_Ion_keys_out, d_Ion_values_in,
                       d_Ion_values_out, d_Ion_offsets);
        setupRadixSort(ifAlpha, d_Alpha_storage, d_Alpha_storage_bytes, d_Alpha_keys_in, d_Alpha_keys_out,
                       d_Alpha_values_in, d_Alpha_values_out, d_Alpha_offsets);
        setupRadixSort(ifBeam, d_Beam_storage, d_Beam_storage_bytes, d_Beam_keys_in, d_Beam_keys_out, d_Beam_values_in,
                       d_Beam_values_out, d_Beam_offsets);

        if (hostId == 0) {
            size_t avail, total, used;
            CUDACHECK(cudaSetDevice(localId * devNums));
            CUDACHECK(cudaMemGetInfo(&avail, &total));
            used = total - avail;
            std::cout << BOLDYELLOW << "Device memory used: " << (double)used / 1024 / 1024 / 1024 << " GB." << RESET
                      << std::endl;
            logDone();
        }
    }
    void memcpyDeviceToHost(const std::string& finalDir) {

        if (hostId == 0)
            std::cout << BOLDYELLOW << "Start: Output results of this simulation." << RESET << std::endl;

        const MPI_Datatype mpiMhdType = std::is_same_v<mhdReal, double> ? MPI_DOUBLE : MPI_FLOAT;
        const size_t hostSlabLen = (size_t)hostNy * gridNxz;
        const size_t devSlabLen = (size_t)devNy * gridNxz;

        /*------------------------------Density / Diffusivity (6 fields)------------------------------*/

        auto reduceField = [&](std::vector<mhdReal>& hostBuf, mhdReal** devBuf) {
            const size_t fieldLen = hostBuf.size();
            std::vector<mhdReal> localSum(fieldLen, 0);
            std::vector<mhdReal> perDev(fieldLen);
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(cudaMemcpy(perDev.data(), devBuf[i], sizeof(mhdReal) * fieldLen, cudaMemcpyDeviceToHost));
                for (size_t k = 0; k < fieldLen; k++)
                    localSum[k] += perDev[k];
            }
            MPICHECK(MPI_Reduce(localSum.data(), hostBuf.data(), fieldLen, mpiMhdType, MPI_SUM, 0, MPI_COMM_WORLD));
        };

        if constexpr (std::is_same_v<::ifDiagDensity, trueType>) {
            if (ifIon)
                reduceField(h_IonDensity, d_IonDensity);
            if (ifAlpha)
                reduceField(h_AlphaDensity, d_AlphaDensity);
            if (ifBeam)
                reduceField(h_BeamDensity, d_BeamDensity);
        }

        if constexpr (std::is_same_v<::ifDiagDiffusivity, trueType>) {
            if (ifIon)
                reduceField(h_IonDiffusivity, d_IonDiffusivity);
            if (ifAlpha)
                reduceField(h_AlphaDiffusivity, d_AlphaDiffusivity);
            if (ifBeam)
                reduceField(h_BeamDiffusivity, d_BeamDiffusivity);
        }

        if constexpr (std::is_same_v<ifOutputPhaceSpaceDeltaF, trueType>) {
            if (ifIon)
                reduceField(h_IonPhaseDeltaF, d_IonPhaseDeltaF);
            if (ifAlpha)
                reduceField(h_AlphaPhaseDeltaF, d_AlphaPhaseDeltaF);
            if (ifBeam)
                reduceField(h_BeamPhaseDeltaF, d_BeamPhaseDeltaF);
        }

        if constexpr (std::is_same_v<ifOutputPhaceSpacePower, trueType>) {
            if (ifIon)
                reduceField(h_IonPhasePower, d_IonPhasePower);
            if (ifAlpha)
                reduceField(h_AlphaPhasePower, d_AlphaPhasePower);
            if (ifBeam)
                reduceField(h_BeamPhasePower, d_BeamPhasePower);
        }

        if constexpr (std::is_same_v<ifOutputPitchSpaceDeltaF, trueType>) {
            if (ifIon)
                reduceField(h_IonPitchDeltaF, d_IonPitchDeltaF);
            if (ifAlpha)
                reduceField(h_AlphaPitchDeltaF, d_AlphaPitchDeltaF);
            if (ifBeam)
                reduceField(h_BeamPitchDeltaF, d_BeamPitchDeltaF);
        }

        if constexpr (std::is_same_v<ifOutputPitchSpacePower, trueType>) {
            if (ifIon)
                reduceField(h_IonPitchPower, d_IonPitchPower);
            if (ifAlpha)
                reduceField(h_AlphaPitchPower, d_AlphaPitchPower);
            if (ifBeam)
                reduceField(h_BeamPitchPower, d_BeamPitchPower);
        }

        /*--------------------------------MHD Perturbation (10 fields)--------------------------------*/

        auto gatherField = [&](mhdReal* hostBase, mhdReal** midlSrc) {
            std::vector<mhdReal> localSlab(hostSlabLen);
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(cudaMemcpy(localSlab.data() + i * devSlabLen, midlSrc[i] + gridGhost * gridNxz,
                                     sizeof(mhdReal) * devSlabLen, cudaMemcpyDeviceToHost));
            }
            MPICHECK(MPI_Gather(localSlab.data(), hostSlabLen, mpiMhdType, hostBase + gridGhost * gridNxz, hostSlabLen,
                                mpiMhdType, 0, MPI_COMM_WORLD));
        };

        gatherField(h_w[0][0], d_w_midl);
        gatherField(h_A[0][0], d_A_midl);
        gatherField(h_dNe[0][0], d_dNe_midl);
        gatherField(h_dTe[0][0], d_dTe_midl);
        gatherField(h_Phi[0][0], d_Phi_midl);
        gatherField(h_dJpB[0][0], d_dJpB_midl);
        gatherField(h_dPe[0][0], d_dPe_midl);
        gatherField(h_globalPi[0][0], d_dPi_midl);
        gatherField(h_globalPa[0][0], d_dPa_midl);
        gatherField(h_globalPb[0][0], d_dPb_midl);

        /*-------------------------------Output Totals (7 conditional)--------------------------------*/

        CUDACHECK(cudaSetDevice(localId * devNums));
        if (ifOutputPhi)
            CUDACHECK(cudaMemcpy(h_totalPhi.data(), d_totalPhi[0], sizeof(mhdReal) * h_totalPhi.size(),
                                 cudaMemcpyDeviceToHost));
        if (ifOutputA)
            CUDACHECK(
                cudaMemcpy(h_totalA.data(), d_totalA[0], sizeof(mhdReal) * h_totalA.size(), cudaMemcpyDeviceToHost));
        if (ifOutputdNe)
            CUDACHECK(cudaMemcpy(h_totaldNe.data(), d_totaldNe[0], sizeof(mhdReal) * h_totaldNe.size(),
                                 cudaMemcpyDeviceToHost));
        if (ifOutputdTe)
            CUDACHECK(cudaMemcpy(h_totaldTe.data(), d_totaldTe[0], sizeof(mhdReal) * h_totaldTe.size(),
                                 cudaMemcpyDeviceToHost));
        if (ifOutputdPi)
            CUDACHECK(cudaMemcpy(h_totaldPi.data(), d_totaldPi[0], sizeof(mhdReal) * h_totaldPi.size(),
                                 cudaMemcpyDeviceToHost));
        if (ifOutputdPa)
            CUDACHECK(cudaMemcpy(h_totaldPa.data(), d_totaldPa[0], sizeof(mhdReal) * h_totaldPa.size(),
                                 cudaMemcpyDeviceToHost));
        if (ifOutputdPb)
            CUDACHECK(cudaMemcpy(h_totaldPb.data(), d_totaldPb[0], sizeof(mhdReal) * h_totaldPb.size(),
                                 cudaMemcpyDeviceToHost));

        /*-------------------Amplitude / Frequency / Mode / Epara (sampling device)-------------------*/

        int localDevIdx;
        if (hostNums == 1)
            localDevIdx = (devNums == 1) ? 0 : devNums / 2;
        else
            localDevIdx = 0;
        CUDACHECK(cudaSetDevice(localId * devNums + localDevIdx));
        if constexpr (std::is_same_v<::ifDiagAmplitude, trueType>) {
            CUDACHECK(cudaMemcpy(h_amplitude.data(), d_amplitude[localDevIdx], sizeof(mhdReal) * h_amplitude.size(),
                                 cudaMemcpyDeviceToHost));
            CUDACHECK(cudaMemcpy(h_modeReal.data(), d_modeReal[localDevIdx], sizeof(mhdReal) * h_modeReal.size(),
                                 cudaMemcpyDeviceToHost));
            CUDACHECK(cudaMemcpy(h_modeImag.data(), d_modeImag[localDevIdx], sizeof(mhdReal) * h_modeImag.size(),
                                 cudaMemcpyDeviceToHost));
        }
        if constexpr (std::is_same_v<::ifDiagFrequency, trueType>)
            CUDACHECK(cudaMemcpy(h_frequency.data(), d_frequency[localDevIdx], sizeof(mhdReal) * h_frequency.size(),
                                 cudaMemcpyDeviceToHost));
        if constexpr (std::is_same_v<::ifDiagEparallel, trueType>) {
            CUDACHECK(cudaMemcpy(h_Epara.data(), d_Epara[localDevIdx], sizeof(mhdReal) * h_Epara.size(),
                                 cudaMemcpyDeviceToHost));
            CUDACHECK(cudaMemcpy(h_EparaES.data(), d_EparaES[localDevIdx], sizeof(mhdReal) * h_EparaES.size(),
                                 cudaMemcpyDeviceToHost));
        }
        if constexpr (std::is_same_v<::ifDiagZFDrive, trueType>) {
            CUDACHECK(cudaMemcpy(h_MaxwellDrive.data(), d_MaxwellDrive[localDevIdx],
                                 sizeof(mhdReal) * h_MaxwellDrive.size(), cudaMemcpyDeviceToHost));
            CUDACHECK(cudaMemcpy(h_ReynoldsDrive.data(), d_ReynoldsDrive[localDevIdx],
                                 sizeof(mhdReal) * h_ReynoldsDrive.size(), cudaMemcpyDeviceToHost));
            CUDACHECK(cudaMemcpy(h_dwdtTotal.data(), d_dwdtTotal[localDevIdx], sizeof(mhdReal) * h_dwdtTotal.size(),
                                 cudaMemcpyDeviceToHost));
        }

        /*------------------------------PIC Keys / Values (per species)-------------------------------*/

        if (ifIon) {
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(
                    cudaMemcpy(h_Ion_keys[i], d_Ion_keys_in[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToHost));
                CUDACHECK(cudaMemcpy(h_Ion_values[i], d_Ion_values_in[i], sizeof(picReal) * picDev * 7,
                                     cudaMemcpyDeviceToHost));
            }
        }
        if (ifAlpha) {
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(
                    cudaMemcpy(h_Alpha_keys[i], d_Alpha_keys_in[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToHost));
                CUDACHECK(cudaMemcpy(h_Alpha_values[i], d_Alpha_values_in[i], sizeof(picReal) * picDev * 7,
                                     cudaMemcpyDeviceToHost));
            }
        }
        if (ifBeam) {
            for (int i = 0; i < devNums; i++) {
                CUDACHECK(cudaSetDevice(localId * devNums + i));
                CUDACHECK(
                    cudaMemcpy(h_Beam_keys[i], d_Beam_keys_in[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToHost));
                CUDACHECK(cudaMemcpy(h_Beam_values[i], d_Beam_values_in[i], sizeof(picReal) * picDev * 7,
                                     cudaMemcpyDeviceToHost));
            }
        }

        /*-----------------------------------PIC Continuation File------------------------------------*/

        const size_t speciesSize = sizeof(picReal) + sizeof(int) * devNums * 8 + sizeof(int) * devNums * picDev * 7 +
                                   sizeof(picReal) * devNums * picDev * 7;
        const int enabledCount = (ifIon ? 1 : 0) + (ifAlpha ? 1 : 0) + (ifBeam ? 1 : 0);
        const size_t rankSize = enabledCount * speciesSize;

        auto writeFinalPicData = [&](picReal* picConst, int** picOffsets, int** picKeys, picReal** picValues,
                                     size_t intraRankOffset) {
            const std::string filename =
                finalDir + "/PICContinue_" + std::to_string(continueSteps + totalSteps) + ".bin";
            const MPI_Datatype picMPIType = std::is_same_v<picReal, double> ? MPI_DOUBLE : MPI_FLOAT;
            const size_t offset = hostId * rankSize + intraRankOffset;

            MPI_File fileHandle;
            MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                          &fileHandle);

            MPI_File_write_at_all(fileHandle, offset, picConst, 1, picMPIType, MPI_STATUS_IGNORE);
            MPI_File_write_at_all(fileHandle, offset + sizeof(picReal), picOffsets[0], devNums * 8, MPI_INT,
                                  MPI_STATUS_IGNORE);
            MPI_File_write_at_all(fileHandle, offset + sizeof(picReal) + sizeof(int) * devNums * 8, picKeys[0],
                                  devNums * picDev * 7, MPI_INT, MPI_STATUS_IGNORE);
            MPI_File_write_at_all(
                fileHandle, offset + sizeof(picReal) + sizeof(int) * devNums * 8 + sizeof(int) * devNums * picDev * 7,
                picValues[0], devNums * picDev * 7, picMPIType, MPI_STATUS_IGNORE);

            MPI_File_close(&fileHandle);
        };

        size_t off = 0;
        if (ifIon) {
            writeFinalPicData(&IonConst, h_Ion_offsets, h_Ion_keys, h_Ion_values, off);
            off += speciesSize;
        }
        if (ifAlpha) {
            writeFinalPicData(&AlphaConst, h_Alpha_offsets, h_Alpha_keys, h_Alpha_values, off);
            off += speciesSize;
        }
        if (ifBeam)
            writeFinalPicData(&BeamConst, h_Beam_offsets, h_Beam_keys, h_Beam_values, off);

        /*--------------------------------------Bin File Output---------------------------------------*/

        auto writeBin = [&](const std::string& path, const void* data, size_t bytes) {
            std::ofstream o(path, std::ios::out | std::ios::binary);
            o.write((const char*)data, bytes);
        };

        if (hostId == hostNums / 2) {
            if (ifDiagAmplitude) {
                writeBin(finalDir + "/amplitude.bin", h_amplitude.data(), sizeof(mhdReal) * h_amplitude.size());
                writeBin(finalDir + "/RealMode.bin", h_modeReal.data(), sizeof(mhdReal) * h_modeReal.size());
                writeBin(finalDir + "/ImagMode.bin", h_modeImag.data(), sizeof(mhdReal) * h_modeImag.size());
            }
            if (ifDiagFrequency)
                writeBin(finalDir + "/frequency.bin", h_frequency.data(), sizeof(mhdReal) * h_frequency.size());
            if (ifDiagEparallel) {
                writeBin(finalDir + "/Epara.bin", h_Epara.data(), sizeof(mhdReal) * h_Epara.size());
                writeBin(finalDir + "/EparaES.bin", h_EparaES.data(), sizeof(mhdReal) * h_EparaES.size());
            }
            if (ifDiagZFDrive) {
                writeBin(finalDir + "/MaxwellDrive.bin", h_MaxwellDrive.data(),
                         sizeof(mhdReal) * h_MaxwellDrive.size());
                writeBin(finalDir + "/ReynoldsDrive.bin", h_ReynoldsDrive.data(),
                         sizeof(mhdReal) * h_ReynoldsDrive.size());
                writeBin(finalDir + "/ZonalDrive.bin", h_dwdtTotal.data(), sizeof(mhdReal) * h_dwdtTotal.size());
            }
        }

        if (hostId == 0) {
            if (ifDiagDensity) {
                if (ifIon)
                    writeBin(finalDir + "/IonDensity.bin", h_IonDensity.data(), sizeof(mhdReal) * h_IonDensity.size());
                if (ifAlpha)
                    writeBin(finalDir + "/AlphaDensity.bin", h_AlphaDensity.data(),
                             sizeof(mhdReal) * h_AlphaDensity.size());
                if (ifBeam)
                    writeBin(finalDir + "/BeamDensity.bin", h_BeamDensity.data(),
                             sizeof(mhdReal) * h_BeamDensity.size());
            }
            if (ifDiagDiffusivity) {
                if (ifIon)
                    writeBin(finalDir + "/IonDiffusivity.bin", h_IonDiffusivity.data(),
                             sizeof(mhdReal) * h_IonDiffusivity.size());
                if (ifAlpha)
                    writeBin(finalDir + "/AlphaDiffusivity.bin", h_AlphaDiffusivity.data(),
                             sizeof(mhdReal) * h_AlphaDiffusivity.size());
                if (ifBeam)
                    writeBin(finalDir + "/BeamDiffusivity.bin", h_BeamDiffusivity.data(),
                             sizeof(mhdReal) * h_BeamDiffusivity.size());
            }

            if (ifOutputPhi)
                writeBin(finalDir + "/totalPhi.bin", h_totalPhi.data(), sizeof(mhdReal) * h_totalPhi.size());
            if (ifOutputA)
                writeBin(finalDir + "/totalA.bin", h_totalA.data(), sizeof(mhdReal) * h_totalA.size());
            if (ifOutputdNe)
                writeBin(finalDir + "/totaldNe.bin", h_totaldNe.data(), sizeof(mhdReal) * h_totaldNe.size());
            if (ifOutputdTe)
                writeBin(finalDir + "/totaldTe.bin", h_totaldTe.data(), sizeof(mhdReal) * h_totaldTe.size());
            if (ifOutputdPi)
                writeBin(finalDir + "/totaldPi.bin", h_totaldPi.data(), sizeof(mhdReal) * h_totaldPi.size());
            if (ifOutputdPa)
                writeBin(finalDir + "/totaldPa.bin", h_totaldPa.data(), sizeof(mhdReal) * h_totaldPa.size());
            if (ifOutputdPb)
                writeBin(finalDir + "/totaldPb.bin", h_totaldPb.data(), sizeof(mhdReal) * h_totaldPb.size());

            if constexpr (std::is_same_v<ifOutputPhaceSpaceDeltaF, trueType>) {

                const size_t snapShots = (size_t)(totalSteps / outputSteps + 1);

                auto updateBoundary = [&](std::vector<mhdReal>& buf) {
                    for (size_t s = 0; s < snapShots; s++) {
                        mhdReal* slice = buf.data() + s * (size_t)gridE * gridPphi * gridLambda;
                        for (int ie = 0; ie < gridE; ie++) {
                            for (int ip = 0; ip < gridPphi; ip++) {
                                for (int il = 0; il < gridLambda; il++) {
                                    int onE = (ie == 0 || ie == gridE - 1);
                                    int onP = (ip == 0 || ip == gridPphi - 1);
                                    int onL = (il == 0 || il == gridLambda - 1);
                                    int hits = onE + onP + onL;
                                    if (hits == 0)
                                        continue;
                                    mhdReal scale = (hits == 3) ? 8 : (hits == 2) ? 4 : 2;
                                    slice[(ie * gridPphi + ip) * gridLambda + il] *= scale;
                                }
                            }
                        }
                    }
                };

                if (ifIon)
                    updateBoundary(h_IonPhaseDeltaF);
                if (ifAlpha)
                    updateBoundary(h_AlphaPhaseDeltaF);
                if (ifBeam)
                    updateBoundary(h_BeamPhaseDeltaF);

                if (ifIon)
                    writeBin(finalDir + "/IonPhaseDeltaF.bin", h_IonPhaseDeltaF.data(),
                             sizeof(mhdReal) * h_IonPhaseDeltaF.size());
                if (ifAlpha)
                    writeBin(finalDir + "/AlphaPhaseDeltaF.bin", h_AlphaPhaseDeltaF.data(),
                             sizeof(mhdReal) * h_AlphaPhaseDeltaF.size());
                if (ifBeam)
                    writeBin(finalDir + "/BeamPhaseDeltaF.bin", h_BeamPhaseDeltaF.data(),
                             sizeof(mhdReal) * h_BeamPhaseDeltaF.size());
            }

            if constexpr (std::is_same_v<ifOutputPhaceSpacePower, trueType>) {

                auto updateBoundary = [&](std::vector<mhdReal>& buf) {
                    for (int snapShot = 0; snapShot < totalSteps / outputSteps + 1; snapShot++) {
                        for (int mode = leftN; mode <= rightN; mode++) {
                            int modeIdx = mode - leftN;
                            mhdReal* slice = buf.data() + ((size_t)snapShot * (rightN - leftN + 1) + modeIdx) * gridE *
                                                              gridPphi * gridLambda;
                            for (int ie = 0; ie < gridE; ie++) {
                                for (int ip = 0; ip < gridPphi; ip++) {
                                    for (int il = 0; il < gridLambda; il++) {
                                        int onE = (ie == 0 || ie == gridE - 1);
                                        int onP = (ip == 0 || ip == gridPphi - 1);
                                        int onL = (il == 0 || il == gridLambda - 1);
                                        int hits = onE + onP + onL;
                                        if (hits == 0)
                                            continue;
                                        mhdReal scale = (hits == 3) ? 8 : (hits == 2) ? 4 : 2;
                                        slice[(ie * gridPphi + ip) * gridLambda + il] *= scale;
                                    }
                                }
                            }
                        }
                    }
                };

                if (ifIon)
                    updateBoundary(h_IonPhasePower);
                if (ifAlpha)
                    updateBoundary(h_AlphaPhasePower);
                if (ifBeam)
                    updateBoundary(h_BeamPhasePower);

                if (ifIon)
                    writeBin(finalDir + "/IonPhasePower.bin", h_IonPhasePower.data(),
                             sizeof(mhdReal) * h_IonPhasePower.size());
                if (ifAlpha)
                    writeBin(finalDir + "/AlphaPhasePower.bin", h_AlphaPhasePower.data(),
                             sizeof(mhdReal) * h_AlphaPhasePower.size());
                if (ifBeam)
                    writeBin(finalDir + "/BeamPhasePower.bin", h_BeamPhasePower.data(),
                             sizeof(mhdReal) * h_BeamPhasePower.size());
            }

            if constexpr (std::is_same_v<ifOutputPitchSpaceDeltaF, trueType>) {

                const size_t snapShots = (size_t)(totalSteps / outputSteps + 1);

                auto updateBoundary = [&](std::vector<mhdReal>& buf) {
                    for (size_t s = 0; s < snapShots; s++) {
                        mhdReal* slice = buf.data() + s * (size_t)gridVpara * gridVperp;
                        for (int iv = 0; iv < gridVpara; iv++) {
                            for (int ip = 0; ip < gridVperp; ip++) {
                                int onV = (iv == 0 || iv == gridVpara - 1);
                                int onP = (ip == 0 || ip == gridVperp - 1);
                                int hits = onV + onP;
                                if (hits == 0)
                                    continue;
                                mhdReal scale = (hits == 2) ? 4 : 2;
                                slice[iv * gridVperp + ip] *= scale;
                            }
                        }
                    }
                };

                if (ifIon)
                    updateBoundary(h_IonPitchDeltaF);
                if (ifAlpha)
                    updateBoundary(h_AlphaPitchDeltaF);
                if (ifBeam)
                    updateBoundary(h_BeamPitchDeltaF);

                if (ifIon)
                    writeBin(finalDir + "/IonPitchDeltaF.bin", h_IonPitchDeltaF.data(),
                             sizeof(mhdReal) * h_IonPitchDeltaF.size());
                if (ifAlpha)
                    writeBin(finalDir + "/AlphaPitchDeltaF.bin", h_AlphaPitchDeltaF.data(),
                             sizeof(mhdReal) * h_AlphaPitchDeltaF.size());
                if (ifBeam)
                    writeBin(finalDir + "/BeamPitchDeltaF.bin", h_BeamPitchDeltaF.data(),
                             sizeof(mhdReal) * h_BeamPitchDeltaF.size());
            }

            if constexpr (std::is_same_v<ifOutputPitchSpacePower, trueType>) {

                auto updateBoundary = [&](std::vector<mhdReal>& buf) {
                    for (int snapShot = 0; snapShot < totalSteps / outputSteps + 1; snapShot++) {
                        for (int mode = leftN; mode <= rightN; mode++) {
                            int modeIdx = mode - leftN;
                            mhdReal* slice = buf.data() + ((size_t)snapShot * (rightN - leftN + 1) + modeIdx) *
                                                              gridVpara * gridVperp;
                            for (int iv = 0; iv < gridVpara; iv++) {
                                for (int ip = 0; ip < gridVperp; ip++) {
                                    int onV = (iv == 0 || iv == gridVpara - 1);
                                    int onP = (ip == 0 || ip == gridVperp - 1);
                                    int hits = onV + onP;
                                    if (hits == 0)
                                        continue;
                                    mhdReal scale = (hits == 2) ? 4 : 2;
                                    slice[iv * gridVperp + ip] *= scale;
                                }
                            }
                        }
                    }
                };

                if (ifIon)
                    updateBoundary(h_IonPitchPower);
                if (ifAlpha)
                    updateBoundary(h_AlphaPitchPower);
                if (ifBeam)
                    updateBoundary(h_BeamPitchPower);

                if (ifIon)
                    writeBin(finalDir + "/IonPitchPower.bin", h_IonPitchPower.data(),
                             sizeof(mhdReal) * h_IonPitchPower.size());
                if (ifAlpha)
                    writeBin(finalDir + "/AlphaPitchPower.bin", h_AlphaPitchPower.data(),
                             sizeof(mhdReal) * h_AlphaPitchPower.size());
                if (ifBeam)
                    writeBin(finalDir + "/BeamPitchPower.bin", h_BeamPitchPower.data(),
                             sizeof(mhdReal) * h_BeamPitchPower.size());
            }

            writeBin(finalDir + "/w.bin", h_w[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/A.bin", h_A[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/dNe.bin", h_dNe[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/dTe.bin", h_dTe[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/Phi.bin", h_Phi[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/dJpB.bin", h_dJpB[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/dPe.bin", h_dPe[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/dPi.bin", h_globalPi[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/dPa.bin", h_globalPa[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);
            writeBin(finalDir + "/dPb.bin", h_globalPb[0][0] + gridGhost * gridNxz, sizeof(mhdReal) * gridNy * gridNxz);

            std::vector<mhdReal> MHDFinalPerturbation(gridNy * gridNxz * 7);

            for (int j = 0; j < gridNy; j++) {
                for (int i = 0; i < gridNx; i++) {
                    for (int k = 0; k < gridNz; k++) {
                        MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 0 * gridNy * gridNxz] =
                            h_Phi[j + gridGhost][i][k];
                        MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 1 * gridNy * gridNxz] =
                            h_A[j + gridGhost][i][k];
                        MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 2 * gridNy * gridNxz] =
                            h_dNe[j + gridGhost][i][k];
                        MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 3 * gridNy * gridNxz] =
                            h_dTe[j + gridGhost][i][k];
                        MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 4 * gridNy * gridNxz] =
                            h_globalPi[j + gridGhost][i][k];
                        MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 5 * gridNy * gridNxz] =
                            h_globalPa[j + gridGhost][i][k];
                        MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 6 * gridNy * gridNxz] =
                            h_globalPb[j + gridGhost][i][k];
                    }
                }
            }

            writeBin(finalDir + "/MHDContinue_" + std::to_string(continueSteps + totalSteps) + ".bin",
                     MHDFinalPerturbation.data(), sizeof(mhdReal) * MHDFinalPerturbation.size());
        }

        if (hostId == 0) {
            std::cout << BOLDGREEN << "Done." << RESET << std::endl;
            std::cout << std::endl;
        }
    }

    void loadMHDEquilibrium(std::string file) {

        logStart("Load MHD equilibrium.");

        std::ifstream input;
        size_t bytes;
        size_t length;

        input.open(file, std::ios::in | std::ios::binary);
        input.seekg(0, std::ios::end);
        bytes = input.tellg();
        length = bytes / sizeof(double);
        std::vector<double> binaryData(length);
        input.seekg(0, std::ios::beg);
        input.read(reinterpret_cast<char*>(&binaryData[0]), bytes);
        input.close();

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNyPlusGhost; j++) {
                h_qtheta[j][i] = binaryData[i * gridNyPlusGhost + j];
            }
        }

        Allocator B2HAllocator;

        B2HAllocator.binaryToHost(gridNx * gridNyPlusGhost, gridNx, gridNyPlusGhost, binaryData, SFAconxx, SFAconxx_px,
                                  SFAconxx_py, SFAconxy, SFAconxy_px, SFAconxy_py, SFAconyy, SFAconyy_px, SFAconyy_py,
                                  SFAconxz, SFAconxz_px, SFAconxz_py, SFAconyz, SFAconyz_px, SFAconyz_py, SFAconzz,
                                  SFAconzz_px, SFAconzz_py, SFAcovxx, SFAcovxx_px, SFAcovxx_py, SFAcovxy, SFAcovxy_px,
                                  SFAcovxy_py, SFAcovyy, SFAcovyy_px, SFAcovyy_py, SFAcovxz, SFAcovxz_px, SFAcovxz_py,
                                  SFAcovyz, SFAcovyz_px, SFAcovyz_py, SFAcovzz, SFAcovzz_px, SFAcovzz_py);

        B2HAllocator.binaryToHost(
            gridNx * gridNyPlusGhost * 37, gridNx, gridNy, binaryData, q, q_px, psip, psip_px, Ni, Ni_px, Ti, Ti_px, Pi,
            Pi_px, Ne, Ne_px, Te, Te_px, Pe, Pe_px, Na, Na_px, Ta, Ta_px, Nb, Nb_px, Tb, Tb_px, gconxx, gconxx_px,
            gconxx_py, gconxy, gconxy_px, gconxy_py, gconyy, gconyy_px, gconyy_py, gconxz, gconxz_px, gconxz_py, gconyz,
            gconyz_px, gconyz_py, gconzz, gconzz_px, gconzz_py, gcovxx, gcovxx_px, gcovxx_py, gcovxy, gcovxy_px,
            gcovxy_py, gcovyy, gcovyy_px, gcovyy_py, gcovxz, gcovxz_px, gcovxz_py, gcovyz, gcovyz_px, gcovyz_py, gcovzz,
            gcovzz_px, gcovzz_py, J, J_px, J_py, Bny, Bny_px, Bny_py, JpB, JpB_px, JpB_py, Rho, Rho_px, Rho_py, Va,
            Va_px, Va_py, R, Z, B, B_px, B_py, B_px2, B_pxy, B_py2);

        auto loadPICProfile = [&](int index, int offset, double**& field) {
            h_pic1d[index][offset + 0] = field[index][0];
            h_pic1d[index][offset + 1] = field[index + 1][0];
        };

        auto loadPICStraight = [&](int index, int offset, int i, int j, double**& field) {
            j = (j - gridGhost + gridNy) % gridNy;
            h_pic2d[index][offset + 0] = field[i][j];
            h_pic2d[index][offset + 1] = field[i + 1][j];
            h_pic2d[index][offset + 2] = field[i][j + 1];
            h_pic2d[index][offset + 3] = field[i + 1][j + 1];
        };

        auto loadPICAligned = [&](int index, int offset, int i, int j, double**& field) {
            h_pic2d[index][offset + 0] = field[i][j];
            h_pic2d[index][offset + 1] = field[i + 1][j];
            h_pic2d[index][offset + 2] = field[i][j + 1];
            h_pic2d[index][offset + 3] = field[i + 1][j + 1];
        };

        for (int index = 0; index < cellNx; index++) {

            loadPICProfile(index, 0, q);
            loadPICProfile(index, 2, q_px);
            loadPICProfile(index, 4, Na);
            loadPICProfile(index, 6, Na_px);
            loadPICProfile(index, 8, Nb);
            loadPICProfile(index, 10, Nb_px);
            loadPICProfile(index, 12, Ni);
            loadPICProfile(index, 14, Ni_px);
            loadPICProfile(index, 16, Ta);
            loadPICProfile(index, 18, Ta_px);
            loadPICProfile(index, 20, Tb);
            loadPICProfile(index, 22, Tb_px);
            loadPICProfile(index, 24, Ti);
            loadPICProfile(index, 26, Ti_px);
            loadPICProfile(index, 28, psip);

            for (int i = 4; i < 16; i++)
                h_pic1d[index][i] *= 1.0e-19;
        }

        for (int index = 0; index < cellNy * cellNx; index++) {

            int i = index % cellNx;
            int j = index / cellNx;

            loadPICStraight(index, 0, i, j, J);
            loadPICStraight(index, 4, i, j, B);
            loadPICStraight(index, 8, i, j, J_px);
            loadPICStraight(index, 12, i, j, J_py);
            loadPICStraight(index, 16, i, j, B_px);
            loadPICStraight(index, 20, i, j, B_py);

            loadPICAligned(index, 24, i, j, SFAcovxy);
            loadPICAligned(index, 28, i, j, SFAcovyy);
            loadPICAligned(index, 32, i, j, SFAcovyz);
            loadPICAligned(index, 36, i, j, SFAcovxy_py);
            loadPICAligned(index, 40, i, j, SFAcovyy_px);
            loadPICAligned(index, 44, i, j, SFAcovyz_px);
            loadPICAligned(index, 48, i, j, SFAcovyz_py);
            loadPICAligned(index, 52, i, j, SFAconxx);
            loadPICAligned(index, 56, i, j, SFAconxy);
            loadPICAligned(index, 60, i, j, SFAconyy);

            loadPICStraight(index, 64, i, j, R);
            loadPICStraight(index, 68, i, j, Z);
        }

        if (hostId == 0) {
            logDone();
        }
    }
    void loadMHDPerturbation(std::string file) {

        logStart("Load MHD perturbation.");

        std::ifstream input;
        size_t bytes;
        size_t length;

        input.open(file, std::ios::in | std::ios::binary);
        input.seekg(0, std::ios::end);
        bytes = input.tellg();
        length = bytes / sizeof(double);
        std::vector<double> binaryData(length);
        input.seekg(0, std::ios::beg);
        input.read(reinterpret_cast<char*>(&binaryData[0]), bytes);
        input.close();

        for (int j = 0; j < gridNy; j++) {
            for (int i = 0; i < gridNx; i++) {
                for (int k = 0; k < gridNz; k++) {
                    h_Phi[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 0 * gridNy * gridNxz];
                    h_A[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 1 * gridNy * gridNxz];
                    h_dNe[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 2 * gridNy * gridNxz];
                    h_dTe[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 3 * gridNy * gridNxz];
                    h_globalPi[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 4 * gridNy * gridNxz];
                    h_globalPa[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 5 * gridNy * gridNxz];
                    h_globalPb[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 6 * gridNy * gridNxz];
                }
            }
        }

        if (hostId == 0) {
            logDone();
        }
    }
    template <int leftN, int rightN>
    void computeMHDPerturbation(int radialIndex, mhdReal width, mhdReal amplitude, std::string initialDir) {

        logStart("Load MHD perturbation.");

        const mhdReal xCenter = (radialIndex - 1) * gridDx;
        const mhdReal sigma2 = 2.0 * width * width;
        const mhdReal qCenter = q[radialIndex - 1][0];
        const mhdReal qRounded = std::round(qCenter * 10.0) / 10.0;

        const int primaryN = rightN * tubes;
        const int primaryM = static_cast<int>(std::round(primaryN * qRounded));

        for (int j = 0; j < gridNy; j++) {
            const mhdReal y = ((j + 0.5) / gridNy) * 2.0 * PI - PI;
            for (int i = 1; i < gridNx - 1; i++) {
                const mhdReal x = i * gridDx;
                const mhdReal gaussian = std::exp(-(x - xCenter) * (x - xCenter) / sigma2);
                for (int k = 0; k < gridNz; k++) {
                    const mhdReal z = ((k + 0.5) / gridNz) * 2.0 * PI / tubes - PI / tubes;

                    mhdReal disturb = amplitude * gaussian * std::cos(primaryM * y - primaryN * z);

                    if constexpr (leftN != rightN) {
                        for (int n = leftN * tubes; n < rightN * tubes; n += tubes) {
                            const int m = static_cast<int>(std::round(n * qRounded));
                            const mhdReal shiftY = 2.0 * PI / primaryM * m;
                            const mhdReal shiftZ = 2.0 * PI / primaryN * n;
                            disturb += amplitude * gaussian * std::cos(m * (y + shiftY) - n * (z + shiftZ));
                        }
                    }

                    h_Phi[j + gridGhost][i][k] = disturb;
                }
            }
        }

        if (hostId == 0) {
            std::vector<mhdReal> buf(gridNy * gridNxz * 7, 0);
            for (int j = 0; j < gridNy; j++)
                for (int i = 0; i < gridNx; i++)
                    for (int k = 0; k < gridNz; k++)
                        buf[j * gridNxz + i * gridNz + k + 0 * gridNy * gridNxz] = h_Phi[j + gridGhost][i][k];

            std::ofstream output(initialDir + "/MHDContinue_0.bin", std::ios::out | std::ios::binary);
            output.write(reinterpret_cast<char*>(buf.data()), sizeof(mhdReal) * buf.size());
            output.close();

            logDone();
        }
    }
    void compressCollocatedCoefficient() {

        logStart("Compress collocated coefficient in shifted metric coordinate.");

        /*---------------------------Linear---------------------------*/

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                // Vorticity

                h_A_w[j][i] =
                    std::pow(B[i][j], -2.0) * std::pow(J[i][j], -1.0) *
                    (B[i][j] * gcovyz[i][j] * (Bny_py[i][j] * JpB_px[i][j] + (-1.0) * Bny_px[i][j] * JpB_py[i][j]) +
                     Bny[i][j] *
                         (gcovyz[i][j] * ((-1.0) * B_py[i][j] * JpB_px[i][j] + B_px[i][j] * JpB_py[i][j]) +
                          B[i][j] * (gcovyz_py[i][j] * JpB_px[i][j] + (-1.0) * gcovyz_px[i][j] * JpB_py[i][j])));
                h_A_px_w[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0) *
                                 JpB_py[i][j];
                h_A_py_w[j][i] =
                    std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0) * JpB_px[i][j];
                h_A_pz_w[j][i] = std::pow(B[i][j], -1.0) * Bny[i][j] * std::pow(J[i][j], -1.0) *
                                 ((-1.0) * gcovyy[i][j] * JpB_px[i][j] + gcovxy[i][j] * JpB_py[i][j]);

                h_dJpB_w[j][i] = 0.0;
                h_dJpB_px_w[j][i] = 0.0;
                h_dJpB_py_w[j][i] = Bny[i][j];
                h_dJpB_pz_w[j][i] = 0.0;

                h_dP_w[j][i] = 0.0;
                h_dP_px_w[j][i] = std::pow(B[i][j], -3.0) *
                                  (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                                   Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                                  std::pow(J[i][j], -1.0);
                h_dP_py_w[j][i] =
                    std::pow(B[i][j], -4.0) *
                    (B[i][j] * Bny[i][j] * B_px[i][j] * gcovyz[i][j] +
                     (-1.0) * std::pow(B[i][j], 2.0) * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j]) +
                     std::pow(Bny[i][j], 3.0) *
                         ((gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]) * gcovyz[i][j] + gcovyy[i][j] * gcovyz_px[i][j] +
                          (-1.0) * gcovxy[i][j] * gcovyz_py[i][j])) *
                    std::pow(J[i][j], -1.0);
                h_dP_pz_w[j][i] = std::pow(B[i][j], -3.0) *
                                  (B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j]) +
                                   Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j] +
                                                B[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) *
                                  std::pow(J[i][j], -1.0);

                h_w_py_w[j][i] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                 std::pow(J[i][j], -1.0) * std::pow(Ni[i][j], -1.0) * Pi_px[i][j] *
                                 std::pow(NormQE, -1.0);
                h_w_pz_w[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                 std::pow(J[i][j], -1.0) * std::pow(Ni[i][j], -1.0) * Pi_px[i][j] *
                                 std::pow(NormQE, -1.0);
                h_w_Phi[j][i] = -std::pow(Rho[i][j], 2.0) * std::pow(Va[i][j], 2.0);

                h_wdPAdJpB2w[j][i][0] = h_w_py_w[j][i];
                h_wdPAdJpB2w[j][i][1] = h_w_pz_w[j][i];
                h_wdPAdJpB2w[j][i][2] = h_dP_px_w[j][i];
                h_wdPAdJpB2w[j][i][3] = h_dP_py_w[j][i];
                h_wdPAdJpB2w[j][i][4] = h_dP_pz_w[j][i];
                h_wdPAdJpB2w[j][i][5] = h_A_w[j][i];
                h_wdPAdJpB2w[j][i][6] = h_A_px_w[j][i];
                h_wdPAdJpB2w[j][i][7] = h_A_py_w[j][i];
                h_wdPAdJpB2w[j][i][8] = h_A_pz_w[j][i];
                h_wdPAdJpB2w[j][i][9] = h_dJpB_py_w[j][i];

                // Poisson Equation

                h_Phi_w[j][i] = 0.0;
                h_Phi_px_w[j][i] =
                    std::pow(Va[i][j], -3.0) * (std::pow(J[i][j], -1.0) *
                                                    ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] +
                                                     gconxx[i][j] * J_px[i][j] + gconxy[i][j] * J_py[i][j]) *
                                                    Va[i][j] +
                                                (-2.0) * (gconxx[i][j] * Va_px[i][j] + gconxy[i][j] * Va_py[i][j]));
                h_Phi_pz_w[j][i] =
                    std::pow(Va[i][j], -3.0) * (std::pow(J[i][j], -1.0) *
                                                    ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] +
                                                     gconxz[i][j] * J_px[i][j] + gconyz[i][j] * J_py[i][j]) *
                                                    Va[i][j] +
                                                (-2.0) * (gconxz[i][j] * Va_px[i][j] + gconyz[i][j] * Va_py[i][j]));
                h_Phi_px2_w[j][i] = gconxx[i][j] * std::pow(Va[i][j], -2.0);
                h_Phi_pxz_w[j][i] = 2.0 * gconxz[i][j] * std::pow(Va[i][j], -2.0);
                h_Phi_pz2_w[j][i] = gconzz[i][j] * std::pow(Va[i][j], -2.0);
                h_Phi2w[j][i][0] = h_Phi_px_w[j][i];
                h_Phi2w[j][i][1] = h_Phi_pz_w[j][i];
                h_Phi2w[j][i][2] = h_Phi_px2_w[j][i];
                h_Phi2w[j][i][3] = h_Phi_pxz_w[j][i];
                h_Phi2w[j][i][4] = h_Phi_pz2_w[j][i];

                // Perturbed Parallel Current

                h_A_dJpB[j][i] = std::pow(B[i][j], -2.0) * std::pow(J[i][j], -1.0) *
                                 ((B_px2[i][j] * gconxx[i][j] + B_px[i][j] * gconxx_px[i][j] +
                                   2.0 * B_pxy[i][j] * gconxy[i][j] + B_px[i][j] * gconxy_py[i][j] +
                                   B_py2[i][j] * gconyy[i][j] + B_py[i][j] * (gconxy_px[i][j] + gconyy_py[i][j])) *
                                      J[i][j] +
                                  B_px[i][j] * gconxx[i][j] * J_px[i][j] + B_py[i][j] * gconxy[i][j] * J_px[i][j] +
                                  B_px[i][j] * gconxy[i][j] * J_py[i][j] + B_py[i][j] * gconyy[i][j] * J_py[i][j]);
                h_A_px_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * std::pow(J[i][j], -1.0) *
                                    ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] + gconxx[i][j] * J_px[i][j] +
                                     gconxy[i][j] * J_py[i][j]);
                h_A_pz_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * std::pow(J[i][j], -1.0) *
                                    ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] + gconxz[i][j] * J_px[i][j] +
                                     gconyz[i][j] * J_py[i][j]);
                h_A_px2_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * gconxx[i][j];
                h_A_pxz_dJpB[j][i] = (-2.0) * std::pow(B[i][j], -1.0) * gconxz[i][j];
                h_A_pz2_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * gconzz[i][j];
                h_A2dJpB[j][i][0] = h_A_dJpB[j][i];
                h_A2dJpB[j][i][1] = h_A_px_dJpB[j][i];
                h_A2dJpB[j][i][2] = h_A_pz_dJpB[j][i];
                h_A2dJpB[j][i][3] = h_A_px2_dJpB[j][i];
                h_A2dJpB[j][i][4] = h_A_pxz_dJpB[j][i];
                h_A2dJpB[j][i][5] = h_A_pz2_dJpB[j][i];

                // Parallel Resistive

                h_A_resistive[j][i] = h_A_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef + 1.0;
                h_A_px_resistive[j][i] = h_A_px_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_pz_resistive[j][i] = h_A_pz_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_px2_resistive[j][i] = h_A_px2_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_pxz_resistive[j][i] = h_A_pxz_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_pz2_resistive[j][i] = h_A_pz2_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;

                // Perpendicular Dissipation

                h_F_perp2[j][i] = 1.0;
                h_F_px_perp2[j][i] = (1.0) * std::pow(J[i][j], -1.0) *
                                     ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] + gconxx[i][j] * J_px[i][j] +
                                      gconxy[i][j] * J_py[i][j]);
                h_F_pz_perp2[j][i] = (1.0) * std::pow(J[i][j], -1.0) *
                                     ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] + gconxz[i][j] * J_px[i][j] +
                                      gconyz[i][j] * J_py[i][j]);
                h_F_px2_perp2[j][i] = (1.0) * gconxx[i][j];
                h_F_pxz_perp2[j][i] = (2.0) * gconxz[i][j];
                h_F_pz2_perp2[j][i] = (1.0) * gconzz[i][j];
                h_F2perp2[j][i][0] = h_F_px_perp2[j][i];
                h_F2perp2[j][i][1] = h_F_pz_perp2[j][i];
                h_F2perp2[j][i][2] = h_F_px2_perp2[j][i];
                h_F2perp2[j][i][3] = h_F_pxz_perp2[j][i];
                h_F2perp2[j][i][4] = h_F_pz2_perp2[j][i];

                // Perturbed Parallel Vector Potential

                h_Phi_A[j][i] = 0.0;
                h_Phi_px_A[j][i] = 0.0;
                h_Phi_py_A[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j];
                h_Phi_pz_A[j][i] = 0.0;

                h_dNe_A[j][i] = 0.0;
                h_dNe_px_A[j][i] = 0.0;
                h_dNe_py_A[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j] * std::pow(NormQE, -1.0);
                h_dNe_pz_A[j][i] = 0.0;

                h_A_A[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                              (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                               Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                              std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                              std::pow(NormQE, -1.0);
                h_A_px_A[j][i] = 0.0;
                h_A_py_A[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                 std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                                 std::pow(NormQE, -1.0);
                h_A_pz_A[j][i] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                 std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                                 std::pow(NormQE, -1.0);

                h_APhidNe2A[j][i][0] = h_A_A[j][i];
                h_APhidNe2A[j][i][1] = h_A_py_A[j][i];
                h_APhidNe2A[j][i][2] = h_A_pz_A[j][i];
                h_APhidNe2A[j][i][3] = h_Phi_py_A[j][i];
                h_APhidNe2A[j][i][4] = h_dNe_py_A[j][i];

                // Perturbed Density

                h_Phi_dNe[j][i] = 0.0;
                h_Phi_px_dNe[j][i] = (-1.0) * std::pow(B[i][j], -3.0) *
                                     (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                                      Bny[i][j] * ((-2.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                                     std::pow(J[i][j], -1.0) * Ne[i][j];
                h_Phi_py_dNe[j][i] =
                    std::pow(B[i][j], -3.0) * std::pow(J[i][j], -1.0) *
                    (B[i][j] * Bny_px[i][j] * gcovyz[i][j] * Ne[i][j] +
                     Bny[i][j] * (B[i][j] * gcovyz_px[i][j] * Ne[i][j] +
                                  gcovyz[i][j] * ((-2.0) * B_px[i][j] * Ne[i][j] + B[i][j] * Ne_px[i][j])));
                h_Phi_pz_dNe[j][i] =
                    std::pow(B[i][j], -3.0) * std::pow(J[i][j], -1.0) *
                    (B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j]) * Ne[i][j] +
                     (-1.0) * Bny[i][j] *
                         (2.0 * B_py[i][j] * gcovxy[i][j] * Ne[i][j] + (-2.0) * B_px[i][j] * gcovyy[i][j] * Ne[i][j] +
                          B[i][j] *
                              (((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]) * Ne[i][j] + gcovyy[i][j] * Ne_px[i][j])));

                h_dPe_dNe[j][i] = 0.0;
                h_dPe_px_dNe[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                                     (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                                      Bny[i][j] * ((-2.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                                     std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_dPe_py_dNe[j][i] =
                    (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                    ((-1.0) * B[i][j] * Bny_px[i][j] * gcovyz[i][j] +
                     Bny[i][j] * (2.0 * B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_px[i][j])) *
                    std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_dPe_pz_dNe[j][i] =
                    (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                    (B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j]) +
                     Bny[i][j] * (2.0 * B_py[i][j] * gcovxy[i][j] + (-2.0) * B_px[i][j] * gcovyy[i][j] +
                                  B[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) *
                    std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);

                h_dJpB_dNe[j][i] = 0.0;
                h_dJpB_px_dNe[j][i] = 0.0;
                h_dJpB_py_dNe[j][i] = Bny[i][j] * std::pow(NormQE, -1.0);
                h_dJpB_pz_dNe[j][i] = 0.0;

                h_A_dNe[j][i] =
                    std::pow(B[i][j], -2.0) * std::pow(J[i][j], -1.0) *
                    (B[i][j] * gcovyz[i][j] * (Bny_py[i][j] * JpB_px[i][j] + (-1.0) * Bny_px[i][j] * JpB_py[i][j]) +
                     Bny[i][j] *
                         (gcovyz[i][j] * ((-1.0) * B_py[i][j] * JpB_px[i][j] + B_px[i][j] * JpB_py[i][j]) +
                          B[i][j] * (gcovyz_py[i][j] * JpB_px[i][j] + (-1.0) * gcovyz_px[i][j] * JpB_py[i][j]))) *
                    std::pow(NormQE, -1.0);
                h_A_px_dNe[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] *
                                   std::pow(J[i][j], -1.0) * JpB_py[i][j] * std::pow(NormQE, -1.0);
                h_A_py_dNe[j][i] = std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0) *
                                   JpB_px[i][j] * std::pow(NormQE, -1.0);
                h_A_pz_dNe[j][i] = std::pow(B[i][j], -1.0) * Bny[i][j] * std::pow(J[i][j], -1.0) *
                                   ((-1.0) * gcovyy[i][j] * JpB_px[i][j] + gcovxy[i][j] * JpB_py[i][j]) *
                                   std::pow(NormQE, -1.0);

                h_dPePhiAdJpB2dNe[j][i][0] = h_dPe_px_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][1] = h_dPe_py_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][2] = h_dPe_pz_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][3] = h_Phi_px_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][4] = h_Phi_py_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][5] = h_Phi_pz_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][6] = h_A_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][7] = h_A_px_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][8] = h_A_py_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][9] = h_A_pz_dNe[j][i];
                h_dPePhiAdJpB2dNe[j][i][10] = h_dJpB_py_dNe[j][i];

                // Perturbed Temperature

                h_Phi_dTe[j][i] = 0.0;
                h_Phi_px_dTe[j][i] = 0.0;
                h_Phi_py_dTe[j][i] =
                    std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0) * Te_px[i][j];
                h_Phi_pz_dTe[j][i] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0) * Te_px[i][j];

                h_dTe_dTe[j][i] = 0.0;
                h_dTe_px_dTe[j][i] = 0.0;
                h_dTe_py_dTe[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                     std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                                     std::pow(NormQE, -1.0);
                h_dTe_pz_dTe[j][i] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                     std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                                     std::pow(NormQE, -1.0);

                h_dNe_dTe[j][i] = 0.0;
                h_dNe_px_dTe[j][i] = 0.0;
                h_dNe_py_dTe[j][i] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                     std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Te[i][j] * Te_px[i][j] *
                                     std::pow(NormQE, -1.0);
                h_dNe_pz_dTe[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                     std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Te[i][j] * Te_px[i][j] *
                                     std::pow(NormQE, -1.0);

                h_PhidTedNe2dTe[j][i][0] = h_Phi_py_dTe[j][i];
                h_PhidTedNe2dTe[j][i][1] = h_Phi_pz_dTe[j][i];
                h_PhidTedNe2dTe[j][i][2] = h_dTe_py_dTe[j][i];
                h_PhidTedNe2dTe[j][i][3] = h_dTe_pz_dTe[j][i];
                h_PhidTedNe2dTe[j][i][4] = h_dNe_py_dTe[j][i];
                h_PhidTedNe2dTe[j][i][5] = h_dNe_pz_dTe[j][i];

                // Equilibrium Density and Temperature

                h_Ne0[j][i] = Ne[i][j];
                h_Te0[j][i] = Te[i][j];
                h_Ne0_px[j][i] = Ne_px[i][j];
                h_Te0_px[j][i] = Te_px[i][j];
                h_Pe0_px[j][i] = Pe_px[i][j];
            }
        }

        /*-------------------------Nonlinear--------------------------*/

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                // Vorticity

                h_wPhi_w[j][i][0] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_wPhi_w[j][i][1] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_wPhi_w[j][i][2] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_wPhi_w[j][i][3] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);
                h_wPhi_w[j][i][4] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_wPhi_w[j][i][5] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);

                h_AdJpB_w[j][i][0] = std::pow(B[i][j], -2.0) *
                                     (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                                      Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                                     std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][1] = std::pow(B[i][j], -2.0) *
                                     (Bny[i][j] * B_px[i][j] * gcovyz[i][j] +
                                      (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) *
                                     std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][2] = std::pow(B[i][j], -2.0) *
                                     (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) +
                                      B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] +
                                                 Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) *
                                     std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][3] =
                    (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][4] = std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][5] = std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][6] =
                    (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][7] =
                    (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_AdJpB_w[j][i][8] = std::pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);

                // Perturbed Parallel Vector Potential

                h_PhiA_A[j][i][0] = std::pow(B[i][j], -3.0) *
                                    (Bny[i][j] * B_py[i][j] * gcovyz[i][j] +
                                     (-1.0) * B[i][j] * (Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_py[i][j])) *
                                    std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][1] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][2] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][3] = std::pow(B[i][j], -3.0) *
                                    (B[i][j] * Bny_px[i][j] * gcovyz[i][j] +
                                     Bny[i][j] * ((-1.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) *
                                    std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][4] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][5] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][6] = std::pow(B[i][j], -3.0) *
                                    (Bny[i][j] * ((-1.0) * B_py[i][j] * gcovxy[i][j] + B_px[i][j] * gcovyy[i][j]) +
                                     B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j] +
                                                Bny[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) *
                                    std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][7] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][8] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);

                h_NeA_A[j][i][0] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                                   (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                                    Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][1] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][2] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][3] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                                   (Bny[i][j] * B_px[i][j] * gcovyz[i][j] +
                                    (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][4] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][5] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][6] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                                   (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) +
                                    B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] +
                                               Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][7] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][8] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];

                // Perturbed Density

                h_AdJpB_dNe[j][i][0] = std::pow(B[i][j], -2.0) *
                                       (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                                        Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                                       std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][1] =
                    std::pow(B[i][j], -2.0) *
                    (Bny[i][j] * B_px[i][j] * gcovyz[i][j] +
                     (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) *
                    std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][2] = std::pow(B[i][j], -2.0) *
                                       (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) +
                                        B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] +
                                                   Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) *
                                       std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][3] = (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] *
                                       std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][4] = std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0) *
                                       std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][5] = std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0) *
                                       std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][6] = (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] *
                                       std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][7] = (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] *
                                       std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0);
                h_AdJpB_dNe[j][i][8] = std::pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0) *
                                       std::pow(NormQE, -1.0);

                h_dNePhi_dNe[j][i][0] =
                    std::pow(B[i][j], -3.0) *
                    ((-1.0) * B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                     Bny[i][j] * (2.0 * B_py[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_py[i][j])) *
                    std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][1] = std::pow(B[i][j], -3.0) *
                                        (B[i][j] * Bny_px[i][j] * gcovyz[i][j] +
                                         Bny[i][j] * ((-2.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) *
                                        std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][2] =
                    std::pow(B[i][j], -3.0) *
                    (B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j]) +
                     Bny[i][j] * ((-2.0) * B_py[i][j] * gcovxy[i][j] + 2.0 * B_px[i][j] * gcovyy[i][j] +
                                  B[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) *
                    std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][3] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][4] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][5] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][6] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][7] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_dNePhi_dNe[j][i][8] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);

                // Perturbed Temperature

                h_PhiTe_dTe[j][i][0] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_PhiTe_dTe[j][i][1] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiTe_dTe[j][i][2] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_PhiTe_dTe[j][i][3] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiTe_dTe[j][i][4] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiTe_dTe[j][i][5] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);

                h_PhiTeA_dTe[j][i][0] =
                    std::pow(B[i][j], -4.0) *
                    (gcovyz[i][j] * (B[i][j] * Bny_px[i][j] * gcovyz[i][j] +
                                     Bny[i][j] * ((-1.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) +
                     gcovxz[i][j] * (Bny[i][j] * B_py[i][j] * gcovyz[i][j] +
                                     (-1.0) * B[i][j] * (Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_py[i][j])) +
                     (Bny[i][j] * ((-1.0) * B_py[i][j] * gcovxy[i][j] + B_px[i][j] * gcovyy[i][j]) +
                      B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j] +
                                 Bny[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) *
                         gcovzz[i][j]) *
                    std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][1] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                        (std::pow(gcovyz[i][j], 2.0) + (-1.0) * gcovyy[i][j] * gcovzz[i][j]) *
                                        std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][2] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                        ((-1.0) * gcovxz[i][j] * gcovyz[i][j] + gcovxy[i][j] * gcovzz[i][j]) *
                                        std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][3] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                        (gcovxz[i][j] * gcovyy[i][j] + (-1.0) * gcovxy[i][j] * gcovyz[i][j]) *
                                        std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][4] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                        (((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]) * gcovyz[i][j] +
                                         (-1.0) * gcovyy[i][j] * gcovyz_px[i][j] + gcovxy[i][j] * gcovyz_py[i][j]) *
                                        std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][5] =
                    std::pow(B[i][j], -4.0) *
                    (Bny[i][j] *
                         (B_px[i][j] * std::pow(gcovyz[i][j], 2.0) + (-1.0) * B[i][j] * gcovyz[i][j] * gcovyz_px[i][j] +
                          gcovxz[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j]) +
                          B_py[i][j] * gcovxy[i][j] * gcovzz[i][j] + (-1.0) * B[i][j] * gcovxy_py[i][j] * gcovzz[i][j] +
                          (-1.0) * B_px[i][j] * gcovyy[i][j] * gcovzz[i][j] +
                          B[i][j] * gcovyy_px[i][j] * gcovzz[i][j]) +
                     B[i][j] * (Bny_py[i][j] * gcovxz[i][j] * gcovyz[i][j] +
                                (-1.0) * Bny_py[i][j] * gcovxy[i][j] * gcovzz[i][j] +
                                Bny_px[i][j] * ((-1.0) * std::pow(gcovyz[i][j], 2.0) + gcovyy[i][j] * gcovzz[i][j]))) *
                    std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][6] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                        ((-1.0) * std::pow(gcovyz[i][j], 2.0) + gcovyy[i][j] * gcovzz[i][j]) *
                                        std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][7] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                        (gcovxz[i][j] * gcovyz[i][j] + (-1.0) * gcovxy[i][j] * gcovzz[i][j]) *
                                        std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][8] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                        ((-1.0) * gcovxz[i][j] * gcovyy[i][j] + gcovxy[i][j] * gcovyz[i][j]) *
                                        std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][9] =
                    std::pow(B[i][j], -4.0) *
                    (B[i][j] * ((-1.0) * Bny_px[i][j] * gcovxz[i][j] * gcovyy[i][j] +
                                (-1.0) * Bny_py[i][j] * gcovxx[i][j] * gcovyz[i][j] +
                                gcovxy[i][j] * (Bny_py[i][j] * gcovxz[i][j] + Bny_px[i][j] * gcovyz[i][j])) +
                     Bny[i][j] * (gcovxz[i][j] * (B_px[i][j] * gcovyy[i][j] +
                                                  B[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j])) +
                                  (-1.0) * gcovxy[i][j] *
                                      (B_py[i][j] * gcovxz[i][j] + B_px[i][j] * gcovyz[i][j] +
                                       (-1.0) * B[i][j] * gcovyz_px[i][j]) +
                                  gcovxx[i][j] * (B_py[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_py[i][j]))) *
                    std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][10] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                         ((-1.0) * gcovxz[i][j] * gcovyy[i][j] + gcovxy[i][j] * gcovyz[i][j]) *
                                         std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][11] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                         (gcovxy[i][j] * gcovxz[i][j] + (-1.0) * gcovxx[i][j] * gcovyz[i][j]) *
                                         std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][12] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                         ((-1.0) * std::pow(gcovxy[i][j], 2.0) + gcovxx[i][j] * gcovyy[i][j]) *
                                         std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][13] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                         ((gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]) * gcovyz[i][j] +
                                          gcovyy[i][j] * gcovyz_px[i][j] + (-1.0) * gcovxy[i][j] * gcovyz_py[i][j]) *
                                         std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][14] =
                    std::pow(B[i][j], -4.0) *
                    (B[i][j] *
                         (Bny_px[i][j] * gcovxz[i][j] * gcovyy[i][j] + Bny_py[i][j] * gcovxx[i][j] * gcovyz[i][j] +
                          (-1.0) * gcovxy[i][j] * (Bny_py[i][j] * gcovxz[i][j] + Bny_px[i][j] * gcovyz[i][j])) +
                     Bny[i][j] *
                         ((-1.0) * gcovxz[i][j] *
                              (B_px[i][j] * gcovyy[i][j] + B[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j])) +
                          gcovxy[i][j] * (B_py[i][j] * gcovxz[i][j] + B_px[i][j] * gcovyz[i][j] +
                                          (-1.0) * B[i][j] * gcovyz_px[i][j]) +
                          gcovxx[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j]))) *
                    std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][15] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                         (gcovxz[i][j] * gcovyy[i][j] + (-1.0) * gcovxy[i][j] * gcovyz[i][j]) *
                                         std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][16] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                         ((-1.0) * gcovxy[i][j] * gcovxz[i][j] + gcovxx[i][j] * gcovyz[i][j]) *
                                         std::pow(J[i][j], -2.0);
                h_PhiTeA_dTe[j][i][17] = std::pow(B[i][j], -3.0) * Bny[i][j] *
                                         (std::pow(gcovxy[i][j], 2.0) + (-1.0) * gcovxx[i][j] * gcovyy[i][j]) *
                                         std::pow(J[i][j], -2.0);
            }
        }

        if (hostId == 0) {
            logDone();
        }
    }
    void compressStaggeredCoefficient() {

        logStart("Compress staggered coefficient in shifted metric coordinate.");

        /*---------------------------Linear---------------------------*/

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                // Perturbed Parallel Current

                h_A_dJpB[j][i] = std::pow(B[i][j], -2.0) * std::pow(J[i][j], -1.0) *
                                 ((B_px2[i][j] * gconxx[i][j] + B_px[i][j] * gconxx_px[i][j] +
                                   2.0 * B_pxy[i][j] * gconxy[i][j] + B_px[i][j] * gconxy_py[i][j] +
                                   B_py2[i][j] * gconyy[i][j] + B_py[i][j] * (gconxy_px[i][j] + gconyy_py[i][j])) *
                                      J[i][j] +
                                  B_px[i][j] * gconxx[i][j] * J_px[i][j] + B_py[i][j] * gconxy[i][j] * J_px[i][j] +
                                  B_px[i][j] * gconxy[i][j] * J_py[i][j] + B_py[i][j] * gconyy[i][j] * J_py[i][j]);
                h_A_px_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * std::pow(J[i][j], -1.0) *
                                    ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] + gconxx[i][j] * J_px[i][j] +
                                     gconxy[i][j] * J_py[i][j]);
                h_A_pz_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * std::pow(J[i][j], -1.0) *
                                    ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] + gconxz[i][j] * J_px[i][j] +
                                     gconyz[i][j] * J_py[i][j]);
                h_A_px2_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * gconxx[i][j];
                h_A_pxz_dJpB[j][i] = (-2.0) * std::pow(B[i][j], -1.0) * gconxz[i][j];
                h_A_pz2_dJpB[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * gconzz[i][j];
                h_A2dJpB[j][i][0] = h_A_dJpB[j][i];
                h_A2dJpB[j][i][1] = h_A_px_dJpB[j][i];
                h_A2dJpB[j][i][2] = h_A_pz_dJpB[j][i];
                h_A2dJpB[j][i][3] = h_A_px2_dJpB[j][i];
                h_A2dJpB[j][i][4] = h_A_pxz_dJpB[j][i];
                h_A2dJpB[j][i][5] = h_A_pz2_dJpB[j][i];

                // Parallel Resistive

                h_A_resistive[j][i] = h_A_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef + 1.0;
                h_A_px_resistive[j][i] = h_A_px_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_pz_resistive[j][i] = h_A_pz_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_px2_resistive[j][i] = h_A_px2_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_pxz_resistive[j][i] = h_A_pxz_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;
                h_A_pz2_resistive[j][i] = h_A_pz2_dJpB[j][i] * B[i][j] * dt * nablaPerp2A.coef;

                // Perturbed Parallel Vector Potential

                h_Phi_A[j][i] = 0.0;
                h_Phi_px_A[j][i] = 0.0;
                h_Phi_py_A[j][i] = (-1.0) * std::pow(B[i][j], -1.0) * Bny[i][j];
                h_Phi_pz_A[j][i] = 0.0;

                h_dNe_A[j][i] = 0.0;
                h_dNe_px_A[j][i] = 0.0;
                h_dNe_py_A[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -1.0) * Bny[i][j] * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j] * std::pow(NormQE, -1.0);
                h_dNe_pz_A[j][i] = 0.0;

                h_A_A[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                              (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                               Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                              std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                              std::pow(NormQE, -1.0);
                h_A_px_A[j][i] = 0.0;
                h_A_py_A[j][i] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                 std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                                 std::pow(NormQE, -1.0);
                h_A_pz_A[j][i] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                 std::pow(J[i][j], -1.0) * std::pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] *
                                 std::pow(NormQE, -1.0);

                h_APhidNe2A[j][i][0] = h_A_A[j][i];
                h_APhidNe2A[j][i][1] = h_A_py_A[j][i];
                h_APhidNe2A[j][i][2] = h_A_pz_A[j][i];
                h_APhidNe2A[j][i][3] = h_Phi_py_A[j][i];
                h_APhidNe2A[j][i][4] = h_dNe_py_A[j][i];
            }
        }

        /*-------------------------Nonlinear--------------------------*/

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                // Perturbed Parallel Vector Potential

                h_PhiA_A[j][i][0] = std::pow(B[i][j], -3.0) *
                                    (Bny[i][j] * B_py[i][j] * gcovyz[i][j] +
                                     (-1.0) * B[i][j] * (Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_py[i][j])) *
                                    std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][1] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][2] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][3] = std::pow(B[i][j], -3.0) *
                                    (B[i][j] * Bny_px[i][j] * gcovyz[i][j] +
                                     Bny[i][j] * ((-1.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) *
                                    std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][4] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][5] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][6] = std::pow(B[i][j], -3.0) *
                                    (Bny[i][j] * ((-1.0) * B_py[i][j] * gcovxy[i][j] + B_px[i][j] * gcovyy[i][j]) +
                                     B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j] +
                                                Bny[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) *
                                    std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][7] =
                    (-1.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * std::pow(J[i][j], -1.0);
                h_PhiA_A[j][i][8] = std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * std::pow(J[i][j], -1.0);

                h_NeA_A[j][i][0] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                                   (B[i][j] * Bny_py[i][j] * gcovyz[i][j] +
                                    Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][1] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][2] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][3] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                                   (Bny[i][j] * B_px[i][j] * gcovyz[i][j] +
                                    (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][4] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][5] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][6] = (1.0 / 2.0) * std::pow(B[i][j], -3.0) *
                                   (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) +
                                    B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] +
                                               Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][7] = (1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
                h_NeA_A[j][i][8] = (-1.0 / 2.0) * std::pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] *
                                   std::pow(J[i][j], -1.0) * std::pow(NormQE, -1.0) * std::pow(Ne[i][j], -1.0) *
                                   Te[i][j];
            }
        }

        if (hostId == 0) {
            logDone();
        }

        logStart("Copy staggered coefficient to device.");

        Allocator H2DAllocator;

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));

            /*---------------------------Linear---------------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      // Perturbed Parallel Vector Potential in Perturbed Parallel Current
                                      d_A_dJpB[i], h_A_dJpB[0], d_A_px_dJpB[i], h_A_px_dJpB[0], d_A_pz_dJpB[i],
                                      h_A_pz_dJpB[0], d_A_px2_dJpB[i], h_A_px2_dJpB[0], d_A_pxz_dJpB[i],
                                      h_A_pxz_dJpB[0], d_A_pz2_dJpB[i], h_A_pz2_dJpB[0]);

            H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
                                      // Perturbed Electric Potential in Perturbed Parallel Vector Potential
                                      d_Phi_A[i], h_Phi_A[0], d_Phi_px_A[i], h_Phi_px_A[0], d_Phi_py_A[i],
                                      h_Phi_py_A[0], d_Phi_pz_A[i], h_Phi_pz_A[0],
                                      // Perturbed Density in Perturbed Parallel Vector Potential
                                      d_dNe_A[i], h_dNe_A[0], d_dNe_px_A[i], h_dNe_px_A[0], d_dNe_py_A[i],
                                      h_dNe_py_A[0], d_dNe_pz_A[i], h_dNe_pz_A[0],
                                      // Perturbed Parallel Vector Potential in Perturbed Parallel Vector Potential
                                      d_A_A[i], h_A_A[0], d_A_px_A[i], h_A_px_A[0], d_A_py_A[i], h_A_py_A[0],
                                      d_A_pz_A[i], h_A_pz_A[0]);

            H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_A2dJpB[i],
                                      h_A2dJpB[0][0]);
            H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_APhidNe2A[i],
                                      h_APhidNe2A[0][0]);

            /*-------------------------Nonlinear--------------------------*/

            H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_PhiA_A[i],
                                      h_PhiA_A[0][0], d_NeA_A[i], h_NeA_A[0][0]);
        }

        if (hostId == 0) {
            logDone();
        }
    }
    template <int matrixType, typename innerDirichlet, typename outerDirichlet>
    void computeSparseMatrix() {

        constexpr const char* matrixName = matrixType == 0   ? "Poisson"
                                           : matrixType == 1 ? "Resistive"
                                           : matrixType == 2 ? "NablaPerp2Phi"
                                           : matrixType == 3 ? "NablaPerp2dNe"
                                           : matrixType == 4 ? "NablaPerp2dTe"
                                           : matrixType == 5 ? "NablaPerp2dPi"
                                           : matrixType == 6 ? "NablaPerp2dPa"
                                                             : "NablaPerp2dPb";

        auto&& [csrR, csrC, csrV] = [&]() {
            if constexpr (matrixType == 0)
                return std::tie(d_laplacianCsrR, d_laplacianCsrC, d_laplacianCsrV);
            else if constexpr (matrixType == 1)
                return std::tie(d_resistiveCsrR, d_resistiveCsrC, d_resistiveCsrV);
            else if constexpr (matrixType == 2)
                return std::tie(d_PhiCsrR, d_PhiCsrC, d_PhiCsrV);
            else if constexpr (matrixType == 3)
                return std::tie(d_dNeCsrR, d_dNeCsrC, d_dNeCsrV);
            else if constexpr (matrixType == 4)
                return std::tie(d_dTeCsrR, d_dTeCsrC, d_dTeCsrV);
            else if constexpr (matrixType == 5)
                return std::tie(d_dPiCsrR, d_dPiCsrC, d_dPiCsrV);
            else if constexpr (matrixType == 6)
                return std::tie(d_dPaCsrR, d_dPaCsrC, d_dPaCsrV);
            else
                return std::tie(d_dPbCsrR, d_dPbCsrC, d_dPbCsrV);
        }();

        auto&& [Configs, Datas, As, Xs, Bs] = [&]() {
            if constexpr (matrixType == 0)
                return std::tie(laplacianConfigs, laplacianDatas, laplacianAs, laplacianXs, laplacianBs);
            else if constexpr (matrixType == 1)
                return std::tie(resistiveConfigs, resistiveDatas, resistiveAs, resistiveXs, resistiveBs);
            else if constexpr (matrixType == 2)
                return std::tie(PhiConfigs, PhiDatas, PhiAs, PhiXs, PhiBs);
            else if constexpr (matrixType == 3)
                return std::tie(dNeConfigs, dNeDatas, dNeAs, dNeXs, dNeBs);
            else if constexpr (matrixType == 4)
                return std::tie(dTeConfigs, dTeDatas, dTeAs, dTeXs, dTeBs);
            else if constexpr (matrixType == 5)
                return std::tie(dPiConfigs, dPiDatas, dPiAs, dPiXs, dPiBs);
            else if constexpr (matrixType == 6)
                return std::tie(dPaConfigs, dPaDatas, dPaAs, dPaXs, dPaBs);
            else
                return std::tie(dPbConfigs, dPbDatas, dPbAs, dPbXs, dPbBs);
        }();

        auto& xField = [&]() -> mhdReal**& {
            if constexpr (matrixType == 0)
                return d_Phi_midl;
            else if constexpr (matrixType == 1)
                return d_A_midr;
            else if constexpr (matrixType == 2)
                return d_Phi_midr;
            else if constexpr (matrixType == 3)
                return d_dNe_midr;
            else if constexpr (matrixType == 4)
                return d_dTe_midr;
            else if constexpr (matrixType == 5)
                return d_dPi_midr;
            else if constexpr (matrixType == 6)
                return d_dPa_midr;
            else
                return d_dPb_midr;
        }();
        auto& bField = [&]() -> mhdReal**& {
            if constexpr (matrixType == 0)
                return d_w_midl;
            else if constexpr (matrixType == 1)
                return d_A_midl;
            else if constexpr (matrixType == 2)
                return d_Phi_midl;
            else if constexpr (matrixType == 3)
                return d_dNe_midl;
            else if constexpr (matrixType == 4)
                return d_dTe_midl;
            else if constexpr (matrixType == 5)
                return d_dPi_midl;
            else if constexpr (matrixType == 6)
                return d_dPa_midl;
            else
                return d_dPb_midl;
        }();

        constexpr auto realType = std::is_same_v<mhdReal, double> ? CUDA_R_64F : CUDA_R_32F;

        auto logPhase = [&](const char* verb, const char* kind) {
            if (hostId == 0) {
                std::string msg = std::string(verb) + " " + matrixName + " " + kind + ".";
                logStart(msg.c_str());
            }
        };

        struct Coef {
            double cc, cx, cz, cx2, cxz, cz2;
        };

        auto loadCoefs = [&](int j, int i) -> Coef {
            if constexpr (matrixType == 0) {
                return {
                    0.0, h_Phi_px_w[j][i], h_Phi_pz_w[j][i], h_Phi_px2_w[j][i], h_Phi_pxz_w[j][i], h_Phi_pz2_w[j][i]};
            } else if constexpr (matrixType == 1) {
                return {h_A_resistive[j][i],     h_A_px_resistive[j][i],  h_A_pz_resistive[j][i],
                        h_A_px2_resistive[j][i], h_A_pxz_resistive[j][i], h_A_pz2_resistive[j][i]};
            } else {
                double coef;
                if constexpr (matrixType == 2)
                    coef = nablaPerp2Phi.coef;
                else if constexpr (matrixType == 3)
                    coef = nablaPerp2dNe.coef;
                else if constexpr (matrixType == 4)
                    coef = nablaPerp2dTe.coef;
                else if constexpr (matrixType == 5)
                    coef = nablaPerp2dPi.coef;
                else if constexpr (matrixType == 6)
                    coef = nablaPerp2dPa.coef;
                else
                    coef = nablaPerp2dPb.coef;
                return {h_F_perp2[j][i],
                        -h_F_px_perp2[j][i] * dt * coef,
                        -h_F_pz_perp2[j][i] * dt * coef,
                        -h_F_px2_perp2[j][i] * dt * coef,
                        -h_F_pxz_perp2[j][i] * dt * coef,
                        -h_F_pz2_perp2[j][i] * dt * coef};
            }
        };

        logPhase("Compute", "COO");

        matrix_i.clear();
        matrix_j.clear();
        matrix_v.clear();

        for (int j = hostId * hostNy; j < (hostId + 1) * hostNy; j++) {

            for (int i = 0; i < gridNx; i++) {
                for (int k = 0; k < gridNz; k++) {

                    int row_index = i * gridNz + k;

                    auto pushCoefs = [&](int i_offset, int k_offset, double value) {
                        int k_index = k + k_offset;
                        if (k_index < 0)
                            k_index += gridNz;
                        else if (k_index >= gridNz)
                            k_index -= gridNz;
                        matrix_i.emplace_back(row_index);
                        matrix_j.emplace_back((i + i_offset) * gridNz + k_index);
                        matrix_v.emplace_back(value);
                    };

                    if (i == 0) {

                        if constexpr (std::is_same_v<innerDirichlet, std::integral_constant<bool, true>>) {

                            pushCoefs(0, 0, 1.0);

                        } else {

                            const double coes[5] = {-25.0 / 12.0 / gridDx, 48.0 / 12.0 / gridDx, -36.0 / 12.0 / gridDx,
                                                    16.0 / 12.0 / gridDx, -3.0 / 12.0 / gridDx};
                            int idx = 0;
                            for (int i_offset = 0; i_offset <= 4; i_offset++) {
                                pushCoefs(i_offset, 0, coes[idx++]);
                            }
                        }

                    } else if (i == gridNx - 1) {

                        if constexpr (std::is_same_v<outerDirichlet, std::integral_constant<bool, true>>) {

                            pushCoefs(0, 0, 1.0);

                        } else {

                            const double coes[5] = {3.0 / 12.0 / gridDx, -16.0 / 12.0 / gridDx, 36.0 / 12.0 / gridDx,
                                                    -48.0 / 12.0 / gridDx, 25.0 / 12.0 / gridDx};
                            int idx = 0;
                            for (int i_offset = -4; i_offset <= 0; i_offset++) {
                                pushCoefs(i_offset, 0, coes[idx++]);
                            }
                        }

                    } else if (i == 1) {

                        auto [cc, cx, cz, cx2, cxz, cz2] = loadCoefs(j, i);
                        double coes[25];

                        coes[0] = (-1.0 / 48.0) * cxz / gridDx / gridDz;
                        coes[1] = (1.0 / 6.0) * cxz / gridDx / gridDz;
                        coes[2] = (1.0 / 12.0) * std::pow(gridDx, -2.0) * (11.0 * cx2 + (-3.0) * cx * gridDx);
                        coes[3] = (-1.0 / 6.0) * cxz / gridDx / gridDz;
                        coes[4] = (1.0 / 48.0) * cxz / gridDx / gridDz;

                        coes[5] = (-1.0 / 72.0) / gridDx * std::pow(gridDz, -2.0) *
                                  (6.0 * cz2 * gridDx + 5.0 * cxz * gridDz + (-6.0) * cz * gridDx * gridDz);
                        coes[6] = (1.0 / 9.0) / gridDx * std::pow(gridDz, -2.0) *
                                  (12.0 * cz2 * gridDx + 5.0 * cxz * gridDz + (-6.0) * cz * gridDx * gridDz);
                        coes[7] = (-5.0 / 6.0) * std::pow(gridDx, -2.0) * (2.0 * cx2 + cx * gridDx) +
                                  (-5.0 / 2.0) * cz2 * std::pow(gridDz, -2.0) + cc;
                        coes[8] = (1.0 / 9.0) / gridDx * std::pow(gridDz, -2.0) *
                                  (12.0 * cz2 * gridDx + (-5.0) * cxz * gridDz + 6.0 * cz * gridDx * gridDz);
                        coes[9] = (-1.0 / 72.0) / gridDx * std::pow(gridDz, -2.0) *
                                  (6.0 * cz2 * gridDx + (-5.0) * cxz * gridDz + 6.0 * cz * gridDx * gridDz);

                        coes[10] = (1.0 / 8.0) * cxz / gridDx / gridDz;
                        coes[11] = (-1.0) * cxz / gridDx / gridDz;
                        coes[12] = (1.0 / 2.0) * std::pow(gridDx, -2.0) * (cx2 + 3.0 * cx * gridDx);
                        coes[13] = cxz / gridDx / gridDz;
                        coes[14] = (-1.0 / 8.0) * cxz / gridDx / gridDz;

                        coes[15] = (-1.0 / 24.0) * cxz / gridDx / gridDz;
                        coes[16] = (1.0 / 3.0) * cxz / gridDx / gridDz;
                        coes[17] = (1.0 / 6.0) * std::pow(gridDx, -2.0) * (2.0 * cx2 + (-3.0) * cx * gridDx);
                        coes[18] = (-1.0 / 3.0) * cxz / gridDx / gridDz;
                        coes[19] = (1.0 / 24.0) * cxz / gridDx / gridDz;

                        coes[20] = (1.0 / 144.0) * cxz / gridDx / gridDz;
                        coes[21] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[22] = (-1.0 / 12.0) * std::pow(gridDx, -2.0) * (cx2 + (-1.0) * cx * gridDx);
                        coes[23] = (1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[24] = (-1.0 / 144.0) * cxz / gridDx / gridDz;

                        int idx = 0;
                        for (int i_offset = -1; i_offset <= 3; i_offset++) {
                            for (int k_offset = -2; k_offset <= 2; k_offset++) {
                                pushCoefs(i_offset, k_offset, coes[idx++]);
                            }
                        }

                    } else if (i == gridNx - 2) {

                        auto [cc, cx, cz, cx2, cxz, cz2] = loadCoefs(j, i);
                        double coes[25];

                        coes[0] = (-1.0 / 144.0) * cxz / gridDx / gridDz;
                        coes[1] = (1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[2] = (-1.0 / 12.0) * std::pow(gridDx, -2.0) * (cx2 + cx * gridDx);
                        coes[3] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[4] = (1.0 / 144.0) * cxz / gridDx / gridDz;

                        coes[5] = (1.0 / 24.0) * cxz / gridDx / gridDz;
                        coes[6] = (-1.0 / 3.0) * cxz / gridDx / gridDz;
                        coes[7] = (1.0 / 6.0) * std::pow(gridDx, -2.0) * (2.0 * cx2 + 3.0 * cx * gridDx);
                        coes[8] = (1.0 / 3.0) * cxz / gridDx / gridDz;
                        coes[9] = (-1.0 / 24.0) * cxz / gridDx / gridDz;

                        coes[10] = (-1.0 / 8.0) * cxz / gridDx / gridDz;
                        coes[11] = cxz / gridDx / gridDz;
                        coes[12] = (1.0 / 2.0) * std::pow(gridDx, -2.0) * (cx2 + (-3.0) * cx * gridDx);
                        coes[13] = (-1.0) * cxz / gridDx / gridDz;
                        coes[14] = (1.0 / 8.0) * cxz / gridDx / gridDz;

                        coes[15] = (1.0 / 72.0) / gridDx * std::pow(gridDz, -2.0) *
                                   ((-6.0) * cz2 * gridDx + 5.0 * cxz * gridDz + 6.0 * cz * gridDx * gridDz);
                        coes[16] = (1.0 / 9.0) / gridDx * std::pow(gridDz, -2.0) *
                                   (12.0 * cz2 * gridDx + (-5.0) * cxz * gridDz + (-6.0) * cz * gridDx * gridDz);
                        coes[17] = (5.0 / 6.0) * ((-2.0) * cx2 * std::pow(gridDx, -2.0) + cx / gridDx +
                                                  (-3.0) * cz2 * std::pow(gridDz, -2.0)) +
                                   cc;
                        coes[18] = (1.0 / 9.0) / gridDx * std::pow(gridDz, -2.0) *
                                   (12.0 * cz2 * gridDx + 5.0 * cxz * gridDz + 6.0 * cz * gridDx * gridDz);
                        coes[19] = (-1.0 / 72.0) / gridDx * std::pow(gridDz, -2.0) *
                                   (6.0 * cz2 * gridDx + 5.0 * cxz * gridDz + 6.0 * cz * gridDx * gridDz);

                        coes[20] = (1.0 / 48.0) * cxz / gridDx / gridDz;
                        coes[21] = (-1.0 / 6.0) * cxz / gridDx / gridDz;
                        coes[22] = (1.0 / 12.0) * std::pow(gridDx, -2.0) * (11.0 * cx2 + 3.0 * cx * gridDx);
                        coes[23] = (1.0 / 6.0) * cxz / gridDx / gridDz;
                        coes[24] = (-1.0 / 48.0) * cxz / gridDx / gridDz;

                        int idx = 0;
                        for (int i_offset = -3; i_offset <= 1; i_offset++) {
                            for (int k_offset = -2; k_offset <= 2; k_offset++) {
                                pushCoefs(i_offset, k_offset, coes[idx++]);
                            }
                        }

                    } else {

                        auto [cc, cx, cz, cx2, cxz, cz2] = loadCoefs(j, i);
                        double coes[25];

                        coes[0] = (1.0 / 144.0) * cxz / gridDx / gridDz;
                        coes[1] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[2] = (-1.0 / 12.0) * std::pow(gridDx, -2.0) * (cx2 + (-1.0) * cx * gridDx);
                        coes[3] = (1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[4] = (-1.0 / 144.0) * cxz / gridDx / gridDz;

                        coes[5] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[6] = (4.0 / 9.0) * cxz / gridDx / gridDz;
                        coes[7] = (1.0 / 3.0) * std::pow(gridDx, -2.0) * (4.0 * cx2 + (-2.0) * cx * gridDx);
                        coes[8] = (-4.0 / 9.0) * cxz / gridDx / gridDz;
                        coes[9] = (1.0 / 18.0) * cxz / gridDx / gridDz;

                        coes[10] = (-1.0 / 12.0) * std::pow(gridDz, -2.0) * (cz2 + (-1.0) * cz * gridDz);
                        coes[11] = (1.0 / 3.0) * std::pow(gridDz, -2.0) * (4.0 * cz2 + (-2.0) * cz * gridDz);
                        coes[12] = (-5.0 / 2.0) * cx2 * std::pow(gridDx, -2.0) +
                                   (-5.0 / 2.0) * cz2 * std::pow(gridDz, -2.0) + cc;
                        coes[13] = (2.0 / 3.0) * std::pow(gridDz, -2.0) * (2.0 * cz2 + cz * gridDz);
                        coes[14] = (-1.0 / 12.0) * std::pow(gridDz, -2.0) * (cz2 + cz * gridDz);

                        coes[15] = (1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[16] = (-4.0 / 9.0) * cxz / gridDx / gridDz;
                        coes[17] = (2.0 / 3.0) * std::pow(gridDx, -2.0) * (2.0 * cx2 + cx * gridDx);
                        coes[18] = (4.0 / 9.0) * cxz / gridDx / gridDz;
                        coes[19] = (-1.0 / 18.0) * cxz / gridDx / gridDz;

                        coes[20] = (-1.0 / 144.0) * cxz / gridDx / gridDz;
                        coes[21] = (1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[22] = (-1.0 / 12.0) * std::pow(gridDx, -2.0) * (cx2 + cx * gridDx);
                        coes[23] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
                        coes[24] = (1.0 / 144.0) * cxz / gridDx / gridDz;

                        int idx = 0;
                        for (int i_offset = -2; i_offset <= 2; i_offset++) {
                            for (int k_offset = -2; k_offset <= 2; k_offset++) {
                                pushCoefs(i_offset, k_offset, coes[idx++]);
                            }
                        }
                    }
                }
            }
        }

        if (hostId == 0) {
            logDone();
        }

        int nnz = matrix_i.size() / hostNy;

        logPhase("Compute", "CSR");

        if (hostId == 0) {
            std::cout << "nnz * hostNy: " << matrix_i.size() << "." << std::endl;
            std::cout << "nnz: " << nnz << "." << std::endl;
        }

        Allocator DeviceAllocator;
        int** d_matrixCooR;
        mhdReal** d_matrixTempV;

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_matrixCooR);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_matrixTempV);

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));
            CUDACHECK(cudaMemcpy(d_matrixCooR[i], matrix_i.data() + i * devNy * nnz, sizeof(int) * devNy * nnz,
                                 cudaMemcpyHostToDevice));
            CUDACHECK(cudaMemcpy(d_matrixTempV[i], matrix_v.data() + i * devNy * nnz, sizeof(mhdReal) * devNy * nnz,
                                 cudaMemcpyHostToDevice));
        }

        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), csrR);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, csrC);
        DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, csrV);

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));
            CUDACHECK(cudaMemcpy(csrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz,
                                 cudaMemcpyHostToDevice));
            CUDACHECK(cudaMemcpy(csrV[i], matrix_v.data() + i * devNy * nnz, sizeof(mhdReal) * devNy * nnz,
                                 cudaMemcpyHostToDevice));
        }

        std::vector<cusparseHandle_t> cusparseHandles(devNums);
        std::vector<std::vector<size_t>> pBufferSize(devNums, std::vector<size_t>(devNy));
        std::vector<std::vector<void*>> pBuffer(devNums, std::vector<void*>(devNy));
        std::vector<std::vector<int*>> permutation(devNums, std::vector<int*>(devNy));
        std::vector<std::vector<cusparseMatDescr_t>> descr(devNums, std::vector<cusparseMatDescr_t>(devNy));
        std::vector<std::vector<cusparseSpVecDescr_t>> sortedSpVec(devNums, std::vector<cusparseSpVecDescr_t>(devNy));
        std::vector<std::vector<cusparseDnVecDescr_t>> unsortedDnVec(devNums, std::vector<cusparseDnVecDescr_t>(devNy));

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));
            CUSPARSECHECK(cusparseCreate(&cusparseHandles[i]));

            for (int j = 0; j < devNy; j++) {

                CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz, nnz, gridNxz,
                                               csrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
                CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
                                                             csrR[i] + j * (gridNxz + 1), csrC[i] + j * nnz,
                                                             &pBufferSize[i][j]));
                CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
                CUDACHECK(cudaMalloc(reinterpret_cast<void**>(&permutation[i][j]), sizeof(int) * nnz));
                CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
                CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
                CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
                                               csrR[i] + j * (gridNxz + 1), csrC[i] + j * nnz, permutation[i][j],
                                               pBuffer[i][j]));
                CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j], csrV[i] + j * nnz,
                                                  CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, realType));
                CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, realType));
                CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));
            }
        }

        DeviceAllocator.releaseDeviceArrays(localId, devNums, d_matrixCooR, d_matrixTempV);

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));
            CUSPARSECHECK(cusparseDestroy(cusparseHandles[i]));

            for (int j = 0; j < devNy; j++) {

                CUDACHECK(cudaFree(pBuffer[i][j]));
                CUDACHECK(cudaFree(permutation[i][j]));
                CUSPARSECHECK(cusparseDestroyMatDescr(descr[i][j]));
                CUSPARSECHECK(cusparseDestroySpVec(sortedSpVec[i][j]));
                CUSPARSECHECK(cusparseDestroyDnVec(unsortedDnVec[i][j]));
            }
        }

        if (hostId == 0) {
            logDone();
        }

        logPhase("Factorize", "CSR");

        if constexpr (matrixType == 0) {
            cudaStreams.resize(devNums, std::vector<cudaStream_t>(devNy));
            cudssHandles.resize(devNums, std::vector<cudssHandle_t>(devNy));
        }
        Configs.resize(devNums, std::vector<cudssConfig_t>(devNy));
        Datas.resize(devNums, std::vector<cudssData_t>(devNy));
        As.resize(devNums, std::vector<cudssMatrix_t>(devNy));
        Xs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
        Bs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));

            for (int j = 0; j < devNy; j++) {

                if constexpr (matrixType == 0) {
                    CUDACHECK(cudaStreamCreate(&cudaStreams[i][j]));
                    CUDSSCHECK(cudssCreate(&cudssHandles[i][j]));
                    CUDSSCHECK(cudssSetStream(cudssHandles[i][j], cudaStreams[i][j]));
                }
                CUDSSCHECK(cudssConfigCreate(&Configs[i][j]));
                CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &Datas[i][j]));

                cudssMatrixType_t mtype = CUDSS_MTYPE_GENERAL;
                cudssMatrixViewType_t mview = CUDSS_MVIEW_FULL;
                cudssIndexBase_t base = CUDSS_BASE_ZERO;

                int64_t nrows = gridNxz;
                int64_t ncols = gridNxz;
                int64_t ld = gridNxz;
                int64_t nrhs = 1;

                CUDSSCHECK(cudssMatrixCreateCsr(&As[i][j], nrows, ncols, nnz, csrR[i] + j * (gridNxz + 1), NULL,
                                                csrC[i] + j * nnz, csrV[i] + j * nnz, CUDA_R_32I, realType, mtype,
                                                mview, base));
                CUDSSCHECK(cudssMatrixCreateDn(&Xs[i][j], nrows, nrhs, ld, xField[i] + (j + gridGhost) * gridNxz,
                                               realType, CUDSS_LAYOUT_COL_MAJOR));
                CUDSSCHECK(cudssMatrixCreateDn(&Bs[i][j], nrows, nrhs, ld, bField[i] + (j + gridGhost) * gridNxz,
                                               realType, CUDSS_LAYOUT_COL_MAJOR));

                CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS, Configs[i][j], Datas[i][j], As[i][j],
                                        Xs[i][j], Bs[i][j]));

                CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION, Configs[i][j], Datas[i][j],
                                        As[i][j], Xs[i][j], Bs[i][j]));
            }
        }

        if (hostId == 0) {
            size_t avail, total, used;
            CUDACHECK(cudaSetDevice(localId * devNums));
            CUDACHECK(cudaMemGetInfo(&avail, &total));
            used = total - avail;
            std::cout << BOLDYELLOW << "Device memory used: " << (double)used / 1024 / 1024 / 1024 << " GB." << RESET
                      << std::endl;
        }

        if (hostId == 0) {
            logDone();
        }
    }

    template <int picType, int disType, int spaceType, int velocityType>
    void loadParticles() {

        if (hostId == 0) {

            const char* label;
            if constexpr (picType == 0)
                label = "thermal ions";
            else if constexpr (picType == 1)
                label = "alpha particles";
            else if constexpr (picType == 2)
                label = "beam particles";

            std::cout << BOLDYELLOW << "Start: Load " << label;
            if constexpr (disType == 0)
                std::cout << BOLDYELLOW << " by Maxwell distribution." << RESET << std::endl;
            else if constexpr (disType == 1)
                std::cout << BOLDYELLOW << " by isotropic slowing-down distribution without erf function." << RESET
                          << std::endl;
            else if constexpr (disType == 2)
                std::cout << BOLDYELLOW << " by isotropic slowing-down distribution with erf function." << RESET
                          << std::endl;
            else if constexpr (disType == 3)
                std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution without erf function." << RESET
                          << std::endl;
            else if constexpr (disType == 4)
                std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution with erf function." << RESET
                          << std::endl;
        }

        double Jmax, Jvmax;
        double Mass, Vmin, Vmax, Vb, DeltaV, Lambda0, DeltaLambda2;

        std::vector<Rand01> xrand(devNums);
        std::vector<Rand01> yrand(devNums);
        std::vector<Rand01> zrand(devNums);
        std::vector<Rand01> Jrand(devNums);

        std::vector<Rand01> vrand(devNums);
        std::vector<Rand01> vprand(devNums);
        std::vector<Rand01> Jvrand(devNums);

        std::vector<RandNormal> vparand(devNums);
        std::vector<Rand01> vperand(devNums);

        std::vector<std::vector<double>> tempJ(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempB(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempN(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempT(gridNx, std::vector<double>(gridNy + 2));

        if constexpr (picType == 0) {

            Mass = IonMass;
            Vmin = IonVmin;
            Vmax = IonVmax;
            Vb = IonVb;
            DeltaV = IonDeltaV;
            Lambda0 = IonLambda0;
            DeltaLambda2 = IonDeltaLambda2;

        } else if constexpr (picType == 1) {

            Mass = AlphaMass;
            Vmin = AlphaVmin;
            Vmax = AlphaVmax;
            Vb = AlphaVb;
            DeltaV = AlphaDeltaV;
            Lambda0 = AlphaLambda0;
            DeltaLambda2 = AlphaDeltaLambda2;

        } else if constexpr (picType == 2) {

            Mass = BeamMass;
            Vmin = BeamVmin;
            Vmax = BeamVmax;
            Vb = BeamVb;
            DeltaV = BeamDeltaV;
            Lambda0 = BeamLambda0;
            DeltaLambda2 = BeamDeltaLambda2;
        }

        int** picOffsets;
        int** picKeys;
        picReal** picValues;

        if constexpr (picType == 0) {
            picOffsets = h_Ion_offsets;
            picKeys = h_Ion_keys;
            picValues = h_Ion_values;
        } else if constexpr (picType == 1) {
            picOffsets = h_Alpha_offsets;
            picKeys = h_Alpha_keys;
            picValues = h_Alpha_values;
        } else if constexpr (picType == 2) {
            picOffsets = h_Beam_offsets;
            picKeys = h_Beam_keys;
            picValues = h_Beam_values;
        }

        Jmax = 0.0;
        Jvmax = std::pow(Vmax, 2.0);

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                tempJ[i][j + 1] = J[i][j];
                tempB[i][j + 1] = B[i][j];

                if constexpr (picType == 0) {
                    tempN[i][j + 1] = Ni[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Ti[i][j];
                } else if constexpr (picType == 1) {
                    tempN[i][j + 1] = Na[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Ta[i][j];
                } else if constexpr (picType == 2) {
                    tempN[i][j + 1] = Nb[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Tb[i][j];
                }

                if constexpr (spaceType == 0) {
                    if (tempJ[i][j + 1] * tempN[i][j + 1] > Jmax)
                        Jmax = tempJ[i][j + 1] * tempN[i][j + 1];
                } else if constexpr (spaceType == 1) {
                    if (tempJ[i][j + 1] > Jmax)
                        Jmax = tempJ[i][j + 1];
                }
            }

            tempJ[i][0] = tempJ[i][gridNy];
            tempJ[i][gridNy + 1] = tempJ[i][1];

            tempB[i][0] = tempB[i][gridNy];
            tempB[i][gridNy + 1] = tempB[i][1];

            tempN[i][0] = tempN[i][gridNy];
            tempN[i][gridNy + 1] = tempN[i][1];

            tempT[i][0] = tempT[i][gridNy];
            tempT[i][gridNy + 1] = tempT[i][1];
        }

        auto interp2d = [&](std::vector<std::vector<double>>& field, double x, double y) {
            double li = (x - x0) / gridDx;
            double lj = (y - y0 + 0.5 * gridDy) / gridDy;
            int i = std::floor(li);
            int j = std::floor(lj);
            double dx = li - i;
            double dy = lj - j;

            double coes[4] = {};
            double cx[4] = {1.0, 1.0, 0.0, 0.0};
            double sx[4] = {-1.0, -1.0, 1.0, 1.0};
            double cy[4] = {1.0, 0.0, 1.0, 0.0};
            double sy[4] = {-1.0, 1.0, -1.0, 1.0};

            double result = 0.0;

            coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
            coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
            coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
            coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

            result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] +
                     field[i + 1][j + 1] * coes[3];

            return result;
        };

#pragma omp parallel for num_threads(devNums)
        for (int devId = 0; devId < devNums; devId++) {
            for (int picId = 0; picId < picDev / gridNz; picId++) {

                int i, j, k;
                double li, lj, lk;
                double J, B, v, Jv, N, T, Lambda;
                double x, y, z, vp, mu, pw;

                if constexpr (spaceType == 0) {
                    do {
                        x = x0 + (x1 - x0) * xrand[devId]();
                        y = y0 + (y1 - y0) * yrand[devId]();
                        z = z0 + gridDz * zrand[devId]();
                        N = interp2d(tempN, x, y);
                        J = interp2d(tempJ, x, y);
                    } while (Jrand[devId]() >= N * J / Jmax);
                } else if constexpr (spaceType == 1) {
                    do {
                        x = x0 + (x1 - x0) * xrand[devId]();
                        y = y0 + (y1 - y0) * yrand[devId]();
                        z = z0 + gridDz * zrand[devId]();
                        J = interp2d(tempJ, x, y);
                    } while (Jrand[devId]() >= J / Jmax);
                }

                B = interp2d(tempB, x, y);
                N = interp2d(tempN, x, y);
                T = interp2d(tempT, x, y);

                if constexpr (velocityType == 0) {

                    if constexpr (disType == 0) {

                        double vth = std::sqrt(T * KEV / (Mass * MP * std::pow(VA0, 2.0)));
                        double vpara, vperp;

                        do {
                            vpara = vth * vparand[devId]();
                            vperp = vth * std::sqrt(-2.0 * std::log(1.0 - vperand[devId]()));
                            if constexpr (picType == 2)
                                vpara = std::fabs(vpara);
                            v = std::sqrt(std::pow(vpara, 2.0) + std::pow(vperp, 2.0));
                        } while (v <= Vmin || v >= Vmax);

                        vp = vpara;
                        mu = 0.5 * Mass * std::pow(vperp, 2.0) / B;
                        Lambda = mu / (0.5 * Mass * std::pow(v, 2.0));

                    } else if constexpr (disType == 1) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);
                        double U = vrand[devId]();
                        v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);

                        double pitch;
                        if constexpr (picType == 2)
                            pitch = vprand[devId]();
                        else
                            pitch = 2.0 * vprand[devId]() - 1.0;

                        vp = v * pitch;
                        Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;

                    } else if constexpr (disType == 2) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);

                        do {
                            double U = vrand[devId]();
                            v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);
                        } while (Jvrand[devId]() >= 0.5 * (1.0 + std::erf((Vb - v) / DeltaV)));

                        double pitch;
                        if constexpr (picType == 2)
                            pitch = vprand[devId]();
                        else
                            pitch = 2.0 * vprand[devId]() - 1.0;

                        vp = v * pitch;
                        Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;

                    } else if constexpr (disType == 3) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);
                        double U = vrand[devId]();
                        v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);

                        double pitch;
                        do {
                            if constexpr (picType == 2)
                                pitch = vprand[devId]();
                            else
                                pitch = 2.0 * vprand[devId]() - 1.0;
                            Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        } while (Jvrand[devId]() >= std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2));

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;

                    } else if constexpr (disType == 4) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);

                        do {
                            double U = vrand[devId]();
                            v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);
                        } while (Jvrand[devId]() >= 0.5 * (1.0 + std::erf((Vb - v) / DeltaV)));

                        double pitch;
                        do {
                            if constexpr (picType == 2)
                                pitch = vprand[devId]();
                            else
                                pitch = 2.0 * vprand[devId]() - 1.0;
                            Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        } while (vperand[devId]() >= std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2));

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;
                    }

                    pw = N;

                } else if constexpr (velocityType == 1) {

                    do {
                        v = Vmax * vrand[devId]();
                        Jv = std::pow(v, 2.0);
                    } while (v <= Vmin || Jvrand[devId]() >= Jv / Jvmax);

                    if constexpr (picType == 2)
                        Lambda = vprand[devId]();
                    else
                        Lambda = 2.0 * vprand[devId]() - 1.0;
                    vp = v * Lambda;
                    mu = 0.5 * Mass * std::pow(v, 2.0) * (1.0 - std::pow(Lambda, 2.0)) / B;
                    Lambda = mu / (0.5 * Mass * std::pow(v, 2.0));

                    if constexpr (disType == 0)
                        pw = N * std::pow(T, -1.5) *
                             std::exp(-0.5 * Mass * std::pow(v, 2.0) * MP * std::pow(VA0, 2.0) / (T * KEV));
                    else if constexpr (disType == 1)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0));
                    else if constexpr (disType == 2)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) * (1.0 + std::erf((Vb - v) / DeltaV));
                    else if constexpr (disType == 3)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) *
                             std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2);
                    else if constexpr (disType == 4)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) *
                             std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2) *
                             (1.0 + std::erf((Vb - v) / DeltaV));
                }

                if constexpr (spaceType == 0)
                    pw /= N;
                else if constexpr (spaceType == 1)
                    pw *= 1.0;

                li = (x - x0PlusGhost) / gridDx;
                lj = (y - y0PlusGhost) / gridDy;
                lk = (z - z0PlusGhost) / gridDz;

                i = std::floor(li);
                j = std::floor(lj);
                k = std::floor(lk);

                for (int repeatId = 0; repeatId < gridNz; repeatId++) {

                    for (int varId = 0; varId < 7; varId++) {
                        picOffsets[devId][varId + 1] = picOffsets[devId][varId] + picDev;
                        picKeys[devId][picId + picDev / gridNz * repeatId + varId * picDev] =
                            j * cellNxz + i * cellNz + k + repeatId;
                    }

                    picValues[devId][picId + picDev / gridNz * repeatId + 0 * picDev] = 0.999999 * x;
                    picValues[devId][picId + picDev / gridNz * repeatId + 1 * picDev] = y;
                    picValues[devId][picId + picDev / gridNz * repeatId + 2 * picDev] = z + repeatId * gridDz;
                    picValues[devId][picId + picDev / gridNz * repeatId + 3 * picDev] = vp;
                    picValues[devId][picId + picDev / gridNz * repeatId + 4 * picDev] = 0.0;
                    picValues[devId][picId + picDev / gridNz * repeatId + 5 * picDev] = pw;
                    picValues[devId][picId + picDev / gridNz * repeatId + 6 * picDev] = mu;
                }
            }
        }
#pragma omp barrier

        if (hostId == 0) {

            logDone();
        }
    }
    void loadParticles(std::string file) {

        logStart("Load particles from PICContinue file.");

        const size_t speciesSize = sizeof(picReal) + sizeof(int) * devNums * 8 + sizeof(int) * devNums * picDev * 7 +
                                   sizeof(picReal) * devNums * picDev * 7;
        const int enabledCount = (ifIon ? 1 : 0) + (ifAlpha ? 1 : 0) + (ifBeam ? 1 : 0);
        const size_t rankSize = enabledCount * speciesSize;

        std::ifstream input(file, std::ios::in | std::ios::binary);
        size_t cursor = hostId * rankSize;

        auto readSection = [&](picReal& Const, int** offsetsH, int** keysH, picReal** valuesH) {
            input.seekg(cursor, std::ios::beg);
            input.read(reinterpret_cast<char*>(&Const), sizeof(picReal));
            input.read(reinterpret_cast<char*>(offsetsH[0]), sizeof(int) * devNums * 8);
            input.read(reinterpret_cast<char*>(keysH[0]), sizeof(int) * devNums * picDev * 7);
            input.read(reinterpret_cast<char*>(valuesH[0]), sizeof(picReal) * devNums * picDev * 7);
            cursor += speciesSize;
        };

        if (ifIon)
            readSection(IonConst, h_Ion_offsets, h_Ion_keys, h_Ion_values);
        if (ifAlpha)
            readSection(AlphaConst, h_Alpha_offsets, h_Alpha_keys, h_Alpha_values);
        if (ifBeam)
            readSection(BeamConst, h_Beam_offsets, h_Beam_keys, h_Beam_values);

        input.close();

        if (hostId == 0) {
            logDone();
        }
    }
    template <int picType>
    void computePhaseSpaceRange() {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Compute phase space range of thermal ions.");
            else if constexpr (picType == 1)
                logStart("Compute phase space range of alpha particles.");
            else if constexpr (picType == 2)
                logStart("Compute phase space range of beam particles.");
        }

        picReal* range;
        double Mass, Char, Vmin, Vmax;

        if constexpr (picType == 0) {
            range = IonEPphiLambda;
            Mass = IonMass;
            Char = IonChar;
            Vmin = IonVmin;
            Vmax = IonVmax;
        } else if constexpr (picType == 1) {
            range = AlphaEPphiLambda;
            Mass = AlphaMass;
            Char = AlphaChar;
            Vmin = AlphaVmin;
            Vmax = AlphaVmax;
        } else if constexpr (picType == 2) {
            range = BeamEPphiLambda;
            Mass = BeamMass;
            Char = BeamChar;
            Vmin = BeamVmin;
            Vmax = BeamVmax;
        }

        const double drho = RHO1 - RHO0;
        const double psitmax = PSITMAX / (B0 * L0 * L0);
        const double cm = VA0 / (L0 * (QE * B0 / MP));

        double minE = 0.5 * Mass * std::pow(Vmin, 2.0);
        double maxE = 0.5 * Mass * std::pow(Vmax, 2.0);

        double minPphi = 20251106.0;
        double maxPphi = -20251106.0;
        double minLambda = 0.0;
        double maxLambda = 0.0;

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                double tempLambda = 1.0 / B[i][j];
                if (tempLambda > maxLambda)
                    maxLambda = tempLambda;

                double base = cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) *
                              SFAcovyz[i][j + gridGhost] / (q[i][j] * J[i][j] * B[i][j]);

                double PphiPlus = base - Char * psip[i][j];
                double PphiMinus = -base - Char * psip[i][j];
                if (PphiPlus > maxPphi)
                    maxPphi = PphiPlus;
                if (PphiMinus < minPphi)
                    minPphi = PphiMinus;
            }
        }

        range[0] = (picReal)minE;
        range[1] = (picReal)maxE;
        range[2] = (picReal)minPphi;
        range[3] = (picReal)maxPphi;
        range[4] = (picReal)minLambda;
        range[5] = (picReal)maxLambda;

        if (hostId == 0) {

            if constexpr (picType == 0)
                std::cout << BOLDYELLOW << "Ion range: ";
            else if constexpr (picType == 1)
                std::cout << BOLDYELLOW << "Alpha range: ";
            else if constexpr (picType == 2)
                std::cout << BOLDYELLOW << "Beam range: ";

            std::cout << std::setprecision(std::numeric_limits<picReal>::digits10) << "E=[" << range[0] << ", "
                      << range[1] << "], "
                      << "Pphi=[" << range[2] << ", " << range[3] << "], "
                      << "Lambda=[" << range[4] << ", " << range[5] << "]." << RESET << std::endl;

            logDone();
        }
    }
    template <int picType>
    void computeEquilibriumPressure(std::string initialDir) {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Compute equilibrium pressure contributed by thermal ions.");
            else if constexpr (picType == 1)
                logStart("Compute equilibrium pressure contributed by alpha particles.");
            else if constexpr (picType == 2)
                logStart("Compute equilibrium pressure contributed by beam particles.");
        }

        picReal* picConst;
        int** picOffsets;
        int** picKeys;
        picReal** picValues;
        picReal Beta;
        size_t intraRankOffset;

        const size_t speciesSize = sizeof(picReal) + sizeof(int) * devNums * 8 + sizeof(int) * devNums * picDev * 7 +
                                   sizeof(picReal) * devNums * picDev * 7;

        if constexpr (picType == 0) {
            picConst = &IonConst;
            picOffsets = h_Ion_offsets;
            picKeys = h_Ion_keys;
            picValues = h_Ion_values;
            Beta = IonBeta;
            intraRankOffset = 0;
        } else if constexpr (picType == 1) {
            picConst = &AlphaConst;
            picOffsets = h_Alpha_offsets;
            picKeys = h_Alpha_keys;
            picValues = h_Alpha_values;
            Beta = AlphaBeta;
            intraRankOffset = ifIon ? speciesSize : 0;
        } else if constexpr (picType == 2) {
            picConst = &BeamConst;
            picOffsets = h_Beam_offsets;
            picKeys = h_Beam_keys;
            picValues = h_Beam_values;
            Beta = BeamBeta;
            intraRankOffset = (ifIon ? speciesSize : 0) + (ifAlpha ? speciesSize : 0);
        }

        picReal coes[8];
        picReal cx[8] = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
        picReal sx[8] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
        picReal cy[8] = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
        picReal sy[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
        picReal cz[8] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        picReal sz[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

        int i, j, k, tileId, cellId;
        picReal li, lj, lk;
        picReal x, y, z, vp, pw, mu;
        picReal dx, dy, dz, J, B, N, P;
        picReal*** localN;
        picReal*** globalN;
        picReal*** localP;
        picReal*** globalP;

        Allocator HostAllocator;
        HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNxPlusGhost, gridNzPlusGhost, localN, globalN);
        HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNxPlusGhost, gridNzPlusGhost, localP, globalP);

        for (int devId = 0; devId < devNums; devId++) {
            for (int picId = 0; picId < picDev; picId++) {

                cellId = picKeys[devId][picId];
                x = picValues[devId][picId + 0 * picDev];
                y = picValues[devId][picId + 1 * picDev];
                z = picValues[devId][picId + 2 * picDev];
                vp = picValues[devId][picId + 3 * picDev];
                pw = picValues[devId][picId + 5 * picDev];
                mu = picValues[devId][picId + 6 * picDev];

                tileId = cellId / cellNz;

                li = (x - x0PlusGhost) / gridDx;
                lj = (y - y0PlusGhost) / gridDy;
                lk = (z - z0PlusGhost) / gridDz;

                i = std::floor(li);
                j = std::floor(lj);
                k = std::floor(lk);

                dx = li - i;
                dy = lj - j;
                dz = lk - k;

                coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
                coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
                coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
                coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

                J = h_pic2d[tileId][0] * coes[0] + h_pic2d[tileId][1] * coes[1] + h_pic2d[tileId][2] * coes[2] +
                    h_pic2d[tileId][3] * coes[3];
                B = h_pic2d[tileId][4] * coes[0] + h_pic2d[tileId][5] * coes[1] + h_pic2d[tileId][6] * coes[2] +
                    h_pic2d[tileId][7] * coes[3];

                coes[4] = coes[0];
                coes[5] = coes[1];
                coes[6] = coes[2];
                coes[7] = coes[3];
                coes[0] *= (cz[0] + sz[0] * dz);
                coes[1] *= (cz[1] + sz[1] * dz);
                coes[2] *= (cz[2] + sz[2] * dz);
                coes[3] *= (cz[3] + sz[3] * dz);
                coes[4] *= (cz[4] + sz[4] * dz);
                coes[5] *= (cz[5] + sz[5] * dz);
                coes[6] *= (cz[6] + sz[6] * dz);
                coes[7] *= (cz[7] + sz[7] * dz);

                N = pw / J * B0 * B0 / 2 / MU0 / (MP * VA0 * VA0);

                if constexpr (picType == 0)
                    P = (IonMass * vp * vp + mu * B) * pw / J / 2;
                else if constexpr (picType == 1)
                    P = (AlphaMass * vp * vp + mu * B) * pw / J / 2;
                else if constexpr (picType == 2)
                    P = (BeamMass * vp * vp + mu * B) * pw / J / 2;

                localN[j][i][k] += N * coes[0];
                localN[j][i + 1][k] += N * coes[1];
                localN[j + 1][i][k] += N * coes[2];
                localN[j + 1][i + 1][k] += N * coes[3];
                localN[j][i][k + 1] += N * coes[4];
                localN[j][i + 1][k + 1] += N * coes[5];
                localN[j + 1][i][k + 1] += N * coes[6];
                localN[j + 1][i + 1][k + 1] += N * coes[7];

                localP[j][i][k] += P * coes[0];
                localP[j][i + 1][k] += P * coes[1];
                localP[j + 1][i][k] += P * coes[2];
                localP[j + 1][i + 1][k] += P * coes[3];
                localP[j][i][k + 1] += P * coes[4];
                localP[j][i + 1][k + 1] += P * coes[5];
                localP[j + 1][i][k + 1] += P * coes[6];
                localP[j + 1][i + 1][k + 1] += P * coes[7];
            }
        }

        if constexpr (std::is_same_v<picReal, double>) {
            MPICHECK(MPI_Allreduce(localN[0][0], globalN[0][0], gridNyPlusGhost * gridNxPlusGhost * gridNzPlusGhost,
                                   MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
        } else {
            MPICHECK(MPI_Allreduce(localN[0][0], globalN[0][0], gridNyPlusGhost * gridNxPlusGhost * gridNzPlusGhost,
                                   MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD));
        }

        if constexpr (std::is_same_v<picReal, double>) {
            MPICHECK(MPI_Allreduce(localP[0][0], globalP[0][0], gridNyPlusGhost * gridNxPlusGhost * gridNzPlusGhost,
                                   MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
        } else {
            MPICHECK(MPI_Allreduce(localP[0][0], globalP[0][0], gridNyPlusGhost * gridNxPlusGhost * gridNzPlusGhost,
                                   MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD));
        }

        for (int j = 0; j < gridNyPlusGhost; j++) {
            for (int k = 0; k < gridNzPlusGhost; k++) {
                globalN[j][0][k] *= 2.0;
                globalN[j][gridNxPlusGhost - 1][k] *= 2.0;
                globalP[j][0][k] *= 2.0;
                globalP[j][gridNxPlusGhost - 1][k] *= 2.0;
            }
            for (int i = 0; i < gridNxPlusGhost; i++) {
                globalN[j][i][gridGhost] += globalN[j][i][gridNz + gridGhost];
                globalN[j][i][gridNz + gridGhost - 1] += globalN[j][i][gridGhost - 1];
                globalP[j][i][gridGhost] += globalP[j][i][gridNz + gridGhost];
                globalP[j][i][gridNz + gridGhost - 1] += globalP[j][i][gridGhost - 1];
            }
        }

        for (int i = 0; i < gridNxPlusGhost; i++) {
            for (int k = 0; k < gridNzPlusGhost; k++) {
                globalN[gridGhost][i][k] += globalN[gridNy + gridGhost][i][k];
                globalN[gridNy + gridGhost - 1][i][k] += globalN[gridGhost - 1][i][k];
                globalP[gridGhost][i][k] += globalP[gridNy + gridGhost][i][k];
                globalP[gridNy + gridGhost - 1][i][k] += globalP[gridGhost - 1][i][k];
            }
        }

        picReal innerP = 0.0;

        for (int j = 0; j < gridNy; j++) {
            for (int k = 0; k < gridNz; k++) {
                innerP += globalP[j + gridGhost][0][k + gridGhost];
            }
        }

        innerP /= gridNy * gridNz;

        *picConst = Beta / innerP;

        std::vector<picReal> density(gridNx);
        std::vector<picReal> pressure(gridNx);

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {
                for (int k = 0; k < gridNz; k++) {

                    globalN[j + gridGhost][i][k + gridGhost] *= *picConst;
                    globalP[j + gridGhost][i][k + gridGhost] *= *picConst * B0 * B0 / 2 / MU0;

                    density[i] += globalN[j + gridGhost][i][k + gridGhost] / (gridNy * gridNz);
                    pressure[i] += globalP[j + gridGhost][i][k + gridGhost] / (gridNy * gridNz);
                }
            }
        }

        HostAllocator.releaseHostArrays(localN, globalN);
        HostAllocator.releaseHostArrays(localP, globalP);

        std::ofstream output;
        std::string fileName;

        if (hostId == 0) {

            if constexpr (picType == 0)
                fileName = initialDir + "/IonDensity.bin";
            else if constexpr (picType == 1)
                fileName = initialDir + "/AlphaDensity.bin";
            else if constexpr (picType == 2)
                fileName = initialDir + "/BeamDensity.bin";

            output.open(fileName.c_str(), std::ios::out | std::ios::binary);
            output.write(reinterpret_cast<char*>(density.data()), sizeof(picReal) * density.size());
            output.close();

            if constexpr (picType == 0)
                fileName = initialDir + "/IonPressure.bin";
            else if constexpr (picType == 1)
                fileName = initialDir + "/AlphaPressure.bin";
            else if constexpr (picType == 2)
                fileName = initialDir + "/BeamPressure.bin";

            output.open(fileName.c_str(), std::ios::out | std::ios::binary);
            output.write(reinterpret_cast<char*>(pressure.data()), sizeof(picReal) * pressure.size());
            output.close();
        }

        const MPI_Datatype picMPIType = std::is_same_v<picReal, double> ? MPI_DOUBLE : MPI_FLOAT;
        const int enabledCount = (ifIon ? 1 : 0) + (ifAlpha ? 1 : 0) + (ifBeam ? 1 : 0);
        const size_t rankSize = enabledCount * speciesSize;
        const size_t offset = hostId * rankSize + intraRankOffset;
        const std::string picFile = initialDir + "/PICContinue_0.bin";

        MPI_File fileHandle;
        MPI_File_open(MPI_COMM_WORLD, picFile.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fileHandle);

        MPI_File_write_at_all(fileHandle, offset, picConst, 1, picMPIType, MPI_STATUS_IGNORE);
        MPI_File_write_at_all(fileHandle, offset + sizeof(picReal), picOffsets[0], devNums * 8, MPI_INT,
                              MPI_STATUS_IGNORE);
        MPI_File_write_at_all(fileHandle, offset + sizeof(picReal) + sizeof(int) * devNums * 8, picKeys[0],
                              devNums * picDev * 7, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_write_at_all(fileHandle,
                              offset + sizeof(picReal) + sizeof(int) * devNums * 8 + sizeof(int) * devNums * picDev * 7,
                              picValues[0], devNums * picDev * 7, picMPIType, MPI_STATUS_IGNORE);

        MPI_File_close(&fileHandle);

        if (hostId == 0) {

            if constexpr (picType == 0)
                std::cout << BOLDYELLOW << "IonConst for computing pressure: " << std::setprecision(10) << IonConst
                          << "." << RESET << std::endl;
            else if constexpr (picType == 1)
                std::cout << BOLDYELLOW << "AlphaConst for computing pressure: " << std::setprecision(10) << AlphaConst
                          << "." << RESET << std::endl;
            else if constexpr (picType == 2)
                std::cout << BOLDYELLOW << "BeamConst for computing pressure: " << std::setprecision(10) << BeamConst
                          << "." << RESET << std::endl;

            logDone();
        }
    }

    template <int picType, int disType, int spaceType, int velocityType, int gridE, int gridPphi, int gridLambda,
              int ppcPhase>
    void computePhaseSpaceF0(std::string initialDir) {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Compute phase space f0 of thermal ions.");
            else if constexpr (picType == 1)
                logStart("Compute phase space f0 of alpha particles.");
            else if constexpr (picType == 2)
                logStart("Compute phase space f0 of beam particles.");
        }

        double minE, maxE, dE;
        double minPphi, maxPphi, dPphi;
        double minLambda, maxLambda, dLambda;
        double Mass, Char, Vmin, Vmax, Vb, DeltaV, Lambda0, DeltaLambda2;

        double drho = RHO1 - RHO0;
        double psitmax = PSITMAX / (B0 * L0 * L0);
        double cm = VA0 / (L0 * (QE * B0 / MP));

        if constexpr (picType == 0) {

            Mass = IonMass;
            Char = IonChar;
            Vmin = IonVmin;
            Vmax = IonVmax;
            Vb = IonVb;
            DeltaV = IonDeltaV;
            Lambda0 = IonLambda0;
            DeltaLambda2 = IonDeltaLambda2;

        } else if constexpr (picType == 1) {

            Mass = AlphaMass;
            Char = AlphaChar;
            Vmin = AlphaVmin;
            Vmax = AlphaVmax;
            Vb = AlphaVb;
            DeltaV = AlphaDeltaV;
            Lambda0 = AlphaLambda0;
            DeltaLambda2 = AlphaDeltaLambda2;

        } else if constexpr (picType == 2) {

            Mass = BeamMass;
            Char = BeamChar;
            Vmin = BeamVmin;
            Vmax = BeamVmax;
            Vb = BeamVb;
            DeltaV = BeamDeltaV;
            Lambda0 = BeamLambda0;
            DeltaLambda2 = BeamDeltaLambda2;
        }

        minE = 0.5 * Mass * std::pow(Vmin, 2.0);
        maxE = 0.5 * Mass * std::pow(Vmax, 2.0);

        minPphi = 20251106;
        maxPphi = -20251106;

        minLambda = 0.0;
        maxLambda = 0.0;

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                double tempPphi, tempLambda;

                tempLambda = 1.0 / B[i][j];
                if (tempLambda > maxLambda)
                    maxLambda = tempLambda;

                tempPphi = cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) *
                               SFAcovyz[i][j + gridGhost] / (q[i][j] * J[i][j] * B[i][j]) -
                           Char * psip[i][j];
                if (tempPphi > maxPphi)
                    maxPphi = tempPphi;

                tempPphi = -cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) *
                               SFAcovyz[i][j + gridGhost] / (q[i][j] * J[i][j] * B[i][j]) -
                           Char * psip[i][j];
                if (tempPphi < minPphi)
                    minPphi = tempPphi;
            }
        }

        dE = (maxE - minE) / (gridE - 1);
        dPphi = (maxPphi - minPphi) / (gridPphi - 1);
        dLambda = (maxLambda - minLambda) / (gridLambda - 1);

        size_t picPhase = (size_t)gridE * gridPphi * gridLambda * ppcPhase / hostNums;

        double*** phaseSpaceF0;
        Allocator HostAllocator;
        HostAllocator.allocateHostArrays(gridE, gridPphi, gridLambda, phaseSpaceF0);

        const size_t f0Len = (size_t)gridE * gridPphi * gridLambda;
        std::vector<std::vector<double>> localF0(devNums, std::vector<double>(f0Len, 0.0));

        double Jmax, Jvmax;
        Jmax = 0.0;
        Jvmax = std::pow(Vmax, 2.0);

        std::vector<std::vector<double>> tempq(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> temppsip(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempJ(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempB(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempN(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempT(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempSFAcovyz(gridNx, std::vector<double>(gridNy + 2));

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                tempq[i][j + 1] = q[i][j];
                temppsip[i][j + 1] = psip[i][j];
                tempJ[i][j + 1] = J[i][j];
                tempB[i][j + 1] = B[i][j];

                if constexpr (picType == 0) {
                    tempN[i][j + 1] = Ni[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Ti[i][j];
                } else if constexpr (picType == 1) {
                    tempN[i][j + 1] = Na[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Ta[i][j];
                } else if constexpr (picType == 2) {
                    tempN[i][j + 1] = Nb[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Tb[i][j];
                }

                if constexpr (spaceType == 0) {
                    if (tempJ[i][j + 1] * tempN[i][j + 1] > Jmax)
                        Jmax = tempJ[i][j + 1] * tempN[i][j + 1];
                } else if constexpr (spaceType == 1) {
                    if (tempJ[i][j + 1] > Jmax)
                        Jmax = tempJ[i][j + 1];
                }
            }

            tempq[i][0] = tempq[i][gridNy];
            tempq[i][gridNy + 1] = tempq[i][1];

            temppsip[i][0] = temppsip[i][gridNy];
            temppsip[i][gridNy + 1] = temppsip[i][1];

            tempJ[i][0] = tempJ[i][gridNy];
            tempJ[i][gridNy + 1] = tempJ[i][1];

            tempB[i][0] = tempB[i][gridNy];
            tempB[i][gridNy + 1] = tempB[i][1];

            tempN[i][0] = tempN[i][gridNy];
            tempN[i][gridNy + 1] = tempN[i][1];

            tempT[i][0] = tempT[i][gridNy];
            tempT[i][gridNy + 1] = tempT[i][1];
        }

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy + 2; j++) {
                tempSFAcovyz[i][j] = SFAcovyz[i][j + gridGhost - 1];
            }
        }

        auto interp2d = [&](std::vector<std::vector<double>>& field, double x, double y) {
            double li = (x - x0) / gridDx;
            double lj = (y - y0 + 0.5 * gridDy) / gridDy;
            int i = std::floor(li);
            int j = std::floor(lj);
            double dx = li - i;
            double dy = lj - j;

            double coes[4] = {};
            double cx[4] = {1.0, 1.0, 0.0, 0.0};
            double sx[4] = {-1.0, -1.0, 1.0, 1.0};
            double cy[4] = {1.0, 0.0, 1.0, 0.0};
            double sy[4] = {-1.0, 1.0, -1.0, 1.0};

            double result = 0.0;

            coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
            coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
            coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
            coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

            result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] +
                     field[i + 1][j + 1] * coes[3];

            return result;
        };

        std::vector<Rand01> xrand(devNums);
        std::vector<Rand01> yrand(devNums);
        std::vector<Rand01> Jrand(devNums);

        std::vector<Rand01> vrand(devNums);
        std::vector<Rand01> vprand(devNums);
        std::vector<Rand01> Jvrand(devNums);

        std::vector<RandNormal> vparand(devNums);
        std::vector<Rand01> vperand(devNums);

        double cx[8] = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
        double sx[8] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
        double cy[8] = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
        double sy[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
        double cz[8] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        double sz[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

#pragma omp parallel for num_threads(devNums)
        for (int devId = 0; devId < devNums; devId++) {

            auto& local = localF0[devId];
            const size_t picNums = picPhase / devNums;

            for (size_t picId = 0; picId < picNums; picId++) {

                if (hostId == 0)
                    if (devId == 0)
                        if (picId % (picNums / 4) == 0)
                            std::cout << BOLDGREEN << 100 * picId / picNums << "%" << RESET << std::endl;

                int i, j, k;
                double li, lj, lk;
                double dx, dy, dz;
                double q, psip, J, B, N, T, SFAcovyz;
                double v, Jv, E, Pphi, Lambda;
                double x, y, vp, mu, pw;
                double coes[8];

                if constexpr (spaceType == 0) {
                    do {
                        x = x0 + (x1 - x0) * xrand[devId]();
                        y = y0 + (y1 - y0) * yrand[devId]();
                        N = interp2d(tempN, x, y);
                        J = interp2d(tempJ, x, y);
                    } while (Jrand[devId]() >= N * J / Jmax);
                } else if constexpr (spaceType == 1) {
                    do {
                        x = x0 + (x1 - x0) * xrand[devId]();
                        y = y0 + (y1 - y0) * yrand[devId]();
                        J = interp2d(tempJ, x, y);
                    } while (Jrand[devId]() >= J / Jmax);
                }

                q = interp2d(tempq, x, y);
                psip = interp2d(temppsip, x, y);
                SFAcovyz = interp2d(tempSFAcovyz, x, y);
                B = interp2d(tempB, x, y);
                N = interp2d(tempN, x, y);
                T = interp2d(tempT, x, y);

                if constexpr (velocityType == 0) {

                    if constexpr (disType == 0) {

                        double vth = std::sqrt(T * KEV / (Mass * MP * std::pow(VA0, 2.0)));
                        double vpara, vperp;

                        do {
                            vpara = vth * vparand[devId]();
                            vperp = vth * std::sqrt(-2.0 * std::log(1.0 - vperand[devId]()));
                            if constexpr (picType == 2)
                                vpara = std::fabs(vpara);
                            v = std::sqrt(std::pow(vpara, 2.0) + std::pow(vperp, 2.0));
                        } while (v <= Vmin || v >= Vmax);

                        vp = vpara;
                        mu = 0.5 * Mass * std::pow(vperp, 2.0) / B;
                        Lambda = mu / (0.5 * Mass * std::pow(v, 2.0));

                    } else if constexpr (disType == 1) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);
                        double U = vrand[devId]();
                        v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);

                        double pitch;
                        if constexpr (picType == 2)
                            pitch = vprand[devId]();
                        else
                            pitch = 2.0 * vprand[devId]() - 1.0;

                        vp = v * pitch;
                        Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;

                    } else if constexpr (disType == 2) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);

                        do {
                            double U = vrand[devId]();
                            v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);
                        } while (Jvrand[devId]() >= 0.5 * (1.0 + std::erf((Vb - v) / DeltaV)));

                        double pitch;
                        if constexpr (picType == 2)
                            pitch = vprand[devId]();
                        else
                            pitch = 2.0 * vprand[devId]() - 1.0;

                        vp = v * pitch;
                        Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;

                    } else if constexpr (disType == 3) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);
                        double U = vrand[devId]();
                        v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);

                        double pitch;
                        do {
                            if constexpr (picType == 2)
                                pitch = vprand[devId]();
                            else
                                pitch = 2.0 * vprand[devId]() - 1.0;
                            Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        } while (Jvrand[devId]() >= std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2));

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;

                    } else if constexpr (disType == 4) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);

                        do {
                            double U = vrand[devId]();
                            v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);
                        } while (Jvrand[devId]() >= 0.5 * (1.0 + std::erf((Vb - v) / DeltaV)));

                        double pitch;
                        do {
                            if constexpr (picType == 2)
                                pitch = vprand[devId]();
                            else
                                pitch = 2.0 * vprand[devId]() - 1.0;
                            Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        } while (vperand[devId]() >= std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2));

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;
                    }

                    pw = N;

                } else if constexpr (velocityType == 1) {

                    do {
                        v = Vmax * vrand[devId]();
                        Jv = std::pow(v, 2.0);
                    } while (v <= Vmin || Jvrand[devId]() >= Jv / Jvmax);

                    if constexpr (picType == 2)
                        Lambda = vprand[devId]();
                    else
                        Lambda = 2.0 * vprand[devId]() - 1.0;
                    vp = v * Lambda;
                    mu = 0.5 * Mass * std::pow(v, 2.0) * (1.0 - std::pow(Lambda, 2.0)) / B;
                    Lambda = (1.0 - std::pow(Lambda, 2.0)) / B;

                    if constexpr (disType == 0)
                        pw = N * std::pow(T, -1.5) *
                             std::exp(-0.5 * Mass * std::pow(v, 2.0) * MP * std::pow(VA0, 2.0) / (T * KEV));
                    else if constexpr (disType == 1)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0));
                    else if constexpr (disType == 2)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) * (1.0 + std::erf((Vb - v) / DeltaV));
                    else if constexpr (disType == 3)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) *
                             std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2);
                    else if constexpr (disType == 4)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) *
                             std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2) *
                             (1.0 + std::erf((Vb - v) / DeltaV));
                }

                E = 0.5 * Mass * std::pow(v, 2.0);
                Pphi = cm * Mass * vp * 2 * psitmax * drho * (RHO0 + x * drho) * SFAcovyz / (q * J * B) - Char * psip;
                Lambda = mu / E;

                if constexpr (spaceType == 0)
                    pw /= N;
                else if constexpr (spaceType == 1)
                    pw *= 1.0;

                li = (E - minE) / dE;
                lj = (Pphi - minPphi) / dPphi;
                lk = (Lambda - minLambda) / dLambda;

                i = std::floor(li);
                j = std::floor(lj);
                k = std::floor(lk);

                dx = li - i;
                dy = lj - j;
                dz = lk - k;

                if (i == gridE - 1) {
                    i--;
                    dx = 1;
                }
                if (j == gridPphi - 1) {
                    j--;
                    dy = 1;
                }
                if (k == gridLambda - 1) {
                    k--;
                    dz = 1;
                }

                for (int index = 0; index < 8; index++)
                    coes[index] =
                        (cx[index] + sx[index] * dx) * (cy[index] + sy[index] * dy) * (cz[index] + sz[index] * dz);

                const size_t strideE = (size_t)gridPphi * gridLambda;
                const size_t strideP = gridLambda;
                const size_t index = (size_t)i * strideE + (size_t)j * strideP + k;

                local[index] += pw * coes[0];
                local[index + strideE] += pw * coes[1];
                local[index + strideP] += pw * coes[2];
                local[index + strideE + strideP] += pw * coes[3];
                local[index + 1] += pw * coes[4];
                local[index + strideE + 1] += pw * coes[5];
                local[index + strideP + 1] += pw * coes[6];
                local[index + strideE + strideP + 1] += pw * coes[7];
            }
        }

        if (hostId == 0)
            std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;

        double* f0 = phaseSpaceF0[0][0];

#pragma omp parallel for num_threads(devNums)
        for (size_t index = 0; index < f0Len; index++) {
            double sum = 0.0;
            for (int devId = 0; devId < devNums; devId++)
                sum += localF0[devId][index];
            f0[index] = sum;
        }

        MPICHECK(MPI_Allreduce(MPI_IN_PLACE, phaseSpaceF0[0][0], gridE * gridPphi * gridLambda, MPI_DOUBLE, MPI_SUM,
                               MPI_COMM_WORLD));

        for (int i = 0; i < gridE; i++) {
            for (int j = 0; j < gridPphi; j++) {
                for (int k = 0; k < gridLambda; k++) {

                    bool isVertex =
                        (i == 0 || i == gridE - 1) && (j == 0 || j == gridPphi - 1) && (k == 0 || k == gridLambda - 1);

                    bool isEdge = ((i == 0 || i == gridE - 1) && (j == 0 || j == gridPphi - 1)) ||
                                  ((i == 0 || i == gridE - 1) && (k == 0 || k == gridLambda - 1)) ||
                                  ((j == 0 || j == gridPphi - 1) && (k == 0 || k == gridLambda - 1));

                    bool isFace =
                        (i == 0 || i == gridE - 1 || j == 0 || j == gridPphi - 1 || k == 0 || k == gridLambda - 1);

                    if (isVertex)
                        phaseSpaceF0[i][j][k] *= 8;
                    else if (isEdge)
                        phaseSpaceF0[i][j][k] *= 4;
                    else if (isFace)
                        phaseSpaceF0[i][j][k] *= 2;
                }
            }
        }

        if (hostId == 0) {

            std::ofstream output;
            std::string fileName;

            if constexpr (picType == 0)
                fileName = initialDir + "/IonPhaseSpaceF0.bin";
            else if constexpr (picType == 1)
                fileName = initialDir + "/AlphaPhaseSpaceF0.bin";
            else if constexpr (picType == 2)
                fileName = initialDir + "/BeamPhaseSpaceF0.bin";

            output.open(fileName.c_str(), std::ios::out | std::ios::binary);

            output.write(reinterpret_cast<char*>(phaseSpaceF0[0][0]), sizeof(double) * gridE * gridPphi * gridLambda);

            output.close();

            logDone();
        }

        HostAllocator.releaseHostArrays(phaseSpaceF0);
    }
    template <int picType, int disType, int spaceType, int velocityType, int gridVpara, int gridVperp, int ppcPitch>
    void computePitchSpaceF0(std::string initialDir) {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Compute pitch space f0 of thermal ions.");
            else if constexpr (picType == 1)
                logStart("Compute pitch space f0 of alpha particles.");
            else if constexpr (picType == 2)
                logStart("Compute pitch space f0 of beam particles.");
        }

        double minVpara, maxVpara, dVpara;
        double minVperp, maxVperp, dVperp;
        double Mass, Vmin, Vmax, Vb, DeltaV, Lambda0, DeltaLambda2;

        if constexpr (picType == 0) {

            Mass = IonMass;
            Vmin = IonVmin;
            Vmax = IonVmax;
            Vb = IonVb;
            DeltaV = IonDeltaV;
            Lambda0 = IonLambda0;
            DeltaLambda2 = IonDeltaLambda2;

        } else if constexpr (picType == 1) {

            Mass = AlphaMass;
            Vmin = AlphaVmin;
            Vmax = AlphaVmax;
            Vb = AlphaVb;
            DeltaV = AlphaDeltaV;
            Lambda0 = AlphaLambda0;
            DeltaLambda2 = AlphaDeltaLambda2;

        } else if constexpr (picType == 2) {

            Mass = BeamMass;
            Vmin = BeamVmin;
            Vmax = BeamVmax;
            Vb = BeamVb;
            DeltaV = BeamDeltaV;
            Lambda0 = BeamLambda0;
            DeltaLambda2 = BeamDeltaLambda2;
        }

        minVpara = (picType == 2) ? 0.0 : -Vmax;
        maxVpara = Vmax;
        minVperp = 0.0;
        maxVperp = Vmax;

        dVpara = (maxVpara - minVpara) / (gridVpara - 1);
        dVperp = (maxVperp - minVperp) / (gridVperp - 1);

        size_t picPitch = (size_t)gridVpara * gridVperp * ppcPitch / hostNums;

        double** pitchSpaceF0;
        Allocator HostAllocator;
        HostAllocator.allocateHostArrays(gridVpara, gridVperp, pitchSpaceF0);

        const size_t f0Len = (size_t)gridVpara * gridVperp;
        std::vector<std::vector<double>> localF0(devNums, std::vector<double>(f0Len, 0.0));

        double Jmax, Jvmax;
        Jmax = 0.0;
        Jvmax = std::pow(Vmax, 2.0);

        std::vector<std::vector<double>> tempJ(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempB(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempN(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempT(gridNx, std::vector<double>(gridNy + 2));

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                tempJ[i][j + 1] = J[i][j];
                tempB[i][j + 1] = B[i][j];

                if constexpr (picType == 0) {
                    tempN[i][j + 1] = Ni[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Ti[i][j];
                } else if constexpr (picType == 1) {
                    tempN[i][j + 1] = Na[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Ta[i][j];
                } else if constexpr (picType == 2) {
                    tempN[i][j + 1] = Nb[i][j] * 1.0e-19;
                    tempT[i][j + 1] = Tb[i][j];
                }

                if constexpr (spaceType == 0) {
                    if (tempJ[i][j + 1] * tempN[i][j + 1] > Jmax)
                        Jmax = tempJ[i][j + 1] * tempN[i][j + 1];
                } else if constexpr (spaceType == 1) {
                    if (tempJ[i][j + 1] > Jmax)
                        Jmax = tempJ[i][j + 1];
                }
            }

            tempJ[i][0] = tempJ[i][gridNy];
            tempJ[i][gridNy + 1] = tempJ[i][1];

            tempB[i][0] = tempB[i][gridNy];
            tempB[i][gridNy + 1] = tempB[i][1];

            tempN[i][0] = tempN[i][gridNy];
            tempN[i][gridNy + 1] = tempN[i][1];

            tempT[i][0] = tempT[i][gridNy];
            tempT[i][gridNy + 1] = tempT[i][1];
        }

        auto interp2d = [&](std::vector<std::vector<double>>& field, double x, double y) {
            double li = (x - x0) / gridDx;
            double lj = (y - y0 + 0.5 * gridDy) / gridDy;
            int i = std::floor(li);
            int j = std::floor(lj);
            double dx = li - i;
            double dy = lj - j;

            double coes[4] = {};
            double cx[4] = {1.0, 1.0, 0.0, 0.0};
            double sx[4] = {-1.0, -1.0, 1.0, 1.0};
            double cy[4] = {1.0, 0.0, 1.0, 0.0};
            double sy[4] = {-1.0, 1.0, -1.0, 1.0};

            double result = 0.0;

            coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
            coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
            coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
            coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

            result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] +
                     field[i + 1][j + 1] * coes[3];

            return result;
        };

        std::vector<Rand01> xrand(devNums);
        std::vector<Rand01> yrand(devNums);
        std::vector<Rand01> Jrand(devNums);

        std::vector<Rand01> vrand(devNums);
        std::vector<Rand01> vprand(devNums);
        std::vector<Rand01> Jvrand(devNums);

        std::vector<RandNormal> vparand(devNums);
        std::vector<Rand01> vperand(devNums);

        double cx[4] = {1.0, 0.0, 1.0, 0.0};
        double sx[4] = {-1.0, 1.0, -1.0, 1.0};
        double cy[4] = {1.0, 1.0, 0.0, 0.0};
        double sy[4] = {-1.0, -1.0, 1.0, 1.0};

#pragma omp parallel for num_threads(devNums)
        for (int devId = 0; devId < devNums; devId++) {

            auto& local = localF0[devId];
            const size_t picNums = picPitch / devNums;

            for (size_t picId = 0; picId < picNums; picId++) {

                if (hostId == 0)
                    if (devId == 0)
                        if (picId % (picNums / 4) == 0)
                            std::cout << BOLDGREEN << 100 * picId / picNums << "%" << RESET << std::endl;

                int i, j;
                double li, lj;
                double dx, dy;
                double J, B, N, T;
                double v, Jv;
                double x, y, vp, mu, pw;
                double vpara, vperp;
                double coes[4];

                if constexpr (spaceType == 0) {
                    do {
                        x = x0 + (x1 - x0) * xrand[devId]();
                        y = y0 + (y1 - y0) * yrand[devId]();
                        N = interp2d(tempN, x, y);
                        J = interp2d(tempJ, x, y);
                    } while (Jrand[devId]() >= N * J / Jmax);
                } else if constexpr (spaceType == 1) {
                    do {
                        x = x0 + (x1 - x0) * xrand[devId]();
                        y = y0 + (y1 - y0) * yrand[devId]();
                        J = interp2d(tempJ, x, y);
                    } while (Jrand[devId]() >= J / Jmax);
                }

                B = interp2d(tempB, x, y);
                N = interp2d(tempN, x, y);
                T = interp2d(tempT, x, y);

                if constexpr (velocityType == 0) {

                    if constexpr (disType == 0) {

                        double vth = std::sqrt(T * KEV / (Mass * MP * std::pow(VA0, 2.0)));

                        do {
                            vpara = vth * vparand[devId]();
                            vperp = vth * std::sqrt(-2.0 * std::log(1.0 - vperand[devId]()));
                            if constexpr (picType == 2)
                                vpara = std::fabs(vpara);
                            v = std::sqrt(std::pow(vpara, 2.0) + std::pow(vperp, 2.0));
                        } while (v <= Vmin || v >= Vmax);

                        vp = vpara;
                        mu = 0.5 * Mass * std::pow(vperp, 2.0) / B;

                    } else if constexpr (disType == 1) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);
                        double U = vrand[devId]();
                        v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);

                        double pitch;
                        if constexpr (picType == 2)
                            pitch = vprand[devId]();
                        else
                            pitch = 2.0 * vprand[devId]() - 1.0;

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * (1.0 - std::pow(pitch, 2.0)) / B;

                    } else if constexpr (disType == 2) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);

                        do {
                            double U = vrand[devId]();
                            v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);
                        } while (Jvrand[devId]() >= 0.5 * (1.0 + std::erf((Vb - v) / DeltaV)));

                        double pitch;
                        if constexpr (picType == 2)
                            pitch = vprand[devId]();
                        else
                            pitch = 2.0 * vprand[devId]() - 1.0;

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * (1.0 - std::pow(pitch, 2.0)) / B;

                    } else if constexpr (disType == 3) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);
                        double U = vrand[devId]();
                        v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);

                        double pitch, Lambda;
                        do {
                            if constexpr (picType == 2)
                                pitch = vprand[devId]();
                            else
                                pitch = 2.0 * vprand[devId]() - 1.0;
                            Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        } while (Jvrand[devId]() >= std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2));

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;

                    } else if constexpr (disType == 4) {

                        double Vmin3 = std::pow(Vmin, 3.0);
                        double Vmax3 = std::pow(Vmax, 3.0);
                        double T3 = std::pow(T, 3.0);

                        do {
                            double U = vrand[devId]();
                            v = std::pow((Vmin3 + T3) * std::pow((Vmax3 + T3) / (Vmin3 + T3), U) - T3, 1.0 / 3.0);
                        } while (Jvrand[devId]() >= 0.5 * (1.0 + std::erf((Vb - v) / DeltaV)));

                        double pitch, Lambda;
                        do {
                            if constexpr (picType == 2)
                                pitch = vprand[devId]();
                            else
                                pitch = 2.0 * vprand[devId]() - 1.0;
                            Lambda = (1.0 - std::pow(pitch, 2.0)) / B;
                        } while (vperand[devId]() >= std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2));

                        vp = v * pitch;
                        mu = 0.5 * Mass * std::pow(v, 2.0) * Lambda;
                    }

                    pw = N;

                } else if constexpr (velocityType == 1) {

                    do {
                        v = Vmax * vrand[devId]();
                        Jv = std::pow(v, 2.0);
                    } while (v <= Vmin || Jvrand[devId]() >= Jv / Jvmax);

                    double pitch;
                    if constexpr (picType == 2)
                        pitch = vprand[devId]();
                    else
                        pitch = 2.0 * vprand[devId]() - 1.0;
                    vp = v * pitch;
                    mu = 0.5 * Mass * std::pow(v, 2.0) * (1.0 - std::pow(pitch, 2.0)) / B;
                    double Lambda = (1.0 - std::pow(pitch, 2.0)) / B;

                    if constexpr (disType == 0)
                        pw = N * std::pow(T, -1.5) *
                             std::exp(-0.5 * Mass * std::pow(v, 2.0) * MP * std::pow(VA0, 2.0) / (T * KEV));
                    else if constexpr (disType == 1)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0));
                    else if constexpr (disType == 2)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) * (1.0 + std::erf((Vb - v) / DeltaV));
                    else if constexpr (disType == 3)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) *
                             std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2);
                    else if constexpr (disType == 4)
                        pw = N / (std::pow(v, 3.0) + std::pow(T, 3.0)) *
                             std::exp(-std::pow(Lambda - Lambda0, 2.0) / DeltaLambda2) *
                             (1.0 + std::erf((Vb - v) / DeltaV));
                }

                vpara = vp;
                vperp = std::sqrt(2.0 * mu * B / Mass);

                if constexpr (spaceType == 0)
                    pw /= N;
                else if constexpr (spaceType == 1)
                    pw *= 1.0;

                li = (vpara - minVpara) / dVpara;
                lj = (vperp - minVperp) / dVperp;

                i = std::floor(li);
                j = std::floor(lj);

                dx = li - i;
                dy = lj - j;

                if (i == gridVpara - 1) {
                    i--;
                    dx = 1;
                }
                if (j == gridVperp - 1) {
                    j--;
                    dy = 1;
                }

                for (int index = 0; index < 4; index++)
                    coes[index] = (cx[index] + sx[index] * dx) * (cy[index] + sy[index] * dy);

                const size_t strideV = gridVperp;
                const size_t index = (size_t)i * strideV + j;

                local[index] += pw * coes[0];
                local[index + strideV] += pw * coes[1];
                local[index + 1] += pw * coes[2];
                local[index + strideV + 1] += pw * coes[3];
            }
        }

        if (hostId == 0)
            std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;

        double* f0 = pitchSpaceF0[0];

#pragma omp parallel for num_threads(devNums)
        for (size_t index = 0; index < f0Len; index++) {
            double sum = 0.0;
            for (int devId = 0; devId < devNums; devId++)
                sum += localF0[devId][index];
            f0[index] = sum;
        }

        MPICHECK(
            MPI_Allreduce(MPI_IN_PLACE, pitchSpaceF0[0], gridVpara * gridVperp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

        for (int i = 0; i < gridVpara; i++) {
            for (int j = 0; j < gridVperp; j++) {

                bool isCorner = (i == 0 || i == gridVpara - 1) && (j == 0 || j == gridVperp - 1);

                bool isEdge = (i == 0 || i == gridVpara - 1 || j == 0 || j == gridVperp - 1);

                if (isCorner)
                    pitchSpaceF0[i][j] *= 4;
                else if (isEdge)
                    pitchSpaceF0[i][j] *= 2;
            }
        }

        if (hostId == 0) {

            std::ofstream output;
            std::string fileName;

            if constexpr (picType == 0)
                fileName = initialDir + "/IonPitchSpaceF0.bin";
            else if constexpr (picType == 1)
                fileName = initialDir + "/AlphaPitchSpaceF0.bin";
            else if constexpr (picType == 2)
                fileName = initialDir + "/BeamPitchSpaceF0.bin";

            output.open(fileName.c_str(), std::ios::out | std::ios::binary);

            output.write(reinterpret_cast<char*>(pitchSpaceF0[0]), sizeof(double) * gridVpara * gridVperp);

            output.close();

            logDone();
        }

        HostAllocator.releaseHostArrays(pitchSpaceF0);
    }
    template <int picType, int gridE, int gridPphi, int gridLambda, int ppcPhase>
    void computePhaseSpaceJacobian(std::string initialDir) {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Compute phase space jacobian of thermal ions.");
            else if constexpr (picType == 1)
                logStart("Compute phase space jacobian of alpha particles.");
            else if constexpr (picType == 2)
                logStart("Compute phase space jacobian of beam particles.");
        }

        double minE, maxE, dE;
        double minPphi, maxPphi, dPphi;
        double minLambda, maxLambda, dLambda;
        double Mass, Char, Vmin, Vmax;

        double drho = RHO1 - RHO0;
        double psitmax = PSITMAX / (B0 * L0 * L0);
        double cm = VA0 / (L0 * (QE * B0 / MP));

        if constexpr (picType == 0) {

            Mass = IonMass;
            Char = IonChar;
            Vmin = IonVmin;
            Vmax = IonVmax;

        } else if constexpr (picType == 1) {

            Mass = AlphaMass;
            Char = AlphaChar;
            Vmin = AlphaVmin;
            Vmax = AlphaVmax;

        } else if constexpr (picType == 2) {

            Mass = BeamMass;
            Char = BeamChar;
            Vmin = BeamVmin;
            Vmax = BeamVmax;
        }

        minE = 0.5 * Mass * std::pow(Vmin, 2.0);
        maxE = 0.5 * Mass * std::pow(Vmax, 2.0);

        minPphi = 20251106;
        maxPphi = -20251106;

        minLambda = 0.0;
        maxLambda = 0.0;

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                double tempPphi, tempLambda;

                tempLambda = 1.0 / B[i][j];
                if (tempLambda > maxLambda)
                    maxLambda = tempLambda;

                tempPphi = cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) *
                               SFAcovyz[i][j + gridGhost] / (q[i][j] * J[i][j] * B[i][j]) -
                           Char * psip[i][j];
                if (tempPphi > maxPphi)
                    maxPphi = tempPphi;

                tempPphi = -cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) *
                               SFAcovyz[i][j + gridGhost] / (q[i][j] * J[i][j] * B[i][j]) -
                           Char * psip[i][j];
                if (tempPphi < minPphi)
                    minPphi = tempPphi;
            }
        }

        dE = (maxE - minE) / (gridE - 1);
        dPphi = (maxPphi - minPphi) / (gridPphi - 1);
        dLambda = (maxLambda - minLambda) / (gridLambda - 1);

        size_t picPhase = (size_t)gridE * gridPphi * gridLambda * ppcPhase / hostNums;

        double*** phaseSpaceJacobian;
        Allocator HostAllocator;
        HostAllocator.allocateHostArrays(gridE, gridPphi, gridLambda, phaseSpaceJacobian);

        const size_t jacobianLen = (size_t)gridE * gridPphi * gridLambda;
        std::vector<std::vector<double>> localJacobian(devNums, std::vector<double>(jacobianLen, 0.0));

        double Jmax, Jvmax;
        Jmax = 0.0;
        Jvmax = std::pow(Vmax, 2.0);

        std::vector<std::vector<double>> tempq(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> temppsip(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempJ(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempB(gridNx, std::vector<double>(gridNy + 2));
        std::vector<std::vector<double>> tempSFAcovyz(gridNx, std::vector<double>(gridNy + 2));

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                tempq[i][j + 1] = q[i][j];
                temppsip[i][j + 1] = psip[i][j];
                tempJ[i][j + 1] = J[i][j];
                tempB[i][j + 1] = B[i][j];

                if (tempJ[i][j + 1] > Jmax)
                    Jmax = tempJ[i][j + 1];
            }

            tempq[i][0] = tempq[i][gridNy];
            tempq[i][gridNy + 1] = tempq[i][1];

            temppsip[i][0] = temppsip[i][gridNy];
            temppsip[i][gridNy + 1] = temppsip[i][1];

            tempJ[i][0] = tempJ[i][gridNy];
            tempJ[i][gridNy + 1] = tempJ[i][1];

            tempB[i][0] = tempB[i][gridNy];
            tempB[i][gridNy + 1] = tempB[i][1];
        }

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy + 2; j++) {
                tempSFAcovyz[i][j] = SFAcovyz[i][j + gridGhost - 1];
            }
        }

        auto interp2d = [&](std::vector<std::vector<double>>& field, double x, double y) {
            double li = (x - x0) / gridDx;
            double lj = (y - y0 + 0.5 * gridDy) / gridDy;
            int i = std::floor(li);
            int j = std::floor(lj);
            double dx = li - i;
            double dy = lj - j;

            double coes[4] = {};
            double cx[4] = {1.0, 1.0, 0.0, 0.0};
            double sx[4] = {-1.0, -1.0, 1.0, 1.0};
            double cy[4] = {1.0, 0.0, 1.0, 0.0};
            double sy[4] = {-1.0, 1.0, -1.0, 1.0};

            double result = 0.0;

            coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
            coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
            coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
            coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

            result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] +
                     field[i + 1][j + 1] * coes[3];

            return result;
        };

        std::vector<Rand01> xrand(devNums);
        std::vector<Rand01> yrand(devNums);
        std::vector<Rand01> Jrand(devNums);

        std::vector<Rand01> vrand(devNums);
        std::vector<Rand01> vprand(devNums);
        std::vector<Rand01> Jvrand(devNums);

        double cx[8] = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
        double sx[8] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
        double cy[8] = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
        double sy[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
        double cz[8] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        double sz[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

#pragma omp parallel for num_threads(devNums)
        for (int devId = 0; devId < devNums; devId++) {

            auto& local = localJacobian[devId];
            const size_t picNums = picPhase / devNums;

            for (size_t picId = 0; picId < picNums; picId++) {

                if (hostId == 0)
                    if (devId == 0)
                        if (picId % (picNums / 4) == 0)
                            std::cout << BOLDGREEN << 100 * picId / picNums << "%" << RESET << std::endl;

                int i, j, k;
                double li, lj, lk;
                double dx, dy, dz;
                double q, psip, J, B, SFAcovyz;
                double v, Jv, E, Pphi, Lambda;
                double x, y, vp, mu;
                double coes[8];

                do {
                    x = x0 + (x1 - x0) * xrand[devId]();
                    y = y0 + (y1 - y0) * yrand[devId]();
                    J = interp2d(tempJ, x, y);
                } while (Jrand[devId]() >= J / Jmax);

                do {
                    v = Vmax * vrand[devId]();
                    Jv = std::pow(v, 2.0);
                } while (v <= Vmin || Jvrand[devId]() >= Jv / Jvmax);

                q = interp2d(tempq, x, y);
                psip = interp2d(temppsip, x, y);
                SFAcovyz = interp2d(tempSFAcovyz, x, y);
                B = interp2d(tempB, x, y);

                if constexpr (picType == 2)
                    Lambda = vprand[devId]();
                else
                    Lambda = 2.0 * vprand[devId]() - 1.0;
                vp = v * Lambda;
                mu = 0.5 * Mass * std::pow(v, 2.0) * (1.0 - std::pow(Lambda, 2.0)) / B;

                E = 0.5 * Mass * std::pow(v, 2.0);
                Pphi = cm * Mass * vp * 2 * psitmax * drho * (RHO0 + x * drho) * SFAcovyz / (q * J * B) - Char * psip;
                Lambda = mu / E;

                li = (E - minE) / dE;
                lj = (Pphi - minPphi) / dPphi;
                lk = (Lambda - minLambda) / dLambda;

                i = std::floor(li);
                j = std::floor(lj);
                k = std::floor(lk);

                dx = li - i;
                dy = lj - j;
                dz = lk - k;

                if (i == gridE - 1) {
                    i--;
                    dx = 1;
                }
                if (j == gridPphi - 1) {
                    j--;
                    dy = 1;
                }
                if (k == gridLambda - 1) {
                    k--;
                    dz = 1;
                }

                for (int index = 0; index < 8; index++)
                    coes[index] =
                        (cx[index] + sx[index] * dx) * (cy[index] + sy[index] * dy) * (cz[index] + sz[index] * dz);

                const size_t strideE = (size_t)gridPphi * gridLambda;
                const size_t strideP = gridLambda;
                const size_t index = (size_t)i * strideE + (size_t)j * strideP + k;

                local[index] += coes[0];
                local[index + strideE] += coes[1];
                local[index + strideP] += coes[2];
                local[index + strideE + strideP] += coes[3];
                local[index + 1] += coes[4];
                local[index + strideE + 1] += coes[5];
                local[index + strideP + 1] += coes[6];
                local[index + strideE + strideP + 1] += coes[7];
            }
        }

        if (hostId == 0)
            std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;

        double* jacobian = phaseSpaceJacobian[0][0];

#pragma omp parallel for num_threads(devNums)
        for (size_t index = 0; index < jacobianLen; index++) {
            double sum = 0.0;
            for (int devId = 0; devId < devNums; devId++)
                sum += localJacobian[devId][index];
            jacobian[index] = sum;
        }

        MPICHECK(MPI_Allreduce(MPI_IN_PLACE, phaseSpaceJacobian[0][0], gridE * gridPphi * gridLambda, MPI_DOUBLE,
                               MPI_SUM, MPI_COMM_WORLD));

        for (int i = 0; i < gridE; i++) {
            for (int j = 0; j < gridPphi; j++) {
                for (int k = 0; k < gridLambda; k++) {

                    bool isVertex =
                        (i == 0 || i == gridE - 1) && (j == 0 || j == gridPphi - 1) && (k == 0 || k == gridLambda - 1);

                    bool isEdge = ((i == 0 || i == gridE - 1) && (j == 0 || j == gridPphi - 1)) ||
                                  ((i == 0 || i == gridE - 1) && (k == 0 || k == gridLambda - 1)) ||
                                  ((j == 0 || j == gridPphi - 1) && (k == 0 || k == gridLambda - 1));

                    bool isFace =
                        (i == 0 || i == gridE - 1 || j == 0 || j == gridPphi - 1 || k == 0 || k == gridLambda - 1);

                    if (isVertex)
                        phaseSpaceJacobian[i][j][k] *= 8;
                    else if (isEdge)
                        phaseSpaceJacobian[i][j][k] *= 4;
                    else if (isFace)
                        phaseSpaceJacobian[i][j][k] *= 2;
                }
            }
        }

        if (hostId == 0) {

            std::ofstream output;
            std::string fileName;

            if constexpr (picType == 0)
                fileName = initialDir + "/IonPhaseSpaceJacobian.bin";
            else if constexpr (picType == 1)
                fileName = initialDir + "/AlphaPhaseSpaceJacobian.bin";
            else if constexpr (picType == 2)
                fileName = initialDir + "/BeamPhaseSpaceJacobian.bin";

            output.open(fileName.c_str(), std::ios::out | std::ios::binary);

            output.write(reinterpret_cast<char*>(phaseSpaceJacobian[0][0]),
                         sizeof(double) * gridE * gridPphi * gridLambda);

            output.close();

            logDone();
        }

        HostAllocator.releaseHostArrays(phaseSpaceJacobian);
    }
    template <int picType, int gridVpara, int gridVperp, int ppcPitch>
    void computePitchSpaceJacobian(std::string initialDir) {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Compute pitch space jacobian of thermal ions.");
            else if constexpr (picType == 1)
                logStart("Compute pitch space jacobian of alpha particles.");
            else if constexpr (picType == 2)
                logStart("Compute pitch space jacobian of beam particles.");
        }

        double minVpara, maxVpara, dVpara;
        double minVperp, maxVperp, dVperp;
        double Vmin, Vmax;

        if constexpr (picType == 0) {

            Vmin = IonVmin;
            Vmax = IonVmax;

        } else if constexpr (picType == 1) {

            Vmin = AlphaVmin;
            Vmax = AlphaVmax;

        } else if constexpr (picType == 2) {

            Vmin = BeamVmin;
            Vmax = BeamVmax;
        }

        minVpara = (picType == 2) ? 0.0 : -Vmax;
        maxVpara = Vmax;
        minVperp = 0.0;
        maxVperp = Vmax;

        dVpara = (maxVpara - minVpara) / (gridVpara - 1);
        dVperp = (maxVperp - minVperp) / (gridVperp - 1);

        size_t picPitch = (size_t)gridVpara * gridVperp * ppcPitch / hostNums;

        double** pitchSpaceJacobian;
        Allocator HostAllocator;
        HostAllocator.allocateHostArrays(gridVpara, gridVperp, pitchSpaceJacobian);

        const size_t jacobianLen = (size_t)gridVpara * gridVperp;
        std::vector<std::vector<double>> localJacobian(devNums, std::vector<double>(jacobianLen, 0.0));

        double Jmax, Jvmax;
        Jmax = 0.0;
        Jvmax = std::pow(Vmax, 2.0);

        std::vector<std::vector<double>> tempJ(gridNx, std::vector<double>(gridNy + 2));

        for (int i = 0; i < gridNx; i++) {
            for (int j = 0; j < gridNy; j++) {

                tempJ[i][j + 1] = J[i][j];

                if (tempJ[i][j + 1] > Jmax)
                    Jmax = tempJ[i][j + 1];
            }

            tempJ[i][0] = tempJ[i][gridNy];
            tempJ[i][gridNy + 1] = tempJ[i][1];
        }

        auto interp2d = [&](std::vector<std::vector<double>>& field, double x, double y) {
            double li = (x - x0) / gridDx;
            double lj = (y - y0 + 0.5 * gridDy) / gridDy;
            int i = std::floor(li);
            int j = std::floor(lj);
            double dx = li - i;
            double dy = lj - j;

            double coes[4] = {};
            double cx[4] = {1.0, 1.0, 0.0, 0.0};
            double sx[4] = {-1.0, -1.0, 1.0, 1.0};
            double cy[4] = {1.0, 0.0, 1.0, 0.0};
            double sy[4] = {-1.0, 1.0, -1.0, 1.0};

            double result = 0.0;

            coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
            coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
            coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
            coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

            result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] +
                     field[i + 1][j + 1] * coes[3];

            return result;
        };

        std::vector<Rand01> xrand(devNums);
        std::vector<Rand01> yrand(devNums);
        std::vector<Rand01> Jrand(devNums);

        std::vector<Rand01> vrand(devNums);
        std::vector<Rand01> vprand(devNums);
        std::vector<Rand01> Jvrand(devNums);

        double cx[4] = {1.0, 0.0, 1.0, 0.0};
        double sx[4] = {-1.0, 1.0, -1.0, 1.0};
        double cy[4] = {1.0, 1.0, 0.0, 0.0};
        double sy[4] = {-1.0, -1.0, 1.0, 1.0};

#pragma omp parallel for num_threads(devNums)
        for (int devId = 0; devId < devNums; devId++) {

            auto& local = localJacobian[devId];
            const size_t picNums = picPitch / devNums;

            for (size_t picId = 0; picId < picNums; picId++) {

                if (hostId == 0)
                    if (devId == 0)
                        if (picId % (picNums / 4) == 0)
                            std::cout << BOLDGREEN << 100 * picId / picNums << "%" << RESET << std::endl;

                int i, j;
                double li, lj;
                double dx, dy;
                double J;
                double v, Jv;
                double x, y;
                double pitch, vpara, vperp;
                double coes[4];

                do {
                    x = x0 + (x1 - x0) * xrand[devId]();
                    y = y0 + (y1 - y0) * yrand[devId]();
                    J = interp2d(tempJ, x, y);
                } while (Jrand[devId]() >= J / Jmax);

                do {
                    v = Vmax * vrand[devId]();
                    Jv = std::pow(v, 2.0);
                } while (v <= Vmin || Jvrand[devId]() >= Jv / Jvmax);

                if constexpr (picType == 2)
                    pitch = vprand[devId]();
                else
                    pitch = 2.0 * vprand[devId]() - 1.0;

                vpara = v * pitch;
                vperp = v * std::sqrt(1.0 - std::pow(pitch, 2.0));

                li = (vpara - minVpara) / dVpara;
                lj = (vperp - minVperp) / dVperp;

                i = std::floor(li);
                j = std::floor(lj);

                dx = li - i;
                dy = lj - j;

                if (i == gridVpara - 1) {
                    i--;
                    dx = 1;
                }
                if (j == gridVperp - 1) {
                    j--;
                    dy = 1;
                }

                for (int index = 0; index < 4; index++)
                    coes[index] = (cx[index] + sx[index] * dx) * (cy[index] + sy[index] * dy);

                const size_t strideV = gridVperp;
                const size_t index = (size_t)i * strideV + j;

                local[index] += coes[0];
                local[index + strideV] += coes[1];
                local[index + 1] += coes[2];
                local[index + strideV + 1] += coes[3];
            }
        }

        if (hostId == 0)
            std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;

        double* jacobian = pitchSpaceJacobian[0];

#pragma omp parallel for num_threads(devNums)
        for (size_t index = 0; index < jacobianLen; index++) {
            double sum = 0.0;
            for (int devId = 0; devId < devNums; devId++)
                sum += localJacobian[devId][index];
            jacobian[index] = sum;
        }

        MPICHECK(MPI_Allreduce(MPI_IN_PLACE, pitchSpaceJacobian[0], gridVpara * gridVperp, MPI_DOUBLE, MPI_SUM,
                               MPI_COMM_WORLD));

        for (int i = 0; i < gridVpara; i++) {
            for (int j = 0; j < gridVperp; j++) {

                bool isCorner = (i == 0 || i == gridVpara - 1) && (j == 0 || j == gridVperp - 1);

                bool isEdge = (i == 0 || i == gridVpara - 1 || j == 0 || j == gridVperp - 1);

                if (isCorner)
                    pitchSpaceJacobian[i][j] *= 4;
                else if (isEdge)
                    pitchSpaceJacobian[i][j] *= 2;
            }
        }

        if (hostId == 0) {

            std::ofstream output;
            std::string fileName;

            if constexpr (picType == 0)
                fileName = initialDir + "/IonPitchSpaceJacobian.bin";
            else if constexpr (picType == 1)
                fileName = initialDir + "/AlphaPitchSpaceJacobian.bin";
            else if constexpr (picType == 2)
                fileName = initialDir + "/BeamPitchSpaceJacobian.bin";

            output.open(fileName.c_str(), std::ios::out | std::ios::binary);

            output.write(reinterpret_cast<char*>(pitchSpaceJacobian[0]), sizeof(double) * gridVpara * gridVperp);

            output.close();

            logDone();
        }

        HostAllocator.releaseHostArrays(pitchSpaceJacobian);
    }
    template <int picType>
    void loadPhaseSpaceMapping(std::string file) {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Load phase space mapping of thermal ions.");
            else if constexpr (picType == 1)
                logStart("Load phase space mapping of alpha particles.");
            else if constexpr (picType == 2)
                logStart("Load phase space mapping of beam particles.");
        }

        std::ifstream input;
        input.open(file, std::ios::in | std::ios::binary);

        input.seekg(0, std::ios::end);
        std::streamsize size = input.tellg();
        input.seekg(0, std::ios::beg);
        int picPhase = size / 28;

        std::vector<int> ids(picPhase);
        std::vector<double> rhos(picPhase);
        std::vector<double> vparas(picPhase);
        std::vector<double> mus(picPhase);

        input.read(reinterpret_cast<char*>(ids.data()), picPhase * sizeof(int));
        input.read(reinterpret_cast<char*>(rhos.data()), picPhase * sizeof(double));
        input.read(reinterpret_cast<char*>(vparas.data()), picPhase * sizeof(double));
        input.read(reinterpret_cast<char*>(mus.data()), picPhase * sizeof(double));
        input.close();

        Allocator HostDeviceAllocator;
        picReal** h_picPhaseMap;
        picReal** d_picPhaseMap;

        if constexpr (picType == 0) {
            HostDeviceAllocator.allocateHostArrays(devNums, (size_t)picPhase * 13, h_IonPhaseSpaceMapping);
            HostDeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picPhase * 13, d_IonPhaseSpaceMapping);
            h_picPhaseMap = h_IonPhaseSpaceMapping;
            d_picPhaseMap = d_IonPhaseSpaceMapping;
        } else if constexpr (picType == 1) {
            HostDeviceAllocator.allocateHostArrays(devNums, (size_t)picPhase * 13, h_AlphaPhaseSpaceMapping);
            HostDeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picPhase * 13, d_AlphaPhaseSpaceMapping);
            h_picPhaseMap = h_AlphaPhaseSpaceMapping;
            d_picPhaseMap = d_AlphaPhaseSpaceMapping;
        } else if constexpr (picType == 2) {
            HostDeviceAllocator.allocateHostArrays(devNums, (size_t)picPhase * 13, h_BeamPhaseSpaceMapping);
            HostDeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picPhase * 13, d_BeamPhaseSpaceMapping);
            h_picPhaseMap = h_BeamPhaseSpaceMapping;
            d_picPhaseMap = d_BeamPhaseSpaceMapping;
        }

        for (int i = 0; i < devNums; i++) {
            for (int picId = 0; picId < picPhase; picId++) {

                // orbit = 0.5 : pad
                // orbit = 1.5 : loss
                // orbit = 2.5 : para
                // orbit = 3.5 : anti
                // orbit = 4.5 : trapped
                // orbit = 5.5 : unknown

                // orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

                if (ids[picId] == 20251106)
                    h_picPhaseMap[i][picId * 13 + 0] = 0.5;
                else
                    h_picPhaseMap[i][picId * 13 + 0] = 5.5;

                h_picPhaseMap[i][picId * 13 + 1] = rhos[picId];
                h_picPhaseMap[i][picId * 13 + 3] = vparas[picId];
                h_picPhaseMap[i][picId * 13 + 4] = mus[picId];
                h_picPhaseMap[i][picId * 13 + 9] = 0.5;
            }
        }

        for (int i = 0; i < devNums; i++) {

            CUDACHECK(cudaSetDevice(localId * devNums + i));
            HostDeviceAllocator.hostToDevice((size_t)picPhase * 13, 0, (size_t)i * picPhase * 13, d_picPhaseMap[i],
                                             h_picPhaseMap[0]);
        }

        if (hostId == 0) {
            logDone();
        }
    }
    template <int picType>
    void computePhaseSpaceFrequency(std::string file, std::string initialDir) {

        if (hostId == 0) {

            if constexpr (picType == 0)
                logStart("Compute phase space orbit frequency of thermal ions.");
            else if constexpr (picType == 1)
                logStart("Compute phase space orbit frequency of alpha particles.");
            else if constexpr (picType == 2)
                logStart("Compute phase space orbit frequency of beam particles.");
        }

        std::ifstream input;
        input.open(file, std::ios::in | std::ios::binary);

        input.seekg(0, std::ios::end);
        std::streamsize size = input.tellg();
        input.seekg(0, std::ios::beg);
        int picPhase = size / 28;

        std::vector<int> ids(picPhase);
        input.read(reinterpret_cast<char*>(ids.data()), picPhase * sizeof(int));
        input.close();

        std::vector<double> phaseSpaceFrequency(picPhase * 9);

        picReal** h_picPhaseMap;
        if constexpr (picType == 0)
            h_picPhaseMap = h_IonPhaseSpaceMapping;
        else if constexpr (picType == 1)
            h_picPhaseMap = h_AlphaPhaseSpaceMapping;
        else if constexpr (picType == 2)
            h_picPhaseMap = h_BeamPhaseSpaceMapping;

        for (int picId = 0; picId < picPhase; picId++) {

            // orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

            phaseSpaceFrequency[picId * 9 + 0] = ids[picId];
            phaseSpaceFrequency[picId * 9 + 1] = h_picPhaseMap[0][picId * 13 + 0];
            phaseSpaceFrequency[picId * 9 + 2] = h_picPhaseMap[0][picId * 13 + 5];
            phaseSpaceFrequency[picId * 9 + 3] = h_picPhaseMap[0][picId * 13 + 6];
            phaseSpaceFrequency[picId * 9 + 4] = h_picPhaseMap[0][picId * 13 + 7];
            phaseSpaceFrequency[picId * 9 + 5] = h_picPhaseMap[0][picId * 13 + 8];
            phaseSpaceFrequency[picId * 9 + 6] = h_picPhaseMap[0][picId * 13 + 10];
            phaseSpaceFrequency[picId * 9 + 7] = h_picPhaseMap[0][picId * 13 + 11];
            phaseSpaceFrequency[picId * 9 + 8] = h_picPhaseMap[0][picId * 13 + 12];
        }

        Allocator HostDeviceAllocator;
        if constexpr (picType == 0) {
            HostDeviceAllocator.releaseHostArrays(h_IonPhaseSpaceMapping);
            HostDeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonPhaseSpaceMapping);
        } else if constexpr (picType == 1) {
            HostDeviceAllocator.releaseHostArrays(h_AlphaPhaseSpaceMapping);
            HostDeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaPhaseSpaceMapping);
        } else if constexpr (picType == 2) {
            HostDeviceAllocator.releaseHostArrays(h_BeamPhaseSpaceMapping);
            HostDeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamPhaseSpaceMapping);
        }

        if (hostId == 0) {

            std::ofstream output;
            std::string fileName;

            if constexpr (picType == 0)
                fileName = initialDir + "/IonPhaseSpaceFrequency.bin";
            else if constexpr (picType == 1)
                fileName = initialDir + "/AlphaPhaseSpaceFrequency.bin";
            else if constexpr (picType == 2)
                fileName = initialDir + "/BeamPhaseSpaceFrequency.bin";

            output.open(fileName.c_str(), std::ios::out | std::ios::binary);
            output.write(reinterpret_cast<char*>(phaseSpaceFrequency.data()), sizeof(double) * picPhase * 9);
            output.close();

            logDone();
        }
    }
};
