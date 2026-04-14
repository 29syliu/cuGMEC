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
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "mpi.h"
#include "nccl.h"
#include "unistd.h"
#include "cufft.h"
#include "cudss.h"
#include "cublas_v2.h"
#include "cub/cub.cuh"
#include "cusparse_v2.h"
#include "curand_kernel.h"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>
#include <math.h>
#include <random>
#include <string>
#include <iomanip>
#include <omp.h>

#define RESET                       "\033[0m"
#define RED                           "\033[31m"
#define BLUE                         "\033[34m"
#define CYAN                        "\033[36m"
#define BLACK                       "\033[30m"
#define WHITE                      "\033[37m"
#define GREEN                      "\033[32m"
#define YELLOW                   "\033[33m"
#define MAGENTA                "\033[35m"
#define BOLDRED           "\033[1m\033[31m"
#define BOLDBLUE          "\033[1m\033[34m"
#define BOLDCYAN        "\033[1m\033[36m"
#define BOLDBLACK       "\033[1m\033[30m"
#define BOLDWHITE       "\033[1m\033[37m"
#define BOLDGREEN      "\033[1m\033[32m"
#define BOLDYELLOW    "\033[1m\033[33m"
#define BOLDMAGENTA "\033[1m\033[35m"

#define MPICHECK(cmd) do {                            \
  int e = cmd;                                                       \
  if( e != MPI_SUCCESS ) {                                  \
    printf("Failed: MPI error %s:%d '%d'\n",              \
        __FILE__,__LINE__, e);                                \
    exit(EXIT_FAILURE);                                      \
  }                                                                        \
} while(0)

#define NCCLCHECK(cmd) do {                         \
  ncclResult_t res = cmd;				            \
  if (res != ncclSuccess) {                                         \
    printf("Failed, NCCL error %s:%d '%s'\n",           \
        __FILE__,__LINE__,ncclGetErrorString(res));  \
    exit(EXIT_FAILURE);                                      \
  }                                                                        \
} while(0)

#define CUDACHECK(cmd) do {                         \
  cudaError_t err = cmd;                                         \
  if (err != cudaSuccess) {                                        \
    printf("Failed: CUDA error %s:%d '%s'\n",           \
        __FILE__,__LINE__,cudaGetErrorString(err)); \
    exit(EXIT_FAILURE);                                      \
  }                                                                        \
} while(0)

#define CUDSSCHECK(cmd) do {                       \
  cudssStatus_t err = cmd;                                      \
  if (err != CUDSS_STATUS_SUCCESS) {            \
    printf("Failed: CUDSS error %s:%d '%d'\n",        \
        __FILE__,__LINE__,err);                              \
    exit(EXIT_FAILURE);                                     \
  }                                                                       \
} while(0)

#define CUSPARSECHECK(cmd) do {                 \
  cusparseStatus_t err = cmd;                                   \
  if (err != CUSPARSE_STATUS_SUCCESS) {      \
    printf("Failed: CUSPARSE error %s:%d '%d'\n",  \
        __FILE__,__LINE__,err);                               \
    exit(EXIT_FAILURE);                                      \
  }                                                                        \
} while(0)

#define CUFFTCHECK(cmd) do {                       \
  cufftResult_t err = cmd;                                      \
  if (err != CUFFT_SUCCESS) {                           \
    printf("Failed: CUFFT error %s:%d '%d'\n",        \
        __FILE__,__LINE__,err);                              \
    exit(EXIT_FAILURE);                                     \
  }                                                                       \
} while(0)

template<typename dataType>
class HybridModel {

public:

	HybridModel(
		std::vector<int> scaleInput,
		std::vector<double> normalizationInput,
		std::tuple<bool, double> perp2AInput,
		std::tuple<bool, double> perp2wInput,
		std::tuple<bool, double> perp2dNeInput,
		std::tuple<bool, double> perp2dTeInput,
		std::tuple<bool, double> perp2dPiInput,
		std::tuple<bool, double> perp2dPaInput,
		std::tuple<bool, double> perp2dPbInput,
		std::tuple<bool, unsigned int> SlowingInput,
		std::tuple<bool, std::vector<double>> IonInput,
		std::tuple<bool, std::vector<double>> AlphaInput,
		std::tuple<bool, std::vector<double>> BeamInput) :

		devNums{ scaleInput[0] },
		hostNums{ scaleInput[1] }, hostId{ scaleInput[2] }, localId{ scaleInput[3] },
		gridNx{ scaleInput[4] }, gridNy{ scaleInput[5] }, gridNz{ scaleInput[6] },
		gridGhost{ scaleInput[7] }, ppcNums{ scaleInput[8] }, tubes{ scaleInput[9] },

		hostNy{ gridNy / hostNums }, devNy{ hostNy / devNums }, gridNxz{ gridNx * gridNz },
		gridNxPlusGhost{ gridNx }, gridNyPlusGhost{ gridNy + 2 * gridGhost }, gridNzPlusGhost{ gridNz + 2 * gridGhost },
		cellNx{ gridNxPlusGhost - 1 }, cellNy{ gridNyPlusGhost - 1 }, cellNz{ gridNzPlusGhost - 1 }, cellNxz{ cellNx * cellNz },
		picHost{ gridNx * gridNy * gridNz / hostNums * ppcNums }, picDev{ picHost / devNums },
		gridDx{ 1.0 / (gridNx - 1) }, gridDy{ 2.0 * PI / gridNy }, gridDz{ 2.0 * PI / tubes / gridNz },
		x0{ 0.0 }, x1{ 1.0 }, x0PlusGhost{ x0 }, x1PlusGhost{ x1 },
		y0{ -PI }, y1{ PI }, y0PlusGhost{ y0 - (gridGhost - 0.5) * gridDy }, y1PlusGhost{ y1 + (gridGhost - 0.5) * gridDy },
		z0{ -PI / tubes }, z1{ PI / tubes }, z0PlusGhost{ z0 - (gridGhost - 0.5) * gridDz }, z1PlusGhost{ z1 + (gridGhost - 0.5) * gridDz },

		B0{ normalizationInput[0] }, L0{ normalizationInput[1] }, VA0{ normalizationInput[2] }, RHO0{ normalizationInput[3] }, RHO1{ normalizationInput[4] },
		PSITMAX{ normalizationInput[5] }, NormQE{ QE / (B0 * L0 * L0 / MU0 / VA0) }, dt{ normalizationInput[6] },

		nablaPerp2A{ perp2AInput }, nablaPerp2w{ perp2wInput }, nablaPerp2dNe{ perp2dNeInput }, nablaPerp2dTe{ perp2dTeInput },
		nablaPerp2dPi{ perp2dPiInput }, nablaPerp2dPa{ perp2dPaInput }, nablaPerp2dPb{ perp2dPbInput },

		ifSlowing{ std::get<0>(SlowingInput) }, randMax{ std::get<1>(SlowingInput) },

		ifIon{ std::get<0>(IonInput) }, IonMass{ std::get<1>(IonInput)[0] }, IonChar{ std::get<1>(IonInput)[1] }, IonBeta{ std::get<1>(IonInput)[2] }, IonVmin{ std::get<1>(IonInput)[3] }, IonVmax{ std::get<1>(IonInput)[4] },
		IonVb{ std::get<1>(IonInput)[5] }, IonDeltaV{ std::get<1>(IonInput)[6] }, IonLambda0{ std::get<1>(IonInput)[7] }, IonDeltaLambda2{ std::get<1>(IonInput)[8] },

		ifAlpha{ std::get<0>(AlphaInput) }, AlphaMass{ std::get<1>(AlphaInput)[0] }, AlphaChar{ std::get<1>(AlphaInput)[1] }, AlphaBeta{ std::get<1>(AlphaInput)[2] }, AlphaVmin{ std::get<1>(AlphaInput)[3] }, AlphaVmax{ std::get<1>(AlphaInput)[4] },
		AlphaVb{ std::get<1>(AlphaInput)[5] }, AlphaDeltaV{ std::get<1>(AlphaInput)[6] }, AlphaLambda0{ std::get<1>(AlphaInput)[7] }, AlphaDeltaLambda2{ std::get<1>(AlphaInput)[8] },

		ifBeam{ std::get<0>(BeamInput) }, BeamMass{ std::get<1>(BeamInput)[0] }, BeamChar{ std::get<1>(BeamInput)[1] }, BeamBeta{ std::get<1>(BeamInput)[2] }, BeamVmin{ std::get<1>(BeamInput)[3] }, BeamVmax{ std::get<1>(BeamInput)[4] },
		BeamVb{ std::get<1>(BeamInput)[5] }, BeamDeltaV{ std::get<1>(BeamInput)[6] }, BeamLambda0{ std::get<1>(BeamInput)[7] }, BeamDeltaLambda2{ std::get<1>(BeamInput)[8] } {

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
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;

		}

	}

	~HybridModel() {}

	/*-------------------------------------------------------Model Parameters-------------------------------------------------------*/

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

	const std::tuple<bool, double> nablaPerp2A;
	const std::tuple<bool, double> nablaPerp2w;
	const std::tuple<bool, double> nablaPerp2dNe;
	const std::tuple<bool, double> nablaPerp2dTe;
	const std::tuple<bool, double> nablaPerp2dPi;
	const std::tuple<bool, double> nablaPerp2dPa;
	const std::tuple<bool, double> nablaPerp2dPb;

	const bool ifSlowing;
	const unsigned int randMax;

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

	/*--------------------------------------------------MHD Equilibrium on CPU--------------------------------------------------*/

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

	/*--------------------------------------------------MHD Perturbation on CPU--------------------------------------------------*/

	dataType** h_qtheta;

	dataType*** h_w;
	dataType*** h_A;
	dataType*** h_dNe;
	dataType*** h_dTe;
	dataType*** h_Phi;
	dataType*** h_dJpB;
	dataType*** h_dPe;

	/*-----------------------------------------------Coefficient Compression on CPU-----------------------------------------------*/

	/*------------------------------------Linear------------------------------------*/

	dataType** h_A_w; dataType** h_A_px_w; dataType** h_A_py_w; dataType** h_A_pz_w;
	dataType** h_dJpB_w; dataType** h_dJpB_px_w; dataType** h_dJpB_py_w; dataType** h_dJpB_pz_w;
	dataType** h_dP_w; dataType** h_dP_px_w; dataType** h_dP_py_w; dataType** h_dP_pz_w;
	dataType** h_w_py_w; dataType** h_w_pz_w; dataType** h_w_Phi;

	double** h_Phi_w;
	double** h_Phi_px_w; double** h_Phi_pz_w;
	double** h_Phi_px2_w;  double** h_Phi_pxz_w; double** h_Phi_pz2_w;

	double** h_A_resistive;
	double** h_A_px_resistive; double** h_A_pz_resistive;
	double** h_A_px2_resistive; double** h_A_pxz_resistive; double** h_A_pz2_resistive;

	double** h_F_perp2;
	double** h_F_px_perp2; double** h_F_pz_perp2;
	double** h_F_px2_perp2; double** h_F_pxz_perp2; double** h_F_pz2_perp2;

	dataType** h_A_dJpB;
	dataType** h_A_px_dJpB;  dataType** h_A_pz_dJpB;
	dataType** h_A_px2_dJpB;  dataType** h_A_pxz_dJpB; dataType** h_A_pz2_dJpB;

	dataType** h_Phi_A; dataType** h_Phi_px_A; dataType** h_Phi_py_A; dataType** h_Phi_pz_A;
	dataType** h_dNe_A; dataType** h_dNe_px_A; dataType** h_dNe_py_A; dataType** h_dNe_pz_A;
	dataType** h_A_A; dataType** h_A_px_A; dataType** h_A_py_A; dataType** h_A_pz_A;

	dataType** h_Phi_dNe; dataType** h_Phi_px_dNe; dataType** h_Phi_py_dNe; dataType** h_Phi_pz_dNe;
	dataType** h_dPe_dNe; dataType** h_dPe_px_dNe; dataType** h_dPe_py_dNe; dataType** h_dPe_pz_dNe;
	dataType** h_dJpB_dNe; dataType** h_dJpB_px_dNe; dataType** h_dJpB_py_dNe; dataType** h_dJpB_pz_dNe;
	dataType** h_A_dNe; dataType** h_A_px_dNe; dataType** h_A_py_dNe; dataType** h_A_pz_dNe;

	dataType** h_Phi_dTe; dataType** h_Phi_px_dTe; dataType** h_Phi_py_dTe; dataType** h_Phi_pz_dTe;
	dataType** h_dTe_dTe; dataType** h_dTe_px_dTe; dataType** h_dTe_py_dTe; dataType** h_dTe_pz_dTe;
	dataType** h_dNe_dTe; dataType** h_dNe_px_dTe; dataType** h_dNe_py_dTe; dataType** h_dNe_pz_dTe;

	dataType** h_Ne0; dataType** h_Te0;
	dataType** h_Ne0_px; dataType** h_Te0_px; dataType** h_Pe0_px;

	dataType*** h_F2perp2;
	dataType*** h_A2dJpB;
	dataType*** h_Phi2w;
	dataType*** h_wdPAdJpB2w;
	dataType*** h_APhidNe2A;
	dataType*** h_dPePhiAdJpB2dNe;
	dataType*** h_PhidTedNe2dTe;

	/*----------------------------------Nonlinear----------------------------------*/

	dataType*** h_wPhi_w; dataType*** h_AdJpB_w;
	dataType*** h_PhiA_A; dataType*** h_NeA_A;
	dataType*** h_AdJpB_dNe; dataType*** h_dNePhi_dNe;
	dataType*** h_PhiTe_dTe; dataType*** h_PhiTeA_dTe;

	/*--------------------------------------------------------Matrix on CPU--------------------------------------------------------*/

	std::vector<int> matrix_i;
	std::vector<int> matrix_j;
	std::vector<dataType> matrix_v;

	/*--------------------------------------------------MHD Perturbation on GPU--------------------------------------------------*/

	dataType** d_qtheta;

	dataType** d_w;
	dataType** d_A;
	dataType** d_dNe;
	dataType** d_dTe;
	dataType** d_Phi;
	dataType** d_dJpB;
	dataType** d_dPe;

	dataType** d_w_beg; dataType** d_w_midl; dataType** d_w_midr; dataType** d_w_end;
	dataType** d_A_beg; dataType** d_A_midl; dataType** d_A_midr; dataType** d_A_end;
	dataType** d_dNe_beg; dataType** d_dNe_midl; dataType** d_dNe_midr; dataType** d_dNe_end;
	dataType** d_dTe_beg; dataType** d_dTe_midl; dataType** d_dTe_midr; dataType** d_dTe_end;
	dataType** d_Phi_midl; dataType** d_Phi_midr;
	dataType** d_dJpB_midl; dataType** d_dJpB_midr;
	dataType** d_dPe_midl; dataType** d_dPe_midr;
	dataType** d_dPi_midl; dataType** d_dPi_midr;
	dataType** d_dPa_midl; dataType** d_dPa_midr;
	dataType** d_dPb_midl; dataType** d_dPb_midr;
	dataType** d_Apt_midl; dataType** d_Apt_midr;

	/*-----------------------------------------------Coefficient Compression on GPU-----------------------------------------------*/

	/*------------------------------------Linear------------------------------------*/

	dataType** d_A_w; dataType** d_A_px_w; dataType** d_A_py_w; dataType** d_A_pz_w;
	dataType** d_dJpB_w; dataType** d_dJpB_px_w; dataType** d_dJpB_py_w; dataType** d_dJpB_pz_w;
	dataType** d_dP_w; dataType** d_dP_px_w; dataType** d_dP_py_w; dataType** d_dP_pz_w;
	dataType** d_w_py_w; dataType** d_w_pz_w; dataType** d_w_Phi;

	dataType** d_F_perp2;
	dataType** d_F_px_perp2; dataType** d_F_pz_perp2;
	dataType** d_F_px2_perp2; dataType** d_F_pxz_perp2; dataType** d_F_pz2_perp2;

	dataType** d_A_dJpB;
	dataType** d_A_px_dJpB; dataType** d_A_pz_dJpB;
	dataType** d_A_px2_dJpB;  dataType** d_A_pxz_dJpB; dataType** d_A_pz2_dJpB;

	dataType** d_Phi_A; dataType** d_Phi_px_A; dataType** d_Phi_py_A; dataType** d_Phi_pz_A;
	dataType** d_dNe_A; dataType** d_dNe_px_A; dataType** d_dNe_py_A; dataType** d_dNe_pz_A;
	dataType** d_A_A; dataType** d_A_px_A; dataType** d_A_py_A; dataType** d_A_pz_A;

	dataType** d_Phi_dNe; dataType** d_Phi_px_dNe; dataType** d_Phi_py_dNe; dataType** d_Phi_pz_dNe;
	dataType** d_dPe_dNe; dataType** d_dPe_px_dNe; dataType** d_dPe_py_dNe; dataType** d_dPe_pz_dNe;
	dataType** d_dJpB_dNe; dataType** d_dJpB_px_dNe; dataType** d_dJpB_py_dNe; dataType** d_dJpB_pz_dNe;
	dataType** d_A_dNe; dataType** d_A_px_dNe; dataType** d_A_py_dNe; dataType** d_A_pz_dNe;

	dataType** d_Phi_dTe; dataType** d_Phi_px_dTe; dataType** d_Phi_py_dTe; dataType** d_Phi_pz_dTe;
	dataType** d_dTe_dTe; dataType** d_dTe_px_dTe; dataType** d_dTe_py_dTe; dataType** d_dTe_pz_dTe;
	dataType** d_dNe_dTe; dataType** d_dNe_px_dTe; dataType** d_dNe_py_dTe; dataType** d_dNe_pz_dTe;

	dataType** d_Ne0; dataType** d_Te0;
	dataType** d_Ne0_px; dataType** d_Te0_px; dataType** d_Pe0_px;

	dataType** d_F2perp2;
	dataType** d_A2dJpB;
	dataType** d_Phi2w;
	dataType** d_wdPAdJpB2w;
	dataType** d_APhidNe2A;
	dataType** d_dPePhiAdJpB2dNe;
	dataType** d_PhidTedNe2dTe;

	/*----------------------------------Nonlinear----------------------------------*/

	dataType** d_wPhi_w; dataType** d_AdJpB_w;
	dataType** d_PhiA_A; dataType** d_NeA_A;
	dataType** d_AdJpB_dNe; dataType** d_dNePhi_dNe;
	dataType** d_PhiTe_dTe; dataType** d_PhiTeA_dTe;

	/*--------------------------------------------------------Matrix on GPU--------------------------------------------------------*/

	int** d_laplacianCsrR;
	int** d_laplacianCsrC;
	dataType** d_laplacianCsrV;

	int** d_resistiveCsrR;
	int** d_resistiveCsrC;
	dataType** d_resistiveCsrV;

	int** d_wCsrR;
	int** d_wCsrC;
	dataType** d_wCsrV;

	int** d_dNeCsrR;
	int** d_dNeCsrC;
	dataType** d_dNeCsrV;

	int** d_dTeCsrR;
	int** d_dTeCsrC;
	dataType** d_dTeCsrV;

	int** d_dPiCsrR;
	int** d_dPiCsrC;
	dataType** d_dPiCsrV;

	int** d_dPaCsrR;
	int** d_dPaCsrC;
	dataType** d_dPaCsrV;

	int** d_dPbCsrR;
	int** d_dPbCsrC;
	dataType** d_dPbCsrV;

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

	std::vector<std::vector<cudssConfig_t>> wConfigs;
	std::vector<std::vector<cudssData_t>> wDatas;
	std::vector<std::vector<cudssMatrix_t>> wAs;
	std::vector<std::vector<cudssMatrix_t>> wXs;
	std::vector<std::vector<cudssMatrix_t>> wBs;

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

	/*---------------------------------------------------------Particle in Cell---------------------------------------------------------*/

	dataType** dpic_Phi_py_A; dataType** dpic_dNe_py_A;
	dataType** dpic_A_A;  dataType** dpic_A_py_A; dataType** dpic_A_pz_A;
	dataType** dpic_APhidNe2A; dataType** dpic_PhiA_A; dataType** dpic_NeA_A;

	dataType IonConst;
	dataType AlphaConst;
	dataType BeamConst;

	dataType*** h_globalPi;
	dataType*** h_globalPa;
	dataType*** h_globalPb;

	dataType** d_globalA; dataType** d_globalPhi; dataType** d_globalApt;
	dataType** d_globalPi; dataType** d_globalPa; dataType** d_globalPb;

	dataType** h_pic1d; dataType** h_pic2d; dataType** h_pic3d;
	dataType** d_pic1d; dataType** d_pic2d; dataType** d_pic3d;

	int** h_Ion_offsets;
	int** h_Ion_keys;
	dataType** h_Ion_values;

	int** h_Alpha_offsets;
	int** h_Alpha_keys;
	dataType** h_Alpha_values;

	int** h_Beam_offsets;
	int** h_Beam_keys;
	dataType** h_Beam_values;

	int** d_Ion_offsets;
	int** d_Ion_keys_in;
	int** d_Ion_keys_out;
	dataType** d_Ion_values_in;
	dataType** d_Ion_values_out;

	int** d_Alpha_offsets;
	int** d_Alpha_keys_in;
	int** d_Alpha_keys_out;
	dataType** d_Alpha_values_in;
	dataType** d_Alpha_values_out;

	int** d_Beam_offsets;
	int** d_Beam_keys_in;
	int** d_Beam_keys_out;
	dataType** d_Beam_values_in;
	dataType** d_Beam_values_out;

	int** h_rand_keys;
	dataType** h_rand_values;
	int** d_rand_keys;
	dataType** d_rand_values;

	dataType** h_IonPhaseSpaceMapping;
	dataType** h_AlphaPhaseSpaceMapping;
	dataType** h_BeamPhaseSpaceMapping;

	dataType** d_IonPhaseSpaceMapping;
	dataType** d_AlphaPhaseSpaceMapping;
	dataType** d_BeamPhaseSpaceMapping;

	/*------------------------------------------------------------Function------------------------------------------------------------*/

	class Rand01 {

	public:

		Rand01() : mt_gen{ std::random_device()() }, rnd_dist{ 0.0,1.0 } {}
		~Rand01() {}
		double operator() () { return rnd_dist(mt_gen); }

		std::mt19937_64 mt_gen;
		std::uniform_real_distribution<double> rnd_dist;

	};
	class Allocator {

	public:

		Allocator() {}
		~Allocator() {}

		/*----------------Allocate and Release Arrays on Host and Device------------------*/

		template<typename T, typename... Ts>
		void allocateHostArrays(size_t dim1, size_t dim2, T**& hostArray, Ts**&... hostArrays) {

			hostArray = new T * [dim1]();
			hostArray[0] = new T[dim1 * dim2]();

			for (int i = 1; i < dim1; i++)
				hostArray[i] = hostArray[i - 1] + dim2;

			if constexpr (sizeof...(hostArrays) > 0)
				allocateHostArrays(dim1, dim2, hostArrays...);

		}
		template<typename T, typename... Ts>
		void allocateHostArrays(size_t dim1, size_t dim2, size_t dim3, T***& hostArray, Ts***&... hostArrays) {

			hostArray = new T * *[dim1]();
			hostArray[0] = new T * [dim1 * dim2]();
			hostArray[0][0] = new T[dim1 * dim2 * dim3]();

			for (int i = 1; i < dim1; i++)
				hostArray[i] = hostArray[i - 1] + dim2;
			for (int i = 1; i < dim1 * dim2; i++)
				hostArray[0][i] = hostArray[0][i - 1] + dim3;

			if constexpr (sizeof...(hostArrays) > 0)
				allocateHostArrays(dim1, dim2, dim3, hostArrays...);

		}
		template<typename T, typename... Ts>
		void releaseHostArrays(T**& hostArray, Ts**&... hostArrays) {

			delete[] hostArray[0];
			delete[] hostArray;

			if constexpr (sizeof...(hostArrays) > 0)
				releaseHostArrays(hostArrays...);

		}
		template<typename T, typename... Ts>
		void releaseHostArrays(T***& hostArray, Ts***&... hostArrays) {

			delete[] hostArray[0][0];
			delete[] hostArray[0];
			delete[] hostArray;

			if constexpr (sizeof...(hostArrays) > 0)
				releaseHostArrays(hostArrays...);

		}

		template<typename T, typename... Ts>
		void allocateDeviceArrays(int localId, int devNums, size_t size, T**& devArray, Ts**&... devArrays) {

			devArray = new T * [devNums]();

			for (int i = 0; i < devNums; i++) {
				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMalloc((void**)&devArray[i], sizeof(T) * size));
			}

			if constexpr (sizeof...(devArrays) > 0)
				allocateDeviceArrays(localId, devNums, size, devArrays...);

		}
		template<typename T, typename... Ts>
		void releaseDeviceArrays(int localId, int devNums, T**& devArray, Ts**&... devArrays) {

			for (int i = 0; i < devNums; i++) {
				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaFree(devArray[i]));
			}

			if constexpr (sizeof...(devArrays) > 0)
				releaseDeviceArrays(localId, devNums, devArrays...);

		}

		/*--------------------Memory Copy Between Host and Device--------------------*/

		template<typename T, typename... Ts>
		void hostToDevice(size_t size, size_t devOffset, size_t hostOffset, T*& devArray, T*& hostArray, Ts*&... devHostArrays) {

			CUDACHECK(cudaMemcpy(devArray + devOffset, hostArray + hostOffset, sizeof(T) * size, cudaMemcpyHostToDevice));
			if constexpr (sizeof...(devHostArrays) > 0)
				hostToDevice(size, devOffset, hostOffset, devHostArrays...);

		}
		template<typename T, typename... Ts>
		void deviceToHost(size_t size, size_t devOffset, size_t hostOffset, T*& devArray, T*& hostArray, Ts*&... devHostArrays) {

			CUDACHECK(cudaMemcpy(hostArray + hostOffset, devArray + devOffset, sizeof(T) * size, cudaMemcpyDeviceToHost));
			if constexpr (sizeof...(devHostArrays) > 0)
				deviceToHost(size, devOffset, hostOffset, devHostArrays...);

		}

		/*----------------------------Help Initialize Equilibrium --------------------------*/

		template<typename T, typename... Ts>
		void binaryToHost(size_t offset, size_t dim1, size_t dim2, std::vector<double>& binary, T**& hostArray, Ts**&... hostArrays) {

			for (int i = 0; i < dim1; i++) {
				for (int j = 0; j < dim2; j++) {
					hostArray[i][j] = binary[offset + i * dim2 + j];
				}
			}

			if constexpr (sizeof...(hostArrays) > 0)
				binaryToHost(offset + dim1 * dim2, dim1, dim2, binary, hostArrays...);

		}

	};
	void allocateHostMemory() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Allocate host memory." << RESET << std::endl;

		Allocator HostAllocator;

		/*----------------------------MHD Equilibrium on CPU----------------------------*/

		HostAllocator.allocateHostArrays(gridNx, gridNy,
			q, q_px, psip, psip_px,
			Ni, Ni_px, Ti, Ti_px, Pi, Pi_px, Ne, Ne_px, Te, Te_px, Pe, Pe_px,
			Na, Na_px, Ta, Ta_px, Nb, Nb_px, Tb, Tb_px,
			B, B_px, B_py, B_px2, B_pxy, B_py2, J, J_px, J_py, Bny, Bny_px, Bny_py,
			Va, Va_px, Va_py, Rho, Rho_px, Rho_py, JpB, JpB_px, JpB_py, R, Z,
			gconxx, gconxx_px, gconxx_py, gconxy, gconxy_px, gconxy_py,
			gconxz, gconxz_px, gconxz_py, gconyy, gconyy_px, gconyy_py,
			gconyz, gconyz_px, gconyz_py, gconzz, gconzz_px, gconzz_py,
			gcovxx, gcovxx_px, gcovxx_py, gcovxy, gcovxy_px, gcovxy_py,
			gcovxz, gcovxz_px, gcovxz_py, gcovyy, gcovyy_px, gcovyy_py,
			gcovyz, gcovyz_px, gcovyz_py, gcovzz, gcovzz_px, gcovzz_py);

		HostAllocator.allocateHostArrays(gridNx, gridNyPlusGhost,
			SFAconxx, SFAconxx_px, SFAconxx_py, SFAconxy, SFAconxy_px, SFAconxy_py,
			SFAconxz, SFAconxz_px, SFAconxz_py, SFAconyy, SFAconyy_px, SFAconyy_py,
			SFAconyz, SFAconyz_px, SFAconyz_py, SFAconzz, SFAconzz_px, SFAconzz_py,
			SFAcovxx, SFAcovxx_px, SFAcovxx_py, SFAcovxy, SFAcovxy_px, SFAcovxy_py,
			SFAcovxz, SFAcovxz_px, SFAcovxz_py, SFAcovyy, SFAcovyy_px, SFAcovyy_py,
			SFAcovyz, SFAcovyz_px, SFAcovyz_py, SFAcovzz, SFAcovzz_px, SFAcovzz_py);

		/*----------------------------MHD Perturbation on CPU----------------------------*/

		HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx, h_qtheta);
		HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx, gridNz, h_w, h_A, h_dNe, h_dTe, h_Phi, h_dJpB, h_dPe);

		/*-------------------------Coefficient Compression on CPU-------------------------*/

		/*-------------------------Linear-------------------------*/

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_A_w, h_A_px_w, h_A_py_w, h_A_pz_w,
			h_dJpB_w, h_dJpB_px_w, h_dJpB_py_w, h_dJpB_pz_w,
			h_dP_w, h_dP_px_w, h_dP_py_w, h_dP_pz_w,
			h_w_py_w, h_w_pz_w, h_w_Phi);

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_Phi_w,
			h_Phi_px_w, h_Phi_pz_w,
			h_Phi_px2_w, h_Phi_pxz_w, h_Phi_pz2_w);

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_A_resistive,
			h_A_px_resistive, h_A_pz_resistive,
			h_A_px2_resistive, h_A_pxz_resistive, h_A_pz2_resistive);

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_F_perp2,
			h_F_px_perp2, h_F_pz_perp2,
			h_F_px2_perp2, h_F_pxz_perp2, h_F_pz2_perp2);

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_A_dJpB,
			h_A_px_dJpB, h_A_pz_dJpB,
			h_A_px2_dJpB, h_A_pxz_dJpB, h_A_pz2_dJpB);

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_Phi_A, h_Phi_px_A, h_Phi_py_A, h_Phi_pz_A,
			h_dNe_A, h_dNe_px_A, h_dNe_py_A, h_dNe_pz_A,
			h_A_A, h_A_px_A, h_A_py_A, h_A_pz_A);

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_Phi_dNe, h_Phi_px_dNe, h_Phi_py_dNe, h_Phi_pz_dNe,
			h_dPe_dNe, h_dPe_px_dNe, h_dPe_py_dNe, h_dPe_pz_dNe,
			h_dJpB_dNe, h_dJpB_px_dNe, h_dJpB_py_dNe, h_dJpB_pz_dNe,
			h_A_dNe, h_A_px_dNe, h_A_py_dNe, h_A_pz_dNe);

		HostAllocator.allocateHostArrays(gridNy, gridNx,
			h_Phi_dTe, h_Phi_px_dTe, h_Phi_py_dTe, h_Phi_pz_dTe,
			h_dTe_dTe, h_dTe_px_dTe, h_dTe_py_dTe, h_dTe_pz_dTe,
			h_dNe_dTe, h_dNe_px_dTe, h_dNe_py_dTe, h_dNe_pz_dTe,
			h_Ne0, h_Te0, h_Ne0_px, h_Te0_px, h_Pe0_px);

		HostAllocator.allocateHostArrays(gridNy, gridNx, 5, h_F2perp2);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_A2dJpB);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 5, h_Phi2w);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 10, h_wdPAdJpB2w);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 5, h_APhidNe2A);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 11, h_dPePhiAdJpB2dNe);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_PhidTedNe2dTe);

		/*-----------------------Nonlinear-----------------------*/

		HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_wPhi_w);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 9, h_AdJpB_w);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 9, h_PhiA_A, h_NeA_A);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 9, h_AdJpB_dNe, h_dNePhi_dNe);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 6, h_PhiTe_dTe);
		HostAllocator.allocateHostArrays(gridNy, gridNx, 18, h_PhiTeA_dTe);

		/*----------------------------------Particle in Cell----------------------------------*/

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

		if (ifSlowing) {

			HostAllocator.allocateHostArrays(devNums, (size_t)randMax, h_rand_keys);
			HostAllocator.allocateHostArrays(devNums, (size_t)randMax * 4, h_rand_values);

		}

		HostAllocator.allocateHostArrays(cellNx, 30, h_pic1d);
		HostAllocator.allocateHostArrays(cellNy * cellNx, 72, h_pic2d);
		HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx * gridNzPlusGhost * 8, h_pic3d);

		HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNx, gridNz, h_globalPi, h_globalPa, h_globalPb);

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void releaseHostMemory() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Release host memory." << RESET << std::endl;

		Allocator HostAllocator;

		/*----------------------------MHD Equilibrium on CPU----------------------------*/

		HostAllocator.releaseHostArrays(
			q, q_px, psip, psip_px,
			Ni, Ni_px, Ti, Ti_px, Pi, Pi_px, Ne, Ne_px, Te, Te_px, Pe, Pe_px,
			Na, Na_px, Ta, Ta_px, Nb, Nb_px, Tb, Tb_px,
			B, B_px, B_py, B_px2, B_pxy, B_py2, J, J_px, J_py, Bny, Bny_px, Bny_py,
			Va, Va_px, Va_py, Rho, Rho_px, Rho_py, JpB, JpB_px, JpB_py, R, Z,
			gconxx, gconxx_px, gconxx_py, gconxy, gconxy_px, gconxy_py,
			gconxz, gconxz_px, gconxz_py, gconyy, gconyy_px, gconyy_py,
			gconyz, gconyz_px, gconyz_py, gconzz, gconzz_px, gconzz_py,
			gcovxx, gcovxx_px, gcovxx_py, gcovxy, gcovxy_px, gcovxy_py,
			gcovxz, gcovxz_px, gcovxz_py, gcovyy, gcovyy_px, gcovyy_py,
			gcovyz, gcovyz_px, gcovyz_py, gcovzz, gcovzz_px, gcovzz_py);

		HostAllocator.releaseHostArrays(
			SFAconxx, SFAconxx_px, SFAconxx_py, SFAconxy, SFAconxy_px, SFAconxy_py,
			SFAconxz, SFAconxz_px, SFAconxz_py, SFAconyy, SFAconyy_px, SFAconyy_py,
			SFAconyz, SFAconyz_px, SFAconyz_py, SFAconzz, SFAconzz_px, SFAconzz_py,
			SFAcovxx, SFAcovxx_px, SFAcovxx_py, SFAcovxy, SFAcovxy_px, SFAcovxy_py,
			SFAcovxz, SFAcovxz_px, SFAcovxz_py, SFAcovyy, SFAcovyy_px, SFAcovyy_py,
			SFAcovyz, SFAcovyz_px, SFAcovyz_py, SFAcovzz, SFAcovzz_px, SFAcovzz_py);

		/*----------------------------MHD Perturbation on CPU----------------------------*/

		HostAllocator.releaseHostArrays(h_qtheta);
		HostAllocator.releaseHostArrays(h_w, h_A, h_dNe, h_dTe, h_Phi, h_dJpB, h_dPe);

		/*-------------------------Coefficient Compression on CPU------------------------*/

		/*-------------------------Linear-------------------------*/

		HostAllocator.releaseHostArrays(
			h_A_w, h_A_px_w, h_A_py_w, h_A_pz_w,
			h_dJpB_w, h_dJpB_px_w, h_dJpB_py_w, h_dJpB_pz_w,
			h_dP_w, h_dP_px_w, h_dP_py_w, h_dP_pz_w,
			h_w_py_w, h_w_pz_w, h_w_Phi);

		HostAllocator.releaseHostArrays(
			h_Phi_w,
			h_Phi_px_w, h_Phi_pz_w,
			h_Phi_px2_w, h_Phi_pxz_w, h_Phi_pz2_w);

		HostAllocator.releaseHostArrays(
			h_A_resistive,
			h_A_px_resistive, h_A_pz_resistive,
			h_A_px2_resistive, h_A_pxz_resistive, h_A_pz2_resistive);

		HostAllocator.releaseHostArrays(
			h_F_perp2,
			h_F_px_perp2, h_F_pz_perp2,
			h_F_px2_perp2, h_F_pxz_perp2, h_F_pz2_perp2);

		HostAllocator.releaseHostArrays(
			h_A_dJpB,
			h_A_px_dJpB, h_A_pz_dJpB,
			h_A_px2_dJpB, h_A_pxz_dJpB, h_A_pz2_dJpB);

		HostAllocator.releaseHostArrays(
			h_Phi_A, h_Phi_px_A, h_Phi_py_A, h_Phi_pz_A,
			h_dNe_A, h_dNe_px_A, h_dNe_py_A, h_dNe_pz_A,
			h_A_A, h_A_px_A, h_A_py_A, h_A_pz_A);

		HostAllocator.releaseHostArrays(
			h_Phi_dNe, h_Phi_px_dNe, h_Phi_py_dNe, h_Phi_pz_dNe,
			h_dPe_dNe, h_dPe_px_dNe, h_dPe_py_dNe, h_dPe_pz_dNe,
			h_dJpB_dNe, h_dJpB_px_dNe, h_dJpB_py_dNe, h_dJpB_pz_dNe,
			h_A_dNe, h_A_px_dNe, h_A_py_dNe, h_A_pz_dNe);

		HostAllocator.releaseHostArrays(
			h_Phi_dTe, h_Phi_px_dTe, h_Phi_py_dTe, h_Phi_pz_dTe,
			h_dTe_dTe, h_dTe_px_dTe, h_dTe_py_dTe, h_dTe_pz_dTe,
			h_dNe_dTe, h_dNe_px_dTe, h_dNe_py_dTe, h_dNe_pz_dTe,
			h_Ne0, h_Te0, h_Ne0_px, h_Te0_px, h_Pe0_px);

		HostAllocator.releaseHostArrays(h_F2perp2);
		HostAllocator.releaseHostArrays(h_A2dJpB);
		HostAllocator.releaseHostArrays(h_Phi2w);
		HostAllocator.releaseHostArrays(h_wdPAdJpB2w);
		HostAllocator.releaseHostArrays(h_APhidNe2A);
		HostAllocator.releaseHostArrays(h_dPePhiAdJpB2dNe);
		HostAllocator.releaseHostArrays(h_PhidTedNe2dTe);

		/*-----------------------Nonlinear-----------------------*/

		HostAllocator.releaseHostArrays(
			h_wPhi_w, h_AdJpB_w,
			h_PhiA_A, h_NeA_A,
			h_AdJpB_dNe, h_dNePhi_dNe,
			h_PhiTe_dTe, h_PhiTeA_dTe);

		/*----------------------------------Particle in Cell----------------------------------*/

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

		if (ifSlowing) {

			HostAllocator.releaseHostArrays(h_rand_keys);
			HostAllocator.releaseHostArrays(h_rand_values);

		}

		HostAllocator.releaseHostArrays(h_pic1d, h_pic2d, h_pic3d);
		HostAllocator.releaseHostArrays(h_globalPi, h_globalPa, h_globalPb);

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void allocateDeviceMemory() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Allocate device memory." << RESET << std::endl;

		Allocator DeviceAllocator;

		/*----------------------------MHD Perturbation on GPU----------------------------*/

		DeviceAllocator.allocateDeviceArrays(localId, devNums, (devNy + 2 * gridGhost) * gridNx, d_qtheta);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, gridNy * gridNxz,
			d_w, d_A, d_dNe, d_dTe, d_Phi, d_dJpB, d_dPe);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, (devNy + 2 * gridGhost) * gridNxz,
			d_w_beg, d_w_midl, d_w_midr, d_w_end,
			d_A_beg, d_A_midl, d_A_midr, d_A_end,
			d_dNe_beg, d_dNe_midl, d_dNe_midr, d_dNe_end,
			d_dTe_beg, d_dTe_midl, d_dTe_midr, d_dTe_end,
			d_Phi_midl, d_Phi_midr,
			d_dJpB_midl, d_dJpB_midr,
			d_dPe_midl, d_dPe_midr,
			d_dPi_midl, d_dPi_midr,
			d_dPa_midl, d_dPa_midr,
			d_dPb_midl, d_dPb_midr,
			d_Apt_midl, d_Apt_midr);

		/*--------------------------Coefficient Compression on GPU--------------------------*/

		/*-------------------------Linear-------------------------*/

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx,
			d_A_w, d_A_px_w, d_A_py_w, d_A_pz_w,
			d_dJpB_w, d_dJpB_px_w, d_dJpB_py_w, d_dJpB_pz_w,
			d_dP_w, d_dP_px_w, d_dP_py_w, d_dP_pz_w,
			d_w_py_w, d_w_pz_w, d_w_Phi);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx,
			d_F_perp2,
			d_F_px_perp2, d_F_pz_perp2,
			d_F_px2_perp2, d_F_pxz_perp2, d_F_pz2_perp2);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx,
			d_A_dJpB,
			d_A_px_dJpB, d_A_pz_dJpB,
			d_A_px2_dJpB, d_A_pxz_dJpB, d_A_pz2_dJpB);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx,
			d_Phi_A, d_Phi_px_A, d_Phi_py_A, d_Phi_pz_A,
			d_dNe_A, d_dNe_px_A, d_dNe_py_A, d_dNe_pz_A,
			d_A_A, d_A_px_A, d_A_py_A, d_A_pz_A);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx,
			d_Phi_dNe, d_Phi_px_dNe, d_Phi_py_dNe, d_Phi_pz_dNe,
			d_dPe_dNe, d_dPe_px_dNe, d_dPe_py_dNe, d_dPe_pz_dNe,
			d_dJpB_dNe, d_dJpB_px_dNe, d_dJpB_py_dNe, d_dJpB_pz_dNe,
			d_A_dNe, d_A_px_dNe, d_A_py_dNe, d_A_pz_dNe);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx,
			d_Phi_dTe, d_Phi_px_dTe, d_Phi_py_dTe, d_Phi_pz_dTe,
			d_dTe_dTe, d_dTe_px_dTe, d_dTe_py_dTe, d_dTe_pz_dTe,
			d_dNe_dTe, d_dNe_px_dTe, d_dNe_py_dTe, d_dNe_pz_dTe,
			d_Ne0, d_Te0, d_Ne0_px, d_Te0_px, d_Pe0_px);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, d_F2perp2);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_A2dJpB);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, d_Phi2w);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 10, d_wdPAdJpB2w);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, d_APhidNe2A);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 11, d_dPePhiAdJpB2dNe);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_PhidTedNe2dTe);

		/*-----------------------Nonlinear-----------------------*/

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_wPhi_w);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, d_AdJpB_w);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, d_PhiA_A, d_NeA_A);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, d_AdJpB_dNe, d_dNePhi_dNe);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 6, d_PhiTe_dTe);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 18, d_PhiTeA_dTe);

		/*----------------------------------Particle in Cell----------------------------------*/

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx,
			dpic_Phi_py_A, dpic_dNe_py_A,
			dpic_A_A, dpic_A_py_A, dpic_A_pz_A);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 5, dpic_APhidNe2A);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * gridNx * 9, dpic_PhiA_A, dpic_NeA_A);

		if (ifIon) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, 8, d_Ion_offsets);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Ion_keys_in, d_Ion_keys_out);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Ion_values_in, d_Ion_values_out);

		}

		if (ifAlpha) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, 8, d_Alpha_offsets);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Alpha_keys_in, d_Alpha_keys_out);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Alpha_values_in, d_Alpha_values_out);

		}

		if (ifBeam) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, 8, d_Beam_offsets);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Beam_keys_in, d_Beam_keys_out);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picDev * 7, d_Beam_values_in, d_Beam_values_out);

		}

		if (ifSlowing) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)randMax, d_rand_keys);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)randMax * 4, d_rand_values);

		}

		DeviceAllocator.allocateDeviceArrays(localId, devNums, cellNx * 30, d_pic1d);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, cellNy * cellNx * 72, d_pic2d);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, gridNyPlusGhost * gridNx * gridNzPlusGhost * 8, d_pic3d);

		DeviceAllocator.allocateDeviceArrays(localId, devNums, gridNyPlusGhost * gridNxz, d_globalA, d_globalPhi, d_globalApt, d_globalPi, d_globalPa, d_globalPb);


		if (hostId == 0) {
			size_t avail, total, used;
			CUDACHECK(cudaSetDevice(localId * devNums));
			CUDACHECK(cudaMemGetInfo(&avail, &total));
			used = total - avail;
			std::cout << BOLDYELLOW << "Device memory used: " << (double)used / 1024 / 1024 / 1024 << " GB." << RESET << std::endl;
		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void releaseDeviceMemory() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Release device memory." << RESET << std::endl;

		Allocator DeviceAllocator;

		/*----------------------------MHD Perturbation on GPU----------------------------*/

		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_qtheta);

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_w, d_A, d_dNe, d_dTe, d_Phi, d_dJpB, d_dPe);

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_w_beg, d_w_midl, d_w_midr, d_w_end,
			d_A_beg, d_A_midl, d_A_midr, d_A_end,
			d_dNe_beg, d_dNe_midl, d_dNe_midr, d_dNe_end,
			d_dTe_beg, d_dTe_midl, d_dTe_midr, d_dTe_end,
			d_Phi_midl, d_Phi_midr,
			d_dJpB_midl, d_dJpB_midr,
			d_dPe_midl, d_dPe_midr,
			d_dPi_midl, d_dPi_midr,
			d_dPa_midl, d_dPa_midr,
			d_dPb_midl, d_dPb_midr,
			d_Apt_midl, d_Apt_midr);

		/*--------------------------Coefficient Compression on GPU--------------------------*/

		/*-------------------------Linear-------------------------*/

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_A_w, d_A_px_w, d_A_py_w, d_A_pz_w,
			d_dJpB_w, d_dJpB_px_w, d_dJpB_py_w, d_dJpB_pz_w,
			d_dP_w, d_dP_px_w, d_dP_py_w, d_dP_pz_w,
			d_w_py_w, d_w_pz_w, d_w_Phi);

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_F_perp2,
			d_F_px_perp2, d_F_pz_perp2,
			d_F_px2_perp2, d_F_pxz_perp2, d_F_pz2_perp2);

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_A_dJpB,
			d_A_px_dJpB, d_A_pz_dJpB,
			d_A_px2_dJpB, d_A_pxz_dJpB, d_A_pz2_dJpB);

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_Phi_A, d_Phi_px_A, d_Phi_py_A, d_Phi_pz_A,
			d_dNe_A, d_dNe_px_A, d_dNe_py_A, d_dNe_pz_A,
			d_A_A, d_A_px_A, d_A_py_A, d_A_pz_A);

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_Phi_dNe, d_Phi_px_dNe, d_Phi_py_dNe, d_Phi_pz_dNe,
			d_dPe_dNe, d_dPe_px_dNe, d_dPe_py_dNe, d_dPe_pz_dNe,
			d_dJpB_dNe, d_dJpB_px_dNe, d_dJpB_py_dNe, d_dJpB_pz_dNe,
			d_A_dNe, d_A_px_dNe, d_A_py_dNe, d_A_pz_dNe);

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_Phi_dTe, d_Phi_px_dTe, d_Phi_py_dTe, d_Phi_pz_dTe,
			d_dTe_dTe, d_dTe_px_dTe, d_dTe_py_dTe, d_dTe_pz_dTe,
			d_dNe_dTe, d_dNe_px_dTe, d_dNe_py_dTe, d_dNe_pz_dTe,
			d_Ne0, d_Te0, d_Ne0_px, d_Te0_px, d_Pe0_px);

		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_F2perp2);
		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_A2dJpB);
		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_Phi2w);
		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_wdPAdJpB2w);
		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_APhidNe2A);
		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPePhiAdJpB2dNe);
		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_PhidTedNe2dTe);

		/*-----------------------Nonlinear-----------------------*/

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			d_wPhi_w, d_AdJpB_w,
			d_PhiA_A, d_NeA_A,
			d_AdJpB_dNe, d_dNePhi_dNe,
			d_PhiTe_dTe, d_PhiTeA_dTe);

		/*------------------------------Inverse Matrix on GPU------------------------------*/

		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_laplacianCsrR, d_laplacianCsrC, d_laplacianCsrV);

		if (std::get<0>(nablaPerp2A))
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_resistiveCsrR, d_resistiveCsrC, d_resistiveCsrV);
		if (std::get<0>(nablaPerp2w))
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_wCsrR, d_wCsrC, d_wCsrV);
		if (std::get<0>(nablaPerp2dNe))
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dNeCsrR, d_dNeCsrC, d_dNeCsrV);
		if (std::get<0>(nablaPerp2dTe))
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dTeCsrR, d_dTeCsrC, d_dTeCsrV);
		if (std::get<0>(nablaPerp2dPi))
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPiCsrR, d_dPiCsrC, d_dPiCsrV);
		if (std::get<0>(nablaPerp2dPa))
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPaCsrR, d_dPaCsrC, d_dPaCsrV);
		if (std::get<0>(nablaPerp2dPb))
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_dPbCsrR, d_dPbCsrC, d_dPbCsrV);

		for (int i = 0; i < devNums; i++) {

			CUDACHECK(cudaSetDevice(localId * devNums + i));

			for (int j = 0; j < devNy; j++) {

				CUDSSCHECK(cudssConfigDestroy(laplacianConfigs[i][j]));
				CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], laplacianDatas[i][j]));
				CUDSSCHECK(cudssMatrixDestroy(laplacianAs[i][j]));
				CUDSSCHECK(cudssMatrixDestroy(laplacianXs[i][j]));
				CUDSSCHECK(cudssMatrixDestroy(laplacianBs[i][j]));

				if (std::get<0>(nablaPerp2A)) {

					CUDSSCHECK(cudssConfigDestroy(resistiveConfigs[i][j]));
					CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], resistiveDatas[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(resistiveAs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(resistiveXs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(resistiveBs[i][j]));

				}

				if (std::get<0>(nablaPerp2w)) {

					CUDSSCHECK(cudssConfigDestroy(wConfigs[i][j]));
					CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], wDatas[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(wAs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(wXs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(wBs[i][j]));

				}

				if (std::get<0>(nablaPerp2dNe)) {

					CUDSSCHECK(cudssConfigDestroy(dNeConfigs[i][j]));
					CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dNeDatas[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dNeAs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dNeXs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dNeBs[i][j]));

				}

				if (std::get<0>(nablaPerp2dTe)) {

					CUDSSCHECK(cudssConfigDestroy(dTeConfigs[i][j]));
					CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dTeDatas[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dTeAs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dTeXs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dTeBs[i][j]));

				}

				if (std::get<0>(nablaPerp2dPi)) {

					CUDSSCHECK(cudssConfigDestroy(dPiConfigs[i][j]));
					CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dPiDatas[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dPiAs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dPiXs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dPiBs[i][j]));

				}

				if (std::get<0>(nablaPerp2dPa)) {

					CUDSSCHECK(cudssConfigDestroy(dPaConfigs[i][j]));
					CUDSSCHECK(cudssDataDestroy(cudssHandles[i][j], dPaDatas[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dPaAs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dPaXs[i][j]));
					CUDSSCHECK(cudssMatrixDestroy(dPaBs[i][j]));

				}

				if (std::get<0>(nablaPerp2dPb)) {

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

		/*----------------------------------Particle in Cell----------------------------------*/

		DeviceAllocator.releaseDeviceArrays(localId, devNums,
			dpic_Phi_py_A, dpic_dNe_py_A,
			dpic_A_A, dpic_A_py_A, dpic_A_pz_A);

		DeviceAllocator.releaseDeviceArrays(localId, devNums, dpic_APhidNe2A);

		DeviceAllocator.releaseDeviceArrays(localId, devNums, dpic_PhiA_A, dpic_NeA_A);

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

		if (ifSlowing) {

			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_rand_keys);
			DeviceAllocator.releaseDeviceArrays(localId, devNums, d_rand_values);

		}

		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_pic1d, d_pic2d, d_pic3d);
		DeviceAllocator.releaseDeviceArrays(localId, devNums, d_globalA, d_globalPhi, d_globalApt, d_globalPi, d_globalPa, d_globalPb);

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void memcpyHostToDevice() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Memory copy from host to device." << RESET << std::endl;

		Allocator H2DAllocator;

		for (int i = 0; i < devNums; i++) {

			CUDACHECK(cudaSetDevice(localId * devNums + i));

			H2DAllocator.hostToDevice((devNy + 2 * gridGhost) * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx, d_qtheta[i], h_qtheta[0]);

			/*------------------------------------------Linear------------------------------------------*/

			/*------------------------------Vorticity------------------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				//Perturbed Parallel Vector Potential in Vorticity
				d_A_w[i], h_A_w[0], d_A_px_w[i], h_A_px_w[0], d_A_py_w[i], h_A_py_w[0], d_A_pz_w[i], h_A_pz_w[0],
				//Perturbed Parallel Current in Vorticity
				d_dJpB_w[i], h_dJpB_w[0], d_dJpB_px_w[i], h_dJpB_px_w[0], d_dJpB_py_w[i], h_dJpB_py_w[0], d_dJpB_pz_w[i], h_dJpB_pz_w[0],
				//Perturbed Pressure in Vorticity
				d_dP_w[i], h_dP_w[0], d_dP_px_w[i], h_dP_px_w[0], d_dP_py_w[i], h_dP_py_w[0], d_dP_pz_w[i], h_dP_pz_w[0],
				//Ion Diamagnetic Drift and Finite Larmor Radius in Vorticity
				d_w_py_w[i], h_w_py_w[0], d_w_pz_w[i], h_w_pz_w[0], d_w_Phi[i], h_w_Phi[0]);

			/*---------------------Perturbed Parallel Current---------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				//Perturbed Parallel Vector Potential in Perturbed Parallel Current
				d_A_dJpB[i], h_A_dJpB[0], d_A_px_dJpB[i], h_A_px_dJpB[0], d_A_pz_dJpB[i], h_A_pz_dJpB[0],
				d_A_px2_dJpB[i], h_A_px2_dJpB[0], d_A_pxz_dJpB[i], h_A_pxz_dJpB[0], d_A_pz2_dJpB[i], h_A_pz2_dJpB[0]);

			/*-----------------Perturbed Parallel Vector Potential-----------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				//Perturbed Electric Potential in Perturbed Parallel Vector Potential
				d_Phi_A[i], h_Phi_A[0], d_Phi_px_A[i], h_Phi_px_A[0], d_Phi_py_A[i], h_Phi_py_A[0], d_Phi_pz_A[i], h_Phi_pz_A[0],
				//Perturbed Density in Perturbed Parallel Vector Potential
				d_dNe_A[i], h_dNe_A[0], d_dNe_px_A[i], h_dNe_px_A[0], d_dNe_py_A[i], h_dNe_py_A[0], d_dNe_pz_A[i], h_dNe_pz_A[0],
				//Perturbed Parallel Vector Potential in Perturbed Parallel Vector Potential
				d_A_A[i], h_A_A[0], d_A_px_A[i], h_A_px_A[0], d_A_py_A[i], h_A_py_A[0], d_A_pz_A[i], h_A_pz_A[0]);

			/*-------------------------Perturbed Density-------------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				//Perturbed Electric Potential in Perturbed Density
				d_Phi_dNe[i], h_Phi_dNe[0], d_Phi_px_dNe[i], h_Phi_px_dNe[0], d_Phi_py_dNe[i], h_Phi_py_dNe[0], d_Phi_pz_dNe[i], h_Phi_pz_dNe[0],
				//Perturbed Pressure in Perturbed Density
				d_dPe_dNe[i], h_dPe_dNe[0], d_dPe_px_dNe[i], h_dPe_px_dNe[0], d_dPe_py_dNe[i], h_dPe_py_dNe[0], d_dPe_pz_dNe[i], h_dPe_pz_dNe[0],
				//Perturbed Parallel Current in Perturbed Density
				d_dJpB_dNe[i], h_dJpB_dNe[0], d_dJpB_px_dNe[i], h_dJpB_px_dNe[0], d_dJpB_py_dNe[i], h_dJpB_py_dNe[0], d_dJpB_pz_dNe[i], h_dJpB_pz_dNe[0],
				//Perturbed Parallel Vector Potential in Perturbed Density
				d_A_dNe[i], h_A_dNe[0], d_A_px_dNe[i], h_A_px_dNe[0], d_A_py_dNe[i], h_A_py_dNe[0], d_A_pz_dNe[i], h_A_pz_dNe[0]);

			/*-----------------------Perturbed Temperature-----------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				//Perturbed Electric Potential in Perturbed Temperature
				d_Phi_dTe[i], h_Phi_dTe[0], d_Phi_px_dTe[i], h_Phi_px_dTe[0], d_Phi_py_dTe[i], h_Phi_py_dTe[0], d_Phi_pz_dTe[i], h_Phi_pz_dTe[0],
				//Perturbed Temperature in Perturbed Temperature
				d_dTe_dTe[i], h_dTe_dTe[0], d_dTe_px_dTe[i], h_dTe_px_dTe[0], d_dTe_py_dTe[i], h_dTe_py_dTe[0], d_dTe_pz_dTe[i], h_dTe_pz_dTe[0],
				//Perturbed Density in Perturbed Temperature
				d_dNe_dTe[i], h_dNe_dTe[0], d_dNe_px_dTe[i], h_dNe_px_dTe[0], d_dNe_py_dTe[i], h_dNe_py_dTe[0], d_dNe_pz_dTe[i], h_dNe_pz_dTe[0]);

			/*---------------Equilibrium Density and Temperature---------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				d_Ne0[i], h_Ne0[0], d_Te0[i], h_Te0[0], d_Ne0_px[i], h_Ne0_px[0], d_Te0_px[i], h_Te0_px[0], d_Pe0_px[i], h_Pe0_px[0]);

			/*-----------------------------wAdNedTe-----------------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_F2perp2[i], h_F2perp2[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_A2dJpB[i], h_A2dJpB[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_Phi2w[i], h_Phi2w[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 10, 0, (hostId * hostNy + i * devNy) * gridNx * 10, d_wdPAdJpB2w[i], h_wdPAdJpB2w[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_APhidNe2A[i], h_APhidNe2A[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 11, 0, (hostId * hostNy + i * devNy) * gridNx * 11, d_dPePhiAdJpB2dNe[i], h_dPePhiAdJpB2dNe[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_PhidTedNe2dTe[i], h_PhidTedNe2dTe[0][0]);

			/*----------------------------------------Nonlinear----------------------------------------*/

			/*------------------------------Vorticity------------------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_wPhi_w[i], h_wPhi_w[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_AdJpB_w[i], h_AdJpB_w[0][0]);

			/*-----------------Perturbed Parallel Vector Potential-----------------*/

			H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_PhiA_A[i], h_PhiA_A[0][0], d_NeA_A[i], h_NeA_A[0][0]);

			/*-------------------------Perturbed Density-------------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_AdJpB_dNe[i], h_AdJpB_dNe[0][0], d_dNePhi_dNe[i], h_dNePhi_dNe[0][0]);

			/*-----------------------Perturbed Temperature-----------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_PhiTe_dTe[i], h_PhiTe_dTe[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 18, 0, (hostId * hostNy + i * devNy) * gridNx * 18, d_PhiTeA_dTe[i], h_PhiTeA_dTe[0][0]);

			/*------------------------------------------MHD------------------------------------------*/

			H2DAllocator.hostToDevice((devNy + 2 * gridGhost) * gridNxz, 0, (hostId * hostNy + i * devNy) * gridNxz,
				d_Phi_midl[i], h_Phi[0][0], d_Phi_midr[i], h_Phi[0][0],
				d_dJpB_midl[i], h_dJpB[0][0], d_dJpB_midr[i], h_dJpB[0][0],
				d_dPe_midl[i], h_dPe[0][0], d_dPe_midr[i], h_dPe[0][0],
				d_dPi_midl[i], h_globalPi[0][0], d_dPi_midr[i], h_globalPi[0][0],
				d_dPa_midl[i], h_globalPa[0][0], d_dPa_midr[i], h_globalPa[0][0],
				d_dPb_midl[i], h_globalPb[0][0], d_dPb_midr[i], h_globalPb[0][0],
				d_Apt_midl[i], h_dPe[0][0], d_Apt_midr[i], h_dPe[0][0],
				d_w_beg[i], h_w[0][0], d_w_midl[i], h_w[0][0], d_w_midr[i], h_dPe[0][0], d_w_end[i], h_dPe[0][0],
				d_A_beg[i], h_A[0][0], d_A_midl[i], h_A[0][0], d_A_midr[i], h_dPe[0][0], d_A_end[i], h_dPe[0][0],
				d_dNe_beg[i], h_dNe[0][0], d_dNe_midl[i], h_dNe[0][0], d_dNe_midr[i], h_dPe[0][0], d_dNe_end[i], h_dPe[0][0],
				d_dTe_beg[i], h_dTe[0][0], d_dTe_midl[i], h_dTe[0][0], d_dTe_midr[i], h_dPe[0][0], d_dTe_end[i], h_dPe[0][0]);

			/*-------------------------------------------PIC-------------------------------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				dpic_Phi_py_A[i], h_Phi_py_A[0], dpic_dNe_py_A[i], h_dNe_py_A[0],
				dpic_A_A[i], h_A_A[0], dpic_A_py_A[i], h_A_py_A[0], dpic_A_pz_A[i], h_A_pz_A[0]);

			H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, dpic_APhidNe2A[i], h_APhidNe2A[0][0]);

			H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, dpic_PhiA_A[i], h_PhiA_A[0][0], dpic_NeA_A[i], h_NeA_A[0][0]);

			if (ifIon) {

				H2DAllocator.hostToDevice(8, 0, i * 8, d_Ion_offsets[i], h_Ion_offsets[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_keys_in[i], h_Ion_keys[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_keys_out[i], h_Ion_keys[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_values_in[i], h_Ion_values[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Ion_values_out[i], h_Ion_values[0]);

			}

			if (ifAlpha) {

				H2DAllocator.hostToDevice(8, 0, i * 8, d_Alpha_offsets[i], h_Alpha_offsets[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_keys_in[i], h_Alpha_keys[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_keys_out[i], h_Alpha_keys[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_values_in[i], h_Alpha_values[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Alpha_values_out[i], h_Alpha_values[0]);

			}

			if (ifBeam) {

				H2DAllocator.hostToDevice(8, 0, i * 8, d_Beam_offsets[i], h_Beam_offsets[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_keys_in[i], h_Beam_keys[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_keys_out[i], h_Beam_keys[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_values_in[i], h_Beam_values[0]);
				H2DAllocator.hostToDevice((size_t)picDev * 7, 0, (size_t)i * picDev * 7, d_Beam_values_out[i], h_Beam_values[0]);

			}

			if (ifSlowing) {

				H2DAllocator.hostToDevice((size_t)randMax, 0, (size_t)i * randMax, d_rand_keys[i], h_rand_keys[0]);
				H2DAllocator.hostToDevice((size_t)randMax * 4, 0, (size_t)i * randMax * 4, d_rand_values[i], h_rand_values[0]);

			}

			H2DAllocator.hostToDevice(cellNx * 30, 0, 0, d_pic1d[i], h_pic1d[0]);
			H2DAllocator.hostToDevice(cellNy * cellNx * 72, 0, 0, d_pic2d[i], h_pic2d[0]);
			H2DAllocator.hostToDevice(gridNyPlusGhost * gridNx * gridNzPlusGhost * 8, 0, 0, d_pic3d[i], h_pic3d[0]);

			H2DAllocator.hostToDevice(gridNyPlusGhost * gridNxz, 0, 0, d_globalPi[i], h_dPe[0][0], d_globalPa[i], h_dPe[0][0], d_globalPb[i], h_dPe[0][0]);

			H2DAllocator.hostToDevice(gridNyPlusGhost * gridNxz, 0, 0, d_globalA[i], h_dPe[0][0], d_globalPhi[i], h_dPe[0][0], d_globalApt[i], h_dPe[0][0]);

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}

	void loadMHDEquilibrium(std::string file) {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Load MHD equilibrium." << RESET << std::endl;

		std::ifstream input;
		size_t bytes;
		size_t length;

		input.open(file, std::ios::in | std::ios::binary);
		input.seekg(0, std::ios::end);
		bytes = input.tellg();
		length = bytes / sizeof(double);
		std::vector<double> binaryData(length);
		input.seekg(0, std::ios::beg);
		input.read((char*)(&binaryData[0]), bytes);
		input.close();

		for (int i = 0; i < gridNx; i++) {
			for (int j = 0; j < gridNyPlusGhost; j++) {
				h_qtheta[j][i] = binaryData[i * gridNyPlusGhost + j];
			}
		}

		Allocator B2HAllocator;

		B2HAllocator.binaryToHost(gridNx * gridNyPlusGhost, gridNx, gridNyPlusGhost, binaryData,
			SFAconxx, SFAconxx_px, SFAconxx_py, SFAconxy, SFAconxy_px, SFAconxy_py,
			SFAconyy, SFAconyy_px, SFAconyy_py, SFAconxz, SFAconxz_px, SFAconxz_py,
			SFAconyz, SFAconyz_px, SFAconyz_py, SFAconzz, SFAconzz_px, SFAconzz_py,
			SFAcovxx, SFAcovxx_px, SFAcovxx_py, SFAcovxy, SFAcovxy_px, SFAcovxy_py,
			SFAcovyy, SFAcovyy_px, SFAcovyy_py, SFAcovxz, SFAcovxz_px, SFAcovxz_py,
			SFAcovyz, SFAcovyz_px, SFAcovyz_py, SFAcovzz, SFAcovzz_px, SFAcovzz_py);

		B2HAllocator.binaryToHost(gridNx * gridNyPlusGhost * 37, gridNx, gridNy, binaryData,
			q, q_px, psip, psip_px,
			Ni, Ni_px, Ti, Ti_px, Pi, Pi_px, Ne, Ne_px, Te, Te_px, Pe, Pe_px,
			Na, Na_px, Ta, Ta_px, Nb, Nb_px, Tb, Tb_px,
			gconxx, gconxx_px, gconxx_py, gconxy, gconxy_px, gconxy_py,
			gconyy, gconyy_px, gconyy_py, gconxz, gconxz_px, gconxz_py,
			gconyz, gconyz_px, gconyz_py, gconzz, gconzz_px, gconzz_py,
			gcovxx, gcovxx_px, gcovxx_py, gcovxy, gcovxy_px, gcovxy_py,
			gcovyy, gcovyy_px, gcovyy_py, gcovxz, gcovxz_px, gcovxz_py,
			gcovyz, gcovyz_px, gcovyz_py, gcovzz, gcovzz_px, gcovzz_py,
			J, J_px, J_py, Bny, Bny_px, Bny_py, JpB, JpB_px, JpB_py,
			Rho, Rho_px, Rho_py, Va, Va_px, Va_py, R, Z,
			B, B_px, B_py, B_px2, B_pxy, B_py2);

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
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void loadMHDPerturbation(std::string file) {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Load MHD perturbation." << RESET << std::endl;

		std::ifstream input;
		size_t bytes;
		size_t length;

		input.open(file, std::ios::in | std::ios::binary);
		input.seekg(0, std::ios::end);
		bytes = input.tellg();
		length = bytes / sizeof(double);
		std::vector<double> binaryData(length);
		input.seekg(0, std::ios::beg);
		input.read((char*)(&binaryData[0]), bytes);
		input.close();

		if (length == gridNy * gridNxz) {

			for (int j = 0; j < gridNy; j++) {
				for (int i = 0; i < gridNx; i++) {
					for (int k = 0; k < gridNz; k++) {
						if (i != 0 && i != gridNx - 1)
							h_Phi[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k];
					}
				}
			}

		}
		else {

			for (int j = 0; j < gridNy; j++) {
				for (int i = 0; i < gridNx; i++) {
					for (int k = 0; k < gridNz; k++) {
						if (i != 0 && i != gridNx - 1) {
							h_w[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 0 * gridNy * gridNxz];
							h_A[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 1 * gridNy * gridNxz];
							h_dNe[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 2 * gridNy * gridNxz];
							h_dTe[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 3 * gridNy * gridNxz];
							h_globalPi[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 4 * gridNy * gridNxz];
							h_globalPa[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 5 * gridNy * gridNxz];
							h_globalPb[j + gridGhost][i][k] = binaryData[j * gridNxz + i * gridNz + k + 6 * gridNy * gridNxz];
						}
					}
				}
			}

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void compressCollocatedCoefficient() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Compress collocated coefficient in shifted metric coordinate." << RESET << std::endl;

		/*-------------------------Linear-------------------------*/

		for (int i = 0; i < gridNx; i++) {
			for (int j = 0; j < gridNy; j++) {

				//Vorticity

				h_A_w[j][i] = pow(B[i][j], -2.0) * pow(J[i][j], -1.0) * (B[i][j] * gcovyz[i][j] * (Bny_py[i][j] * JpB_px[i][j] + (-1.0) * Bny_px[i][j] * JpB_py[i][j]) + Bny[i][j] * (gcovyz[i][j] * ((-1.0) * B_py[i][j] * JpB_px[i][j] + B_px[i][j] * JpB_py[i][j]) + B[i][j] * (gcovyz_py[i][j] * JpB_px[i][j] + (-1.0) * gcovyz_px[i][j] * JpB_py[i][j])));
				h_A_px_w[j][i] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * JpB_py[i][j];
				h_A_py_w[j][i] = pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * JpB_px[i][j];
				h_A_pz_w[j][i] = pow(B[i][j], -1.0) * Bny[i][j] * pow(J[i][j], -1.0) * ((-1.0) * gcovyy[i][j] * JpB_px[i][j] + gcovxy[i][j] * JpB_py[i][j]);

				h_dJpB_w[j][i] = 0.0;
				h_dJpB_px_w[j][i] = 0.0;
				h_dJpB_py_w[j][i] = Bny[i][j];
				h_dJpB_pz_w[j][i] = 0.0;

				h_dP_w[j][i] = 0.0;
				h_dP_px_w[j][i] = pow(B[i][j], -3.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0);
				h_dP_py_w[j][i] = pow(B[i][j], -4.0) * (B[i][j] * Bny[i][j] * B_px[i][j] * gcovyz[i][j] + (-1.0) * pow(B[i][j], 2.0) * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j]) + pow(Bny[i][j], 3.0) * ((gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]) * gcovyz[i][j] + gcovyy[i][j] * gcovyz_px[i][j] + (-1.0) * gcovxy[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0);
				h_dP_pz_w[j][i] = pow(B[i][j], -3.0) * (B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j]) + Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j] + B[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) * pow(J[i][j], -1.0);

				h_w_py_w[j][i] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(Ni[i][j], -1.0) * Pi_px[i][j] * pow(NormQE, -1.0);
				h_w_pz_w[j][i] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(Ni[i][j], -1.0) * Pi_px[i][j] * pow(NormQE, -1.0);
				h_w_Phi[j][i] = -pow(Rho[i][j], 2.0) * pow(Va[i][j], 2.0);

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

				//Poisson Equation

				h_Phi_w[j][i] = 0.0;
				h_Phi_px_w[j][i] = pow(Va[i][j], -3.0) * (pow(J[i][j], -1.0) * ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] + gconxx[i][j] * J_px[i][j] + gconxy[i][j] * J_py[i][j]) * Va[i][j] + (-2.0) * (gconxx[i][j] * Va_px[i][j] + gconxy[i][j] * Va_py[i][j]));
				h_Phi_pz_w[j][i] = pow(Va[i][j], -3.0) * (pow(J[i][j], -1.0) * ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] + gconxz[i][j] * J_px[i][j] + gconyz[i][j] * J_py[i][j]) * Va[i][j] + (-2.0) * (gconxz[i][j] * Va_px[i][j] + gconyz[i][j] * Va_py[i][j]));
				h_Phi_px2_w[j][i] = gconxx[i][j] * pow(Va[i][j], -2.0);
				h_Phi_pxz_w[j][i] = 2.0 * gconxz[i][j] * pow(Va[i][j], -2.0);
				h_Phi_pz2_w[j][i] = gconzz[i][j] * pow(Va[i][j], -2.0);
				h_Phi2w[j][i][0] = h_Phi_px_w[j][i];
				h_Phi2w[j][i][1] = h_Phi_pz_w[j][i];
				h_Phi2w[j][i][2] = h_Phi_px2_w[j][i];
				h_Phi2w[j][i][3] = h_Phi_pxz_w[j][i];
				h_Phi2w[j][i][4] = h_Phi_pz2_w[j][i];

				//Perturbed Parallel Current

				h_A_dJpB[j][i] = pow(B[i][j], -2.0) * pow(J[i][j], -1.0) * ((B_px2[i][j] * gconxx[i][j] + B_px[i][j] * gconxx_px[i][j] + 2.0 * B_pxy[i][j] * gconxy[i][j] + B_px[i][j] * gconxy_py[i][j] + B_py2[i][j] * gconyy[i][j] + B_py[i][j] * (gconxy_px[i][j] + gconyy_py[i][j])) * J[i][j] + B_px[i][j] * gconxx[i][j] * J_px[i][j] + B_py[i][j] * gconxy[i][j] * J_px[i][j] + B_px[i][j] * gconxy[i][j] * J_py[i][j] + B_py[i][j] * gconyy[i][j] * J_py[i][j]);
				h_A_px_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * pow(J[i][j], -1.0) * ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] + gconxx[i][j] * J_px[i][j] + gconxy[i][j] * J_py[i][j]);
				h_A_pz_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * pow(J[i][j], -1.0) * ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] + gconxz[i][j] * J_px[i][j] + gconyz[i][j] * J_py[i][j]);
				h_A_px2_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * gconxx[i][j];
				h_A_pxz_dJpB[j][i] = (-2.0) * pow(B[i][j], -1.0) * gconxz[i][j];
				h_A_pz2_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * gconzz[i][j];
				h_A2dJpB[j][i][0] = h_A_dJpB[j][i];
				h_A2dJpB[j][i][1] = h_A_px_dJpB[j][i];
				h_A2dJpB[j][i][2] = h_A_pz_dJpB[j][i];
				h_A2dJpB[j][i][3] = h_A_px2_dJpB[j][i];
				h_A2dJpB[j][i][4] = h_A_pxz_dJpB[j][i];
				h_A2dJpB[j][i][5] = h_A_pz2_dJpB[j][i];

				//Parallel Resistive

				h_A_resistive[j][i] = h_A_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A) + 1.0;
				h_A_px_resistive[j][i] = h_A_px_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_pz_resistive[j][i] = h_A_pz_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_px2_resistive[j][i] = h_A_px2_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_pxz_resistive[j][i] = h_A_pxz_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_pz2_resistive[j][i] = h_A_pz2_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);

				//Perpendicular Dissipation

				h_F_perp2[j][i] = 1.0;
				h_F_px_perp2[j][i] = (1.0) * pow(J[i][j], -1.0) * ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] + gconxx[i][j] * J_px[i][j] + gconxy[i][j] * J_py[i][j]);
				h_F_pz_perp2[j][i] = (1.0) * pow(J[i][j], -1.0) * ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] + gconxz[i][j] * J_px[i][j] + gconyz[i][j] * J_py[i][j]);
				h_F_px2_perp2[j][i] = (1.0) * gconxx[i][j];
				h_F_pxz_perp2[j][i] = (2.0) * gconxz[i][j];
				h_F_pz2_perp2[j][i] = (1.0) * gconzz[i][j];
				h_F2perp2[j][i][0] = h_F_px_perp2[j][i];
				h_F2perp2[j][i][1] = h_F_pz_perp2[j][i];
				h_F2perp2[j][i][2] = h_F_px2_perp2[j][i];
				h_F2perp2[j][i][3] = h_F_pxz_perp2[j][i];
				h_F2perp2[j][i][4] = h_F_pz2_perp2[j][i];

				//Perturbed Parallel Vector Potential

				h_Phi_A[j][i] = 0.0;
				h_Phi_px_A[j][i] = 0.0;
				h_Phi_py_A[j][i] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j];
				h_Phi_pz_A[j][i] = 0.0;

				h_dNe_A[j][i] = 0.0;
				h_dNe_px_A[j][i] = 0.0;
				h_dNe_py_A[j][i] = (1.0 / 2.0) * pow(B[i][j], -1.0) * Bny[i][j] * pow(Ne[i][j], -1.0) * Te[i][j] * pow(NormQE, -1.0);
				h_dNe_pz_A[j][i] = 0.0;

				h_A_A[j][i] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);
				h_A_px_A[j][i] = 0.0;
				h_A_py_A[j][i] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);
				h_A_pz_A[j][i] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);

				h_APhidNe2A[j][i][0] = h_A_A[j][i];
				h_APhidNe2A[j][i][1] = h_A_py_A[j][i];
				h_APhidNe2A[j][i][2] = h_A_pz_A[j][i];
				h_APhidNe2A[j][i][3] = h_Phi_py_A[j][i];
				h_APhidNe2A[j][i][4] = h_dNe_py_A[j][i];

				//Perturbed Density

				h_Phi_dNe[j][i] = 0.0;
				h_Phi_px_dNe[j][i] = (-1.0) * pow(B[i][j], -3.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-2.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0) * Ne[i][j];
				h_Phi_py_dNe[j][i] = pow(B[i][j], -3.0) * pow(J[i][j], -1.0) * (B[i][j] * Bny_px[i][j] * gcovyz[i][j] * Ne[i][j] + Bny[i][j] * (B[i][j] * gcovyz_px[i][j] * Ne[i][j] + gcovyz[i][j] * ((-2.0) * B_px[i][j] * Ne[i][j] + B[i][j] * Ne_px[i][j])));
				h_Phi_pz_dNe[j][i] = pow(B[i][j], -3.0) * pow(J[i][j], -1.0) * (B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j]) * Ne[i][j] + (-1.0) * Bny[i][j] * (2.0 * B_py[i][j] * gcovxy[i][j] * Ne[i][j] + (-2.0) * B_px[i][j] * gcovyy[i][j] * Ne[i][j] + B[i][j] * (((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]) * Ne[i][j] + gcovyy[i][j] * Ne_px[i][j])));

				h_dPe_dNe[j][i] = 0.0;
				h_dPe_px_dNe[j][i] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-2.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_dPe_py_dNe[j][i] = (1.0 / 2.0) * pow(B[i][j], -3.0) * ((-1.0) * B[i][j] * Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * (2.0 * B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_dPe_pz_dNe[j][i] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j]) + Bny[i][j] * (2.0 * B_py[i][j] * gcovxy[i][j] + (-2.0) * B_px[i][j] * gcovyy[i][j] + B[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) * pow(J[i][j], -1.0) * pow(NormQE, -1.0);

				h_dJpB_dNe[j][i] = 0.0;
				h_dJpB_px_dNe[j][i] = 0.0;
				h_dJpB_py_dNe[j][i] = Bny[i][j] * pow(NormQE, -1.0);
				h_dJpB_pz_dNe[j][i] = 0.0;

				h_A_dNe[j][i] = pow(B[i][j], -2.0) * pow(J[i][j], -1.0) * (B[i][j] * gcovyz[i][j] * (Bny_py[i][j] * JpB_px[i][j] + (-1.0) * Bny_px[i][j] * JpB_py[i][j]) + Bny[i][j] * (gcovyz[i][j] * ((-1.0) * B_py[i][j] * JpB_px[i][j] + B_px[i][j] * JpB_py[i][j]) + B[i][j] * (gcovyz_py[i][j] * JpB_px[i][j] + (-1.0) * gcovyz_px[i][j] * JpB_py[i][j]))) * pow(NormQE, -1.0);
				h_A_px_dNe[j][i] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * JpB_py[i][j] * pow(NormQE, -1.0);
				h_A_py_dNe[j][i] = pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * JpB_px[i][j] * pow(NormQE, -1.0);
				h_A_pz_dNe[j][i] = pow(B[i][j], -1.0) * Bny[i][j] * pow(J[i][j], -1.0) * ((-1.0) * gcovyy[i][j] * JpB_px[i][j] + gcovxy[i][j] * JpB_py[i][j]) * pow(NormQE, -1.0);

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

				//Perturbed Temperature

				h_Phi_dTe[j][i] = 0.0;
				h_Phi_px_dTe[j][i] = 0.0;
				h_Phi_py_dTe[j][i] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * Te_px[i][j];
				h_Phi_pz_dTe[j][i] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * Te_px[i][j];

				h_dTe_dTe[j][i] = 0.0;
				h_dTe_px_dTe[j][i] = 0.0;
				h_dTe_py_dTe[j][i] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);
				h_dTe_pz_dTe[j][i] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);

				h_dNe_dTe[j][i] = 0.0;
				h_dNe_px_dTe[j][i] = 0.0;
				h_dNe_py_dTe[j][i] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Te[i][j] * Te_px[i][j] * pow(NormQE, -1.0);
				h_dNe_pz_dTe[j][i] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Te[i][j] * Te_px[i][j] * pow(NormQE, -1.0);

				h_PhidTedNe2dTe[j][i][0] = h_Phi_py_dTe[j][i];
				h_PhidTedNe2dTe[j][i][1] = h_Phi_pz_dTe[j][i];
				h_PhidTedNe2dTe[j][i][2] = h_dTe_py_dTe[j][i];
				h_PhidTedNe2dTe[j][i][3] = h_dTe_pz_dTe[j][i];
				h_PhidTedNe2dTe[j][i][4] = h_dNe_py_dTe[j][i];
				h_PhidTedNe2dTe[j][i][5] = h_dNe_pz_dTe[j][i];

				//Equilibrium Density and Temperature

				h_Ne0[j][i] = Ne[i][j];
				h_Te0[j][i] = Te[i][j];
				h_Ne0_px[j][i] = Ne_px[i][j];
				h_Te0_px[j][i] = Te_px[i][j];
				h_Pe0_px[j][i] = Pe_px[i][j];

			}
		}

		/*-----------------------Nonlinear-----------------------*/

		for (int i = 0; i < gridNx; i++) {
			for (int j = 0; j < gridNy; j++) {

				//Vorticity

				h_wPhi_w[j][i][0] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_wPhi_w[j][i][1] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_wPhi_w[j][i][2] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_wPhi_w[j][i][3] = pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);
				h_wPhi_w[j][i][4] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_wPhi_w[j][i][5] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);

				h_AdJpB_w[j][i][0] = pow(B[i][j], -2.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][1] = pow(B[i][j], -2.0) * (Bny[i][j] * B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][2] = pow(B[i][j], -2.0) * (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) + B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] + Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][3] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][4] = pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][5] = pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][6] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][7] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_AdJpB_w[j][i][8] = pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);

				//Perturbed Parallel Vector Potential

				h_PhiA_A[j][i][0] = pow(B[i][j], -3.0) * (Bny[i][j] * B_py[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * (Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][1] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][2] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][3] = pow(B[i][j], -3.0) * (B[i][j] * Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][4] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][5] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][6] = pow(B[i][j], -3.0) * (Bny[i][j] * ((-1.0) * B_py[i][j] * gcovxy[i][j] + B_px[i][j] * gcovyy[i][j]) + B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j] + Bny[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][7] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][8] = pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);

				h_NeA_A[j][i][0] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][1] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][2] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][3] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (Bny[i][j] * B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][4] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][5] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][6] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) + B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] + Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][7] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][8] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];

				//Perturbed Density

				h_AdJpB_dNe[j][i][0] = pow(B[i][j], -2.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][1] = pow(B[i][j], -2.0) * (Bny[i][j] * B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][2] = pow(B[i][j], -2.0) * (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) + B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] + Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][3] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][4] = pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][5] = pow(B[i][j], -1.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][6] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][7] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0);
				h_AdJpB_dNe[j][i][8] = pow(B[i][j], -1.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0);

				h_dNePhi_dNe[j][i][0] = pow(B[i][j], -3.0) * ((-1.0) * B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * (2.0 * B_py[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][1] = pow(B[i][j], -3.0) * (B[i][j] * Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * ((-2.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][2] = pow(B[i][j], -3.0) * (B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j]) + Bny[i][j] * ((-2.0) * B_py[i][j] * gcovxy[i][j] + 2.0 * B_px[i][j] * gcovyy[i][j] + B[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][3] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][4] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][5] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][6] = pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][7] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_dNePhi_dNe[j][i][8] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);

				//Perturbed Temperature

				h_PhiTe_dTe[j][i][0] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_PhiTe_dTe[j][i][1] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_PhiTe_dTe[j][i][2] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_PhiTe_dTe[j][i][3] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);
				h_PhiTe_dTe[j][i][4] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_PhiTe_dTe[j][i][5] = pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);

				h_PhiTeA_dTe[j][i][0] = pow(B[i][j], -4.0) * (gcovyz[i][j] * (B[i][j] * Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) + gcovxz[i][j] * (Bny[i][j] * B_py[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * (Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_py[i][j])) + (Bny[i][j] * ((-1.0) * B_py[i][j] * gcovxy[i][j] + B_px[i][j] * gcovyy[i][j]) + B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j] + Bny[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) * gcovzz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][1] = pow(B[i][j], -3.0) * Bny[i][j] * (pow(gcovyz[i][j], 2.0) + (-1.0) * gcovyy[i][j] * gcovzz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][2] = pow(B[i][j], -3.0) * Bny[i][j] * ((-1.0) * gcovxz[i][j] * gcovyz[i][j] + gcovxy[i][j] * gcovzz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][3] = pow(B[i][j], -3.0) * Bny[i][j] * (gcovxz[i][j] * gcovyy[i][j] + (-1.0) * gcovxy[i][j] * gcovyz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][4] = pow(B[i][j], -3.0) * Bny[i][j] * (((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]) * gcovyz[i][j] + (-1.0) * gcovyy[i][j] * gcovyz_px[i][j] + gcovxy[i][j] * gcovyz_py[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][5] = pow(B[i][j], -4.0) * (Bny[i][j] * (B_px[i][j] * pow(gcovyz[i][j], 2.0) + (-1.0) * B[i][j] * gcovyz[i][j] * gcovyz_px[i][j] + gcovxz[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j]) + B_py[i][j] * gcovxy[i][j] * gcovzz[i][j] + (-1.0) * B[i][j] * gcovxy_py[i][j] * gcovzz[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j] * gcovzz[i][j] + B[i][j] * gcovyy_px[i][j] * gcovzz[i][j]) + B[i][j] * (Bny_py[i][j] * gcovxz[i][j] * gcovyz[i][j] + (-1.0) * Bny_py[i][j] * gcovxy[i][j] * gcovzz[i][j] + Bny_px[i][j] * ((-1.0) * pow(gcovyz[i][j], 2.0) + gcovyy[i][j] * gcovzz[i][j]))) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][6] = pow(B[i][j], -3.0) * Bny[i][j] * ((-1.0) * pow(gcovyz[i][j], 2.0) + gcovyy[i][j] * gcovzz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][7] = pow(B[i][j], -3.0) * Bny[i][j] * (gcovxz[i][j] * gcovyz[i][j] + (-1.0) * gcovxy[i][j] * gcovzz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][8] = pow(B[i][j], -3.0) * Bny[i][j] * ((-1.0) * gcovxz[i][j] * gcovyy[i][j] + gcovxy[i][j] * gcovyz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][9] = pow(B[i][j], -4.0) * (B[i][j] * ((-1.0) * Bny_px[i][j] * gcovxz[i][j] * gcovyy[i][j] + (-1.0) * Bny_py[i][j] * gcovxx[i][j] * gcovyz[i][j] + gcovxy[i][j] * (Bny_py[i][j] * gcovxz[i][j] + Bny_px[i][j] * gcovyz[i][j])) + Bny[i][j] * (gcovxz[i][j] * (B_px[i][j] * gcovyy[i][j] + B[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j])) + (-1.0) * gcovxy[i][j] * (B_py[i][j] * gcovxz[i][j] + B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_px[i][j]) + gcovxx[i][j] * (B_py[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_py[i][j]))) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][10] = pow(B[i][j], -3.0) * Bny[i][j] * ((-1.0) * gcovxz[i][j] * gcovyy[i][j] + gcovxy[i][j] * gcovyz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][11] = pow(B[i][j], -3.0) * Bny[i][j] * (gcovxy[i][j] * gcovxz[i][j] + (-1.0) * gcovxx[i][j] * gcovyz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][12] = pow(B[i][j], -3.0) * Bny[i][j] * ((-1.0) * pow(gcovxy[i][j], 2.0) + gcovxx[i][j] * gcovyy[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][13] = pow(B[i][j], -3.0) * Bny[i][j] * ((gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]) * gcovyz[i][j] + gcovyy[i][j] * gcovyz_px[i][j] + (-1.0) * gcovxy[i][j] * gcovyz_py[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][14] = pow(B[i][j], -4.0) * (B[i][j] * (Bny_px[i][j] * gcovxz[i][j] * gcovyy[i][j] + Bny_py[i][j] * gcovxx[i][j] * gcovyz[i][j] + (-1.0) * gcovxy[i][j] * (Bny_py[i][j] * gcovxz[i][j] + Bny_px[i][j] * gcovyz[i][j])) + Bny[i][j] * ((-1.0) * gcovxz[i][j] * (B_px[i][j] * gcovyy[i][j] + B[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j])) + gcovxy[i][j] * (B_py[i][j] * gcovxz[i][j] + B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * gcovyz_px[i][j]) + gcovxx[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j]))) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][15] = pow(B[i][j], -3.0) * Bny[i][j] * (gcovxz[i][j] * gcovyy[i][j] + (-1.0) * gcovxy[i][j] * gcovyz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][16] = pow(B[i][j], -3.0) * Bny[i][j] * ((-1.0) * gcovxy[i][j] * gcovxz[i][j] + gcovxx[i][j] * gcovyz[i][j]) * pow(J[i][j], -2.0);
				h_PhiTeA_dTe[j][i][17] = pow(B[i][j], -3.0) * Bny[i][j] * (pow(gcovxy[i][j], 2.0) + (-1.0) * gcovxx[i][j] * gcovyy[i][j]) * pow(J[i][j], -2.0);

			}
		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void compressStaggeredCoefficient() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Compress staggered coefficient in shifted metric coordinate." << RESET << std::endl;

		/*-------------------------Linear-------------------------*/

		for (int i = 0; i < gridNx; i++) {
			for (int j = 0; j < gridNy; j++) {

				//Perturbed Parallel Current

				h_A_dJpB[j][i] = pow(B[i][j], -2.0) * pow(J[i][j], -1.0) * ((B_px2[i][j] * gconxx[i][j] + B_px[i][j] * gconxx_px[i][j] + 2.0 * B_pxy[i][j] * gconxy[i][j] + B_px[i][j] * gconxy_py[i][j] + B_py2[i][j] * gconyy[i][j] + B_py[i][j] * (gconxy_px[i][j] + gconyy_py[i][j])) * J[i][j] + B_px[i][j] * gconxx[i][j] * J_px[i][j] + B_py[i][j] * gconxy[i][j] * J_px[i][j] + B_px[i][j] * gconxy[i][j] * J_py[i][j] + B_py[i][j] * gconyy[i][j] * J_py[i][j]);
				h_A_px_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * pow(J[i][j], -1.0) * ((gconxx_px[i][j] + gconxy_py[i][j]) * J[i][j] + gconxx[i][j] * J_px[i][j] + gconxy[i][j] * J_py[i][j]);
				h_A_pz_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * pow(J[i][j], -1.0) * ((gconxz_px[i][j] + gconyz_py[i][j]) * J[i][j] + gconxz[i][j] * J_px[i][j] + gconyz[i][j] * J_py[i][j]);
				h_A_px2_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * gconxx[i][j];
				h_A_pxz_dJpB[j][i] = (-2.0) * pow(B[i][j], -1.0) * gconxz[i][j];
				h_A_pz2_dJpB[j][i] = (-1.0) * pow(B[i][j], -1.0) * gconzz[i][j];
				h_A2dJpB[j][i][0] = h_A_dJpB[j][i];
				h_A2dJpB[j][i][1] = h_A_px_dJpB[j][i];
				h_A2dJpB[j][i][2] = h_A_pz_dJpB[j][i];
				h_A2dJpB[j][i][3] = h_A_px2_dJpB[j][i];
				h_A2dJpB[j][i][4] = h_A_pxz_dJpB[j][i];
				h_A2dJpB[j][i][5] = h_A_pz2_dJpB[j][i];

				//Parallel Resistive

				h_A_resistive[j][i] = h_A_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A) + 1.0;
				h_A_px_resistive[j][i] = h_A_px_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_pz_resistive[j][i] = h_A_pz_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_px2_resistive[j][i] = h_A_px2_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_pxz_resistive[j][i] = h_A_pxz_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);
				h_A_pz2_resistive[j][i] = h_A_pz2_dJpB[j][i] * B[i][j] * dt * std::get<1>(nablaPerp2A);

				//Perturbed Parallel Vector Potential

				h_Phi_A[j][i] = 0.0;
				h_Phi_px_A[j][i] = 0.0;
				h_Phi_py_A[j][i] = (-1.0) * pow(B[i][j], -1.0) * Bny[i][j];
				h_Phi_pz_A[j][i] = 0.0;

				h_dNe_A[j][i] = 0.0;
				h_dNe_px_A[j][i] = 0.0;
				h_dNe_py_A[j][i] = (1.0 / 2.0) * pow(B[i][j], -1.0) * Bny[i][j] * pow(Ne[i][j], -1.0) * Te[i][j] * pow(NormQE, -1.0);
				h_dNe_pz_A[j][i] = 0.0;

				h_A_A[j][i] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);
				h_A_px_A[j][i] = 0.0;
				h_A_py_A[j][i] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);
				h_A_pz_A[j][i] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(Ne[i][j], -1.0) * Ne_px[i][j] * Te[i][j] * pow(NormQE, -1.0);

				h_APhidNe2A[j][i][0] = h_A_A[j][i];
				h_APhidNe2A[j][i][1] = h_A_py_A[j][i];
				h_APhidNe2A[j][i][2] = h_A_pz_A[j][i];
				h_APhidNe2A[j][i][3] = h_Phi_py_A[j][i];
				h_APhidNe2A[j][i][4] = h_dNe_py_A[j][i];

			}
		}

		/*-----------------------Nonlinear-----------------------*/

		for (int i = 0; i < gridNx; i++) {
			for (int j = 0; j < gridNy; j++) {

				//Perturbed Parallel Vector Potential

				h_PhiA_A[j][i][0] = pow(B[i][j], -3.0) * (Bny[i][j] * B_py[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * (Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][1] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][2] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][3] = pow(B[i][j], -3.0) * (B[i][j] * Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_px[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][4] = pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][5] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][6] = pow(B[i][j], -3.0) * (Bny[i][j] * ((-1.0) * B_py[i][j] * gcovxy[i][j] + B_px[i][j] * gcovyy[i][j]) + B[i][j] * (Bny_py[i][j] * gcovxy[i][j] + (-1.0) * Bny_px[i][j] * gcovyy[i][j] + Bny[i][j] * (gcovxy_py[i][j] + (-1.0) * gcovyy_px[i][j]))) * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][7] = (-1.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0);
				h_PhiA_A[j][i][8] = pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0);

				h_NeA_A[j][i][0] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (B[i][j] * Bny_py[i][j] * gcovyz[i][j] + Bny[i][j] * ((-1.0) * B_py[i][j] * gcovyz[i][j] + B[i][j] * gcovyz_py[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][1] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][2] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][3] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (Bny[i][j] * B_px[i][j] * gcovyz[i][j] + (-1.0) * B[i][j] * (Bny_px[i][j] * gcovyz[i][j] + Bny[i][j] * gcovyz_px[i][j])) * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][4] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyz[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][5] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][6] = (1.0 / 2.0) * pow(B[i][j], -3.0) * (Bny[i][j] * (B_py[i][j] * gcovxy[i][j] + (-1.0) * B_px[i][j] * gcovyy[i][j]) + B[i][j] * ((-1.0) * Bny_py[i][j] * gcovxy[i][j] + Bny_px[i][j] * gcovyy[i][j] + Bny[i][j] * ((-1.0) * gcovxy_py[i][j] + gcovyy_px[i][j]))) * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][7] = (1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovyy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];
				h_NeA_A[j][i][8] = (-1.0 / 2.0) * pow(B[i][j], -2.0) * Bny[i][j] * gcovxy[i][j] * pow(J[i][j], -1.0) * pow(NormQE, -1.0) * pow(Ne[i][j], -1.0) * Te[i][j];

			}
		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Copy staggered coefficient to device." << RESET << std::endl;

		Allocator H2DAllocator;

		for (int i = 0; i < devNums; i++) {

			CUDACHECK(cudaSetDevice(localId * devNums + i));

			/*-------------------------Linear-------------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				//Perturbed Parallel Vector Potential in Perturbed Parallel Current
				d_A_dJpB[i], h_A_dJpB[0], d_A_px_dJpB[i], h_A_px_dJpB[0], d_A_pz_dJpB[i], h_A_pz_dJpB[0],
				d_A_px2_dJpB[i], h_A_px2_dJpB[0], d_A_pxz_dJpB[i], h_A_pxz_dJpB[0], d_A_pz2_dJpB[i], h_A_pz2_dJpB[0]);

			H2DAllocator.hostToDevice(devNy * gridNx, 0, (hostId * hostNy + i * devNy) * gridNx,
				//Perturbed Electric Potential in Perturbed Parallel Vector Potential
				d_Phi_A[i], h_Phi_A[0], d_Phi_px_A[i], h_Phi_px_A[0], d_Phi_py_A[i], h_Phi_py_A[0], d_Phi_pz_A[i], h_Phi_pz_A[0],
				//Perturbed Density in Perturbed Parallel Vector Potential
				d_dNe_A[i], h_dNe_A[0], d_dNe_px_A[i], h_dNe_px_A[0], d_dNe_py_A[i], h_dNe_py_A[0], d_dNe_pz_A[i], h_dNe_pz_A[0],
				//Perturbed Parallel Vector Potential in Perturbed Parallel Vector Potential
				d_A_A[i], h_A_A[0], d_A_px_A[i], h_A_px_A[0], d_A_py_A[i], h_A_py_A[0], d_A_pz_A[i], h_A_pz_A[0]);

			H2DAllocator.hostToDevice(devNy * gridNx * 6, 0, (hostId * hostNy + i * devNy) * gridNx * 6, d_A2dJpB[i], h_A2dJpB[0][0]);
			H2DAllocator.hostToDevice(devNy * gridNx * 5, 0, (hostId * hostNy + i * devNy) * gridNx * 5, d_APhidNe2A[i], h_APhidNe2A[0][0]);

			/*-----------------------Nonlinear-----------------------*/

			H2DAllocator.hostToDevice(devNy * gridNx * 9, 0, (hostId * hostNy + i * devNy) * gridNx * 9, d_PhiA_A[i], h_PhiA_A[0][0], d_NeA_A[i], h_NeA_A[0][0]);

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	template<int matrixType, typename innerDirichlet, typename outerDirichlet>
	void computeSparseMatrix() {

		if (hostId == 0) {

			if constexpr (matrixType == 0)
				std::cout << BOLDYELLOW << "Start: Compute Poisson COO." << RESET << std::endl;
			else if constexpr (matrixType == 1)
				std::cout << BOLDYELLOW << "Start: Compute Resistive COO." << RESET << std::endl;
			else if constexpr (matrixType == 2)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2Phi COO." << RESET << std::endl;
			else if constexpr (matrixType == 3)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dNe COO." << RESET << std::endl;
			else if constexpr (matrixType == 4)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dTe COO." << RESET << std::endl;
			else if constexpr (matrixType == 5)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dPi COO." << RESET << std::endl;
			else if constexpr (matrixType == 6)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dPa COO." << RESET << std::endl;
			else if constexpr (matrixType == 7)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dPb COO." << RESET << std::endl;

		}

		matrix_i.clear();
		matrix_j.clear();
		matrix_v.clear();

		for (int j = hostId * hostNy; j < (hostId + 1) * hostNy; j++) {

			for (int i = 0; i < gridNx; i++) {
				for (int k = 0; k < gridNz; k++) {

					int row_index = i * gridNz + k;
					std::vector<double> coes(25, 0.0);
					double cc, cx, cz, cx2, cxz, cz2;

					if constexpr (matrixType == 0) {

						cc = 0.0;
						cx = h_Phi_px_w[j][i];
						cz = h_Phi_pz_w[j][i];
						cx2 = h_Phi_px2_w[j][i];
						cxz = h_Phi_pxz_w[j][i];
						cz2 = h_Phi_pz2_w[j][i];

					}
					else if constexpr (matrixType == 1) {

						cc = h_A_resistive[j][i];
						cx = h_A_px_resistive[j][i];
						cz = h_A_pz_resistive[j][i];
						cx2 = h_A_px2_resistive[j][i];
						cxz = h_A_pxz_resistive[j][i];
						cz2 = h_A_pz2_resistive[j][i];

					}
					else if constexpr (matrixType == 2) {

						cc = h_F_perp2[j][i];
						cx = -h_F_px_perp2[j][i] * dt * std::get<1>(nablaPerp2w);
						cz = -h_F_pz_perp2[j][i] * dt * std::get<1>(nablaPerp2w);
						cx2 = -h_F_px2_perp2[j][i] * dt * std::get<1>(nablaPerp2w);
						cxz = -h_F_pxz_perp2[j][i] * dt * std::get<1>(nablaPerp2w);
						cz2 = -h_F_pz2_perp2[j][i] * dt * std::get<1>(nablaPerp2w);

					}
					else if constexpr (matrixType == 3) {

						cc = h_F_perp2[j][i];
						cx = -h_F_px_perp2[j][i] * dt * std::get<1>(nablaPerp2dNe);
						cz = -h_F_pz_perp2[j][i] * dt * std::get<1>(nablaPerp2dNe);
						cx2 = -h_F_px2_perp2[j][i] * dt * std::get<1>(nablaPerp2dNe);
						cxz = -h_F_pxz_perp2[j][i] * dt * std::get<1>(nablaPerp2dNe);
						cz2 = -h_F_pz2_perp2[j][i] * dt * std::get<1>(nablaPerp2dNe);

					}
					else if constexpr (matrixType == 4) {

						cc = h_F_perp2[j][i];
						cx = -h_F_px_perp2[j][i] * dt * std::get<1>(nablaPerp2dTe);
						cz = -h_F_pz_perp2[j][i] * dt * std::get<1>(nablaPerp2dTe);
						cx2 = -h_F_px2_perp2[j][i] * dt * std::get<1>(nablaPerp2dTe);
						cxz = -h_F_pxz_perp2[j][i] * dt * std::get<1>(nablaPerp2dTe);
						cz2 = -h_F_pz2_perp2[j][i] * dt * std::get<1>(nablaPerp2dTe);

					}
					else if constexpr (matrixType == 5) {

						cc = h_F_perp2[j][i];
						cx = -h_F_px_perp2[j][i] * dt * std::get<1>(nablaPerp2dPi);
						cz = -h_F_pz_perp2[j][i] * dt * std::get<1>(nablaPerp2dPi);
						cx2 = -h_F_px2_perp2[j][i] * dt * std::get<1>(nablaPerp2dPi);
						cxz = -h_F_pxz_perp2[j][i] * dt * std::get<1>(nablaPerp2dPi);
						cz2 = -h_F_pz2_perp2[j][i] * dt * std::get<1>(nablaPerp2dPi);

					}
					else if constexpr (matrixType == 6) {

						cc = h_F_perp2[j][i];
						cx = -h_F_px_perp2[j][i] * dt * std::get<1>(nablaPerp2dPa);
						cz = -h_F_pz_perp2[j][i] * dt * std::get<1>(nablaPerp2dPa);
						cx2 = -h_F_px2_perp2[j][i] * dt * std::get<1>(nablaPerp2dPa);
						cxz = -h_F_pxz_perp2[j][i] * dt * std::get<1>(nablaPerp2dPa);
						cz2 = -h_F_pz2_perp2[j][i] * dt * std::get<1>(nablaPerp2dPa);

					}
					else if constexpr (matrixType == 7) {

						cc = h_F_perp2[j][i];
						cx = -h_F_px_perp2[j][i] * dt * std::get<1>(nablaPerp2dPb);
						cz = -h_F_pz_perp2[j][i] * dt * std::get<1>(nablaPerp2dPb);
						cx2 = -h_F_px2_perp2[j][i] * dt * std::get<1>(nablaPerp2dPb);
						cxz = -h_F_pxz_perp2[j][i] * dt * std::get<1>(nablaPerp2dPb);
						cz2 = -h_F_pz2_perp2[j][i] * dt * std::get<1>(nablaPerp2dPb);

					}

					if (i == 0) {

						if constexpr (std::is_same_v<innerDirichlet, std::integral_constant<bool, true>>) {

							int i_index = i;
							int k_index = k;

							int col_index = i_index * gridNz + k_index;

							matrix_i.emplace_back(row_index);
							matrix_j.emplace_back(col_index);
							matrix_v.emplace_back(1);

						}
						else {

							coes[0] = -25.0 / 12.0 / gridDx;
							coes[1] = 48.0 / 12.0 / gridDx;
							coes[2] = -36.0 / 12.0 / gridDx;
							coes[3] = 16.0 / 12.0 / gridDx;
							coes[4] = -3.0 / 12.0 / gridDx;

							int coes_index = -1;

							for (int i_offset = 0; i_offset <= 4; i_offset++) {

								coes_index++;
								int i_index = i + i_offset;
								int k_index = k;

								int col_index = i_index * gridNz + k_index;

								matrix_i.emplace_back(row_index);
								matrix_j.emplace_back(col_index);
								matrix_v.emplace_back(coes[coes_index]);

							}

						}

					}
					else if (i == gridNx - 1) {

						if constexpr (std::is_same_v<outerDirichlet, std::integral_constant<bool, true>>) {

							int i_index = i;
							int k_index = k;

							int col_index = i_index * gridNz + k_index;

							matrix_i.emplace_back(row_index);
							matrix_j.emplace_back(col_index);
							matrix_v.emplace_back(1);

						}
						else {

							coes[0] = 3.0 / 12.0 / gridDx;
							coes[1] = -16.0 / 12.0 / gridDx;
							coes[2] = 36.0 / 12.0 / gridDx;
							coes[3] = -48.0 / 12.0 / gridDx;
							coes[4] = 25.0 / 12.0 / gridDx;

							int coes_index = -1;

							for (int i_offset = -4; i_offset <= 0; i_offset++) {

								coes_index++;
								int i_index = i + i_offset;
								int k_index = k;

								int col_index = i_index * gridNz + k_index;

								matrix_i.emplace_back(row_index);
								matrix_j.emplace_back(col_index);
								matrix_v.emplace_back(coes[coes_index]);

							}

						}

					}
					else if (i == 1) {

						coes[0] = (-1.0 / 48.0) * cxz / gridDx / gridDz;
						coes[1] = (1.0 / 6.0) * cxz / gridDx / gridDz;
						coes[2] = (1.0 / 12.0) * pow(gridDx, -2.0) * (11.0 * cx2 + (-3.0) * cx * gridDx);
						coes[3] = (-1.0 / 6.0) * cxz / gridDx / gridDz;
						coes[4] = (1.0 / 48.0) * cxz / gridDx / gridDz;

						coes[5] = (-1.0 / 72.0) / gridDx * pow(gridDz, -2.0) * (6.0 * cz2 * gridDx + 5.0 * cxz * gridDz + (-6.0) * cz * gridDx * gridDz);
						coes[6] = (1.0 / 9.0) / gridDx * pow(gridDz, -2.0) * (12.0 * cz2 * gridDx + 5.0 * cxz * gridDz + (-6.0) * cz * gridDx * gridDz);
						coes[7] = (-5.0 / 6.0) * pow(gridDx, -2.0) * (2.0 * cx2 + cx * gridDx) + (-5.0 / 2.0) * cz2 * pow(gridDz, -2.0) + cc;
						coes[8] = (1.0 / 9.0) / gridDx * pow(gridDz, -2.0) * (12.0 * cz2 * gridDx + (-5.0) * cxz * gridDz + 6.0 * cz * gridDx * gridDz);
						coes[9] = (-1.0 / 72.0) / gridDx * pow(gridDz, -2.0) * (6.0 * cz2 * gridDx + (-5.0) * cxz * gridDz + 6.0 * cz * gridDx * gridDz);

						coes[10] = (1.0 / 8.0) * cxz / gridDx / gridDz;
						coes[11] = (-1.0) * cxz / gridDx / gridDz;
						coes[12] = (1.0 / 2.0) * pow(gridDx, -2.0) * (cx2 + 3.0 * cx * gridDx);
						coes[13] = cxz / gridDx / gridDz;
						coes[14] = (-1.0 / 8.0) * cxz / gridDx / gridDz;

						coes[15] = (-1.0 / 24.0) * cxz / gridDx / gridDz;
						coes[16] = (1.0 / 3.0) * cxz / gridDx / gridDz;
						coes[17] = (1.0 / 6.0) * pow(gridDx, -2.0) * (2.0 * cx2 + (-3.0) * cx * gridDx);
						coes[18] = (-1.0 / 3.0) * cxz / gridDx / gridDz;
						coes[19] = (1.0 / 24.0) * cxz / gridDx / gridDz;

						coes[20] = (1.0 / 144.0) * cxz / gridDx / gridDz;
						coes[21] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[22] = (-1.0 / 12.0) * pow(gridDx, -2.0) * (cx2 + (-1.0) * cx * gridDx);
						coes[23] = (1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[24] = (-1.0 / 144.0) * cxz / gridDx / gridDz;

						int coes_index = -1;

						for (int i_offset = -1; i_offset <= 3; i_offset++) {
							for (int k_offset = -2; k_offset <= 2; k_offset++) {

								coes_index++;
								int i_index = i + i_offset;
								int k_index = k + k_offset;

								if (k_index < 0)
									k_index = k_index + gridNz;
								else if (k_index > gridNz - 1)
									k_index = k_index - gridNz;

								int col_index = i_index * gridNz + k_index;

								matrix_i.emplace_back(row_index);
								matrix_j.emplace_back(col_index);
								matrix_v.emplace_back(coes[coes_index]);

							}
						}

					}
					else if (i == gridNx - 2) {

						coes[0] = (-1.0 / 144.0) * cxz / gridDx / gridDz;
						coes[1] = (1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[2] = (-1.0 / 12.0) * pow(gridDx, -2.0) * (cx2 + cx * gridDx);
						coes[3] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[4] = (1.0 / 144.0) * cxz / gridDx / gridDz;

						coes[5] = (1.0 / 24.0) * cxz / gridDx / gridDz;
						coes[6] = (-1.0 / 3.0) * cxz / gridDx / gridDz;
						coes[7] = (1.0 / 6.0) * pow(gridDx, -2.0) * (2.0 * cx2 + 3.0 * cx * gridDx);
						coes[8] = (1.0 / 3.0) * cxz / gridDx / gridDz;
						coes[9] = (-1.0 / 24.0) * cxz / gridDx / gridDz;

						coes[10] = (-1.0 / 8.0) * cxz / gridDx / gridDz;
						coes[11] = cxz / gridDx / gridDz;
						coes[12] = (1.0 / 2.0) * pow(gridDx, -2.0) * (cx2 + (-3.0) * cx * gridDx);
						coes[13] = (-1.0) * cxz / gridDx / gridDz;
						coes[14] = (1.0 / 8.0) * cxz / gridDx / gridDz;

						coes[15] = (1.0 / 72.0) / gridDx * pow(gridDz, -2.0) * ((-6.0) * cz2 * gridDx + 5.0 * cxz * gridDz + 6.0 * cz * gridDx * gridDz);
						coes[16] = (1.0 / 9.0) / gridDx * pow(gridDz, -2.0) * (12.0 * cz2 * gridDx + (-5.0) * cxz * gridDz + (-6.0) * cz * gridDx * gridDz);
						coes[17] = (5.0 / 6.0) * ((-2.0) * cx2 * pow(gridDx, -2.0) + cx / gridDx + (-3.0) * cz2 * pow(gridDz, -2.0)) + cc;
						coes[18] = (1.0 / 9.0) / gridDx * pow(gridDz, -2.0) * (12.0 * cz2 * gridDx + 5.0 * cxz * gridDz + 6.0 * cz * gridDx * gridDz);
						coes[19] = (-1.0 / 72.0) / gridDx * pow(gridDz, -2.0) * (6.0 * cz2 * gridDx + 5.0 * cxz * gridDz + 6.0 * cz * gridDx * gridDz);

						coes[20] = (1.0 / 48.0) * cxz / gridDx / gridDz;
						coes[21] = (-1.0 / 6.0) * cxz / gridDx / gridDz;
						coes[22] = (1.0 / 12.0) * pow(gridDx, -2.0) * (11.0 * cx2 + 3.0 * cx * gridDx);
						coes[23] = (1.0 / 6.0) * cxz / gridDx / gridDz;
						coes[24] = (-1.0 / 48.0) * cxz / gridDx / gridDz;

						int coes_index = -1;

						for (int i_offset = -3; i_offset <= 1; i_offset++) {
							for (int k_offset = -2; k_offset <= 2; k_offset++) {

								coes_index++;
								int i_index = i + i_offset;
								int k_index = k + k_offset;

								if (k_index < 0)
									k_index = k_index + gridNz;
								else if (k_index > gridNz - 1)
									k_index = k_index - gridNz;

								int col_index = i_index * gridNz + k_index;

								matrix_i.emplace_back(row_index);
								matrix_j.emplace_back(col_index);
								matrix_v.emplace_back(coes[coes_index]);

							}
						}

					}
					else {

						coes[0] = (1.0 / 144.0) * cxz / gridDx / gridDz;
						coes[1] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[2] = (-1.0 / 12.0) * pow(gridDx, -2.0) * (cx2 + (-1.0) * cx * gridDx);
						coes[3] = (1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[4] = (-1.0 / 144.0) * cxz / gridDx / gridDz;

						coes[5] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[6] = (4.0 / 9.0) * cxz / gridDx / gridDz;
						coes[7] = (1.0 / 3.0) * pow(gridDx, -2.0) * (4.0 * cx2 + (-2.0) * cx * gridDx);
						coes[8] = (-4.0 / 9.0) * cxz / gridDx / gridDz;
						coes[9] = (1.0 / 18.0) * cxz / gridDx / gridDz;

						coes[10] = (-1.0 / 12.0) * pow(gridDz, -2.0) * (cz2 + (-1.0) * cz * gridDz);
						coes[11] = (1.0 / 3.0) * pow(gridDz, -2.0) * (4.0 * cz2 + (-2.0) * cz * gridDz);
						coes[12] = (-5.0 / 2.0) * cx2 * pow(gridDx, -2.0) + (-5.0 / 2.0) * cz2 * pow(gridDz, -2.0) + cc;
						coes[13] = (2.0 / 3.0) * pow(gridDz, -2.0) * (2.0 * cz2 + cz * gridDz);
						coes[14] = (-1.0 / 12.0) * pow(gridDz, -2.0) * (cz2 + cz * gridDz);

						coes[15] = (1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[16] = (-4.0 / 9.0) * cxz / gridDx / gridDz;
						coes[17] = (2.0 / 3.0) * pow(gridDx, -2.0) * (2.0 * cx2 + cx * gridDx);
						coes[18] = (4.0 / 9.0) * cxz / gridDx / gridDz;
						coes[19] = (-1.0 / 18.0) * cxz / gridDx / gridDz;

						coes[20] = (-1.0 / 144.0) * cxz / gridDx / gridDz;
						coes[21] = (1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[22] = (-1.0 / 12.0) * pow(gridDx, -2.0) * (cx2 + cx * gridDx);
						coes[23] = (-1.0 / 18.0) * cxz / gridDx / gridDz;
						coes[24] = (1.0 / 144.0) * cxz / gridDx / gridDz;

						int coes_index = -1;

						for (int i_offset = -2; i_offset <= 2; i_offset++) {
							for (int k_offset = -2; k_offset <= 2; k_offset++) {

								coes_index++;
								int i_index = i + i_offset;
								int k_index = k + k_offset;

								if (k_index < 0)
									k_index = k_index + gridNz;
								else if (k_index > gridNz - 1)
									k_index = k_index - gridNz;

								int col_index = i_index * gridNz + k_index;

								matrix_i.emplace_back(row_index);
								matrix_j.emplace_back(col_index);
								matrix_v.emplace_back(coes[coes_index]);

							}
						}

					}

				}
			}

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

		int nnz = matrix_i.size() / hostNy;

		if (hostId == 0) {

			if constexpr (matrixType == 0)
				std::cout << BOLDYELLOW << "Start: Compute Poisson CSR." << std::endl;
			else if constexpr (matrixType == 1)
				std::cout << BOLDYELLOW << "Start: Compute Resistive CSR." << std::endl;
			else if constexpr (matrixType == 2)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2Phi CSR." << std::endl;
			else if constexpr (matrixType == 3)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dNe CSR." << std::endl;
			else if constexpr (matrixType == 4)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dTe CSR." << std::endl;
			else if constexpr (matrixType == 5)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dPi CSR." << std::endl;
			else if constexpr (matrixType == 6)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dPa CSR." << std::endl;
			else if constexpr (matrixType == 7)
				std::cout << BOLDYELLOW << "Start: Compute NablaPerp2dPb CSR." << std::endl;

		}

		if (hostId == 0) {
			std::cout << "nnz * hostNy: " << matrix_i.size() << "." << std::endl;
			std::cout << "nnz: " << nnz << "." << std::endl;
		}

		Allocator DeviceAllocator;
		int** d_matrixCooR;
		dataType** d_matrixTempV;

		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_matrixCooR);
		DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_matrixTempV);

		for (int i = 0; i < devNums; i++) {

			CUDACHECK(cudaSetDevice(localId * devNums + i));
			CUDACHECK(cudaMemcpy(d_matrixCooR[i], matrix_i.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
			CUDACHECK(cudaMemcpy(d_matrixTempV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

		}

		if constexpr (matrixType == 0) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_laplacianCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_laplacianCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_laplacianCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_laplacianCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_laplacianCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

		}
		else if constexpr (matrixType == 1) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_resistiveCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_resistiveCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_resistiveCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_resistiveCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_resistiveCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

		}
		else if constexpr (matrixType == 2) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_wCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_wCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_wCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_wCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_wCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

		}
		else if constexpr (matrixType == 3) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_dNeCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dNeCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dNeCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_dNeCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_dNeCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

		}
		else if constexpr (matrixType == 4) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_dTeCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dTeCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dTeCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_dTeCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_dTeCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

		}
		else if constexpr (matrixType == 5) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_dPiCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dPiCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dPiCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_dPiCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_dPiCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

		}
		else if constexpr (matrixType == 6) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_dPaCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dPaCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dPaCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_dPaCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_dPaCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

		}
		else if constexpr (matrixType == 7) {

			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * (gridNxz + 1), d_dPbCsrR);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dPbCsrC);
			DeviceAllocator.allocateDeviceArrays(localId, devNums, devNy * nnz, d_dPbCsrV);

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				CUDACHECK(cudaMemcpy(d_dPbCsrC[i], matrix_j.data() + i * devNy * nnz, sizeof(int) * devNy * nnz, cudaMemcpyHostToDevice));
				CUDACHECK(cudaMemcpy(d_dPbCsrV[i], matrix_v.data() + i * devNy * nnz, sizeof(dataType) * devNy * nnz, cudaMemcpyHostToDevice));

			}

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

				if constexpr (matrixType == 0) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_laplacianCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_laplacianCsrR[i] + j * (gridNxz + 1), d_laplacianCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_laplacianCsrR[i] + j * (gridNxz + 1), d_laplacianCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_laplacianCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_laplacianCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}
				else if constexpr (matrixType == 1) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_resistiveCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_resistiveCsrR[i] + j * (gridNxz + 1), d_resistiveCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_resistiveCsrR[i] + j * (gridNxz + 1), d_resistiveCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_resistiveCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_resistiveCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}
				else if constexpr (matrixType == 2) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_wCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_wCsrR[i] + j * (gridNxz + 1), d_wCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_wCsrR[i] + j * (gridNxz + 1), d_wCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_wCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_wCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}
				else if constexpr (matrixType == 3) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_dNeCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_dNeCsrR[i] + j * (gridNxz + 1), d_dNeCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_dNeCsrR[i] + j * (gridNxz + 1), d_dNeCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dNeCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dNeCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}
				else if constexpr (matrixType == 4) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_dTeCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_dTeCsrR[i] + j * (gridNxz + 1), d_dTeCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_dTeCsrR[i] + j * (gridNxz + 1), d_dTeCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dTeCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dTeCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}
				else if constexpr (matrixType == 5) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_dPiCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_dPiCsrR[i] + j * (gridNxz + 1), d_dPiCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_dPiCsrR[i] + j * (gridNxz + 1), d_dPiCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dPiCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dPiCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}
				else if constexpr (matrixType == 6) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_dPaCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_dPaCsrR[i] + j * (gridNxz + 1), d_dPaCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_dPaCsrR[i] + j * (gridNxz + 1), d_dPaCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dPaCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dPaCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}
				else if constexpr (matrixType == 7) {

					CUSPARSECHECK(cusparseXcoo2csr(cusparseHandles[i], d_matrixCooR[i] + j * nnz,
						nnz, gridNxz, d_dPbCsrR[i] + j * (gridNxz + 1), CUSPARSE_INDEX_BASE_ZERO));
					CUSPARSECHECK(cusparseXcsrsort_bufferSizeExt(cusparseHandles[i], gridNxz, gridNxz, nnz,
						d_dPbCsrR[i] + j * (gridNxz + 1), d_dPbCsrC[i] + j * nnz, &pBufferSize[i][j]));
					CUDACHECK(cudaMalloc(&pBuffer[i][j], sizeof(char) * pBufferSize[i][j]));
					CUDACHECK(cudaMalloc((void**)&permutation[i][j], sizeof(int) * nnz));
					CUSPARSECHECK(cusparseCreateIdentityPermutation(cusparseHandles[i], nnz, permutation[i][j]));
					CUSPARSECHECK(cusparseCreateMatDescr(&descr[i][j]));
					CUSPARSECHECK(cusparseXcsrsort(cusparseHandles[i], gridNxz, gridNxz, nnz, descr[i][j],
						d_dPbCsrR[i] + j * (gridNxz + 1), d_dPbCsrC[i] + j * nnz, permutation[i][j], pBuffer[i][j]));
					if constexpr (std::is_same_v<dataType, double>) {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dPbCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_64F));
					}
					else {
						CUSPARSECHECK(cusparseCreateSpVec(&sortedSpVec[i][j], nnz, nnz, permutation[i][j],
							d_dPbCsrV[i] + j * nnz, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
						CUSPARSECHECK(cusparseCreateDnVec(&unsortedDnVec[i][j], nnz, d_matrixTempV[i] + j * nnz, CUDA_R_32F));
					}
					CUSPARSECHECK(cusparseGather(cusparseHandles[i], unsortedDnVec[i][j], sortedSpVec[i][j]));

				}

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
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

		if (hostId == 0) {

			if constexpr (matrixType == 0)
				std::cout << BOLDYELLOW << "Start: Factorize Poisson CSR." << RESET << std::endl;
			else if constexpr (matrixType == 1)
				std::cout << BOLDYELLOW << "Start: Factorize Resistive CSR." << RESET << std::endl;
			else if constexpr (matrixType == 2)
				std::cout << BOLDYELLOW << "Start: Factorize NablaPerp2Phi CSR." << RESET << std::endl;
			else if constexpr (matrixType == 3)
				std::cout << BOLDYELLOW << "Start: Factorize NablaPerp2dNe CSR." << RESET << std::endl;
			else if constexpr (matrixType == 4)
				std::cout << BOLDYELLOW << "Start: Factorize NablaPerp2dTe CSR." << RESET << std::endl;
			else if constexpr (matrixType == 5)
				std::cout << BOLDYELLOW << "Start: Factorize NablaPerp2dPi CSR." << RESET << std::endl;
			else if constexpr (matrixType == 6)
				std::cout << BOLDYELLOW << "Start: Factorize NablaPerp2dPa CSR." << RESET << std::endl;
			else if constexpr (matrixType == 7)
				std::cout << BOLDYELLOW << "Start: Factorize NablaPerp2dPb CSR." << RESET << std::endl;

		}

		if constexpr (matrixType == 0) {

			cudaStreams.resize(devNums, std::vector<cudaStream_t>(devNy));
			cudssHandles.resize(devNums, std::vector<cudssHandle_t>(devNy));

			laplacianConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			laplacianDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			laplacianAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			laplacianXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			laplacianBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}
		else if constexpr (matrixType == 1) {

			resistiveConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			resistiveDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			resistiveAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			resistiveXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			resistiveBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}
		else if constexpr (matrixType == 2) {

			wConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			wDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			wAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			wXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			wBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}
		else if constexpr (matrixType == 3) {

			dNeConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			dNeDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			dNeAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dNeXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dNeBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}
		else if constexpr (matrixType == 4) {

			dTeConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			dTeDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			dTeAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dTeXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dTeBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}
		else if constexpr (matrixType == 5) {

			dPiConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			dPiDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			dPiAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dPiXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dPiBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}
		else if constexpr (matrixType == 6) {

			dPaConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			dPaDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			dPaAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dPaXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dPaBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}
		else if constexpr (matrixType == 7) {

			dPbConfigs.resize(devNums, std::vector<cudssConfig_t>(devNy));
			dPbDatas.resize(devNums, std::vector<cudssData_t>(devNy));
			dPbAs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dPbXs.resize(devNums, std::vector<cudssMatrix_t>(devNy));
			dPbBs.resize(devNums, std::vector<cudssMatrix_t>(devNy));

		}

		for (int i = 0; i < devNums; i++) {

			CUDACHECK(cudaSetDevice(localId * devNums + i));

			for (int j = 0; j < devNy; j++) {

				if constexpr (matrixType == 0) {

					CUDACHECK(cudaStreamCreate(&cudaStreams[i][j]));
					CUDSSCHECK(cudssCreate(&cudssHandles[i][j]));
					CUDSSCHECK(cudssSetStream(cudssHandles[i][j], cudaStreams[i][j]));

					CUDSSCHECK(cudssConfigCreate(&laplacianConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &laplacianDatas[i][j]));

				}
				else if constexpr (matrixType == 1) {

					CUDSSCHECK(cudssConfigCreate(&resistiveConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &resistiveDatas[i][j]));

				}
				else if constexpr (matrixType == 2) {

					CUDSSCHECK(cudssConfigCreate(&wConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &wDatas[i][j]));

				}
				else if constexpr (matrixType == 3) {

					CUDSSCHECK(cudssConfigCreate(&dNeConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &dNeDatas[i][j]));

				}
				else if constexpr (matrixType == 4) {

					CUDSSCHECK(cudssConfigCreate(&dTeConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &dTeDatas[i][j]));

				}
				else if constexpr (matrixType == 5) {

					CUDSSCHECK(cudssConfigCreate(&dPiConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &dPiDatas[i][j]));

				}
				else if constexpr (matrixType == 6) {

					CUDSSCHECK(cudssConfigCreate(&dPaConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &dPaDatas[i][j]));

				}
				else if constexpr (matrixType == 7) {

					CUDSSCHECK(cudssConfigCreate(&dPbConfigs[i][j]));
					CUDSSCHECK(cudssDataCreate(cudssHandles[i][j], &dPbDatas[i][j]));

				}

				cudssMatrixType_t mtype = CUDSS_MTYPE_GENERAL;
				cudssMatrixViewType_t mview = CUDSS_MVIEW_FULL;
				cudssIndexBase_t base = CUDSS_BASE_ZERO;

				int64_t nrows = gridNxz;
				int64_t ncols = gridNxz;
				int64_t ld = gridNxz;
				int64_t nrhs = 1;

				if constexpr (matrixType == 0) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&laplacianAs[i][j], nrows, ncols, nnz, d_laplacianCsrR[i] + j * (gridNxz + 1), NULL,
							d_laplacianCsrC[i] + j * nnz, d_laplacianCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&laplacianXs[i][j], nrows, nrhs, ld, d_Phi_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&laplacianBs[i][j], nrows, nrhs, ld, d_w_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&laplacianAs[i][j], nrows, ncols, nnz, d_laplacianCsrR[i] + j * (gridNxz + 1), NULL,
							d_laplacianCsrC[i] + j * nnz, d_laplacianCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&laplacianXs[i][j], nrows, nrhs, ld, d_Phi_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&laplacianBs[i][j], nrows, nrhs, ld, d_w_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						laplacianConfigs[i][j], laplacianDatas[i][j], laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						laplacianConfigs[i][j], laplacianDatas[i][j], laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]));

				}
				else if constexpr (matrixType == 1) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&resistiveAs[i][j], nrows, ncols, nnz, d_resistiveCsrR[i] + j * (gridNxz + 1), NULL,
							d_resistiveCsrC[i] + j * nnz, d_resistiveCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&resistiveXs[i][j], nrows, nrhs, ld, d_A_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&resistiveBs[i][j], nrows, nrhs, ld, d_A_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&resistiveAs[i][j], nrows, ncols, nnz, d_resistiveCsrR[i] + j * (gridNxz + 1), NULL,
							d_resistiveCsrC[i] + j * nnz, d_resistiveCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&resistiveXs[i][j], nrows, nrhs, ld, d_A_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&resistiveBs[i][j], nrows, nrhs, ld, d_A_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						resistiveConfigs[i][j], resistiveDatas[i][j], resistiveAs[i][j], resistiveXs[i][j], resistiveBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						resistiveConfigs[i][j], resistiveDatas[i][j], resistiveAs[i][j], resistiveXs[i][j], resistiveBs[i][j]));

				}
				else if constexpr (matrixType == 2) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&wAs[i][j], nrows, ncols, nnz, d_wCsrR[i] + j * (gridNxz + 1), NULL,
							d_wCsrC[i] + j * nnz, d_wCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&wXs[i][j], nrows, nrhs, ld, d_Phi_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&wBs[i][j], nrows, nrhs, ld, d_Phi_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&wAs[i][j], nrows, ncols, nnz, d_wCsrR[i] + j * (gridNxz + 1), NULL,
							d_wCsrC[i] + j * nnz, d_wCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&wXs[i][j], nrows, nrhs, ld, d_Phi_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&wBs[i][j], nrows, nrhs, ld, d_Phi_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						wConfigs[i][j], wDatas[i][j], wAs[i][j], wXs[i][j], wBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						wConfigs[i][j], wDatas[i][j], wAs[i][j], wXs[i][j], wBs[i][j]));

				}
				else if constexpr (matrixType == 3) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&dNeAs[i][j], nrows, ncols, nnz, d_dNeCsrR[i] + j * (gridNxz + 1), NULL,
							d_dNeCsrC[i] + j * nnz, d_dNeCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dNeXs[i][j], nrows, nrhs, ld, d_dNe_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dNeBs[i][j], nrows, nrhs, ld, d_dNe_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&dNeAs[i][j], nrows, ncols, nnz, d_dNeCsrR[i] + j * (gridNxz + 1), NULL,
							d_dNeCsrC[i] + j * nnz, d_dNeCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dNeXs[i][j], nrows, nrhs, ld, d_dNe_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dNeBs[i][j], nrows, nrhs, ld, d_dNe_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						dNeConfigs[i][j], dNeDatas[i][j], dNeAs[i][j], dNeXs[i][j], dNeBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						dNeConfigs[i][j], dNeDatas[i][j], dNeAs[i][j], dNeXs[i][j], dNeBs[i][j]));

				}
				else if constexpr (matrixType == 4) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&dTeAs[i][j], nrows, ncols, nnz, d_dTeCsrR[i] + j * (gridNxz + 1), NULL,
							d_dTeCsrC[i] + j * nnz, d_dTeCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dTeXs[i][j], nrows, nrhs, ld, d_dTe_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dTeBs[i][j], nrows, nrhs, ld, d_dTe_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&dTeAs[i][j], nrows, ncols, nnz, d_dTeCsrR[i] + j * (gridNxz + 1), NULL,
							d_dTeCsrC[i] + j * nnz, d_dTeCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dTeXs[i][j], nrows, nrhs, ld, d_dTe_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dTeBs[i][j], nrows, nrhs, ld, d_dTe_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						dTeConfigs[i][j], dTeDatas[i][j], dTeAs[i][j], dTeXs[i][j], dTeBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						dTeConfigs[i][j], dTeDatas[i][j], dTeAs[i][j], dTeXs[i][j], dTeBs[i][j]));

				}
				else if constexpr (matrixType == 5) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&dPiAs[i][j], nrows, ncols, nnz, d_dPiCsrR[i] + j * (gridNxz + 1), NULL,
							d_dPiCsrC[i] + j * nnz, d_dPiCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dPiXs[i][j], nrows, nrhs, ld, d_dPi_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dPiBs[i][j], nrows, nrhs, ld, d_dPi_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&dPiAs[i][j], nrows, ncols, nnz, d_dPiCsrR[i] + j * (gridNxz + 1), NULL,
							d_dPiCsrC[i] + j * nnz, d_dPiCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dPiXs[i][j], nrows, nrhs, ld, d_dPi_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dPiBs[i][j], nrows, nrhs, ld, d_dPi_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						dPiConfigs[i][j], dPiDatas[i][j], dPiAs[i][j], dPiXs[i][j], dPiBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						dPiConfigs[i][j], dPiDatas[i][j], dPiAs[i][j], dPiXs[i][j], dPiBs[i][j]));

				}
				else if constexpr (matrixType == 6) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&dPaAs[i][j], nrows, ncols, nnz, d_dPaCsrR[i] + j * (gridNxz + 1), NULL,
							d_dPaCsrC[i] + j * nnz, d_dPaCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dPaXs[i][j], nrows, nrhs, ld, d_dPa_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dPaBs[i][j], nrows, nrhs, ld, d_dPa_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&dPaAs[i][j], nrows, ncols, nnz, d_dPaCsrR[i] + j * (gridNxz + 1), NULL,
							d_dPaCsrC[i] + j * nnz, d_dPaCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dPaXs[i][j], nrows, nrhs, ld, d_dPa_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dPaBs[i][j], nrows, nrhs, ld, d_dPa_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						dPaConfigs[i][j], dPaDatas[i][j], dPaAs[i][j], dPaXs[i][j], dPaBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						dPaConfigs[i][j], dPaDatas[i][j], dPaAs[i][j], dPaXs[i][j], dPaBs[i][j]));

				}
				else if constexpr (matrixType == 7) {

					if constexpr (std::is_same_v<dataType, double>) {

						CUDSSCHECK(cudssMatrixCreateCsr(&dPbAs[i][j], nrows, ncols, nnz, d_dPbCsrR[i] + j * (gridNxz + 1), NULL,
							d_dPbCsrC[i] + j * nnz, d_dPbCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_64F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dPbXs[i][j], nrows, nrhs, ld, d_dPb_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dPbBs[i][j], nrows, nrhs, ld, d_dPb_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_64F,
							CUDSS_LAYOUT_COL_MAJOR));

					}
					else {

						CUDSSCHECK(cudssMatrixCreateCsr(&dPbAs[i][j], nrows, ncols, nnz, d_dPbCsrR[i] + j * (gridNxz + 1), NULL,
							d_dPbCsrC[i] + j * nnz, d_dPbCsrV[i] + j * nnz, CUDA_R_32I, CUDA_R_32F, mtype, mview, base));
						CUDSSCHECK(cudssMatrixCreateDn(&dPbXs[i][j], nrows, nrhs, ld, d_dPb_midr[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));
						CUDSSCHECK(cudssMatrixCreateDn(&dPbBs[i][j], nrows, nrhs, ld, d_dPb_midl[i] + (j + gridGhost) * gridNxz, CUDA_R_32F,
							CUDSS_LAYOUT_COL_MAJOR));

					}

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_ANALYSIS,
						dPbConfigs[i][j], dPbDatas[i][j], dPbAs[i][j], dPbXs[i][j], dPbBs[i][j]));

					CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_FACTORIZATION,
						dPbConfigs[i][j], dPbDatas[i][j], dPbAs[i][j], dPbXs[i][j], dPbBs[i][j]));

				}

			}

		}

		if (hostId == 0) {
			size_t avail, total, used;
			CUDACHECK(cudaSetDevice(localId * devNums));
			CUDACHECK(cudaMemGetInfo(&avail, &total));
			used = total - avail;
			std::cout << BOLDYELLOW << "Device memory used: " << (double)used / 1024 / 1024 / 1024 << " GB." << RESET << std::endl;
		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	template<int picType, int disType, int markerType>
	void loadParticles() {

		if (hostId == 0) {

			if constexpr (picType == 0) {

				std::cout << BOLDYELLOW << "Start: Load thermal ions";
				if constexpr (disType == 0)
					std::cout << BOLDYELLOW << " by Maxwell distribution." << RESET << std::endl;
				else if constexpr (disType == 1)
					std::cout << BOLDYELLOW << " by isotropic slowing-down distribution without erf function." << RESET << std::endl;
				else if constexpr (disType == 2)
					std::cout << BOLDYELLOW << " by isotropic slowing-down distribution with erf function." << RESET << std::endl;
				else if constexpr (disType == 3)
					std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution without erf function." << RESET << std::endl;
				else if constexpr (disType == 4)
					std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution with erf function." << RESET << std::endl;

			}
			else if constexpr (picType == 1) {

				std::cout << BOLDYELLOW << "Start: Load alpha particles";
				if constexpr (disType == 0)
					std::cout << BOLDYELLOW << " by Maxwell distribution." << RESET << std::endl;
				else if constexpr (disType == 1)
					std::cout << BOLDYELLOW << " by isotropic slowing-down distribution without erf function." << RESET << std::endl;
				else if constexpr (disType == 2)
					std::cout << BOLDYELLOW << " by isotropic slowing-down distribution with erf function." << RESET << std::endl;
				else if constexpr (disType == 3)
					std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution without erf function." << RESET << std::endl;
				else if constexpr (disType == 4)
					std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution with erf function." << RESET << std::endl;

			}
			else if constexpr (picType == 2) {

				std::cout << BOLDYELLOW << "Start: Load beam particles";
				if constexpr (disType == 0)
					std::cout << BOLDYELLOW << " by Maxwell distribution." << RESET << std::endl;
				else if constexpr (disType == 1)
					std::cout << BOLDYELLOW << " by isotropic slowing-down distribution without erf function." << RESET << std::endl;
				else if constexpr (disType == 2)
					std::cout << BOLDYELLOW << " by isotropic slowing-down distribution with erf function." << RESET << std::endl;
				else if constexpr (disType == 3)
					std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution without erf function." << RESET << std::endl;
				else if constexpr (disType == 4)
					std::cout << BOLDYELLOW << " by anisotropic slowing-down distribution with erf function." << RESET << std::endl;

			}

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

		}
		else if constexpr (picType == 1) {

			Mass = AlphaMass;
			Vmin = AlphaVmin;
			Vmax = AlphaVmax;
			Vb = AlphaVb;
			DeltaV = AlphaDeltaV;
			Lambda0 = AlphaLambda0;
			DeltaLambda2 = AlphaDeltaLambda2;

		}
		else if constexpr (picType == 2) {

			Mass = BeamMass;
			Vmin = BeamVmin;
			Vmax = BeamVmax;
			Vb = BeamVb;
			DeltaV = BeamDeltaV;
			Lambda0 = BeamLambda0;
			DeltaLambda2 = BeamDeltaLambda2;

		}

		Jmax = 0.0;
		Jvmax = pow(Vmax, 2.0);

		for (int i = 0; i < gridNx; i++) {
			for (int j = 0; j < gridNy; j++) {

				tempJ[i][j + 1] = J[i][j];
				tempB[i][j + 1] = B[i][j];

				if constexpr (picType == 0) {
					tempN[i][j + 1] = Ni[i][j] * 1.0e-19;
					tempT[i][j + 1] = Ti[i][j];
				}
				else if constexpr (picType == 1) {
					tempN[i][j + 1] = Na[i][j] * 1.0e-19;
					tempT[i][j + 1] = Ta[i][j];
				}
				else if constexpr (picType == 2) {
					tempN[i][j + 1] = Nb[i][j] * 1.0e-19;
					tempT[i][j + 1] = Tb[i][j];
				}

				if constexpr (markerType == 0) {
					if (tempJ[i][j + 1] * tempN[i][j + 1] > Jmax)
						Jmax = tempJ[i][j + 1] * tempN[i][j + 1];
				}
				else if constexpr (markerType == 1) {
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
			int i = floor(li);
			int j = floor(lj);
			double dx = li - i;
			double dy = lj - j;

			double coes[4] = {};
			double cx[4] = { 1.0, 1.0, 0.0, 0.0 };
			double sx[4] = { -1.0, -1.0, 1.0, 1.0 };
			double cy[4] = { 1.0, 0.0, 1.0, 0.0 };
			double sy[4] = { -1.0, 1.0, -1.0, 1.0 };

			double result = 0.0;

			coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
			coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
			coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
			coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

			result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] + field[i + 1][j + 1] * coes[3];

			return result;

		};

#pragma omp parallel for num_threads(devNums)
		for (int devId = 0; devId < devNums; devId++) {
			for (int picId = 0; picId < picDev / gridNz; picId++) {

				int i, j, k;
				double li, lj, lk;
				double J, B, v, Jv, N, T, Lambda;
				double x, y, z, vp, mu, pw;

				if constexpr (markerType == 0) {
					do {
						x = x0 + (x1 - x0) * xrand[devId]();
						y = y0 + (y1 - y0) * yrand[devId]();
						z = z0 + gridDz * zrand[devId]();
						N = interp2d(tempN, x, y);
						J = interp2d(tempJ, x, y);
					} while (Jrand[devId]() >= N * J / Jmax);
				}
				else if constexpr (markerType == 1) {
					do {
						x = x0 + (x1 - x0) * xrand[devId]();
						y = y0 + (y1 - y0) * yrand[devId]();
						z = z0 + gridDz * zrand[devId]();
						J = interp2d(tempJ, x, y);
					} while (Jrand[devId]() >= J / Jmax);
				}
				else if constexpr (markerType == 2) {
					x = x0 + (x1 - x0) * xrand[devId]();
					y = y0 + (y1 - y0) * yrand[devId]();
					z = z0 + gridDz * zrand[devId]();
					J = interp2d(tempJ, x, y);
				}

				do {
					v = Vmax * vrand[devId]();
					Jv = pow(v, 2.0);
				} while (v <= Vmin || Jvrand[devId]() >= Jv / Jvmax);

				B = interp2d(tempB, x, y);
				N = interp2d(tempN, x, y);
				T = interp2d(tempT, x, y);

				if constexpr (picType == 2)
					Lambda = vprand[devId]();
				else
					Lambda = 2.0 * vprand[devId]() - 1.0;
				vp = v * Lambda;
				mu = 0.5 * Mass * pow(v, 2.0) * (1.0 - pow(Lambda, 2.0)) / B;
				Lambda = mu / (0.5 * Mass * pow(v, 2.0));

				if constexpr (disType == 0)
					pw = N * pow(T, -1.5) * exp(-0.5 * Mass * pow(v, 2.0) * MP * pow(VA0, 2.0) / (T * KEV));
				else if constexpr (disType == 1)
					pw = N / (pow(v, 3.0) + pow(T, 3.0));
				else if constexpr (disType == 2)
					pw = N / (pow(v, 3.0) + pow(T, 3.0)) * (1.0 + erf((Vb - v) / DeltaV));
				else if constexpr (disType == 3)
					pw = N / (pow(v, 3.0) + pow(T, 3.0)) * exp(-pow(Lambda - Lambda0, 2.0) / DeltaLambda2);
				else if constexpr (disType == 4)
					pw = N / (pow(v, 3.0) + pow(T, 3.0)) * exp(-pow(Lambda - Lambda0, 2.0) / DeltaLambda2) * (1.0 + erf((Vb - v) / DeltaV));

				if constexpr (markerType == 0)
					pw /= N;
				else if constexpr (markerType == 1)
					pw *= 1.0;
				else if constexpr (markerType == 2)
					pw *= J;

				li = (x - x0PlusGhost) / gridDx;
				lj = (y - y0PlusGhost) / gridDy;
				lk = (z - z0PlusGhost) / gridDz;

				i = floor(li);
				j = floor(lj);
				k = floor(lk);

				if constexpr (picType == 0) {

					for (int repeatId = 0; repeatId < gridNz; repeatId++) {

						for (int varId = 0; varId < 7; varId++) {
							h_Ion_offsets[devId][varId + 1] = h_Ion_offsets[devId][varId] + picDev;
							h_Ion_keys[devId][picId + picDev / gridNz * repeatId + varId * picDev] = j * cellNxz + i * cellNz + k + repeatId;
						}

						h_Ion_values[devId][picId + picDev / gridNz * repeatId + 0 * picDev] = 0.999999 * x;
						h_Ion_values[devId][picId + picDev / gridNz * repeatId + 1 * picDev] = y;
						h_Ion_values[devId][picId + picDev / gridNz * repeatId + 2 * picDev] = z + repeatId * gridDz;
						h_Ion_values[devId][picId + picDev / gridNz * repeatId + 3 * picDev] = vp;
						h_Ion_values[devId][picId + picDev / gridNz * repeatId + 4 * picDev] = 0.0;
						h_Ion_values[devId][picId + picDev / gridNz * repeatId + 5 * picDev] = pw;
						h_Ion_values[devId][picId + picDev / gridNz * repeatId + 6 * picDev] = mu;

					}

				}
				else if constexpr (picType == 1) {

					for (int repeatId = 0; repeatId < gridNz; repeatId++) {

						for (int varId = 0; varId < 7; varId++) {
							h_Alpha_offsets[devId][varId + 1] = h_Alpha_offsets[devId][varId] + picDev;
							h_Alpha_keys[devId][picId + picDev / gridNz * repeatId + varId * picDev] = j * cellNxz + i * cellNz + k + repeatId;
						}

						h_Alpha_values[devId][picId + picDev / gridNz * repeatId + 0 * picDev] = 0.999999 * x;
						h_Alpha_values[devId][picId + picDev / gridNz * repeatId + 1 * picDev] = y;
						h_Alpha_values[devId][picId + picDev / gridNz * repeatId + 2 * picDev] = z + repeatId * gridDz;
						h_Alpha_values[devId][picId + picDev / gridNz * repeatId + 3 * picDev] = vp;
						h_Alpha_values[devId][picId + picDev / gridNz * repeatId + 4 * picDev] = 0.0;
						h_Alpha_values[devId][picId + picDev / gridNz * repeatId + 5 * picDev] = pw;
						h_Alpha_values[devId][picId + picDev / gridNz * repeatId + 6 * picDev] = mu;

					}

				}
				else if constexpr (picType == 2) {

					for (int repeatId = 0; repeatId < gridNz; repeatId++) {

						for (int varId = 0; varId < 7; varId++) {
							h_Beam_offsets[devId][varId + 1] = h_Beam_offsets[devId][varId] + picDev;
							h_Beam_keys[devId][picId + picDev / gridNz * repeatId + varId * picDev] = j * cellNxz + i * cellNz + k + repeatId;
						}

						h_Beam_values[devId][picId + picDev / gridNz * repeatId + 0 * picDev] = 0.999999 * x;
						h_Beam_values[devId][picId + picDev / gridNz * repeatId + 1 * picDev] = y;
						h_Beam_values[devId][picId + picDev / gridNz * repeatId + 2 * picDev] = z + repeatId * gridDz;
						h_Beam_values[devId][picId + picDev / gridNz * repeatId + 3 * picDev] = vp;
						h_Beam_values[devId][picId + picDev / gridNz * repeatId + 4 * picDev] = 0.0;
						h_Beam_values[devId][picId + picDev / gridNz * repeatId + 5 * picDev] = pw;
						h_Beam_values[devId][picId + picDev / gridNz * repeatId + 6 * picDev] = mu;

					}

				}

			}

		}
#pragma omp barrier

		std::ofstream output;
		std::string fileName;

		if constexpr (picType == 0) {

			fileName = "IonOffsets_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Ion_offsets[0]), sizeof(int) * devNums * 8);
			output.close();

			fileName = "IonKeys_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Ion_keys[0]), sizeof(int) * devNums * picDev * 7);
			output.close();

			fileName = "IonValues_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Ion_values[0]), sizeof(dataType) * devNums * picDev * 7);
			output.close();

		}
		else if constexpr (picType == 1) {

			fileName = "AlphaOffsets_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Alpha_offsets[0]), sizeof(int) * devNums * 8);
			output.close();

			fileName = "AlphaKeys_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Alpha_keys[0]), sizeof(int) * devNums * picDev * 7);
			output.close();

			fileName = "AlphaValues_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Alpha_values[0]), sizeof(dataType) * devNums * picDev * 7);
			output.close();

		}
		else if constexpr (picType == 2) {

			fileName = "BeamOffsets_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Beam_offsets[0]), sizeof(int) * devNums * 8);
			output.close();

			fileName = "BeamKeys_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Beam_keys[0]), sizeof(int) * devNums * picDev * 7);
			output.close();

			fileName = "BeamValues_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(h_Beam_values[0]), sizeof(dataType) * devNums * picDev * 7);
			output.close();

		}

		if (hostId == 0) {

			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;

		}

	}
	template<int picType>
	void loadParticles(std::string file1, std::string file2, std::string file3, std::string file4) {

		if (hostId == 0) {

			if constexpr (picType == 0)
				std::cout << BOLDYELLOW << "Start: Load thermal ions from existed files." << RESET << std::endl;
			else if constexpr (picType == 1)
				std::cout << BOLDYELLOW << "Start: Load alpha particles from existed files." << RESET << std::endl;
			else if constexpr (picType == 2)
				std::cout << BOLDYELLOW << "Start: Load beam particles from existed files." << RESET << std::endl;

		}

		std::ifstream input;

		if constexpr (picType == 0) {

			input.open(file1, std::ios::in | std::ios::binary);
			input.read((char*)(&IonConst), sizeof(dataType));
			input.close();

			input.open(file2, std::ios::in | std::ios::binary);
			input.read((char*)(h_Ion_offsets[0]), sizeof(int) * devNums * 8);
			input.close();

			input.open(file3, std::ios::in | std::ios::binary);
			input.read((char*)(h_Ion_keys[0]), sizeof(int) * devNums * picDev * 7);
			input.close();

			input.open(file4, std::ios::in | std::ios::binary);
			input.read((char*)(h_Ion_values[0]), sizeof(dataType) * devNums * picDev * 7);
			input.close();

		}
		else if constexpr (picType == 1) {

			input.open(file1, std::ios::in | std::ios::binary);
			input.read((char*)(&AlphaConst), sizeof(dataType));
			input.close();

			input.open(file2, std::ios::in | std::ios::binary);
			input.read((char*)(h_Alpha_offsets[0]), sizeof(int) * devNums * 8);
			input.close();

			input.open(file3, std::ios::in | std::ios::binary);
			input.read((char*)(h_Alpha_keys[0]), sizeof(int) * devNums * picDev * 7);
			input.close();

			input.open(file4, std::ios::in | std::ios::binary);
			input.read((char*)(h_Alpha_values[0]), sizeof(dataType) * devNums * picDev * 7);
			input.close();

		}
		else if constexpr (picType == 2) {

			input.open(file1, std::ios::in | std::ios::binary);
			input.read((char*)(&BeamConst), sizeof(dataType));
			input.close();

			input.open(file2, std::ios::in | std::ios::binary);
			input.read((char*)(h_Beam_offsets[0]), sizeof(int) * devNums * 8);
			input.close();

			input.open(file3, std::ios::in | std::ios::binary);
			input.read((char*)(h_Beam_keys[0]), sizeof(int) * devNums * picDev * 7);
			input.close();

			input.open(file4, std::ios::in | std::ios::binary);
			input.read((char*)(h_Beam_values[0]), sizeof(dataType) * devNums * picDev * 7);
			input.close();

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	void loadRandom() {

		if (hostId == 0)
			std::cout << BOLDYELLOW << "Start: Load random coordinates for slowing down operation." << RESET << std::endl;

		double Jmax = 0.0;
		std::vector<Rand01> xrand(devNums);
		std::vector<Rand01> yrand(devNums);
		std::vector<Rand01> zrand(devNums);
		std::vector<Rand01> Jrand(devNums);
		std::vector<Rand01> vprand(devNums);

		std::vector<std::vector<double>> tempJ(gridNx, std::vector<double>(gridNy + 2));

		for (int i = 0; i < gridNx; i++) {

			for (int j = 0; j < gridNy; j++) {
				tempJ[i][j + 1] = J[i][j];
				if (J[i][j] > Jmax)
					Jmax = J[i][j];
			}

			tempJ[i][0] = tempJ[i][gridNy];
			tempJ[i][gridNy + 1] = tempJ[i][1];

		}

		auto interp2d = [&](std::vector<std::vector<double>>& field, double x, double y) {

			double li = (x - x0) / gridDx;
			double lj = (y - y0 + 0.5 * gridDy) / gridDy;
			int i = floor(li);
			int j = floor(lj);
			double dx = li - i;
			double dy = lj - j;

			double coes[4] = {};
			double cx[4] = { 1.0, 1.0, 0.0, 0.0 };
			double sx[4] = { -1.0, -1.0, 1.0, 1.0 };
			double cy[4] = { 1.0, 0.0, 1.0, 0.0 };
			double sy[4] = { -1.0, 1.0, -1.0, 1.0 };

			double result = 0.0;

			coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
			coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
			coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
			coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

			result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] + field[i + 1][j + 1] * coes[3];

			return result;

		};

#pragma omp parallel for num_threads(devNums)
		for (int devId = 0; devId < devNums; devId++) {
			for (int picId = 0; picId < randMax; picId++) {

				int i, j, k;
				double li, lj, lk;
				double x, y, z, J;

				do {
					x = x0 + (x1 - x0) * xrand[devId]();
					y = y0 + (y1 - y0) * yrand[devId]();
					z = z0 + (z1 - z0) * zrand[devId]();
					J = interp2d(tempJ, x, y);
				} while (Jrand[devId]() >= J / Jmax);

				li = (x - x0PlusGhost) / gridDx;
				lj = (y - y0PlusGhost) / gridDy;
				lk = (z - z0PlusGhost) / gridDz;

				i = floor(li);
				j = floor(lj);
				k = floor(lk);

				h_rand_keys[devId][picId + 0 * randMax] = j * cellNxz + i * cellNz + k;
				h_rand_values[devId][picId + 0 * randMax] = 0.999999 * x;
				h_rand_values[devId][picId + 1 * randMax] = y;
				h_rand_values[devId][picId + 2 * randMax] = z;
				h_rand_values[devId][picId + 3 * randMax] = 2.0 * vprand[devId]() - 1.0;

			}
		}
#pragma omp barrier

		if (hostId == 0) {

			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;

		}

	}
	template<int picType>
	void computeEquilibriumPressure() {

		if (hostId == 0) {

			if constexpr (picType == 0)
				std::cout << BOLDYELLOW << "Start: Compute equilibrium pressure contributed by thermal ions." << RESET << std::endl;
			else if constexpr (picType == 1)
				std::cout << BOLDYELLOW << "Start: Compute equilibrium pressure contributed by alpha particles." << RESET << std::endl;
			else if constexpr (picType == 2)
				std::cout << BOLDYELLOW << "Start: Compute equilibrium pressure contributed by beam particles." << RESET << std::endl;

		}

		dataType coes[8];
		dataType cx[8] = { 1.0, 0.0, 1.0, 0.0,1.0, 0.0, 1.0, 0.0 };
		dataType sx[8] = { -1.0, 1.0, -1.0, 1.0,-1.0, 1.0, -1.0, 1.0 };
		dataType cy[8] = { 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0 };
		dataType sy[8] = { -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0 };
		dataType cz[8] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
		dataType sz[8] = { -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 };

		int i, j, k, tileId, cellId;
		dataType li, lj, lk;
		dataType x, y, z, vp, pw, mu;
		dataType dx, dy, dz, J, B, P;
		dataType*** localP;
		dataType*** globalP;

		Allocator HostAllocator;
		HostAllocator.allocateHostArrays(gridNyPlusGhost, gridNxPlusGhost, gridNzPlusGhost, localP, globalP);

		for (int devId = 0; devId < devNums; devId++) {
			for (int picId = 0; picId < picDev; picId++) {

				if constexpr (picType == 0) {
					cellId = h_Ion_keys[devId][picId];
					x = h_Ion_values[devId][picId + 0 * picDev];
					y = h_Ion_values[devId][picId + 1 * picDev];
					z = h_Ion_values[devId][picId + 2 * picDev];
					vp = h_Ion_values[devId][picId + 3 * picDev];
					pw = h_Ion_values[devId][picId + 5 * picDev];
					mu = h_Ion_values[devId][picId + 6 * picDev];
				}
				else if constexpr (picType == 1) {
					cellId = h_Alpha_keys[devId][picId];
					x = h_Alpha_values[devId][picId + 0 * picDev];
					y = h_Alpha_values[devId][picId + 1 * picDev];
					z = h_Alpha_values[devId][picId + 2 * picDev];
					vp = h_Alpha_values[devId][picId + 3 * picDev];
					pw = h_Alpha_values[devId][picId + 5 * picDev];
					mu = h_Alpha_values[devId][picId + 6 * picDev];
				}
				else if constexpr (picType == 2) {
					cellId = h_Beam_keys[devId][picId];
					x = h_Beam_values[devId][picId + 0 * picDev];
					y = h_Beam_values[devId][picId + 1 * picDev];
					z = h_Beam_values[devId][picId + 2 * picDev];
					vp = h_Beam_values[devId][picId + 3 * picDev];
					pw = h_Beam_values[devId][picId + 5 * picDev];
					mu = h_Beam_values[devId][picId + 6 * picDev];
				}

				tileId = cellId / cellNz;

				li = (x - x0PlusGhost) / gridDx;
				lj = (y - y0PlusGhost) / gridDy;
				lk = (z - z0PlusGhost) / gridDz;

				i = floor(li);
				j = floor(lj);
				k = floor(lk);

				dx = li - i;
				dy = lj - j;
				dz = lk - k;

				coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
				coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
				coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
				coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

				J = h_pic2d[tileId][0] * coes[0] + h_pic2d[tileId][1] * coes[1] + h_pic2d[tileId][2] * coes[2] + h_pic2d[tileId][3] * coes[3];
				B = h_pic2d[tileId][4] * coes[0] + h_pic2d[tileId][5] * coes[1] + h_pic2d[tileId][6] * coes[2] + h_pic2d[tileId][7] * coes[3];

				coes[4] = coes[0]; coes[5] = coes[1]; coes[6] = coes[2]; coes[7] = coes[3];
				coes[0] *= (cz[0] + sz[0] * dz); coes[1] *= (cz[1] + sz[1] * dz);
				coes[2] *= (cz[2] + sz[2] * dz); coes[3] *= (cz[3] + sz[3] * dz);
				coes[4] *= (cz[4] + sz[4] * dz); coes[5] *= (cz[5] + sz[5] * dz);
				coes[6] *= (cz[6] + sz[6] * dz); coes[7] *= (cz[7] + sz[7] * dz);

				if constexpr (picType == 0)
					P = IonMass * vp * vp * pw / J + mu * B * pw / J;
				else if constexpr (picType == 1)
					P = AlphaMass * vp * vp * pw / J + mu * B * pw / J;
				else if constexpr (picType == 2)
					P = BeamMass * vp * vp * pw / J + mu * B * pw / J;

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

		if constexpr (std::is_same_v<dataType, double>) {
			MPICHECK(MPI_Allreduce(localP[0][0], globalP[0][0],
				gridNyPlusGhost * gridNxPlusGhost * gridNzPlusGhost,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
		}
		else {
			MPICHECK(MPI_Allreduce(localP[0][0], globalP[0][0],
				gridNyPlusGhost * gridNxPlusGhost * gridNzPlusGhost,
				MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD));
		}

		for (int j = 0; j < gridNyPlusGhost; j++) {
			for (int k = 0; k < gridNzPlusGhost; k++) {
				globalP[j][0][k] *= 2.0;
				globalP[j][gridNxPlusGhost - 1][k] *= 2.0;
			}
			for (int i = 0; i < gridNxPlusGhost; i++) {
				globalP[j][i][gridGhost] += globalP[j][i][gridNz + gridGhost];
				globalP[j][i][gridNz + gridGhost - 1] += globalP[j][i][gridGhost - 1];
			}
		}

		for (int i = 0; i < gridNxPlusGhost; i++) {
			for (int k = 0; k < gridNzPlusGhost; k++) {
				globalP[gridGhost][i][k] += globalP[gridNy + gridGhost][i][k];
				globalP[gridNy + gridGhost - 1][i][k] += globalP[gridGhost - 1][i][k];
			}
		}

		dataType innerP = 0.0;

		for (int j = 0; j < gridNy; j++) {
			for (int k = 0; k < gridNz; k++) {
				innerP += globalP[j + gridGhost][0][k + gridGhost];
			}
		}

		innerP /= (2 * gridNy * gridNz);

		if constexpr (picType == 0)
			IonConst = IonBeta / innerP;
		else if constexpr (picType == 1)
			AlphaConst = AlphaBeta / innerP;
		else if constexpr (picType == 2)
			BeamConst = BeamBeta / innerP;

		HostAllocator.releaseHostArrays(localP, globalP);

		std::ofstream output;
		std::string fileName;

		if constexpr (picType == 0) {

			fileName = "IonConst_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(&IonConst), sizeof(dataType));
			output.close();

		}
		else if constexpr (picType == 1) {

			fileName = "AlphaConst_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(&AlphaConst), sizeof(dataType));
			output.close();

		}
		else if constexpr (picType == 2) {

			fileName = "BeamConst_" + std::to_string(hostId) + "_0" + ".bin";
			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(&BeamConst), sizeof(dataType));
			output.close();

		}

		if (hostId == 0) {

			if constexpr (picType == 0)
				std::cout << BOLDYELLOW << "IonConst for computing pressure: " << std::setprecision(10) << IonConst << "." << RESET << std::endl;
			else if constexpr (picType == 1)
				std::cout << BOLDYELLOW << "AlphaConst for computing pressure: " << std::setprecision(10) << AlphaConst << "." << RESET << std::endl;
			else if constexpr (picType == 2)
				std::cout << BOLDYELLOW << "BeamConst for computing pressure: " << std::setprecision(10) << BeamConst << "." << RESET << std::endl;

			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;

		}

	}

	template<int picType, int disType, int markerType, int gridE, int gridPphi, int gridLambda, int ppcPhase>
	void computePhaseSpaceF0() {

		if (hostId == 0) {

			if constexpr (picType == 0)
				std::cout << BOLDYELLOW << "Start: Compute phase space f0 of thermal ions." << RESET << std::endl;
			else if constexpr (picType == 1)
				std::cout << BOLDYELLOW << "Start: Compute phase space f0 of alpha particles." << RESET << std::endl;
			else if constexpr (picType == 2)
				std::cout << BOLDYELLOW << "Start: Compute phase space f0 of beam particles." << RESET << std::endl;

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

		}
		else if constexpr (picType == 1) {

			Mass = AlphaMass;
			Char = AlphaChar;
			Vmin = AlphaVmin;
			Vmax = AlphaVmax;
			Vb = AlphaVb;
			DeltaV = AlphaDeltaV;
			Lambda0 = AlphaLambda0;
			DeltaLambda2 = AlphaDeltaLambda2;

		}
		else if constexpr (picType == 2) {

			Mass = BeamMass;
			Char = BeamChar;
			Vmin = BeamVmin;
			Vmax = BeamVmax;
			Vb = BeamVb;
			DeltaV = BeamDeltaV;
			Lambda0 = BeamLambda0;
			DeltaLambda2 = BeamDeltaLambda2;

		}

		minE = 0.5 * Mass * pow(Vmin, 2.0);
		maxE = 0.5 * Mass * pow(Vmax, 2.0);

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

				tempPphi = cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) * SFAcovyz[i][j+gridGhost] / (q[i][j] * J[i][j] * B[i][j]) - Char * psip[i][j];
				if (tempPphi > maxPphi)
					maxPphi = tempPphi;

				tempPphi = -cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) * SFAcovyz[i][j+gridGhost] / (q[i][j] * J[i][j] * B[i][j]) - Char * psip[i][j];
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

		double Jmax, Jvmax;
		Jmax = 0.0;
		Jvmax = pow(Vmax, 2.0);

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
				}
				else if constexpr (picType == 1) {
					tempN[i][j + 1] = Na[i][j] * 1.0e-19;
					tempT[i][j + 1] = Ta[i][j];
				}
				else if constexpr (picType == 2) {
					tempN[i][j + 1] = Nb[i][j] * 1.0e-19;
					tempT[i][j + 1] = Tb[i][j];
				}

				if constexpr (markerType == 0) {
					if (tempJ[i][j + 1] * tempN[i][j + 1] > Jmax)
						Jmax = tempJ[i][j + 1] * tempN[i][j + 1];
				}
				else if constexpr (markerType == 1) {
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
			int i = floor(li);
			int j = floor(lj);
			double dx = li - i;
			double dy = lj - j;

			double coes[4] = {};
			double cx[4] = { 1.0, 1.0, 0.0, 0.0 };
			double sx[4] = { -1.0, -1.0, 1.0, 1.0 };
			double cy[4] = { 1.0, 0.0, 1.0, 0.0 };
			double sy[4] = { -1.0, 1.0, -1.0, 1.0 };

			double result = 0.0;

			coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
			coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
			coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
			coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

			result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] + field[i + 1][j + 1] * coes[3];

			return result;

		};

		Rand01 xrand;
		Rand01 yrand;
		Rand01 Jrand;

		Rand01 vrand;
		Rand01 vprand;
		Rand01 Jvrand;

		double cx[8] = { 1.0, 0.0, 1.0, 0.0,1.0, 0.0, 1.0, 0.0 };
		double sx[8] = { -1.0, 1.0, -1.0, 1.0,-1.0, 1.0, -1.0, 1.0 };
		double cy[8] = { 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0 };
		double sy[8] = { -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0 };
		double cz[8] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
		double sz[8] = { -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 };

		for (size_t picId = 0; picId < picPhase; picId++) {

			if (picId % (picPhase / 10) == 0)
				if (hostId == 0)
					std::cout << BOLDGREEN << 100 * picId / picPhase << "%" << RESET << std::endl;

			int i, j, k;
			double li, lj, lk;
			double dx, dy, dz;
			double q, psip, J, B, N, T, SFAcovyz;
			double v, Jv, E, Pphi, Lambda;
			double x, y, z, vp, mu, pw;
			double coes[8];

			if constexpr (markerType == 0) {
				do {
					x = x0 + (x1 - x0) * xrand();
					y = y0 + (y1 - y0) * yrand();
					N = interp2d(tempN, x, y);
					J = interp2d(tempJ, x, y);
				} while (Jrand() >= N * J / Jmax);
			}
			else if constexpr (markerType == 1) {
				do {
					x = x0 + (x1 - x0) * xrand();
					y = y0 + (y1 - y0) * yrand();
					J = interp2d(tempJ, x, y);
				} while (Jrand() >= J / Jmax);
			}
			else if constexpr (markerType == 2) {
				x = x0 + (x1 - x0) * xrand();
				y = y0 + (y1 - y0) * yrand();
				J = interp2d(tempJ, x, y);
			}

			do {
				v = Vmax * vrand();
				Jv = pow(v, 2.0);
			} while (v <= Vmin || Jvrand() >= Jv / Jvmax);

			q = interp2d(tempq, x, y);
			psip = interp2d(temppsip, x, y);
			SFAcovyz = interp2d(tempSFAcovyz, x, y);
			B = interp2d(tempB, x, y);
			N = interp2d(tempN, x, y);
			T = interp2d(tempT, x, y);

			if constexpr (picType == 2)
				Lambda = vprand();
			else
				Lambda = 2.0 * vprand() - 1.0;
			vp = v * Lambda;
			mu = 0.5 * Mass * pow(v, 2.0) * (1.0 - pow(Lambda, 2.0)) / B;

			E = 0.5 * Mass * pow(v, 2.0);
			Pphi = cm * Mass * vp * 2 * psitmax * drho * (RHO0 + x * drho) * SFAcovyz / (q * J * B) - Char * psip;
			Lambda = mu / E;

			if constexpr (disType == 0)
				pw = N * pow(T, -1.5) * exp(-0.5 * Mass * pow(v, 2.0) * MP * pow(VA0, 2.0) / (T * KEV));
			else if constexpr (disType == 1)
				pw = N / (pow(v, 3.0) + pow(T, 3.0));
			else if constexpr (disType == 2)
				pw = N / (pow(v, 3.0) + pow(T, 3.0)) * (1.0 + erf((Vb - v) / DeltaV));
			else if constexpr (disType == 3)
				pw = N / (pow(v, 3.0) + pow(T, 3.0)) * exp(-pow(Lambda - Lambda0, 2.0) / DeltaLambda2);
			else if constexpr (disType == 4)
				pw = N / (pow(v, 3.0) + pow(T, 3.0)) * exp(-pow(Lambda - Lambda0, 2.0) / DeltaLambda2) * (1.0 + erf((Vb - v) / DeltaV));

			if constexpr (markerType == 0)
				pw /= N;
			else if constexpr (markerType == 1)
				pw *= 1.0;
			else if constexpr (markerType == 2)
				pw *= J;

			li = (E - minE) / dE;
			lj = (Pphi - minPphi) / dPphi;
			lk = (Lambda - minLambda) / dLambda;

			i = floor(li);
			j = floor(lj);
			k = floor(lk);

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
				coes[index] = (cx[index] + sx[index] * dx) * (cy[index] + sy[index] * dy) * (cz[index] + sz[index] * dz);

			phaseSpaceF0[i][j][k] += pw * coes[0];
			phaseSpaceF0[i + 1][j][k] += pw * coes[1];
			phaseSpaceF0[i][j + 1][k] += pw * coes[2];
			phaseSpaceF0[i + 1][j + 1][k] += pw * coes[3];
			phaseSpaceF0[i][j][k + 1] += pw * coes[4];
			phaseSpaceF0[i + 1][j][k + 1] += pw * coes[5];
			phaseSpaceF0[i][j + 1][k + 1] += pw * coes[6];
			phaseSpaceF0[i + 1][j + 1][k + 1] += pw * coes[7];

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;
			std::cout << std::endl;
		}

		MPICHECK(MPI_Allreduce(MPI_IN_PLACE, phaseSpaceF0[0][0], gridE * gridPphi * gridLambda, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

		for (int i = 0; i < gridE; i++) {
			for (int j = 0; j < gridPphi; j++) {
				for (int k = 0; k < gridLambda; k++) {

					bool isVertex = (i == 0 || i == gridE - 1) &&
						(j == 0 || j == gridPphi - 1) &&
						(k == 0 || k == gridLambda - 1);

					bool isEdge = ((i == 0 || i == gridE - 1) && (j == 0 || j == gridPphi - 1)) ||
						((i == 0 || i == gridE - 1) && (k == 0 || k == gridLambda - 1)) ||
						((j == 0 || j == gridPphi - 1) && (k == 0 || k == gridLambda - 1));

					bool isFace = (i == 0 || i == gridE - 1 || j == 0 || j == gridPphi - 1 || k == 0 || k == gridLambda - 1);

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
				fileName = "IonPhaseSpaceF0.bin";
			else if constexpr (picType == 1)
				fileName = "AlphaPhaseSpaceF0.bin";
			else if constexpr (picType == 2)
				fileName = "BeamPhaseSpaceF0.bin";

			output.open(fileName.c_str(), std::ios::out | std::ios::binary);

			int tempGridE = gridE, tempGridPphi = gridPphi, tempGridLambda = gridLambda;
			output.write((char*)(&tempGridE), sizeof(int));
			output.write((char*)(&tempGridPphi), sizeof(int));
			output.write((char*)(&tempGridLambda), sizeof(int));

			output.write((char*)(&minE), sizeof(double));
			output.write((char*)(&maxE), sizeof(double));

			output.write((char*)(&minPphi), sizeof(double));
			output.write((char*)(&maxPphi), sizeof(double));

			output.write((char*)(&minLambda), sizeof(double));
			output.write((char*)(&maxLambda), sizeof(double));

			output.write((char*)(phaseSpaceF0[0][0]), sizeof(double) * gridE * gridPphi * gridLambda);

			output.close();

			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;

		}

		HostAllocator.releaseHostArrays(phaseSpaceF0);

	}
	template<int picType, int gridE, int gridPphi, int gridLambda, int ppcPhase>
	void computePhaseSpaceJacobian() {

		if (hostId == 0) {

			if constexpr (picType == 0)
				std::cout << BOLDYELLOW << "Start: Compute phase space jacobian of thermal ions." << RESET << std::endl;
			else if constexpr (picType == 1)
				std::cout << BOLDYELLOW << "Start: Compute phase space jacobian of alpha particles." << RESET << std::endl;
			else if constexpr (picType == 2)
				std::cout << BOLDYELLOW << "Start: Compute phase space jacobian of beam particles." << RESET << std::endl;

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

		}
		else if constexpr (picType == 1) {

			Mass = AlphaMass;
			Char = AlphaChar;
			Vmin = AlphaVmin;
			Vmax = AlphaVmax;

		}
		else if constexpr (picType == 2) {

			Mass = BeamMass;
			Char = BeamChar;
			Vmin = BeamVmin;
			Vmax = BeamVmax;

		}

		minE = 0.5 * Mass * pow(Vmin, 2.0);
		maxE = 0.5 * Mass * pow(Vmax, 2.0);

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

				tempPphi = cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) * SFAcovyz[i][j + gridGhost] / (q[i][j] * J[i][j] * B[i][j]) - Char * psip[i][j];
				if (tempPphi > maxPphi)
					maxPphi = tempPphi;

				tempPphi = -cm * Mass * Vmax * 2 * psitmax * drho * (RHO0 + i * gridDx * drho) * SFAcovyz[i][j + gridGhost] / (q[i][j] * J[i][j] * B[i][j]) - Char * psip[i][j];
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

		double Jmax, Jvmax;
		Jmax = 0.0;
		Jvmax = pow(Vmax, 2.0);

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
			int i = floor(li);
			int j = floor(lj);
			double dx = li - i;
			double dy = lj - j;

			double coes[4] = {};
			double cx[4] = { 1.0, 1.0, 0.0, 0.0 };
			double sx[4] = { -1.0, -1.0, 1.0, 1.0 };
			double cy[4] = { 1.0, 0.0, 1.0, 0.0 };
			double sy[4] = { -1.0, 1.0, -1.0, 1.0 };

			double result = 0.0;

			coes[0] = (cx[0] + sx[0] * dx) * (cy[0] + sy[0] * dy);
			coes[1] = (cx[1] + sx[1] * dx) * (cy[1] + sy[1] * dy);
			coes[2] = (cx[2] + sx[2] * dx) * (cy[2] + sy[2] * dy);
			coes[3] = (cx[3] + sx[3] * dx) * (cy[3] + sy[3] * dy);

			result = field[i][j] * coes[0] + field[i][j + 1] * coes[1] + field[i + 1][j] * coes[2] + field[i + 1][j + 1] * coes[3];

			return result;

		};

		Rand01 xrand;
		Rand01 yrand;
		Rand01 Jrand;

		Rand01 vrand;
		Rand01 vprand;
		Rand01 Jvrand;

		double cx[8] = { 1.0, 0.0, 1.0, 0.0,1.0, 0.0, 1.0, 0.0 };
		double sx[8] = { -1.0, 1.0, -1.0, 1.0,-1.0, 1.0, -1.0, 1.0 };
		double cy[8] = { 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0 };
		double sy[8] = { -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0 };
		double cz[8] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
		double sz[8] = { -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 };

		for (size_t picId = 0; picId < picPhase; picId++) {

			if (picId % (picPhase / 10) == 0)
				if (hostId == 0)
					std::cout << BOLDGREEN << 100 * picId / picPhase << "%" << RESET << std::endl;

			int i, j, k;
			double li, lj, lk;
			double dx, dy, dz;
			double q, psip, J, B, SFAcovyz;
			double v, Jv, E, Pphi, Lambda;
			double x, y, z, vp, mu;
			double coes[8];

			do {
				x = x0 + (x1 - x0) * xrand();
				y = y0 + (y1 - y0) * yrand();
				J = interp2d(tempJ, x, y);
			} while (Jrand() >= J / Jmax);

			do {
				v = Vmax * vrand();
				Jv = pow(v, 2.0);
			} while (v <= Vmin || Jvrand() >= Jv / Jvmax);

			q = interp2d(tempq, x, y);
			psip = interp2d(temppsip, x, y);
			SFAcovyz = interp2d(tempSFAcovyz, x, y);
			B = interp2d(tempB, x, y);

			if constexpr (picType == 2)
				Lambda = vprand();
			else
				Lambda = 2.0 * vprand() - 1.0;
			vp = v * Lambda;
			mu = 0.5 * Mass * pow(v, 2.0) * (1.0 - pow(Lambda, 2.0)) / B;

			E = 0.5 * Mass * pow(v, 2.0);
			Pphi = cm * Mass * vp * 2 * psitmax * drho * (RHO0 + x * drho) * SFAcovyz / (q * J * B) - Char * psip;
			Lambda = mu / E;

			li = (E - minE) / dE;
			lj = (Pphi - minPphi) / dPphi;
			lk = (Lambda - minLambda) / dLambda;

			i = floor(li);
			j = floor(lj);
			k = floor(lk);

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
				coes[index] = (cx[index] + sx[index] * dx) * (cy[index] + sy[index] * dy) * (cz[index] + sz[index] * dz);

			phaseSpaceJacobian[i][j][k] += coes[0];
			phaseSpaceJacobian[i + 1][j][k] += coes[1];
			phaseSpaceJacobian[i][j + 1][k] += coes[2];
			phaseSpaceJacobian[i + 1][j + 1][k] += coes[3];
			phaseSpaceJacobian[i][j][k + 1] += coes[4];
			phaseSpaceJacobian[i + 1][j][k + 1] += coes[5];
			phaseSpaceJacobian[i][j + 1][k + 1] += coes[6];
			phaseSpaceJacobian[i + 1][j + 1][k + 1] += coes[7];

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;
			std::cout << std::endl;
		}

		MPICHECK(MPI_Allreduce(MPI_IN_PLACE, phaseSpaceJacobian[0][0], gridE * gridPphi * gridLambda, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

		for (int i = 0; i < gridE; i++) {
			for (int j = 0; j < gridPphi; j++) {
				for (int k = 0; k < gridLambda; k++) {

					bool isVertex = (i == 0 || i == gridE - 1) &&
						(j == 0 || j == gridPphi - 1) &&
						(k == 0 || k == gridLambda - 1);

					bool isEdge = ((i == 0 || i == gridE - 1) && (j == 0 || j == gridPphi - 1)) ||
						((i == 0 || i == gridE - 1) && (k == 0 || k == gridLambda - 1)) ||
						((j == 0 || j == gridPphi - 1) && (k == 0 || k == gridLambda - 1));

					bool isFace = (i == 0 || i == gridE - 1 || j == 0 || j == gridPphi - 1 || k == 0 || k == gridLambda - 1);

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
				fileName = "IonPhaseSpaceJacobian.bin";
			else if constexpr (picType == 1)
				fileName = "AlphaPhaseSpaceJacobian.bin";
			else if constexpr (picType == 2)
				fileName = "BeamPhaseSpaceJacobian.bin";

			output.open(fileName.c_str(), std::ios::out | std::ios::binary);

			int tempGridE = gridE, tempGridPphi = gridPphi, tempGridLambda = gridLambda;
			output.write((char*)(&tempGridE), sizeof(int));
			output.write((char*)(&tempGridPphi), sizeof(int));
			output.write((char*)(&tempGridLambda), sizeof(int));

			output.write((char*)(&minE), sizeof(double));
			output.write((char*)(&maxE), sizeof(double));

			output.write((char*)(&minPphi), sizeof(double));
			output.write((char*)(&maxPphi), sizeof(double));

			output.write((char*)(&minLambda), sizeof(double));
			output.write((char*)(&maxLambda), sizeof(double));

			output.write((char*)(phaseSpaceJacobian[0][0]), sizeof(double) * gridE * gridPphi * gridLambda);

			output.close();

			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;

		}

		HostAllocator.releaseHostArrays(phaseSpaceJacobian);

	}
	template<int picType>
	void loadPhaseSpaceMapping(std::string file) {

		if (hostId == 0) {

			if constexpr (picType == 0)
				std::cout << BOLDYELLOW << "Start: Load phase space mapping of thermal ions." << RESET << std::endl;
			else if constexpr (picType == 1)
				std::cout << BOLDYELLOW << "Start: Load phase space mapping of alpha particles." << RESET << std::endl;
			else if constexpr (picType == 2)
				std::cout << BOLDYELLOW << "Start: Load phase space mapping of beam particles." << RESET << std::endl;

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

		input.read((char*)(ids.data()), picPhase * sizeof(int));
		input.read((char*)(rhos.data()), picPhase * sizeof(double));
		input.read((char*)(vparas.data()), picPhase * sizeof(double));
		input.read((char*)(mus.data()), picPhase * sizeof(double));
		input.close();

		if constexpr (picType == 0) {

			Allocator HostDeviceAllocator;
			HostDeviceAllocator.allocateHostArrays(devNums, (size_t)picPhase * 13, h_IonPhaseSpaceMapping);
			HostDeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picPhase * 13, d_IonPhaseSpaceMapping);

			for (int i = 0; i < devNums; i++) {
				for (int picId = 0; picId < picPhase; picId++) {

					//orbit = 0.5 : pad
					//orbit = 1.5 : loss
					//orbit = 2.5 : para
					//orbit = 3.5 : anti
					//orbit = 4.5 : trapped
					//orbit = 5.5 : unknown

					//orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

					if (ids[picId] == 20251106)
						h_IonPhaseSpaceMapping[i][picId * 13 + 0] = 0.5;
					else
						h_IonPhaseSpaceMapping[i][picId * 13 + 0] = 5.5;

					h_IonPhaseSpaceMapping[i][picId * 13 + 1] = rhos[picId];
					h_IonPhaseSpaceMapping[i][picId * 13 + 3] = vparas[picId];
					h_IonPhaseSpaceMapping[i][picId * 13 + 4] = mus[picId];
					h_IonPhaseSpaceMapping[i][picId * 13 + 9] = 0.5;

				}
			}

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				HostDeviceAllocator.hostToDevice((size_t)picPhase * 13, 0, (size_t)i * picPhase * 13, d_IonPhaseSpaceMapping[i], h_IonPhaseSpaceMapping[0]);

			}

		}
		else if constexpr (picType == 1) {

			Allocator HostDeviceAllocator;
			HostDeviceAllocator.allocateHostArrays(devNums, (size_t)picPhase * 13, h_AlphaPhaseSpaceMapping);
			HostDeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picPhase * 13, d_AlphaPhaseSpaceMapping);

			for (int i = 0; i < devNums; i++) {
				for (int picId = 0; picId < picPhase; picId++) {

					//orbit = 0.5 : pad
					//orbit = 1.5 : loss
					//orbit = 2.5 : para
					//orbit = 3.5 : anti
					//orbit = 4.5 : trapped
					//orbit = 5.5 : unknown

					//orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

					if (ids[picId] == 20251106)
						h_AlphaPhaseSpaceMapping[i][picId * 13 + 0] = 0.5;
					else
						h_AlphaPhaseSpaceMapping[i][picId * 13 + 0] = 5.5;

					h_AlphaPhaseSpaceMapping[i][picId * 13 + 1] = rhos[picId];
					h_AlphaPhaseSpaceMapping[i][picId * 13 + 3] = vparas[picId];
					h_AlphaPhaseSpaceMapping[i][picId * 13 + 4] = mus[picId];
					h_AlphaPhaseSpaceMapping[i][picId * 13 + 9] = 0.5;

					//h_AlphaPhaseSpaceMapping[i][picId * 13 + 0] = 5.5;
					//h_AlphaPhaseSpaceMapping[i][picId * 13 + 1] = 0.292253437648246;
					//h_AlphaPhaseSpaceMapping[i][picId * 13 + 2] = 0.0;
					//h_AlphaPhaseSpaceMapping[i][picId * 13 + 3] = -0.568243955750533;
					//h_AlphaPhaseSpaceMapping[i][picId * 13 + 4] = 2.355984934257279;
					//h_AlphaPhaseSpaceMapping[i][picId * 13 + 9] = 0.5;

				}
			}

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				HostDeviceAllocator.hostToDevice((size_t)picPhase * 13, 0, (size_t)i * picPhase * 13, d_AlphaPhaseSpaceMapping[i], h_AlphaPhaseSpaceMapping[0]);

			}

		}
		else if constexpr (picType == 2) {

			Allocator HostDeviceAllocator;
			HostDeviceAllocator.allocateHostArrays(devNums, (size_t)picPhase * 13, h_BeamPhaseSpaceMapping);
			HostDeviceAllocator.allocateDeviceArrays(localId, devNums, (size_t)picPhase * 13, d_BeamPhaseSpaceMapping);

			for (int i = 0; i < devNums; i++) {
				for (int picId = 0; picId < picPhase; picId++) {

					//orbit = 0.5 : pad
					//orbit = 1.5 : loss
					//orbit = 2.5 : para
					//orbit = 3.5 : anti
					//orbit = 4.5 : trapped
					//orbit = 5.5 : unknown

					//orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

					if (ids[picId] == 20251106)
						h_BeamPhaseSpaceMapping[i][picId * 13 + 0] = 0.5;
					else
						h_BeamPhaseSpaceMapping[i][picId * 13 + 0] = 5.5;

					h_BeamPhaseSpaceMapping[i][picId * 13 + 1] = rhos[picId];
					h_BeamPhaseSpaceMapping[i][picId * 13 + 3] = vparas[picId];
					h_BeamPhaseSpaceMapping[i][picId * 13 + 4] = mus[picId];
					h_BeamPhaseSpaceMapping[i][picId * 13 + 9] = 0.5;

				}
			}

			for (int i = 0; i < devNums; i++) {

				CUDACHECK(cudaSetDevice(localId * devNums + i));
				HostDeviceAllocator.hostToDevice((size_t)picPhase * 13, 0, (size_t)i * picPhase * 13, d_BeamPhaseSpaceMapping[i], h_BeamPhaseSpaceMapping[0]);

			}

		}

		if (hostId == 0) {
			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;
		}

	}
	template<int picType>
	void computePhaseSpaceOrbit(std::string file) {

		if (hostId == 0) {

			if constexpr (picType == 0)
				std::cout << BOLDYELLOW << "Start: Compute phase space orbit frequency of thermal ions." << RESET << std::endl;
			else if constexpr (picType == 1)
				std::cout << BOLDYELLOW << "Start: Compute phase space orbit frequency of alpha particles." << RESET << std::endl;
			else if constexpr (picType == 2)
				std::cout << BOLDYELLOW << "Start: Compute phase space orbit frequency of beam particles." << RESET << std::endl;

		}

		std::ifstream input;
		input.open(file, std::ios::in | std::ios::binary);

		input.seekg(0, std::ios::end);
		std::streamsize size = input.tellg();
		input.seekg(0, std::ios::beg);
		int picPhase = size / 28;

		std::vector<int> ids(picPhase);
		input.read((char*)(ids.data()), picPhase * sizeof(int));
		input.close();

		std::vector<double> PhaseSpaceOrbit(picPhase * 9);

		if constexpr (picType == 0) {

			for (int picId = 0; picId < picPhase; picId++) {

				//orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

				PhaseSpaceOrbit[picId * 9 + 0] = ids[picId];
				PhaseSpaceOrbit[picId * 9 + 1] = h_IonPhaseSpaceMapping[0][picId * 13 + 0];
				PhaseSpaceOrbit[picId * 9 + 2] = h_IonPhaseSpaceMapping[0][picId * 13 + 5];
				PhaseSpaceOrbit[picId * 9 + 3] = h_IonPhaseSpaceMapping[0][picId * 13 + 6];
				PhaseSpaceOrbit[picId * 9 + 4] = h_IonPhaseSpaceMapping[0][picId * 13 + 7];
				PhaseSpaceOrbit[picId * 9 + 5] = h_IonPhaseSpaceMapping[0][picId * 13 + 8];
				PhaseSpaceOrbit[picId * 9 + 6] = h_IonPhaseSpaceMapping[0][picId * 13 + 10];
				PhaseSpaceOrbit[picId * 9 + 7] = h_IonPhaseSpaceMapping[0][picId * 13 + 11];
				PhaseSpaceOrbit[picId * 9 + 8] = h_IonPhaseSpaceMapping[0][picId * 13 + 12];

			}

			Allocator HostDeviceAllocator;
			HostDeviceAllocator.releaseHostArrays(h_IonPhaseSpaceMapping);
			HostDeviceAllocator.releaseDeviceArrays(localId, devNums, d_IonPhaseSpaceMapping);

		}
		else if constexpr (picType == 1) {

			for (int picId = 0; picId < picPhase; picId++) {

				//orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

				PhaseSpaceOrbit[picId * 9 + 0] = ids[picId];
				PhaseSpaceOrbit[picId * 9 + 1] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 0];
				PhaseSpaceOrbit[picId * 9 + 2] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 5];
				PhaseSpaceOrbit[picId * 9 + 3] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 6];
				PhaseSpaceOrbit[picId * 9 + 4] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 7];
				PhaseSpaceOrbit[picId * 9 + 5] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 8];
				PhaseSpaceOrbit[picId * 9 + 6] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 10];
				PhaseSpaceOrbit[picId * 9 + 7] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 11];
				PhaseSpaceOrbit[picId * 9 + 8] = h_AlphaPhaseSpaceMapping[0][picId * 13 + 12];

			}

			Allocator HostDeviceAllocator;
			HostDeviceAllocator.releaseHostArrays(h_AlphaPhaseSpaceMapping);
			HostDeviceAllocator.releaseDeviceArrays(localId, devNums, d_AlphaPhaseSpaceMapping);

		}
		else if constexpr (picType == 2) {

			for (int picId = 0; picId < picPhase; picId++) {

				//orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

				PhaseSpaceOrbit[picId * 9 + 0] = ids[picId];
				PhaseSpaceOrbit[picId * 9 + 1] = h_BeamPhaseSpaceMapping[0][picId * 13 + 0];
				PhaseSpaceOrbit[picId * 9 + 2] = h_BeamPhaseSpaceMapping[0][picId * 13 + 5];
				PhaseSpaceOrbit[picId * 9 + 3] = h_BeamPhaseSpaceMapping[0][picId * 13 + 6];
				PhaseSpaceOrbit[picId * 9 + 4] = h_BeamPhaseSpaceMapping[0][picId * 13 + 7];
				PhaseSpaceOrbit[picId * 9 + 5] = h_BeamPhaseSpaceMapping[0][picId * 13 + 8];
				PhaseSpaceOrbit[picId * 9 + 6] = h_BeamPhaseSpaceMapping[0][picId * 13 + 10];
				PhaseSpaceOrbit[picId * 9 + 7] = h_BeamPhaseSpaceMapping[0][picId * 13 + 11];
				PhaseSpaceOrbit[picId * 9 + 8] = h_BeamPhaseSpaceMapping[0][picId * 13 + 12];

			}

			Allocator HostDeviceAllocator;
			HostDeviceAllocator.releaseHostArrays(h_BeamPhaseSpaceMapping);
			HostDeviceAllocator.releaseDeviceArrays(localId, devNums, d_BeamPhaseSpaceMapping);

		}

		if (hostId == 0) {

			std::ofstream output;
			std::string fileName;

			if constexpr (picType == 0)
				fileName = "IonPhaseSpaceOrbit.bin";
			else if constexpr (picType == 1)
				fileName = "AlphaPhaseSpaceOrbit.bin";
			else if constexpr (picType == 2)
				fileName = "BeamPhaseSpaceOrbit.bin";

			output.open(fileName.c_str(), std::ios::out | std::ios::binary);
			output.write((char*)(PhaseSpaceOrbit.data()), sizeof(double) * picPhase * 9);
			output.close();

			std::cout << BOLDGREEN << "Done." << RESET << std::endl;
			std::cout << std::endl;

		}

	}

};