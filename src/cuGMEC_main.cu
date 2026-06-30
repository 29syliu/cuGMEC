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

#include "cuGMEC_model.h"
#include "cuGMEC_const.h"
#include "cuGMEC_param.h"
#include "cuGMEC_setup.h"

#include "kernels/utils.cuh"
#include "kernels/mhd.cuh"
#include "kernels/pic.cuh"

int main(int argc, char* argv[]) {

    /*--------------------------------Phase 1: MPI Initialization---------------------------------*/

    int myRank, nRanks, localRank = 0;

    MPICHECK(MPI_Init(&argc, &argv));
    MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
    MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));

    if (myRank == 0)
        std::cout << BOLDYELLOW << "Start: MPI Initialization." << RESET << std::endl;

    if (nRanks != hostNums) {
        std::cout << "Error: nRanks != hostNums." << std::endl;
        return 0;
    }

    if (myRank == 0) {
        std::cout << BOLDGREEN << "Done." << RESET << std::endl;
        std::cout << std::endl;
    }

    /*--------------------------------Phase 2: NCCL Initialization--------------------------------*/

    if (myRank == 0)
        std::cout << BOLDYELLOW << "Start: NCCL Initialization." << RESET << std::endl;

    uint64_t hostHashs[nRanks];
    char hostname[1024];
    getHostName(hostname, 1024);
    hostHashs[myRank] = getHostHash(hostname);
    MPICHECK(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hostHashs, sizeof(uint64_t), MPI_BYTE, MPI_COMM_WORLD));
    for (int p = 0; p < nRanks; p++) {
        if (p == myRank)
            break;
        if (hostHashs[p] == hostHashs[myRank])
            localRank++;
    }

    ncclUniqueId id;
    ncclComm_t comms[devNums];
    if (myRank == 0)
        ncclGetUniqueId(&id);
    MPICHECK(MPI_Bcast((void*)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));

    NCCLCHECK(ncclGroupStart());
    for (int i = 0; i < devNums; i++) {
        CUDACHECK(cudaSetDevice(localRank * devNums + i));
        NCCLCHECK(ncclCommInitRank(comms + i, nRanks * devNums, id, myRank * devNums + i));
    }
    NCCLCHECK(ncclGroupEnd());

    if (myRank == 0) {
        std::cout << BOLDGREEN << "Done." << RESET << std::endl;
        std::cout << std::endl;
    }

    std::vector<int> ncclLeftNei(devNums);
    std::vector<int> ncclRightNei(devNums);
    int ncclLeftSend;
    int ncclLeftRecv;
    int ncclRightSend;
    int ncclRightRecv;

    for (int i = 0; i < devNums; i++) {
        ncclLeftNei[i] = (myRank * devNums + i - 1 + nRanks * devNums) % (nRanks * devNums);
        ncclRightNei[i] = (myRank * devNums + i + 1) % (nRanks * devNums);
        ncclLeftSend = gridGhost * gridNxz;
        ncclLeftRecv = 0;
        ncclRightSend = devNy * gridNxz;
        ncclRightRecv = (devNy + gridGhost) * gridNxz;
    }

    ncclDataType_t ncclType;
    if constexpr (std::is_same_v<mhdReal, double>)
        ncclType = ncclDouble;
    else
        ncclType = ncclFloat;

    /*----------------------------Phase 3: HybridModel Initialization-----------------------------*/

    HybridModel<mhdReal, picReal> cuGMEC(HybridModelConfig{
        .scale = {.devNums = devNums,
                  .hostNums = nRanks,
                  .hostId = myRank,
                  .localId = localRank,
                  .gridNx = gridNx,
                  .gridNy = gridNy,
                  .gridNz = gridNz,
                  .NFP = NFP,
                  .gridGhost = gridGhost,
                  .ppcNums = ppcNums,
                  .tubes = tubes},
        .norm = {.B0 = B0, .L0 = L0, .VA0 = VA0, .RHO0 = RHO0, .RHO1 = RHO1, .PSITMAX = PSITMAX, .dt = dt},
        .hyperDiff = {.A = {toBool<ifNablaPerp2A>, perp2A},
                      .Phi = {toBool<ifNablaPerp2Phi>, perp2Phi},
                      .dNe = {toBool<ifNablaPerp2dNe>, perp2dNe},
                      .dTe = {toBool<ifNablaPerp2dTe>, perp2dTe},
                      .dP = {toBool<ifNablaPerp2dP>, perp2dP}},
        .ion = {toBool<ifIon>, IonMass, IonChar, IonBeta, IonVmin, IonVmax, IonVb, IonDeltaV, IonLambda0,
                IonDeltaLambda2},
        .alpha = {toBool<ifAlpha>, AlphaMass, AlphaChar, AlphaBeta, AlphaVmin, AlphaVmax, AlphaVb, AlphaDeltaV,
                  AlphaLambda0, AlphaDeltaLambda2},
        .beam = {toBool<ifBeam>, BeamMass, BeamChar, BeamBeta, BeamVmin, BeamVmax, BeamVb, BeamDeltaV, BeamLambda0,
                 BeamDeltaLambda2},
        .time = {.total = totalSteps, .diag = diagSteps, .output = outputSteps},
        .filter = {.leftN = leftN, .rightN = rightN},
        .nFFT = {.time = nFFTTimeSize, .batch = nFFTBatchSize, .freq = nFFTFreqSize},
        .mFFT = {.time = mFFTTimeSize},
        .diag = {.amplitude = toBool<ifDiagAmplitude>,
                 .frequency = toBool<ifDiagFrequency>,
                 .Eparallel = toBool<ifDiagEparallel>,
                 .density = toBool<ifDiagDensity>,
                 .diffusivity = toBool<ifDiagDiffusivity>,
                 .ZFDrive = toBool<ifDiagZFDrive>,
                 .checkNAN = toBool<ifCheckNAN>},
        .output = {.Phi = toBool<ifOutputPhi>,
                   .A = toBool<ifOutputA>,
                   .dNe = toBool<ifOutputdNe>,
                   .dTe = toBool<ifOutputdTe>,
                   .dP = toBool<ifOutputdP>,
                   .dPi = toBool<ifOutputdPi>,
                   .dPa = toBool<ifOutputdPa>,
                   .dPb = toBool<ifOutputdPb>}});

    cuGMEC.setup();

    /*---------------------------Phase 4: Kernel Launch Configurations----------------------------*/

    dim3 MRK4GridSize(MRK4GridDimx, MRK4GridDimy, MRK4GridDimz);
    dim3 MRK4BlockSize(MRK4BlockDimx, MRK4BlockDimy, MRK4BlockDimz);

    dim3 GhostGridSize(GhostGridDimx, GhostGridDimy, GhostGridDimz);
    dim3 GhostBlockSize(GhostBlockDimx, GhostBlockDimy, GhostBlockDimz);

    dim3 M2PGridSize(M2PGridDimx, M2PGridDimy, M2PGridDimz);
    dim3 M2PBlockSize(M2PBlockDimx, M2PBlockDimy, M2PBlockDimz);

    dim3 LocalNMGridSize(LocalNMGridDimx);
    dim3 GhostNMGridSize(GhostNMGridDimx);
    dim3 RefinedNMGridSize(RefinedNMGridDimx);
    dim3 NMBlockSize(NMBlockDimx);

    dim3 PICGridSize(PICGridDimx);
    dim3 PICBlockSize(PICBlockDimx);

    dim3 MergeGridSize(MergeGridDimx);
    dim3 MergeBlockSize(MergeBlockDimx);

    /*------------------------------Phase 5: cuGMEC Resource Binding------------------------------*/

    int currentStep = 0;
    int diagIndex = 0;
    int outputIndex = 0;
    auto& start = cuGMEC.startEvents;
    auto& end = cuGMEC.endEvents;
    auto& time = cuGMEC.elapsedTime;
    auto& cudssHandles = cuGMEC.cudssHandles;

#define CUDSSBIND(P)                                                                                                   \
    auto& P##Configs = cuGMEC.P##Configs;                                                                              \
    auto& P##Datas = cuGMEC.P##Datas;                                                                                  \
    auto& P##As = cuGMEC.P##As;                                                                                        \
    auto& P##Xs = cuGMEC.P##Xs;                                                                                        \
    auto& P##Bs = cuGMEC.P##Bs;

    CUDSSBIND(laplacian)
    CUDSSBIND(resistive)
    CUDSSBIND(Phi)
    CUDSSBIND(dNe)
    CUDSSBIND(dTe)
    CUDSSBIND(dP)

#undef CUDSSBIND

    // clang-format off
#define REF(x) auto& x = cuGMEC.d_##x

    REF(qtheta);
    REF(w_beg);   REF(w_midl);   REF(w_midr);   REF(w_end);
    REF(A_beg);   REF(A_midl);   REF(A_midr);   REF(A_end);
    REF(dNe_beg); REF(dNe_midl); REF(dNe_midr); REF(dNe_end);
    REF(dTe_beg); REF(dTe_midl); REF(dTe_midr); REF(dTe_end);
    REF(dP_beg);  REF(dP_midl);  REF(dP_midr);  REF(dP_end);

    REF(Phi_midl);  REF(Phi_midr);
    REF(dJpB_midl); REF(dJpB_midr);
    REF(dPe_midl);  REF(dPe_midr);
    REF(Apt_midl);  REF(Apt_midr);
    REF(Ne0); REF(Te0); REF(Ne0_px); REF(Te0_px);

    REF(w2Phi); REF(A2dJpB); REF(Phi2w);
    REF(wdPAdJpB2w); REF(APhidNe2A); REF(dPePhiAdJpB2dNe); REF(PhidTedNe2dTe); REF(Phi2dP);
    REF(wPhi_w); REF(AdJpB_w); REF(PhiA_A); REF(NeA_A);
    REF(dNePhi_dNe); REF(PhiTe_dTe); REF(PhiTeA_dTe); REF(dPPhi_dP);

    REF(pic1d); REF(pic2d); REF(pic3d);
    REF(globalA);  REF(globalPhi); REF(globalApt);
    REF(globalPa); REF(globalPi);  REF(globalPb);
    REF(globalNa); REF(globalNi);  REF(globalNb);
    REF(dPa_midl); REF(dPa_midr); REF(dPi_midl); REF(dPi_midr); REF(dPb_midl); REF(dPb_midr);
    REF(dNa_midl); REF(dNa_midr); REF(dNi_midl); REF(dNi_midr); REF(dNb_midl); REF(dNb_midr);

    REF(Alpha_keys_in); REF(Alpha_keys_out); REF(Alpha_sort_ids_in); REF(Alpha_sort_ids_out);
    REF(Ion_keys_in);   REF(Ion_keys_out);   REF(Ion_sort_ids_in);   REF(Ion_sort_ids_out);
    REF(Beam_keys_in);  REF(Beam_keys_out);  REF(Beam_sort_ids_in);  REF(Beam_sort_ids_out);
    REF(Alpha_values_in); REF(Alpha_values_out);
    REF(Ion_values_in);   REF(Ion_values_out);
    REF(Beam_values_in);  REF(Beam_values_out);

    REF(amplitude); REF(frequency); REF(modeReal); REF(modeImag); REF(Epara); REF(EparaES);
    REF(MaxwellDrive); REF(ReynoldsDrive); REF(dwdtTotal);
    REF(Maxwell); REF(Reynolds);
    REF(IonDensity);     REF(AlphaDensity);     REF(BeamDensity);
    REF(IonDiffusivity); REF(AlphaDiffusivity); REF(BeamDiffusivity);
    REF(IonPhaseDeltaF); REF(AlphaPhaseDeltaF); REF(BeamPhaseDeltaF);
    REF(IonPitchDeltaF); REF(AlphaPitchDeltaF); REF(BeamPitchDeltaF);
    REF(IonPhasePower);  REF(AlphaPhasePower);  REF(BeamPhasePower);
    REF(IonPitchPower);  REF(AlphaPitchPower);  REF(BeamPitchPower);
    REF(NANFlag);
    REF(totalPhi); REF(totalA); REF(totaldNe); REF(totaldTe); REF(totaldP);
    REF(totaldPi); REF(totaldPa); REF(totaldPb);

    REF(nmLocal); REF(nmGlobal); REF(nmRefined);

    REF(nPlanR2Cs); REF(nPlanC2Rs); REF(nFreqd); REF(nFreqf);
    REF(mPlanC2Cs);

    REF(Ion_storage);   REF(Ion_storage_bytes);
    REF(Alpha_storage); REF(Alpha_storage_bytes);
    REF(Beam_storage);  REF(Beam_storage_bytes);

#undef REF
    // clang-format on

    /*--------------------------------Phase 6: Lambda Definitions---------------------------------*/

    auto forEachDev = [&](auto&& body) {
        for (int i = 0; i < devNums; i++) {
            cudaSetDevice(localRank * devNums + i);
            body(i);
        }
    };

    auto haloExchange = [&]<typename... Guards>(auto... fields) {
        if constexpr (allTrue<Guards...>) {
            if constexpr (devNy >= gridGhost) {
                for (int i = 0; i < devNums; i++) {
                    cudaSetDevice(localRank * devNums + i);
                    (
                        [&](mhdReal** field) {
                            ncclSend(field[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i],
                                     0);
                            ncclRecv(field[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i],
                                     0);
                            ncclSend(field[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i],
                                     0);
                            ncclRecv(field[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i],
                                     0);
                        }(fields),
                        ...);
                }
            } else {
                for (int g = 0; g < gridGhost; g++) {
                    for (int i = 0; i < devNums; i++) {
                        cudaSetDevice(localRank * devNums + i);

                        const int ncclSelf = myRank * devNums + i;
                        const int ncclTotal = nRanks * devNums;
                        const int ncclLayerShift = g / devNy + 1;
                        const int ncclLayerId = g % devNy;

                        const int ncclLeftLayerNei = (ncclSelf - ncclLayerShift + ncclTotal) % ncclTotal;
                        const int ncclRightLayerNei = (ncclSelf + ncclLayerShift) % ncclTotal;
                        const int ncclLeftLayerSend = (gridGhost + ncclLayerId) * gridNxz;
                        const int ncclRightLayerSend = (gridGhost + devNy - 1 - ncclLayerId) * gridNxz;
                        const int ncclLeftLayerRecv = (gridGhost - 1 - g) * gridNxz;
                        const int ncclRightLayerRecv = (gridGhost + devNy + g) * gridNxz;

                        (
                            [&](mhdReal** field) {
                                ncclSend(field[i] + ncclLeftLayerSend, gridNxz, ncclType, ncclLeftLayerNei, comms[i],
                                         0);
                                ncclRecv(field[i] + ncclRightLayerRecv, gridNxz, ncclType, ncclRightLayerNei, comms[i],
                                         0);
                                ncclSend(field[i] + ncclRightLayerSend, gridNxz, ncclType, ncclRightLayerNei, comms[i],
                                         0);
                                ncclRecv(field[i] + ncclLeftLayerRecv, gridNxz, ncclType, ncclLeftLayerNei, comms[i],
                                         0);
                            }(fields),
                            ...);
                    }
                }
            }
        }
    };

    auto loadMHDPIC = [&]() {
        forEachDev([&](int i) {
            MHD2w<ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(Phi_midl[i], w_midl[i], Phi2w[i]);
            if constexpr (std::is_same_v<ifQNeutrality, trueType>) {
                MHDQNeutrality2dNe<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(w_midl[i], dNi_midl[i], dNa_midl[i],
                                                                          dNb_midl[i], dNe_midl[i]);
            }
            MHD2dJpBdPePhi<ifNonlinearMHD, ifLocal, ifFLRMHD, ifReducedMHD><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                A_midl[i], dJpB_midl[i], A2dJpB[i], w_midl[i], Phi_midl[i], w2Phi[i], dNe_midl[i], dTe_midl[i],
                dP_midl[i], dPe_midl[i], Ne0[i], Te0[i]);
            cudaMemcpyAsync(w_beg[i], w_midl[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                            cudaMemcpyDeviceToDevice, 0);
        });

        ncclGroupStart();
        haloExchange(w_midl, A_midl, dNe_midl, dTe_midl, dP_midl, Phi_midl, dJpB_midl, dPe_midl);
        haloExchange.template operator()<ifIon>(dPi_midl);
        haloExchange.template operator()<ifAlpha>(dPa_midl);
        haloExchange.template operator()<ifBeam>(dPb_midl);
        haloExchange.template operator()<ifIon, ifQNeutrality>(dNi_midl);
        haloExchange.template operator()<ifAlpha, ifQNeutrality>(dNa_midl);
        haloExchange.template operator()<ifBeam, ifQNeutrality>(dNb_midl);
        ncclGroupEnd();
    };

    auto runMHDRK4Stage = [&]<int stage>(mhdReal** w_in = nullptr, mhdReal** A_in = nullptr, mhdReal** dNe_in = nullptr,
                                         mhdReal** dTe_in = nullptr, mhdReal** dP_in = nullptr,
                                         mhdReal** w_out = nullptr, mhdReal** A_out = nullptr,
                                         mhdReal** dNe_out = nullptr, mhdReal** dTe_out = nullptr,
                                         mhdReal** dP_out = nullptr) {
        if constexpr (stage != 5) {
            forEachDev([&](int i) {
                MHDLinearRK4<stage, ifNonlinearMHD, ifLocal, ifEparallel, ifReducedMHD>
                    <<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                        qtheta[i], w_beg[i], w_in[i], w_out[i], w_end[i], A_beg[i], A_in[i], A_out[i], A_end[i],
                        dNe_beg[i], dNe_in[i], dNe_out[i], dNe_end[i], dTe_beg[i], dTe_in[i], dTe_out[i], dTe_end[i],
                        dP_beg[i], dP_in[i], dP_out[i], dP_end[i], Phi_midl[i], dJpB_midl[i], dPe_midl[i], dPi_midl[i],
                        dPa_midl[i], dPb_midl[i], wdPAdJpB2w[i], APhidNe2A[i], dPePhiAdJpB2dNe[i], PhidTedNe2dTe[i],
                        Phi2dP[i]);
                if constexpr (std::is_same_v<ifNonlinearMHD, trueType>) {
                    MHDNonlinearRK4<stage, ifMaxwellStress, ifReynoldsStress, ifDiagZFDrive, ifLocal, ifEparallel,
                                    ifReducedMHD><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                        qtheta[i], w_in[i], w_out[i], w_end[i], A_in[i], A_out[i], A_end[i], dNe_in[i], dNe_out[i],
                        dNe_end[i], dTe_in[i], dTe_out[i], dTe_end[i], dP_in[i], dP_out[i], dP_end[i], Phi_midl[i],
                        dJpB_midl[i], dPe_midl[i], Ne0[i], Te0[i], Ne0_px[i], Te0_px[i], APhidNe2A[i], wPhi_w[i],
                        AdJpB_w[i], PhiA_A[i], NeA_A[i], dNePhi_dNe[i], PhiTe_dTe[i], PhiTeA_dTe[i], dPPhi_dP[i],
                        Maxwell[i], Reynolds[i]);
                }
                for (int j = 0; j < devNy; j++) {
                    cudssMatrixSetValues(laplacianBs[i][j], w_out[i] + (j + gridGhost) * gridNxz);
                    cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, laplacianConfigs[i][j], laplacianDatas[i][j],
                                 laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]);
                }
                if constexpr (std::is_same_v<ifQNeutrality, trueType>) {
                    MHDQNeutrality2dNe<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(w_out[i], dNi_midl[i], dNa_midl[i],
                                                                              dNb_midl[i], dNe_out[i]);
                }
                if constexpr (stage != 4) {
                    MHD2dJpBdPePhi<ifNonlinearMHD, ifLocal, ifFLRMHD, ifReducedMHD>
                        <<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(A_out[i], dJpB_midl[i], A2dJpB[i], w_out[i],
                                                                Phi_midl[i], w2Phi[i], dNe_out[i], dTe_out[i],
                                                                dP_out[i], dPe_midl[i], Ne0[i], Te0[i]);
                }
                if constexpr (stage == 4 && std::is_same_v<ifDiagZFDrive, trueType>) {
                    cudaMemcpyAsync(w_in[i], w_beg[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                    cudaMemcpyDeviceToDevice, 0);
                }
            });
            if constexpr (stage != 4) {
                ncclGroupStart();
                haloExchange(w_out, A_out, dNe_out, dTe_out, dP_out, Phi_midl, dJpB_midl, dPe_midl);
                ncclGroupEnd();
            }
        } else {
            forEachDev([&](int i) {
                MHD2w<ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(Phi_midl[i], w_midl[i], Phi2w[i]);
                if constexpr (std::is_same_v<ifQNeutrality, trueType>) {
                    MHDQNeutrality2dNe<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(w_midl[i], dNi_midl[i], dNa_midl[i],
                                                                              dNb_midl[i], dNe_midl[i]);
                }
                MHD2dJpBdPePhi<ifNonlinearMHD, ifLocal, ifFLRMHD, ifReducedMHD><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                    A_midl[i], dJpB_midl[i], A2dJpB[i], w_midl[i], Phi_midl[i], w2Phi[i], dNe_midl[i], dTe_midl[i],
                    dP_midl[i], dPe_midl[i], Ne0[i], Te0[i]);
                cudaMemcpyAsync(w_beg[i], w_midl[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                cudaMemcpyDeviceToDevice, 0);
                cudaMemcpyAsync(A_beg[i], A_midl[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                cudaMemcpyDeviceToDevice, 0);
                cudaMemcpyAsync(dNe_beg[i], dNe_midl[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                cudaMemcpyDeviceToDevice, 0);
                cudaMemcpyAsync(dTe_beg[i], dTe_midl[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                cudaMemcpyDeviceToDevice, 0);
                cudaMemcpyAsync(dP_beg[i], dP_midl[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                cudaMemcpyDeviceToDevice, 0);
            });
            ncclGroupStart();
            haloExchange(w_midl, A_midl, dNe_midl, dTe_midl, dP_midl, Phi_midl, dJpB_midl, dPe_midl);
            ncclGroupEnd();
        }
    };

    auto nablaPerp2 = [&]<typename... Guards>(auto& Configs, auto& Datas, auto& As, auto& Xs, auto& Bs, mhdReal** midl,
                                              mhdReal** midr) {
        if constexpr (allTrue<Guards...>) {
            forEachDev([&](int i) {
                for (int j = 0; j < devNy; j++) {
                    CUDSSCHECK(cudssMatrixSetValues(Xs[i][j], midr[i] + (j + gridGhost) * gridNxz));
                    CUDSSCHECK(cudssMatrixSetValues(Bs[i][j], midl[i] + (j + gridGhost) * gridNxz));
                    CUDSSCHECK(cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, Configs[i][j], Datas[i][j], As[i][j],
                                            Xs[i][j], Bs[i][j]));
                }
                cudaMemcpyAsync(midl[i], midr[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                cudaMemcpyDeviceToDevice, 0);
            });
        }
    };

    auto nablaPara4 = [&]<int Tag, typename... Guards>(mhdReal** midl, mhdReal** midr) {
        if constexpr (allTrue<Guards...>) {
            ncclGroupStart();
            haloExchange(midl);
            ncclGroupEnd();

            forEachDev([&](int i) {
                MHDNablaPara2<ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(qtheta[i], midl[i], midr[i]);
            });

            ncclGroupStart();
            haloExchange(midr);
            ncclGroupEnd();

            forEachDev([&](int i) {
                MHDNablaPara4<Tag, ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(qtheta[i], midl[i], midr[i]);
            });
        }
    };

    auto filterModeN = [&]<typename... Guards>(mhdReal** midl, int left, int right) {
        if constexpr (allTrue<Guards...>) {
            forEachDev([&](int i) {
                if constexpr (std::is_same_v<mhdReal, double>) {
                    cufftExecD2Z(nPlanR2Cs[i], (double*)midl[i] + gridGhost * gridNxz, nFreqd[i]);
                    MHDFilterModeN<<<devNy, gridNx, 0, 0>>>(nFreqd[i], left, right);
                    cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)midl[i] + gridGhost * gridNxz);
                } else {
                    cufftExecR2C(nPlanR2Cs[i], (float*)midl[i] + gridGhost * gridNxz, nFreqf[i]);
                    MHDFilterModeN<<<devNy, gridNx, 0, 0>>>(nFreqf[i], left, right);
                    cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)midl[i] + gridGhost * gridNxz);
                }
                MHDFilterResizeN<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(midl[i]);
            });
        }
    };

    auto removeModeN = [&]<const auto & removeN, typename... Guards>(mhdReal** midl, mhdReal** midr) {
        if constexpr (allTrue<Guards...> && removeN.size() > 0) {
            for (int tor : removeN) {
                forEachDev([&](int i) {
                    if constexpr (std::is_same_v<mhdReal, double>) {
                        cufftExecD2Z(nPlanR2Cs[i], (double*)midl[i] + gridGhost * gridNxz, nFreqd[i]);
                        MHDFilterModeN<<<devNy, gridNx, 0, 0>>>(nFreqd[i], tor);
                        cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)midr[i] + gridGhost * gridNxz);
                    } else {
                        cufftExecR2C(nPlanR2Cs[i], (float*)midl[i] + gridGhost * gridNxz, nFreqf[i]);
                        MHDFilterModeN<<<devNy, gridNx, 0, 0>>>(nFreqf[i], tor);
                        cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)midr[i] + gridGhost * gridNxz);
                    }
                    MHDFilterResizeN<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(midr[i]);
                    MHDSubtractMode<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(midr[i], midl[i]);
                });
            }
        }
    };

    auto selectModeNM = [&]<const auto & selectNM, typename... Guards>(mhdReal** midl) {
        if constexpr (allTrue<Guards...> && selectNM.size() > 0) {
            for (const auto& [toroidal, poloidalLeft, poloidalRight] : selectNM) {
                forEachDev([&](int i) {
                    if constexpr (std::is_same_v<mhdReal, double>) {
                        cufftExecD2Z(nPlanR2Cs[i], (double*)midl[i] + gridGhost * gridNxz, nFreqd[i]);
                        MHDSelectNMSubtractMode<<<LocalNMGridSize, NMBlockSize, 0, 0>>>(nFreqd[i], midl[i], nmLocal[i],
                                                                                        toroidal);
                    } else {
                        cufftExecR2C(nPlanR2Cs[i], (float*)midl[i] + gridGhost * gridNxz, nFreqf[i]);
                        MHDSelectNMSubtractMode<<<LocalNMGridSize, NMBlockSize, 0, 0>>>(nFreqf[i], midl[i], nmLocal[i],
                                                                                        toroidal);
                    }

                    MHDSelectNMShifted2A<<<LocalNMGridSize, NMBlockSize, 0, 0>>>(qtheta[i], nmLocal[i], nmLocal[i],
                                                                                 toroidal);
                });

                ncclGroupStart();
                for (int i = 0; i < devNums; i++) {
                    cudaSetDevice(localRank * devNums + i);
                    ncclAllGather(nmLocal[i], nmGlobal[i] + gridGhost * gridNx, 2 * devNy * gridNx, ncclType, comms[i],
                                  0);
                }
                ncclGroupEnd();

                forEachDev([&](int i) {
                    int localIndex = myRank * hostNy + i * devNy;

                    MHDSelectNMAlignedGhost<<<GhostNMGridSize, NMBlockSize, 0, 0>>>(qtheta[i], nmGlobal[i], toroidal);
                    MHDSelectNMRefineY<mhdComplex, mhdReal>
                        <<<RefinedNMGridSize, NMBlockSize, 0, 0>>>(nmGlobal[i], nmRefined[i]);
                    MHDSelectNMAligned2S<<<RefinedNMGridSize, NMBlockSize, 0, 0>>>(qtheta[i], nmRefined[i],
                                                                                   nmRefined[i], localIndex, toroidal);

                    if constexpr (std::is_same_v<mhdReal, double>) {
                        cufftExecZ2Z(mPlanC2Cs[i], reinterpret_cast<cufftDoubleComplex*>(nmRefined[i]),
                                     reinterpret_cast<cufftDoubleComplex*>(nmRefined[i]), CUFFT_FORWARD);
                    } else {
                        cufftExecC2C(mPlanC2Cs[i], reinterpret_cast<cufftComplex*>(nmRefined[i]),
                                     reinterpret_cast<cufftComplex*>(nmRefined[i]), CUFFT_FORWARD);
                    }

                    MHDSelectNMFilterM<<<RefinedNMGridSize, NMBlockSize, 0, 0>>>(nmRefined[i], poloidalLeft,
                                                                                 poloidalRight);

                    if constexpr (std::is_same_v<mhdReal, double>) {
                        cufftExecZ2Z(mPlanC2Cs[i], reinterpret_cast<cufftDoubleComplex*>(nmRefined[i]),
                                     reinterpret_cast<cufftDoubleComplex*>(nmRefined[i]), CUFFT_INVERSE);
                    } else {
                        cufftExecC2C(mPlanC2Cs[i], reinterpret_cast<cufftComplex*>(nmRefined[i]),
                                     reinterpret_cast<cufftComplex*>(nmRefined[i]), CUFFT_INVERSE);
                    }

                    MHDSelectNMAddMode<<<LocalNMGridSize, NMBlockSize, 0, 0>>>(nmRefined[i], midl[i], localIndex,
                                                                               toroidal);
                });
            }
        }
    };

    auto mhdToPICMode = [&]<typename... Guards>(int mode) {
        if constexpr (allTrue<Guards...>) {
            auto filterModeN = [&](mhdReal** src, mhdReal** dst, int leftMode, int rightMode) {
                forEachDev([&](int i) {
                    if constexpr (std::is_same_v<mhdReal, double>) {
                        cufftExecD2Z(nPlanR2Cs[i], (double*)src[i] + gridGhost * gridNxz, nFreqd[i]);
                        MHDFilterModeN<<<devNy, gridNx, 0, 0>>>(nFreqd[i], leftMode, rightMode);
                        cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dst[i] + gridGhost * gridNxz);
                    } else {
                        cufftExecR2C(nPlanR2Cs[i], (float*)src[i] + gridGhost * gridNxz, nFreqf[i]);
                        MHDFilterModeN<<<devNy, gridNx, 0, 0>>>(nFreqf[i], leftMode, rightMode);
                        cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dst[i] + gridGhost * gridNxz);
                    }
                    MHDFilterResizeN<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(dst[i]);
                });
            };

            forEachDev([&](int i) {
                MHD2Apt<ifNonlinearMHD, ifLocal, ifEparallel><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                    qtheta[i], A_midl[i], dNe_midl[i], dTe_midl[i], Phi_midl[i], Ne0[i], Te0[i], Ne0_px[i],
                    APhidNe2A[i], PhiA_A[i], NeA_A[i], Apt_midl[i]);
            });

            filterModeN(Phi_midl, w_midr, mode, mode);
            filterModeN(A_midl, dJpB_midr, mode, mode);
            filterModeN(Apt_midl, dPe_midr, mode, mode);

            forEachDev([&](int i) {
                MHDShifted2A<0, ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                    qtheta[i], dJpB_midr[i], A_midr[i], w_midr[i], Phi_midr[i], dPe_midr[i], Apt_midr[i]);
            });

            ncclGroupStart();
            for (int i = 0; i < devNums; i++) {
                cudaSetDevice(localRank * devNums + i);
                ncclAllGather(A_midr[i] + gridGhost * gridNxz, globalA[i] + gridGhost * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
                ncclAllGather(Phi_midr[i] + gridGhost * gridNxz, globalPhi[i] + gridGhost * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
                ncclAllGather(Apt_midr[i] + gridGhost * gridNxz, globalApt[i] + gridGhost * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            }
            ncclGroupEnd();

            forEachDev([&](int i) {
                MHDAlignedGhost<ifLocal>
                    <<<GhostGridSize, GhostBlockSize, 0, 0>>>(qtheta[i], globalA[i], globalPhi[i], globalApt[i]);
                MHD2PIC<mhdReal, picReal>
                    <<<M2PGridSize, M2PBlockSize, 0, 0>>>(pic3d[i], globalA[i], globalPhi[i], globalApt[i]);
            });
        }
    };

    auto runDiagnostics = [&](int diagIdx) {
        if (diagIdx % diagSteps != 0)
            return;

        if constexpr (std::is_same_v<ifDiagAmplitude, trueType>) {
            forEachDev([&](int i) {
                if constexpr (std::is_same_v<mhdReal, double>) {
                    cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
                    MHDDiagAmplitude<<<1, gridNx, 0, 0>>>(
                        nFreqd[i], amplitude[i] + diagIdx / diagSteps * (rightN - leftN + 1) * gridNx,
                        modeReal[i] + diagIdx / diagSteps * (rightN - leftN + 1) * gridNx,
                        modeImag[i] + diagIdx / diagSteps * (rightN - leftN + 1) * gridNx);
                } else {
                    cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
                    MHDDiagAmplitude<<<1, gridNx, 0, 0>>>(
                        nFreqf[i], amplitude[i] + diagIdx / diagSteps * (rightN - leftN + 1) * gridNx,
                        modeReal[i] + diagIdx / diagSteps * (rightN - leftN + 1) * gridNx,
                        modeImag[i] + diagIdx / diagSteps * (rightN - leftN + 1) * gridNx);
                }
            });
        }

        if constexpr (std::is_same_v<ifDiagFrequency, trueType>) {
            forEachDev([&](int i) {
                MHDDiagFrequency<<<1, gridNx, 0, 0>>>(Phi_midl[i], frequency[i] + diagIdx / diagSteps * gridNx);
            });
        }

        if constexpr (std::is_same_v<ifDiagEparallel, trueType>) {
            forEachDev([&](int i) {
                MHDDiagEparallel<ifNonlinearMHD, ifEparallel><<<8, gridNx / 8, 0, 0>>>(
                    qtheta[i], A_midl[i], dNe_midl[i], dTe_midl[i], Phi_midl[i], Ne0[i], Te0[i], Ne0_px[i],
                    APhidNe2A[i], PhiA_A[i], NeA_A[i], Epara[i] + diagIdx / diagSteps * gridNx,
                    EparaES[i] + diagIdx / diagSteps * gridNx);
            });
        }

        if constexpr (std::is_same_v<ifDiagZFDrive, trueType>) {
            filterModeN.template operator()<ifDiagZFDrive>(Maxwell, 0, 0);
            filterModeN.template operator()<ifDiagZFDrive>(Reynolds, 0, 0);
            filterModeN.template operator()<ifDiagZFDrive>(w_midr, 0, 0);
            filterModeN.template operator()<ifDiagZFDrive>(w_end, 0, 0);

            forEachDev([&](int i) {
                MHDDiagZFDrive<<<1, gridNx, 0, 0>>>(
                    w_midr[i], w_end[i], Maxwell[i], Reynolds[i], MaxwellDrive[i] + diagIdx / diagSteps * gridNx,
                    ReynoldsDrive[i] + diagIdx / diagSteps * gridNx, dwdtTotal[i] + diagIdx / diagSteps * gridNx);
            });
        }

        if constexpr (std::is_same_v<ifDiagDensity, trueType>) {
            forEachDev([&](int i) {
                if constexpr (std::is_same_v<ifIon, trueType>)
                    PICDiagDensity<Ion, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic2d[i], Ion_keys_in[i], Ion_values_in[i], IonDensity[i] + diagIdx / diagSteps * gridNx);
                if constexpr (std::is_same_v<ifAlpha, trueType>)
                    PICDiagDensity<Alpha, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic2d[i], Alpha_keys_in[i], Alpha_values_in[i], AlphaDensity[i] + diagIdx / diagSteps * gridNx);
                if constexpr (std::is_same_v<ifBeam, trueType>)
                    PICDiagDensity<Beam, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic2d[i], Beam_keys_in[i], Beam_values_in[i], BeamDensity[i] + diagIdx / diagSteps * gridNx);
            });
        }

        if constexpr (std::is_same_v<ifDiagDiffusivity, trueType>) {
            for (int mode = leftN; mode <= rightN; mode++) {
                int modeIdx = mode - leftN;
                mhdToPICMode.template operator()<ifPIC>(mode);

                forEachDev([&](int i) {
                    const size_t offset = ((size_t)diagIdx / diagSteps * (rightN - leftN + 1) + modeIdx) * gridNx;

                    if constexpr (std::is_same_v<ifIon, trueType>)
                        PICDiagDiffusivity<Ion, mhdReal, picReal, ifFLRPIC, gyroNums>
                            <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Ion_keys_in[i],
                                                                  Ion_values_in[i], IonDiffusivity[i] + offset);
                    if constexpr (std::is_same_v<ifAlpha, trueType>)
                        PICDiagDiffusivity<Alpha, mhdReal, picReal, ifFLRPIC, gyroNums>
                            <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Alpha_keys_in[i],
                                                                  Alpha_values_in[i], AlphaDiffusivity[i] + offset);
                    if constexpr (std::is_same_v<ifBeam, trueType>)
                        PICDiagDiffusivity<Beam, mhdReal, picReal, ifFLRPIC, gyroNums>
                            <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Beam_keys_in[i],
                                                                  Beam_values_in[i], BeamDiffusivity[i] + offset);
                });
            }
        }
    };

    auto runOutput = [&](int outputIdx) {
        if (outputIdx % outputSteps != 0)
            return;

        ncclGroupStart();
        for (int i = 0; i < devNums; i++) {
            cudaSetDevice(localRank * devNums + i);
            if constexpr (std::is_same_v<ifOutputPhi, trueType>)
                ncclAllGather(Phi_midl[i] + gridGhost * gridNxz,
                              totalPhi[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            if constexpr (std::is_same_v<ifOutputA, trueType>)
                ncclAllGather(A_midl[i] + gridGhost * gridNxz,
                              totalA[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz, ncclType,
                              comms[i], 0);
            if constexpr (std::is_same_v<ifOutputdNe, trueType>)
                ncclAllGather(dNe_midl[i] + gridGhost * gridNxz,
                              totaldNe[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            if constexpr (std::is_same_v<ifOutputdTe, trueType>)
                ncclAllGather(dTe_midl[i] + gridGhost * gridNxz,
                              totaldTe[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            if constexpr (std::is_same_v<ifOutputdP, trueType>)
                ncclAllGather(dP_midl[i] + gridGhost * gridNxz,
                              totaldP[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            if constexpr (std::is_same_v<ifOutputdPi, trueType>)
                ncclAllGather(dPi_midl[i] + gridGhost * gridNxz,
                              totaldPi[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            if constexpr (std::is_same_v<ifOutputdPa, trueType>)
                ncclAllGather(dPa_midl[i] + gridGhost * gridNxz,
                              totaldPa[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            if constexpr (std::is_same_v<ifOutputdPb, trueType>)
                ncclAllGather(dPb_midl[i] + gridGhost * gridNxz,
                              totaldPb[i] + (size_t)outputIdx / outputSteps * gridNy * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
        }
        ncclGroupEnd();

        if constexpr (std::is_same_v<ifOutputPhaseSpaceDeltaF, trueType>) {
            forEachDev([&](int i) {
                if constexpr (std::is_same_v<ifIon, trueType>)
                    PICDiagPhaseDeltaF<Ion, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic1d[i], pic2d[i], Ion_values_in[i],
                        IonPhaseDeltaF[i] + (size_t)outputIdx / outputSteps * gridE * gridPphi * gridLambda);
                if constexpr (std::is_same_v<ifAlpha, trueType>)
                    PICDiagPhaseDeltaF<Alpha, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic1d[i], pic2d[i], Alpha_values_in[i],
                        AlphaPhaseDeltaF[i] + (size_t)outputIdx / outputSteps * gridE * gridPphi * gridLambda);
                if constexpr (std::is_same_v<ifBeam, trueType>)
                    PICDiagPhaseDeltaF<Beam, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic1d[i], pic2d[i], Beam_values_in[i],
                        BeamPhaseDeltaF[i] + (size_t)outputIdx / outputSteps * gridE * gridPphi * gridLambda);
            });
        }

        if constexpr (std::is_same_v<ifOutputPitchSpaceDeltaF, trueType>) {
            forEachDev([&](int i) {
                if constexpr (std::is_same_v<ifIon, trueType>)
                    PICDiagPitchDeltaF<Ion, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic1d[i], pic2d[i], Ion_values_in[i],
                        IonPitchDeltaF[i] + (size_t)outputIdx / outputSteps * gridVpara * gridVperp);
                if constexpr (std::is_same_v<ifAlpha, trueType>)
                    PICDiagPitchDeltaF<Alpha, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic1d[i], pic2d[i], Alpha_values_in[i],
                        AlphaPitchDeltaF[i] + (size_t)outputIdx / outputSteps * gridVpara * gridVperp);
                if constexpr (std::is_same_v<ifBeam, trueType>)
                    PICDiagPitchDeltaF<Beam, mhdReal, picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic1d[i], pic2d[i], Beam_values_in[i],
                        BeamPitchDeltaF[i] + (size_t)outputIdx / outputSteps * gridVpara * gridVperp);
            });
        }

        if constexpr (std::is_same_v<ifOutputPhaseSpacePower, trueType> ||
                      std::is_same_v<ifOutputPitchSpacePower, trueType>) {
            for (int mode = leftN; mode <= rightN; mode++) {
                int modeIdx = mode - leftN;
                mhdToPICMode.template operator()<ifPIC>(mode);

                forEachDev([&](int i) {
                    const size_t phaseOffset = ((size_t)outputIdx / outputSteps * (rightN - leftN + 1) + modeIdx) *
                                               gridE * gridPphi * gridLambda;
                    const size_t pitchOffset =
                        ((size_t)outputIdx / outputSteps * (rightN - leftN + 1) + modeIdx) * gridVpara * gridVperp;

                    if constexpr (std::is_same_v<ifOutputPhaseSpacePower, trueType>) {
                        if constexpr (std::is_same_v<ifIon, trueType>)
                            PICDiagPhasePower<Ion, mhdReal, picReal, ifFLRPIC, gyroNums>
                                <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Ion_keys_in[i],
                                                                      Ion_values_in[i], IonPhasePower[i] + phaseOffset);
                        if constexpr (std::is_same_v<ifAlpha, trueType>)
                            PICDiagPhasePower<Alpha, mhdReal, picReal, ifFLRPIC, gyroNums>
                                <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Alpha_keys_in[i],
                                                                      Alpha_values_in[i],
                                                                      AlphaPhasePower[i] + phaseOffset);
                        if constexpr (std::is_same_v<ifBeam, trueType>)
                            PICDiagPhasePower<Beam, mhdReal, picReal, ifFLRPIC, gyroNums>
                                <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Beam_keys_in[i],
                                                                      Beam_values_in[i],
                                                                      BeamPhasePower[i] + phaseOffset);
                    }

                    if constexpr (std::is_same_v<ifOutputPitchSpacePower, trueType>) {
                        if constexpr (std::is_same_v<ifIon, trueType>)
                            PICDiagPitchPower<Ion, mhdReal, picReal, ifFLRPIC, gyroNums>
                                <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Ion_keys_in[i],
                                                                      Ion_values_in[i], IonPitchPower[i] + pitchOffset);
                        if constexpr (std::is_same_v<ifAlpha, trueType>)
                            PICDiagPitchPower<Alpha, mhdReal, picReal, ifFLRPIC, gyroNums>
                                <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Alpha_keys_in[i],
                                                                      Alpha_values_in[i],
                                                                      AlphaPitchPower[i] + pitchOffset);
                        if constexpr (std::is_same_v<ifBeam, trueType>)
                            PICDiagPitchPower<Beam, mhdReal, picReal, ifFLRPIC, gyroNums>
                                <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], Beam_keys_in[i],
                                                                      Beam_values_in[i],
                                                                      BeamPitchPower[i] + pitchOffset);
                    }
                });
            }
        }
    };

    auto mhdToPIC = [&]<typename... Guards>() {
        if constexpr (allTrue<Guards...>) {
            forEachDev([&](int i) {
                MHD2Apt<ifNonlinearMHD, ifLocal, ifEparallel><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                    qtheta[i], A_midl[i], dNe_midl[i], dTe_midl[i], Phi_midl[i], Ne0[i], Te0[i], Ne0_px[i],
                    APhidNe2A[i], PhiA_A[i], NeA_A[i], Apt_midl[i]);
                MHDShifted2A<0, ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                    qtheta[i], A_midl[i], A_midr[i], Phi_midl[i], Phi_midr[i], Apt_midl[i], Apt_midr[i]);
            });

            ncclGroupStart();
            for (int i = 0; i < devNums; i++) {
                cudaSetDevice(localRank * devNums + i);
                ncclAllGather(A_midr[i] + gridGhost * gridNxz, globalA[i] + gridGhost * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
                ncclAllGather(Phi_midr[i] + gridGhost * gridNxz, globalPhi[i] + gridGhost * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
                ncclAllGather(Apt_midr[i] + gridGhost * gridNxz, globalApt[i] + gridGhost * gridNxz, devNy * gridNxz,
                              ncclType, comms[i], 0);
            }
            ncclGroupEnd();

            forEachDev([&](int i) {
                MHDAlignedGhost<ifLocal>
                    <<<GhostGridSize, GhostBlockSize, 0, 0>>>(qtheta[i], globalA[i], globalPhi[i], globalApt[i]);
                MHD2PIC<mhdReal, picReal>
                    <<<M2PGridSize, M2PBlockSize, 0, 0>>>(pic3d[i], globalA[i], globalPhi[i], globalApt[i]);
            });
        }
    };

    auto runPICRK4 = [&]<picType particle, disType distribution, typename... Guards>(
                         int** keys_in, picReal** values_in, mhdReal** globalP, mhdReal** globalN) {
        if constexpr (allTrue<Guards...>) {
            forEachDev([&](int i) {
                if constexpr (std::is_same_v<ifFLRPIC, trueType>)
                    GyroAlignedRK4<ratioDt, gyroNums, particle, distribution, ifNonlinearPIC, ifQNeutrality, mhdReal,
                                   picReal><<<PICGridSize, PICBlockSize, 0, 0>>>(
                        pic1d[i], pic2d[i], pic3d[i], keys_in[i], values_in[i], globalP[i], globalN[i]);
                else
                    DriftAlignedRK4<ratioDt, particle, distribution, ifNonlinearPIC, ifQNeutrality, mhdReal, picReal>
                        <<<PICGridSize, PICBlockSize, 0, 0>>>(pic1d[i], pic2d[i], pic3d[i], keys_in[i], values_in[i],
                                                              globalP[i], globalN[i]);
            });
        }
    };

    auto mergePICBuffers = [&]<typename... Guards>(mhdReal** globalP) {
        if constexpr (allTrue<Guards...>) {
            if constexpr (depositBufferNums > 1) {
                forEachDev([&](int i) {
                    PICMergeBuffers<depositBufferNums, mhdReal><<<MergeGridSize, MergeBlockSize, 0, 0>>>(globalP[i]);
                });
            }
        }
    };

    auto allReducePressure = [&]<typename... Guards>(mhdReal** globalP) {
        if constexpr (allTrue<Guards...>) {
            for (int i = 0; i < devNums; i++) {
                cudaSetDevice(localRank * devNums + i);
                ncclAllReduce(globalP[i], globalP[i], (gridNy + 2 * gridGhost) * gridNxz, ncclType, ncclSum, comms[i],
                              0);
            }
        }
    };

    auto updateAlignedGhost = [&]<typename... Guards>(mhdReal** globalP, mhdReal** dP_midl) {
        if constexpr (allTrue<Guards...>) {
            forEachDev([&](int i) {
                PICAlignedGhost<ifLocal><<<GhostGridSize, GhostBlockSize, 0, 0>>>(qtheta[i], globalP[i]);
                cudaMemcpyAsync(dP_midl[i], globalP[i] + (myRank * hostNy + i * devNy) * gridNxz,
                                sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
                cudaMemsetAsync(globalP[i], 0, sizeof(mhdReal) * (gridNy + 2 * gridGhost) * gridNxz, 0);
            });
        }
    };

    auto picToMHD = [&]<typename... Guards>() {
        if constexpr (allTrue<Guards...>) {
            forEachDev([&](int i) {
                MHDShifted2A<1, ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                    qtheta[i], dPi_midl[i], dPi_midr[i], dPa_midl[i], dPa_midr[i], dPb_midl[i], dPb_midr[i]);
                if constexpr (std::is_same_v<ifQNeutrality, trueType>) {
                    MHDShifted2A<1, ifLocal><<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(
                        qtheta[i], dNi_midl[i], dNi_midr[i], dNa_midl[i], dNa_midr[i], dNb_midl[i], dNb_midr[i]);
                }
            });

            auto right2Left = [&]<typename... SubGuards>(mhdReal** midl, mhdReal** midr) {
                if constexpr (allTrue<SubGuards...>) {
                    forEachDev([&](int i) {
                        cudaMemcpyAsync(midl[i], midr[i], sizeof(mhdReal) * (devNy + 2 * gridGhost) * gridNxz,
                                        cudaMemcpyDeviceToDevice, 0);
                    });
                }
            };

            right2Left.template operator()<ifIon>(dPi_midl, dPi_midr);
            right2Left.template operator()<ifAlpha>(dPa_midl, dPa_midr);
            right2Left.template operator()<ifBeam>(dPb_midl, dPb_midr);
            right2Left.template operator()<ifIon, ifQNeutrality>(dNi_midl, dNi_midr);
            right2Left.template operator()<ifAlpha, ifQNeutrality>(dNa_midl, dNa_midr);
            right2Left.template operator()<ifBeam, ifQNeutrality>(dNb_midl, dNb_midr);
        }
    };

    auto sortParticles = [&]<typename... Guards>(void** storage, std::vector<size_t>& storage_bytes, int** keys_in,
                                                 int** keys_out, int** sort_ids_in, int** sort_ids_out,
                                                 picReal** values_in, picReal** values_out) {
        if constexpr (allTrue<Guards...>) {
            forEachDev([&](int i) {
                cub::DeviceRadixSort::SortPairs(storage[i], storage_bytes[i], keys_in[i], keys_out[i], sort_ids_in[i],
                                                sort_ids_out[i], picDev);
                PICReorderValues<<<SortGridDimx, SortBlockDimx, 0, 0>>>(sort_ids_out[i], values_in[i], values_out[i]);

                int* temp_keys = keys_in[i];
                keys_in[i] = keys_out[i];
                keys_out[i] = temp_keys;

                picReal* temp_values = values_in[i];
                values_in[i] = values_out[i];
                values_out[i] = temp_values;
            });
        }
    };

    auto checkNAN = [&]<typename... Guards>() -> bool {
        if constexpr (allTrue<Guards...>) {
            int local = 0;
            forEachDev([&](int i) {
                cudaMemsetAsync(NANFlag[i], 0, sizeof(int), 0);
                MHDCheckNAN<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(NANFlag[i], Phi_midl[i], A_midl[i], dNe_midl[i],
                                                                   dTe_midl[i], dP_midl[i]);
            });
            forEachDev([&](int i) {
                int devLocal;
                cudaMemcpy(&devLocal, NANFlag[i], sizeof(int), cudaMemcpyDeviceToHost);
                local |= devLocal;
            });
            int global;
            MPICHECK(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));
            return global != 0;
        } else {
            return false;
        }
    };

    /*-------------------------------Phase 7: Initial State Loading-------------------------------*/

    loadMHDPIC();
    runDiagnostics(diagIndex);
    runOutput(outputIndex);

    /*------------------------------Phase 8: Main Time-Stepping Loop------------------------------*/

    if (myRank == 0) {
        std::cout << BOLDGREEN << "Gyrokinetic-MHD hybrid simulation is running." << RESET << std::endl;
        std::cout << std::endl;
    }

    forEachDev([&](int i) { cudaEventRecord(start[i]); });

    for (int outerLoop = 0; outerLoop < outerLoopMax; outerLoop++) {

        if constexpr (outerLoopMax / 10 != 0) {

            if (outerLoop % (outerLoopMax / 10) == 0) {
                if (myRank == 0) {
                    std::cout << BOLDGREEN << 100 * outerLoop / outerLoopMax << "%" << RESET << std::endl;
                }
            }

        } else {

            if (myRank == 0) {
                std::cout << BOLDGREEN << 100 * outerLoop / outerLoopMax << "%" << RESET << std::endl;
            }
        }

        for (int innerLoop = 0; innerLoop < innerLoopMax; innerLoop++) {

            for (int MHDLoop = 0; MHDLoop < MHDLoopMax; MHDLoop++) {

                currentStep++;

                /*----------------------------------MHD RK4-----------------------------------*/

                runMHDRK4Stage.operator()<1>(w_midl, A_midl, dNe_midl, dTe_midl, dP_midl, w_midr, A_midr, dNe_midr,
                                             dTe_midr, dP_midr);
                runMHDRK4Stage.operator()<2>(w_midr, A_midr, dNe_midr, dTe_midr, dP_midr, w_midl, A_midl, dNe_midl,
                                             dTe_midl, dP_midl);
                runMHDRK4Stage.operator()<3>(w_midl, A_midl, dNe_midl, dTe_midl, dP_midl, w_midr, A_midr, dNe_midr,
                                             dTe_midr, dP_midr);
                runMHDRK4Stage.operator()<4>(w_midr, A_midr, dNe_midr, dTe_midr, dP_midr, w_midl, A_midl, dNe_midl,
                                             dTe_midl, dP_midl);

                nablaPerp2.template operator()<ifNablaPerp2Phi>(PhiConfigs, PhiDatas, PhiAs, PhiXs, PhiBs, Phi_midl,
                                                                Phi_midr);
                nablaPerp2.template operator()<ifNablaPerp2A>(resistiveConfigs, resistiveDatas, resistiveAs,
                                                              resistiveXs, resistiveBs, A_midl, A_midr);
                nablaPerp2.template operator()<ifNablaPerp2dNe>(dNeConfigs, dNeDatas, dNeAs, dNeXs, dNeBs, dNe_midl,
                                                                dNe_midr);
                nablaPerp2.template operator()<ifNablaPerp2dTe>(dTeConfigs, dTeDatas, dTeAs, dTeXs, dTeBs, dTe_midl,
                                                                dTe_midr);
                nablaPerp2.template operator()<ifNablaPerp2dP, ifReducedMHD>(dPConfigs, dPDatas, dPAs, dPXs, dPBs,
                                                                             dP_midl, dP_midr);

                nablaPara4.template operator()<0, ifNablaPara4Phi>(Phi_midl, Phi_midr);
                nablaPara4.template operator()<1, ifNablaPara4A>(A_midl, A_midr);
                nablaPara4.template operator()<2, ifNablaPara4dNe>(dNe_midl, dNe_midr);
                nablaPara4.template operator()<3, ifNablaPara4dTe>(dTe_midl, dTe_midr);
                nablaPara4.template operator()<4, ifNablaPara4dP, ifReducedMHD>(dP_midl, dP_midr);

                filterModeN.template operator()<ifFilterN_Phi>(Phi_midl, leftN, rightN);
                filterModeN.template operator()<ifFilterN_A>(A_midl, leftN, rightN);
                filterModeN.template operator()<ifFilterN_dNe>(dNe_midl, leftN, rightN);
                filterModeN.template operator()<ifFilterN_dTe>(dTe_midl, leftN, rightN);
                filterModeN.template operator()<ifFilterN_dP, ifReducedMHD>(dP_midl, leftN, rightN);

                removeModeN.template operator()<removeN_Phi>(Phi_midl, Phi_midr);
                removeModeN.template operator()<removeN_A>(A_midl, A_midr);
                removeModeN.template operator()<removeN_dNe>(dNe_midl, dNe_midr);
                removeModeN.template operator()<removeN_dTe>(dTe_midl, dTe_midr);
                removeModeN.template operator()<removeN_dP, ifReducedMHD>(dP_midl, dP_midr);

                selectModeNM.template operator()<selectNM_Phi>(Phi_midl);
                selectModeNM.template operator()<selectNM_A>(A_midl);
                selectModeNM.template operator()<selectNM_dNe>(dNe_midl);
                selectModeNM.template operator()<selectNM_dTe>(dTe_midl);
                selectModeNM.template operator()<selectNM_dP, ifReducedMHD>(dP_midl);

                runMHDRK4Stage.operator()<5>();

                if (checkNAN.template operator()<ifCheckNAN>()) {
                    if (myRank == 0) {
                        std::cout << BOLDRED << "NAN detected at " << currentStep
                                  << " steps. All accumulated data will be dumped and the program will exit now."
                                  << RESET << std::endl;
                        std::cout << std::endl;
                    }
                    goto finalize;
                }

                /*----------------------------Diagnostics & Output----------------------------*/

                diagIndex++;
                outputIndex++;

                runDiagnostics(diagIndex);
                runOutput(outputIndex);
            }

            /*----------------------------------PIC RK4-----------------------------------*/

            mhdToPIC.template operator()<ifPIC>();

            runPICRK4.template operator()<Ion, IonType, ifIon>(Ion_keys_in, Ion_values_in, globalPi, globalNi);
            runPICRK4.template operator()<Alpha, AlphaType, ifAlpha>(Alpha_keys_in, Alpha_values_in, globalPa,
                                                                     globalNa);
            runPICRK4.template operator()<Beam, BeamType, ifBeam>(Beam_keys_in, Beam_values_in, globalPb, globalNb);

            mergePICBuffers.template operator()<ifIon>(globalPi);
            mergePICBuffers.template operator()<ifAlpha>(globalPa);
            mergePICBuffers.template operator()<ifBeam>(globalPb);
            mergePICBuffers.template operator()<ifIon, ifQNeutrality>(globalNi);
            mergePICBuffers.template operator()<ifAlpha, ifQNeutrality>(globalNa);
            mergePICBuffers.template operator()<ifBeam, ifQNeutrality>(globalNb);

            ncclGroupStart();
            allReducePressure.template operator()<ifIon>(globalPi);
            allReducePressure.template operator()<ifAlpha>(globalPa);
            allReducePressure.template operator()<ifBeam>(globalPb);
            allReducePressure.template operator()<ifIon, ifQNeutrality>(globalNi);
            allReducePressure.template operator()<ifAlpha, ifQNeutrality>(globalNa);
            allReducePressure.template operator()<ifBeam, ifQNeutrality>(globalNb);
            ncclGroupEnd();

            updateAlignedGhost.template operator()<ifIon>(globalPi, dPi_midl);
            updateAlignedGhost.template operator()<ifAlpha>(globalPa, dPa_midl);
            updateAlignedGhost.template operator()<ifBeam>(globalPb, dPb_midl);
            updateAlignedGhost.template operator()<ifIon, ifQNeutrality>(globalNi, dNi_midl);
            updateAlignedGhost.template operator()<ifAlpha, ifQNeutrality>(globalNa, dNa_midl);
            updateAlignedGhost.template operator()<ifBeam, ifQNeutrality>(globalNb, dNb_midl);

            picToMHD.template operator()<ifPIC>();

            nablaPerp2.template operator()<ifNablaPerp2dP, ifIon>(dPConfigs, dPDatas, dPAs, dPXs, dPBs, dPi_midl,
                                                                  dPi_midr);
            nablaPerp2.template operator()<ifNablaPerp2dP, ifAlpha>(dPConfigs, dPDatas, dPAs, dPXs, dPBs, dPa_midl,
                                                                    dPa_midr);
            nablaPerp2.template operator()<ifNablaPerp2dP, ifBeam>(dPConfigs, dPDatas, dPAs, dPXs, dPBs, dPb_midl,
                                                                   dPb_midr);
            nablaPerp2.template operator()<ifNablaPerp2dNe, ifIon, ifQNeutrality>(dNeConfigs, dNeDatas, dNeAs, dNeXs,
                                                                                  dNeBs, dNi_midl, dNi_midr);
            nablaPerp2.template operator()<ifNablaPerp2dNe, ifAlpha, ifQNeutrality>(dNeConfigs, dNeDatas, dNeAs, dNeXs,
                                                                                    dNeBs, dNa_midl, dNa_midr);
            nablaPerp2.template operator()<ifNablaPerp2dNe, ifBeam, ifQNeutrality>(dNeConfigs, dNeDatas, dNeAs, dNeXs,
                                                                                   dNeBs, dNb_midl, dNb_midr);

            nablaPara4.template operator()<4, ifNablaPara4dP, ifIon>(dPi_midl, dPi_midr);
            nablaPara4.template operator()<4, ifNablaPara4dP, ifAlpha>(dPa_midl, dPa_midr);
            nablaPara4.template operator()<4, ifNablaPara4dP, ifBeam>(dPb_midl, dPb_midr);
            nablaPara4.template operator()<2, ifNablaPara4dNe, ifIon, ifQNeutrality>(dNi_midl, dNi_midr);
            nablaPara4.template operator()<2, ifNablaPara4dNe, ifAlpha, ifQNeutrality>(dNa_midl, dNa_midr);
            nablaPara4.template operator()<2, ifNablaPara4dNe, ifBeam, ifQNeutrality>(dNb_midl, dNb_midr);

            filterModeN.template operator()<ifFilterN_dP, ifIon>(dPi_midl, leftN, rightN);
            filterModeN.template operator()<ifFilterN_dP, ifAlpha>(dPa_midl, leftN, rightN);
            filterModeN.template operator()<ifFilterN_dP, ifBeam>(dPb_midl, leftN, rightN);
            filterModeN.template operator()<ifFilterN_dNe, ifIon, ifQNeutrality>(dNi_midl, leftN, rightN);
            filterModeN.template operator()<ifFilterN_dNe, ifAlpha, ifQNeutrality>(dNa_midl, leftN, rightN);
            filterModeN.template operator()<ifFilterN_dNe, ifBeam, ifQNeutrality>(dNb_midl, leftN, rightN);

            removeModeN.template operator()<removeN_dP, ifIon>(dPi_midl, dPi_midr);
            removeModeN.template operator()<removeN_dP, ifAlpha>(dPa_midl, dPa_midr);
            removeModeN.template operator()<removeN_dP, ifBeam>(dPb_midl, dPb_midr);
            removeModeN.template operator()<removeN_dNe, ifIon, ifQNeutrality>(dNi_midl, dNi_midr);
            removeModeN.template operator()<removeN_dNe, ifAlpha, ifQNeutrality>(dNa_midl, dNa_midr);
            removeModeN.template operator()<removeN_dNe, ifBeam, ifQNeutrality>(dNb_midl, dNb_midr);

            selectModeNM.template operator()<selectNM_dP, ifIon>(dPi_midl);
            selectModeNM.template operator()<selectNM_dP, ifAlpha>(dPa_midl);
            selectModeNM.template operator()<selectNM_dP, ifBeam>(dPb_midl);
            selectModeNM.template operator()<selectNM_dNe, ifIon, ifQNeutrality>(dNi_midl);
            selectModeNM.template operator()<selectNM_dNe, ifAlpha, ifQNeutrality>(dNa_midl);
            selectModeNM.template operator()<selectNM_dNe, ifBeam, ifQNeutrality>(dNb_midl);

            ncclGroupStart();
            haloExchange.template operator()<ifIon>(dPi_midl);
            haloExchange.template operator()<ifAlpha>(dPa_midl);
            haloExchange.template operator()<ifBeam>(dPb_midl);
            haloExchange.template operator()<ifIon, ifQNeutrality>(dNi_midl);
            haloExchange.template operator()<ifAlpha, ifQNeutrality>(dNa_midl);
            haloExchange.template operator()<ifBeam, ifQNeutrality>(dNb_midl);
            ncclGroupEnd();

            if constexpr (std::is_same_v<ifQNeutrality, trueType>) {
                forEachDev([&](int i) {
                    MHDQNeutrality2dNe<<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(w_midl[i], dNi_midl[i], dNa_midl[i],
                                                                              dNb_midl[i], dNe_midl[i]);
                    MHD2dJpBdPePhi<ifNonlinearMHD, ifLocal, falseType, ifReducedMHD>
                        <<<MRK4GridSize, MRK4BlockSize, 0, 0>>>(A_midl[i], dJpB_midl[i], A2dJpB[i], w_midl[i],
                                                                Phi_midl[i], w2Phi[i], dNe_midl[i], dTe_midl[i],
                                                                dP_midl[i], dPe_midl[i], Ne0[i], Te0[i]);
                });
                ncclGroupStart();
                haloExchange(dNe_midl, dJpB_midl, dPe_midl);
                ncclGroupEnd();
            }
        }

        sortParticles.template operator()<ifIon>(Ion_storage, Ion_storage_bytes, Ion_keys_in, Ion_keys_out,
                                                 Ion_sort_ids_in, Ion_sort_ids_out, Ion_values_in, Ion_values_out);
        sortParticles.template operator()<ifAlpha>(Alpha_storage, Alpha_storage_bytes, Alpha_keys_in, Alpha_keys_out,
                                                   Alpha_sort_ids_in, Alpha_sort_ids_out, Alpha_values_in,
                                                   Alpha_values_out);
        sortParticles.template operator()<ifBeam>(Beam_storage, Beam_storage_bytes, Beam_keys_in, Beam_keys_out,
                                                  Beam_sort_ids_in, Beam_sort_ids_out, Beam_values_in, Beam_values_out);
    }

    /*-----------------------------Phase 9: Finalize Timing & Memory------------------------------*/

    if (myRank == 0) {
        std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;
        std::cout << std::endl;
    }

finalize:

    forEachDev([&](int i) {
        cudaEventRecord(end[i]);
        cudaEventSynchronize(end[i]);
        cudaEventElapsedTime(&time[i], start[i], end[i]);
        CUDACHECK(cudaGetLastError());
    });

    if (myRank == 0) {
        for (int i = 1; i < devNums; i++)
            time[0] += time[i];
        time[0] /= devNums;
        if (time[0] > 1000)
            std::cout << BOLDGREEN << "Time used: " << std::setprecision(10) << time[0] / 1000 << "s." << RESET
                      << std::endl;
        else
            std::cout << BOLDGREEN << "Time used: " << std::setprecision(10) << time[0] << "ms." << RESET << std::endl;
        std::cout << std::endl;
    }

    cuGMEC.memcpyDeviceToHost(finalDir);

    cuGMEC.releaseDeviceMemory();
    cuGMEC.releaseHostMemory();

    /*---------------------------Phase 10: NCCL Destroy & MPI Finalize----------------------------*/

    // finalizing NCCL
    for (int i = 0; i < devNums; i++) {
        cudaSetDevice(localRank * devNums + i);
        ncclCommDestroy(comms[i]);
    }

    if (myRank == 0) {
        std::cout << BOLDYELLOW << "Start: Exit the program." << RESET << std::endl;
        std::cout << BOLDGREEN << "Done." << RESET << std::endl;
        std::cout << std::endl;
    }

    // finalizing MPI
    MPICHECK(MPI_Finalize());

    return 0;
}
