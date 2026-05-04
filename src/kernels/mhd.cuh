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

/*----------------------------------------MHD Kernels-----------------------------------------*/

template <int rk4, typename nonlinear, typename local, typename staggered, typename Eparallel, typename type>
__global__ void MHDLinearRK4(type* __restrict__ d_qtheta, type* __restrict__ w_beg, type* __restrict__ w_midl,
                             type* __restrict__ w_midr, type* __restrict__ w_end, type* __restrict__ A_beg,
                             type* __restrict__ A_midl, type* __restrict__ A_midr, type* __restrict__ A_end,
                             type* __restrict__ dNe_beg, type* __restrict__ dNe_midl, type* __restrict__ dNe_midr,
                             type* __restrict__ dNe_end, type* __restrict__ dTe_beg, type* __restrict__ dTe_midl,
                             type* __restrict__ dTe_midr, type* __restrict__ dTe_end, type* __restrict__ Phi_mid,
                             type* __restrict__ dJpB_mid, type* __restrict__ dPe_mid, type* __restrict__ dPi_mid,
                             type* __restrict__ dPa_mid, type* __restrict__ dPb_mid, type* __restrict__ d_wdPAdJpB2w,
                             type* __restrict__ d_APhidNe2A, type* __restrict__ d_dPePhiAdJpB2dNe,
                             type* __restrict__ d_PhidTedNe2dTe) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = (j + gridGhost) * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    type qtheta;
    type qtheta_lr[4];
    int shift_k;
    type shift_lk;
    type shift_dk;

    /*----------------------------------Field and Derivative----------------------------------*/

    type field;
    type field_px, field_py, field_pz;
    type field_du[4];
    type field_lr[4];

    /*---------------------------------Compressed Coefficient---------------------------------*/

    type compcoes[6];

    /*-------------------------------------RK4 Variables--------------------------------------*/

    type w_begin, dwdt;
    type A_begin, dAdt;
    type dNe_begin, dNedt;
    type dTe_begin, dTedt;

    /*---------------------------------------Initialize---------------------------------------*/

    w_begin = w_beg[offset3d];
    A_begin = A_beg[offset3d];
    dNe_begin = dNe_beg[offset3d];
    dTe_begin = dTe_beg[offset3d];

    qtheta = d_qtheta[offset2d];
    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    offset2d = j * gridNx + i;

    dwdt = 0;
    dAdt = 0;
    dNedt = 0;
    dTedt = 0;

    /*------------------Electron Pressure in Vorticity and Electron Density-------------------*/

    field = dPe_mid[offset3d];

    PartialZ<local>(k, offset3d, lane_id, dPe_mid, field, field_du, field_pz);
    PartialX(0, i, k, offset3d, dPe_mid, field, field_lr, field_px);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPe_mid, field_du, field_lr,
                    field_py);

    for (int index = 0; index < 3; index++)
        compcoes[index] = d_dPePhiAdJpB2dNe[offset2d * 11 + index];

    dNedt += compcoes[0] * field_px + compcoes[1] * field_py + compcoes[2] * field_pz;

    for (int index = 0; index < 5; index++)
        compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index];

    dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;

    /*---------------------------Ion Diamagnetic Drift in Vorticity---------------------------*/

    field = w_midl[offset3d];

    PartialZ<local>(k, offset3d, lane_id, w_midl, field, field_du, field_pz);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, w_midl, field_du, field_lr,
                    field_py);

    dwdt += compcoes[0] * field_py + compcoes[1] * field_pz;

    /*-------------------------------Ion Pressure in Vorticity--------------------------------*/

    if constexpr (std::is_same_v<ifIon, trueType>) {

        field = dPi_mid[offset3d];

        PartialZ<local>(k, offset3d, lane_id, dPi_mid, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, dPi_mid, field, field_lr, field_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPi_mid, field_du,
                        field_lr, field_py);

        dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;
    }

    /*------------------------------Alpha Pressure in Vorticity-------------------------------*/

    if constexpr (std::is_same_v<ifAlpha, trueType>) {

        field = dPa_mid[offset3d];

        PartialZ<local>(k, offset3d, lane_id, dPa_mid, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, dPa_mid, field, field_lr, field_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPa_mid, field_du,
                        field_lr, field_py);

        dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;
    }

    /*-------------------------------Beam Pressure in Vorticity-------------------------------*/

    if constexpr (std::is_same_v<ifBeam, trueType>) {

        field = dPb_mid[offset3d];

        PartialZ<local>(k, offset3d, lane_id, dPb_mid, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, dPb_mid, field, field_lr, field_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPb_mid, field_du,
                        field_lr, field_py);

        dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;
    }

    /*-----------------Electric Potential in Electron Density and Temperature-----------------*/

    field = Phi_mid[offset3d];

    PartialZ<local>(k, offset3d, lane_id, Phi_mid, field, field_du, field_pz);
    PartialX(0, i, k, offset3d, Phi_mid, field, field_lr, field_px);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr,
                    field_py);

    for (int index = 0; index < 3; index++)
        compcoes[index] = d_dPePhiAdJpB2dNe[offset2d * 11 + index + 3];

    dNedt += compcoes[0] * field_px + compcoes[1] * field_py + compcoes[2] * field_pz;

    if constexpr (std::is_same_v<nonlinear, falseType>) {

        for (int index = 0; index < 6; index++)
            compcoes[index] = d_PhidTedNe2dTe[offset2d * 6 + index];

        dTedt += compcoes[0] * field_py + compcoes[1] * field_pz;
    }

    /*----------------------Electron Temperature in Electron Temperature----------------------*/

    if constexpr (std::is_same_v<nonlinear, falseType>) {

        field = dTe_midl[offset3d];

        PartialZ<local>(k, offset3d, lane_id, dTe_midl, field, field_du, field_pz);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dTe_midl, field_du,
                        field_lr, field_py);

        dTedt += compcoes[2] * field_py + compcoes[3] * field_pz;
    }

    /*------------------------Electron Density in Electron Temperature------------------------*/

    if constexpr (std::is_same_v<nonlinear, falseType>) {

        field = dNe_midl[offset3d];

        PartialZ<local>(k, offset3d, lane_id, dNe_midl, field, field_du, field_pz);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl, field_du,
                        field_lr, field_py);

        dTedt += compcoes[4] * field_py + compcoes[5] * field_pz;
    }

    /*-----------------Parallel Vector Potential in Parallel Vector Potential-----------------*/

    for (int index = 0; index < 5; index++)
        compcoes[index] = d_APhidNe2A[offset2d * 5 + index];

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        if constexpr (std::is_same_v<nonlinear, falseType>) {

            field = A_midl[offset3d];
            PartialZ<local>(k, offset3d, lane_id, A_midl, field, field_du, field_pz);
            PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du,
                            field_lr, field_py);
            dAdt += compcoes[0] * field + compcoes[1] * field_py + compcoes[2] * field_pz;
        }
    }

    /*--------------------Electric Potential in Parallel Vector Potential---------------------*/

    if constexpr (std::is_same_v<staggered, trueType>) {

        Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, Phi_mid, field_du, field_lr, field);
        field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * mhdGridDy);

    } else {

        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du,
                        field_lr, field_py);
    }

    dAdt += compcoes[3] * field_py;

    /*---------------------Electron Density in Parallel Vector Potential----------------------*/

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        if constexpr (std::is_same_v<nonlinear, falseType>) {

            if constexpr (std::is_same_v<staggered, trueType>) {

                Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta,
                                    qtheta, qtheta_lr, dNe_midl, field_du, field_lr, field);
                field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * mhdGridDy);

            } else {

                PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl,
                                field_du, field_lr, field_py);
            }

            dAdt += compcoes[4] * field_py;
        }
    }

    /*--------------Parallel Vector Potential in Vorticity and Electron Density---------------*/

    if constexpr (std::is_same_v<staggered, trueType>) {

        S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                             qtheta_lr, A_midl, field_du, field_lr, field, field_px, field_py, field_pz);

        for (int index = 0; index < 4; index++)
            compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 5];

        dwdt += compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz;
        dNedt +=
            (compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz) / NormQE;

    } else {

        field = A_midl[offset3d];

        PartialZ<local>(k, offset3d, lane_id, A_midl, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, A_midl, field, field_lr, field_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du,
                        field_lr, field_py);

        for (int index = 0; index < 4; index++)
            compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 5];

        dwdt += compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz;
        dNedt +=
            (compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz) / NormQE;
    }

    /*-----------------------------Parallel Current in Vorticity------------------------------*/

    if constexpr (std::is_same_v<staggered, trueType>) {

        Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, dJpB_mid, field_du, field_lr, field);
        field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * mhdGridDy);

        for (int index = 0; index < 1; index++)
            compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 9];

        dwdt += compcoes[0] * field_py;
        dNedt += compcoes[0] * field_py / NormQE;

    } else {

        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dJpB_mid, field_du,
                        field_lr, field_py);

        for (int index = 0; index < 1; index++)
            compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 9];

        dwdt += compcoes[0] * field_py;
        dNedt += compcoes[0] * field_py / NormQE;
    }

    /*------------------------------------------RK4-------------------------------------------*/

    if constexpr (rk4 == 1) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] = w_begin + dwdt * mhdGridDt / 2;
            w_end[offset3d] = w_begin + dwdt * mhdGridDt / 6;

            A_midr[offset3d] = A_begin + dAdt * mhdGridDt / 2;
            A_end[offset3d] = A_begin + dAdt * mhdGridDt / 6;

            dNe_midr[offset3d] = dNe_begin + dNedt * mhdGridDt / 2;
            dNe_end[offset3d] = dNe_begin + dNedt * mhdGridDt / 6;

            dTe_midr[offset3d] = dTe_begin + dTedt * mhdGridDt / 2;
            dTe_end[offset3d] = dTe_begin + dTedt * mhdGridDt / 6;
        }

    } else if constexpr (rk4 == 2) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] = w_begin + dwdt * mhdGridDt / 2;
            w_end[offset3d] += dwdt * mhdGridDt / 3;

            A_midr[offset3d] = A_begin + dAdt * mhdGridDt / 2;
            A_end[offset3d] += dAdt * mhdGridDt / 3;

            dNe_midr[offset3d] = dNe_begin + dNedt * mhdGridDt / 2;
            dNe_end[offset3d] += dNedt * mhdGridDt / 3;

            dTe_midr[offset3d] = dTe_begin + dTedt * mhdGridDt / 2;
            dTe_end[offset3d] += dTedt * mhdGridDt / 3;
        }

    } else if constexpr (rk4 == 3) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] = w_begin + dwdt * mhdGridDt;
            w_end[offset3d] += dwdt * mhdGridDt / 3;

            A_midr[offset3d] = A_begin + dAdt * mhdGridDt;
            A_end[offset3d] += dAdt * mhdGridDt / 3;

            dNe_midr[offset3d] = dNe_begin + dNedt * mhdGridDt;
            dNe_end[offset3d] += dNedt * mhdGridDt / 3;

            dTe_midr[offset3d] = dTe_begin + dTedt * mhdGridDt;
            dTe_end[offset3d] += dTedt * mhdGridDt / 3;
        }

    } else if constexpr (rk4 == 4) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] = w_end[offset3d] + dwdt * mhdGridDt / 6;
            w_end[offset3d] += dwdt * mhdGridDt / 6;

            A_midr[offset3d] = A_end[offset3d] + dAdt * mhdGridDt / 6;
            A_end[offset3d] += dAdt * mhdGridDt / 6;

            dNe_midr[offset3d] = dNe_end[offset3d] + dNedt * mhdGridDt / 6;
            dNe_end[offset3d] += dNedt * mhdGridDt / 6;

            dTe_midr[offset3d] = dTe_end[offset3d] + dTedt * mhdGridDt / 6;
            dTe_end[offset3d] += dTedt * mhdGridDt / 6;
        }
    }
}

template <int rk4, typename MaxwellStress, typename ReynoldsStress, typename diagZFDrive, typename local,
          typename staggered, typename Eparallel, typename type>
__global__ void MHDNonlinearRK4(type* __restrict__ d_qtheta, type* __restrict__ w_midl, type* __restrict__ w_midr,
                                type* __restrict__ w_end, type* __restrict__ A_midl, type* __restrict__ A_midr,
                                type* __restrict__ A_end, type* __restrict__ dNe_midl, type* __restrict__ dNe_midr,
                                type* __restrict__ dNe_end, type* __restrict__ dTe_midl, type* __restrict__ dTe_midr,
                                type* __restrict__ dTe_end, type* __restrict__ Phi_mid, type* __restrict__ dJpB_mid,
                                type* __restrict__ dPe_mid, type* __restrict__ d_Ne0, type* __restrict__ d_Te0,
                                type* __restrict__ d_Ne0_px, type* __restrict__ d_Te0_px,
                                type* __restrict__ d_APhidNe2A, type* __restrict__ d_wPhi_w,
                                type* __restrict__ d_AdJpB_w, type* __restrict__ d_PhiA_A, type* __restrict__ d_NeA_A,
                                type* __restrict__ d_dNePhi_dNe, type* __restrict__ d_PhiTe_dTe,
                                type* __restrict__ d_PhiTeA_dTe, type* __restrict__ d_Maxwell,
                                type* __restrict__ d_Reynolds) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = j * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    type qtheta;
    type qtheta_lr[4];
    int shift_k;
    type shift_lk;
    type shift_dk;

    /*----------------------------------Field and Derivative----------------------------------*/

    type Ne0, Te0, Ne0_px, Te0_px, Pe0_px, Ne, Te;
    type A, A_px, A_py, A_pz;
    type Phi, Phi_px, Phi_py, Phi_pz;
    type field, field_px, field_py, field_pz;
    type field_du[4];
    type field_lr[4];

    /*---------------------------------Compressed Coefficient---------------------------------*/

    type compcoes[9];

    /*-------------------------------------RK4 Variables--------------------------------------*/

    type dwdt, dAdt, dNedt, dTedt;
    type MaxwellDrive, ReynoldsDrive, dwdtBefore;

    /*---------------------------------------Initialize---------------------------------------*/

    Ne0 = d_Ne0[offset2d];
    Te0 = d_Te0[offset2d];
    Ne0_px = d_Ne0_px[offset2d];
    Te0_px = d_Te0_px[offset2d];
    Pe0_px = Ne0_px * Te0 + Ne0 * Te0_px;

    offset2d = (j + gridGhost) * gridNx + i;

    qtheta = d_qtheta[offset2d];
    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    offset2d = j * gridNx + i;

    dwdt = 0;
    dAdt = 0;
    dNedt = 0;
    dTedt = 0;
    MaxwellDrive = 0;
    ReynoldsDrive = 0;

    Phi = Phi_mid[offset3d];
    PartialZ<local>(k, offset3d, lane_id, Phi_mid, Phi, field_du, Phi_pz);
    PartialX(0, i, k, offset3d, Phi_mid, Phi, field_lr, Phi_px);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr,
                    Phi_py);

    if constexpr (std::is_same_v<staggered, falseType>) {

        A = A_midl[offset3d];
        PartialZ<local>(k, offset3d, lane_id, A_midl, A, field_du, A_pz);
        PartialX(0, i, k, offset3d, A_midl, A, field_lr, A_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du,
                        field_lr, A_py);

        field = dJpB_mid[offset3d];
        PartialZ<local>(k, offset3d, lane_id, dJpB_mid, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, dJpB_mid, field, field_lr, field_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dJpB_mid, field_du,
                        field_lr, field_py);

    } else {

        S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                             qtheta_lr, A_midl, field_du, field_lr, A, A_px, A_py, A_pz);

        S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                             qtheta_lr, dJpB_mid, field_du, field_lr, field, field_px, field_py, field_pz);

        offset2d = (j + gridGhost) * gridNx + i;

        qtheta = d_qtheta[offset2d];
        qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
        qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
        qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
        qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

        offset2d = j * gridNx + i;
    }

    /*-----------------------------------AdJpB in Vorticity-----------------------------------*/

    for (int index = 0; index < 9; index++)
        compcoes[index] = d_AdJpB_w[offset2d * 9 + index];

    dwdtBefore = dwdt;

    if constexpr (std::is_same_v<MaxwellStress, trueType>)
        dwdt += MaxwellStressCoef *
                (compcoes[0] * A * field_px + compcoes[1] * A * field_py + compcoes[2] * A * field_pz +
                 compcoes[3] * A_px * field_py + compcoes[4] * A_px * field_pz + compcoes[5] * A_py * field_px +
                 compcoes[6] * A_py * field_pz + compcoes[7] * A_pz * field_px + compcoes[8] * A_pz * field_py);

    MaxwellDrive = dwdt - dwdtBefore;

    /*-------------------------------AdJpB in Electron Density--------------------------------*/

    dNedt += (compcoes[0] * A * field_px + compcoes[1] * A * field_py + compcoes[2] * A * field_pz +
              compcoes[3] * A_px * field_py + compcoes[4] * A_px * field_pz + compcoes[5] * A_py * field_px +
              compcoes[6] * A_py * field_pz + compcoes[7] * A_pz * field_px + compcoes[8] * A_pz * field_py) /
             NormQE;

    /*-----------------------------------wPhi in Vorticity------------------------------------*/

    for (int index = 0; index < 6; index++)
        compcoes[index] = d_wPhi_w[offset2d * 6 + index];

    field = w_midl[offset3d];
    PartialZ<local>(k, offset3d, lane_id, w_midl, field, field_du, field_pz);
    PartialX(0, i, k, offset3d, w_midl, field, field_lr, field_px);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, w_midl, field_du, field_lr,
                    field_py);

    dwdtBefore = dwdt;

    if constexpr (std::is_same_v<ReynoldsStress, trueType>)
        dwdt += ReynoldsStressCoef *
                (compcoes[0] * field_px * Phi_py + compcoes[1] * field_px * Phi_pz + compcoes[2] * field_py * Phi_px +
                 compcoes[3] * field_py * Phi_pz + compcoes[4] * field_pz * Phi_px + compcoes[5] * field_pz * Phi_py);

    ReynoldsDrive = dwdt - dwdtBefore;

    /*-------------------------------dNePhi in Electron Density-------------------------------*/

    for (int index = 0; index < 9; index++)
        compcoes[index] = d_dNePhi_dNe[offset2d * 9 + index];

    field = dNe_midl[offset3d];
    PartialZ<local>(k, offset3d, lane_id, dNe_midl, field, field_du, field_pz);
    PartialX(0, i, k, offset3d, dNe_midl, field, field_lr, field_px);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl, field_du, field_lr,
                    field_py);

    dNedt += compcoes[0] * field * Phi_px + compcoes[1] * field * Phi_py + compcoes[2] * field * Phi_pz +
             compcoes[3] * field_px * Phi_py + compcoes[4] * field_px * Phi_pz + compcoes[5] * field_py * Phi_px +
             compcoes[6] * field_py * Phi_pz + compcoes[7] * field_pz * Phi_px + compcoes[8] * field_pz * Phi_py;

    Ne = Ne0 + field;

    /*-----------------------------PhiTeA in Electron Temperature-----------------------------*/

    field = dTe_midl[offset3d];
    PartialZ<local>(k, offset3d, lane_id, dTe_midl, field, field_du, field_pz);
    PartialX(0, i, k, offset3d, dTe_midl, field, field_lr, field_px);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dTe_midl, field_du, field_lr,
                    field_py);
    field_px += Te0_px;

    Te = Te0 + field;

    for (int index = 0; index < 6; index++)
        compcoes[index] = d_PhiTe_dTe[offset2d * 6 + index];

    dTedt += compcoes[0] * Phi_px * field_py + compcoes[1] * Phi_px * field_pz + compcoes[2] * Phi_py * field_px +
             compcoes[3] * Phi_py * field_pz + compcoes[4] * Phi_pz * field_px + compcoes[5] * Phi_pz * field_py;

    for (int index = 0; index < 9; index++)
        compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index];

    dTedt += compcoes[0] * Phi_px * field_py * A + compcoes[1] * Phi_px * field_py * A_px +
             compcoes[2] * Phi_px * field_py * A_py + compcoes[3] * Phi_px * field_py * A_pz;
    dTedt += compcoes[4] * Phi_px * field_pz * A;
    dTedt += compcoes[5] * Phi_py * field_px * A + compcoes[6] * Phi_py * field_px * A_px +
             compcoes[7] * Phi_py * field_px * A_py + compcoes[8] * Phi_py * field_px * A_pz;

    for (int index = 0; index < 9; index++)
        compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index + 9];

    dTedt += compcoes[0] * Phi_py * field_pz * A + compcoes[1] * Phi_py * field_pz * A_px +
             compcoes[2] * Phi_py * field_pz * A_py + compcoes[3] * Phi_py * field_pz * A_pz;
    dTedt += compcoes[4] * Phi_pz * field_px * A;
    dTedt += compcoes[5] * Phi_pz * field_py * A + compcoes[6] * Phi_pz * field_py * A_px +
             compcoes[7] * Phi_pz * field_py * A_py + compcoes[8] * Phi_pz * field_py * A_pz;

    /*-----------------------------PeTeA in Electron Temperature------------------------------*/

    Phi = dPe_mid[offset3d];
    PartialZ<local>(k, offset3d, lane_id, dPe_mid, Phi, field_du, Phi_pz);
    PartialX(0, i, k, offset3d, dPe_mid, Phi, field_lr, Phi_px);
    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPe_mid, field_du, field_lr,
                    Phi_py);
    Phi_px += Pe0_px;

    for (int index = 0; index < 6; index++)
        compcoes[index] = d_PhiTe_dTe[offset2d * 6 + index];

    dTedt += -(compcoes[0] * Phi_px * field_py + compcoes[1] * Phi_px * field_pz + compcoes[2] * Phi_py * field_px +
               compcoes[3] * Phi_py * field_pz + compcoes[4] * Phi_pz * field_px + compcoes[5] * Phi_pz * field_py) /
             (2 * NormQE * Ne);

    for (int index = 0; index < 9; index++)
        compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index];

    dTedt += -(compcoes[0] * Phi_px * field_py * A + compcoes[1] * Phi_px * field_py * A_px +
               compcoes[2] * Phi_px * field_py * A_py + compcoes[3] * Phi_px * field_py * A_pz) /
             (2 * NormQE * Ne);
    dTedt += -(compcoes[4] * Phi_px * field_pz * A) / (2 * NormQE * Ne);
    dTedt += -(compcoes[5] * Phi_py * field_px * A + compcoes[6] * Phi_py * field_px * A_px +
               compcoes[7] * Phi_py * field_px * A_py + compcoes[8] * Phi_py * field_px * A_pz) /
             (2 * NormQE * Ne);

    for (int index = 0; index < 9; index++)
        compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index + 9];

    dTedt += -(compcoes[0] * Phi_py * field_pz * A + compcoes[1] * Phi_py * field_pz * A_px +
               compcoes[2] * Phi_py * field_pz * A_py + compcoes[3] * Phi_py * field_pz * A_pz) /
             (2 * NormQE * Ne);
    dTedt += -(compcoes[4] * Phi_pz * field_px * A) / (2 * NormQE * Ne);
    dTedt += -(compcoes[5] * Phi_pz * field_py * A + compcoes[6] * Phi_pz * field_py * A_px +
               compcoes[7] * Phi_pz * field_py * A_py + compcoes[8] * Phi_pz * field_py * A_pz) /
             (2 * NormQE * Ne);

    /*---------------------------PhiA in Parallel Vector Potential----------------------------*/

    for (int index = 0; index < 9; index++)
        compcoes[index] = d_PhiA_A[offset2d * 9 + index];

    if constexpr (std::is_same_v<staggered, trueType>) {

        A = A_midl[offset3d];
        PartialZ<local>(k, offset3d, lane_id, A_midl, A, field_du, A_pz);
        PartialX(0, i, k, offset3d, A_midl, A, field_lr, A_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du,
                        field_lr, A_py);

        C2SPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                             qtheta_lr, Phi_mid, field_du, field_lr, Phi, Phi_px, Phi_py, Phi_pz);

    } else {

        Phi = Phi_mid[offset3d];
        PartialZ<local>(k, offset3d, lane_id, Phi_mid, Phi, field_du, Phi_pz);
        PartialX(0, i, k, offset3d, Phi_mid, Phi, field_lr, Phi_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du,
                        field_lr, Phi_py);
    }

    dAdt += compcoes[0] * Phi_px * A + compcoes[1] * Phi_px * A_py + compcoes[2] * Phi_px * A_pz +
            compcoes[3] * Phi_py * A + compcoes[4] * Phi_py * A_px + compcoes[5] * Phi_py * A_pz +
            compcoes[6] * Phi_pz * A + compcoes[7] * Phi_pz * A_px + compcoes[8] * Phi_pz * A_py;

    /*----------------------------dNe in Parallel Vector Potential----------------------------*/

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        compcoes[0] = d_APhidNe2A[offset2d * 5 + 4];

        if constexpr (std::is_same_v<staggered, falseType>) {

            field = dNe_midl[offset3d];
            PartialZ<local>(k, offset3d, lane_id, dNe_midl, field, field_du, field_pz);
            PartialX(0, i, k, offset3d, dNe_midl, field, field_lr, field_px);
            PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl, field_du,
                            field_lr, field_py);

        } else {

            C2SPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                                 qtheta_lr, dNe_midl, field_du, field_lr, field, field_px, field_py, field_pz);

            Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta,
                                qtheta, qtheta_lr, dTe_midl, field_du, field_lr, Te);

            Ne = Ne0 + field;
            Te = Te0 + Te;
        }

        dAdt += compcoes[0] * (Te / Te0 * Ne0 / Ne) * field_py;
    }

    /*----------------------------NeA in Parallel Vector Potential----------------------------*/

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        for (int index = 0; index < 9; index++)
            compcoes[index] = d_NeA_A[offset2d * 9 + index];

        field_px += Ne0_px;

        dAdt += (compcoes[0] * field_px * A + compcoes[1] * field_px * A_py + compcoes[2] * field_px * A_pz +
                 compcoes[3] * field_py * A + compcoes[4] * field_py * A_px + compcoes[5] * field_py * A_pz +
                 compcoes[6] * field_pz * A + compcoes[7] * field_pz * A_px + compcoes[8] * field_pz * A_py) *
                (Te / Te0 * Ne0 / Ne);
    }

    /*------------------------------------------RK4-------------------------------------------*/

    if constexpr (rk4 == 1) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] += dwdt * mhdGridDt / 2;
            w_end[offset3d] += dwdt * mhdGridDt / 6;

            A_midr[offset3d] += dAdt * mhdGridDt / 2;
            A_end[offset3d] += dAdt * mhdGridDt / 6;

            dNe_midr[offset3d] += dNedt * mhdGridDt / 2;
            dNe_end[offset3d] += dNedt * mhdGridDt / 6;

            dTe_midr[offset3d] += dTedt * mhdGridDt / 2;
            dTe_end[offset3d] += dTedt * mhdGridDt / 6;
        }

    } else if constexpr (rk4 == 2) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] += dwdt * mhdGridDt / 2;
            w_end[offset3d] += dwdt * mhdGridDt / 3;

            A_midr[offset3d] += dAdt * mhdGridDt / 2;
            A_end[offset3d] += dAdt * mhdGridDt / 3;

            dNe_midr[offset3d] += dNedt * mhdGridDt / 2;
            dNe_end[offset3d] += dNedt * mhdGridDt / 3;

            dTe_midr[offset3d] += dTedt * mhdGridDt / 2;
            dTe_end[offset3d] += dTedt * mhdGridDt / 3;
        }

    } else if constexpr (rk4 == 3) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] += dwdt * mhdGridDt;
            w_end[offset3d] += dwdt * mhdGridDt / 3;

            A_midr[offset3d] += dAdt * mhdGridDt;
            A_end[offset3d] += dAdt * mhdGridDt / 3;

            dNe_midr[offset3d] += dNedt * mhdGridDt;
            dNe_end[offset3d] += dNedt * mhdGridDt / 3;

            dTe_midr[offset3d] += dTedt * mhdGridDt;
            dTe_end[offset3d] += dTedt * mhdGridDt / 3;
        }

    } else if constexpr (rk4 == 4) {

        if (i != 0 && i != gridNx - 1) {
            w_midr[offset3d] += dwdt * mhdGridDt / 6;
            w_end[offset3d] += dwdt * mhdGridDt / 6;

            A_midr[offset3d] += dAdt * mhdGridDt / 6;
            A_end[offset3d] += dAdt * mhdGridDt / 6;

            dNe_midr[offset3d] += dNedt * mhdGridDt / 6;
            dNe_end[offset3d] += dNedt * mhdGridDt / 6;

            dTe_midr[offset3d] += dTedt * mhdGridDt / 6;
            dTe_end[offset3d] += dTedt * mhdGridDt / 6;
        }
    }

    if constexpr (std::is_same_v<diagZFDrive, trueType>) {
        if (i != 0 && i != gridNx - 1) {
            if constexpr (rk4 == 1) {
                d_Maxwell[offset3d] = MaxwellDrive * mhdGridDt / 6;
                d_Reynolds[offset3d] = ReynoldsDrive * mhdGridDt / 6;
            } else if constexpr (rk4 == 2) {
                d_Maxwell[offset3d] += MaxwellDrive * mhdGridDt / 3;
                d_Reynolds[offset3d] += ReynoldsDrive * mhdGridDt / 3;
            } else if constexpr (rk4 == 3) {
                d_Maxwell[offset3d] += MaxwellDrive * mhdGridDt / 3;
                d_Reynolds[offset3d] += ReynoldsDrive * mhdGridDt / 3;
            } else if constexpr (rk4 == 4) {
                d_Maxwell[offset3d] += MaxwellDrive * mhdGridDt / 6;
                d_Reynolds[offset3d] += ReynoldsDrive * mhdGridDt / 6;
            }
        }
    }
}

template <typename nonlinear, typename local, typename FLRMHD, typename type>
__global__ void MHD2dJpBdPePhi(type* __restrict__ A_mid, type* __restrict__ dJpB_mid, type* __restrict__ d_A2dJpB,
                               type* __restrict__ w_mid, type* __restrict__ Phi_mid, type* __restrict__ d_w2Phi,
                               type* __restrict__ dNe_mid, type* __restrict__ dTe_mid, type* __restrict__ dPe_mid,
                               type* __restrict__ d_Ne0, type* __restrict__ d_Te0) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = j * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------Field and Derivative----------------------------------*/

    type dJpB, Ne0, Te0, dNe, dTe;

    type field;
    type field_px, field_pz, field_px2, field_pxz, field_pz2;
    type field_du[4];
    type field_lr[4];

    /*---------------------------------Compressed Coefficient---------------------------------*/

    type compcoes[6];

    /*---------------------------------------Initialize---------------------------------------*/

    for (int index = 0; index < 6; index++)
        compcoes[index] = d_A2dJpB[offset2d * 6 + index];

    Ne0 = d_Ne0[offset2d];
    Te0 = d_Te0[offset2d];
    dNe = dNe_mid[offset3d];
    dTe = dTe_mid[offset3d];

    /*---------------------Parallel Vector Potential in Parallel Current----------------------*/

    field = A_mid[offset3d];

    PartialZ2<local>(k, offset3d, lane_id, A_mid, field, field_du, field_pz, field_pz2);
    PartialX2(i, offset3d, A_mid, field, field_lr, field_px, field_px2);

    dJpB = compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_pz + compcoes[3] * field_px2 +
           compcoes[5] * field_pz2;

    if constexpr (std::is_same_v<local, trueType>) {

        PartialZ<local>(k, offset3d, lane_id, A_mid, field_px, field_du, field_pxz);

    } else {

        field = A_mid[offset3d - k + (k - 2 + gridNz) % gridNz];
        PartialX(-2, i, k, offset3d, A_mid, field, field_lr, field_px);
        field_pxz = field_px;

        field = A_mid[offset3d - k + (k - 1 + gridNz) % gridNz];
        PartialX(-1, i, k, offset3d, A_mid, field, field_lr, field_px);
        field_pxz += -8 * field_px;

        field = A_mid[offset3d - k + (k + 1 + gridNz) % gridNz];
        PartialX(1, i, k, offset3d, A_mid, field, field_lr, field_px);
        field_pxz += 8 * field_px;

        field = A_mid[offset3d - k + (k + 2 + gridNz) % gridNz];
        PartialX(2, i, k, offset3d, A_mid, field, field_lr, field_px);
        field_pxz += -field_px;

        field_pxz /= (12 * mhdGridDz);
    }

    dJpB += compcoes[4] * field_pxz;

    if (i != 0 && i != gridNx - 1)
        dJpB_mid[offset3d] = dJpB;

    /*---------------------------------------FLR Effect---------------------------------------*/

    if constexpr (std::is_same_v<FLRMHD, trueType>)
        Phi_mid[offset3d] += d_w2Phi[offset2d] * w_mid[offset3d];

    /*-----------------------Electron Density, Temperature and Pressure-----------------------*/

    if constexpr (std::is_same_v<nonlinear, falseType>)
        dPe_mid[offset3d] = dNe * Te0 + Ne0 * dTe;
    else
        dPe_mid[offset3d] = dNe * Te0 + Ne0 * dTe + dNe * dTe;
}

template <typename local, typename type>
__global__ void MHD2w(type* __restrict__ Phi_mid, type* __restrict__ w_mid, type* __restrict__ d_Phi2w) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = j * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------Field and Derivative----------------------------------*/

    type w;

    type field;
    type field_px, field_pz, field_px2, field_pxz, field_pz2;
    type field_du[4];
    type field_lr[4];

    /*---------------------------------Compressed Coefficient---------------------------------*/

    type compcoes[5];

    /*---------------------------------------Initialize---------------------------------------*/

    for (int index = 0; index < 5; index++)
        compcoes[index] = d_Phi2w[offset2d * 5 + index];

    /*----------------------------Electric Potential in Vorticity-----------------------------*/

    field = Phi_mid[offset3d];

    PartialZ2<local>(k, offset3d, lane_id, Phi_mid, field, field_du, field_pz, field_pz2);
    PartialX2(i, offset3d, Phi_mid, field, field_lr, field_px, field_px2);

    w = compcoes[0] * field_px + compcoes[1] * field_pz + compcoes[2] * field_px2 + compcoes[4] * field_pz2;

    if constexpr (std::is_same_v<local, trueType>) {

        PartialZ<local>(k, offset3d, lane_id, Phi_mid, field_px, field_du, field_pxz);

    } else {

        field = Phi_mid[offset3d - k + (k - 2 + gridNz) % gridNz];
        PartialX(-2, i, k, offset3d, Phi_mid, field, field_lr, field_px);
        field_pxz = field_px;

        field = Phi_mid[offset3d - k + (k - 1 + gridNz) % gridNz];
        PartialX(-1, i, k, offset3d, Phi_mid, field, field_lr, field_px);
        field_pxz += -8 * field_px;

        field = Phi_mid[offset3d - k + (k + 1 + gridNz) % gridNz];
        PartialX(1, i, k, offset3d, Phi_mid, field, field_lr, field_px);
        field_pxz += 8 * field_px;

        field = Phi_mid[offset3d - k + (k + 2 + gridNz) % gridNz];
        PartialX(2, i, k, offset3d, Phi_mid, field, field_lr, field_px);
        field_pxz += -field_px;

        field_pxz /= (12 * mhdGridDz);
    }

    w += compcoes[3] * field_pxz;

    if (i != 0 && i != gridNx - 1)
        w_mid[offset3d] = w;
}

template <typename dirichlet0, typename dirichlet1, typename... types>
__global__ void MHDBoundary(types* __restrict__... fields) {

    Boundary<dirichlet0, dirichlet1>(fields...);
}

template <typename local, typename type, typename... types>
__global__ void MHDAlignedGhost(type* __restrict__ d_qtheta, types* __restrict__... fields) {

    AlignedGhost<local>(d_qtheta, fields...);
}

template <typename local, typename type, typename... types>
__global__ void MHDStaggered2C(type* __restrict__ d_qtheta, types* __restrict__... fields) {

    Staggered2C<local>(d_qtheta, fields...);
}

template <int dir, typename local, typename type, typename... types>
__global__ void MHDShifted2A(type* __restrict__ d_qtheta, types* __restrict__... fields) {

    Shifted2A<dir, local>(d_qtheta, fields...);
}

template <typename local, typename type>
__global__ void MHDNablaPara2(type* __restrict__ d_qtheta, type* __restrict__ d_field,
                              type* __restrict__ d_nablaPara2) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = (j + gridGhost) * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    type qtheta;
    type qtheta_lr[4];
    int shift_k;
    type shift_lk;
    type shift_dk;

    /*----------------------------------Field and Derivative----------------------------------*/

    type field;
    type field_py;
    type field_du[4];
    type field_lr[4];

    type nablaPara2;

    /*---------------------------------------Initialize---------------------------------------*/

    qtheta = d_qtheta[offset2d];
    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    /*---------------------------------------NablaPara2---------------------------------------*/

    field = d_field[offset3d];

    PartialY2<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, d_field, field, field_du,
                     field_lr, field_py, nablaPara2);

    d_nablaPara2[offset3d] = nablaPara2;
}

template <int F, typename local, typename type>
__global__ void MHDNablaPara4(type* __restrict__ d_qtheta, type* __restrict__ d_field,
                              type* __restrict__ d_nablaPara2) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = (j + gridGhost) * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    type qtheta;
    type qtheta_lr[4];
    int shift_k;
    type shift_lk;
    type shift_dk;

    /*----------------------------------Field and Derivative----------------------------------*/

    type field;
    type field_py;
    type field_du[4];
    type field_lr[4];

    type para4;
    type nablaPara4;

    /*---------------------------------------Initialize---------------------------------------*/

    qtheta = d_qtheta[offset2d];
    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    if constexpr (F == 0)
        para4 = para4Phi;
    else if constexpr (F == 1)
        para4 = para4A;
    else if constexpr (F == 2)
        para4 = para4dNe;
    else if constexpr (F == 3)
        para4 = para4dTe;
    else if constexpr (F == 4)
        para4 = para4dPi;
    else if constexpr (F == 5)
        para4 = para4dPa;
    else if constexpr (F == 6)
        para4 = para4dPb;

    /*---------------------------------------NablaPara4---------------------------------------*/

    field = d_nablaPara2[offset3d];

    PartialY2<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, d_nablaPara2, field,
                     field_du, field_lr, field_py, nablaPara4);

    d_field[offset3d] -= mhdGridDt * para4 * nablaPara4;
}

template <typename type>
__global__ void MHDGaussianAverage(type* __restrict__ field_in, type* __restrict__ field_out, type w0, type w1,
                                   type w2) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    /*------------------------------------5pt Gaussian in y-----------------------------------*/

    field_out[offset3d] = w0 * field_in[offset3d] +
                          w1 * (field_in[offset3d - 1 * gridNxz] + field_in[offset3d + 1 * gridNxz]) +
                          w2 * (field_in[offset3d - 2 * gridNxz] + field_in[offset3d + 2 * gridNxz]);
}

template <typename cufftType>
__global__ void MHDFilterModeN(cufftType* __restrict__ d_freq, int modeNumber) {

    /*-------------------------------------Related Index--------------------------------------*/

    int offset = (blockIdx.x * blockDim.x + threadIdx.x) * nFFTFreqSize;

    /*-----------------------------------------Filter-----------------------------------------*/

    for (int mode = 0; mode < nFFTFreqSize; mode++) {
        if (mode != modeNumber) {

            d_freq[offset + mode].x = 0;
            d_freq[offset + mode].y = 0;
        }
    }
}

template <typename cufftType>
__global__ void MHDFilterModeN(cufftType* __restrict__ d_freq, int leftModeNumber, int rightModeNumber) {

    /*-------------------------------------Related Index--------------------------------------*/

    int offset = (blockIdx.x * blockDim.x + threadIdx.x) * nFFTFreqSize;

    /*-----------------------------------------Filter-----------------------------------------*/

    for (int mode = 0; mode < nFFTFreqSize; mode++) {
        if (mode < leftModeNumber || mode > rightModeNumber) {

            d_freq[offset + mode].x = 0;
            d_freq[offset + mode].y = 0;
        }
    }
}

template <typename type>
__global__ void MHDFilterResizeN(type* __restrict__ d_field) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    /*-----------------------------------------Resize-----------------------------------------*/

    d_field[offset3d] /= gridNz;
}

template <typename cufftType>
__global__ void MHDFilterModeM(cufftType* __restrict__ d_freq, int modeNumber) {

    /*-------------------------------------Related Index--------------------------------------*/

    int offset = (blockIdx.x * blockDim.x + threadIdx.x) * mFFTFreqSize;

    /*-----------------------------------------Filter-----------------------------------------*/

    for (int mode = 0; mode < mFFTFreqSize; mode++) {
        if (mode != modeNumber) {

            d_freq[offset + mode].x = 0;
            d_freq[offset + mode].y = 0;
        }
    }
}

template <typename cufftType>
__global__ void MHDFilterModeM(cufftType* __restrict__ d_freq, int leftModeNumber, int rightModeNumber) {

    /*-------------------------------------Related Index--------------------------------------*/

    int offset = (blockIdx.x * blockDim.x + threadIdx.x) * mFFTFreqSize;

    /*-----------------------------------------Filter-----------------------------------------*/

    for (int mode = 0; mode < mFFTFreqSize; mode++) {
        if (mode < leftModeNumber || mode > rightModeNumber) {

            d_freq[offset + mode].x = 0;
            d_freq[offset + mode].y = 0;
        }
    }
}

template <typename type>
__global__ void MHDFilterResizeM(type* __restrict__ d_field) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    /*-----------------------------------------Resize-----------------------------------------*/

    d_field[offset3d] /= gridNy;
}

template <typename type>
__global__ void MHDTransposeLeft(type* __restrict__ d_yxzField, type* __restrict__ d_xzyField) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = threadIdx.x;
    int j = blockIdx.x;

    /*-------------------------------------Transpose Left-------------------------------------*/

    for (int k = 0; k < gridNz; k++)
        d_xzyField[i * gridNz * gridNy + k * gridNy + j] = d_yxzField[j * gridNxz + i * gridNz + k];
}

template <typename type>
__global__ void MHDTransposeRight(type* __restrict__ d_xzyField, type* __restrict__ d_yxzField) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x;
    int k = threadIdx.x;

    /*------------------------------------Transpose Right-------------------------------------*/

    for (int j = 0; j < gridNy; j++)
        d_yxzField[j * gridNxz + i * gridNz + k] = d_xzyField[i * gridNz * gridNy + k * gridNy + j];
}

template <typename type>
__global__ void MHDAddMode(type* __restrict__ d_Addend, type* __restrict__ d_Augend) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    /*-----------------------------------------Add N------------------------------------------*/

    d_Augend[offset3d] += d_Addend[offset3d];
}

template <typename type>
__global__ void MHDSubtractMode(type* __restrict__ d_Subtrahend, type* __restrict__ d_Minuend) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    /*---------------------------------------Subtract N---------------------------------------*/

    d_Minuend[offset3d] -= d_Subtrahend[offset3d];
}

template <typename nonlinear, typename local, typename staggered, typename Eparallel, typename type>
__global__ void MHD2Apt(type* __restrict__ d_qtheta, type* __restrict__ A_mid, type* __restrict__ dNe_mid,
                        type* __restrict__ dTe_mid, type* __restrict__ Phi_mid, type* __restrict__ d_Ne0,
                        type* __restrict__ d_Te0, type* __restrict__ d_Ne0_px, type* __restrict__ d_APhidNe2A,
                        type* __restrict__ d_PhiA_A, type* __restrict__ d_NeA_A, type* __restrict__ d_A_pt) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = (j + gridGhost) * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    type qtheta;
    type qtheta_lr[4];
    int shift_k;
    type shift_lk;
    type shift_dk;

    /*----------------------------------Field and Derivative----------------------------------*/

    type field;
    type field_px, field_py, field_pz;
    type field_du[4];
    type field_lr[4];

    /*---------------------------------Compressed Coefficient---------------------------------*/

    type compcoes[9];

    /*-------------------------------------RK4 Variables--------------------------------------*/

    type dAdt;

    /*---------------------------------------Initialize---------------------------------------*/

    qtheta = d_qtheta[offset2d];
    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    offset2d = j * gridNx + i;

    dAdt = 0;

    /*--------------------Electric Potential in Parallel Vector Potential---------------------*/

    for (int index = 0; index < 5; index++)
        compcoes[index] = d_APhidNe2A[offset2d * 5 + index];

    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr,
                    field_py);
    dAdt += compcoes[3] * field_py;

    /*---------------------Electron Density in Parallel Vector Potential----------------------*/

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        if constexpr (std::is_same_v<nonlinear, falseType>) {

            PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid, field_du,
                            field_lr, field_py);
            dAdt += compcoes[4] * field_py;
        }
    }

    /*-----------------Parallel Vector Potential in Parallel Vector Potential-----------------*/

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        if constexpr (std::is_same_v<nonlinear, falseType>) {

            if constexpr (std::is_same_v<staggered, falseType>) {

                field = A_mid[offset3d];
                PartialZ<local>(k, offset3d, lane_id, A_mid, field, field_du, field_pz);
                PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du,
                                field_lr, field_py);

            } else {

                S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta,
                                     qtheta, qtheta_lr, A_mid, field_du, field_lr, field, field_px, field_py, field_pz);
            }

            dAdt += compcoes[0] * field + compcoes[1] * field_py + compcoes[2] * field_pz;
        }
    }

    if constexpr (std::is_same_v<nonlinear, trueType>) {

        type Ne0, Te0, Ne0_px, Ne, Te;
        type A, A_px, A_py, A_pz;
        type Phi, Phi_px, Phi_py, Phi_pz;

        Ne0 = d_Ne0[offset2d];
        Te0 = d_Te0[offset2d];
        Ne0_px = d_Ne0_px[offset2d];
        Te = Te0 + dTe_mid[offset3d];

        Phi = Phi_mid[offset3d];
        PartialZ<local>(k, offset3d, lane_id, Phi_mid, Phi, field_du, Phi_pz);
        PartialX(0, i, k, offset3d, Phi_mid, Phi, field_lr, Phi_px);
        PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du,
                        field_lr, Phi_py);

        if constexpr (std::is_same_v<staggered, falseType>) {

            A = A_mid[offset3d];
            PartialZ<local>(k, offset3d, lane_id, A_mid, A, field_du, A_pz);
            PartialX(0, i, k, offset3d, A_mid, A, field_lr, A_px);
            PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du,
                            field_lr, A_py);

        } else {

            S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                                 qtheta_lr, A_mid, field_du, field_lr, A, A_px, A_py, A_pz);

            offset2d = (j + gridGhost) * gridNx + i;

            qtheta = d_qtheta[offset2d];
            qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
            qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
            qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
            qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

            offset2d = j * gridNx + i;
        }

        /*-------------------------PhiA in Parallel Vector Potential--------------------------*/

        for (int index = 0; index < 9; index++)
            compcoes[index] = d_PhiA_A[offset2d * 9 + index];

        dAdt += compcoes[0] * Phi_px * A + compcoes[1] * Phi_px * A_py + compcoes[2] * Phi_px * A_pz +
                compcoes[3] * Phi_py * A + compcoes[4] * Phi_py * A_px + compcoes[5] * Phi_py * A_pz +
                compcoes[6] * Phi_pz * A + compcoes[7] * Phi_pz * A_px + compcoes[8] * Phi_pz * A_py;

        /*--------------------------dNe in Parallel Vector Potential--------------------------*/

        if constexpr (std::is_same_v<Eparallel, trueType>) {

            compcoes[0] = d_APhidNe2A[offset2d * 5 + 4];

            field = dNe_mid[offset3d];
            PartialZ<local>(k, offset3d, lane_id, dNe_mid, field, field_du, field_pz);
            PartialX(0, i, k, offset3d, dNe_mid, field, field_lr, field_px);
            PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid, field_du,
                            field_lr, field_py);

            Ne = Ne0 + field;

            dAdt += compcoes[0] * (Te / Te0 * Ne0 / Ne) * field_py;
        }

        /*--------------------------NeA in Parallel Vector Potential--------------------------*/

        if constexpr (std::is_same_v<Eparallel, trueType>) {

            for (int index = 0; index < 9; index++)
                compcoes[index] = d_NeA_A[offset2d * 9 + index];

            field_px += Ne0_px;

            dAdt += (compcoes[0] * field_px * A + compcoes[1] * field_px * A_py + compcoes[2] * field_px * A_pz +
                     compcoes[3] * field_py * A + compcoes[4] * field_py * A_px + compcoes[5] * field_py * A_pz +
                     compcoes[6] * field_pz * A + compcoes[7] * field_pz * A_px + compcoes[8] * field_pz * A_py) *
                    (Te / Te0 * Ne0 / Ne);
        }
    }

    d_A_pt[offset3d] = dAdt;
}

template <typename local, typename mhdReal, typename picReal>
__global__ void MHD2PIC(picReal* __restrict__ pic3d, mhdReal* __restrict__ globalA, mhdReal* __restrict__ globalPhi,
                        mhdReal* __restrict__ globalApt) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;

    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    /*----------------------------------Field and Derivative----------------------------------*/

    mhdReal field, field_px, field_py, field_pz;
    mhdReal field_du[4];

    /*----------------------------------------MHD2PIC-----------------------------------------*/

    field = globalA[offset3d];

    PartialZ<local>(k, offset3d, offset3d, globalA, field, field_du, field_pz);
    PartialX(0, i, k, offset3d, globalA, field, field_du, field_px);

    field_du[0] = globalA[offset3d - 2 * gridNxz];
    field_du[1] = globalA[offset3d - 1 * gridNxz];
    field_du[2] = globalA[offset3d + 1 * gridNxz];
    field_du[3] = globalA[offset3d + 2 * gridNxz];
    field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * mhdGridDy);

    offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

    pic3d[offset3d + 0] = static_cast<picReal>(field);
    pic3d[offset3d + 1] = static_cast<picReal>(field_px);
    pic3d[offset3d + 2] = static_cast<picReal>(field_py);
    pic3d[offset3d + 3] = static_cast<picReal>(field_pz);

    if (k < gridGhost) {

        offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

        pic3d[offset3d + 0] = static_cast<picReal>(field);
        pic3d[offset3d + 1] = static_cast<picReal>(field_px);
        pic3d[offset3d + 2] = static_cast<picReal>(field_py);
        pic3d[offset3d + 3] = static_cast<picReal>(field_pz);

    } else if (k > gridNz - gridGhost - 1) {

        offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

        pic3d[offset3d + 0] = static_cast<picReal>(field);
        pic3d[offset3d + 1] = static_cast<picReal>(field_px);
        pic3d[offset3d + 2] = static_cast<picReal>(field_py);
        pic3d[offset3d + 3] = static_cast<picReal>(field_pz);
    }

    offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    field = globalPhi[offset3d];

    PartialZ<local>(k, offset3d, offset3d, globalPhi, field, field_du, field_pz);
    PartialX(0, i, k, offset3d, globalPhi, field, field_du, field_px);

    field_du[0] = globalPhi[offset3d - 2 * gridNxz];
    field_du[1] = globalPhi[offset3d - 1 * gridNxz];
    field_du[2] = globalPhi[offset3d + 1 * gridNxz];
    field_du[3] = globalPhi[offset3d + 2 * gridNxz];
    field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * mhdGridDy);

    offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

    pic3d[offset3d + 4] = static_cast<picReal>(field_px);
    pic3d[offset3d + 5] = static_cast<picReal>(field_py);
    pic3d[offset3d + 6] = static_cast<picReal>(field_pz);

    if (k < gridGhost) {

        offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

        pic3d[offset3d + 4] = static_cast<picReal>(field_px);
        pic3d[offset3d + 5] = static_cast<picReal>(field_py);
        pic3d[offset3d + 6] = static_cast<picReal>(field_pz);

    } else if (k > gridNz - gridGhost - 1) {

        offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

        pic3d[offset3d + 4] = static_cast<picReal>(field_px);
        pic3d[offset3d + 5] = static_cast<picReal>(field_py);
        pic3d[offset3d + 6] = static_cast<picReal>(field_pz);
    }

    offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    field = globalApt[offset3d];

    offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

    pic3d[offset3d + 7] = static_cast<picReal>(field);

    if (k < gridGhost) {

        offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

        pic3d[offset3d + 7] = static_cast<picReal>(field);

    } else if (k > gridNz - gridGhost - 1) {

        offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

        pic3d[offset3d + 7] = static_cast<picReal>(field);
    }

    if (j < gridGhost) {

        offset3d = j * gridNxz + i * gridNz + k;

        field = globalA[offset3d];

        PartialZ<local>(k, offset3d, offset3d, globalA, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, globalA, field, field_du, field_px);

        if (j == 0) {
            field_du[0] = globalA[offset3d + 1 * gridNxz];
            field_du[1] = globalA[offset3d + 2 * gridNxz];
            field_du[2] = globalA[offset3d + 3 * gridNxz];
            field_du[3] = globalA[offset3d + 4 * gridNxz];
            field_py = (-25 * field + 48 * field_du[0] - 36 * field_du[1] + 16 * field_du[2] - 3 * field_du[3]) /
                       (12 * mhdGridDy);
        } else if (j == 1) {
            field_du[0] = globalA[offset3d - 1 * gridNxz];
            field_du[1] = globalA[offset3d + 1 * gridNxz];
            field_du[2] = globalA[offset3d + 2 * gridNxz];
            field_du[3] = globalA[offset3d + 3 * gridNxz];
            field_py =
                (-3 * field_du[0] - 10 * field + 18 * field_du[1] - 6 * field_du[2] + field_du[3]) / (12 * mhdGridDy);
        } else {
            field_du[0] = globalA[offset3d - 2 * gridNxz];
            field_du[1] = globalA[offset3d - 1 * gridNxz];
            field_du[2] = globalA[offset3d + 1 * gridNxz];
            field_du[3] = globalA[offset3d + 2 * gridNxz];
            field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * mhdGridDy);
        }

        offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

        pic3d[offset3d + 0] = static_cast<picReal>(field);
        pic3d[offset3d + 1] = static_cast<picReal>(field_px);
        pic3d[offset3d + 2] = static_cast<picReal>(field_py);
        pic3d[offset3d + 3] = static_cast<picReal>(field_pz);

        if (k < gridGhost) {

            offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

            pic3d[offset3d + 0] = static_cast<picReal>(field);
            pic3d[offset3d + 1] = static_cast<picReal>(field_px);
            pic3d[offset3d + 2] = static_cast<picReal>(field_py);
            pic3d[offset3d + 3] = static_cast<picReal>(field_pz);

        } else if (k > gridNz - gridGhost - 1) {

            offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

            pic3d[offset3d + 0] = static_cast<picReal>(field);
            pic3d[offset3d + 1] = static_cast<picReal>(field_px);
            pic3d[offset3d + 2] = static_cast<picReal>(field_py);
            pic3d[offset3d + 3] = static_cast<picReal>(field_pz);
        }

        offset3d = j * gridNxz + i * gridNz + k;

        field = globalPhi[offset3d];

        PartialZ<local>(k, offset3d, offset3d, globalPhi, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, globalPhi, field, field_du, field_px);

        if (j == 0) {
            field_du[0] = globalPhi[offset3d + 1 * gridNxz];
            field_du[1] = globalPhi[offset3d + 2 * gridNxz];
            field_du[2] = globalPhi[offset3d + 3 * gridNxz];
            field_du[3] = globalPhi[offset3d + 4 * gridNxz];
            field_py = (-25 * field + 48 * field_du[0] - 36 * field_du[1] + 16 * field_du[2] - 3 * field_du[3]) /
                       (12 * mhdGridDy);
        } else if (j == 1) {
            field_du[0] = globalPhi[offset3d - 1 * gridNxz];
            field_du[1] = globalPhi[offset3d + 1 * gridNxz];
            field_du[2] = globalPhi[offset3d + 2 * gridNxz];
            field_du[3] = globalPhi[offset3d + 3 * gridNxz];
            field_py =
                (-3 * field_du[0] - 10 * field + 18 * field_du[1] - 6 * field_du[2] + field_du[3]) / (12 * mhdGridDy);
        } else {
            field_du[0] = globalPhi[offset3d - 2 * gridNxz];
            field_du[1] = globalPhi[offset3d - 1 * gridNxz];
            field_du[2] = globalPhi[offset3d + 1 * gridNxz];
            field_du[3] = globalPhi[offset3d + 2 * gridNxz];
            field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * mhdGridDy);
        }

        offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

        pic3d[offset3d + 4] = static_cast<picReal>(field_px);
        pic3d[offset3d + 5] = static_cast<picReal>(field_py);
        pic3d[offset3d + 6] = static_cast<picReal>(field_pz);

        if (k < gridGhost) {

            offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

            pic3d[offset3d + 4] = static_cast<picReal>(field_px);
            pic3d[offset3d + 5] = static_cast<picReal>(field_py);
            pic3d[offset3d + 6] = static_cast<picReal>(field_pz);

        } else if (k > gridNz - gridGhost - 1) {

            offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

            pic3d[offset3d + 4] = static_cast<picReal>(field_px);
            pic3d[offset3d + 5] = static_cast<picReal>(field_py);
            pic3d[offset3d + 6] = static_cast<picReal>(field_pz);
        }

        offset3d = j * gridNxz + i * gridNz + k;

        field = globalApt[offset3d];

        offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

        pic3d[offset3d + 7] = static_cast<picReal>(field);

        if (k < gridGhost) {

            offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

            pic3d[offset3d + 7] = static_cast<picReal>(field);

        } else if (k > gridNz - gridGhost - 1) {

            offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

            pic3d[offset3d + 7] = static_cast<picReal>(field);
        }

    } else if (j > gridNy - gridGhost - 1) {

        offset3d = (j + 2 * gridGhost) * gridNxz + i * gridNz + k;

        field = globalA[offset3d];

        PartialZ<local>(k, offset3d, offset3d, globalA, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, globalA, field, field_du, field_px);

        if (j == gridNy - 1) {
            field_du[0] = globalA[offset3d - 4 * gridNxz];
            field_du[1] = globalA[offset3d - 3 * gridNxz];
            field_du[2] = globalA[offset3d - 2 * gridNxz];
            field_du[3] = globalA[offset3d - 1 * gridNxz];
            field_py = (3 * field_du[0] - 16 * field_du[1] + 36 * field_du[2] - 48 * field_du[3] + 25 * field) /
                       (12 * mhdGridDy);
        } else if (j == gridNy - 2) {
            field_du[0] = globalA[offset3d - 3 * gridNxz];
            field_du[1] = globalA[offset3d - 2 * gridNxz];
            field_du[2] = globalA[offset3d - 1 * gridNxz];
            field_du[3] = globalA[offset3d + 1 * gridNxz];
            field_py =
                (-field_du[0] + 6 * field_du[1] - 18 * field_du[2] + 10 * field + 3 * field_du[3]) / (12 * mhdGridDy);
        } else {
            field_du[0] = globalA[offset3d - 2 * gridNxz];
            field_du[1] = globalA[offset3d - 1 * gridNxz];
            field_du[2] = globalA[offset3d + 1 * gridNxz];
            field_du[3] = globalA[offset3d + 2 * gridNxz];
            field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * mhdGridDy);
        }

        offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

        pic3d[offset3d + 0] = static_cast<picReal>(field);
        pic3d[offset3d + 1] = static_cast<picReal>(field_px);
        pic3d[offset3d + 2] = static_cast<picReal>(field_py);
        pic3d[offset3d + 3] = static_cast<picReal>(field_pz);

        if (k < gridGhost) {

            offset3d =
                ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

            pic3d[offset3d + 0] = static_cast<picReal>(field);
            pic3d[offset3d + 1] = static_cast<picReal>(field_px);
            pic3d[offset3d + 2] = static_cast<picReal>(field_py);
            pic3d[offset3d + 3] = static_cast<picReal>(field_pz);

        } else if (k > gridNz - gridGhost - 1) {

            offset3d =
                ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

            pic3d[offset3d + 0] = static_cast<picReal>(field);
            pic3d[offset3d + 1] = static_cast<picReal>(field_px);
            pic3d[offset3d + 2] = static_cast<picReal>(field_py);
            pic3d[offset3d + 3] = static_cast<picReal>(field_pz);
        }

        offset3d = (j + 2 * gridGhost) * gridNxz + i * gridNz + k;

        field = globalPhi[offset3d];

        PartialZ<local>(k, offset3d, offset3d, globalPhi, field, field_du, field_pz);
        PartialX(0, i, k, offset3d, globalPhi, field, field_du, field_px);

        if (j == gridNy - 1) {
            field_du[0] = globalPhi[offset3d - 4 * gridNxz];
            field_du[1] = globalPhi[offset3d - 3 * gridNxz];
            field_du[2] = globalPhi[offset3d - 2 * gridNxz];
            field_du[3] = globalPhi[offset3d - 1 * gridNxz];
            field_py = (3 * field_du[0] - 16 * field_du[1] + 36 * field_du[2] - 48 * field_du[3] + 25 * field) /
                       (12 * mhdGridDy);
        } else if (j == gridNy - 2) {
            field_du[0] = globalPhi[offset3d - 3 * gridNxz];
            field_du[1] = globalPhi[offset3d - 2 * gridNxz];
            field_du[2] = globalPhi[offset3d - 1 * gridNxz];
            field_du[3] = globalPhi[offset3d + 1 * gridNxz];
            field_py =
                (-field_du[0] + 6 * field_du[1] - 18 * field_du[2] + 10 * field + 3 * field_du[3]) / (12 * mhdGridDy);
        } else {
            field_du[0] = globalPhi[offset3d - 2 * gridNxz];
            field_du[1] = globalPhi[offset3d - 1 * gridNxz];
            field_du[2] = globalPhi[offset3d + 1 * gridNxz];
            field_du[3] = globalPhi[offset3d + 2 * gridNxz];
            field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * mhdGridDy);
        }

        offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

        pic3d[offset3d + 4] = static_cast<picReal>(field_px);
        pic3d[offset3d + 5] = static_cast<picReal>(field_py);
        pic3d[offset3d + 6] = static_cast<picReal>(field_pz);

        if (k < gridGhost) {

            offset3d =
                ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

            pic3d[offset3d + 4] = static_cast<picReal>(field_px);
            pic3d[offset3d + 5] = static_cast<picReal>(field_py);
            pic3d[offset3d + 6] = static_cast<picReal>(field_pz);

        } else if (k > gridNz - gridGhost - 1) {

            offset3d =
                ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

            pic3d[offset3d + 4] = static_cast<picReal>(field_px);
            pic3d[offset3d + 5] = static_cast<picReal>(field_py);
            pic3d[offset3d + 6] = static_cast<picReal>(field_pz);
        }

        offset3d = (j + 2 * gridGhost) * gridNxz + i * gridNz + k;

        field = globalApt[offset3d];

        offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

        pic3d[offset3d + 7] = static_cast<picReal>(field);

        if (k < gridGhost) {

            offset3d =
                ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

            pic3d[offset3d + 7] = static_cast<picReal>(field);

        } else if (k > gridNz - gridGhost - 1) {

            offset3d =
                ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

            pic3d[offset3d + 7] = static_cast<picReal>(field);
        }
    }
}

template <typename cufftType, typename type>
__global__ void MHDDiagAmplitude(cufftType* __restrict__ d_freq, type* __restrict__ d_amplitude,
                                 type* __restrict__ d_modeReal, type* __restrict__ d_modeImag) {

    /*-------------------------------------Related Index--------------------------------------*/

    int offset = (diagY * gridNx + threadIdx.x) * nFFTFreqSize;
    type real, imag, factor;

    /*---------------------------------------Amplitude----------------------------------------*/

    for (int mode = leftN; mode <= rightN; mode++) {

        real = d_freq[offset + mode].x;
        imag = d_freq[offset + mode].y;
        factor = (mode == 0 || mode == gridNz / 2) ? 1 : 2;

        if constexpr (std::is_same_v<type, double>)
            d_amplitude[threadIdx.x * (rightN - leftN + 1) + mode - leftN] =
                sqrt(real * real + imag * imag) / gridNz * factor;
        else
            d_amplitude[threadIdx.x * (rightN - leftN + 1) + mode - leftN] =
                sqrtf(real * real + imag * imag) / gridNz * factor;

        d_modeReal[threadIdx.x * (rightN - leftN + 1) + mode - leftN] = real;
        d_modeImag[threadIdx.x * (rightN - leftN + 1) + mode - leftN] = imag;
    }
}

template <typename type>
__global__ void MHDDiagFrequency(type* __restrict__ Phi_mid, type* __restrict__ d_frequency) {

    /*-------------------------------------Related Index--------------------------------------*/

    int offset3d = (diagY + gridGhost) * gridNxz + threadIdx.x * gridNz;

    /*---------------------------------------Frequency----------------------------------------*/

    d_frequency[threadIdx.x] = Phi_mid[offset3d];
}

template <typename type>
__global__ void MHDDiagZFDrive(type* __restrict__ w_beg, type* __restrict__ w_end, type* __restrict__ Maxwell,
                               type* __restrict__ Reynolds, type* __restrict__ MaxwellDrive,
                               type* __restrict__ ReynoldsDrive, type* __restrict__ dwdtTotal) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = threadIdx.x;
    int offset3d = (diagY + gridGhost) * gridNxz + i * gridNz;

    /*----------------------------------------Diagnose----------------------------------------*/

    MaxwellDrive[i] = Maxwell[offset3d] / mhdGridDt;
    ReynoldsDrive[i] = Reynolds[offset3d] / mhdGridDt;
    dwdtTotal[i] = (w_end[offset3d] - w_beg[offset3d]) / mhdGridDt;
}

template <typename nonlinear, typename staggered, typename Eparallel, typename type>
__global__ void MHDDiagEparallel(type* __restrict__ d_qtheta, type* __restrict__ A_mid, type* __restrict__ dNe_mid,
                                 type* __restrict__ dTe_mid, type* __restrict__ Phi_mid, type* __restrict__ d_Ne0,
                                 type* __restrict__ d_Te0, type* __restrict__ d_Ne0_px, type* __restrict__ d_APhidNe2A,
                                 type* __restrict__ d_PhiA_A, type* __restrict__ d_NeA_A, type* __restrict__ d_Epara,
                                 type* __restrict__ d_EparaES) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = diagY;
    int k = 0;
    int offset2d = (j + gridGhost) * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    type qtheta;
    type qtheta_lr[4];
    int shift_k;
    type shift_lk;
    type shift_dk;

    /*----------------------------------Field and Derivative----------------------------------*/

    type field;
    type field_px, field_py, field_pz;
    type field_du[4];
    type field_lr[4];

    /*---------------------------------Compressed Coefficient---------------------------------*/

    type compcoes[9];

    /*----------------------------------------Diagnose----------------------------------------*/

    type Epara, EparaES;

    /*---------------------------------------Initialize---------------------------------------*/

    qtheta = d_qtheta[offset2d];
    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    offset2d = j * gridNx + i;

    Epara = 0;
    EparaES = 0;

    /*-----------------Parallel Vector Potential in Parallel Vector Potential-----------------*/

    for (int index = 0; index < 5; index++)
        compcoes[index] = d_APhidNe2A[offset2d * 5 + index];

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        if constexpr (std::is_same_v<nonlinear, falseType>) {

            field = A_mid[offset3d];
            PartialZ<falseType>(k, offset3d, lane_id, A_mid, field, field_du, field_pz);
            PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du,
                                field_lr, field_py);
            Epara += compcoes[0] * field + compcoes[1] * field_py + compcoes[2] * field_pz;
        }
    }

    /*--------------------Electric Potential in Parallel Vector Potential---------------------*/

    if constexpr (std::is_same_v<staggered, trueType>) {

        Collocated2S<falseType>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta,
                                qtheta, qtheta_lr, Phi_mid, field_du, field_lr, field);
        field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * mhdGridDy);

    } else {

        PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du,
                            field_lr, field_py);
    }

    EparaES += compcoes[3] * field_py;

    /*---------------------Electron Density in Parallel Vector Potential----------------------*/

    if constexpr (std::is_same_v<Eparallel, trueType>) {

        if constexpr (std::is_same_v<nonlinear, falseType>) {

            if constexpr (std::is_same_v<staggered, trueType>) {

                Collocated2S<falseType>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
                                        d_qtheta, qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field);
                field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * mhdGridDy);

            } else {

                PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid,
                                    field_du, field_lr, field_py);
            }

            Epara += compcoes[4] * field_py;
        }
    }

    if constexpr (std::is_same_v<nonlinear, trueType>) {

        type Ne0, Te0, Ne0_px, Ne, Te;
        type A, A_px, A_py, A_pz;
        type Phi, Phi_px, Phi_py, Phi_pz;

        Ne0 = d_Ne0[offset2d];
        Te0 = d_Te0[offset2d];
        Ne0_px = d_Ne0_px[offset2d];
        Te = Te0 + dTe_mid[offset3d];

        offset2d = (j + gridGhost) * gridNx + i;

        qtheta = d_qtheta[offset2d];
        qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
        qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
        qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
        qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

        offset2d = j * gridNx + i;

        A = A_mid[offset3d];
        PartialZ<falseType>(k, offset3d, lane_id, A_mid, A, field_du, A_pz);
        PartialX(0, i, k, offset3d, A_mid, A, field_lr, A_px);
        PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du,
                            field_lr, A_py);

        if constexpr (std::is_same_v<staggered, falseType>) {

            Phi = Phi_mid[offset3d];
            PartialZ<falseType>(k, offset3d, lane_id, Phi_mid, Phi, field_du, Phi_pz);
            PartialX(0, i, k, offset3d, Phi_mid, Phi, field_lr, Phi_px);
            PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du,
                                field_lr, Phi_py);

        } else {

            C2SPartialXYZ<falseType>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta,
                                     qtheta, qtheta_lr, Phi_mid, field_du, field_lr, Phi, Phi_px, Phi_py, Phi_pz);
        }

        /*-------------------------PhiA in Parallel Vector Potential--------------------------*/

        for (int index = 0; index < 9; index++)
            compcoes[index] = d_PhiA_A[offset2d * 9 + index];

        EparaES += compcoes[0] * Phi_px * A + compcoes[1] * Phi_px * A_py + compcoes[2] * Phi_px * A_pz +
                   compcoes[3] * Phi_py * A + compcoes[4] * Phi_py * A_px + compcoes[5] * Phi_py * A_pz +
                   compcoes[6] * Phi_pz * A + compcoes[7] * Phi_pz * A_px + compcoes[8] * Phi_pz * A_py;

        /*--------------------------dNe in Parallel Vector Potential--------------------------*/

        if constexpr (std::is_same_v<Eparallel, trueType>) {

            compcoes[0] = d_APhidNe2A[offset2d * 5 + 4];

            if constexpr (std::is_same_v<staggered, falseType>) {

                field = dNe_mid[offset3d];
                PartialZ<falseType>(k, offset3d, lane_id, dNe_mid, field, field_du, field_pz);
                PartialX(0, i, k, offset3d, dNe_mid, field, field_lr, field_px);
                PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid,
                                    field_du, field_lr, field_py);

                Ne = Ne0 + field;

            } else {

                C2SPartialXYZ<falseType>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta,
                                         qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field, field_px, field_py,
                                         field_pz);

                Collocated2S<falseType>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
                                        d_qtheta, qtheta, qtheta_lr, dTe_mid, field_du, field_lr, Te);

                Ne = Ne0 + field;
                Te = Te0 + Te;
            }

            Epara += compcoes[0] * (Te / Te0 * Ne0 / Ne) * field_py;
        }

        /*--------------------------NeA in Parallel Vector Potential--------------------------*/

        if constexpr (std::is_same_v<Eparallel, trueType>) {

            for (int index = 0; index < 9; index++)
                compcoes[index] = d_NeA_A[offset2d * 9 + index];

            field_px += Ne0_px;

            Epara += (compcoes[0] * field_px * A + compcoes[1] * field_px * A_py + compcoes[2] * field_px * A_pz +
                      compcoes[3] * field_py * A + compcoes[4] * field_py * A_px + compcoes[5] * field_py * A_pz +
                      compcoes[6] * field_pz * A + compcoes[7] * field_pz * A_px + compcoes[8] * field_pz * A_py) *
                     (Te / Te0 * Ne0 / Ne);
        }
    }

    d_Epara[blockIdx.x * blockDim.x + threadIdx.x] = -Epara;
    d_EparaES[blockIdx.x * blockDim.x + threadIdx.x] = EparaES;
}

template <typename type, typename... types>
__global__ void MHDCheckNAN(int* __restrict__ flag, type* __restrict__ first, types* __restrict__... second) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    /*------------------------------------------NaN-------------------------------------------*/

    bool bad = isnan(first[offset3d]);
    if constexpr (sizeof...(second) > 0)
        ((bad |= isnan(second[offset3d])), ...);

    if (bad)
        atomicOr(flag, 1);
}
