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

/*----------------------------------------PIC Kernels-----------------------------------------*/

template <typename local, typename type>
__global__ void PICAlignedGhost(type* __restrict__ d_qtheta, type* __restrict__ dP_mid) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = j * gridNx + i;
    int offset3d = j * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    int shift_k;
    type shift_lk;
    type shift_dk;
    type qtheta;
    type field;
    type field_du[4];

    qtheta = (d_qtheta[offset2d + gridNx] - d_qtheta[offset2d]) * gridNy;

    /*---------------------------------------Left Ghost---------------------------------------*/

    shift_lk = qtheta / mhdGridDz;
    if constexpr (std::is_same_v<type, double>)
        shift_k = __double2int_rd(shift_lk);
    else
        shift_k = __float2int_rd(shift_lk);
    shift_dk = shift_lk - shift_k;

    field = dP_mid[offset3d];

    if constexpr (std::is_same_v<local, trueType>) {

        field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
        field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
        field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
        field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

    } else {

        field_du[0] = dP_mid[offset3d - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
        field_du[1] = dP_mid[offset3d - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
        field_du[2] = dP_mid[offset3d - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
        field_du[3] = dP_mid[offset3d - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];
    }

    field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
            (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
            shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
            shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

    dP_mid[offset3d + gridNy * gridNxz] += field;

    /*--------------------------------------Right Ghost---------------------------------------*/

    shift_lk = -qtheta / mhdGridDz;
    if constexpr (std::is_same_v<type, double>)
        shift_k = __double2int_rd(shift_lk);
    else
        shift_k = __float2int_rd(shift_lk);
    shift_dk = shift_lk - shift_k;

    field = dP_mid[offset3d + (gridNy + gridGhost) * gridNxz];

    if constexpr (std::is_same_v<local, trueType>) {

        field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
        field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
        field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
        field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

    } else {

        field_du[0] =
            dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
        field_du[1] =
            dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
        field_du[2] =
            dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
        field_du[3] =
            dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];
    }

    field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
            (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
            shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
            shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

    dP_mid[offset3d + gridGhost * gridNxz] += field;
}

template <int ratioDt, picType particle, disType distribution, typename nonlinear, typename QNeutrality,
          typename mhdReal, typename picReal>
__global__ void DriftAlignedRK4(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d, picReal* __restrict__ pic3d,
                                int* __restrict__ pic_keys_in, picReal* __restrict__ pic_values_in,
                                mhdReal* __restrict__ dP_mid, mhdReal* __restrict__ dN_mid) {

    int illegal;
    int i, j, k;
    int qId, tileId, cellId, picId;

    picReal flag;
    picReal li, lj, lk;
    picReal coes[8] = {};

    picReal dx, dy, dz, disP, mu;
    picReal ddt[5] = {};
    picReal vec0[5] = {};
    picReal vec1[5] = {};
    picReal vec2[5] = {};

    picReal q, q_px, J, J_px, J_py, B, B_px, B_py;
    picReal gcovxy, gcovyy, gcovyz;
    picReal gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
    picReal APhiApt[8] = {};

    picReal bx, by, bz;
    picReal rho, bcony;
    picReal cx, cy, cz;
    picReal m2e, mu2e;
    picReal dxy, dxz, dyz;
    picReal dxB, dyB, dzB;
    picReal Bstarx, Bstary, Bstarz, Bstar;
    picReal invJ, invB, invQ, invRho, invBstar, invM2e, bconyOverJ;
    picReal na, na_px, nb, nb_px, ni, ni_px;
    picReal ta, ta_px, tb, tb_px, ti, ti_px;
    picReal V, E, cdwdt;

    picReal dvdt1, dxdt1, dydt1;
    picReal dxPhi, dyPhi, dzPhi;
    picReal cxdxA, cydyA, czdzA;

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partChar = (particle == Ion) ? IonChar : (particle == Alpha) ? AlphaChar : BeamChar;
    const picReal partVb = (particle == Ion) ? IonVb : (particle == Alpha) ? AlphaVb : BeamVb;
    const picReal partDeltaV = (particle == Ion) ? IonDeltaV : (particle == Alpha) ? AlphaDeltaV : BeamDeltaV;
    const picReal partLambda0 = (particle == Ion) ? IonLambda0 : (particle == Alpha) ? AlphaLambda0 : BeamLambda0;
    const picReal partDeltaLambda2 = (particle == Ion)     ? IonDeltaLambda2
                                     : (particle == Alpha) ? AlphaDeltaLambda2
                                                           : BeamDeltaLambda2;
    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        cellId = pic_keys_in[picId];
        vec0[0] = pic_values_in[picId + 0 * picDev];
        vec0[1] = pic_values_in[picId + 1 * picDev];
        vec0[2] = pic_values_in[picId + 2 * picDev];
        vec0[3] = pic_values_in[picId + 3 * picDev];
        vec0[4] = pic_values_in[picId + 4 * picDev];
        dis = pic_values_in[picId + 5 * picDev];
        mu = pic_values_in[picId + 6 * picDev];

        li = (vec0[0] - xbeg) / picGridDx;
        lj = (vec0[1] - ybeg) / picGridDy;
        lk = (vec0[2] - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dx = li - i;
        dy = lj - j;
        dz = lk - k;

        for (int index = 0; index < 5; index++) {
            vec1[index] = vec0[index];
            vec2[index] = vec0[index];
        }

        qId = i * qStride;
        tileId = (j * cellNx + i) * tileStride;
        cellId *= cellStride;
        illegal = 0;

        auto dfdt_XVpara = [&]() {
            picReal& n = (particle == Ion ? ni : (particle == Alpha ? na : nb));
            picReal& n_px = (particle == Ion ? ni_px : (particle == Alpha ? na_px : nb_px));
            picReal& t = (particle == Ion ? ti : (particle == Alpha ? ta : tb));
            picReal& t_px = (particle == Ion ? ti_px : (particle == Alpha ? ta_px : tb_px));

            if constexpr (distribution == Maxwell) {

                cdwdt = mp * va * va / (t * kev);

                ddt[4] = (-n_px / n + 3 * t_px / (2 * t) + (mu * B_px - t_px / t * E) * cdwdt) * dxdt1;
                ddt[4] += mu * B_py * cdwdt * dydt1;
                ddt[4] += partMass * vec1[3] * cdwdt * dvdt1;

            } else {

                cdwdt = 3 * V / (2 * E * V + partMass * t * t * t);

                ddt[4] = (-n_px / n + mu * B_px * cdwdt) * dxdt1;
                ddt[4] += mu * B_py * cdwdt * dydt1;
                ddt[4] += partMass * vec1[3] * cdwdt * dvdt1;

                ddt[4] += 3 * t * t / (V * V * V + t * t * t) * t_px * dxdt1;

                if constexpr (distribution == Slowing0) {

                    cdwdt = 0;

                } else if constexpr (distribution == Slowing1) {

                    if constexpr (std::is_same_v<picReal, double>)
                        cdwdt = 2 * exp(-pow((partVb - V) / partDeltaV, 2.0)) /
                                (partMass * V * partDeltaV * sqrt(pi) * (1 + erf((partVb - V) / partDeltaV)));
                    else
                        cdwdt = 2 * expf(-powf((partVb - V) / partDeltaV, 2.0f)) /
                                (partMass * V * partDeltaV * sqrtf(pi) * (1 + erff((partVb - V) / partDeltaV)));

                } else if constexpr (distribution == Slowing2) {

                    cdwdt = 2 * mu * (partLambda0 * E - mu) / (partDeltaLambda2 * E * E * E);

                } else if constexpr (distribution == Slowing3) {

                    cdwdt = 2 * mu * (partLambda0 * E - mu) / (partDeltaLambda2 * E * E * E);

                    if constexpr (std::is_same_v<picReal, double>)
                        cdwdt += 2 * exp(-pow((partVb - V) / partDeltaV, 2.0)) /
                                 (partMass * V * partDeltaV * sqrt(pi) * (1 + erf((partVb - V) / partDeltaV)));
                    else
                        cdwdt += 2 * expf(-powf((partVb - V) / partDeltaV, 2.0f)) /
                                 (partMass * V * partDeltaV * sqrtf(pi) * (1 + erff((partVb - V) / partDeltaV)));
                }

                ddt[4] += mu * B_px * cdwdt * dxdt1;
                ddt[4] += mu * B_py * cdwdt * dydt1;
                ddt[4] += partMass * vec1[3] * cdwdt * dvdt1;
            }

            if constexpr (std::is_same_v<nonlinear, trueType>)
                ddt[4] *= (dis - vec1[4]);
            else
                ddt[4] *= dis;
        };

        auto interpRK4 = [&]() {
            for (int index = 0; index < 2; index++)
                coes[index] = (hx[index] + sx[index] * dx);
            FieldGather1d2d<2>(qId, coes, pic1d, q, q_px, na, na_px, nb, nb_px, ni, ni_px, ta, ta_px, tb, tb_px, ti,
                               ti_px);

            for (int index = 0; index < 2; index++)
                coes[index + 2] = coes[index];

            for (int index = 0; index < 4; index++)
                coes[index] *= (hy[index] + sy[index] * dy);
            FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py, gcovxy, gcovyy, gcovyz, gcovxy_py,
                               gcovyy_px, gcovyz_px, gcovyz_py);

            for (int index = 0; index < 4; index++)
                coes[index + 4] = coes[index];

            for (int index = 0; index < 8; index++)
                coes[index] *= (hz[index] + sz[index] * dz);

            for (int index = 0; index < 8; index++)
                APhiApt[index] = 0;

            FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

            m2e = cm * partMass / partChar;
            mu2e = cm * mu / partChar;
            E = partMass * vec1[3] * vec1[3] / 2 + mu * B;
            if constexpr (std::is_same_v<picReal, double>)
                V = sqrt(2.0 * E / partMass);
            else
                V = sqrtf(2.0f * E / partMass);

            rho = rho0 + vec1[0] * drho;
            bcony = 2 * psitmax * drho * rho / (q * J * B);

            invJ = 1 / J;
            invB = 1 / B;
            invQ = 1 / q;
            invRho = 1 / rho;
            invM2e = 1 / m2e;
            bconyOverJ = bcony * invJ;

            bx = bcony * gcovxy;
            by = bcony * gcovyy;
            bz = bcony * gcovyz;

            cx = bconyOverJ * (gcovyz_py - gcovyz * (J_py * invJ + B_py * invB));
            cy = -bconyOverJ *
                 (gcovyz_px - gcovyz * (J_px * invJ + B_px * invB) + gcovyz * (drho * invRho - q_px * invQ));
            cz = bconyOverJ * (gcovyy_px - gcovyy * (J_px * invJ + B_px * invB) - gcovxy_py +
                               gcovxy * (J_py * invJ + B_py * invB) + gcovyy * (drho * invRho - q_px * invQ));

            dxy = bconyOverJ * gcovyz;
            dxz = bconyOverJ * gcovyy;
            dyz = bconyOverJ * gcovxy;

            dxB = dxy * B_py;
            dyB = -dxy * B_px;
            dzB = dxz * B_px - dyz * B_py;

            Bstarx = cx * m2e * vec1[3];
            Bstary = cy * m2e * vec1[3] + B * bcony;
            Bstarz = cz * m2e * vec1[3];

            cxdxA = cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
            cydyA = cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
            czdzA = cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

            dxPhi = dxy * APhiApt[5] - dxz * APhiApt[6];
            dyPhi = dyz * APhiApt[6] - dxy * APhiApt[4];
            dzPhi = dxz * APhiApt[4] - dyz * APhiApt[5];

            if constexpr (std::is_same_v<nonlinear, trueType>) {

                Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
                invBstar = 1 / Bstar;
                dxdt1 = invBstar * (vec1[3] * Bstarx - mu2e * dxB);
                dydt1 = invBstar * (vec1[3] * Bstary - mu2e * dyB);
                dvdt1 = -invM2e * invBstar * mu2e * (Bstarx * B_px + Bstary * B_py);

                Bstarx += cxdxA;
                Bstary += cydyA;
                Bstarz += czdzA;

                Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
                invBstar = 1 / Bstar;

                ddt[0] = invBstar * (vec1[3] * Bstarx - dxPhi - mu2e * dxB);
                ddt[1] = invBstar * (vec1[3] * Bstary - dyPhi - mu2e * dyB);
                ddt[2] = invBstar * (vec1[3] * Bstarz - dzPhi - mu2e * dzB);
                ddt[3] =
                    -invM2e * (APhiApt[7] + invBstar * (Bstarx * (APhiApt[4] + mu2e * B_px) +
                                                        Bstary * (APhiApt[5] + mu2e * B_py) + Bstarz * APhiApt[6]));

                dxdt1 = ddt[0] - dxdt1;
                dydt1 = ddt[1] - dydt1;
                dvdt1 = ddt[3] - dvdt1;
            } else {
                Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
                invBstar = 1 / Bstar;

                ddt[0] = invBstar * (vec1[3] * Bstarx - mu2e * dxB);
                ddt[1] = invBstar * (vec1[3] * Bstary - mu2e * dyB);
                ddt[2] = invBstar * (vec1[3] * Bstarz - mu2e * dzB);
                ddt[3] = -invM2e * invBstar * mu2e * (Bstarx * B_px + Bstary * B_py);

                dxdt1 = invBstar * (vec1[3] * cxdxA - dxPhi);
                dydt1 = invBstar * (vec1[3] * cydyA - dyPhi);
                dvdt1 =
                    -invM2e * (APhiApt[7] + invBstar * (Bstarx * APhiApt[4] + Bstary * APhiApt[5] +
                                                        Bstarz * APhiApt[6] + mu2e * (cxdxA * B_px + cydyA * B_py)));
            }

            dfdt_XVpara();
        };

        auto advanceRK4 = [&](picReal bRK, picReal cRK) {
            for (int index = 0; index < 5; index++)
                vec2[index] += ddt[index] * picGridDt * ratioDt * bRK;

            flag = vec0[0] + ddt[0] * picGridDt * ratioDt * cRK;
            if (flag < xbeg || flag >= xend)
                illegal = 1;
            if (!illegal)
                for (int index = 0; index < 5; index++)
                    vec1[index] = vec0[index] + ddt[index] * picGridDt * ratioDt * cRK;
            else
                for (int index = 0; index < 5; index++)
                    vec1[index] = vec0[index];

            li = (vec1[0] - xbeg) / picGridDx;
            if constexpr (std::is_same_v<picReal, double>)
                i = __double2int_rd(li);
            else
                i = __float2int_rd(li);
            dx = li - i;
            qId = i * qStride;
            coes[0] = hx[0] + sx[0] * dx;
            coes[1] = hx[1] + sx[1] * dx;
            FieldGather1d2d<2>(qId, coes, pic1d, q);

            if constexpr (std::is_same_v<picReal, double>) {
                vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
                vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
                vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
            } else {
                vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
                vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
                vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
            }

            lj = (vec1[1] - ybeg) / picGridDy;
            lk = (vec1[2] - zbeg) / picGridDz;

            if constexpr (std::is_same_v<picReal, double>) {
                j = __double2int_rd(lj);
                k = __double2int_rd(lk);
            } else {
                j = __float2int_rd(lj);
                k = __float2int_rd(lk);
            }

            dy = lj - j;
            dz = lk - k;
            tileId = (j * cellNx + i) * tileStride;
            cellId = (j * cellNxz + i * cellNz + k) * cellStride;
        };

        /*--------------------------------------1st RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 6, picReal(1) / 2);

        /*--------------------------------------2nd RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 3, picReal(1) / 2);

        /*--------------------------------------3rd RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 3, picReal(1));

        /*--------------------------------------4th RK4---------------------------------------*/

        interpRK4();

        for (int index = 0; index < 5; index++)
            vec2[index] += ddt[index] * picGridDt * ratioDt / 6;

        if (vec2[0] < xbeg || vec2[0] >= xend)
            illegal = 1;
        if (illegal) {
            for (int index = 0; index < 5; index++)
                vec2[index] = vec0[index];
            vec2[1] = -vec0[1];
            vec2[4] = 0;
        }

        li = (vec2[0] - xbeg) / picGridDx;
        if constexpr (std::is_same_v<picReal, double>)
            i = __double2int_rd(li);
        else
            i = __float2int_rd(li);
        dx = li - i;
        qId = i * qStride;
        coes[0] = hx[0] + sx[0] * dx;
        coes[1] = hx[1] + sx[1] * dx;
        FieldGather1d2d<2>(qId, coes, pic1d, q);

        if constexpr (std::is_same_v<picReal, double>) {
            vec2[2] = vec2[2] + q * floor((vec2[1] - yori) / yrange) * yrange;
            vec2[1] = vec2[1] - floor((vec2[1] - yori) / yrange) * yrange;
            vec2[2] = vec2[2] - floor((vec2[2] - zori) / zrange) * zrange;
        } else {
            vec2[2] = vec2[2] + q * floorf((vec2[1] - yori) / yrange) * yrange;
            vec2[1] = vec2[1] - floorf((vec2[1] - yori) / yrange) * yrange;
            vec2[2] = vec2[2] - floorf((vec2[2] - zori) / zrange) * zrange;
        }

        lj = (vec2[1] - ybeg) / picGridDy;
        lk = (vec2[2] - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dy = lj - j;
        dz = lk - k;
        tileId = (j * cellNx + i) * tileStride;
        cellId = j * cellNxz + i * cellNz + k;

        for (int index = 0; index < 2; index++)
            coes[index + 2] = coes[index];

        for (int index = 0; index < 4; index++)
            coes[index] *= (hy[index] + sy[index] * dy);

        FieldGather1d2d<4>(tileId, coes, pic2d, J, B);

        disP = (partMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * partConst / 2;

        for (int index = 0; index < 4; index++)
            coes[index + 4] = coes[index];

        for (int index = 0; index < 8; index++)
            coes[index] *= (hz[index] + sz[index] * dz);

        if (i == 0) {

            for (int index = 0; index < 4; index++)
                coes[2 * index] = 0;
        } else if (i == gridNx - 2) {

            for (int index = 0; index < 4; index++)
                coes[2 * index + 1] = 0;
        }

        k = (k - gridGhost + gridNz) % gridNz;
        qId = j * gridNxz + i * gridNz + k;

        i = gridNz;
        j = gridNxz;
        k = (k + 1 + gridNz) % gridNz - k;

        atomicAdd(&dP_mid[qId], static_cast<mhdReal>(coes[0] * disP));
        atomicAdd(&dP_mid[qId + i], static_cast<mhdReal>(coes[1] * disP));
        atomicAdd(&dP_mid[qId + j], static_cast<mhdReal>(coes[2] * disP));
        atomicAdd(&dP_mid[qId + i + j], static_cast<mhdReal>(coes[3] * disP));
        atomicAdd(&dP_mid[qId + k], static_cast<mhdReal>(coes[4] * disP));
        atomicAdd(&dP_mid[qId + i + k], static_cast<mhdReal>(coes[5] * disP));
        atomicAdd(&dP_mid[qId + j + k], static_cast<mhdReal>(coes[6] * disP));
        atomicAdd(&dP_mid[qId + i + j + k], static_cast<mhdReal>(coes[7] * disP));

        if constexpr (std::is_same_v<QNeutrality, trueType>) {
            picReal disN = vec2[4] / J * partConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) * l0 * l0 * l0;

            atomicAdd(&dN_mid[qId], static_cast<mhdReal>(coes[0] * disN));
            atomicAdd(&dN_mid[qId + i], static_cast<mhdReal>(coes[1] * disN));
            atomicAdd(&dN_mid[qId + j], static_cast<mhdReal>(coes[2] * disN));
            atomicAdd(&dN_mid[qId + i + j], static_cast<mhdReal>(coes[3] * disN));
            atomicAdd(&dN_mid[qId + k], static_cast<mhdReal>(coes[4] * disN));
            atomicAdd(&dN_mid[qId + i + k], static_cast<mhdReal>(coes[5] * disN));
            atomicAdd(&dN_mid[qId + j + k], static_cast<mhdReal>(coes[6] * disN));
            atomicAdd(&dN_mid[qId + i + j + k], static_cast<mhdReal>(coes[7] * disN));
        }

        for (int index = 0; index < 7; index++)
            pic_keys_in[picId + index * picDev] = cellId;
        pic_values_in[picId + 0 * picDev] = vec2[0];
        pic_values_in[picId + 1 * picDev] = vec2[1];
        pic_values_in[picId + 2 * picDev] = vec2[2];
        pic_values_in[picId + 3 * picDev] = vec2[3];
        pic_values_in[picId + 4 * picDev] = vec2[4];
    }
}

template <int ratioDt, int gyroNums, picType particle, disType distribution, typename nonlinear, typename QNeutrality,
          typename mhdReal, typename picReal>
__global__ void GyroAlignedRK4(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d, picReal* __restrict__ pic3d,
                               int* __restrict__ pic_keys_in, picReal* __restrict__ pic_values_in,
                               mhdReal* __restrict__ dP_mid, mhdReal* __restrict__ dN_mid) {

    int illegal;
    int i, j, k;
    int qId, tileId, cellId, picId;

    picReal flag;
    picReal li, lj, lk;
    picReal coes[8] = {};

    picReal dx, dy, dz, disP, mu;
    picReal ddt[5] = {};
    picReal vec0[5] = {};
    picReal vec1[5] = {};
    picReal vec2[5] = {};

    picReal q, q_px, J, J_px, J_py, B, B_px, B_py;
    picReal gcovxy, gcovyy, gcovyz;
    picReal gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
    picReal APhiApt[8] = {};

    picReal bx, by, bz;
    picReal rho, bcony;
    picReal cx, cy, cz;
    picReal m2e, mu2e;
    picReal dxy, dxz, dyz;
    picReal dxB, dyB, dzB;
    picReal Bstarx, Bstary, Bstarz, Bstar;
    picReal invJ, invB, invQ, invRho, invBstar, invM2e, bconyOverJ;
    picReal na, na_px, nb, nb_px, ni, ni_px;
    picReal ta, ta_px, tb, tb_px, ti, ti_px;
    picReal V, E, cdwdt;

    picReal gyroDx, gyroDy;
    picReal gyroX, gyroY, gyroZ;
    picReal gconxx, gconxy, gconyy;
    picReal R0, Z0, R1, Z1, angle, radius;
    picReal halfAngle, sinAngle, sinHalfAngle, cosHalfAngle, invSinAngle, gyroTheta, sinA, cosA;
    picReal avecxdxA, avecydyA, aveczdzA;
    picReal avedxPhi, avedyPhi, avedzPhi;
    picReal avePhipx, avePhipy, avePhipz;
    picReal aveAptbx, aveAptby, aveAptbz;
    picReal dvdt1, dxdt1, dydt1;

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partChar = (particle == Ion) ? IonChar : (particle == Alpha) ? AlphaChar : BeamChar;
    const picReal partVb = (particle == Ion) ? IonVb : (particle == Alpha) ? AlphaVb : BeamVb;
    const picReal partDeltaV = (particle == Ion) ? IonDeltaV : (particle == Alpha) ? AlphaDeltaV : BeamDeltaV;
    const picReal partLambda0 = (particle == Ion) ? IonLambda0 : (particle == Alpha) ? AlphaLambda0 : BeamLambda0;
    const picReal partDeltaLambda2 = (particle == Ion)     ? IonDeltaLambda2
                                     : (particle == Alpha) ? AlphaDeltaLambda2
                                                           : BeamDeltaLambda2;
    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        cellId = pic_keys_in[picId];
        vec0[0] = pic_values_in[picId + 0 * picDev];
        vec0[1] = pic_values_in[picId + 1 * picDev];
        vec0[2] = pic_values_in[picId + 2 * picDev];
        vec0[3] = pic_values_in[picId + 3 * picDev];
        vec0[4] = pic_values_in[picId + 4 * picDev];
        dis = pic_values_in[picId + 5 * picDev];
        mu = pic_values_in[picId + 6 * picDev];

        li = (vec0[0] - xbeg) / picGridDx;
        lj = (vec0[1] - ybeg) / picGridDy;
        lk = (vec0[2] - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dx = li - i;
        dy = lj - j;
        dz = lk - k;

        for (int index = 0; index < 5; index++) {
            vec1[index] = vec0[index];
            vec2[index] = vec0[index];
        }

        qId = i * qStride;
        tileId = (j * cellNx + i) * tileStride;
        cellId *= cellStride;
        illegal = 0;

        auto dfdt_XVpara = [&]() {
            picReal& n = (particle == Ion ? ni : (particle == Alpha ? na : nb));
            picReal& n_px = (particle == Ion ? ni_px : (particle == Alpha ? na_px : nb_px));
            picReal& t = (particle == Ion ? ti : (particle == Alpha ? ta : tb));
            picReal& t_px = (particle == Ion ? ti_px : (particle == Alpha ? ta_px : tb_px));

            if constexpr (distribution == Maxwell) {

                cdwdt = mp * va * va / (t * kev);

                ddt[4] = (-n_px / n + 3 * t_px / (2 * t) + (mu * B_px - t_px / t * E) * cdwdt) * dxdt1;
                ddt[4] += mu * B_py * cdwdt * dydt1;
                ddt[4] += partMass * vec1[3] * cdwdt * dvdt1;

            } else {

                cdwdt = 3 * V / (2 * E * V + partMass * t * t * t);

                ddt[4] = (-n_px / n + mu * B_px * cdwdt) * dxdt1;
                ddt[4] += mu * B_py * cdwdt * dydt1;
                ddt[4] += partMass * vec1[3] * cdwdt * dvdt1;

                ddt[4] += 3 * t * t / (V * V * V + t * t * t) * t_px * dxdt1;

                if constexpr (distribution == Slowing0) {

                    cdwdt = 0;

                } else if constexpr (distribution == Slowing1) {

                    if constexpr (std::is_same_v<picReal, double>)
                        cdwdt = 2 * exp(-pow((partVb - V) / partDeltaV, 2.0)) /
                                (partMass * V * partDeltaV * sqrt(pi) * (1 + erf((partVb - V) / partDeltaV)));
                    else
                        cdwdt = 2 * expf(-powf((partVb - V) / partDeltaV, 2.0f)) /
                                (partMass * V * partDeltaV * sqrtf(pi) * (1 + erff((partVb - V) / partDeltaV)));

                } else if constexpr (distribution == Slowing2) {

                    cdwdt = 2 * mu * (partLambda0 * E - mu) / (partDeltaLambda2 * E * E * E);

                } else if constexpr (distribution == Slowing3) {

                    cdwdt = 2 * mu * (partLambda0 * E - mu) / (partDeltaLambda2 * E * E * E);

                    if constexpr (std::is_same_v<picReal, double>)
                        cdwdt += 2 * exp(-pow((partVb - V) / partDeltaV, 2.0)) /
                                 (partMass * V * partDeltaV * sqrt(pi) * (1 + erf((partVb - V) / partDeltaV)));
                    else
                        cdwdt += 2 * expf(-powf((partVb - V) / partDeltaV, 2.0f)) /
                                 (partMass * V * partDeltaV * sqrtf(pi) * (1 + erff((partVb - V) / partDeltaV)));
                }

                ddt[4] += mu * B_px * cdwdt * dxdt1;
                ddt[4] += mu * B_py * cdwdt * dydt1;
                ddt[4] += partMass * vec1[3] * cdwdt * dvdt1;
            }

            if constexpr (std::is_same_v<nonlinear, trueType>)
                ddt[4] *= (dis - vec1[4]);
            else
                ddt[4] *= dis;
        };

        auto interpRK4 = [&]() {
            for (int index = 0; index < 2; index++)
                coes[index] = (hx[index] + sx[index] * dx);
            FieldGather1d2d<2>(qId, coes, pic1d, q, q_px, na, na_px, nb, nb_px, ni, ni_px, ta, ta_px, tb, tb_px, ti,
                               ti_px);

            for (int index = 0; index < 2; index++)
                coes[index + 2] = coes[index];

            for (int index = 0; index < 4; index++)
                coes[index] *= (hy[index] + sy[index] * dy);
            FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py, gcovxy, gcovyy, gcovyz, gcovxy_py,
                               gcovyy_px, gcovyz_px, gcovyz_py, gconxx, gconxy, gconyy, R0, Z0);

            m2e = cm * partMass / partChar;
            mu2e = cm * mu / partChar;
            E = partMass * vec1[3] * vec1[3] / 2 + mu * B;
            if constexpr (std::is_same_v<picReal, double>) {
                radius = cm / partChar * sqrt(2.0 * mu * partMass / B);
                V = sqrt(2.0 * E / partMass);
            } else {
                radius = cm / partChar * sqrtf(2.0f * mu * partMass / B);
                V = sqrtf(2.0f * E / partMass);
            }

            rho = rho0 + vec1[0] * drho;
            bcony = 2 * psitmax * drho * rho / (q * J * B);

            invJ = 1 / J;
            invB = 1 / B;
            invQ = 1 / q;
            invRho = 1 / rho;
            invM2e = 1 / m2e;
            bconyOverJ = bcony * invJ;

            bx = bcony * gcovxy;
            by = bcony * gcovyy;
            bz = bcony * gcovyz;

            cx = bconyOverJ * (gcovyz_py - gcovyz * (J_py * invJ + B_py * invB));
            cy = -bconyOverJ *
                 (gcovyz_px - gcovyz * (J_px * invJ + B_px * invB) + gcovyz * (drho * invRho - q_px * invQ));
            cz = bconyOverJ * (gcovyy_px - gcovyy * (J_px * invJ + B_px * invB) - gcovxy_py +
                               gcovxy * (J_py * invJ + B_py * invB) + gcovyy * (drho * invRho - q_px * invQ));

            dxy = bconyOverJ * gcovyz;
            dxz = bconyOverJ * gcovyy;
            dyz = bconyOverJ * gcovxy;

            dxB = dxy * B_py;
            dyB = -dxy * B_px;
            dzB = dxz * B_px - dyz * B_py;

            Bstarx = cx * m2e * vec1[3];
            Bstary = cy * m2e * vec1[3] + B * bcony;
            Bstarz = cz * m2e * vec1[3];

            Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;

            if constexpr (std::is_same_v<nonlinear, trueType>) {

                invBstar = 1 / Bstar;
                dxdt1 = invBstar * (vec1[3] * Bstarx - mu2e * dxB);
                dydt1 = invBstar * (vec1[3] * Bstary - mu2e * dyB);
                dvdt1 = -invM2e * invBstar * mu2e * (Bstarx * B_px + Bstary * B_py);
            }

            if constexpr (std::is_same_v<picReal, double>)
                angle = acos(gconxy / sqrt(gconxx * gconyy));
            else
                angle = acosf(gconxy / sqrtf(gconxx * gconyy));

            if (i == gridNx - 2) {
                tileId = (j * cellNx + i - 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            } else {
                tileId = (j * cellNx + i + 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            }

            tileId = ((j + 1) * cellNx + i) * tileStride + 64;
            FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
            if constexpr (std::is_same_v<picReal, double>)
                gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;
            else
                gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;

            R0 = B_px;
            Z0 = B_py;
            avecxdxA = 0;
            avecydyA = 0;
            aveczdzA = 0;
            avedxPhi = 0;
            avedyPhi = 0;
            avedzPhi = 0;
            avePhipx = 0;
            avePhipy = 0;
            avePhipz = 0;
            aveAptbx = 0;
            aveAptby = 0;
            aveAptbz = 0;

            halfAngle = angle / 2;
            if constexpr (std::is_same_v<picReal, double>) {
                sincos(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sin(angle);
            } else {
                sincosf(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sinf(angle);
            }
            invSinAngle = 1 / sinAngle;
            gyroTheta = 2 * pi / gyroNums;

            for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

                if constexpr (std::is_same_v<picReal, double>)
                    sincos(gyroId * gyroTheta, &sinA, &cosA);
                else
                    sincosf(gyroId * gyroTheta, &sinA, &cosA);

                gyroX = vec1[0] + (sinA * cosHalfAngle + cosA * sinHalfAngle) * invSinAngle * gyroDx;
                gyroY = vec1[1] + (sinA * cosHalfAngle - cosA * sinHalfAngle) * invSinAngle * gyroDy;

                if (gyroX < 0 || gyroX >= 1)
                    continue;

                li = (gyroX - xbeg) / picGridDx;
                if constexpr (std::is_same_v<picReal, double>)
                    i = __double2int_rd(li);
                else
                    i = __float2int_rd(li);
                dx = li - i;
                qId = i * qStride;
                coes[0] = hx[0] + sx[0] * dx;
                coes[1] = hx[1] + sx[1] * dx;
                FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

                gyroZ = vec1[2] - q * (gyroY - vec1[1]) - vec1[1] * (gyroX - q) - (gyroY - vec1[1]) * (gyroX - q);

                if constexpr (std::is_same_v<picReal, double>) {
                    gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
                } else {
                    gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
                }

                lj = (gyroY - ybeg) / picGridDy;
                lk = (gyroZ - zbeg) / picGridDz;

                if constexpr (std::is_same_v<picReal, double>) {
                    j = __double2int_rd(lj);
                    k = __double2int_rd(lk);
                } else {
                    j = __float2int_rd(lj);
                    k = __float2int_rd(lk);
                }

                dy = lj - j;
                dz = lk - k;
                tileId = (j * cellNx + i) * tileStride;
                cellId = (j * cellNxz + i * cellNz + k) * cellStride;

                for (int index = 0; index < 2; index++)
                    coes[index + 2] = coes[index];

                for (int index = 0; index < 4; index++)
                    coes[index] *= (hy[index] + sy[index] * dy);

                for (int index = 0; index < 4; index++)
                    coes[index + 4] = coes[index];

                for (int index = 0; index < 8; index++)
                    coes[index] *= (hz[index] + sz[index] * dz);

                for (int index = 0; index < 8; index++)
                    APhiApt[index] = 0;

                FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

                avecxdxA += cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
                avecydyA += cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
                aveczdzA += cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

                avedxPhi += dxy * APhiApt[5] - dxz * APhiApt[6];
                avedyPhi += dyz * APhiApt[6] - dxy * APhiApt[4];
                avedzPhi += dxz * APhiApt[4] - dyz * APhiApt[5];

                avePhipx += APhiApt[4];
                avePhipy += APhiApt[5];
                avePhipz += APhiApt[6];

                aveAptbx += APhiApt[7] * bx;
                aveAptby += APhiApt[7] * by;
                aveAptbz += APhiApt[7] * bz;
            }

            avecxdxA /= gyroNums;
            avecydyA /= gyroNums;
            aveczdzA /= gyroNums;

            avedxPhi /= gyroNums;
            avedyPhi /= gyroNums;
            avedzPhi /= gyroNums;

            avePhipx /= gyroNums;
            avePhipy /= gyroNums;
            avePhipz /= gyroNums;

            aveAptbx /= gyroNums;
            aveAptby /= gyroNums;
            aveAptbz /= gyroNums;

            if constexpr (std::is_same_v<nonlinear, trueType>) {
                Bstarx += avecxdxA;
                Bstary += avecydyA;
                Bstarz += aveczdzA;

                Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
                invBstar = 1 / Bstar;

                ddt[0] = invBstar * (vec1[3] * Bstarx - avedxPhi - mu2e * dxB);
                ddt[1] = invBstar * (vec1[3] * Bstary - avedyPhi - mu2e * dyB);
                ddt[2] = invBstar * (vec1[3] * Bstarz - avedzPhi - mu2e * dzB);
                ddt[3] = -invM2e * invBstar *
                         (Bstarx * (avePhipx + aveAptbx + mu2e * R0) + Bstary * (avePhipy + aveAptby + mu2e * Z0) +
                          Bstarz * (avePhipz + aveAptbz));

                dxdt1 = ddt[0] - dxdt1;
                dydt1 = ddt[1] - dydt1;
                dvdt1 = ddt[3] - dvdt1;
            } else {
                Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
                invBstar = 1 / Bstar;

                ddt[0] = invBstar * (vec1[3] * Bstarx - mu2e * dxB);
                ddt[1] = invBstar * (vec1[3] * Bstary - mu2e * dyB);
                ddt[2] = invBstar * (vec1[3] * Bstarz - mu2e * dzB);
                ddt[3] = -invM2e * invBstar * mu2e * (Bstarx * R0 + Bstary * Z0);

                dxdt1 = invBstar * (vec1[3] * avecxdxA - avedxPhi);
                dydt1 = invBstar * (vec1[3] * avecydyA - avedyPhi);
                dvdt1 = -invM2e * invBstar *
                        (Bstarx * (avePhipx + aveAptbx) + Bstary * (avePhipy + aveAptby) +
                         Bstarz * (avePhipz + aveAptbz) + mu2e * (avecxdxA * R0 + avecydyA * Z0));
            }

            dfdt_XVpara();
        };

        auto advanceRK4 = [&](picReal bRK, picReal cRK) {
            for (int index = 0; index < 5; index++)
                vec2[index] += ddt[index] * picGridDt * ratioDt * bRK;

            flag = vec0[0] + ddt[0] * picGridDt * ratioDt * cRK;
            if (flag < xbeg || flag >= xend)
                illegal = 1;
            if (!illegal)
                for (int index = 0; index < 5; index++)
                    vec1[index] = vec0[index] + ddt[index] * picGridDt * ratioDt * cRK;
            else
                for (int index = 0; index < 5; index++)
                    vec1[index] = vec0[index];

            li = (vec1[0] - xbeg) / picGridDx;
            if constexpr (std::is_same_v<picReal, double>)
                i = __double2int_rd(li);
            else
                i = __float2int_rd(li);
            dx = li - i;
            qId = i * qStride;
            coes[0] = hx[0] + sx[0] * dx;
            coes[1] = hx[1] + sx[1] * dx;
            FieldGather1d2d<2>(qId, coes, pic1d, q);

            if constexpr (std::is_same_v<picReal, double>) {
                vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
                vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
                vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
            } else {
                vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
                vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
                vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
            }

            lj = (vec1[1] - ybeg) / picGridDy;
            lk = (vec1[2] - zbeg) / picGridDz;

            if constexpr (std::is_same_v<picReal, double>) {
                j = __double2int_rd(lj);
                k = __double2int_rd(lk);
            } else {
                j = __float2int_rd(lj);
                k = __float2int_rd(lk);
            }

            dy = lj - j;
            dz = lk - k;
            tileId = (j * cellNx + i) * tileStride;
            cellId = (j * cellNxz + i * cellNz + k) * cellStride;
        };

        /*--------------------------------------1st RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 6, picReal(1) / 2);

        /*--------------------------------------2nd RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 3, picReal(1) / 2);

        /*--------------------------------------3rd RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 3, picReal(1));

        /*--------------------------------------4th RK4---------------------------------------*/

        interpRK4();

        for (int index = 0; index < 5; index++)
            vec2[index] += ddt[index] * picGridDt * ratioDt / 6;

        if (vec2[0] < xbeg || vec2[0] >= xend)
            illegal = 1;
        if (illegal) {
            for (int index = 0; index < 5; index++)
                vec2[index] = vec0[index];
            vec2[1] = -vec0[1];
            vec2[4] = 0;
        }

        li = (vec2[0] - xbeg) / picGridDx;
        if constexpr (std::is_same_v<picReal, double>)
            i = __double2int_rd(li);
        else
            i = __float2int_rd(li);
        dx = li - i;
        qId = i * qStride;
        coes[0] = hx[0] + sx[0] * dx;
        coes[1] = hx[1] + sx[1] * dx;
        FieldGather1d2d<2>(qId, coes, pic1d, q);

        if constexpr (std::is_same_v<picReal, double>) {
            vec2[2] = vec2[2] + q * floor((vec2[1] - yori) / yrange) * yrange;
            vec2[1] = vec2[1] - floor((vec2[1] - yori) / yrange) * yrange;
            vec2[2] = vec2[2] - floor((vec2[2] - zori) / zrange) * zrange;
        } else {
            vec2[2] = vec2[2] + q * floorf((vec2[1] - yori) / yrange) * yrange;
            vec2[1] = vec2[1] - floorf((vec2[1] - yori) / yrange) * yrange;
            vec2[2] = vec2[2] - floorf((vec2[2] - zori) / zrange) * zrange;
        }

        lj = (vec2[1] - ybeg) / picGridDy;
        lk = (vec2[2] - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dy = lj - j;
        dz = lk - k;
        cellId = j * cellNxz + i * cellNz + k;

        for (int index = 0; index < 2; index++)
            coes[index + 2] = coes[index];

        for (int index = 0; index < 4; index++)
            coes[index] *= (hy[index] + sy[index] * dy);

        tileId = (j * cellNx + i) * tileStride + 4;
        FieldGather1d2d<4>(tileId, coes, pic2d, B);
        tileId = (j * cellNx + i) * tileStride + 52;
        FieldGather1d2d<4>(tileId, coes, pic2d, gconxx, gconxy, gconyy, R0, Z0);

        if constexpr (std::is_same_v<picReal, double>)
            radius = cm / partChar * sqrt(2.0 * mu * partMass / B);
        else
            radius = cm / partChar * sqrtf(2.0f * mu * partMass / B);

        if constexpr (std::is_same_v<picReal, double>)
            angle = acos(gconxy / sqrt(gconxx * gconyy));
        else
            angle = acosf(gconxy / sqrtf(gconxx * gconyy));

        if (i == gridNx - 2) {
            tileId = (j * cellNx + i - 1) * tileStride + 64;
            FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
            if constexpr (std::is_same_v<picReal, double>)
                gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            else
                gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
        } else {
            tileId = (j * cellNx + i + 1) * tileStride + 64;
            FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
            if constexpr (std::is_same_v<picReal, double>)
                gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            else
                gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
        }

        tileId = ((j + 1) * cellNx + i) * tileStride + 64;
        FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
        if constexpr (std::is_same_v<picReal, double>)
            gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;
        else
            gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;

        halfAngle = angle / 2;
        if constexpr (std::is_same_v<picReal, double>) {
            sincos(halfAngle, &sinHalfAngle, &cosHalfAngle);
            sinAngle = sin(angle);
        } else {
            sincosf(halfAngle, &sinHalfAngle, &cosHalfAngle);
            sinAngle = sinf(angle);
        }
        invSinAngle = 1 / sinAngle;
        gyroTheta = 2 * pi / gyroNums;

        for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

            if constexpr (std::is_same_v<picReal, double>)
                sincos(gyroId * gyroTheta, &sinA, &cosA);
            else
                sincosf(gyroId * gyroTheta, &sinA, &cosA);

            gyroX = vec2[0] + (sinA * cosHalfAngle + cosA * sinHalfAngle) * invSinAngle * gyroDx;
            gyroY = vec2[1] + (sinA * cosHalfAngle - cosA * sinHalfAngle) * invSinAngle * gyroDy;

            if (gyroX < 0 || gyroX >= 1)
                continue;

            li = (gyroX - xbeg) / picGridDx;
            if constexpr (std::is_same_v<picReal, double>)
                i = __double2int_rd(li);
            else
                i = __float2int_rd(li);
            dx = li - i;
            qId = i * qStride;
            coes[0] = hx[0] + sx[0] * dx;
            coes[1] = hx[1] + sx[1] * dx;
            FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

            gyroZ = vec2[2] - q * (gyroY - vec2[1]) - vec2[1] * (gyroX - q);

            if constexpr (std::is_same_v<picReal, double>) {
                gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
                gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
                gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
            } else {
                gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
                gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
                gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
            }

            lj = (gyroY - ybeg) / picGridDy;
            lk = (gyroZ - zbeg) / picGridDz;

            if constexpr (std::is_same_v<picReal, double>) {
                j = __double2int_rd(lj);
                k = __double2int_rd(lk);
            } else {
                j = __float2int_rd(lj);
                k = __float2int_rd(lk);
            }

            dy = lj - j;
            dz = lk - k;
            tileId = (j * cellNx + i) * tileStride;

            for (int index = 0; index < 2; index++)
                coes[index + 2] = coes[index];

            for (int index = 0; index < 4; index++)
                coes[index] *= (hy[index] + sy[index] * dy);
            FieldGather1d2d<4>(tileId, coes, pic2d, J, B);

            disP = (partMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * partConst / 2 / gyroNums;

            for (int index = 0; index < 4; index++)
                coes[index + 4] = coes[index];

            for (int index = 0; index < 8; index++)
                coes[index] *= (hz[index] + sz[index] * dz);

            if (i == 0) {

                for (int index = 0; index < 4; index++)
                    coes[2 * index] = 0;
            } else if (i == gridNx - 2) {

                for (int index = 0; index < 4; index++)
                    coes[2 * index + 1] = 0;
            }

            k = (k - gridGhost + gridNz) % gridNz;
            qId = j * gridNxz + i * gridNz + k;

            i = gridNz;
            j = gridNxz;
            k = (k + 1 + gridNz) % gridNz - k;

            atomicAdd(&dP_mid[qId], static_cast<mhdReal>(coes[0] * disP));
            atomicAdd(&dP_mid[qId + i], static_cast<mhdReal>(coes[1] * disP));
            atomicAdd(&dP_mid[qId + j], static_cast<mhdReal>(coes[2] * disP));
            atomicAdd(&dP_mid[qId + i + j], static_cast<mhdReal>(coes[3] * disP));
            atomicAdd(&dP_mid[qId + k], static_cast<mhdReal>(coes[4] * disP));
            atomicAdd(&dP_mid[qId + i + k], static_cast<mhdReal>(coes[5] * disP));
            atomicAdd(&dP_mid[qId + j + k], static_cast<mhdReal>(coes[6] * disP));
            atomicAdd(&dP_mid[qId + i + j + k], static_cast<mhdReal>(coes[7] * disP));

            if constexpr (std::is_same_v<QNeutrality, trueType>) {
                picReal disN =
                    vec2[4] / J * partConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) * l0 * l0 * l0 / gyroNums;

                atomicAdd(&dN_mid[qId], static_cast<mhdReal>(coes[0] * disN));
                atomicAdd(&dN_mid[qId + i], static_cast<mhdReal>(coes[1] * disN));
                atomicAdd(&dN_mid[qId + j], static_cast<mhdReal>(coes[2] * disN));
                atomicAdd(&dN_mid[qId + i + j], static_cast<mhdReal>(coes[3] * disN));
                atomicAdd(&dN_mid[qId + k], static_cast<mhdReal>(coes[4] * disN));
                atomicAdd(&dN_mid[qId + i + k], static_cast<mhdReal>(coes[5] * disN));
                atomicAdd(&dN_mid[qId + j + k], static_cast<mhdReal>(coes[6] * disN));
                atomicAdd(&dN_mid[qId + i + j + k], static_cast<mhdReal>(coes[7] * disN));
            }
        }

        for (int index = 0; index < 7; index++)
            pic_keys_in[picId + index * picDev] = cellId;
        pic_values_in[picId + 0 * picDev] = vec2[0];
        pic_values_in[picId + 1 * picDev] = vec2[1];
        pic_values_in[picId + 2 * picDev] = vec2[2];
        pic_values_in[picId + 3 * picDev] = vec2[3];
        pic_values_in[picId + 4 * picDev] = vec2[4];
    }
}

template <picType particle, typename mhdReal, typename picReal>
__global__ void PICDiagDensity(picReal* __restrict__ pic2d, int* __restrict__ pic_keys_in,
                               picReal* __restrict__ pic_values_in, mhdReal* __restrict__ pic_density) {

    int i, j;
    int tileId, picId;

    picReal J, li, lj, dx, dy, dis;
    picReal vec0[3] = {};
    picReal coes[4] = {};

    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        vec0[0] = pic_values_in[picId + 0 * picDev];
        vec0[1] = pic_values_in[picId + 1 * picDev];
        vec0[2] = pic_values_in[picId + 4 * picDev];

        li = (vec0[0] - xbeg) / picGridDx;
        lj = (vec0[1] - ybeg) / picGridDy;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
        }

        dx = li - i;
        dy = lj - j;

        tileId = (j * cellNx + i) * tileStride;

        /*------------------------------------Diag Density------------------------------------*/

        for (int index = 0; index < 4; index++)
            coes[index] = (hx[index] + sx[index] * dx) * (hy[index] + sy[index] * dy);

        FieldGather1d2d<4>(tileId, coes, pic2d, J);

        dis = vec0[2] / J * partConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz);

        coes[0] = hx[0] + sx[0] * dx;
        coes[1] = hx[1] + sx[1] * dx;

        if (i == 0)
            coes[0] *= 2;
        else if (i == gridNx - 2)
            coes[1] *= 2;

        atomicAdd(&pic_density[i], static_cast<mhdReal>(coes[0] * dis));
        atomicAdd(&pic_density[i + 1], static_cast<mhdReal>(coes[1] * dis));
    }
}

template <picType particle, typename mhdReal, typename picReal>
__global__ void PICDiagPhaseDeltaF(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d,
                                   picReal* __restrict__ pic_values_in, mhdReal* __restrict__ phaseDeltaF) {

    int i, j, ie, ip, il;
    int qId, tileId, picId;

    picReal x, y, vp, pw, mu;
    picReal li, lj, dx, dy;
    picReal q, psip, J, B, gcovyz;
    picReal E, Pphi, Lambda, rho;
    picReal le, lp, ll, de, dp, dl;
    picReal coes[8] = {};

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partChar = (particle == Ion) ? IonChar : (particle == Alpha) ? AlphaChar : BeamChar;
    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;

    const picReal minE = (particle == Ion)     ? IonEPphiLambda[0]
                         : (particle == Alpha) ? AlphaEPphiLambda[0]
                                               : BeamEPphiLambda[0];
    const picReal maxE = (particle == Ion)     ? IonEPphiLambda[1]
                         : (particle == Alpha) ? AlphaEPphiLambda[1]
                                               : BeamEPphiLambda[1];
    const picReal minPphi = (particle == Ion)     ? IonEPphiLambda[2]
                            : (particle == Alpha) ? AlphaEPphiLambda[2]
                                                  : BeamEPphiLambda[2];
    const picReal maxPphi = (particle == Ion)     ? IonEPphiLambda[3]
                            : (particle == Alpha) ? AlphaEPphiLambda[3]
                                                  : BeamEPphiLambda[3];
    const picReal minLambda = (particle == Ion)     ? IonEPphiLambda[4]
                              : (particle == Alpha) ? AlphaEPphiLambda[4]
                                                    : BeamEPphiLambda[4];
    const picReal maxLambda = (particle == Ion)     ? IonEPphiLambda[5]
                              : (particle == Alpha) ? AlphaEPphiLambda[5]
                                                    : BeamEPphiLambda[5];

    const picReal dE = (maxE - minE) / (gridE - 1);
    const picReal dPphi = (maxPphi - minPphi) / (gridPphi - 1);
    const picReal dLambda = (maxLambda - minLambda) / (gridLambda - 1);

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        x = pic_values_in[picId + 0 * picDev];
        y = pic_values_in[picId + 1 * picDev];
        vp = pic_values_in[picId + 3 * picDev];
        pw = pic_values_in[picId + 4 * picDev];
        mu = pic_values_in[picId + 6 * picDev];

        li = (x - xbeg) / picGridDx;
        lj = (y - ybeg) / picGridDy;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
        }

        dx = li - i;
        dy = lj - j;

        /*--------------------------------1D gather: q, psip----------------------------------*/

        coes[0] = hx[0] + sx[0] * dx;
        coes[1] = hx[1] + sx[1] * dx;

        qId = i * qStride;
        FieldGather1d2d<2>(qId, coes, pic1d, q);
        qId = i * qStride + 28;
        FieldGather1d2d<2>(qId, coes, pic1d, psip);

        /*------------------------------2D gather: J, B, gcovyz-------------------------------*/

        for (int index = 0; index < 2; index++)
            coes[index + 2] = coes[index];
        for (int index = 0; index < 4; index++)
            coes[index] *= (hy[index] + sy[index] * dy);

        tileId = (j * cellNx + i) * tileStride;
        FieldGather1d2d<4>(tileId, coes, pic2d, J, B);
        tileId = (j * cellNx + i) * tileStride + 32;
        FieldGather1d2d<4>(tileId, coes, pic2d, gcovyz);

        /*----------------------------(E, Pphi, Lambda) of marker-----------------------------*/

        E = partMass * vp * vp / 2 + mu * B;
        rho = rho0 + x * drho;
        Pphi = cm * partMass * vp * 2 * psitmax * drho * rho * gcovyz / (q * J * B) - partChar * psip;
        Lambda = mu / E;

        /*----------------------------phase-space index + clamping----------------------------*/

        le = (E - minE) / dE;
        lp = (Pphi - minPphi) / dPphi;
        ll = (Lambda - minLambda) / dLambda;

        if constexpr (std::is_same_v<picReal, double>) {
            ie = __double2int_rd(le);
            ip = __double2int_rd(lp);
            il = __double2int_rd(ll);
        } else {
            ie = __float2int_rd(le);
            ip = __float2int_rd(lp);
            il = __float2int_rd(ll);
        }

        if (ie < 0 || ie >= gridE - 1)
            continue;
        if (ip < 0 || ip >= gridPphi - 1)
            continue;
        if (il < 0 || il >= gridLambda - 1)
            continue;

        de = le - ie;
        dp = lp - ip;
        dl = ll - il;

        /*-----------------------phase-space 8-corner trilinear weights-----------------------*/

        for (int index = 0; index < 8; index++)
            coes[index] = (hx[index] + sx[index] * de) * (hy[index] + sy[index] * dp) * (hz[index] + sz[index] * dl);

        /*---------------------------8-corner trilinear deposition----------------------------*/

        const int strideE = gridPphi * gridLambda;
        const int strideP = gridLambda;
        const int index = ie * strideE + ip * strideP + il;

        atomicAdd(&phaseDeltaF[index], (mhdReal)(coes[0] * pw * partConst * tubes));
        atomicAdd(&phaseDeltaF[index + strideE], (mhdReal)(coes[1] * pw * partConst * tubes));
        atomicAdd(&phaseDeltaF[index + strideP], (mhdReal)(coes[2] * pw * partConst * tubes));
        atomicAdd(&phaseDeltaF[index + strideE + strideP], (mhdReal)(coes[3] * pw * partConst * tubes));
        atomicAdd(&phaseDeltaF[index + 1], (mhdReal)(coes[4] * pw * partConst * tubes));
        atomicAdd(&phaseDeltaF[index + strideE + 1], (mhdReal)(coes[5] * pw * partConst * tubes));
        atomicAdd(&phaseDeltaF[index + strideP + 1], (mhdReal)(coes[6] * pw * partConst * tubes));
        atomicAdd(&phaseDeltaF[index + strideE + strideP + 1], (mhdReal)(coes[7] * pw * partConst * tubes));
    }
}

template <picType particle, typename mhdReal, typename picReal>
__global__ void PICDiagPitchDeltaF(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d,
                                   picReal* __restrict__ pic_values_in, mhdReal* __restrict__ pitchDeltaF) {

    int i, j, iv, ip;
    int tileId, picId;

    picReal x, y, vp, pw, mu;
    picReal li, lj, dx, dy;
    picReal J, B;
    picReal vpara, vperp;
    picReal lvp, lvr, dvp, dvr;
    picReal coes[4] = {};

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partVmax = (particle == Ion) ? IonVmax : (particle == Alpha) ? AlphaVmax : BeamVmax;
    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;

    const picReal minVpara = (particle == Beam) ? 0 : -partVmax;
    const picReal maxVpara = partVmax;
    const picReal minVperp = 0;
    const picReal maxVperp = partVmax;

    const picReal dVpara = (maxVpara - minVpara) / (gridVpara - 1);
    const picReal dVperp = (maxVperp - minVperp) / (gridVperp - 1);

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        x = pic_values_in[picId + 0 * picDev];
        y = pic_values_in[picId + 1 * picDev];
        vp = pic_values_in[picId + 3 * picDev];
        pw = pic_values_in[picId + 4 * picDev];
        mu = pic_values_in[picId + 6 * picDev];

        li = (x - xbeg) / picGridDx;
        lj = (y - ybeg) / picGridDy;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
        }

        dx = li - i;
        dy = lj - j;

        /*-------------------------------------2D gather: J, B--------------------------------*/

        for (int index = 0; index < 4; index++)
            coes[index] = (hx[index] + sx[index] * dx) * (hy[index] + sy[index] * dy);

        tileId = (j * cellNx + i) * tileStride;
        FieldGather1d2d<4>(tileId, coes, pic2d, J, B);

        /*-----------------------------(vpara, vperp) of marker-------------------------------*/

        vpara = vp;
        if constexpr (std::is_same_v<picReal, double>)
            vperp = sqrt(2.0 * mu * B / partMass);
        else
            vperp = sqrtf(2.0f * mu * B / partMass);

        /*----------------------------pitch-space index + clamping----------------------------*/

        lvp = (vpara - minVpara) / dVpara;
        lvr = (vperp - minVperp) / dVperp;

        if constexpr (std::is_same_v<picReal, double>) {
            iv = __double2int_rd(lvp);
            ip = __double2int_rd(lvr);
        } else {
            iv = __float2int_rd(lvp);
            ip = __float2int_rd(lvr);
        }

        if (iv < 0 || iv >= gridVpara - 1)
            continue;
        if (ip < 0 || ip >= gridVperp - 1)
            continue;

        dvp = lvp - iv;
        dvr = lvr - ip;

        /*-----------------------pitch-space 4-corner bilinear weights------------------------*/

        for (int index = 0; index < 4; index++)
            coes[index] = (hx[index] + sx[index] * dvp) * (hy[index] + sy[index] * dvr);

        /*---------------------------4-corner bilinear deposition-----------------------------*/

        const int strideV = gridVperp;
        const int index = iv * strideV + ip;

        atomicAdd(&pitchDeltaF[index], (mhdReal)(coes[0] * pw * partConst * tubes));
        atomicAdd(&pitchDeltaF[index + strideV], (mhdReal)(coes[1] * pw * partConst * tubes));
        atomicAdd(&pitchDeltaF[index + 1], (mhdReal)(coes[2] * pw * partConst * tubes));
        atomicAdd(&pitchDeltaF[index + strideV + 1], (mhdReal)(coes[3] * pw * partConst * tubes));
    }
}

template <picType particle, typename mhdReal, typename picReal, typename FLRPIC, int gyroNums>
__global__ void PICDiagPhasePower(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d, picReal* __restrict__ pic3d,
                                  int* __restrict__ pic_keys_in, picReal* __restrict__ pic_values_in,
                                  mhdReal* __restrict__ phasePower) {

    int i, j, k, ie, ip, il;
    int qId, tileId, picId, cellId;

    picReal x, y, z, vp, pw, mu;
    picReal li, lj, lk, dx, dy, dz;
    picReal q, q_px, psip, J, J_px, J_py, B, B_px, B_py, gcovyz;
    picReal gcovxy, gcovyy;
    picReal gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
    picReal gconxx, gconxy, gconyy;
    picReal E, Pphi, Lambda, rho, power;
    picReal le, lp, ll, de, dp, dl;
    picReal APhiApt[8] = {};
    picReal coes[8] = {};
    picReal bx, by, bz;
    picReal bcony;
    picReal cx, cy, cz;
    picReal m2e, mu2e;
    picReal dxy, dxz, dyz;
    picReal dxB, dyB, dzB;
    picReal Bstarx, Bstary, Bstarz, Bstar;
    picReal invJ, invB, invQ, invRho, invBstar, bconyOverJ;
    picReal gyroDx, gyroDy;
    picReal gyroX, gyroY, gyroZ;
    picReal R0, Z0, R1, Z1, angle, radius;
    picReal halfAngle, sinAngle, sinHalfAngle, cosHalfAngle, invSinAngle, gyroTheta, sinA, cosA;
    picReal avecxdxA, avecydyA, aveczdzA;
    picReal avedxPhi, avedyPhi, avedzPhi;
    picReal avePhipx, avePhipy, avePhipz;
    picReal aveAptbx, aveAptby, aveAptbz;
    picReal dxdt, dydt, dzdt;
    picReal singleParticlePower, deltaParticleNumber;

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partChar = (particle == Ion) ? IonChar : (particle == Alpha) ? AlphaChar : BeamChar;
    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;

    const picReal minE = (particle == Ion)     ? IonEPphiLambda[0]
                         : (particle == Alpha) ? AlphaEPphiLambda[0]
                                               : BeamEPphiLambda[0];
    const picReal maxE = (particle == Ion)     ? IonEPphiLambda[1]
                         : (particle == Alpha) ? AlphaEPphiLambda[1]
                                               : BeamEPphiLambda[1];
    const picReal minPphi = (particle == Ion)     ? IonEPphiLambda[2]
                            : (particle == Alpha) ? AlphaEPphiLambda[2]
                                                  : BeamEPphiLambda[2];
    const picReal maxPphi = (particle == Ion)     ? IonEPphiLambda[3]
                            : (particle == Alpha) ? AlphaEPphiLambda[3]
                                                  : BeamEPphiLambda[3];
    const picReal minLambda = (particle == Ion)     ? IonEPphiLambda[4]
                              : (particle == Alpha) ? AlphaEPphiLambda[4]
                                                    : BeamEPphiLambda[4];
    const picReal maxLambda = (particle == Ion)     ? IonEPphiLambda[5]
                              : (particle == Alpha) ? AlphaEPphiLambda[5]
                                                    : BeamEPphiLambda[5];

    const picReal dE = (maxE - minE) / (gridE - 1);
    const picReal dPphi = (maxPphi - minPphi) / (gridPphi - 1);
    const picReal dLambda = (maxLambda - minLambda) / (gridLambda - 1);

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        cellId = pic_keys_in[picId];
        x = pic_values_in[picId + 0 * picDev];
        y = pic_values_in[picId + 1 * picDev];
        z = pic_values_in[picId + 2 * picDev];
        vp = pic_values_in[picId + 3 * picDev];
        pw = pic_values_in[picId + 4 * picDev];
        mu = pic_values_in[picId + 6 * picDev];

        li = (x - xbeg) / picGridDx;
        lj = (y - ybeg) / picGridDy;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
        }

        dx = li - i;
        dy = lj - j;

        for (int index = 0; index < 2; index++)
            coes[index] = hx[index] + sx[index] * dx;

        qId = i * qStride;
        FieldGather1d2d<2>(qId, coes, pic1d, q);
        qId = i * qStride + 28;
        FieldGather1d2d<2>(qId, coes, pic1d, psip);

        for (int index = 0; index < 2; index++)
            coes[index + 2] = coes[index];
        for (int index = 0; index < 4; index++)
            coes[index] *= (hy[index] + sy[index] * dy);

        tileId = (j * cellNx + i) * tileStride;
        FieldGather1d2d<4>(tileId, coes, pic2d, J, B);
        tileId = (j * cellNx + i) * tileStride + 32;
        FieldGather1d2d<4>(tileId, coes, pic2d, gcovyz);

        E = partMass * vp * vp / 2 + mu * B;
        rho = rho0 + x * drho;
        Pphi = cm * partMass * vp * 2 * psitmax * drho * rho * gcovyz / (q * J * B) - partChar * psip;
        Lambda = mu / E;

        le = (E - minE) / dE;
        lp = (Pphi - minPphi) / dPphi;
        ll = (Lambda - minLambda) / dLambda;

        if constexpr (std::is_same_v<picReal, double>) {
            ie = __double2int_rd(le);
            ip = __double2int_rd(lp);
            il = __double2int_rd(ll);
        } else {
            ie = __float2int_rd(le);
            ip = __float2int_rd(lp);
            il = __float2int_rd(ll);
        }

        if (ie < 0 || ie >= gridE - 1)
            continue;
        if (ip < 0 || ip >= gridPphi - 1)
            continue;
        if (il < 0 || il >= gridLambda - 1)
            continue;

        /*-------------------------------------Diag Power-------------------------------------*/

        li = (x - xbeg) / picGridDx;
        lj = (y - ybeg) / picGridDy;
        lk = (z - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dx = li - i;
        dy = lj - j;
        dz = lk - k;

        qId = i * qStride;
        tileId = (j * cellNx + i) * tileStride;
        cellId *= cellStride;

        for (int index = 0; index < 2; index++)
            coes[index] = (hx[index] + sx[index] * dx);
        FieldGather1d2d<2>(qId, coes, pic1d, q, q_px);

        for (int index = 0; index < 2; index++)
            coes[index + 2] = coes[index];

        for (int index = 0; index < 4; index++)
            coes[index] *= (hy[index] + sy[index] * dy);
        FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py, gcovxy, gcovyy, gcovyz, gcovxy_py,
                           gcovyy_px, gcovyz_px, gcovyz_py, gconxx, gconxy, gconyy, R0, Z0);

        m2e = cm * partMass / partChar;
        mu2e = cm * mu / partChar;

        if constexpr (std::is_same_v<FLRPIC, trueType>) {
            if constexpr (std::is_same_v<picReal, double>)
                radius = cm / partChar * sqrt(2.0 * mu * partMass / B);
            else
                radius = cm / partChar * sqrtf(2.0f * mu * partMass / B);
        }

        rho = rho0 + x * drho;
        bcony = 2 * psitmax * drho * rho / (q * J * B);

        invJ = 1 / J;
        invB = 1 / B;
        invQ = 1 / q;
        invRho = 1 / rho;
        bconyOverJ = bcony * invJ;

        bx = bcony * gcovxy;
        by = bcony * gcovyy;
        bz = bcony * gcovyz;

        cx = bconyOverJ * (gcovyz_py - gcovyz * (J_py * invJ + B_py * invB));
        cy = -bconyOverJ * (gcovyz_px - gcovyz * (J_px * invJ + B_px * invB) + gcovyz * (drho * invRho - q_px * invQ));
        cz = bconyOverJ * (gcovyy_px - gcovyy * (J_px * invJ + B_px * invB) - gcovxy_py +
                           gcovxy * (J_py * invJ + B_py * invB) + gcovyy * (drho * invRho - q_px * invQ));

        dxy = bconyOverJ * gcovyz;
        dxz = bconyOverJ * gcovyy;
        dyz = bconyOverJ * gcovxy;

        dxB = dxy * B_py;
        dyB = -dxy * B_px;
        dzB = dxz * B_px - dyz * B_py;

        Bstarx = cx * m2e * vp;
        Bstary = cy * m2e * vp + B * bcony;
        Bstarz = cz * m2e * vp;

        if constexpr (std::is_same_v<FLRPIC, falseType>) {
            for (int index = 0; index < 4; index++)
                coes[index + 4] = coes[index];

            for (int index = 0; index < 8; index++)
                coes[index] *= (hz[index] + sz[index] * dz);

            for (int index = 0; index < 8; index++)
                APhiApt[index] = 0;

            FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

            avecxdxA = cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
            avecydyA = cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
            aveczdzA = cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

            avedxPhi = dxy * APhiApt[5] - dxz * APhiApt[6];
            avedyPhi = dyz * APhiApt[6] - dxy * APhiApt[4];
            avedzPhi = dxz * APhiApt[4] - dyz * APhiApt[5];

            avePhipx = APhiApt[4];
            avePhipy = APhiApt[5];
            avePhipz = APhiApt[6];

            aveAptbx = APhiApt[7] * bx;
            aveAptby = APhiApt[7] * by;
            aveAptbz = APhiApt[7] * bz;

        } else {
            if constexpr (std::is_same_v<picReal, double>)
                angle = acos(gconxy / sqrt(gconxx * gconyy));
            else
                angle = acosf(gconxy / sqrtf(gconxx * gconyy));

            if (i == gridNx - 2) {
                tileId = (j * cellNx + i - 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            } else {
                tileId = (j * cellNx + i + 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            }

            tileId = ((j + 1) * cellNx + i) * tileStride + 64;
            FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
            if constexpr (std::is_same_v<picReal, double>)
                gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;
            else
                gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;

            avecxdxA = 0;
            avecydyA = 0;
            aveczdzA = 0;
            avedxPhi = 0;
            avedyPhi = 0;
            avedzPhi = 0;
            avePhipx = 0;
            avePhipy = 0;
            avePhipz = 0;
            aveAptbx = 0;
            aveAptby = 0;
            aveAptbz = 0;

            halfAngle = angle / picReal(2);
            if constexpr (std::is_same_v<picReal, double>) {
                sincos(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sin(angle);
            } else {
                sincosf(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sinf(angle);
            }
            invSinAngle = 1 / sinAngle;
            gyroTheta = 2 * pi / gyroNums;

            for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

                if constexpr (std::is_same_v<picReal, double>)
                    sincos(gyroId * gyroTheta, &sinA, &cosA);
                else
                    sincosf(gyroId * gyroTheta, &sinA, &cosA);

                gyroX = x + (sinA * cosHalfAngle + cosA * sinHalfAngle) * invSinAngle * gyroDx;
                gyroY = y + (sinA * cosHalfAngle - cosA * sinHalfAngle) * invSinAngle * gyroDy;

                if (gyroX < 0 || gyroX >= 1)
                    continue;

                li = (gyroX - xbeg) / picGridDx;
                if constexpr (std::is_same_v<picReal, double>)
                    i = __double2int_rd(li);
                else
                    i = __float2int_rd(li);
                dx = li - i;
                qId = i * qStride;
                coes[0] = hx[0] + sx[0] * dx;
                coes[1] = hx[1] + sx[1] * dx;
                FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

                gyroZ = z - q * (gyroY - y) - y * (gyroX - q) - (gyroY - y) * (gyroX - q);

                if constexpr (std::is_same_v<picReal, double>) {
                    gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
                } else {
                    gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
                }

                lj = (gyroY - ybeg) / picGridDy;
                lk = (gyroZ - zbeg) / picGridDz;

                if constexpr (std::is_same_v<picReal, double>) {
                    j = __double2int_rd(lj);
                    k = __double2int_rd(lk);
                } else {
                    j = __float2int_rd(lj);
                    k = __float2int_rd(lk);
                }

                dy = lj - j;
                dz = lk - k;
                tileId = (j * cellNx + i) * tileStride;
                cellId = (j * cellNxz + i * cellNz + k) * cellStride;

                for (int index = 0; index < 2; index++)
                    coes[index + 2] = coes[index];

                for (int index = 0; index < 4; index++)
                    coes[index] *= (hy[index] + sy[index] * dy);

                for (int index = 0; index < 4; index++)
                    coes[index + 4] = coes[index];

                for (int index = 0; index < 8; index++)
                    coes[index] *= (hz[index] + sz[index] * dz);

                for (int index = 0; index < 8; index++)
                    APhiApt[index] = 0;

                FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

                avecxdxA += cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
                avecydyA += cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
                aveczdzA += cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

                avedxPhi += dxy * APhiApt[5] - dxz * APhiApt[6];
                avedyPhi += dyz * APhiApt[6] - dxy * APhiApt[4];
                avedzPhi += dxz * APhiApt[4] - dyz * APhiApt[5];

                avePhipx += APhiApt[4];
                avePhipy += APhiApt[5];
                avePhipz += APhiApt[6];

                aveAptbx += APhiApt[7] * bx;
                aveAptby += APhiApt[7] * by;
                aveAptbz += APhiApt[7] * bz;
            }

            avecxdxA /= gyroNums;
            avecydyA /= gyroNums;
            aveczdzA /= gyroNums;

            avedxPhi /= gyroNums;
            avedyPhi /= gyroNums;
            avedzPhi /= gyroNums;

            avePhipx /= gyroNums;
            avePhipy /= gyroNums;
            avePhipz /= gyroNums;

            aveAptbx /= gyroNums;
            aveAptby /= gyroNums;
            aveAptbz /= gyroNums;
        }

        Bstarx += avecxdxA;
        Bstary += avecydyA;
        Bstarz += aveczdzA;
        Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
        invBstar = 1 / Bstar;

        dxdt = invBstar * (vp * Bstarx - avedxPhi - mu2e * dxB);
        dydt = invBstar * (vp * Bstary - avedyPhi - mu2e * dyB);
        dzdt = invBstar * (vp * Bstarz - avedzPhi - mu2e * dzB);

        singleParticlePower =
            partChar / cm *
            ((avePhipx + aveAptbx) * dxdt + (avePhipy + aveAptby) * dydt + (avePhipz + aveAptbz) * dzdt) *
            (mp * va * va * va / l0);
        deltaParticleNumber = (pw / J * partConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va)) *
                              (J * picGridDx * picGridDy * picGridDz * l0 * l0 * l0);
        power = deltaParticleNumber * singleParticlePower;

        de = le - ie;
        dp = lp - ip;
        dl = ll - il;

        for (int index = 0; index < 8; index++)
            coes[index] = (hx[index] + sx[index] * de) * (hy[index] + sy[index] * dp) * (hz[index] + sz[index] * dl);

        const int strideE = gridPphi * gridLambda;
        const int strideP = gridLambda;
        const int index = ie * strideE + ip * strideP + il;

        atomicAdd(&phasePower[index], (mhdReal)(coes[0] * power * tubes));
        atomicAdd(&phasePower[index + strideE], (mhdReal)(coes[1] * power * tubes));
        atomicAdd(&phasePower[index + strideP], (mhdReal)(coes[2] * power * tubes));
        atomicAdd(&phasePower[index + strideE + strideP], (mhdReal)(coes[3] * power * tubes));
        atomicAdd(&phasePower[index + 1], (mhdReal)(coes[4] * power * tubes));
        atomicAdd(&phasePower[index + strideE + 1], (mhdReal)(coes[5] * power * tubes));
        atomicAdd(&phasePower[index + strideP + 1], (mhdReal)(coes[6] * power * tubes));
        atomicAdd(&phasePower[index + strideE + strideP + 1], (mhdReal)(coes[7] * power * tubes));
    }
}

template <picType particle, typename mhdReal, typename picReal, typename FLRPIC, int gyroNums>
__global__ void PICDiagPitchPower(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d, picReal* __restrict__ pic3d,
                                  int* __restrict__ pic_keys_in, picReal* __restrict__ pic_values_in,
                                  mhdReal* __restrict__ pitchPower) {

    int i, j, k, iv, ip;
    int qId, tileId, picId, cellId;

    picReal x, y, z, vp, pw, mu;
    picReal li, lj, lk, dx, dy, dz;
    picReal q, q_px, J, J_px, J_py, B, B_px, B_py;
    picReal gcovxy, gcovyy, gcovyz;
    picReal gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
    picReal gconxx, gconxy, gconyy;
    picReal vpara, vperp, power;
    picReal lvp, lvr, dvp, dvr;
    picReal APhiApt[8] = {};
    picReal coes[8] = {};
    picReal bx, by, bz;
    picReal rho, bcony;
    picReal cx, cy, cz;
    picReal m2e, mu2e;
    picReal dxy, dxz, dyz;
    picReal dxB, dyB, dzB;
    picReal Bstarx, Bstary, Bstarz, Bstar;
    picReal invJ, invB, invQ, invRho, invBstar, bconyOverJ;
    picReal gyroDx, gyroDy;
    picReal gyroX, gyroY, gyroZ;
    picReal R0, Z0, R1, Z1, angle, radius;
    picReal halfAngle, sinAngle, sinHalfAngle, cosHalfAngle, invSinAngle, gyroTheta, sinA, cosA;
    picReal avecxdxA, avecydyA, aveczdzA;
    picReal avedxPhi, avedyPhi, avedzPhi;
    picReal avePhipx, avePhipy, avePhipz;
    picReal aveAptbx, aveAptby, aveAptbz;
    picReal dxdt, dydt, dzdt;
    picReal singleParticlePower, deltaParticleNumber;

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partChar = (particle == Ion) ? IonChar : (particle == Alpha) ? AlphaChar : BeamChar;
    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;
    const picReal partVmax = (particle == Ion) ? IonVmax : (particle == Alpha) ? AlphaVmax : BeamVmax;

    const picReal minVpara = (particle == Beam) ? 0 : -partVmax;
    const picReal maxVpara = partVmax;
    const picReal minVperp = 0;
    const picReal maxVperp = partVmax;

    const picReal dVpara = (maxVpara - minVpara) / (gridVpara - 1);
    const picReal dVperp = (maxVperp - minVperp) / (gridVperp - 1);

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        cellId = pic_keys_in[picId];
        x = pic_values_in[picId + 0 * picDev];
        y = pic_values_in[picId + 1 * picDev];
        z = pic_values_in[picId + 2 * picDev];
        vp = pic_values_in[picId + 3 * picDev];
        pw = pic_values_in[picId + 4 * picDev];
        mu = pic_values_in[picId + 6 * picDev];

        li = (x - xbeg) / picGridDx;
        lj = (y - ybeg) / picGridDy;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
        }

        dx = li - i;
        dy = lj - j;

        for (int index = 0; index < 4; index++)
            coes[index] = (hx[index] + sx[index] * dx) * (hy[index] + sy[index] * dy);

        tileId = (j * cellNx + i) * tileStride;
        FieldGather1d2d<4>(tileId, coes, pic2d, J, B);

        vpara = vp;
        if constexpr (std::is_same_v<picReal, double>)
            vperp = sqrt(2.0 * mu * B / partMass);
        else
            vperp = sqrtf(2.0f * mu * B / partMass);

        lvp = (vpara - minVpara) / dVpara;
        lvr = (vperp - minVperp) / dVperp;

        if constexpr (std::is_same_v<picReal, double>) {
            iv = __double2int_rd(lvp);
            ip = __double2int_rd(lvr);
        } else {
            iv = __float2int_rd(lvp);
            ip = __float2int_rd(lvr);
        }

        if (iv < 0 || iv >= gridVpara - 1)
            continue;
        if (ip < 0 || ip >= gridVperp - 1)
            continue;

        /*-------------------------------------Diag Power-------------------------------------*/

        li = (x - xbeg) / picGridDx;
        lj = (y - ybeg) / picGridDy;
        lk = (z - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dx = li - i;
        dy = lj - j;
        dz = lk - k;

        qId = i * qStride;
        tileId = (j * cellNx + i) * tileStride;
        cellId *= cellStride;

        for (int index = 0; index < 2; index++)
            coes[index] = (hx[index] + sx[index] * dx);
        FieldGather1d2d<2>(qId, coes, pic1d, q, q_px);

        for (int index = 0; index < 2; index++)
            coes[index + 2] = coes[index];

        for (int index = 0; index < 4; index++)
            coes[index] *= (hy[index] + sy[index] * dy);
        FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py, gcovxy, gcovyy, gcovyz, gcovxy_py,
                           gcovyy_px, gcovyz_px, gcovyz_py, gconxx, gconxy, gconyy, R0, Z0);

        m2e = cm * partMass / partChar;
        mu2e = cm * mu / partChar;

        if constexpr (std::is_same_v<FLRPIC, trueType>) {
            if constexpr (std::is_same_v<picReal, double>)
                radius = cm / partChar * sqrt(2.0 * mu * partMass / B);
            else
                radius = cm / partChar * sqrtf(2.0f * mu * partMass / B);
        }

        rho = rho0 + x * drho;
        bcony = 2 * psitmax * drho * rho / (q * J * B);

        invJ = 1 / J;
        invB = 1 / B;
        invQ = 1 / q;
        invRho = 1 / rho;
        bconyOverJ = bcony * invJ;

        bx = bcony * gcovxy;
        by = bcony * gcovyy;
        bz = bcony * gcovyz;

        cx = bconyOverJ * (gcovyz_py - gcovyz * (J_py * invJ + B_py * invB));
        cy = -bconyOverJ * (gcovyz_px - gcovyz * (J_px * invJ + B_px * invB) + gcovyz * (drho * invRho - q_px * invQ));
        cz = bconyOverJ * (gcovyy_px - gcovyy * (J_px * invJ + B_px * invB) - gcovxy_py +
                           gcovxy * (J_py * invJ + B_py * invB) + gcovyy * (drho * invRho - q_px * invQ));

        dxy = bconyOverJ * gcovyz;
        dxz = bconyOverJ * gcovyy;
        dyz = bconyOverJ * gcovxy;

        dxB = dxy * B_py;
        dyB = -dxy * B_px;
        dzB = dxz * B_px - dyz * B_py;

        Bstarx = cx * m2e * vp;
        Bstary = cy * m2e * vp + B * bcony;
        Bstarz = cz * m2e * vp;

        if constexpr (std::is_same_v<FLRPIC, falseType>) {
            for (int index = 0; index < 4; index++)
                coes[index + 4] = coes[index];

            for (int index = 0; index < 8; index++)
                coes[index] *= (hz[index] + sz[index] * dz);

            for (int index = 0; index < 8; index++)
                APhiApt[index] = 0;

            FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

            avecxdxA = cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
            avecydyA = cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
            aveczdzA = cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

            avedxPhi = dxy * APhiApt[5] - dxz * APhiApt[6];
            avedyPhi = dyz * APhiApt[6] - dxy * APhiApt[4];
            avedzPhi = dxz * APhiApt[4] - dyz * APhiApt[5];

            avePhipx = APhiApt[4];
            avePhipy = APhiApt[5];
            avePhipz = APhiApt[6];

            aveAptbx = APhiApt[7] * bx;
            aveAptby = APhiApt[7] * by;
            aveAptbz = APhiApt[7] * bz;

        } else {
            if constexpr (std::is_same_v<picReal, double>)
                angle = acos(gconxy / sqrt(gconxx * gconyy));
            else
                angle = acosf(gconxy / sqrtf(gconxx * gconyy));

            if (i == gridNx - 2) {
                tileId = (j * cellNx + i - 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            } else {
                tileId = (j * cellNx + i + 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            }

            tileId = ((j + 1) * cellNx + i) * tileStride + 64;
            FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
            if constexpr (std::is_same_v<picReal, double>)
                gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;
            else
                gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;

            avecxdxA = 0;
            avecydyA = 0;
            aveczdzA = 0;
            avedxPhi = 0;
            avedyPhi = 0;
            avedzPhi = 0;
            avePhipx = 0;
            avePhipy = 0;
            avePhipz = 0;
            aveAptbx = 0;
            aveAptby = 0;
            aveAptbz = 0;

            halfAngle = angle / picReal(2);
            if constexpr (std::is_same_v<picReal, double>) {
                sincos(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sin(angle);
            } else {
                sincosf(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sinf(angle);
            }
            invSinAngle = 1 / sinAngle;
            gyroTheta = 2 * pi / gyroNums;

            for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

                if constexpr (std::is_same_v<picReal, double>)
                    sincos(gyroId * gyroTheta, &sinA, &cosA);
                else
                    sincosf(gyroId * gyroTheta, &sinA, &cosA);

                gyroX = x + (sinA * cosHalfAngle + cosA * sinHalfAngle) * invSinAngle * gyroDx;
                gyroY = y + (sinA * cosHalfAngle - cosA * sinHalfAngle) * invSinAngle * gyroDy;

                if (gyroX < 0 || gyroX >= 1)
                    continue;

                li = (gyroX - xbeg) / picGridDx;
                if constexpr (std::is_same_v<picReal, double>)
                    i = __double2int_rd(li);
                else
                    i = __float2int_rd(li);
                dx = li - i;
                qId = i * qStride;
                coes[0] = hx[0] + sx[0] * dx;
                coes[1] = hx[1] + sx[1] * dx;
                FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

                gyroZ = z - q * (gyroY - y) - y * (gyroX - q) - (gyroY - y) * (gyroX - q);

                if constexpr (std::is_same_v<picReal, double>) {
                    gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
                } else {
                    gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
                }

                lj = (gyroY - ybeg) / picGridDy;
                lk = (gyroZ - zbeg) / picGridDz;

                if constexpr (std::is_same_v<picReal, double>) {
                    j = __double2int_rd(lj);
                    k = __double2int_rd(lk);
                } else {
                    j = __float2int_rd(lj);
                    k = __float2int_rd(lk);
                }

                dy = lj - j;
                dz = lk - k;
                tileId = (j * cellNx + i) * tileStride;
                cellId = (j * cellNxz + i * cellNz + k) * cellStride;

                for (int index = 0; index < 2; index++)
                    coes[index + 2] = coes[index];

                for (int index = 0; index < 4; index++)
                    coes[index] *= (hy[index] + sy[index] * dy);

                for (int index = 0; index < 4; index++)
                    coes[index + 4] = coes[index];

                for (int index = 0; index < 8; index++)
                    coes[index] *= (hz[index] + sz[index] * dz);

                for (int index = 0; index < 8; index++)
                    APhiApt[index] = 0;

                FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

                avecxdxA += cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
                avecydyA += cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
                aveczdzA += cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

                avedxPhi += dxy * APhiApt[5] - dxz * APhiApt[6];
                avedyPhi += dyz * APhiApt[6] - dxy * APhiApt[4];
                avedzPhi += dxz * APhiApt[4] - dyz * APhiApt[5];

                avePhipx += APhiApt[4];
                avePhipy += APhiApt[5];
                avePhipz += APhiApt[6];

                aveAptbx += APhiApt[7] * bx;
                aveAptby += APhiApt[7] * by;
                aveAptbz += APhiApt[7] * bz;
            }

            avecxdxA /= gyroNums;
            avecydyA /= gyroNums;
            aveczdzA /= gyroNums;

            avedxPhi /= gyroNums;
            avedyPhi /= gyroNums;
            avedzPhi /= gyroNums;

            avePhipx /= gyroNums;
            avePhipy /= gyroNums;
            avePhipz /= gyroNums;

            aveAptbx /= gyroNums;
            aveAptby /= gyroNums;
            aveAptbz /= gyroNums;
        }

        Bstarx += avecxdxA;
        Bstary += avecydyA;
        Bstarz += aveczdzA;
        Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
        invBstar = 1 / Bstar;

        dxdt = invBstar * (vp * Bstarx - avedxPhi - mu2e * dxB);
        dydt = invBstar * (vp * Bstary - avedyPhi - mu2e * dyB);
        dzdt = invBstar * (vp * Bstarz - avedzPhi - mu2e * dzB);

        singleParticlePower =
            partChar / cm *
            ((avePhipx + aveAptbx) * dxdt + (avePhipy + aveAptby) * dydt + (avePhipz + aveAptbz) * dzdt) *
            (mp * va * va * va / l0);
        deltaParticleNumber = (pw / J * partConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va)) *
                              (J * picGridDx * picGridDy * picGridDz * l0 * l0 * l0);
        power = deltaParticleNumber * singleParticlePower;

        dvp = lvp - iv;
        dvr = lvr - ip;

        for (int index = 0; index < 4; index++)
            coes[index] = (hx[index] + sx[index] * dvp) * (hy[index] + sy[index] * dvr);

        const int strideV = gridVperp;
        const int index = iv * strideV + ip;

        atomicAdd(&pitchPower[index], (mhdReal)(coes[0] * power * tubes));
        atomicAdd(&pitchPower[index + strideV], (mhdReal)(coes[1] * power * tubes));
        atomicAdd(&pitchPower[index + 1], (mhdReal)(coes[2] * power * tubes));
        atomicAdd(&pitchPower[index + strideV + 1], (mhdReal)(coes[3] * power * tubes));
    }
}

template <picType particle, typename mhdReal, typename picReal, typename FLRPIC, int gyroNums>
__global__ void PICDiagDiffusivity(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d,
                                   picReal* __restrict__ pic3d, int* __restrict__ pic_keys_in,
                                   picReal* __restrict__ pic_values_in, mhdReal* __restrict__ pic_diffusivity) {

    int i, j, k;
    int qId, tileId, cellId, picId;

    picReal li, lj, lk;
    picReal coes[8] = {};

    picReal dx, dy, dz, dis, mu, dxdt;
    picReal vec0[5] = {};

    picReal q, q_px, J, J_px, J_py, B, B_px, B_py;
    picReal gcovxy, gcovyy, gcovyz;
    picReal gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
    picReal APhiApt[8] = {};

    picReal bx, by, bz;
    picReal rho, bcony;
    picReal cx, cy, cz;
    picReal m2e;
    picReal dxy, dxz, dyz;
    picReal Bstarx, Bstary, Bstarz, Bstar;
    picReal invJ, invB, invQ, invRho, invBstar, bconyOverJ;
    picReal na, na_px, nb, nb_px, ni, ni_px;

    picReal gyroDx, gyroDy;
    picReal gyroX, gyroY, gyroZ;
    picReal gconxx, gconxy, gconyy;
    picReal R0, Z0, R1, Z1, angle, radius;
    picReal halfAngle, sinAngle, sinHalfAngle, cosHalfAngle, invSinAngle, gyroTheta, sinA, cosA;
    picReal avecxdxA, avecydyA, aveczdzA;
    picReal avedxPhi, avedyPhi, avedzPhi;
    picReal avePhipx, avePhipy, avePhipz;
    picReal aveAptbx, aveAptby, aveAptbz;

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partChar = (particle == Ion) ? IonChar : (particle == Alpha) ? AlphaChar : BeamChar;
    const picReal partConst = (particle == Ion) ? IonConst : (particle == Alpha) ? AlphaConst : BeamConst;
    picReal& part_n_px = (particle == Ion ? ni_px : (particle == Alpha ? na_px : nb_px));

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
        cellId = pic_keys_in[picId];
        vec0[0] = pic_values_in[picId + 0 * picDev];
        vec0[1] = pic_values_in[picId + 1 * picDev];
        vec0[2] = pic_values_in[picId + 2 * picDev];
        vec0[3] = pic_values_in[picId + 3 * picDev];
        vec0[4] = pic_values_in[picId + 4 * picDev];
        mu = pic_values_in[picId + 6 * picDev];

        li = (vec0[0] - xbeg) / picGridDx;
        lj = (vec0[1] - ybeg) / picGridDy;
        lk = (vec0[2] - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dx = li - i;
        dy = lj - j;
        dz = lk - k;

        qId = i * qStride;
        tileId = (j * cellNx + i) * tileStride;
        cellId *= cellStride;

        auto interpRK4 = [&]() {
            for (int index = 0; index < 2; index++)
                coes[index] = (hx[index] + sx[index] * dx);
            FieldGather1d2d<2>(qId, coes, pic1d, q, q_px, na, na_px, nb, nb_px, ni, ni_px);

            for (int index = 0; index < 2; index++)
                coes[index + 2] = coes[index];

            for (int index = 0; index < 4; index++)
                coes[index] *= (hy[index] + sy[index] * dy);
            FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py, gcovxy, gcovyy, gcovyz, gcovxy_py,
                               gcovyy_px, gcovyz_px, gcovyz_py, gconxx, gconxy, gconyy, R0, Z0);

            m2e = cm * partMass / partChar;
            if constexpr (std::is_same_v<FLRPIC, trueType>) {
                if constexpr (std::is_same_v<picReal, double>)
                    radius = cm / partChar * sqrt(2.0 * mu * partMass / B);
                else
                    radius = cm / partChar * sqrtf(2.0f * mu * partMass / B);
            }

            rho = rho0 + vec0[0] * drho;
            bcony = 2 * psitmax * drho * rho / (q * J * B);

            invJ = 1 / J;
            invB = 1 / B;
            invQ = 1 / q;
            invRho = 1 / rho;
            bconyOverJ = bcony * invJ;

            bx = bcony * gcovxy;
            by = bcony * gcovyy;
            bz = bcony * gcovyz;

            cx = bconyOverJ * (gcovyz_py - gcovyz * (J_py * invJ + B_py * invB));
            cy = -bconyOverJ *
                 (gcovyz_px - gcovyz * (J_px * invJ + B_px * invB) + gcovyz * (drho * invRho - q_px * invQ));
            cz = bconyOverJ * (gcovyy_px - gcovyy * (J_px * invJ + B_px * invB) - gcovxy_py +
                               gcovxy * (J_py * invJ + B_py * invB) + gcovyy * (drho * invRho - q_px * invQ));

            dxy = bconyOverJ * gcovyz;
            dxz = bconyOverJ * gcovyy;
            dyz = bconyOverJ * gcovxy;

            Bstarx = cx * m2e * vec0[3];
            Bstary = cy * m2e * vec0[3] + B * bcony;
            Bstarz = cz * m2e * vec0[3];

            Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;

            if constexpr (std::is_same_v<FLRPIC, falseType>) {
                for (int index = 0; index < 4; index++)
                    coes[index + 4] = coes[index];

                for (int index = 0; index < 8; index++)
                    coes[index] *= (hz[index] + sz[index] * dz);

                for (int index = 0; index < 8; index++)
                    APhiApt[index] = 0;

                FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

                avecxdxA = cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
                avecydyA = cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
                aveczdzA = cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

                avedxPhi = dxy * APhiApt[5] - dxz * APhiApt[6];

                invBstar = 1 / Bstar;
                dxdt = invBstar * (vec0[3] * Bstarx - cm * mu / partChar * dxy * B_py);

                Bstarx += avecxdxA;
                Bstary += avecydyA;
                Bstarz += aveczdzA;
                Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
                invBstar = 1 / Bstar;

                dxdt = invBstar * (vec0[3] * Bstarx - avedxPhi - cm * mu / partChar * dxy * B_py) - dxdt;
                return;
            }

            if constexpr (std::is_same_v<picReal, double>)
                angle = acos(gconxy / sqrt(gconxx * gconyy));
            else
                angle = acosf(gconxy / sqrtf(gconxx * gconyy));

            if (i == gridNx - 2) {
                tileId = (j * cellNx + i - 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            } else {
                tileId = (j * cellNx + i + 1) * tileStride + 64;
                FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
                if constexpr (std::is_same_v<picReal, double>)
                    gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
                else
                    gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDx;
            }

            tileId = ((j + 1) * cellNx + i) * tileStride + 64;
            FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
            if constexpr (std::is_same_v<picReal, double>)
                gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;
            else
                gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * picGridDy;

            R0 = B_px;
            Z0 = B_py;
            avecxdxA = 0;
            avecydyA = 0;
            aveczdzA = 0;
            avedxPhi = 0;
            avedyPhi = 0;
            avedzPhi = 0;
            avePhipx = 0;
            avePhipy = 0;
            avePhipz = 0;
            aveAptbx = 0;
            aveAptby = 0;
            aveAptbz = 0;

            halfAngle = angle / picReal(2);
            if constexpr (std::is_same_v<picReal, double>) {
                sincos(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sin(angle);
            } else {
                sincosf(halfAngle, &sinHalfAngle, &cosHalfAngle);
                sinAngle = sinf(angle);
            }
            invSinAngle = 1 / sinAngle;
            gyroTheta = 2 * pi / gyroNums;

            for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

                if constexpr (std::is_same_v<picReal, double>)
                    sincos(gyroId * gyroTheta, &sinA, &cosA);
                else
                    sincosf(gyroId * gyroTheta, &sinA, &cosA);

                gyroX = vec0[0] + (sinA * cosHalfAngle + cosA * sinHalfAngle) * invSinAngle * gyroDx;
                gyroY = vec0[1] + (sinA * cosHalfAngle - cosA * sinHalfAngle) * invSinAngle * gyroDy;

                if (gyroX < 0 || gyroX >= 1)
                    continue;

                li = (gyroX - xbeg) / picGridDx;
                if constexpr (std::is_same_v<picReal, double>)
                    i = __double2int_rd(li);
                else
                    i = __float2int_rd(li);
                dx = li - i;
                qId = i * qStride;
                coes[0] = hx[0] + sx[0] * dx;
                coes[1] = hx[1] + sx[1] * dx;
                FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

                gyroZ = vec0[2] - q * (gyroY - vec0[1]) - vec0[1] * (gyroX - q) - (gyroY - vec0[1]) * (gyroX - q);

                if constexpr (std::is_same_v<picReal, double>) {
                    gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
                } else {
                    gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
                    gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
                    gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
                }

                lj = (gyroY - ybeg) / picGridDy;
                lk = (gyroZ - zbeg) / picGridDz;

                if constexpr (std::is_same_v<picReal, double>) {
                    j = __double2int_rd(lj);
                    k = __double2int_rd(lk);
                } else {
                    j = __float2int_rd(lj);
                    k = __float2int_rd(lk);
                }

                dy = lj - j;
                dz = lk - k;
                tileId = (j * cellNx + i) * tileStride;
                cellId = (j * cellNxz + i * cellNz + k) * cellStride;

                for (int index = 0; index < 2; index++)
                    coes[index + 2] = coes[index];

                for (int index = 0; index < 4; index++)
                    coes[index] *= (hy[index] + sy[index] * dy);

                for (int index = 0; index < 4; index++)
                    coes[index + 4] = coes[index];

                for (int index = 0; index < 8; index++)
                    coes[index] *= (hz[index] + sz[index] * dz);

                for (int index = 0; index < 8; index++)
                    APhiApt[index] = 0;

                FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

                avecxdxA += cx * APhiApt[0] + dxy * APhiApt[2] - dxz * APhiApt[3];
                avecydyA += cy * APhiApt[0] + dyz * APhiApt[3] - dxy * APhiApt[1];
                aveczdzA += cz * APhiApt[0] + dxz * APhiApt[1] - dyz * APhiApt[2];

                avedxPhi += dxy * APhiApt[5] - dxz * APhiApt[6];
                avedyPhi += dyz * APhiApt[6] - dxy * APhiApt[4];
                avedzPhi += dxz * APhiApt[4] - dyz * APhiApt[5];

                avePhipx += APhiApt[4];
                avePhipy += APhiApt[5];
                avePhipz += APhiApt[6];

                aveAptbx += APhiApt[7] * bx;
                aveAptby += APhiApt[7] * by;
                aveAptbz += APhiApt[7] * bz;
            }

            avecxdxA /= gyroNums;
            avecydyA /= gyroNums;
            aveczdzA /= gyroNums;

            avedxPhi /= gyroNums;
            avedyPhi /= gyroNums;
            avedzPhi /= gyroNums;

            avePhipx /= gyroNums;
            avePhipy /= gyroNums;
            avePhipz /= gyroNums;

            aveAptbx /= gyroNums;
            aveAptby /= gyroNums;
            aveAptbz /= gyroNums;

            Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
            invBstar = 1 / Bstar;
            dxdt = invBstar * (vec0[3] * Bstarx - cm * mu / partChar * dxy * B_py);

            Bstarx += avecxdxA;
            Bstary += avecydyA;
            Bstarz += aveczdzA;
            Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
            invBstar = 1 / Bstar;

            dxdt = invBstar * (vec0[3] * Bstarx - avedxPhi - cm * mu / partChar * dxy * B_py) - dxdt;
        };

        interpRK4();

        li = (vec0[0] - xbeg) / picGridDx;
        lj = (vec0[1] - ybeg) / picGridDy;
        lk = (vec0[2] - zbeg) / picGridDz;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
            k = __double2int_rd(lk);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
            k = __float2int_rd(lk);
        }

        dx = li - i;
        dy = lj - j;
        dz = lk - k;

        /*----------------------------------Diag Diffusivity----------------------------------*/

        for (int index = 0; index < 4; index++)
            coes[index] = (hx[index] + sx[index] * dx) * (hy[index] + sy[index] * dy);

        tileId = (j * cellNx + i) * tileStride;
        FieldGather1d2d<4>(tileId, coes, pic2d, J);
        tileId = (j * cellNx + i) * tileStride + 52;
        FieldGather1d2d<4>(tileId, coes, pic2d, gconxx);

        dis = -vec0[4] / J * partConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz) * dxdt * va /
              (part_n_px * gconxx * l4);

        coes[0] = hx[0] + sx[0] * dx;
        coes[1] = hx[1] + sx[1] * dx;

        if (i == 0)
            coes[0] *= 2;
        else if (i == gridNx - 2)
            coes[1] *= 2;

        atomicAdd(&pic_diffusivity[i], static_cast<mhdReal>(coes[0] * dis));
        atomicAdd(&pic_diffusivity[i + 1], static_cast<mhdReal>(coes[1] * dis));
    }
}

template <int ratioDt, picType particle, typename picReal>
__global__ void PICDiagOrbit(picReal* __restrict__ pic1d, picReal* __restrict__ pic2d,
                             picReal* __restrict__ phaseSpaceMapping) {

    int illegal;
    int i, j;
    int qId, tileId, picId;

    picReal flag;
    picReal li, lj, dx, dy;
    picReal coes[4] = {};
    picReal ddt[4] = {};
    picReal vec0[4] = {};
    picReal vec1[4] = {};
    picReal vec2[4] = {};

    picReal q, q_px, J, J_px, J_py, B, B_px, B_py;
    picReal gcovxy, gcovyy, gcovyz;
    picReal gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;

    picReal bx, by, bz;
    picReal rho, bcony;
    picReal cx, cy, cz;
    picReal m2e, mu2e;
    picReal dxy, dxz, dyz;
    picReal dxB, dyB, dzB;
    picReal Bstarx, Bstary, Bstarz, Bstar;
    picReal invJ, invB, invQ, invRho, invBstar, invM2e, bconyOverJ;

    picReal orbit, mu, q0, frac, psip, E, Pphi, Lambda;
    picReal dtheta, dphiTotal, dphiVpara, dT, bounce;

    const picReal partMass = (particle == Ion) ? IonMass : (particle == Alpha) ? AlphaMass : BeamMass;
    const picReal partChar = (particle == Ion) ? IonChar : (particle == Alpha) ? AlphaChar : BeamChar;

    for (int id = 0; id < pptNums; id++) {

        picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;

        // orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

        orbit = phaseSpaceMapping[picId * 13 + 0];
        vec0[0] = phaseSpaceMapping[picId * 13 + 1];
        vec0[1] = phaseSpaceMapping[picId * 13 + 2];
        vec0[3] = phaseSpaceMapping[picId * 13 + 3];
        mu = phaseSpaceMapping[picId * 13 + 4];
        dtheta = phaseSpaceMapping[picId * 13 + 5];
        dphiTotal = phaseSpaceMapping[picId * 13 + 6];
        dphiVpara = phaseSpaceMapping[picId * 13 + 7];
        dT = phaseSpaceMapping[picId * 13 + 8];
        bounce = phaseSpaceMapping[picId * 13 + 9];

        // orbit = 0.5 : pad
        // orbit = 1.5 : loss
        // orbit = 2.5 : para
        // orbit = 3.5 : anti
        // orbit = 4.5 : trapped
        // orbit = 5.5 : unknown

        if (orbit < 5)
            continue;

        vec0[2] = 0;
        for (int index = 0; index < 4; index++) {
            vec1[index] = vec0[index];
            vec2[index] = vec0[index];
        }

        li = (vec0[0] - xbeg) / picGridDx;
        lj = (vec0[1] - ybeg) / picGridDy;

        if constexpr (std::is_same_v<picReal, double>) {
            i = __double2int_rd(li);
            j = __double2int_rd(lj);
        } else {
            i = __float2int_rd(li);
            j = __float2int_rd(lj);
        }

        dx = li - i;
        dy = lj - j;

        qId = i * qStride;
        tileId = (j * cellNx + i) * tileStride;
        illegal = 0;

        coes[0] = hx[0] + sx[0] * dx;
        coes[1] = hx[1] + sx[1] * dx;
        FieldGather1d2d<2>(qId, coes, pic1d, q0);

        auto interpRK4 = [&]() {
            for (int index = 0; index < 2; index++)
                coes[index] = (hx[index] + sx[index] * dx);
            FieldGather1d2d<2>(qId, coes, pic1d, q, q_px);

            for (int index = 0; index < 2; index++)
                coes[index + 2] = coes[index];

            for (int index = 0; index < 4; index++)
                coes[index] *= (hy[index] + sy[index] * dy);
            FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py, gcovxy, gcovyy, gcovyz, gcovxy_py,
                               gcovyy_px, gcovyz_px, gcovyz_py);

            m2e = cm * partMass / partChar;
            mu2e = cm * mu / partChar;

            rho = rho0 + vec1[0] * drho;
            bcony = 2 * psitmax * drho * rho / (q * J * B);

            invJ = 1 / J;
            invB = 1 / B;
            invQ = 1 / q;
            invRho = 1 / rho;
            invM2e = 1 / m2e;
            bconyOverJ = bcony * invJ;

            bx = bcony * gcovxy;
            by = bcony * gcovyy;
            bz = bcony * gcovyz;

            cx = bconyOverJ * (gcovyz_py - gcovyz * (J_py * invJ + B_py * invB));
            cy = -bconyOverJ *
                 (gcovyz_px - gcovyz * (J_px * invJ + B_px * invB) + gcovyz * (drho * invRho - q_px * invQ));
            cz = bconyOverJ * (gcovyy_px - gcovyy * (J_px * invJ + B_px * invB) - gcovxy_py +
                               gcovxy * (J_py * invJ + B_py * invB) + gcovyy * (drho * invRho - q_px * invQ));

            dxy = bconyOverJ * gcovyz;
            dxz = bconyOverJ * gcovyy;
            dyz = bconyOverJ * gcovxy;

            dxB = dxy * B_py;
            dyB = -dxy * B_px;
            dzB = dxz * B_px - dyz * B_py;

            Bstarx = cx * m2e * vec1[3];
            Bstary = cy * m2e * vec1[3] + B * bcony;
            Bstarz = cz * m2e * vec1[3];

            Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;
            invBstar = 1 / Bstar;

            ddt[0] = invBstar * (vec1[3] * Bstarx - mu2e * dxB);
            ddt[1] = invBstar * (vec1[3] * Bstary - mu2e * dyB);
            ddt[2] = invBstar * (vec1[3] * Bstarz - mu2e * dzB);
            ddt[3] = -invM2e * invBstar * mu2e * (Bstarx * B_px + Bstary * B_py);
        };

        auto advanceRK4 = [&](picReal bRK, picReal cRK) {
            for (int index = 0; index < 4; index++)
                vec2[index] += ddt[index] * picGridDt * ratioDt * bRK;

            flag = vec0[0] + ddt[0] * picGridDt * ratioDt * cRK;
            if (flag < xbeg || flag >= xend)
                illegal = 1;
            if (!illegal)
                for (int index = 0; index < 4; index++)
                    vec1[index] = vec0[index] + ddt[index] * picGridDt * ratioDt * cRK;
            else
                for (int index = 0; index < 4; index++)
                    vec1[index] = vec0[index];

            li = (vec1[0] - xbeg) / picGridDx;

            if constexpr (std::is_same_v<picReal, double>) {
                i = __double2int_rd(li);
                vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
                lj = (vec1[1] - ybeg) / picGridDy;
                j = __double2int_rd(lj);
            } else {
                i = __float2int_rd(li);
                vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
                lj = (vec1[1] - ybeg) / picGridDy;
                j = __float2int_rd(lj);
            }

            dx = li - i;
            dy = lj - j;
            qId = i * qStride;
            tileId = (j * cellNx + i) * tileStride;
        };

        /*--------------------------------------1st RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 6, picReal(1) / 2);

        /*--------------------------------------2nd RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 3, picReal(1) / 2);

        /*--------------------------------------3rd RK4---------------------------------------*/

        interpRK4();
        advanceRK4(picReal(1) / 3, picReal(1));

        /*--------------------------------------4th RK4---------------------------------------*/

        interpRK4();

        for (int index = 0; index < 4; index++)
            vec2[index] += ddt[index] * picGridDt * ratioDt / 6;

        if (illegal || vec2[0] < xbeg || vec2[0] >= xend) {
            orbit = 1.5;
            for (int index = 0; index < 4; index++)
                vec2[index] = vec0[index];
        }

        frac = 1;

        if (orbit > 5 && bounce < 1 &&
            ((dtheta + vec2[1] - vec0[1]) < -2 * pi || (dtheta + vec2[1] - vec0[1]) > 2 * pi)) {

            if (vec2[3] > 0)
                orbit = 2.5;
            else
                orbit = 3.5;

            if ((dtheta + vec2[1] - vec0[1]) < -2 * pi)
                frac = (-2 * pi - dtheta) / (vec2[1] - vec0[1]);
            else
                frac = (2 * pi - dtheta) / (vec2[1] - vec0[1]);
        }

        if (orbit > 5 && (vec0[3] * vec2[3] < 0)) {
            if (bounce < 1) {
                dT = 0;
                dtheta = 0;
                dphiTotal = 0;
                dphiVpara = 0;
                frac = vec2[3] / (vec2[3] - vec0[3]);
            } else if (bounce > 2) {
                orbit = 4.5;
                frac = 1 - vec2[3] / (vec2[3] - vec0[3]);
            }
            bounce = bounce + 1;
        }

        li = (vec2[0] - xbeg) / picGridDx;
        if constexpr (std::is_same_v<picReal, double>)
            i = __double2int_rd(li);
        else
            i = __float2int_rd(li);
        dx = li - i;
        qId = i * qStride;
        coes[0] = hx[0] + sx[0] * dx;
        coes[1] = hx[1] + sx[1] * dx;
        FieldGather1d2d<2>(qId, coes, pic1d, q);

        dT = dT + picGridDt * ratioDt * frac;
        dtheta = dtheta + (vec2[1] - vec0[1]) * frac;
        dphiTotal = dphiTotal + (vec2[2] + q * vec2[1] - (vec0[2] + q0 * vec0[1])) * frac;
        dphiVpara = dphiVpara + (q + q0) * (vec2[1] - vec0[1]) * frac / 2;

        if constexpr (std::is_same_v<picReal, double>) {
            vec2[1] = vec2[1] - floor((vec2[1] - yori) / yrange) * yrange;
            lj = (vec2[1] - ybeg) / picGridDy;
            j = __double2int_rd(lj);
        } else {
            vec2[1] = vec2[1] - floorf((vec2[1] - yori) / yrange) * yrange;
            lj = (vec2[1] - ybeg) / picGridDy;
            j = __float2int_rd(lj);
        }

        dy = lj - j;
        qId = i * qStride + 28;
        tileId = (j * cellNx + i) * tileStride;

        FieldGather1d2d<2>(qId, coes, pic1d, psip);

        for (int index = 0; index < 2; index++)
            coes[index + 2] = coes[index];
        for (int index = 0; index < 4; index++)
            coes[index] *= (hy[index] + sy[index] * dy);

        FieldGather1d2d<4>(tileId, coes, pic2d, J, B);
        tileId = (j * cellNx + i) * tileStride + 32;
        FieldGather1d2d<4>(tileId, coes, pic2d, gcovyz);

        E = partMass * vec2[3] * vec2[3] / 2 + mu * B;
        Pphi = cm * partMass * vec2[3] * 2 * psitmax * drho * (RHO0 + vec2[0] * drho) * gcovyz / (q * J * B) -
               partChar * psip;

        Lambda = mu / E;

        // orbit x y vp mu dtheta dphiTotal dphiVpara dT bounce E Pphi Lambda

        phaseSpaceMapping[picId * 13 + 0] = orbit;
        phaseSpaceMapping[picId * 13 + 1] = vec2[0];
        phaseSpaceMapping[picId * 13 + 2] = vec2[1];
        phaseSpaceMapping[picId * 13 + 3] = vec2[3];
        phaseSpaceMapping[picId * 13 + 5] = dtheta;
        phaseSpaceMapping[picId * 13 + 6] = dphiTotal;
        phaseSpaceMapping[picId * 13 + 7] = dphiVpara;
        phaseSpaceMapping[picId * 13 + 8] = dT;
        phaseSpaceMapping[picId * 13 + 9] = bounce;
        phaseSpaceMapping[picId * 13 + 10] = E;
        phaseSpaceMapping[picId * 13 + 11] = Pphi;
        phaseSpaceMapping[picId * 13 + 12] = Lambda;
    }
}
