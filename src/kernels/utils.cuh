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

/*----------------------------------Device Helper Functions-----------------------------------*/

template <typename local, typename type>
__device__ void Staggered2C(int offsetx, int offsetz, int& i, int& j, int& k, int& offset2d, int& offset3d,
                            int& lane_id, int& shift_k, type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta,
                            type qtheta_lr[4], type* address, type field_du[4], type field_lr[4], type& field) {

    offset2d = (j + gridGhost) * gridNx + i + offsetx;

    qtheta = d_qtheta[offset2d];
    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    qtheta_lr[0] = (qtheta_lr[0] + qtheta_lr[1]) / 2;
    qtheta_lr[3] = (qtheta_lr[2] + qtheta_lr[3]) / 2;
    qtheta_lr[1] = (qtheta_lr[1] + qtheta) / 2;
    qtheta_lr[2] = (qtheta_lr[2] + qtheta) / 2;

    if constexpr (std::is_same_v<local, trueType>) {

        field_lr[0] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz - 1 * gridNxz];
        field_lr[1] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz + 0 * gridNxz];
        field_lr[2] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz + 1 * gridNxz];
        field_lr[3] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz + 2 * gridNxz];

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
            field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
            field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
            field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }

    } else {

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k +
                                  ((k + offsetz + shift_k - 1) % gridNz + gridNz) % gridNz];
            field_du[1] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k +
                                  ((k + offsetz + shift_k + 0) % gridNz + gridNz) % gridNz];
            field_du[2] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k +
                                  ((k + offsetz + shift_k + 1) % gridNz + gridNz) % gridNz];
            field_du[3] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k +
                                  ((k + offsetz + shift_k + 2) % gridNz + gridNz) % gridNz];

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }
    }

    field = ((field_lr[1] + field_lr[2]) * 9 - (field_lr[0] + field_lr[3])) / 16;

    offset2d = j * gridNx + i;
}

template <typename local, typename type>
__device__ void Collocated2S(int offsetx, int offsetz, int& i, int& j, int& k, int& offset2d, int& offset3d,
                             int& lane_id, int& shift_k, type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta,
                             type qtheta_lr[4], type* address, type field_du[4], type field_lr[4], type& field) {

    offset2d = (j + gridGhost) * gridNx + i + offsetx;

    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 0 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 1 * gridNx];
    qtheta = (qtheta_lr[1] + qtheta_lr[2]) / 2;

    if constexpr (std::is_same_v<local, trueType>) {

        field_lr[0] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz - 2 * gridNxz];
        field_lr[1] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz - 1 * gridNxz];
        field_lr[2] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz + 0 * gridNxz];
        field_lr[3] = address[offset3d + offsetx * gridNz - k + (k + offsetz + gridNz) % gridNz + 1 * gridNxz];

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
            field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
            field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
            field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }

    } else {

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k +
                                  ((k + offsetz + shift_k - 1) % gridNz + gridNz) % gridNz];
            field_du[1] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k +
                                  ((k + offsetz + shift_k + 0) % gridNz + gridNz) % gridNz];
            field_du[2] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k +
                                  ((k + offsetz + shift_k + 1) % gridNz + gridNz) % gridNz];
            field_du[3] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k +
                                  ((k + offsetz + shift_k + 2) % gridNz + gridNz) % gridNz];

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }
    }

    field = ((field_lr[1] + field_lr[2]) * 9 - (field_lr[0] + field_lr[3])) / 16;

    offset2d = j * gridNx + i;
}

template <typename type>
__device__ void PartialX(int offsetz, int& i, int& k, int& offset3d, type* address, type& field, type field_lr[4],
                         type& field_px) {

    offsetz = -k + (k + offsetz + gridNz) % gridNz;

    if (i == 0) {
        field_lr[0] = address[offset3d + offsetz + 1 * gridNz];
        field_lr[1] = address[offset3d + offsetz + 2 * gridNz];
        field_lr[2] = address[offset3d + offsetz + 3 * gridNz];
        field_lr[3] = address[offset3d + offsetz + 4 * gridNz];
        field_px =
            (-25 * field + 48 * field_lr[0] - 36 * field_lr[1] + 16 * field_lr[2] - 3 * field_lr[3]) / (12 * mhdGridDx);
    } else if (i == gridNx - 1) {
        field_lr[0] = address[offset3d + offsetz - 4 * gridNz];
        field_lr[1] = address[offset3d + offsetz - 3 * gridNz];
        field_lr[2] = address[offset3d + offsetz - 2 * gridNz];
        field_lr[3] = address[offset3d + offsetz - 1 * gridNz];
        field_px =
            (3 * field_lr[0] - 16 * field_lr[1] + 36 * field_lr[2] - 48 * field_lr[3] + 25 * field) / (12 * mhdGridDx);
    } else if (i == 1) {
        field_lr[0] = address[offset3d + offsetz - 1 * gridNz];
        field_lr[1] = address[offset3d + offsetz + 1 * gridNz];
        field_lr[2] = address[offset3d + offsetz + 2 * gridNz];
        field_lr[3] = address[offset3d + offsetz + 3 * gridNz];
        field_px =
            (-3 * field_lr[0] - 10 * field + 18 * field_lr[1] - 6 * field_lr[2] + field_lr[3]) / (12 * mhdGridDx);
    } else if (i == gridNx - 2) {
        field_lr[0] = address[offset3d + offsetz - 3 * gridNz];
        field_lr[1] = address[offset3d + offsetz - 2 * gridNz];
        field_lr[2] = address[offset3d + offsetz - 1 * gridNz];
        field_lr[3] = address[offset3d + offsetz + 1 * gridNz];
        field_px =
            (-field_lr[0] + 6 * field_lr[1] - 18 * field_lr[2] + 10 * field + 3 * field_lr[3]) / (12 * mhdGridDx);
    } else {
        field_lr[0] = address[offset3d + offsetz - 2 * gridNz];
        field_lr[1] = address[offset3d + offsetz - 1 * gridNz];
        field_lr[2] = address[offset3d + offsetz + 1 * gridNz];
        field_lr[3] = address[offset3d + offsetz + 2 * gridNz];
        field_px = (field_lr[0] - 8 * field_lr[1] + 8 * field_lr[2] - field_lr[3]) / (12 * mhdGridDx);
    }
}

template <typename type>
__device__ void PartialX2(int& i, int& offset3d, type* address, type& field, type field_lr[4], type& field_px,
                          type& field_px2) {

    if (i == 0) {
        field_lr[0] = address[offset3d + 1 * gridNz];
        field_lr[1] = address[offset3d + 2 * gridNz];
        field_lr[2] = address[offset3d + 3 * gridNz];
        field_lr[3] = address[offset3d + 4 * gridNz];
        field_px =
            (-25 * field + 48 * field_lr[0] - 36 * field_lr[1] + 16 * field_lr[2] - 3 * field_lr[3]) / (12 * mhdGridDx);
        field_px2 = (35 * field - 104 * field_lr[0] + 114 * field_lr[1] - 56 * field_lr[2] + 11 * field_lr[3]) /
                    (12 * mhdGridDx * mhdGridDx);
    } else if (i == gridNx - 1) {
        field_lr[0] = address[offset3d - 4 * gridNz];
        field_lr[1] = address[offset3d - 3 * gridNz];
        field_lr[2] = address[offset3d - 2 * gridNz];
        field_lr[3] = address[offset3d - 1 * gridNz];
        field_px =
            (3 * field_lr[0] - 16 * field_lr[1] + 36 * field_lr[2] - 48 * field_lr[3] + 25 * field) / (12 * mhdGridDx);
        field_px2 = (11 * field_lr[0] - 56 * field_lr[1] + 114 * field_lr[2] - 104 * field_lr[3] + 35 * field) /
                    (12 * mhdGridDx * mhdGridDx);
    } else if (i == 1) {
        field_lr[0] = address[offset3d - 1 * gridNz];
        field_lr[1] = address[offset3d + 1 * gridNz];
        field_lr[2] = address[offset3d + 2 * gridNz];
        field_lr[3] = address[offset3d + 3 * gridNz];
        field_px =
            (-3 * field_lr[0] - 10 * field + 18 * field_lr[1] - 6 * field_lr[2] + field_lr[3]) / (12 * mhdGridDx);
        field_px2 = (11 * field_lr[0] - 20 * field + 6 * field_lr[1] + 4 * field_lr[2] - field_lr[3]) /
                    (12 * mhdGridDx * mhdGridDx);
    } else if (i == gridNx - 2) {
        field_lr[0] = address[offset3d - 3 * gridNz];
        field_lr[1] = address[offset3d - 2 * gridNz];
        field_lr[2] = address[offset3d - 1 * gridNz];
        field_lr[3] = address[offset3d + 1 * gridNz];
        field_px =
            (-field_lr[0] + 6 * field_lr[1] - 18 * field_lr[2] + 10 * field + 3 * field_lr[3]) / (12 * mhdGridDx);
        field_px2 = (-field_lr[0] + 4 * field_lr[1] + 6 * field_lr[2] - 20 * field + 11 * field_lr[3]) /
                    (12 * mhdGridDx * mhdGridDx);
    } else {
        field_lr[0] = address[offset3d - 2 * gridNz];
        field_lr[1] = address[offset3d - 1 * gridNz];
        field_lr[2] = address[offset3d + 1 * gridNz];
        field_lr[3] = address[offset3d + 2 * gridNz];
        field_px = (field_lr[0] - 8 * field_lr[1] + 8 * field_lr[2] - field_lr[3]) / (12 * mhdGridDx);
        field_px2 = (-field_lr[0] + 16 * field_lr[1] - 30 * field + 16 * field_lr[2] - field_lr[3]) /
                    (12 * mhdGridDx * mhdGridDx);
    }
}

template <typename local, typename type>
__device__ void PartialY(int& k, int& offset3d, int& lane_id, int& shift_k, type& shift_lk, type& shift_dk,
                         type& qtheta, type qtheta_lr[4], type* address, type field_du[4], type field_lr[4],
                         type& field_py) {

    if constexpr (std::is_same_v<local, trueType>) {

        field_lr[0] = address[offset3d - 2 * gridNxz];
        field_lr[1] = address[offset3d - 1 * gridNxz];
        field_lr[2] = address[offset3d + 1 * gridNxz];
        field_lr[3] = address[offset3d + 2 * gridNxz];

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
            field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
            field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
            field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }

    } else {

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] = address[offset3d + (index - 2 + index / 2) * gridNxz - k +
                                  ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
            field_du[1] = address[offset3d + (index - 2 + index / 2) * gridNxz - k +
                                  ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
            field_du[2] = address[offset3d + (index - 2 + index / 2) * gridNxz - k +
                                  ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
            field_du[3] = address[offset3d + (index - 2 + index / 2) * gridNxz - k +
                                  ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }
    }

    field_py = (field_lr[0] - 8 * field_lr[1] + 8 * field_lr[2] - field_lr[3]) / (12 * mhdGridDy);
}

template <typename local, typename type>
__device__ void PartialY2(int& k, int& offset3d, int& lane_id, int& shift_k, type& shift_lk, type& shift_dk,
                          type& qtheta, type qtheta_lr[4], type* address, type& field, type field_du[4],
                          type field_lr[4], type& field_py, type& field_py2) {

    PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, address, field_du, field_lr,
                    field_py);

    field_py2 =
        (-field_lr[0] + 16 * field_lr[1] - 30 * field + 16 * field_lr[2] - field_lr[3]) / (12 * mhdGridDy * mhdGridDy);
}

template <typename local, typename type>
__device__ void PartialZ(int& k, int& offset3d, int& lane_id, type* address, type& field, type field_du[4],
                         type& field_pz) {

    if constexpr (std::is_same_v<local, trueType>) {

        field_du[0] = __shfl_sync(0xffffffff, field, lane_id - 2, gridNz);
        field_du[1] = __shfl_sync(0xffffffff, field, lane_id - 1, gridNz);
        field_du[2] = __shfl_sync(0xffffffff, field, lane_id + 1, gridNz);
        field_du[3] = __shfl_sync(0xffffffff, field, lane_id + 2, gridNz);

    } else {

        field_du[0] = address[offset3d - k + (k - 2 + gridNz) % gridNz];
        field_du[1] = address[offset3d - k + (k - 1 + gridNz) % gridNz];
        field_du[2] = address[offset3d - k + (k + 1 + gridNz) % gridNz];
        field_du[3] = address[offset3d - k + (k + 2 + gridNz) % gridNz];
    }

    field_pz = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * mhdGridDz);
}

template <typename local, typename type>
__device__ void PartialZ2(int& k, int& offset3d, int& lane_id, type* address, type& field, type field_du[4],
                          type& field_pz, type& field_pz2) {

    PartialZ<local>(k, offset3d, lane_id, address, field, field_du, field_pz);

    field_pz2 =
        (-field_du[0] + 16 * field_du[1] - 30 * field + 16 * field_du[2] - field_du[3]) / (12 * mhdGridDz * mhdGridDz);
}

template <typename local, typename type>
__device__ void S2CPartialXYZ(int& i, int& j, int& k, int& offset2d, int& offset3d, int& lane_id, int& shift_k,
                              type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta, type qtheta_lr[4],
                              type* address, type field_du[4], type field_lr[4], type& field, type& field_px,
                              type& field_py, type& field_pz) {

    if (i == 0) {

        //+0
        Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px = -25 * field;

        //+1
        Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 48 * field;

        //+2
        Staggered2C<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -36 * field;

        //+3
        Staggered2C<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 16 * field;

        //+4
        Staggered2C<local>(4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -3 * field;

        field_px /= (12 * mhdGridDx);

    } else if (i == gridNx - 1) {

        //+0
        Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px = 25 * field;

        //-4
        Staggered2C<local>(-4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 3 * field;

        //-3
        Staggered2C<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -16 * field;

        //-2
        Staggered2C<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 36 * field;

        //-1
        Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -48 * field;

        field_px /= (12 * mhdGridDx);

    } else if (i == 1) {

        //+0
        Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px = -10 * field;

        //-1
        Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -3 * field;

        //+1
        Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 18 * field;

        //+2
        Staggered2C<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -6 * field;

        //+3
        Staggered2C<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += field;

        field_px /= (12 * mhdGridDx);

    } else if (i == gridNx - 2) {

        //+0
        Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px = 10 * field;

        //-3
        Staggered2C<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -field;

        //-2
        Staggered2C<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 6 * field;

        //-1
        Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -18 * field;

        //+1
        Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 3 * field;

        field_px /= (12 * mhdGridDx);

    } else {

        //-2
        Staggered2C<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px = field;

        //-1
        Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -8 * field;

        //+1
        Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += 8 * field;

        //+2
        Staggered2C<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_px += -field;

        field_px /= (12 * mhdGridDx);
    }

    Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                       qtheta_lr, address, field_du, field_lr, field);

    field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * mhdGridDy);

    if constexpr (std::is_same_v<local, trueType>) {

        PartialZ<local>(k, offset3d, lane_id, address, field, field_du, field_pz);

    } else {

        Staggered2C<local>(0, -2, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_pz = field;

        Staggered2C<local>(0, -1, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_pz += -8 * field;

        Staggered2C<local>(0, 1, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_pz += 8 * field;

        Staggered2C<local>(0, 2, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
        field_pz += -field;

        field_pz /= (12 * mhdGridDz);

        Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                           qtheta_lr, address, field_du, field_lr, field);
    }
}

template <typename local, typename type>
__device__ void C2SPartialXYZ(int& i, int& j, int& k, int& offset2d, int& offset3d, int& lane_id, int& shift_k,
                              type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta, type qtheta_lr[4],
                              type* address, type field_du[4], type field_lr[4], type& field, type& field_px,
                              type& field_py, type& field_pz) {

    if (i == 0) {

        //+0
        Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px = -25 * field;

        //+1
        Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 48 * field;

        //+2
        Collocated2S<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -36 * field;

        //+3
        Collocated2S<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 16 * field;

        //+4
        Collocated2S<local>(4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -3 * field;

        field_px /= (12 * mhdGridDx);

    } else if (i == gridNx - 1) {

        //+0
        Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px = 25 * field;

        //-4
        Collocated2S<local>(-4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 3 * field;

        //-3
        Collocated2S<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -16 * field;

        //-2
        Collocated2S<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 36 * field;

        //-1
        Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -48 * field;

        field_px /= (12 * mhdGridDx);

    } else if (i == 1) {

        //+0
        Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px = -10 * field;

        //-1
        Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -3 * field;

        //+1
        Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 18 * field;

        //+2
        Collocated2S<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -6 * field;

        //+3
        Collocated2S<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += field;

        field_px /= (12 * mhdGridDx);

    } else if (i == gridNx - 2) {

        //+0
        Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px = 10 * field;

        //-3
        Collocated2S<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -field;

        //-2
        Collocated2S<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 6 * field;

        //-1
        Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -18 * field;

        //+1
        Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 3 * field;

        field_px /= (12 * mhdGridDx);

    } else {

        //-2
        Collocated2S<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px = field;

        //-1
        Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -8 * field;

        //+1
        Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += 8 * field;

        //+2
        Collocated2S<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_px += -field;

        field_px /= (12 * mhdGridDx);
    }

    Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                        qtheta_lr, address, field_du, field_lr, field);

    field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * mhdGridDy);

    if constexpr (std::is_same_v<local, trueType>) {

        PartialZ<local>(k, offset3d, lane_id, address, field, field_du, field_pz);

    } else {

        Collocated2S<local>(0, -2, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_pz = field;

        Collocated2S<local>(0, -1, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_pz += -8 * field;

        Collocated2S<local>(0, 1, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_pz += 8 * field;

        Collocated2S<local>(0, 2, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
        field_pz += -field;

        field_pz /= (12 * mhdGridDz);

        Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk, d_qtheta, qtheta,
                            qtheta_lr, address, field_du, field_lr, field);
    }
}

template <typename dirichlet0, typename dirichlet1, typename type, typename... types>
__device__ void Boundary(type* first, types*... second) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i;
    int j = blockIdx.x;
    int k = threadIdx.x;
    int offset3d;
    type field_lr[4];

    /*-------------------------------------Inner Boundary-------------------------------------*/

    i = 0;
    offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    if constexpr (std::is_same_v<dirichlet0, trueType>) {
        first[offset3d] = 0;
    } else {
        field_lr[0] = first[offset3d + 1 * gridNz];
        field_lr[1] = first[offset3d + 2 * gridNz];
        field_lr[2] = first[offset3d + 3 * gridNz];
        field_lr[3] = first[offset3d + 4 * gridNz];
        first[offset3d] = (48 * field_lr[0] - 36 * field_lr[1] + 16 * field_lr[2] - 3 * field_lr[3]) / 25;
    }

    /*-------------------------------------Outer Boundary-------------------------------------*/

    i = gridNx - 1;
    offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

    if constexpr (std::is_same_v<dirichlet1, trueType>) {
        first[offset3d] = 0;
    } else {
        field_lr[0] = first[offset3d - 4 * gridNz];
        field_lr[1] = first[offset3d - 3 * gridNz];
        field_lr[2] = first[offset3d - 2 * gridNz];
        field_lr[3] = first[offset3d - 1 * gridNz];
        first[offset3d] = (-3 * field_lr[0] + 16 * field_lr[1] - 36 * field_lr[2] + 48 * field_lr[3]) / 25;
    }

    if constexpr (sizeof...(second) > 0)
        Boundary<dirichlet0, dirichlet1>(second...);
}

template <typename local, typename type, typename... types>
__device__ void AlignedGhost(type* d_qtheta, type* first, types*... second) {

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

    shift_lk = -qtheta / mhdGridDz;
    if constexpr (std::is_same_v<type, double>)
        shift_k = __double2int_rd(shift_lk);
    else
        shift_k = __float2int_rd(shift_lk);
    shift_dk = shift_lk - shift_k;

    if constexpr (std::is_same_v<local, trueType>) {

        field = first[offset3d + gridNy * gridNxz];

        field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
        field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
        field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
        field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

    } else {

        field_du[0] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
        field_du[1] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
        field_du[2] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
        field_du[3] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];
    }

    field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
            (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
            shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
            shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

    first[offset3d] = field;

    /*--------------------------------------Right Ghost---------------------------------------*/

    shift_lk = qtheta / mhdGridDz;
    if constexpr (std::is_same_v<type, double>)
        shift_k = __double2int_rd(shift_lk);
    else
        shift_k = __float2int_rd(shift_lk);
    shift_dk = shift_lk - shift_k;

    if constexpr (std::is_same_v<local, trueType>) {

        field = first[offset3d + gridGhost * gridNxz];

        field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
        field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
        field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
        field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

    } else {

        field_du[0] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
        field_du[1] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
        field_du[2] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
        field_du[3] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];
    }

    field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
            (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
            shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
            shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

    first[offset3d + (gridNy + gridGhost) * gridNxz] = field;

    if constexpr (sizeof...(second) > 0)
        AlignedGhost<local>(d_qtheta, second...);
}

template <typename local, typename type, typename... types>
__device__ void Staggered2C(type* d_qtheta, type* staggered, type* collocated, types*... fields) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = (j + gridGhost) * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    int shift_k;
    type shift_lk;
    type shift_dk;
    type qtheta;
    type qtheta_lr[4];
    type field_du[4];
    type field_lr[4];

    /*----------------------------------Staggered2Collocated----------------------------------*/

    qtheta = d_qtheta[offset2d];

    qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
    qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
    qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
    qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

    qtheta_lr[0] = (qtheta_lr[0] + qtheta_lr[1]) / 2;
    qtheta_lr[3] = (qtheta_lr[2] + qtheta_lr[3]) / 2;
    qtheta_lr[1] = (qtheta_lr[1] + qtheta) / 2;
    qtheta_lr[2] = (qtheta_lr[2] + qtheta) / 2;

    if constexpr (std::is_same_v<local, trueType>) {

        field_lr[0] = staggered[offset3d - 1 * gridNxz];
        field_lr[1] = staggered[offset3d + 0 * gridNxz];
        field_lr[2] = staggered[offset3d + 1 * gridNxz];
        field_lr[3] = staggered[offset3d + 2 * gridNxz];

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
            field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
            field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
            field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }

    } else {

        for (int index = 0; index < 4; index++) {

            shift_lk = (qtheta_lr[index] - qtheta) / mhdGridDz;
            if constexpr (std::is_same_v<type, double>)
                shift_k = __double2int_rd(shift_lk);
            else
                shift_k = __float2int_rd(shift_lk);
            shift_dk = shift_lk - shift_k;

            field_du[0] =
                staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
            field_du[1] =
                staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
            field_du[2] =
                staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
            field_du[3] =
                staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

            field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
                              (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
                              shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
                              shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];
        }
    }

    collocated[offset3d] = ((field_lr[1] + field_lr[2]) * 9 - (field_lr[0] + field_lr[3])) / 16;

    if constexpr (sizeof...(fields) > 0)
        Staggered2C<local>(d_qtheta, fields...);
}

template <int dir, typename local, typename type, typename... types>
__device__ void Shifted2A(type* d_qtheta, type* shifted, type* aligned, types*... fields) {

    /*-------------------------------------Related Index--------------------------------------*/

    int i = blockIdx.x * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.x + threadIdx.x;
    int offset2d = (j + gridGhost) * gridNx + i;
    int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
    int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

    /*----------------------------------------Shifted-----------------------------------------*/

    int shift_k;
    type shift_lk;
    type shift_dk;
    type qtheta;
    type field;
    type field_du[4];

    /*---------------------------------Shifted2Aligned(dir=0)---------------------------------*/
    /*---------------------------------Aligned2Shifted(dir=1)---------------------------------*/

    qtheta = d_qtheta[offset2d];

    if constexpr (dir == 0)
        shift_lk = qtheta / mhdGridDz;
    else
        shift_lk = -qtheta / mhdGridDz;

    if constexpr (std::is_same_v<type, double>)
        shift_k = __double2int_rd(shift_lk);
    else
        shift_k = __float2int_rd(shift_lk);

    shift_dk = shift_lk - shift_k;

    if constexpr (std::is_same_v<local, trueType>) {

        field = shifted[offset3d];

        field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
        field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
        field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
        field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

    } else {

        field_du[0] = shifted[offset3d - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
        field_du[1] = shifted[offset3d - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
        field_du[2] = shifted[offset3d - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
        field_du[3] = shifted[offset3d - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];
    }

    field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0] +
            (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1] -
            shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2] +
            shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

    aligned[offset3d] = field;

    if constexpr (sizeof...(fields) > 0)
        Shifted2A<dir, local>(d_qtheta, fields...);
}

template <int size, typename type, typename... types>
__device__ void FieldGather1d2d(int& address, type* coes, type* redundant, type& field, types&... fields) {

    field = 0;

    for (int index = 0; index < size; index++)
        field += redundant[address + index] * coes[index];

    if constexpr (sizeof...(fields) > 0) {
        address += size;
        FieldGather1d2d<size>(address, coes, redundant, fields...);
    }
}

template <typename type>
__device__ void FieldGather3d(int& i, int& j, int& k, int& offset, type* coes, type* redundant, type* fields) {

    offset = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[0] * redundant[offset + index];
    }

    offset = (j * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[1] * redundant[offset + index];
    }

    offset = ((j + 1) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[2] * redundant[offset + index];
    }

    offset = ((j + 1) * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[3] * redundant[offset + index];
    }

    offset = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + 1) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[4] * redundant[offset + index];
    }

    offset = (j * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k + 1) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[5] * redundant[offset + index];
    }

    offset = ((j + 1) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + 1) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[6] * redundant[offset + index];
    }

    offset = ((j + 1) * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k + 1) * 8;

    for (int index = 0; index < 8; index++) {
        fields[index] += coes[7] * redundant[offset + index];
    }
}

