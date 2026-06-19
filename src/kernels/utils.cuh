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
__device__ void FieldGather1d2d(int& address, type* coes, type* pic1d2d, type& field, types&... fields) {

    field = 0;

    for (int index = 0; index < size; index++)
        field += pic1d2d[address + index] * coes[index];

    if constexpr (sizeof...(fields) > 0) {
        address += size;
        FieldGather1d2d<size>(address, coes, pic1d2d, fields...);
    }
}

template <typename type>
__device__ void FieldGather3d(int& i, int& j, int& k, int& offset, type* coes, type* pic3d, type* fields) {

    const int gridPointCount = gridNyPlusGhost * gridNx * gridNzPlusGhost;
    const int xStride = gridNzPlusGhost;
    const int yStride = gridNx * gridNzPlusGhost;
    const int fieldGroup0Base = 0 * gridPointCount;
    const int fieldGroup1Base = 1 * gridPointCount;
    const int fieldGroup2Base = 2 * gridPointCount;
    const int fieldGroup3Base = 3 * gridPointCount;
    const int cornerX0Y0Z0Id = j * yStride + i * xStride + k;
    const int cornerX1Y0Z0Id = cornerX0Y0Z0Id + xStride;
    const int cornerX0Y1Z0Id = cornerX0Y0Z0Id + yStride;
    const int cornerX1Y1Z0Id = cornerX0Y1Z0Id + xStride;

    offset = (cornerX1Y1Z0Id + 1) * 8;

    if constexpr (std::is_same_v<type, float>) {
        const float4* pic3dFloat4 = reinterpret_cast<const float4*>(pic3d);
        float4 fieldGroup0 = make_float4(fields[0], fields[1], fields[2], fields[3]);
        float4 fieldGroup1 = make_float4(fields[4], fields[5], fields[6], fields[7]);
        float4 cornerFieldGroup0;
        float4 cornerFieldGroup1;

        cornerFieldGroup0 = pic3dFloat4[cornerX0Y0Z0Id];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX0Y0Z0Id];
        fieldGroup0.x += coes[0] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[0] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[0] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[0] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[0] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[0] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[0] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[0] * cornerFieldGroup1.w;

        cornerFieldGroup0 = pic3dFloat4[cornerX1Y0Z0Id];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX1Y0Z0Id];
        fieldGroup0.x += coes[1] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[1] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[1] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[1] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[1] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[1] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[1] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[1] * cornerFieldGroup1.w;

        cornerFieldGroup0 = pic3dFloat4[cornerX0Y1Z0Id];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX0Y1Z0Id];
        fieldGroup0.x += coes[2] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[2] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[2] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[2] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[2] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[2] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[2] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[2] * cornerFieldGroup1.w;

        cornerFieldGroup0 = pic3dFloat4[cornerX1Y1Z0Id];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX1Y1Z0Id];
        fieldGroup0.x += coes[3] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[3] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[3] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[3] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[3] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[3] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[3] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[3] * cornerFieldGroup1.w;

        cornerFieldGroup0 = pic3dFloat4[cornerX0Y0Z0Id + 1];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX0Y0Z0Id + 1];
        fieldGroup0.x += coes[4] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[4] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[4] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[4] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[4] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[4] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[4] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[4] * cornerFieldGroup1.w;

        cornerFieldGroup0 = pic3dFloat4[cornerX1Y0Z0Id + 1];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX1Y0Z0Id + 1];
        fieldGroup0.x += coes[5] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[5] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[5] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[5] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[5] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[5] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[5] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[5] * cornerFieldGroup1.w;

        cornerFieldGroup0 = pic3dFloat4[cornerX0Y1Z0Id + 1];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX0Y1Z0Id + 1];
        fieldGroup0.x += coes[6] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[6] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[6] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[6] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[6] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[6] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[6] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[6] * cornerFieldGroup1.w;

        cornerFieldGroup0 = pic3dFloat4[cornerX1Y1Z0Id + 1];
        cornerFieldGroup1 = pic3dFloat4[fieldGroup1Base + cornerX1Y1Z0Id + 1];
        fieldGroup0.x += coes[7] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[7] * cornerFieldGroup0.y;
        fieldGroup0.z += coes[7] * cornerFieldGroup0.z;
        fieldGroup0.w += coes[7] * cornerFieldGroup0.w;
        fieldGroup1.x += coes[7] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[7] * cornerFieldGroup1.y;
        fieldGroup1.z += coes[7] * cornerFieldGroup1.z;
        fieldGroup1.w += coes[7] * cornerFieldGroup1.w;

        fields[0] = fieldGroup0.x;
        fields[1] = fieldGroup0.y;
        fields[2] = fieldGroup0.z;
        fields[3] = fieldGroup0.w;
        fields[4] = fieldGroup1.x;
        fields[5] = fieldGroup1.y;
        fields[6] = fieldGroup1.z;
        fields[7] = fieldGroup1.w;
    } else {
        const double2* pic3dDouble2 = reinterpret_cast<const double2*>(pic3d);
        double2 fieldGroup0 = make_double2(fields[0], fields[1]);
        double2 fieldGroup1 = make_double2(fields[2], fields[3]);
        double2 fieldGroup2 = make_double2(fields[4], fields[5]);
        double2 fieldGroup3 = make_double2(fields[6], fields[7]);
        double2 cornerFieldGroup0;
        double2 cornerFieldGroup1;
        double2 cornerFieldGroup2;
        double2 cornerFieldGroup3;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX0Y0Z0Id];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX0Y0Z0Id];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX0Y0Z0Id];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX0Y0Z0Id];
        fieldGroup0.x += coes[0] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[0] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[0] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[0] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[0] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[0] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[0] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[0] * cornerFieldGroup3.y;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX1Y0Z0Id];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX1Y0Z0Id];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX1Y0Z0Id];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX1Y0Z0Id];
        fieldGroup0.x += coes[1] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[1] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[1] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[1] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[1] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[1] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[1] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[1] * cornerFieldGroup3.y;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX0Y1Z0Id];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX0Y1Z0Id];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX0Y1Z0Id];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX0Y1Z0Id];
        fieldGroup0.x += coes[2] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[2] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[2] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[2] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[2] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[2] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[2] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[2] * cornerFieldGroup3.y;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX1Y1Z0Id];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX1Y1Z0Id];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX1Y1Z0Id];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX1Y1Z0Id];
        fieldGroup0.x += coes[3] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[3] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[3] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[3] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[3] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[3] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[3] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[3] * cornerFieldGroup3.y;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX0Y0Z0Id + 1];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX0Y0Z0Id + 1];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX0Y0Z0Id + 1];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX0Y0Z0Id + 1];
        fieldGroup0.x += coes[4] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[4] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[4] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[4] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[4] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[4] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[4] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[4] * cornerFieldGroup3.y;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX1Y0Z0Id + 1];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX1Y0Z0Id + 1];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX1Y0Z0Id + 1];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX1Y0Z0Id + 1];
        fieldGroup0.x += coes[5] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[5] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[5] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[5] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[5] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[5] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[5] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[5] * cornerFieldGroup3.y;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX0Y1Z0Id + 1];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX0Y1Z0Id + 1];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX0Y1Z0Id + 1];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX0Y1Z0Id + 1];
        fieldGroup0.x += coes[6] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[6] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[6] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[6] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[6] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[6] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[6] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[6] * cornerFieldGroup3.y;

        cornerFieldGroup0 = pic3dDouble2[fieldGroup0Base + cornerX1Y1Z0Id + 1];
        cornerFieldGroup1 = pic3dDouble2[fieldGroup1Base + cornerX1Y1Z0Id + 1];
        cornerFieldGroup2 = pic3dDouble2[fieldGroup2Base + cornerX1Y1Z0Id + 1];
        cornerFieldGroup3 = pic3dDouble2[fieldGroup3Base + cornerX1Y1Z0Id + 1];
        fieldGroup0.x += coes[7] * cornerFieldGroup0.x;
        fieldGroup0.y += coes[7] * cornerFieldGroup0.y;
        fieldGroup1.x += coes[7] * cornerFieldGroup1.x;
        fieldGroup1.y += coes[7] * cornerFieldGroup1.y;
        fieldGroup2.x += coes[7] * cornerFieldGroup2.x;
        fieldGroup2.y += coes[7] * cornerFieldGroup2.y;
        fieldGroup3.x += coes[7] * cornerFieldGroup3.x;
        fieldGroup3.y += coes[7] * cornerFieldGroup3.y;

        fields[0] = fieldGroup0.x;
        fields[1] = fieldGroup0.y;
        fields[2] = fieldGroup1.x;
        fields[3] = fieldGroup1.y;
        fields[4] = fieldGroup2.x;
        fields[5] = fieldGroup2.y;
        fields[6] = fieldGroup3.x;
        fields[7] = fieldGroup3.y;
    }
}
