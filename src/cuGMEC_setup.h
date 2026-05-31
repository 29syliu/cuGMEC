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
#include "cuGMEC_param.h"

/*-------------------------------I/O Paths------------------------------*/

const std::string initialDir = outputDir + "/initial";
const std::string finalDir = outputDir + "/final";

/*----------------------------MHD Parameter----------------------------*/

const int gridGhost = 2;
const int diagY = ((hostNums * devNums == 1) ? gridNy / 2 : 0);
const int outerLoopMax = totalSteps / sortSteps / ratioDt;
const int innerLoopMax = sortSteps;
const int MHDLoopMax = ratioDt;

const int cellNx = gridNx - 1;
const int cellNy = gridNy + 2 * gridGhost - 1;
const int cellNz = gridNz + 2 * gridGhost - 1;
const int cellNxz = cellNx * cellNz;
const int gridNxz = gridNx * gridNz;
const int gridNyPlusGhost = gridNy + 2 * gridGhost;
const int gridNzPlusGhost = gridNz + 2 * gridGhost;
const int hostNy = gridNy / hostNums;
const int devNy = hostNy / devNums;
const int refinedNy = gridNy * refinedTimes;

const mhdReal NormQE = QE / (B0 * L0 * L0 / MU0 / VA0);
const mhdReal NormMP = MP / (B0 * B0 * L0 * L0 * L0 / MU0 / VA0 / VA0);
const mhdReal mhdGridDx = 1.0 / (gridNx - 1);
const mhdReal mhdGridDy = 2.0 * PI / gridNy;
const mhdReal mhdGridDz = 2.0 * PI / tubes / gridNz;
const mhdReal mhdGridDt = dt;
const picReal picGridDx = 1.0 / (gridNx - 1);
const picReal picGridDy = 2.0 * PI / gridNy;
const picReal picGridDz = 2.0 * PI / tubes / gridNz;
const picReal picGridDt = dt;

/*-----------------------------PIC Parameter-----------------------------*/

const int qStride = 30;
const int tileStride = 72;
const int cellStride = 64;

const int picHost = gridNx * gridNy * gridNz / hostNums * ppcNums;
const int picDev = gridNx * gridNy * gridNz / hostNums / devNums * ppcNums;

const picReal rho0 = RHO0;
const picReal drho = RHO1 - RHO0;
const picReal psitmax = PSITMAX / (B0 * L0 * L0);
const picReal xbeg = 0.0;
const picReal xend = 1.0;
const picReal ybeg = -PI - (gridGhost - 0.5) * picGridDy;
const picReal zbeg = -PI / tubes - (gridGhost - 0.5) * picGridDz;
const picReal yori = -PI;
const picReal zori = -PI / tubes;
const picReal yrange = 2.0 * PI;
const picReal zrange = 2.0 * PI / tubes;
const picReal pi = PI;
const picReal mp = MP;
const picReal mu0 = MU0;
const picReal pitchB0 = B0;
const picReal kev = KEV;
const picReal va = VA0;
const picReal l0 = L0;
const picReal l4 = 1e19 / (L0 * L0 * L0 * L0);
const picReal cm = VA0 / (L0 * (QE * B0 / MP));

__constant__ picReal IonConst;
__constant__ picReal AlphaConst;
__constant__ picReal BeamConst;
__constant__ picReal IonEPphiLambda[6];
__constant__ picReal AlphaEPphiLambda[6];
__constant__ picReal BeamEPphiLambda[6];
__constant__ picReal hx[8] = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
__constant__ picReal sx[8] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
__constant__ picReal hy[8] = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
__constant__ picReal sy[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
__constant__ picReal hz[8] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
__constant__ picReal sz[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

using ifPIC = std::conditional<std::is_same_v<ifIon, trueType> || std::is_same_v<ifAlpha, trueType> ||
                                   std::is_same_v<ifBeam, trueType>,
                               trueType, falseType>::type;

/*-------------------------Thread Block Setting-------------------------*/

using ifLocal = std::conditional<gridNz == 8 || gridNz == 16 || gridNz == 32, trueType, falseType>::type;

const int MRK4BlockDimx = (gridNz == 8) ? 8 : (gridNz == 16) ? 16 : (gridNz == 32) ? 32 : 16;
const int MRK4BlockDimy = (devNy >= 2) ? 2 : 1;
const int MRK4BlockDimz = 4;
const int MRK4GridDimx = gridNx / MRK4BlockDimz;
const int MRK4GridDimy = devNy / MRK4BlockDimy;
const int MRK4GridDimz = gridNz / MRK4BlockDimx;

const int GhostBlockDimx = (gridNz == 8) ? 8 : (gridNz == 16) ? 16 : (gridNz == 32) ? 32 : 16;
const int GhostBlockDimy = 2;
const int GhostBlockDimz = ((gridNz == 32) ? 4 : 8);
const int GhostGridDimx = gridNx / GhostBlockDimz;
const int GhostGridDimy = gridGhost / GhostBlockDimy;
const int GhostGridDimz = gridNz / GhostBlockDimx;

const int M2PBlockDimx = (gridNz == 8) ? 8 : (gridNz == 16) ? 16 : (gridNz == 32) ? 32 : 16;
const int M2PBlockDimy = ((gridNz == 32) ? 2 : 4);
const int M2PBlockDimz = 4;
const int M2PGridDimx = gridNx / M2PBlockDimz;
const int M2PGridDimy = gridNy / M2PBlockDimy;
const int M2PGridDimz = gridNz / M2PBlockDimx;

const int NMBlockDimx = 256;
const int LocalNMGridDimx = devNy * gridNx / NMBlockDimx;
const int GhostNMGridDimx = gridGhost * gridNx / NMBlockDimx;
const int RefinedNMGridDimx = refinedNy * gridNx / NMBlockDimx;

const int pptNums = 32;
const int PICBlockDimx = 128;
const int PICGridDimx = picDev / pptNums / PICBlockDimx;

const int nFFTBatchSize = devNy * gridNx;
const int nFFTTimeSize = gridNz;
const int nFFTFreqSize = gridNz / 2 + 1;

const int mFFTTimeSize = refinedNy;
