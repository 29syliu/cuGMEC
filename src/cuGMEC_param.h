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
#include "cuGMEC_const.h"

/*-----------------------Data Type and File Address-----------------------*/

using mhdReal = double;
using picReal = float;

const std::string inputDir = "/home-ssd/Users/nsgm_pjj/ITER/0.08-0.90";
const std::string outputDir = "/home-ssd/Users/nsgm_pjj/ITER/0.08-0.90/6_1-6_400";
const std::string MHDCollocated = "MHDCollocated_256_32.bin";
const std::string MHDStaggered = "MHDStaggered_256_32.bin";
const std::string IonPhaseSpaceMapping = "IonPhaseSpaceMapping.bin";
const std::string BeamPhaseSpaceMapping = "BeamPhaseSpaceMapping.bin";
const std::string AlphaPhaseSpaceMapping = "AlphaPhaseSpaceMapping.bin";

/*--------------------------Normalization Setting-------------------------*/

const double B0 = 4.921751144619735;
const double L0 = 6.595629295925759;
const double VA0 = 8.864164667700194e+06;
const double RHO0 = 0.08;
const double RHO1 = 0.90;
const double PSITMAX = 18.868213504765762;

/*------------------------------Scale Setting-----------------------------*/

const int hostNums = 4;
const int devNums = 4;
const int gridNx = 256;
const int gridNy = 32;
const int gridNz = 96;
const int ppcNums = 256;

/*-------------------------------MHD Setting------------------------------*/

const int tubes = 6;
const int leftN = 1;
const int rightN = 6;

const int perturbLeftN = 1;
const int perturbRightN = 6;
const int perturbRadialIndex = 85;
const mhdReal perturbWidth = 0.04;
const mhdReal perturbAmplitude = 2.5e-9;

using ifFilterN_Phi = trueType;
using ifFilterN_A = trueType;
using ifFilterN_dNe = trueType;
using ifFilterN_dTe = trueType;
using ifFilterN_dP = trueType;
constexpr std::array<int, 0> removeN_Phi = {};
constexpr std::array<int, 0> removeN_A = {};
constexpr std::array<int, 0> removeN_dNe = {};
constexpr std::array<int, 0> removeN_dTe = {};
constexpr std::array<int, 0> removeN_dP = {};
constexpr std::array<std::tuple<int, int, int>, 1> selectNM_Phi = {{{0, 0, 0}}};
constexpr std::array<std::tuple<int, int, int>, 1> selectNM_A = {{{0, 0, 1}}};
constexpr std::array<std::tuple<int, int, int>, 1> selectNM_dNe = {{{0, 0, 1}}};
constexpr std::array<std::tuple<int, int, int>, 1> selectNM_dTe = {{{0, 0, 1}}};
constexpr std::array<std::tuple<int, int, int>, 1> selectNM_dP = {{{0, 0, 1}}};

using ifStaggered = falseType;
using ifNonlinearMHD = trueType;
using ifEparallel = trueType;
using ifFLRMHD = falseType;
using ifMaxwellStress = trueType;
using ifReynoldsStress = trueType;
const mhdReal MaxwellStressCoef = 1.0;
const mhdReal ReynoldsStressCoef = 1.0;

using ifNablaPerp2Phi = trueType;
using ifNablaPara4Phi = trueType;
const mhdReal perp2Phi = 2.0e-7;
const mhdReal para4Phi = 1.0e-7;

using ifNablaPerp2A = trueType;
using ifNablaPara4A = trueType;
const mhdReal perp2A = 2.0e-7;
const mhdReal para4A = 1.0e-7;

using ifNablaPerp2dNe = trueType;
using ifNablaPara4dNe = trueType;
const mhdReal perp2dNe = 2.0e-7;
const mhdReal para4dNe = 1.0e-7;

using ifNablaPerp2dTe = trueType;
using ifNablaPara4dTe = trueType;
const mhdReal perp2dTe = 2.0e-7;
const mhdReal para4dTe = 1.0e-7;

using ifNablaPerp2dPi = trueType;
using ifNablaPara4dPi = falseType;
const mhdReal perp2dPi = 1.0e-6;
const mhdReal para4dPi = 0.0;

using ifNablaPerp2dPa = trueType;
using ifNablaPara4dPa = falseType;
const mhdReal perp2dPa = 1.0e-6;
const mhdReal para4dPa = 0.0;

using ifNablaPerp2dPb = trueType;
using ifNablaPara4dPb = falseType;
const mhdReal perp2dPb = 1.0e-6;
const mhdReal para4dPb = 0.0;

using ifConvolveAligned = trueType;
const int convolveTimes = 1;
const mhdReal convolveSigmaMax = 1.5;
const mhdReal convolveSigmaMin = 0.05;
const mhdReal convolveTSwitch = 200.0;
const mhdReal convolveDtSwitch = 25.0;

/*-------------------------------PIC Setting------------------------------*/

using ifFLRPIC = trueType;
using ifNonlinearPIC = trueType;
const int gyroNums = 4;

using ifIon = trueType;
const disType IonType = Maxwell;
const spaceType IonSpace = spaceReal;
const velocityType IonVelocity = velocityUniform;
const picReal IonMass = 2.5;
const picReal IonChar = 1.0;
const picReal IonBeta = 0.037793721898356;
const picReal IonVmin = 0.0135;
const picReal IonVmax = 0.54;
const picReal IonVb = 0.0;
const picReal IonDeltaV = 0.0;
const picReal IonLambda0 = 0.0;
const picReal IonDeltaLambda2 = 0.0;

using ifAlpha = trueType;
const disType AlphaType = Slowing0;
const spaceType AlphaSpace = spaceReal;
const velocityType AlphaVelocity = velocityUniform;
const picReal AlphaMass = 4.0;
const picReal AlphaChar = 2.0;
const picReal AlphaBeta = 0.018554328058860;
const picReal AlphaVmin = 0.0733;
const picReal AlphaVmax = 1.466;
const picReal AlphaVb = 0.0;
const picReal AlphaDeltaV = 0.0;
const picReal AlphaLambda0 = 0.0;
const picReal AlphaDeltaLambda2 = 0.0;

using ifBeam = trueType;
const disType BeamType = Slowing2;
const spaceType BeamSpace = spaceReal;
const velocityType BeamVelocity = velocityUniform;
const picReal BeamMass = 2.0;
const picReal BeamChar = 1.0;
const picReal BeamBeta = 0.010583527513998;
const picReal BeamVmin = 0.0552;
const picReal BeamVmax = 1.104;
const picReal BeamVb = 0.0;
const picReal BeamDeltaV = 0.0;
const picReal BeamLambda0 = 0.4;
const picReal BeamDeltaLambda2 = 1.0 / (4.5 * 4.5);

/*-------------------------------Run Setting------------------------------*/

using ifDiagAmplitude = trueType;
using ifDiagFrequency = trueType;
using ifDiagEparallel = trueType;
using ifDiagDensity = trueType;
using ifDiagDiffusivity = trueType;
using ifDiagZFDrive = falseType;
using ifCheckNAN = trueType;

using ifOutputPhi = trueType;
using ifOutputA = falseType;
using ifOutputdNe = falseType;
using ifOutputdTe = falseType;
using ifOutputdPi = falseType;
using ifOutputdPa = falseType;
using ifOutputdPb = falseType;

const int gridE = 96;
const int gridPphi = 128;
const int gridLambda = 48;
const int ppcPhase = 2048;
using ifOutputPhaseSpaceJacobian = trueType;
using ifOutputPhaseSpaceFrequency = falseType;
using ifOutputPhaseSpaceF0 = falseType;
using ifOutputPhaseSpaceDeltaF = falseType;
using ifOutputPhaseSpacePower = falseType;

const int gridVpara = 128;
const int gridVperp = 64;
const int ppcPitch = 2048;
using ifOutputPitchSpaceJacobian = trueType;
using ifOutputPitchSpaceF0 = falseType;
using ifOutputPitchSpaceDeltaF = falseType;
using ifOutputPitchSpacePower = falseType;

using ifContinue = falseType;
const int continueSteps = 0;
const double dt = 0.02;
const int totalSteps = 20000;
const int ratioDt = 1;
const int sortSteps = 25;
const int diagSteps = 1;
const int outputSteps = 2500;
