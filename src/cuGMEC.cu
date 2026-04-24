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

#include "cuGMEC.h"
#include <nlohmann/json.hpp>
#include <stdexcept>

using json = nlohmann::json;

static uint64_t getHostHash(const char* string) {
	// Based on DJB2a, result = result * 33 ^ char
	uint64_t result = 5381;
	for (int c = 0; string[c] != '\0'; c++) {
		result = ((result << 5) + result) ^ string[c];
	}
	return result;
}

static void getHostName(char* hostname, int maxlen) {
	gethostname(hostname, maxlen);
	for (int i = 0; i < maxlen; i++) {
		if (hostname[i] == '.') {
			hostname[i] = '\0';
			return;
		}
	}
}

/*--------------------------------Constant--------------------------------*/

const double QE = 1.6021766208e-19;
const double MP = 1.672621637e-27;
const double PI = 3.1415926535897932;
const double MU0 = 4.0 * PI * 1.0e-7;
const double KEV = 1000.0 * QE;
using trueType = std::integral_constant<bool, true>;
using falseType = std::integral_constant<bool, false>;
enum picType { Ion, Alpha, Beam };
enum disType { Maxwell, Slowing0, Slowing1, Slowing2, Slowing3 };
enum markerType { physicalReal, physicalUniform, numericalUniform };
enum matrixType { Laplacian, Resistive, Perp2Phi, Perp2dNe, Perp2dTe, Perp2dPi, Perp2dPa, Perp2dPb };

/*-----------------------Data Type and File Address-----------------------*/

using dataType = double;

std::string h_MHDCollocated = "../MHDCollocated_512_64.bin";
std::string h_MHDStaggered = "../MHDStaggered_512_64.bin";
std::string h_MHDPerturbation = "../MHDPerturbation_512_64_16_30_1.bin";

/*-------------------------Normalization Setting-------------------------*/

dataType h_B0 = 4.921751144619735;
dataType h_L0 = 6.595629295925759;
dataType h_VA0 = 8.864164667700194e+06;
dataType h_RHO0 = 0.08;
dataType h_RHO1 = 0.95;
dataType h_PSITMAX = 18.868213504765762;
dataType h_OMEGACI = QE * h_B0 / MP;

__constant__ dataType B0;
__constant__ dataType L0;
__constant__ dataType VA0;
__constant__ dataType RHO0;
__constant__ dataType RHO1;
__constant__ dataType PSITMAX;
__constant__ dataType OMEGACI;

/*------------------------------Scale Setting------------------------------*/

const int hostNums = 2;
const int devNums = 4;
const int gridNx = 512;
const int gridNy = 64;
const int gridNz = 16;
const int gridGhost = 2;
int h_ppcNums = 128;
__constant__ int ppcNums;

/*-----------------------------MHD Setting-----------------------------*/

int h_tubes = 30;
__constant__ int tubes;
int h_leftN = 1;
int h_rightN = 1;
__constant__ int leftN;
__constant__ int rightN;
using ifFilterN_Phi = trueType;
using ifFilterN_A = trueType;
using ifFilterN_dNe = trueType;
using ifFilterN_dTe = trueType;
using ifFilterN_dP = trueType;
constexpr std::array<int, 0> removeN_Phi = { };
constexpr std::array<int, 0> removeN_A = { };
constexpr std::array<int, 0> removeN_dNe = { };
constexpr std::array<int, 0> removeN_dTe = { };
constexpr std::array<int, 0> removeN_dP = { };
constexpr std::array<std::tuple<int, int, int>, 0> selectNM_Phi = { };
constexpr std::array<std::tuple<int, int, int>, 0> selectNM_A = { };
constexpr std::array<std::tuple<int, int, int>, 0> selectNM_dNe = { };
constexpr std::array<std::tuple<int, int, int>, 0> selectNM_dTe = { };
constexpr std::array<std::tuple<int, int, int>, 0> selectNM_dP = { };

using ifStaggered = falseType;
using ifNonlinear = falseType;
using ifEparallel = trueType;
using ifFLRMHD = falseType;

using ifNablaPerp2Phi = falseType;
using ifNablaPara4Phi = falseType;
dataType h_perp2Phi = 0.0;
dataType h_para4Phi = 0.0;

using ifNablaPerp2A = falseType;
using ifNablaPara4A = falseType;
dataType h_perp2A = 0.0;
dataType h_para4A = 0.0;

using ifNablaPerp2dNe = falseType;
using ifNablaPara4dNe = falseType;
dataType h_perp2dNe = 0.0;
dataType h_para4dNe = 0.0;

using ifNablaPerp2dTe = falseType;
using ifNablaPara4dTe = falseType;
dataType h_perp2dTe = 0.0;
dataType h_para4dTe = 0.0;

using ifNablaPerp2dPi = falseType;
using ifNablaPara4dPi = falseType;
dataType h_perp2dPi = 0.0;
dataType h_para4dPi = 0.0;

using ifNablaPerp2dPa = falseType;
using ifNablaPara4dPa = falseType;
dataType h_perp2dPa = 0.0;
dataType h_para4dPa = 0.0;

using ifNablaPerp2dPb = falseType;
using ifNablaPara4dPb = falseType;
dataType h_perp2dPb = 0.0;
dataType h_para4dPb = 0.0;

/*------------------------------PIC Setting------------------------------*/

using randomState = curandStateXORWOW_t;
const unsigned int randMax = 1 << 24;

using ifFLRPIC = trueType;
const int gyroNums = 4;

using ifIon = trueType;
const disType IonType = Maxwell;
const markerType IonMarker = physicalReal;
using ifIonSlowing = falseType;
dataType h_IonMass = 2.5;
dataType h_IonChar = 1.0;
dataType h_IonBeta = 0.037793721898356;
dataType h_IonVmin = 0.0135;
dataType h_IonVmax = 0.54;
dataType h_IonVb = 0.0;
dataType h_IonDeltaV = 0.0;
dataType h_IonLambda0 = 0.0;
dataType h_IonDeltaLambda2 = 0.0;
__constant__ dataType IonMass;
__constant__ dataType IonChar;
__constant__ dataType IonBeta;
__constant__ dataType IonVmin;
__constant__ dataType IonVmax;
__constant__ dataType IonVb;
__constant__ dataType IonDeltaV;
__constant__ dataType IonLambda0;
__constant__ dataType IonDeltaLambda2;
const dataType IonDragRate = 0.0;

using ifAlpha = trueType;
const disType AlphaType = Slowing0;
const markerType AlphaMarker = physicalReal;
using ifAlphaSlowing = falseType;
dataType h_AlphaMass = 4.0;
dataType h_AlphaChar = 2.0;
dataType h_AlphaBeta = 0.018554328058860;
dataType h_AlphaVmin = 0.0733;
dataType h_AlphaVmax = 1.466;
dataType h_AlphaVb = 0.0;
dataType h_AlphaDeltaV = 0.0;
dataType h_AlphaLambda0 = 0.0;
dataType h_AlphaDeltaLambda2 = 0.0;
__constant__ dataType AlphaMass;
__constant__ dataType AlphaChar;
__constant__ dataType AlphaBeta;
__constant__ dataType AlphaVmin;
__constant__ dataType AlphaVmax;
__constant__ dataType AlphaVb;
__constant__ dataType AlphaDeltaV;
__constant__ dataType AlphaLambda0;
__constant__ dataType AlphaDeltaLambda2;
const dataType AlphaDragRate = 0.0;

using ifBeam = trueType;
const disType BeamType = Slowing2;
const markerType BeamMarker = physicalReal;
using ifBeamSlowing = falseType;
dataType h_BeamMass = 2.0;
dataType h_BeamChar = 1.0;
dataType h_BeamBeta = 0.010583527513998;
dataType h_BeamVmin = 0.0552;
dataType h_BeamVmax = 1.104;
dataType h_BeamVb = 0.0;
dataType h_BeamDeltaV = 0.0;
dataType h_BeamLambda0 = 0.4;
dataType h_BeamDeltaLambda2 = 1.0 / (4.5 * 4.5);
__constant__ dataType BeamMass;
__constant__ dataType BeamChar;
__constant__ dataType BeamBeta;
__constant__ dataType BeamVmin;
__constant__ dataType BeamVmax;
__constant__ dataType BeamVb;
__constant__ dataType BeamDeltaV;
__constant__ dataType BeamLambda0;
__constant__ dataType BeamDeltaLambda2;
const dataType BeamDragRate = 0.0;

/*------------------------------Run Setting------------------------------*/

using ifContinue = falseType;
using ifDiagAmplitude = trueType;
using ifDiagFrequency = trueType;
using ifDiagEparallel = trueType;
using ifDiagDensity = falseType;
using ifDiagDiffusivity = falseType;

using ifOutputw = falseType;
using ifOutputA = falseType;
using ifOutputdNe = falseType;
using ifOutputdTe = falseType;
using ifOutputPhi = falseType;
using ifOutputdPi = falseType;
using ifOutputdPa = falseType;
using ifOutputdPb = falseType;

//const int gridE = 64;
//const int gridPphi = 64;
//const int gridLambda = 64;
//const int ppcPhase = 4096;
//using ifOutputPhaceSpaceF0 = falseType;
//using ifOutputPhaceSpaceOribit = falseType;
//using ifOutputPhaceSpaceDeltaF = falseType;
//using ifOutputPhaceSpaceJacobian = falseType;
//const std::string IonPhaseSpaceMapping = "../MHDCollocated_512_64.bin";
//const std::string AlphaPhaseSpaceMapping = "../MHDCollocated_512_64.bin";
//const std::string BeamPhaseSpaceMapping = "../MHDCollocated_512_64.bin";

dataType h_dt = 0.02;
int h_continueSteps = 0;
const int ratioDt = 1;
int h_totalSteps = 4000;
int h_sortSteps = 25;
int h_diagSteps = 1;
int h_outputSteps = 5000;
int h_diagLeftX = 0;
int h_diagRightX = gridNx - 1;

int h_diagY = ((hostNums * devNums == 1) ? gridNy / 2 : 0);
__constant__ int diagLeftX;
__constant__ int diagRightX;
__constant__ int diagY;
int h_outerLoopMax = h_totalSteps / h_sortSteps / ratioDt;
int h_innerLoopMax = h_sortSteps;
int h_MHDLoopMax = ratioDt;

/*----------------------------MHD Parameter----------------------------*/

const int cellNx = gridNx - 1;
const int cellNy = gridNy + 2 * gridGhost - 1;
const int cellNz = gridNz + 2 * gridGhost - 1;
const int cellNxz = cellNx * cellNz;
const int gridNxz = gridNx * gridNz;
const int gridNyPlusGhost = gridNy + 2 * gridGhost;
const int gridNzPlusGhost = gridNz + 2 * gridGhost;
const int hostNy = gridNy / hostNums;
const int devNy = hostNy / devNums;

__constant__ dataType NormQE;
__constant__ dataType NormMP;
const dataType gridDx = 1.0 / (gridNx - 1);
const dataType gridDy = 2.0 * PI / gridNy;
__constant__ dataType gridDz;
__constant__ dataType gridDt;

/*-----------------------------PIC Parameter-----------------------------*/

int h_picHost = gridNx * gridNy * gridNz / hostNums * h_ppcNums;
int h_picDev = gridNx * gridNy * gridNz / hostNums / devNums * h_ppcNums;
__constant__ int picHost;
__constant__ int picDev;

__constant__ dataType rho0;
__constant__ dataType drho;
__constant__ dataType psitmax;
const dataType xbeg = 0.0;
const dataType xend = 1.0;
const dataType ybeg = -PI - (gridGhost - 0.5) * gridDy;
__constant__ dataType zbeg;
const dataType yori = -PI;
__constant__ dataType zori;
const dataType yrange = 2.0 * PI;
__constant__ dataType zrange;
const dataType pi = PI;
const dataType mp = MP;
const dataType mu0 = MU0;
__constant__ dataType pitchB0;
const dataType kev = KEV;
__constant__ dataType va;
__constant__ dataType l3;
__constant__ dataType cm;

__constant__ dataType perp2Phi;
__constant__ dataType para4Phi;
__constant__ dataType perp2A;
__constant__ dataType para4A;
__constant__ dataType perp2dNe;
__constant__ dataType para4dNe;
__constant__ dataType perp2dTe;
__constant__ dataType para4dTe;
__constant__ dataType perp2dPi;
__constant__ dataType para4dPi;
__constant__ dataType perp2dPa;
__constant__ dataType para4dPa;
__constant__ dataType perp2dPb;
__constant__ dataType para4dPb;

static void loadConfigFromJson(const std::string& configPath) {
	std::ifstream ifs(configPath);
	if (!ifs.is_open()) {
		return;
	}

	json j;
	ifs >> j;

	h_MHDCollocated = j.value("MHDCollocated", h_MHDCollocated);
	h_MHDStaggered = j.value("MHDStaggered", h_MHDStaggered);
	h_MHDPerturbation = j.value("MHDPerturbation", h_MHDPerturbation);

	h_B0 = j.value("B0", h_B0);
	h_L0 = j.value("L0", h_L0);
	h_VA0 = j.value("VA0", h_VA0);
	h_dt = j.value("dt", h_dt);
	h_continueSteps = j.value("continueSteps", h_continueSteps);
	h_totalSteps = j.value("totalSteps", h_totalSteps);
	h_sortSteps = j.value("sortSteps", h_sortSteps);
	h_diagSteps = j.value("diagSteps", h_diagSteps);
	h_outputSteps = j.value("outputSteps", h_outputSteps);
	int cfgRatioDt = j.value("ratioDt", ratioDt);
	if (cfgRatioDt != ratioDt) {
		throw std::runtime_error("ratioDt is a compile-time parameter and must match compiled values.");
	}
	if (h_sortSteps <= 0 || h_diagSteps <= 0) {
		throw std::runtime_error("sortSteps and diagSteps must be positive.");
	}
	h_outerLoopMax = h_totalSteps / h_sortSteps / ratioDt;
	h_innerLoopMax = h_sortSteps;
	h_MHDLoopMax = ratioDt;
	h_RHO0 = j.value("RHO0", h_RHO0);
	h_RHO1 = j.value("RHO1", h_RHO1);
	h_PSITMAX = j.value("PSITMAX", h_PSITMAX);
	h_tubes = j.value("tubes", h_tubes);
	int cfgGridNx = j.value("gridNx", gridNx);
	int cfgGridNy = j.value("gridNy", gridNy);
	int cfgGridNz = j.value("gridNz", gridNz);
	if (cfgGridNx != gridNx || cfgGridNy != gridNy || cfgGridNz != gridNz) {
		throw std::runtime_error("gridNx/gridNy/gridNz are compile-time parameters and must match compiled values.");
	}
	h_ppcNums = j.value("ppcNums", h_ppcNums);
	h_picHost = gridNx * gridNy * gridNz / hostNums * h_ppcNums;
	h_picDev = gridNx * gridNy * gridNz / hostNums / devNums * h_ppcNums;
	h_leftN = j.value("leftN", h_leftN);
	h_rightN = j.value("rightN", h_rightN);
	h_diagLeftX = j.value("diagLeftX", h_diagLeftX);
	h_diagRightX = j.value("diagRightX", h_diagRightX);
	h_diagY = j.value("diagY", h_diagY);
	h_OMEGACI = QE * h_B0 / MP;

	h_perp2Phi = j.value("perp2Phi", h_perp2Phi);
	h_para4Phi = j.value("para4Phi", h_para4Phi);
	h_perp2A = j.value("perp2A", h_perp2A);
	h_para4A = j.value("para4A", h_para4A);
	h_perp2dNe = j.value("perp2dNe", h_perp2dNe);
	h_para4dNe = j.value("para4dNe", h_para4dNe);
	h_perp2dTe = j.value("perp2dTe", h_perp2dTe);
	h_para4dTe = j.value("para4dTe", h_para4dTe);
	h_perp2dPi = j.value("perp2dPi", h_perp2dPi);
	h_para4dPi = j.value("para4dPi", h_para4dPi);
	h_perp2dPa = j.value("perp2dPa", h_perp2dPa);
	h_para4dPa = j.value("para4dPa", h_para4dPa);
	h_perp2dPb = j.value("perp2dPb", h_perp2dPb);
	h_para4dPb = j.value("para4dPb", h_para4dPb);

	h_IonMass = j.value("IonMass", h_IonMass);
	h_IonChar = j.value("IonChar", h_IonChar);
	h_IonBeta = j.value("IonBeta", h_IonBeta);
	h_IonVmin = j.value("IonVmin", h_IonVmin);
	h_IonVmax = j.value("IonVmax", h_IonVmax);
	h_IonVb = j.value("IonVb", h_IonVb);
	h_IonDeltaV = j.value("IonDeltaV", h_IonDeltaV);
	h_IonLambda0 = j.value("IonLambda0", h_IonLambda0);
	h_IonDeltaLambda2 = j.value("IonDeltaLambda2", h_IonDeltaLambda2);

	h_AlphaMass = j.value("AlphaMass", h_AlphaMass);
	h_AlphaChar = j.value("AlphaChar", h_AlphaChar);
	h_AlphaBeta = j.value("AlphaBeta", h_AlphaBeta);
	h_AlphaVmin = j.value("AlphaVmin", h_AlphaVmin);
	h_AlphaVmax = j.value("AlphaVmax", h_AlphaVmax);
	h_AlphaVb = j.value("AlphaVb", h_AlphaVb);
	h_AlphaDeltaV = j.value("AlphaDeltaV", h_AlphaDeltaV);
	h_AlphaLambda0 = j.value("AlphaLambda0", h_AlphaLambda0);
	h_AlphaDeltaLambda2 = j.value("AlphaDeltaLambda2", h_AlphaDeltaLambda2);

	h_BeamMass = j.value("BeamMass", h_BeamMass);
	h_BeamChar = j.value("BeamChar", h_BeamChar);
	h_BeamBeta = j.value("BeamBeta", h_BeamBeta);
	h_BeamVmin = j.value("BeamVmin", h_BeamVmin);
	h_BeamVmax = j.value("BeamVmax", h_BeamVmax);
	h_BeamVb = j.value("BeamVb", h_BeamVb);
	h_BeamDeltaV = j.value("BeamDeltaV", h_BeamDeltaV);
	h_BeamLambda0 = j.value("BeamLambda0", h_BeamLambda0);
	h_BeamDeltaLambda2 = j.value("BeamDeltaLambda2", h_BeamDeltaLambda2);
}

static void syncConfigToDeviceConstants() {
	const dataType normQE = QE / (h_B0 * h_L0 * h_L0 / MU0 / h_VA0);
	const dataType normMP = MP / (h_B0 * h_B0 * h_L0 * h_L0 * h_L0 / MU0 / h_VA0 / h_VA0);
	const dataType gridDzVal = 2.0 * PI / h_tubes / gridNz;
	const dataType gridDtVal = h_dt;
	const dataType zbegVal = -PI / h_tubes - (gridGhost - 0.5) * gridDzVal;
	const dataType zoriVal = -PI / h_tubes;
	const dataType zrangeVal = 2.0 * PI / h_tubes;
	const dataType rho0Val = h_RHO0;
	const dataType drhoVal = h_RHO1 - h_RHO0;
	const dataType psitmaxVal = h_PSITMAX / (h_B0 * h_L0 * h_L0);
	const dataType pitchB0Val = h_B0;
	const dataType vaVal = h_VA0;
	const dataType l3Val = 1e19 / (h_L0 * h_L0 * h_L0);
	const dataType cmVal = h_VA0 / (h_L0 * h_OMEGACI);

	CUDACHECK(cudaMemcpyToSymbol(B0, &h_B0, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(L0, &h_L0, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(VA0, &h_VA0, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(RHO0, &h_RHO0, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(RHO1, &h_RHO1, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(PSITMAX, &h_PSITMAX, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(OMEGACI, &h_OMEGACI, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(ppcNums, &h_ppcNums, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(picHost, &h_picHost, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(picDev, &h_picDev, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(tubes, &h_tubes, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(leftN, &h_leftN, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(rightN, &h_rightN, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(diagLeftX, &h_diagLeftX, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(diagRightX, &h_diagRightX, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(diagY, &h_diagY, sizeof(int)));
	CUDACHECK(cudaMemcpyToSymbol(gridDt, &gridDtVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(gridDz, &gridDzVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(zbeg, &zbegVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(zori, &zoriVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(zrange, &zrangeVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(NormQE, &normQE, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(NormMP, &normMP, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(rho0, &rho0Val, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(drho, &drhoVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(psitmax, &psitmaxVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(pitchB0, &pitchB0Val, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(va, &vaVal, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(l3, &l3Val, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(cm, &cmVal, sizeof(dataType)));

	CUDACHECK(cudaMemcpyToSymbol(perp2Phi, &h_perp2Phi, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(para4Phi, &h_para4Phi, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(perp2A, &h_perp2A, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(para4A, &h_para4A, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(perp2dNe, &h_perp2dNe, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(para4dNe, &h_para4dNe, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(perp2dTe, &h_perp2dTe, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(para4dTe, &h_para4dTe, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(perp2dPi, &h_perp2dPi, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(para4dPi, &h_para4dPi, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(perp2dPa, &h_perp2dPa, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(para4dPa, &h_para4dPa, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(perp2dPb, &h_perp2dPb, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(para4dPb, &h_para4dPb, sizeof(dataType)));

	CUDACHECK(cudaMemcpyToSymbol(IonMass, &h_IonMass, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonChar, &h_IonChar, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonBeta, &h_IonBeta, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonVmin, &h_IonVmin, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonVmax, &h_IonVmax, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonVb, &h_IonVb, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonDeltaV, &h_IonDeltaV, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonLambda0, &h_IonLambda0, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(IonDeltaLambda2, &h_IonDeltaLambda2, sizeof(dataType)));

	CUDACHECK(cudaMemcpyToSymbol(AlphaMass, &h_AlphaMass, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaChar, &h_AlphaChar, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaBeta, &h_AlphaBeta, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaVmin, &h_AlphaVmin, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaVmax, &h_AlphaVmax, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaVb, &h_AlphaVb, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaDeltaV, &h_AlphaDeltaV, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaLambda0, &h_AlphaLambda0, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(AlphaDeltaLambda2, &h_AlphaDeltaLambda2, sizeof(dataType)));

	CUDACHECK(cudaMemcpyToSymbol(BeamMass, &h_BeamMass, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamChar, &h_BeamChar, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamBeta, &h_BeamBeta, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamVmin, &h_BeamVmin, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamVmax, &h_BeamVmax, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamVb, &h_BeamVb, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamDeltaV, &h_BeamDeltaV, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamLambda0, &h_BeamLambda0, sizeof(dataType)));
	CUDACHECK(cudaMemcpyToSymbol(BeamDeltaLambda2, &h_BeamDeltaLambda2, sizeof(dataType)));
}

__constant__ dataType IonConst;
__constant__ dataType AlphaConst;
__constant__ dataType BeamConst;
__constant__ dataType hx[8] = { 1.0, 0.0, 1.0, 0.0,1.0, 0.0, 1.0, 0.0 };
__constant__ dataType sx[8] = { -1.0, 1.0, -1.0, 1.0,-1.0, 1.0, -1.0, 1.0 };
__constant__ dataType hy[8] = { 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0 };
__constant__ dataType sy[8] = { -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0 };
__constant__ dataType hz[8] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
__constant__ dataType sz[8] = { -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 };

/*-------------------------Thread Block Setting-------------------------*/

using ifLocal = std::conditional<gridNz == 8 || gridNz == 16 || gridNz == 32, trueType, falseType>::type;

const int MRK4BlockDimx = (gridNz == 8) ? 8 : (gridNz == 16) ? 16 : (gridNz == 32) ? 32 : 16;
const int MRK4BlockDimy = 2;
const int MRK4BlockDimz = 4;
const int MRK4GridDimx = gridNx / MRK4BlockDimz;
const int MRK4GridDimy = gridNy / hostNums / devNums / MRK4BlockDimy;
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

const int pptNums = 32;
const int PICBlockDimx = 256;
int h_PICGridDimx = h_picDev / pptNums / PICBlockDimx;

const int nFFTBatchSize = devNy * gridNx;
const int nFFTTimeSize = gridNz;
const int nFFTFreqSize = gridNz / 2 + 1;

const int mFFTBatchSize = gridNx * gridNz;
const int mFFTTimeSize = gridNy;
const int mFFTFreqSize = gridNy / 2 + 1;

/*-----------------------------------------------------------------Device Function-----------------------------------------------------------------*/

/*------------------------------------------------------MHD RK4------------------------------------------------------*/

template<typename local, typename type>
__device__ void Staggered2C(int offsetx, int offsetz, int& i, int& j, int& k, int& offset2d, int& offset3d, int& lane_id,
	int& shift_k, type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta, type qtheta_lr[4],
	type* address, type field_du[4], type field_lr[4], type& field) {

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

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
			field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
			field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
			field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}
	else {

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k + ((k + offsetz + shift_k - 1) % gridNz + gridNz) % gridNz];
			field_du[1] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k + ((k + offsetz + shift_k + 0) % gridNz + gridNz) % gridNz];
			field_du[2] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k + ((k + offsetz + shift_k + 1) % gridNz + gridNz) % gridNz];
			field_du[3] = address[offset3d + offsetx * gridNz + (index - 1) * gridNxz - k + ((k + offsetz + shift_k + 2) % gridNz + gridNz) % gridNz];

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}

	field = ((field_lr[1] + field_lr[2]) * 9 - (field_lr[0] + field_lr[3])) / 16;

	offset2d = j * gridNx + i;

}

template<typename local, typename type>
__device__ void Collocated2S(int offsetx, int offsetz, int& i, int& j, int& k, int& offset2d, int& offset3d, int& lane_id,
	int& shift_k, type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta, type qtheta_lr[4],
	type* address, type field_du[4], type field_lr[4], type& field) {

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

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
			field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
			field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
			field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}
	else {

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k + ((k + offsetz + shift_k - 1) % gridNz + gridNz) % gridNz];
			field_du[1] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k + ((k + offsetz + shift_k + 0) % gridNz + gridNz) % gridNz];
			field_du[2] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k + ((k + offsetz + shift_k + 1) % gridNz + gridNz) % gridNz];
			field_du[3] = address[offset3d + offsetx * gridNz + (index - 2) * gridNxz - k + ((k + offsetz + shift_k + 2) % gridNz + gridNz) % gridNz];

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}

	field = ((field_lr[1] + field_lr[2]) * 9 - (field_lr[0] + field_lr[3])) / 16;

	offset2d = j * gridNx + i;

}

template<typename type>
__device__ void PartialX(int offsetz, int& i, int& k, int& offset3d, type* address, type& field, type field_lr[4], type& field_px) {

	offsetz = -k + (k + offsetz + gridNz) % gridNz;

	if (i == 0) {
		field_lr[0] = address[offset3d + offsetz + 1 * gridNz];
		field_lr[1] = address[offset3d + offsetz + 2 * gridNz];
		field_lr[2] = address[offset3d + offsetz + 3 * gridNz];
		field_lr[3] = address[offset3d + offsetz + 4 * gridNz];
		field_px = (-25 * field + 48 * field_lr[0] - 36 * field_lr[1] + 16 * field_lr[2] - 3 * field_lr[3]) / (12 * gridDx);
	}
	else if (i == gridNx - 1) {
		field_lr[0] = address[offset3d + offsetz - 4 * gridNz];
		field_lr[1] = address[offset3d + offsetz - 3 * gridNz];
		field_lr[2] = address[offset3d + offsetz - 2 * gridNz];
		field_lr[3] = address[offset3d + offsetz - 1 * gridNz];
		field_px = (3 * field_lr[0] - 16 * field_lr[1] + 36 * field_lr[2] - 48 * field_lr[3] + 25 * field) / (12 * gridDx);
	}
	else if (i == 1) {
		field_lr[0] = address[offset3d + offsetz - 1 * gridNz];
		field_lr[1] = address[offset3d + offsetz + 1 * gridNz];
		field_lr[2] = address[offset3d + offsetz + 2 * gridNz];
		field_lr[3] = address[offset3d + offsetz + 3 * gridNz];
		field_px = (-3 * field_lr[0] - 10 * field + 18 * field_lr[1] - 6 * field_lr[2] + field_lr[3]) / (12 * gridDx);
	}
	else if (i == gridNx - 2) {
		field_lr[0] = address[offset3d + offsetz - 3 * gridNz];
		field_lr[1] = address[offset3d + offsetz - 2 * gridNz];
		field_lr[2] = address[offset3d + offsetz - 1 * gridNz];
		field_lr[3] = address[offset3d + offsetz + 1 * gridNz];
		field_px = (-field_lr[0] + 6 * field_lr[1] - 18 * field_lr[2] + 10 * field + 3 * field_lr[3]) / (12 * gridDx);
	}
	else {
		field_lr[0] = address[offset3d + offsetz - 2 * gridNz];
		field_lr[1] = address[offset3d + offsetz - 1 * gridNz];
		field_lr[2] = address[offset3d + offsetz + 1 * gridNz];
		field_lr[3] = address[offset3d + offsetz + 2 * gridNz];
		field_px = (field_lr[0] - 8 * field_lr[1] + 8 * field_lr[2] - field_lr[3]) / (12 * gridDx);
	}

}

template<typename type>
__device__ void PartialX2(int& i, int& offset3d, type* address, type& field, type field_lr[4], type& field_px, type& field_px2) {

	if (i == 0) {
		field_lr[0] = address[offset3d + 1 * gridNz];
		field_lr[1] = address[offset3d + 2 * gridNz];
		field_lr[2] = address[offset3d + 3 * gridNz];
		field_lr[3] = address[offset3d + 4 * gridNz];
		field_px = (-25 * field + 48 * field_lr[0] - 36 * field_lr[1] + 16 * field_lr[2] - 3 * field_lr[3]) / (12 * gridDx);
		field_px2 = (35 * field - 104 * field_lr[0] + 114 * field_lr[1] - 56 * field_lr[2] + 11 * field_lr[3]) / (12 * gridDx * gridDx);
	}
	else if (i == gridNx - 1) {
		field_lr[0] = address[offset3d - 4 * gridNz];
		field_lr[1] = address[offset3d - 3 * gridNz];
		field_lr[2] = address[offset3d - 2 * gridNz];
		field_lr[3] = address[offset3d - 1 * gridNz];
		field_px = (3 * field_lr[0] - 16 * field_lr[1] + 36 * field_lr[2] - 48 * field_lr[3] + 25 * field) / (12 * gridDx);
		field_px2 = (11 * field_lr[0] - 56 * field_lr[1] + 114 * field_lr[2] - 104 * field_lr[3] + 35 * field) / (12 * gridDx * gridDx);
	}
	else if (i == 1) {
		field_lr[0] = address[offset3d - 1 * gridNz];
		field_lr[1] = address[offset3d + 1 * gridNz];
		field_lr[2] = address[offset3d + 2 * gridNz];
		field_lr[3] = address[offset3d + 3 * gridNz];
		field_px = (-3 * field_lr[0] - 10 * field + 18 * field_lr[1] - 6 * field_lr[2] + field_lr[3]) / (12 * gridDx);
		field_px2 = (11 * field_lr[0] - 20 * field + 6 * field_lr[1] + 4 * field_lr[2] - field_lr[3]) / (12 * gridDx * gridDx);
	}
	else if (i == gridNx - 2) {
		field_lr[0] = address[offset3d - 3 * gridNz];
		field_lr[1] = address[offset3d - 2 * gridNz];
		field_lr[2] = address[offset3d - 1 * gridNz];
		field_lr[3] = address[offset3d + 1 * gridNz];
		field_px = (-field_lr[0] + 6 * field_lr[1] - 18 * field_lr[2] + 10 * field + 3 * field_lr[3]) / (12 * gridDx);
		field_px2 = (-field_lr[0] + 4 * field_lr[1] + 6 * field_lr[2] - 20 * field + 11 * field_lr[3]) / (12 * gridDx * gridDx);
	}
	else {
		field_lr[0] = address[offset3d - 2 * gridNz];
		field_lr[1] = address[offset3d - 1 * gridNz];
		field_lr[2] = address[offset3d + 1 * gridNz];
		field_lr[3] = address[offset3d + 2 * gridNz];
		field_px = (field_lr[0] - 8 * field_lr[1] + 8 * field_lr[2] - field_lr[3]) / (12 * gridDx);
		field_px2 = (-field_lr[0] + 16 * field_lr[1] - 30 * field + 16 * field_lr[2] - field_lr[3]) / (12 * gridDx * gridDx);
	}

}

template<typename local, typename type>
__device__ void PartialY(int& k, int& offset3d, int& lane_id, int& shift_k, type& shift_lk, type& shift_dk,
	type& qtheta, type qtheta_lr[4], type* address, type field_du[4], type field_lr[4], type& field_py) {

	if constexpr (std::is_same_v<local, trueType>) {

		field_lr[0] = address[offset3d - 2 * gridNxz];
		field_lr[1] = address[offset3d - 1 * gridNxz];
		field_lr[2] = address[offset3d + 1 * gridNxz];
		field_lr[3] = address[offset3d + 2 * gridNxz];

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
			field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
			field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
			field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}
	else {

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = address[offset3d + (index - 2 + index / 2) * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
			field_du[1] = address[offset3d + (index - 2 + index / 2) * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
			field_du[2] = address[offset3d + (index - 2 + index / 2) * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
			field_du[3] = address[offset3d + (index - 2 + index / 2) * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}

	field_py = (field_lr[0] - 8 * field_lr[1] + 8 * field_lr[2] - field_lr[3]) / (12 * gridDy);

}

template<typename local, typename type>
__device__ void PartialY2(int& k, int& offset3d, int& lane_id, int& shift_k, type& shift_lk, type& shift_dk,
	type& qtheta, type qtheta_lr[4], type* address, type& field, type field_du[4], type field_lr[4], type& field_py, type& field_py2) {

	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk,
		qtheta, qtheta_lr, address, field_du, field_lr, field_py);

	field_py2 = (-field_lr[0] + 16 * field_lr[1] - 30 * field + 16 * field_lr[2] - field_lr[3]) / (12 * gridDy * gridDy);

}

template<typename local, typename type>
__device__ void PartialZ(int& k, int& offset3d, int& lane_id, type* address, type& field, type field_du[4], type& field_pz) {

	if constexpr (std::is_same_v<local, trueType>) {

		field_du[0] = __shfl_sync(0xffffffff, field, lane_id - 2, gridNz);
		field_du[1] = __shfl_sync(0xffffffff, field, lane_id - 1, gridNz);
		field_du[2] = __shfl_sync(0xffffffff, field, lane_id + 1, gridNz);
		field_du[3] = __shfl_sync(0xffffffff, field, lane_id + 2, gridNz);

	}
	else {

		field_du[0] = address[offset3d - k + (k - 2 + gridNz) % gridNz];
		field_du[1] = address[offset3d - k + (k - 1 + gridNz) % gridNz];
		field_du[2] = address[offset3d - k + (k + 1 + gridNz) % gridNz];
		field_du[3] = address[offset3d - k + (k + 2 + gridNz) % gridNz];

	}

	field_pz = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * gridDz);

}

template<typename local, typename type>
__device__ void PartialZ2(int& k, int& offset3d, int& lane_id, type* address, type& field, type field_du[4], type& field_pz, type& field_pz2) {

	PartialZ<local>(k, offset3d, lane_id, address, field, field_du, field_pz);

	field_pz2 = (-field_du[0] + 16 * field_du[1] - 30 * field + 16 * field_du[2] - field_du[3]) / (12 * gridDz * gridDz);

}

template<typename local, typename type>
__device__ void S2CPartialXYZ(int& i, int& j, int& k, int& offset2d, int& offset3d, int& lane_id,
	int& shift_k, type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta, type qtheta_lr[4],
	type* address, type field_du[4], type field_lr[4], type& field, type& field_px, type& field_py, type& field_pz) {

	if (i == 0) {

		//+0
		Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = -25 * field;

		//+1
		Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 48 * field;

		//+2
		Staggered2C<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -36 * field;

		//+3
		Staggered2C<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 16 * field;

		//+4
		Staggered2C<local>(4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -3 * field;

		field_px /= (12 * gridDx);

	}
	else if (i == gridNx - 1) {

		//+0
		Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = 25 * field;

		//-4
		Staggered2C<local>(-4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 3 * field;

		//-3
		Staggered2C<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -16 * field;

		//-2
		Staggered2C<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 36 * field;

		//-1
		Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -48 * field;

		field_px /= (12 * gridDx);

	}
	else if (i == 1) {

		//+0
		Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = -10 * field;

		//-1
		Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -3 * field;

		//+1
		Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 18 * field;

		//+2
		Staggered2C<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -6 * field;

		//+3
		Staggered2C<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += field;

		field_px /= (12 * gridDx);

	}
	else if (i == gridNx - 2) {

		//+0
		Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = 10 * field;

		//-3
		Staggered2C<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -field;

		//-2
		Staggered2C<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 6 * field;

		//-1
		Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -18 * field;

		//+1
		Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 3 * field;

		field_px /= (12 * gridDx);

	}
	else {

		//-2
		Staggered2C<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = field;

		//-1
		Staggered2C<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -8 * field;

		//+1
		Staggered2C<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 8 * field;

		//+2
		Staggered2C<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -field;

		field_px /= (12 * gridDx);

	}

	Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
		shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);

	field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * gridDy);

	if constexpr (std::is_same_v<local, trueType>) {

		PartialZ<local>(k, offset3d, lane_id, address, field, field_du, field_pz);

	}
	else {

		Staggered2C<local>(0, -2, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz = field;

		Staggered2C<local>(0, -1, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz += -8 * field;

		Staggered2C<local>(0, 1, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz += 8 * field;

		Staggered2C<local>(0, 2, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz += -field;

		field_pz /= (12 * gridDz);

		Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);

	}

}

template<typename local, typename type>
__device__ void C2SPartialXYZ(int& i, int& j, int& k, int& offset2d, int& offset3d, int& lane_id,
	int& shift_k, type& shift_lk, type& shift_dk, type* d_qtheta, type& qtheta, type qtheta_lr[4],
	type* address, type field_du[4], type field_lr[4], type& field, type& field_px, type& field_py, type& field_pz) {

	if (i == 0) {

		//+0
		Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = -25 * field;

		//+1
		Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 48 * field;

		//+2
		Collocated2S<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -36 * field;

		//+3
		Collocated2S<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 16 * field;

		//+4
		Collocated2S<local>(4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -3 * field;

		field_px /= (12 * gridDx);

	}
	else if (i == gridNx - 1) {

		//+0
		Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = 25 * field;

		//-4
		Collocated2S<local>(-4, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 3 * field;

		//-3
		Collocated2S<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -16 * field;

		//-2
		Collocated2S<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 36 * field;

		//-1
		Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -48 * field;

		field_px /= (12 * gridDx);

	}
	else if (i == 1) {

		//+0
		Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = -10 * field;

		//-1
		Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -3 * field;

		//+1
		Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 18 * field;

		//+2
		Collocated2S<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -6 * field;

		//+3
		Collocated2S<local>(3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += field;

		field_px /= (12 * gridDx);

	}
	else if (i == gridNx - 2) {

		//+0
		Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = 10 * field;

		//-3
		Collocated2S<local>(-3, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -field;

		//-2
		Collocated2S<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 6 * field;

		//-1
		Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -18 * field;

		//+1
		Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 3 * field;

		field_px /= (12 * gridDx);

	}
	else {

		//-2
		Collocated2S<local>(-2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px = field;

		//-1
		Collocated2S<local>(-1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -8 * field;

		//+1
		Collocated2S<local>(1, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += 8 * field;

		//+2
		Collocated2S<local>(2, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_px += -field;

		field_px /= (12 * gridDx);

	}

	Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
		shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);

	field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * gridDy);

	if constexpr (std::is_same_v<local, trueType>) {

		PartialZ<local>(k, offset3d, lane_id, address, field, field_du, field_pz);

	}
	else {

		Collocated2S<local>(0, -2, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz = field;

		Collocated2S<local>(0, -1, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz += -8 * field;

		Collocated2S<local>(0, 1, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz += 8 * field;

		Collocated2S<local>(0, 2, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);
		field_pz += -field;

		field_pz /= (12 * gridDz);

		Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, address, field_du, field_lr, field);

	}

}

template<int rk4, typename nonlinear, typename local, typename staggered, typename Eparallel, typename type>
__global__ void MHDLinearRK4(type* __restrict__  d_qtheta,
	type* __restrict__ w_beg, type* __restrict__ w_midl, type* __restrict__ w_midr, type* __restrict__ w_end,
	type* __restrict__ A_beg, type* __restrict__ A_midl, type* __restrict__ A_midr, type* __restrict__ A_end,
	type* __restrict__ dNe_beg, type* __restrict__ dNe_midl, type* __restrict__ dNe_midr, type* __restrict__ dNe_end,
	type* __restrict__ dTe_beg, type* __restrict__ dTe_midl, type* __restrict__ dTe_midr, type* __restrict__ dTe_end,
	type* __restrict__ Phi_mid, type* __restrict__ dJpB_mid, type* __restrict__ dPe_mid, type* __restrict__ dPi_mid, type* __restrict__ dPa_mid, type* __restrict__ dPb_mid,
	type* __restrict__ d_wdPAdJpB2w, type* __restrict__ d_APhidNe2A, type* __restrict__ d_dPePhiAdJpB2dNe, type* __restrict__ d_PhidTedNe2dTe) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = (j + gridGhost) * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	type qtheta;
	type qtheta_lr[4];
	int shift_k;
	type shift_lk;
	type shift_dk;

	/*-----------------------Field and Derivative-----------------------*/

	type field;
	type field_px, field_py, field_pz;
	type field_du[4];
	type field_lr[4];

	/*---------------------Compressed Coefficient---------------------*/

	type compcoes[6];

	/*--------------------------RK4 Variables--------------------------*/

	type w_begin, dwdt;
	type A_begin, dAdt;
	type dNe_begin, dNedt;
	type dTe_begin, dTedt;

	/*-----------------------------Initialize-----------------------------*/

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

	dwdt = 0; dAdt = 0; dNedt = 0; dTedt = 0;

	/*--------Electron Pressure in Vorticity and Electron Density--------*/

	field = dPe_mid[offset3d];

	PartialZ<local>(k, offset3d, lane_id, dPe_mid, field, field_du, field_pz);
	PartialX(0, i, k, offset3d, dPe_mid, field, field_lr, field_px);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPe_mid, field_du, field_lr, field_py);

	for (int index = 0; index < 3; index++)
		compcoes[index] = d_dPePhiAdJpB2dNe[offset2d * 11 + index];

	dNedt += compcoes[0] * field_px + compcoes[1] * field_py + compcoes[2] * field_pz;

	for (int index = 0; index < 5; index++)
		compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index];

	dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;

	/*----------------Ion Diamagnetic Drift in Vorticity----------------*/

	field = w_midl[offset3d];

	PartialZ<local>(k, offset3d, lane_id, w_midl, field, field_du, field_pz);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, w_midl, field_du, field_lr, field_py);

	dwdt += compcoes[0] * field_py + compcoes[1] * field_pz;

	/*---------------------Ion Pressure in Vorticity---------------------*/

	if constexpr (std::is_same_v<ifIon, trueType>) {

		field = dPi_mid[offset3d];

		PartialZ<local>(k, offset3d, lane_id, dPi_mid, field, field_du, field_pz);
		PartialX(0, i, k, offset3d, dPi_mid, field, field_lr, field_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPi_mid, field_du, field_lr, field_py);

		dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;

	}

	/*--------------------Alpha Pressure in Vorticity--------------------*/

	if constexpr (std::is_same_v<ifAlpha, trueType>) {

		field = dPa_mid[offset3d];

		PartialZ<local>(k, offset3d, lane_id, dPa_mid, field, field_du, field_pz);
		PartialX(0, i, k, offset3d, dPa_mid, field, field_lr, field_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPa_mid, field_du, field_lr, field_py);

		dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;

	}

	/*--------------------Beam Pressure in Vorticity--------------------*/

	if constexpr (std::is_same_v<ifBeam, trueType>) {

		field = dPb_mid[offset3d];

		PartialZ<local>(k, offset3d, lane_id, dPb_mid, field, field_du, field_pz);
		PartialX(0, i, k, offset3d, dPb_mid, field, field_lr, field_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPb_mid, field_du, field_lr, field_py);

		dwdt += compcoes[2] * field_px + compcoes[3] * field_py + compcoes[4] * field_pz;

	}

	/*------Electric Potential in Electron Density and Temperature------*/

	field = Phi_mid[offset3d];

	PartialZ<local>(k, offset3d, lane_id, Phi_mid, field, field_du, field_pz);
	PartialX(0, i, k, offset3d, Phi_mid, field, field_lr, field_px);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, field_py);

	for (int index = 0; index < 3; index++)
		compcoes[index] = d_dPePhiAdJpB2dNe[offset2d * 11 + index + 3];

	dNedt += compcoes[0] * field_px + compcoes[1] * field_py + compcoes[2] * field_pz;

	if constexpr (std::is_same_v<nonlinear, falseType>) {

		for (int index = 0; index < 6; index++)
			compcoes[index] = d_PhidTedNe2dTe[offset2d * 6 + index];

		dTedt += compcoes[0] * field_py + compcoes[1] * field_pz;

	}

	/*---------Electron Temperature in Electron Temperature---------*/

	if constexpr (std::is_same_v<nonlinear, falseType>) {

		field = dTe_midl[offset3d];

		PartialZ<local>(k, offset3d, lane_id, dTe_midl, field, field_du, field_pz);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dTe_midl, field_du, field_lr, field_py);

		dTedt += compcoes[2] * field_py + compcoes[3] * field_pz;

	}

	/*------------Electron Density in Electron Temperature------------*/

	if constexpr (std::is_same_v<nonlinear, falseType>) {

		field = dNe_midl[offset3d];

		PartialZ<local>(k, offset3d, lane_id, dNe_midl, field, field_du, field_pz);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl, field_du, field_lr, field_py);

		dTedt += compcoes[4] * field_py + compcoes[5] * field_pz;

	}

	/*-------Parallel Vector Potential in Parallel Vector Potential-------*/

	for (int index = 0; index < 5; index++)
		compcoes[index] = d_APhidNe2A[offset2d * 5 + index];

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		if constexpr (std::is_same_v<nonlinear, falseType>) {

			field = A_midl[offset3d];
			PartialZ<local>(k, offset3d, lane_id, A_midl, field, field_du, field_pz);
			PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du, field_lr, field_py);
			dAdt += compcoes[0] * field + compcoes[1] * field_py + compcoes[2] * field_pz;

		}

	}

	/*-----------Electric Potential in Parallel Vector Potential-----------*/

	if constexpr (std::is_same_v<staggered, trueType>) {

		Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, field);
		field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * gridDy);

	}
	else {

		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, field_py);

	}

	dAdt += compcoes[3] * field_py;

	/*-----------Electron Density in Parallel Vector Potential-----------*/

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		if constexpr (std::is_same_v<nonlinear, falseType>) {

			if constexpr (std::is_same_v<staggered, trueType>) {

				Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
					shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, dNe_midl, field_du, field_lr, field);
				field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * gridDy);

			}
			else {

				PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl, field_du, field_lr, field_py);

			}

			dAdt += compcoes[4] * field_py;

		}

	}

	/*----Parallel Vector Potential in Vorticity and Electron Density----*/

	if constexpr (std::is_same_v<staggered, trueType>) {

		S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
			d_qtheta, qtheta, qtheta_lr, A_midl, field_du, field_lr, field, field_px, field_py, field_pz);

		for (int index = 0; index < 4; index++)
			compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 5];

		dwdt += compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz;
		dNedt += (compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz) / NormQE;

	}
	else {

		field = A_midl[offset3d];

		PartialZ<local>(k, offset3d, lane_id, A_midl, field, field_du, field_pz);
		PartialX(0, i, k, offset3d, A_midl, field, field_lr, field_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du, field_lr, field_py);

		for (int index = 0; index < 4; index++)
			compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 5];

		dwdt += compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz;
		dNedt += (compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_py + compcoes[3] * field_pz) / NormQE;

	}

	/*-------------------Parallel Current in Vorticity-------------------*/

	if constexpr (std::is_same_v<staggered, trueType>) {

		Staggered2C<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, dJpB_mid, field_du, field_lr, field);
		field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * gridDy);

		for (int index = 0; index < 1; index++)
			compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 9];

		dwdt += compcoes[0] * field_py;
		dNedt += compcoes[0] * field_py / NormQE;

	}
	else {

		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dJpB_mid, field_du, field_lr, field_py);

		for (int index = 0; index < 1; index++)
			compcoes[index] = d_wdPAdJpB2w[offset2d * 10 + index + 9];

		dwdt += compcoes[0] * field_py;
		dNedt += compcoes[0] * field_py / NormQE;

	}

	/*-------------------------------RK4-------------------------------*/

	if constexpr (rk4 == 1) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] = w_begin + dwdt * gridDt / 2;
			w_end[offset3d] = w_begin + dwdt * gridDt / 6;

			A_midr[offset3d] = A_begin + dAdt * gridDt / 2;
			A_end[offset3d] = A_begin + dAdt * gridDt / 6;

			dNe_midr[offset3d] = dNe_begin + dNedt * gridDt / 2;
			dNe_end[offset3d] = dNe_begin + dNedt * gridDt / 6;

			dTe_midr[offset3d] = dTe_begin + dTedt * gridDt / 2;
			dTe_end[offset3d] = dTe_begin + dTedt * gridDt / 6;
		}

	}
	else if constexpr (rk4 == 2) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] = w_begin + dwdt * gridDt / 2;
			w_end[offset3d] += dwdt * gridDt / 3;

			A_midr[offset3d] = A_begin + dAdt * gridDt / 2;
			A_end[offset3d] += dAdt * gridDt / 3;

			dNe_midr[offset3d] = dNe_begin + dNedt * gridDt / 2;
			dNe_end[offset3d] += dNedt * gridDt / 3;

			dTe_midr[offset3d] = dTe_begin + dTedt * gridDt / 2;
			dTe_end[offset3d] += dTedt * gridDt / 3;
		}

	}
	else if constexpr (rk4 == 3) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] = w_begin + dwdt * gridDt;
			w_end[offset3d] += dwdt * gridDt / 3;

			A_midr[offset3d] = A_begin + dAdt * gridDt;
			A_end[offset3d] += dAdt * gridDt / 3;

			dNe_midr[offset3d] = dNe_begin + dNedt * gridDt;
			dNe_end[offset3d] += dNedt * gridDt / 3;

			dTe_midr[offset3d] = dTe_begin + dTedt * gridDt;
			dTe_end[offset3d] += dTedt * gridDt / 3;
		}

	}
	else if constexpr (rk4 == 4) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] = w_end[offset3d] + dwdt * gridDt / 6;
			w_end[offset3d] += dwdt * gridDt / 6;

			A_midr[offset3d] = A_end[offset3d] + dAdt * gridDt / 6;
			A_end[offset3d] += dAdt * gridDt / 6;

			dNe_midr[offset3d] = dNe_end[offset3d] + dNedt * gridDt / 6;
			dNe_end[offset3d] += dNedt * gridDt / 6;

			dTe_midr[offset3d] = dTe_end[offset3d] + dTedt * gridDt / 6;
			dTe_end[offset3d] += dTedt * gridDt / 6;
		}

	}

}

template<int rk4, typename local, typename staggered, typename Eparallel, typename type>
__global__ void MHDNonlinearRK4(type* __restrict__  d_qtheta,
	type* __restrict__ w_midl, type* __restrict__ w_midr, type* __restrict__ w_end,
	type* __restrict__ A_midl, type* __restrict__ A_midr, type* __restrict__ A_end,
	type* __restrict__ dNe_midl, type* __restrict__ dNe_midr, type* __restrict__ dNe_end,
	type* __restrict__ dTe_midl, type* __restrict__ dTe_midr, type* __restrict__ dTe_end,
	type* __restrict__ Phi_mid, type* __restrict__ dJpB_mid, type* __restrict__ dPe_mid,
	type* __restrict__ d_Ne0, type* __restrict__ d_Te0, type* __restrict__ d_Ne0_px, type* __restrict__ d_Te0_px, type* __restrict__ d_APhidNe2A,
	type* __restrict__ d_wPhi_w, type* __restrict__ d_AdJpB_w, type* __restrict__ d_PhiA_A, type* __restrict__ d_NeA_A,
	type* __restrict__ d_dNePhi_dNe, type* __restrict__ d_PhiTe_dTe, type* __restrict__ d_PhiTeA_dTe) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = j * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	type qtheta;
	type qtheta_lr[4];
	int shift_k;
	type shift_lk;
	type shift_dk;

	/*-----------------------Field and Derivative-----------------------*/

	type Ne0, Te0, Ne0_px, Te0_px, Pe0_px, Ne, Te;
	type A, A_px, A_py, A_pz;
	type Phi, Phi_px, Phi_py, Phi_pz;
	type field, field_px, field_py, field_pz;
	type field_du[4];
	type field_lr[4];

	/*---------------------Compressed Coefficient---------------------*/

	type compcoes[9];

	/*--------------------------RK4 Variables--------------------------*/

	type dwdt, dAdt, dNedt, dTedt;

	/*-----------------------------Initialize-----------------------------*/

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

	dwdt = 0; dAdt = 0; dNedt = 0; dTedt = 0;

	Phi = Phi_mid[offset3d];
	PartialZ<local>(k, offset3d, lane_id, Phi_mid, Phi, field_du, Phi_pz);
	PartialX(0, i, k, offset3d, Phi_mid, Phi, field_lr, Phi_px);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, Phi_py);

	if constexpr (std::is_same_v<staggered, falseType>) {

		A = A_midl[offset3d];
		PartialZ<local>(k, offset3d, lane_id, A_midl, A, field_du, A_pz);
		PartialX(0, i, k, offset3d, A_midl, A, field_lr, A_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du, field_lr, A_py);

		field = dJpB_mid[offset3d];
		PartialZ<local>(k, offset3d, lane_id, dJpB_mid, field, field_du, field_pz);
		PartialX(0, i, k, offset3d, dJpB_mid, field, field_lr, field_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dJpB_mid, field_du, field_lr, field_py);

	}
	else {

		S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
			d_qtheta, qtheta, qtheta_lr, A_midl, field_du, field_lr, A, A_px, A_py, A_pz);

		S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
			d_qtheta, qtheta, qtheta_lr, dJpB_mid, field_du, field_lr, field, field_px, field_py, field_pz);

		offset2d = (j + gridGhost) * gridNx + i;

		qtheta = d_qtheta[offset2d];
		qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
		qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
		qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
		qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

		offset2d = j * gridNx + i;

	}

	/*------------------------AdJpB in Vorticity------------------------*/

	for (int index = 0; index < 9; index++)
		compcoes[index] = d_AdJpB_w[offset2d * 9 + index];

	dwdt += compcoes[0] * A * field_px + compcoes[1] * A * field_py + compcoes[2] * A * field_pz
		+ compcoes[3] * A_px * field_py + compcoes[4] * A_px * field_pz
		+ compcoes[5] * A_py * field_px + compcoes[6] * A_py * field_pz
		+ compcoes[7] * A_pz * field_px + compcoes[8] * A_pz * field_py;

	/*--------------------AdJpB in Electron Density--------------------*/

	dNedt += (compcoes[0] * A * field_px + compcoes[1] * A * field_py + compcoes[2] * A * field_pz
		+ compcoes[3] * A_px * field_py + compcoes[4] * A_px * field_pz
		+ compcoes[5] * A_py * field_px + compcoes[6] * A_py * field_pz
		+ compcoes[7] * A_pz * field_px + compcoes[8] * A_pz * field_py) / NormQE;

	/*-------------------------wPhi in Vorticity-------------------------*/

	for (int index = 0; index < 6; index++)
		compcoes[index] = d_wPhi_w[offset2d * 6 + index];

	field = w_midl[offset3d];
	PartialZ<local>(k, offset3d, lane_id, w_midl, field, field_du, field_pz);
	PartialX(0, i, k, offset3d, w_midl, field, field_lr, field_px);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, w_midl, field_du, field_lr, field_py);

	dwdt += compcoes[0] * field_px * Phi_py + compcoes[1] * field_px * Phi_pz
		+ compcoes[2] * field_py * Phi_px + compcoes[3] * field_py * Phi_pz
		+ compcoes[4] * field_pz * Phi_px + compcoes[5] * field_pz * Phi_py;

	/*-------------------dNePhi in Electron Density-------------------*/

	for (int index = 0; index < 9; index++)
		compcoes[index] = d_dNePhi_dNe[offset2d * 9 + index];

	field = dNe_midl[offset3d];
	PartialZ<local>(k, offset3d, lane_id, dNe_midl, field, field_du, field_pz);
	PartialX(0, i, k, offset3d, dNe_midl, field, field_lr, field_px);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl, field_du, field_lr, field_py);

	dNedt += compcoes[0] * field * Phi_px + compcoes[1] * field * Phi_py + compcoes[2] * field * Phi_pz
		+ compcoes[3] * field_px * Phi_py + compcoes[4] * field_px * Phi_pz
		+ compcoes[5] * field_py * Phi_px + compcoes[6] * field_py * Phi_pz
		+ compcoes[7] * field_pz * Phi_px + compcoes[8] * field_pz * Phi_py;

	Ne = Ne0 + field;

	/*-----------------PhiTeA in Electron Temperature-----------------*/

	field = dTe_midl[offset3d];
	PartialZ<local>(k, offset3d, lane_id, dTe_midl, field, field_du, field_pz);
	PartialX(0, i, k, offset3d, dTe_midl, field, field_lr, field_px);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dTe_midl, field_du, field_lr, field_py);
	field_px += Te0_px;

	Te = Te0 + field;

	for (int index = 0; index < 6; index++)
		compcoes[index] = d_PhiTe_dTe[offset2d * 6 + index];

	dTedt += compcoes[0] * Phi_px * field_py + compcoes[1] * Phi_px * field_pz
		+ compcoes[2] * Phi_py * field_px + compcoes[3] * Phi_py * field_pz
		+ compcoes[4] * Phi_pz * field_px + compcoes[5] * Phi_pz * field_py;

	for (int index = 0; index < 9; index++)
		compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index];

	dTedt += compcoes[0] * Phi_px * field_py * A + compcoes[1] * Phi_px * field_py * A_px
		+ compcoes[2] * Phi_px * field_py * A_py + compcoes[3] * Phi_px * field_py * A_pz;
	dTedt += compcoes[4] * Phi_px * field_pz * A;
	dTedt += compcoes[5] * Phi_py * field_px * A + compcoes[6] * Phi_py * field_px * A_px
		+ compcoes[7] * Phi_py * field_px * A_py + compcoes[8] * Phi_py * field_px * A_pz;

	for (int index = 0; index < 9; index++)
		compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index + 9];

	dTedt += compcoes[0] * Phi_py * field_pz * A + compcoes[1] * Phi_py * field_pz * A_px
		+ compcoes[2] * Phi_py * field_pz * A_py + compcoes[3] * Phi_py * field_pz * A_pz;
	dTedt += compcoes[4] * Phi_pz * field_px * A;
	dTedt += compcoes[5] * Phi_pz * field_py * A + compcoes[6] * Phi_pz * field_py * A_px
		+ compcoes[7] * Phi_pz * field_py * A_py + compcoes[8] * Phi_pz * field_py * A_pz;

	/*-----------------PeTeA in Electron Temperature-----------------*/

	Phi = dPe_mid[offset3d];
	PartialZ<local>(k, offset3d, lane_id, dPe_mid, Phi, field_du, Phi_pz);
	PartialX(0, i, k, offset3d, dPe_mid, Phi, field_lr, Phi_px);
	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dPe_mid, field_du, field_lr, Phi_py);
	Phi_px += Pe0_px;

	for (int index = 0; index < 6; index++)
		compcoes[index] = d_PhiTe_dTe[offset2d * 6 + index];

	dTedt += -(compcoes[0] * Phi_px * field_py + compcoes[1] * Phi_px * field_pz
		+ compcoes[2] * Phi_py * field_px + compcoes[3] * Phi_py * field_pz
		+ compcoes[4] * Phi_pz * field_px + compcoes[5] * Phi_pz * field_py) / (2 * NormQE * Ne);

	for (int index = 0; index < 9; index++)
		compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index];

	dTedt += -(compcoes[0] * Phi_px * field_py * A + compcoes[1] * Phi_px * field_py * A_px
		+ compcoes[2] * Phi_px * field_py * A_py + compcoes[3] * Phi_px * field_py * A_pz) / (2 * NormQE * Ne);
	dTedt += -(compcoes[4] * Phi_px * field_pz * A) / (2 * NormQE * Ne);
	dTedt += -(compcoes[5] * Phi_py * field_px * A + compcoes[6] * Phi_py * field_px * A_px
		+ compcoes[7] * Phi_py * field_px * A_py + compcoes[8] * Phi_py * field_px * A_pz) / (2 * NormQE * Ne);

	for (int index = 0; index < 9; index++)
		compcoes[index] = d_PhiTeA_dTe[offset2d * 18 + index + 9];

	dTedt += -(compcoes[0] * Phi_py * field_pz * A + compcoes[1] * Phi_py * field_pz * A_px
		+ compcoes[2] * Phi_py * field_pz * A_py + compcoes[3] * Phi_py * field_pz * A_pz) / (2 * NormQE * Ne);
	dTedt += -(compcoes[4] * Phi_pz * field_px * A) / (2 * NormQE * Ne);
	dTedt += -(compcoes[5] * Phi_pz * field_py * A + compcoes[6] * Phi_pz * field_py * A_px
		+ compcoes[7] * Phi_pz * field_py * A_py + compcoes[8] * Phi_pz * field_py * A_pz) / (2 * NormQE * Ne);

	/*-----------------PhiA in Parallel Vector Potential-----------------*/

	for (int index = 0; index < 9; index++)
		compcoes[index] = d_PhiA_A[offset2d * 9 + index];

	if constexpr (std::is_same_v<staggered, trueType>) {

		A = A_midl[offset3d];
		PartialZ<local>(k, offset3d, lane_id, A_midl, A, field_du, A_pz);
		PartialX(0, i, k, offset3d, A_midl, A, field_lr, A_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_midl, field_du, field_lr, A_py);

		C2SPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
			d_qtheta, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, Phi, Phi_px, Phi_py, Phi_pz);

	}
	else {

		Phi = Phi_mid[offset3d];
		PartialZ<local>(k, offset3d, lane_id, Phi_mid, Phi, field_du, Phi_pz);
		PartialX(0, i, k, offset3d, Phi_mid, Phi, field_lr, Phi_px);
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, Phi_py);

	}

	dAdt += compcoes[0] * Phi_px * A + compcoes[1] * Phi_px * A_py + compcoes[2] * Phi_px * A_pz
		+ compcoes[3] * Phi_py * A + compcoes[4] * Phi_py * A_px + compcoes[5] * Phi_py * A_pz
		+ compcoes[6] * Phi_pz * A + compcoes[7] * Phi_pz * A_px + compcoes[8] * Phi_pz * A_py;

	/*-----------------dNe in Parallel Vector Potential-----------------*/

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		compcoes[0] = d_APhidNe2A[offset2d * 5 + 4];

		if constexpr (std::is_same_v<staggered, falseType>) {

			field = dNe_midl[offset3d];
			PartialZ<local>(k, offset3d, lane_id, dNe_midl, field, field_du, field_pz);
			PartialX(0, i, k, offset3d, dNe_midl, field, field_lr, field_px);
			PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_midl, field_du, field_lr, field_py);

		}
		else {

			C2SPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
				d_qtheta, qtheta, qtheta_lr, dNe_midl, field_du, field_lr, field, field_px, field_py, field_pz);

			Collocated2S<local>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
				shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, dTe_midl, field_du, field_lr, Te);

			Ne = Ne0 + field;
			Te = Te0 + Te;

		}

		dAdt += compcoes[0] * (Te / Te0 * Ne0 / Ne) * field_py;

	}

	/*-----------------NeA in Parallel Vector Potential-----------------*/

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		for (int index = 0; index < 9; index++)
			compcoes[index] = d_NeA_A[offset2d * 9 + index];

		field_px += Ne0_px;

		dAdt += (compcoes[0] * field_px * A + compcoes[1] * field_px * A_py + compcoes[2] * field_px * A_pz
			+ compcoes[3] * field_py * A + compcoes[4] * field_py * A_px + compcoes[5] * field_py * A_pz
			+ compcoes[6] * field_pz * A + compcoes[7] * field_pz * A_px + compcoes[8] * field_pz * A_py) * (Te / Te0 * Ne0 / Ne);

	}

	/*-------------------------------RK4-------------------------------*/

	if constexpr (rk4 == 1) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] += dwdt * gridDt / 2;
			w_end[offset3d] += dwdt * gridDt / 6;

			A_midr[offset3d] += dAdt * gridDt / 2;
			A_end[offset3d] += dAdt * gridDt / 6;

			dNe_midr[offset3d] += dNedt * gridDt / 2;
			dNe_end[offset3d] += dNedt * gridDt / 6;

			dTe_midr[offset3d] += dTedt * gridDt / 2;
			dTe_end[offset3d] += dTedt * gridDt / 6;
		}

	}
	else if constexpr (rk4 == 2) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] += dwdt * gridDt / 2;
			w_end[offset3d] += dwdt * gridDt / 3;

			A_midr[offset3d] += dAdt * gridDt / 2;
			A_end[offset3d] += dAdt * gridDt / 3;

			dNe_midr[offset3d] += dNedt * gridDt / 2;
			dNe_end[offset3d] += dNedt * gridDt / 3;

			dTe_midr[offset3d] += dTedt * gridDt / 2;
			dTe_end[offset3d] += dTedt * gridDt / 3;
		}

	}
	else if constexpr (rk4 == 3) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] += dwdt * gridDt;
			w_end[offset3d] += dwdt * gridDt / 3;

			A_midr[offset3d] += dAdt * gridDt;
			A_end[offset3d] += dAdt * gridDt / 3;

			dNe_midr[offset3d] += dNedt * gridDt;
			dNe_end[offset3d] += dNedt * gridDt / 3;

			dTe_midr[offset3d] += dTedt * gridDt;
			dTe_end[offset3d] += dTedt * gridDt / 3;
		}

	}
	else if constexpr (rk4 == 4) {

		if (i != 0 && i != gridNx - 1) {
			w_midr[offset3d] += dwdt * gridDt / 6;
			w_end[offset3d] += dwdt * gridDt / 6;

			A_midr[offset3d] += dAdt * gridDt / 6;
			A_end[offset3d] += dAdt * gridDt / 6;

			dNe_midr[offset3d] += dNedt * gridDt / 6;
			dNe_end[offset3d] += dNedt * gridDt / 6;

			dTe_midr[offset3d] += dTedt * gridDt / 6;
			dTe_end[offset3d] += dTedt * gridDt / 6;
		}

	}

}

template<typename nonlinear, typename local, typename FLRMHD, typename type>
__global__ void MHD2dJpBdPePhi(
	type* __restrict__ A_mid, type* __restrict__ dJpB_mid, type* __restrict__ d_A2dJpB,
	type* __restrict__ w_mid, type* __restrict__ Phi_mid, type* __restrict__ d_w2Phi,
	type* __restrict__ dNe_mid, type* __restrict__ dTe_mid, type* __restrict__ dPe_mid,
	type* __restrict__ d_Ne0, type* __restrict__ d_Te0) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = j * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*-----------------------Field and Derivative-----------------------*/

	type dJpB, Ne0, Te0, dNe, dTe;

	type field;
	type field_px, field_pz, field_px2, field_pxz, field_pz2;
	type field_du[4];
	type field_lr[4];

	/*---------------------Compressed Coefficient---------------------*/

	type compcoes[6];

	/*-----------------------------Initialize-----------------------------*/

	for (int index = 0; index < 6; index++)
		compcoes[index] = d_A2dJpB[offset2d * 6 + index];

	Ne0 = d_Ne0[offset2d];
	Te0 = d_Te0[offset2d];
	dNe = dNe_mid[offset3d];
	dTe = dTe_mid[offset3d];

	/*-----------Parallel Vector Potential in Parallel Current-----------*/

	field = A_mid[offset3d];

	PartialZ2<local>(k, offset3d, lane_id, A_mid, field, field_du, field_pz, field_pz2);
	PartialX2(i, offset3d, A_mid, field, field_lr, field_px, field_px2);

	dJpB = compcoes[0] * field + compcoes[1] * field_px + compcoes[2] * field_pz
		+ compcoes[3] * field_px2 + compcoes[5] * field_pz2;

	if constexpr (std::is_same_v<local, trueType>) {

		PartialZ<local>(k, offset3d, lane_id, A_mid, field_px, field_du, field_pxz);

	}
	else {

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

		field_pxz /= (12 * gridDz);

	}

	dJpB += compcoes[4] * field_pxz;

	if (i != 0 && i != gridNx - 1)
		dJpB_mid[offset3d] = dJpB;

	/*---------------------------FLR Effect---------------------------*/

	if constexpr (std::is_same_v<FLRMHD, trueType>)
		Phi_mid[offset3d] += d_w2Phi[offset2d] * w_mid[offset3d];

	/*-----------Electron Density, Temperature and Pressure-----------*/

	if constexpr (std::is_same_v<nonlinear, falseType>)
		dPe_mid[offset3d] = dNe * Te0 + Ne0 * dTe;
	else
		dPe_mid[offset3d] = dNe * Te0 + Ne0 * dTe + dNe * dTe;

}

template< typename local, typename type>
__global__ void MHD2w(
	type* __restrict__ Phi_mid, type* __restrict__ w_mid, type* __restrict__ d_Phi2w) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = j * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*-----------------------Field and Derivative-----------------------*/

	type w;

	type field;
	type field_px, field_pz, field_px2, field_pxz, field_pz2;
	type field_du[4];
	type field_lr[4];

	/*---------------------Compressed Coefficient---------------------*/

	type compcoes[5];

	/*-----------------------------Initialize-----------------------------*/

	for (int index = 0; index < 5; index++)
		compcoes[index] = d_Phi2w[offset2d * 5 + index];

	/*------------------Electric Potential in Vorticity------------------*/

	field = Phi_mid[offset3d];

	PartialZ2<local>(k, offset3d, lane_id, Phi_mid, field, field_du, field_pz, field_pz2);
	PartialX2(i, offset3d, Phi_mid, field, field_lr, field_px, field_px2);

	w = compcoes[0] * field_px + compcoes[1] * field_pz + compcoes[2] * field_px2 + compcoes[4] * field_pz2;

	if constexpr (std::is_same_v<local, trueType>) {

		PartialZ<local>(k, offset3d, lane_id, Phi_mid, field_px, field_du, field_pxz);

	}
	else {

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

		field_pxz /= (12 * gridDz);

	}

	w += compcoes[3] * field_pxz;

	if (i != 0 && i != gridNx - 1)
		w_mid[offset3d] = w;

}

/*------------------------------------------------MHD Boundary Ghost------------------------------------------------*/

template<typename dirichlet0, typename dirichlet1, typename type, typename... types>
__device__ void Boundary(type* first, types*... second) {

	/*--------------------------Related Index--------------------------*/

	int i;
	int j = blockIdx.x;
	int k = threadIdx.x;
	int offset3d;
	type field_lr[4];

	/*-------------------------Inner Boundary-------------------------*/

	i = 0;
	offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	if constexpr (std::is_same_v<dirichlet0, trueType>) {
		first[offset3d] = 0;
	}
	else {
		field_lr[0] = first[offset3d + 1 * gridNz];
		field_lr[1] = first[offset3d + 2 * gridNz];
		field_lr[2] = first[offset3d + 3 * gridNz];
		field_lr[3] = first[offset3d + 4 * gridNz];
		first[offset3d] = (48 * field_lr[0] - 36 * field_lr[1] + 16 * field_lr[2] - 3 * field_lr[3]) / 25;
	}

	/*-------------------------Outer Boundary-------------------------*/

	i = gridNx - 1;
	offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	if constexpr (std::is_same_v<dirichlet1, trueType>) {
		first[offset3d] = 0;
	}
	else {
		field_lr[0] = first[offset3d - 4 * gridNz];
		field_lr[1] = first[offset3d - 3 * gridNz];
		field_lr[2] = first[offset3d - 2 * gridNz];
		field_lr[3] = first[offset3d - 1 * gridNz];
		first[offset3d] = (-3 * field_lr[0] + 16 * field_lr[1] - 36 * field_lr[2] + 48 * field_lr[3]) / 25;
	}

	if constexpr (sizeof...(second) > 0)
		Boundary<dirichlet0, dirichlet1>(second...);

}

template<typename dirichlet0, typename dirichlet1, typename... types>
__global__ void MHDBoundary(types* __restrict__ ... fields) {

	Boundary<dirichlet0, dirichlet1>(fields...);

}

template<typename local, typename type, typename... types>
__device__ void AlignedGhost(type* d_qtheta, type* first, types*... second) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = j * gridNx + i;
	int offset3d = j * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	int shift_k; type shift_lk; type shift_dk;
	type qtheta; type field; type field_du[4];

	qtheta = (d_qtheta[offset2d + gridNx] - d_qtheta[offset2d]) * gridNy;

	/*----------------------------Left Ghost----------------------------*/

	shift_lk = -qtheta / gridDz;
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

	}
	else {

		field_du[0] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
		field_du[1] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
		field_du[2] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
		field_du[3] = first[offset3d + gridNy * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

	}

	field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
		+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
		- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
		+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

	first[offset3d] = field;

	/*---------------------------Right Ghost---------------------------*/

	shift_lk = qtheta / gridDz;
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

	}
	else {

		field_du[0] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
		field_du[1] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
		field_du[2] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
		field_du[3] = first[offset3d + gridGhost * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

	}

	field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
		+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
		- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
		+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

	first[offset3d + (gridNy + gridGhost) * gridNxz] = field;

	if constexpr (sizeof...(second) > 0)
		AlignedGhost<local>(d_qtheta, second...);

}

template<typename local, typename type, typename... types>
__global__ void MHDAlignedGhost(type* __restrict__ d_qtheta, types* __restrict__ ... fields) {

	AlignedGhost<local>(d_qtheta, fields...);

}

/*-------------------------------------MHD Staggered2Collocated Shifted2Aligned-------------------------------------*/

template<typename local, typename type, typename... types>
__device__ void Staggered2C(type* d_qtheta, type* staggered, type* collocated, types*... fields) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = (j + gridGhost) * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	int shift_k; type shift_lk; type shift_dk;
	type qtheta; type qtheta_lr[4];
	type field_du[4]; type field_lr[4];

	/*-----------------------Staggered2Collocated-----------------------*/

	qtheta = d_qtheta[offset2d];

	qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
	qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
	qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
	qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

	qtheta_lr[0] = (qtheta_lr[0] + qtheta_lr[1]) / 2;
	qtheta_lr[3] = (qtheta_lr[2] + qtheta_lr[3]) / 2;
	qtheta_lr[1] = (qtheta_lr[1] + qtheta) / 2;
	qtheta_lr[2] = (qtheta_lr[2] + qtheta) / 2;

	if  constexpr (std::is_same_v<local, trueType>) {

		field_lr[0] = staggered[offset3d - 1 * gridNxz];
		field_lr[1] = staggered[offset3d + 0 * gridNxz];
		field_lr[2] = staggered[offset3d + 1 * gridNxz];
		field_lr[3] = staggered[offset3d + 2 * gridNxz];

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k - 1, gridNz);
			field_du[1] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 0, gridNz);
			field_du[2] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 1, gridNz);
			field_du[3] = __shfl_sync(0xffffffff, field_lr[index], lane_id + shift_k + 2, gridNz);

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}
	else {

#pragma unroll
		for (int index = 0; index < 4; index++) {

			shift_lk = (qtheta_lr[index] - qtheta) / gridDz;
			if constexpr (std::is_same_v<type, double>)
				shift_k = __double2int_rd(shift_lk);
			else
				shift_k = __float2int_rd(shift_lk);
			shift_dk = shift_lk - shift_k;

			field_du[0] = staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
			field_du[1] = staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
			field_du[2] = staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
			field_du[3] = staggered[offset3d + (index - 1) * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

			field_lr[index] = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
				+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
				- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
				+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

		}

	}

	collocated[offset3d] = ((field_lr[1] + field_lr[2]) * 9 - (field_lr[0] + field_lr[3])) / 16;

	if constexpr (sizeof...(fields) > 0)
		Staggered2C<local>(d_qtheta, fields...);

}

template<typename local, typename type, typename... types>
__global__ void MHDStaggered2C(type* __restrict__ d_qtheta, types* __restrict__ ... fields) {

	Staggered2C<local>(d_qtheta, fields...);

}

template<int dir, typename local, typename type, typename... types>
__device__ void Shifted2A(type* d_qtheta, type* shifted, type* aligned, types*... fields) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = (j + gridGhost) * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	int shift_k; type shift_lk; type shift_dk;
	type qtheta; type field; type field_du[4];

	/*----------------------Shifted2Aligned(dir=0)----------------------*/
	/*----------------------Aligned2Shifted(dir=1)----------------------*/

	qtheta = d_qtheta[offset2d];

	if constexpr (dir == 0)
		shift_lk = qtheta / gridDz;
	else
		shift_lk = -qtheta / gridDz;

	if constexpr (std::is_same_v<type, double>)
		shift_k = __double2int_rd(shift_lk);
	else
		shift_k = __float2int_rd(shift_lk);

	shift_dk = shift_lk - shift_k;

	if  constexpr (std::is_same_v<local, trueType>) {

		field = shifted[offset3d];

		field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
		field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
		field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
		field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

	}
	else {

		field_du[0] = shifted[offset3d - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
		field_du[1] = shifted[offset3d - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
		field_du[2] = shifted[offset3d - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
		field_du[3] = shifted[offset3d - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

	}

	field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
		+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
		- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
		+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

	aligned[offset3d] = field;

	if constexpr (sizeof...(fields) > 0)
		Shifted2A<dir, local>(d_qtheta, fields...);

}

template<int dir, typename local, typename type, typename... types>
__global__ void MHDShifted2A(type* __restrict__ d_qtheta, types* __restrict__ ... fields) {

	Shifted2A<dir, local>(d_qtheta, fields...);

}

/*--------------------------------------------------MHD NablaPara2, 4--------------------------------------------------*/

template<typename local, typename type>
__global__ void MHDNablaPara2(type* __restrict__ d_qtheta, type* __restrict__ d_field, type* __restrict__ d_nablaPara2) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = (j + gridGhost) * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	type qtheta;
	type qtheta_lr[4];
	int shift_k;
	type shift_lk;
	type shift_dk;

	/*-----------------------Field and Derivative-----------------------*/

	type field;
	type field_py;
	type field_du[4];
	type field_lr[4];

	type nablaPara2;

	/*-----------------------------Initialize-----------------------------*/

	qtheta = d_qtheta[offset2d];
	qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
	qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
	qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
	qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

	/*---------------------------NablaPara2---------------------------*/

	field = d_field[offset3d];

	PartialY2<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, d_field, field, field_du, field_lr, field_py, nablaPara2);

	d_nablaPara2[offset3d] = nablaPara2;

}

template<int F, typename local, typename type>
__global__ void MHDNablaPara4(type* __restrict__ d_qtheta, type* __restrict__ d_field, type* __restrict__ d_nablaPara2) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = (j + gridGhost) * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	type qtheta;
	type qtheta_lr[4];
	int shift_k;
	type shift_lk;
	type shift_dk;

	/*-----------------------Field and Derivative-----------------------*/

	type field;
	type field_py;
	type field_du[4];
	type field_lr[4];

	type para4;
	type nablaPara4;

	/*-----------------------------Initialize-----------------------------*/

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

	/*---------------------------NablaPara4---------------------------*/

	field = d_nablaPara2[offset3d];

	PartialY2<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, d_nablaPara2, field, field_du, field_lr, field_py, nablaPara4);

	d_field[offset3d] -= gridDt * para4 * nablaPara4;

}

/*----------------------------------------------MHD Filter Mode Number----------------------------------------------*/

template<typename cufftType>
__global__ void MHDFilterModeN(cufftType* __restrict__ d_freq, int modeNumber) {

	/*--------------------------Related Index--------------------------*/

	int offset = (blockIdx.x * blockDim.x + threadIdx.x) * nFFTFreqSize;

	/*-------------------------------Filter-------------------------------*/

	for (int mode = 0; mode < nFFTFreqSize; mode++) {
		if (mode != modeNumber) {

			d_freq[offset + mode].x = 0;
			d_freq[offset + mode].y = 0;

		}
	}

}

template<typename cufftType>
__global__ void MHDFilterModeN(cufftType* __restrict__ d_freq, int leftModeNumber, int rightModeNumber) {

	/*--------------------------Related Index--------------------------*/

	int offset = (blockIdx.x * blockDim.x + threadIdx.x) * nFFTFreqSize;

	/*-------------------------------Filter-------------------------------*/

	for (int mode = 0; mode < nFFTFreqSize; mode++) {
		if (mode < leftModeNumber || mode > rightModeNumber) {

			d_freq[offset + mode].x = 0;
			d_freq[offset + mode].y = 0;

		}
	}

}

template<typename type>
__global__ void MHDFilterResizeN(type* __restrict__ d_field) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	/*------------------------------Resize------------------------------*/

	d_field[offset3d] /= gridNz;

}

template<typename cufftType>
__global__ void MHDFilterModeM(cufftType* __restrict__ d_freq, int modeNumber) {

	/*--------------------------Related Index--------------------------*/

	int offset = (blockIdx.x * blockDim.x + threadIdx.x) * mFFTFreqSize;

	/*-------------------------------Filter-------------------------------*/

	for (int mode = 0; mode < mFFTFreqSize; mode++) {
		if (mode != modeNumber) {

			d_freq[offset + mode].x = 0;
			d_freq[offset + mode].y = 0;

		}
	}

}

template<typename cufftType>
__global__ void MHDFilterModeM(cufftType* __restrict__ d_freq, int leftModeNumber, int rightModeNumber) {

	/*--------------------------Related Index--------------------------*/

	int offset = (blockIdx.x * blockDim.x + threadIdx.x) * mFFTFreqSize;

	/*-------------------------------Filter-------------------------------*/

	for (int mode = 0; mode < mFFTFreqSize; mode++) {
		if (mode < leftModeNumber || mode > rightModeNumber) {

			d_freq[offset + mode].x = 0;
			d_freq[offset + mode].y = 0;

		}
	}

}

template<typename type>
__global__ void MHDFilterResizeM(type* __restrict__ d_field) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	/*------------------------------Resize------------------------------*/

	d_field[offset3d] /= gridNy;

}

template<typename type>
__global__ void MHDTransposeLeft(type* __restrict__ d_yxzField, type* __restrict__ d_xzyField) {

	/*--------------------------Related Index--------------------------*/

	int i = threadIdx.x;
	int j = blockIdx.x;

	/*-------------------------Transpose Left-------------------------*/

	for (int k = 0; k < gridNz; k++)
		d_xzyField[i * gridNz * gridNy + k * gridNy + j] = d_yxzField[j * gridNxz + i * gridNz + k];

}

template<typename type>
__global__ void MHDTransposeRight(type* __restrict__ d_xzyField, type* __restrict__ d_yxzField) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x;
	int k = threadIdx.x;

	/*-------------------------Transpose Right-------------------------*/

	for (int j = 0; j < gridNy; j++)
		d_yxzField[j * gridNxz + i * gridNz + k] = d_xzyField[i * gridNz * gridNy + k * gridNy + j];

}

template<typename type>
__global__ void MHDAddMode(type* __restrict__ d_Addend, type* __restrict__ d_Augend) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	/*-----------------------------Add N-----------------------------*/

	d_Augend[offset3d] += d_Addend[offset3d];

}

template<typename type>
__global__ void MHDSubtractMode(type* __restrict__ d_Subtrahend, type* __restrict__ d_Minuend) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	/*---------------------------Subtract N---------------------------*/

	d_Minuend[offset3d] -= d_Subtrahend[offset3d];

}

/*----------------------------------------------------MHD Diagnose----------------------------------------------------*/

template<typename cufftType, typename type>
__global__ void MHDDiagAmplitude(cufftType* __restrict__ d_freq, type* __restrict__  d_amplitude, type* __restrict__  d_modeReal, type* __restrict__  d_modeImag) {

	/*--------------------------Related Index--------------------------*/

	int offset = (diagY * gridNx + diagLeftX + threadIdx.x) * nFFTFreqSize;
	type real, imag;

	/*----------------------------Amplitude----------------------------*/

	for (int mode = leftN; mode <= rightN; mode++) {

		real = d_freq[offset + mode].x;
		imag = d_freq[offset + mode].y;

		if constexpr (std::is_same_v<type, double>)
			d_amplitude[threadIdx.x * (rightN - leftN + 1) + mode - leftN] = sqrt(real * real + imag * imag) / gridNz * 2;
		else
			d_amplitude[threadIdx.x * (rightN - leftN + 1) + mode - leftN] = sqrtf(real * real + imag * imag) / gridNz * 2;

		d_modeReal[threadIdx.x * (rightN - leftN + 1) + mode - leftN] = real / gridNz;
		d_modeImag[threadIdx.x * (rightN - leftN + 1) + mode - leftN] = imag / gridNz;

	}

}

template<typename type>
__global__ void MHDDiagFrequency(type* __restrict__ Phi_mid, type* __restrict__ d_frequency) {

	/*--------------------------Related Index--------------------------*/

	int offset3d = (diagY + gridGhost) * gridNxz + (diagLeftX + threadIdx.x) * gridNz;

	/*----------------------------Frequency----------------------------*/

	d_frequency[threadIdx.x] = Phi_mid[offset3d];

}

template<typename nonlinear, typename staggered, typename Eparallel, typename type>
__global__ void MHDDiagEparallel(type* __restrict__  d_qtheta,
	type* __restrict__ A_mid, type* __restrict__ dNe_mid, type* __restrict__ dTe_mid, type* __restrict__ Phi_mid,
	type* __restrict__ d_Ne0, type* __restrict__ d_Te0, type* __restrict__ d_Ne0_px, type* __restrict__ d_APhidNe2A,
	type* __restrict__ d_PhiA_A, type* __restrict__ d_NeA_A, type* __restrict__ d_Epara, type* __restrict__ d_EparaES) {

	/*--------------------------Related Index--------------------------*/

	int i = diagLeftX + blockIdx.x * blockDim.x + threadIdx.x;
	int j = diagY;
	int k = 0;
	int offset2d = (j + gridGhost) * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	type qtheta;
	type qtheta_lr[4];
	int shift_k;
	type shift_lk;
	type shift_dk;

	/*-----------------------Field and Derivative-----------------------*/

	type field;
	type field_px, field_py, field_pz;
	type field_du[4];
	type field_lr[4];

	/*---------------------Compressed Coefficient---------------------*/

	type compcoes[9];

	/*-----------------------------Diagnose-----------------------------*/

	type Epara, EparaES;

	/*-----------------------------Initialize-----------------------------*/

	qtheta = d_qtheta[offset2d];
	qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
	qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
	qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
	qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

	offset2d = j * gridNx + i;

	Epara = 0; EparaES = 0;

	/*-------Parallel Vector Potential in Parallel Vector Potential-------*/

	for (int index = 0; index < 5; index++)
		compcoes[index] = d_APhidNe2A[offset2d * 5 + index];

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		if constexpr (std::is_same_v<nonlinear, falseType>) {

			field = A_mid[offset3d];
			PartialZ<falseType>(k, offset3d, lane_id, A_mid, field, field_du, field_pz);
			PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du, field_lr, field_py);
			Epara += compcoes[0] * field + compcoes[1] * field_py + compcoes[2] * field_pz;

		}

	}

	/*-----------Electric Potential in Parallel Vector Potential-----------*/

	if constexpr (std::is_same_v<staggered, trueType>) {

		Collocated2S<falseType>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
			shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, field);
		field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * gridDy);

	}
	else {

		PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, field_py);

	}

	EparaES += compcoes[3] * field_py;

	/*-----------Electron Density in Parallel Vector Potential-----------*/

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		if constexpr (std::is_same_v<nonlinear, falseType>) {

			if constexpr (std::is_same_v<staggered, trueType>) {

				Collocated2S<falseType>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
					shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field);
				field_py = (field_lr[0] - 27 * field_lr[1] + 27 * field_lr[2] - field_lr[3]) / (24 * gridDy);

			}
			else {

				PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field_py);

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
		PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du, field_lr, A_py);

		if constexpr (std::is_same_v<staggered, falseType>) {

			Phi = Phi_mid[offset3d];
			PartialZ<falseType>(k, offset3d, lane_id, Phi_mid, Phi, field_du, Phi_pz);
			PartialX(0, i, k, offset3d, Phi_mid, Phi, field_lr, Phi_px);
			PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, Phi_py);

		}
		else {

			C2SPartialXYZ<falseType>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
				d_qtheta, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, Phi, Phi_px, Phi_py, Phi_pz);

		}

		/*-----------------PhiA in Parallel Vector Potential-----------------*/

		for (int index = 0; index < 9; index++)
			compcoes[index] = d_PhiA_A[offset2d * 9 + index];

		EparaES += compcoes[0] * Phi_px * A + compcoes[1] * Phi_px * A_py + compcoes[2] * Phi_px * A_pz
			+ compcoes[3] * Phi_py * A + compcoes[4] * Phi_py * A_px + compcoes[5] * Phi_py * A_pz
			+ compcoes[6] * Phi_pz * A + compcoes[7] * Phi_pz * A_px + compcoes[8] * Phi_pz * A_py;

		/*-----------------dNe in Parallel Vector Potential-----------------*/

		if constexpr (std::is_same_v<Eparallel, trueType>) {

			compcoes[0] = d_APhidNe2A[offset2d * 5 + 4];

			if constexpr (std::is_same_v<staggered, falseType>) {

				field = dNe_mid[offset3d];
				PartialZ<falseType>(k, offset3d, lane_id, dNe_mid, field, field_du, field_pz);
				PartialX(0, i, k, offset3d, dNe_mid, field, field_lr, field_px);
				PartialY<falseType>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field_py);

				Ne = Ne0 + field;

			}
			else {

				C2SPartialXYZ<falseType>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
					d_qtheta, qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field, field_px, field_py, field_pz);

				Collocated2S<falseType>(0, 0, i, j, k, offset2d, offset3d, lane_id, shift_k,
					shift_lk, shift_dk, d_qtheta, qtheta, qtheta_lr, dTe_mid, field_du, field_lr, Te);

				Ne = Ne0 + field;
				Te = Te0 + Te;

			}

			Epara += compcoes[0] * (Te / Te0 * Ne0 / Ne) * field_py;

		}

		/*-----------------NeA in Parallel Vector Potential-----------------*/

		if constexpr (std::is_same_v<Eparallel, trueType>) {

			for (int index = 0; index < 9; index++)
				compcoes[index] = d_NeA_A[offset2d * 9 + index];

			field_px += Ne0_px;

			Epara += (compcoes[0] * field_px * A + compcoes[1] * field_px * A_py + compcoes[2] * field_px * A_pz
				+ compcoes[3] * field_py * A + compcoes[4] * field_py * A_px + compcoes[5] * field_py * A_pz
				+ compcoes[6] * field_pz * A + compcoes[7] * field_pz * A_px + compcoes[8] * field_pz * A_py) * (Te / Te0 * Ne0 / Ne);

		}

	}

	d_Epara[blockIdx.x * blockDim.x + threadIdx.x] = -Epara;
	d_EparaES[blockIdx.x * blockDim.x + threadIdx.x] = EparaES;

}

/*----------------------------------------------------------PIC----------------------------------------------------------*/

template<typename local, typename type>
__global__ void PICAlignedGhost(type* __restrict__ d_qtheta, type* __restrict__ dP_mid) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = j * gridNx + i;
	int offset3d = j * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted------------------------------*/

	int shift_k; type shift_lk; type shift_dk;
	type qtheta; type field; type field_du[4];

	qtheta = (d_qtheta[offset2d + gridNx] - d_qtheta[offset2d]) * gridNy;

	/*----------------------------Left Ghost----------------------------*/

	shift_lk = qtheta / gridDz;
	if constexpr (std::is_same_v<type, double>)
		shift_k = __double2int_rd(shift_lk);
	else
		shift_k = __float2int_rd(shift_lk);
	shift_dk = shift_lk - shift_k;

	field = dP_mid[offset3d];

	if  constexpr (std::is_same_v<local, trueType>) {

		field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
		field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
		field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
		field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

	}
	else {

		field_du[0] = dP_mid[offset3d - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
		field_du[1] = dP_mid[offset3d - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
		field_du[2] = dP_mid[offset3d - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
		field_du[3] = dP_mid[offset3d - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

	}

	field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
		+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
		- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
		+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

	dP_mid[offset3d + gridNy * gridNxz] += field;

	/*---------------------------Right Ghost---------------------------*/

	shift_lk = -qtheta / gridDz;
	if constexpr (std::is_same_v<type, double>)
		shift_k = __double2int_rd(shift_lk);
	else
		shift_k = __float2int_rd(shift_lk);
	shift_dk = shift_lk - shift_k;

	field = dP_mid[offset3d + (gridNy + gridGhost) * gridNxz];

	if  constexpr (std::is_same_v<local, trueType>) {

		field_du[0] = __shfl_sync(0xffffffff, field, lane_id + shift_k - 1, gridNz);
		field_du[1] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 0, gridNz);
		field_du[2] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 1, gridNz);
		field_du[3] = __shfl_sync(0xffffffff, field, lane_id + shift_k + 2, gridNz);

	}
	else {

		field_du[0] = dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k - 1) % gridNz + gridNz) % gridNz];
		field_du[1] = dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k + 0) % gridNz + gridNz) % gridNz];
		field_du[2] = dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k + 1) % gridNz + gridNz) % gridNz];
		field_du[3] = dP_mid[offset3d + (gridNy + gridGhost) * gridNxz - k + ((k + shift_k + 2) % gridNz + gridNz) % gridNz];

	}

	field = -shift_dk * (-1 + shift_dk) * (-2 + shift_dk) / 6 * field_du[0]
		+ (1 + shift_dk) * (-1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[1]
		- shift_dk * (1 + shift_dk) * (-2 + shift_dk) / 2 * field_du[2]
		+ shift_dk * (1 + shift_dk) * (-1 + shift_dk) / 6 * field_du[3];

	dP_mid[offset3d + gridGhost * gridNxz] += field;

}

template<typename nonlinear, typename local, typename staggered, typename Eparallel, typename type>
__global__ void MHD2Apt(type* __restrict__  d_qtheta,
	type* __restrict__ A_mid, type* __restrict__ dNe_mid, type* __restrict__ dTe_mid, type* __restrict__ Phi_mid,
	type* __restrict__ d_Ne0, type* __restrict__ d_Te0, type* __restrict__ d_Ne0_px, type* __restrict__ d_APhidNe2A,
	type* __restrict__ d_PhiA_A, type* __restrict__ d_NeA_A, type* __restrict__ d_A_pt) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;
	int offset2d = (j + gridGhost) * gridNx + i;
	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;
	int lane_id = (threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x) % 32;

	/*------------------------------Shifted-----------------------------*/

	type qtheta;
	type qtheta_lr[4];
	int shift_k;
	type shift_lk;
	type shift_dk;

	/*-----------------------Field and Derivative-----------------------*/

	type field;
	type field_px, field_py, field_pz;
	type field_du[4];
	type field_lr[4];

	/*---------------------Compressed Coefficient---------------------*/

	type compcoes[9];

	/*--------------------------RK4 Variables-------------------------*/

	type dAdt;

	/*-----------------------------Initialize----------------------------*/

	qtheta = d_qtheta[offset2d];
	qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
	qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
	qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
	qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

	offset2d = j * gridNx + i;

	dAdt = 0;

	/*-----------Electric Potential in Parallel Vector Potential----------*/

	for (int index = 0; index < 5; index++)
		compcoes[index] = d_APhidNe2A[offset2d * 5 + index];

	PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, field_py);
	dAdt += compcoes[3] * field_py;

	/*-----------Electron Density in Parallel Vector Potential----------*/

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		if constexpr (std::is_same_v<nonlinear, falseType>) {

			PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field_py);
			dAdt += compcoes[4] * field_py;

		}

	}

	/*--------Parallel Vector Potential in Parallel Vector Potential------*/

	if constexpr (std::is_same_v<Eparallel, trueType>) {

		if constexpr (std::is_same_v<nonlinear, falseType>) {

			if constexpr (std::is_same_v<staggered, falseType>) {

				field = A_mid[offset3d];
				PartialZ<local>(k, offset3d, lane_id, A_mid, field, field_du, field_pz);
				PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du, field_lr, field_py);

			}
			else {

				S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
					d_qtheta, qtheta, qtheta_lr, A_mid, field_du, field_lr, field, field_px, field_py, field_pz);

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
		PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, Phi_mid, field_du, field_lr, Phi_py);

		if constexpr (std::is_same_v<staggered, falseType>) {

			A = A_mid[offset3d];
			PartialZ<local>(k, offset3d, lane_id, A_mid, A, field_du, A_pz);
			PartialX(0, i, k, offset3d, A_mid, A, field_lr, A_px);
			PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, A_mid, field_du, field_lr, A_py);

		}
		else {

			S2CPartialXYZ<local>(i, j, k, offset2d, offset3d, lane_id, shift_k, shift_lk, shift_dk,
				d_qtheta, qtheta, qtheta_lr, A_mid, field_du, field_lr, A, A_px, A_py, A_pz);

			offset2d = (j + gridGhost) * gridNx + i;

			qtheta = d_qtheta[offset2d];
			qtheta_lr[0] = d_qtheta[offset2d - 2 * gridNx];
			qtheta_lr[1] = d_qtheta[offset2d - 1 * gridNx];
			qtheta_lr[2] = d_qtheta[offset2d + 1 * gridNx];
			qtheta_lr[3] = d_qtheta[offset2d + 2 * gridNx];

			offset2d = j * gridNx + i;

		}

		/*-------------------------PhiA in Parallel Vector Potential---------------------*/

		for (int index = 0; index < 9; index++)
			compcoes[index] = d_PhiA_A[offset2d * 9 + index];

		dAdt += compcoes[0] * Phi_px * A + compcoes[1] * Phi_px * A_py + compcoes[2] * Phi_px * A_pz
			+ compcoes[3] * Phi_py * A + compcoes[4] * Phi_py * A_px + compcoes[5] * Phi_py * A_pz
			+ compcoes[6] * Phi_pz * A + compcoes[7] * Phi_pz * A_px + compcoes[8] * Phi_pz * A_py;

		/*-------------------------dNe in Parallel Vector Potential----------------------*/

		if constexpr (std::is_same_v<Eparallel, trueType>) {

			compcoes[0] = d_APhidNe2A[offset2d * 5 + 4];

			field = dNe_mid[offset3d];
			PartialZ<local>(k, offset3d, lane_id, dNe_mid, field, field_du, field_pz);
			PartialX(0, i, k, offset3d, dNe_mid, field, field_lr, field_px);
			PartialY<local>(k, offset3d, lane_id, shift_k, shift_lk, shift_dk, qtheta, qtheta_lr, dNe_mid, field_du, field_lr, field_py);

			Ne = Ne0 + field;

			dAdt += compcoes[0] * (Te / Te0 * Ne0 / Ne) * field_py;

		}

		/*-------------------------NeA in Parallel Vector Potential---------------------*/

		if constexpr (std::is_same_v<Eparallel, trueType>) {

			for (int index = 0; index < 9; index++)
				compcoes[index] = d_NeA_A[offset2d * 9 + index];

			field_px += Ne0_px;

			dAdt += (compcoes[0] * field_px * A + compcoes[1] * field_px * A_py + compcoes[2] * field_px * A_pz
				+ compcoes[3] * field_py * A + compcoes[4] * field_py * A_px + compcoes[5] * field_py * A_pz
				+ compcoes[6] * field_pz * A + compcoes[7] * field_pz * A_px + compcoes[8] * field_pz * A_py) * (Te / Te0 * Ne0 / Ne);

		}

	}

	d_A_pt[offset3d] = dAdt;

}

template<typename local, typename type>
__global__ void MHD2PIC(type* __restrict__ pic3d, type* __restrict__ globalA, type* __restrict__ globalPhi, type* __restrict__ globalApt) {

	/*--------------------------Related Index--------------------------*/

	int i = blockIdx.x * blockDim.z + threadIdx.z;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.x + threadIdx.x;

	int offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	/*-----------------------Field and Derivative-----------------------*/

	type field, field_px, field_py, field_pz;
	type field_du[4];

	/*----------------------------MHD2PIC----------------------------*/

	field = globalA[offset3d];

	PartialZ<local>(k, offset3d, offset3d, globalA, field, field_du, field_pz);
	PartialX(0, i, k, offset3d, globalA, field, field_du, field_px);

	field_du[0] = globalA[offset3d - 2 * gridNxz];
	field_du[1] = globalA[offset3d - 1 * gridNxz];
	field_du[2] = globalA[offset3d + 1 * gridNxz];
	field_du[3] = globalA[offset3d + 2 * gridNxz];
	field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * gridDy);

	offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

	pic3d[offset3d + 0] = field;
	pic3d[offset3d + 1] = field_px;
	pic3d[offset3d + 2] = field_py;
	pic3d[offset3d + 3] = field_pz;

	if (k < gridGhost) {

		offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

		pic3d[offset3d + 0] = field;
		pic3d[offset3d + 1] = field_px;
		pic3d[offset3d + 2] = field_py;
		pic3d[offset3d + 3] = field_pz;

	}
	else if (k > gridNz - gridGhost - 1) {

		offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

		pic3d[offset3d + 0] = field;
		pic3d[offset3d + 1] = field_px;
		pic3d[offset3d + 2] = field_py;
		pic3d[offset3d + 3] = field_pz;

	}

	offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	field = globalPhi[offset3d];

	PartialZ<local>(k, offset3d, offset3d, globalPhi, field, field_du, field_pz);
	PartialX(0, i, k, offset3d, globalPhi, field, field_du, field_px);

	field_du[0] = globalPhi[offset3d - 2 * gridNxz];
	field_du[1] = globalPhi[offset3d - 1 * gridNxz];
	field_du[2] = globalPhi[offset3d + 1 * gridNxz];
	field_du[3] = globalPhi[offset3d + 2 * gridNxz];
	field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * gridDy);

	offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

	pic3d[offset3d + 4] = field_px;
	pic3d[offset3d + 5] = field_py;
	pic3d[offset3d + 6] = field_pz;

	if (k < gridGhost) {

		offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

		pic3d[offset3d + 4] = field_px;
		pic3d[offset3d + 5] = field_py;
		pic3d[offset3d + 6] = field_pz;

	}
	else if (k > gridNz - gridGhost - 1) {

		offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

		pic3d[offset3d + 4] = field_px;
		pic3d[offset3d + 5] = field_py;
		pic3d[offset3d + 6] = field_pz;

	}

	offset3d = (j + gridGhost) * gridNxz + i * gridNz + k;

	field = globalApt[offset3d];

	offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

	pic3d[offset3d + 7] = field;

	if (k < gridGhost) {

		offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

		pic3d[offset3d + 7] = field;

	}
	else if (k > gridNz - gridGhost - 1) {

		offset3d = ((j + gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

		pic3d[offset3d + 7] = field;

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
			field_py = (-25 * field + 48 * field_du[0] - 36 * field_du[1] + 16 * field_du[2] - 3 * field_du[3]) / (12 * gridDy);
		}
		else if (j == 1) {
			field_du[0] = globalA[offset3d - 1 * gridNxz];
			field_du[1] = globalA[offset3d + 1 * gridNxz];
			field_du[2] = globalA[offset3d + 2 * gridNxz];
			field_du[3] = globalA[offset3d + 3 * gridNxz];
			field_py = (-3 * field_du[0] - 10 * field + 18 * field_du[1] - 6 * field_du[2] + field_du[3]) / (12 * gridDy);
		}
		else {
			field_du[0] = globalA[offset3d - 2 * gridNxz];
			field_du[1] = globalA[offset3d - 1 * gridNxz];
			field_du[2] = globalA[offset3d + 1 * gridNxz];
			field_du[3] = globalA[offset3d + 2 * gridNxz];
			field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * gridDy);
		}

		offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

		pic3d[offset3d + 0] = field;
		pic3d[offset3d + 1] = field_px;
		pic3d[offset3d + 2] = field_py;
		pic3d[offset3d + 3] = field_pz;

		if (k < gridGhost) {

			offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

			pic3d[offset3d + 0] = field;
			pic3d[offset3d + 1] = field_px;
			pic3d[offset3d + 2] = field_py;
			pic3d[offset3d + 3] = field_pz;

		}
		else if (k > gridNz - gridGhost - 1) {

			offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

			pic3d[offset3d + 0] = field;
			pic3d[offset3d + 1] = field_px;
			pic3d[offset3d + 2] = field_py;
			pic3d[offset3d + 3] = field_pz;

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
			field_py = (-25 * field + 48 * field_du[0] - 36 * field_du[1] + 16 * field_du[2] - 3 * field_du[3]) / (12 * gridDy);
		}
		else if (j == 1) {
			field_du[0] = globalPhi[offset3d - 1 * gridNxz];
			field_du[1] = globalPhi[offset3d + 1 * gridNxz];
			field_du[2] = globalPhi[offset3d + 2 * gridNxz];
			field_du[3] = globalPhi[offset3d + 3 * gridNxz];
			field_py = (-3 * field_du[0] - 10 * field + 18 * field_du[1] - 6 * field_du[2] + field_du[3]) / (12 * gridDy);
		}
		else {
			field_du[0] = globalPhi[offset3d - 2 * gridNxz];
			field_du[1] = globalPhi[offset3d - 1 * gridNxz];
			field_du[2] = globalPhi[offset3d + 1 * gridNxz];
			field_du[3] = globalPhi[offset3d + 2 * gridNxz];
			field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * gridDy);
		}

		offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

		pic3d[offset3d + 4] = field_px;
		pic3d[offset3d + 5] = field_py;
		pic3d[offset3d + 6] = field_pz;

		if (k < gridGhost) {

			offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

			pic3d[offset3d + 4] = field_px;
			pic3d[offset3d + 5] = field_py;
			pic3d[offset3d + 6] = field_pz;

		}
		else if (k > gridNz - gridGhost - 1) {

			offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

			pic3d[offset3d + 4] = field_px;
			pic3d[offset3d + 5] = field_py;
			pic3d[offset3d + 6] = field_pz;

		}

		offset3d = j * gridNxz + i * gridNz + k;

		field = globalApt[offset3d];

		offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

		pic3d[offset3d + 7] = field;

		if (k < gridGhost) {

			offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

			pic3d[offset3d + 7] = field;

		}
		else if (k > gridNz - gridGhost - 1) {

			offset3d = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

			pic3d[offset3d + 7] = field;

		}

	}
	else if (j > gridNy - gridGhost - 1) {

		offset3d = (j + 2 * gridGhost) * gridNxz + i * gridNz + k;

		field = globalA[offset3d];

		PartialZ<local>(k, offset3d, offset3d, globalA, field, field_du, field_pz);
		PartialX(0, i, k, offset3d, globalA, field, field_du, field_px);

		if (j == gridNy - 1) {
			field_du[0] = globalA[offset3d - 4 * gridNxz];
			field_du[1] = globalA[offset3d - 3 * gridNxz];
			field_du[2] = globalA[offset3d - 2 * gridNxz];
			field_du[3] = globalA[offset3d - 1 * gridNxz];
			field_py = (3 * field_du[0] - 16 * field_du[1] + 36 * field_du[2] - 48 * field_du[3] + 25 * field) / (12 * gridDy);
		}
		else if (j == gridNy - 2) {
			field_du[0] = globalA[offset3d - 3 * gridNxz];
			field_du[1] = globalA[offset3d - 2 * gridNxz];
			field_du[2] = globalA[offset3d - 1 * gridNxz];
			field_du[3] = globalA[offset3d + 1 * gridNxz];
			field_py = (-field_du[0] + 6 * field_du[1] - 18 * field_du[2] + 10 * field + 3 * field_du[3]) / (12 * gridDy);
		}
		else {
			field_du[0] = globalA[offset3d - 2 * gridNxz];
			field_du[1] = globalA[offset3d - 1 * gridNxz];
			field_du[2] = globalA[offset3d + 1 * gridNxz];
			field_du[3] = globalA[offset3d + 2 * gridNxz];
			field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * gridDy);
		}

		offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

		pic3d[offset3d + 0] = field;
		pic3d[offset3d + 1] = field_px;
		pic3d[offset3d + 2] = field_py;
		pic3d[offset3d + 3] = field_pz;

		if (k < gridGhost) {

			offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

			pic3d[offset3d + 0] = field;
			pic3d[offset3d + 1] = field_px;
			pic3d[offset3d + 2] = field_py;
			pic3d[offset3d + 3] = field_pz;

		}
		else if (k > gridNz - gridGhost - 1) {

			offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

			pic3d[offset3d + 0] = field;
			pic3d[offset3d + 1] = field_px;
			pic3d[offset3d + 2] = field_py;
			pic3d[offset3d + 3] = field_pz;

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
			field_py = (3 * field_du[0] - 16 * field_du[1] + 36 * field_du[2] - 48 * field_du[3] + 25 * field) / (12 * gridDy);
		}
		else if (j == gridNy - 2) {
			field_du[0] = globalPhi[offset3d - 3 * gridNxz];
			field_du[1] = globalPhi[offset3d - 2 * gridNxz];
			field_du[2] = globalPhi[offset3d - 1 * gridNxz];
			field_du[3] = globalPhi[offset3d + 1 * gridNxz];
			field_py = (-field_du[0] + 6 * field_du[1] - 18 * field_du[2] + 10 * field + 3 * field_du[3]) / (12 * gridDy);
		}
		else {
			field_du[0] = globalPhi[offset3d - 2 * gridNxz];
			field_du[1] = globalPhi[offset3d - 1 * gridNxz];
			field_du[2] = globalPhi[offset3d + 1 * gridNxz];
			field_du[3] = globalPhi[offset3d + 2 * gridNxz];
			field_py = (field_du[0] - 8 * field_du[1] + 8 * field_du[2] - field_du[3]) / (12 * gridDy);
		}

		offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

		pic3d[offset3d + 4] = field_px;
		pic3d[offset3d + 5] = field_py;
		pic3d[offset3d + 6] = field_pz;

		if (k < gridGhost) {

			offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

			pic3d[offset3d + 4] = field_px;
			pic3d[offset3d + 5] = field_py;
			pic3d[offset3d + 6] = field_pz;

		}
		else if (k > gridNz - gridGhost - 1) {

			offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

			pic3d[offset3d + 4] = field_px;
			pic3d[offset3d + 5] = field_py;
			pic3d[offset3d + 6] = field_pz;

		}

		offset3d = (j + 2 * gridGhost) * gridNxz + i * gridNz + k;

		field = globalApt[offset3d];

		offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost) * 8;

		pic3d[offset3d + 7] = field;

		if (k < gridGhost) {

			offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost + gridNz) * 8;

			pic3d[offset3d + 7] = field;

		}
		else if (k > gridNz - gridGhost - 1) {

			offset3d = ((j + 2 * gridGhost) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + gridGhost - gridNz) * 8;

			pic3d[offset3d + 7] = field;

		}

	}

}

template<int size, typename type, typename... types>
__device__ void FieldGather1d2d(int& address, type* coes, type* redundant, type& field, types&... fields) {

	field = 0;
#pragma unroll
	for (int index = 0; index < size; index++)
		field += redundant[address + index] * coes[index];

	if constexpr (sizeof...(fields) > 0) {
		address += size;
		FieldGather1d2d<size>(address, coes, redundant, fields...);
	}

}

template<typename type>
__device__ void FieldGather3d(int& i, int& j, int& k, int& offset, type* coes, type* redundant, type* fields) {

	offset = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[0] * redundant[offset + index];
	}

	offset = (j * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[1] * redundant[offset + index];
	}

	offset = ((j + 1) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[2] * redundant[offset + index];
	}

	offset = ((j + 1) * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[3] * redundant[offset + index];
	}

	offset = (j * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + 1) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[4] * redundant[offset + index];
	}

	offset = (j * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k + 1) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[5] * redundant[offset + index];
	}

	offset = ((j + 1) * gridNx * gridNzPlusGhost + i * gridNzPlusGhost + k + 1) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[6] * redundant[offset + index];
	}

	offset = ((j + 1) * gridNx * gridNzPlusGhost + (i + 1) * gridNzPlusGhost + k + 1) * 8;
#pragma unroll
	for (int index = 0; index < 8; index++) {
		fields[index] += coes[7] * redundant[offset + index];
	}

}

template<int ratioDt, picType particle, disType distribution, typename type>
__global__ void DriftAlignedRK4(type* __restrict__ pic1d, type* __restrict__ pic2d, type* __restrict__ pic3d, int* __restrict__ pic_keys_in, type* __restrict__ pic_values_in, type* __restrict__ dP_mid) {

	int illegal;
	int i, j, k;
	int qId, tileId, cellId, picId;

	type flag;
	type li, lj, lk;
	type coes[8] = {};

	type dx, dy, dz, dis, mu;
	type ddt[5] = {};
	type vec0[5] = {};
	type vec1[5] = {};
	type vec2[5] = {};

	type q, q_px, J, J_px, J_py, B, B_px, B_py;
	type gcovxy, gcovyy, gcovyz;
	type gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
	type APhiApt[8] = {};

	type bx, by, bz;
	type rho, bcony;
	type cx, cy, cz;
	type m2e, mu2e;
	type dxy, dxz, dyz;
	type dxB, dyB, dzB;
	type Bstarx, Bstary, Bstarz, Bstar;
	type na, na_px, nb, nb_px, ni, ni_px;
	type ta, ta_px, tb, tb_px, ti, ti_px;
	type V, E, cdwdt;

	type dvdt1, dxdt1, dydt1;
	type dxPhi, dyPhi, dzPhi;
	type cxdxA, cydyA, czdzA;

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

		li = (vec0[0] - xbeg) / gridDx;
		lj = (vec0[1] - ybeg) / gridDy;
		lk = (vec0[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			i = __double2int_rd(li);
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			i = __float2int_rd(li);
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dx = li - i;
		dy = lj - j;
		dz = lk - k;

#pragma unroll
		for (int index = 0; index < 5; index++) {
			vec1[index] = vec0[index];
			vec2[index] = vec0[index];
		}

		qId = i * 28;
		tileId = (j * cellNx + i) * 72;
		cellId *= 64;
		illegal = 0;

		auto dfdt_XVpara = [&]() {

			dxdt1 = 1 / Bstar * (vec1[3] * cxdxA - dxPhi);
			dydt1 = 1 / Bstar * (vec1[3] * cydyA - dyPhi);
			dvdt1 = -1 / m2e * (APhiApt[7] + 1 / Bstar * (Bstarx * APhiApt[4] + Bstary * APhiApt[5] + Bstarz * APhiApt[6] + mu2e * (cxdxA * B_px + cydyA * B_py)));

			if constexpr (particle == Ion) {

				if constexpr (distribution == Maxwell) {

					cdwdt = mp * va * va / (ti * kev);

					ddt[4] = (-ni_px / ni + 3 * ti_px / (2 * ti) + (mu * B_px - ti_px / ti * E) * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += IonMass * vec1[3] * cdwdt * dvdt1;

				}
				else {

					cdwdt = 3 * V / (2 * E * V + IonMass * ti * ti * ti);

					ddt[4] = (-ni_px / ni + mu * B_px * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += IonMass * vec1[3] * cdwdt * dvdt1;

					ddt[4] += 3 * ti * ti / (V * V * V + ti * ti * ti) * ti_px * dxdt1;

					if constexpr (distribution == Slowing0) {

						cdwdt = 0;

					}
					else if constexpr (distribution == Slowing1) {

						if constexpr (std::is_same_v<type, double>)
							cdwdt = 2 * exp(-pow((IonVb - V) / IonDeltaV, 2.0)) / (IonMass * V * IonDeltaV * sqrt(pi) * (1 + erf((IonVb - V) / IonDeltaV)));
						else
							cdwdt = 2 * expf(-powf((IonVb - V) / IonDeltaV, 2.0f)) / (IonMass * V * IonDeltaV * sqrtf(pi) * (1 + erff((IonVb - V) / IonDeltaV)));

					}
					else if constexpr (distribution == Slowing2) {

						cdwdt = 2 * mu * (IonLambda0 * E - mu) / (IonDeltaLambda2 * E * E * E);

					}
					else if constexpr (distribution == Slowing3) {

						cdwdt = 2 * mu * (IonLambda0 * E - mu) / (IonDeltaLambda2 * E * E * E);

						if constexpr (std::is_same_v<type, double>)
							cdwdt += 2 * exp(-pow((IonVb - V) / IonDeltaV, 2.0)) / (IonMass * V * IonDeltaV * sqrt(pi) * (1 + erf((IonVb - V) / IonDeltaV)));
						else
							cdwdt += 2 * expf(-powf((IonVb - V) / IonDeltaV, 2.0f)) / (IonMass * V * IonDeltaV * sqrtf(pi) * (1 + erff((IonVb - V) / IonDeltaV)));

					}

					ddt[4] += mu * B_px * cdwdt * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += IonMass * vec1[3] * cdwdt * dvdt1;

				}

			}
			else if constexpr (particle == Alpha) {

				if constexpr (distribution == Maxwell) {

					cdwdt = mp * va * va / (ta * kev);

					ddt[4] = (-na_px / na + 3 * ta_px / (2 * ta) + (mu * B_px - ta_px / ta * E) * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += AlphaMass * vec1[3] * cdwdt * dvdt1;

				}
				else {

					cdwdt = 3 * V / (2 * E * V + AlphaMass * ta * ta * ta);

					ddt[4] = (-na_px / na + mu * B_px * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += AlphaMass * vec1[3] * cdwdt * dvdt1;

					ddt[4] += 3 * ta * ta / (V * V * V + ta * ta * ta) * ta_px * dxdt1;

					if constexpr (distribution == Slowing0) {

						cdwdt = 0;

					}
					else if constexpr (distribution == Slowing1) {

						if constexpr (std::is_same_v<type, double>)
							cdwdt = 2 * exp(-pow((AlphaVb - V) / AlphaDeltaV, 2.0)) / (AlphaMass * V * AlphaDeltaV * sqrt(pi) * (1 + erf((AlphaVb - V) / AlphaDeltaV)));
						else
							cdwdt = 2 * expf(-powf((AlphaVb - V) / AlphaDeltaV, 2.0f)) / (AlphaMass * V * AlphaDeltaV * sqrtf(pi) * (1 + erff((AlphaVb - V) / AlphaDeltaV)));

					}
					else if constexpr (distribution == Slowing2) {

						cdwdt = 2 * mu * (AlphaLambda0 * E - mu) / (AlphaDeltaLambda2 * E * E * E);

					}
					else if constexpr (distribution == Slowing3) {

						cdwdt = 2 * mu * (AlphaLambda0 * E - mu) / (AlphaDeltaLambda2 * E * E * E);

						if constexpr (std::is_same_v<type, double>)
							cdwdt += 2 * exp(-pow((AlphaVb - V) / AlphaDeltaV, 2.0)) / (AlphaMass * V * AlphaDeltaV * sqrt(pi) * (1 + erf((AlphaVb - V) / AlphaDeltaV)));
						else
							cdwdt += 2 * expf(-powf((AlphaVb - V) / AlphaDeltaV, 2.0f)) / (AlphaMass * V * AlphaDeltaV * sqrtf(pi) * (1 + erff((AlphaVb - V) / AlphaDeltaV)));

					}

					ddt[4] += mu * B_px * cdwdt * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += AlphaMass * vec1[3] * cdwdt * dvdt1;

				}

			}
			else if constexpr (particle == Beam) {

				if constexpr (distribution == Maxwell) {

					cdwdt = mp * va * va / (tb * kev);

					ddt[4] = (-nb_px / nb + 3 * tb_px / (2 * tb) + (mu * B_px - tb_px / tb * E) * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += BeamMass * vec1[3] * cdwdt * dvdt1;

				}
				else {

					cdwdt = 3 * V / (2 * E * V + BeamMass * tb * tb * tb);

					ddt[4] = (-nb_px / nb + mu * B_px * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += BeamMass * vec1[3] * cdwdt * dvdt1;

					ddt[4] += 3 * tb * tb / (V * V * V + tb * tb * tb) * tb_px * dxdt1;

					if constexpr (distribution == Slowing0) {

						cdwdt = 0;

					}
					else if constexpr (distribution == Slowing1) {

						if constexpr (std::is_same_v<type, double>)
							cdwdt = 2 * exp(-pow((BeamVb - V) / BeamDeltaV, 2.0)) / (BeamMass * V * BeamDeltaV * sqrt(pi) * (1 + erf((BeamVb - V) / BeamDeltaV)));
						else
							cdwdt = 2 * expf(-powf((BeamVb - V) / BeamDeltaV, 2.0f)) / (BeamMass * V * BeamDeltaV * sqrtf(pi) * (1 + erff((BeamVb - V) / BeamDeltaV)));

					}
					else if constexpr (distribution == Slowing2) {

						cdwdt = 2 * mu * (BeamLambda0 * E - mu) / (BeamDeltaLambda2 * E * E * E);

					}
					else if constexpr (distribution == Slowing3) {

						cdwdt = 2 * mu * (BeamLambda0 * E - mu) / (BeamDeltaLambda2 * E * E * E);

						if constexpr (std::is_same_v<type, double>)
							cdwdt += 2 * exp(-pow((BeamVb - V) / BeamDeltaV, 2.0)) / (BeamMass * V * BeamDeltaV * sqrt(pi) * (1 + erf((BeamVb - V) / BeamDeltaV)));
						else
							cdwdt += 2 * expf(-powf((BeamVb - V) / BeamDeltaV, 2.0f)) / (BeamMass * V * BeamDeltaV * sqrtf(pi) * (1 + erff((BeamVb - V) / BeamDeltaV)));

					}

					ddt[4] += mu * B_px * cdwdt * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += BeamMass * vec1[3] * cdwdt * dvdt1;

				}

			}

			ddt[4] *= (dis - vec1[4]);

		};

		auto interpRK4 = [&]() {

#pragma unroll
			for (int index = 0; index < 2; index++)
				coes[index] = (hx[index] + sx[index] * dx);
			FieldGather1d2d<2>(qId, coes, pic1d, q, q_px, na, na_px, nb, nb_px, ni, ni_px, ta, ta_px, tb, tb_px, ti, ti_px);

#pragma unroll
			for (int index = 0; index < 2; index++)
				coes[index + 2] = coes[index];
#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[index] *= (hy[index] + sy[index] * dy);
			FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py,
				gcovxy, gcovyy, gcovyz, gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py);

#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[index + 4] = coes[index];
#pragma unroll
			for (int index = 0; index < 8; index++)
				coes[index] *= (hz[index] + sz[index] * dz);
#pragma unroll
			for (int index = 0; index < 8; index++)
				APhiApt[index] = 0;

			FieldGather3d(i, j, k, cellId, coes, pic3d, APhiApt);

			if constexpr (particle == Ion) {

				m2e = cm * IonMass / IonChar;
				mu2e = cm * mu / IonChar;
				E = IonMass * vec1[3] * vec1[3] / 2 + mu * B;
				if constexpr (std::is_same_v<type, double>)
					V = sqrt(2.0 * E / IonMass);
				else
					V = sqrtf(2.0f * E / IonMass);

			}
			else if constexpr (particle == Alpha) {

				m2e = cm * AlphaMass / AlphaChar;
				mu2e = cm * mu / AlphaChar;
				E = AlphaMass * vec1[3] * vec1[3] / 2 + mu * B;
				if constexpr (std::is_same_v<type, double>)
					V = sqrt(2.0 * E / AlphaMass);
				else
					V = sqrtf(2.0f * E / AlphaMass);

			}
			else if constexpr (particle == Beam) {

				m2e = cm * BeamMass / BeamChar;
				mu2e = cm * mu / BeamChar;
				E = BeamMass * vec1[3] * vec1[3] / 2 + mu * B;
				if constexpr (std::is_same_v<type, double>)
					V = sqrt(2.0 * E / BeamMass);
				else
					V = sqrtf(2.0f * E / BeamMass);

			}

			rho = rho0 + vec1[0] * drho;
			bcony = 2 * psitmax * drho * rho / (q * J * B);

			bx = bcony * gcovxy;
			by = bcony * gcovyy;
			bz = bcony * gcovyz;

			cx = bcony / J * (gcovyz_py - gcovyz * (J_py / J + B_py / B));
			cy = -bcony / J * (gcovyz_px - gcovyz * (J_px / J + B_px / B) + gcovyz * (drho / rho - q_px / q));
			cz = bcony / J * (gcovyy_px - gcovyy * (J_px / J + B_px / B) - gcovxy_py + gcovxy * (J_py / J + B_py / B) + gcovyy * (drho / rho - q_px / q));

			dxy = bcony / J * gcovyz;
			dxz = bcony / J * gcovyy;
			dyz = bcony / J * gcovxy;

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

			Bstarx += cxdxA;
			Bstary += cydyA;
			Bstarz += czdzA;
			Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;

			ddt[0] = 1 / Bstar * (vec1[3] * Bstarx - dxPhi - mu2e * dxB);
			ddt[1] = 1 / Bstar * (vec1[3] * Bstary - dyPhi - mu2e * dyB);
			ddt[2] = 1 / Bstar * (vec1[3] * Bstarz - dzPhi - mu2e * dzB);
			ddt[3] = -1 / m2e * (APhiApt[7] + 1 / Bstar * (Bstarx * (APhiApt[4] + mu2e * B_px) + Bstary * (APhiApt[5] + mu2e * B_py) + Bstarz * APhiApt[6]));

			dfdt_XVpara();

		};

		/*-----------------------------1st RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 6;

		flag = vec0[0] + ddt[0] * gridDt * ratioDt / 2;
		if (flag < xbeg || flag >= xend)
			illegal = 1;
		if (!illegal)
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index] + ddt[index] * gridDt * ratioDt / 2;
		else
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index];

		li = (vec1[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
		}
		else {
			vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
		}

		lj = (vec1[1] - ybeg) / gridDy;
		lk = (vec1[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		tileId = (j * cellNx + i) * 72;
		cellId = (j * cellNxz + i * cellNz + k) * 64;

		/*-----------------------------2nd RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 3;

		flag = vec0[0] + ddt[0] * gridDt * ratioDt / 2;
		if (flag < xbeg || flag >= xend)
			illegal = 1;
		if (!illegal)
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index] + ddt[index] * gridDt * ratioDt / 2;
		else
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index];

		li = (vec1[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
		}
		else {
			vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
		}

		lj = (vec1[1] - ybeg) / gridDy;
		lk = (vec1[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		tileId = (j * cellNx + i) * 72;
		cellId = (j * cellNxz + i * cellNz + k) * 64;

		/*-----------------------------3rd RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 3;

		flag = vec0[0] + ddt[0] * gridDt * ratioDt;
		if (flag < xbeg || flag >= xend)
			illegal = 1;
		if (!illegal)
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index] + ddt[index] * gridDt * ratioDt;
		else
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index];

		li = (vec1[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
		}
		else {
			vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
		}

		lj = (vec1[1] - ybeg) / gridDy;
		lk = (vec1[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		tileId = (j * cellNx + i) * 72;
		cellId = (j * cellNxz + i * cellNz + k) * 64;

		/*-----------------------------4th RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 6;

		if (vec2[0] < xbeg || vec2[0] >= xend)
			illegal = 1;
		if (illegal) {
			for (int index = 0; index < 5; index++)
				vec2[index] = vec0[index];
			vec2[1] = -vec0[1];
			vec2[4] = 0;
		}

		li = (vec2[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec2[2] = vec2[2] + q * floor((vec2[1] - yori) / yrange) * yrange;
			vec2[1] = vec2[1] - floor((vec2[1] - yori) / yrange) * yrange;
			vec2[2] = vec2[2] - floor((vec2[2] - zori) / zrange) * zrange;
		}
		else {
			vec2[2] = vec2[2] + q * floorf((vec2[1] - yori) / yrange) * yrange;
			vec2[1] = vec2[1] - floorf((vec2[1] - yori) / yrange) * yrange;
			vec2[2] = vec2[2] - floorf((vec2[2] - zori) / zrange) * zrange;
		}

		lj = (vec2[1] - ybeg) / gridDy;
		lk = (vec2[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		tileId = (j * cellNx + i) * 72;
		cellId = j * cellNxz + i * cellNz + k;

#pragma unroll
		for (int index = 0; index < 2; index++)
			coes[index + 2] = coes[index];
#pragma unroll
		for (int index = 0; index < 4; index++)
			coes[index] *= (hy[index] + sy[index] * dy);

		FieldGather1d2d<4>(tileId, coes, pic2d, J, B);

		if constexpr (particle == Ion)
			dis = (IonMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * IonConst / 2;
		else if constexpr (particle == Alpha)
			dis = (AlphaMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * AlphaConst / 2;
		else if constexpr (particle == Beam)
			dis = (BeamMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * BeamConst / 2;

#pragma unroll
		for (int index = 0; index < 4; index++)
			coes[index + 4] = coes[index];
#pragma unroll
		for (int index = 0; index < 8; index++)
			coes[index] *= (hz[index] + sz[index] * dz);

		if (i == 0) {
#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[2 * index] = 0;
		}
		else if (i == gridNx - 2) {
#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[2 * index + 1] = 0;
		}

		k = (k - gridGhost + gridNz) % gridNz;
		qId = j * gridNxz + i * gridNz + k;

		i = gridNz;
		j = gridNxz;
		k = (k + 1 + gridNz) % gridNz - k;

		atomicAdd(&dP_mid[qId], coes[0] * dis);
		atomicAdd(&dP_mid[qId + i], coes[1] * dis);
		atomicAdd(&dP_mid[qId + j], coes[2] * dis);
		atomicAdd(&dP_mid[qId + i + j], coes[3] * dis);
		atomicAdd(&dP_mid[qId + k], coes[4] * dis);
		atomicAdd(&dP_mid[qId + i + k], coes[5] * dis);
		atomicAdd(&dP_mid[qId + j + k], coes[6] * dis);
		atomicAdd(&dP_mid[qId + i + j + k], coes[7] * dis);

#pragma unroll
		for (int index = 0; index < 7; index++)
			pic_keys_in[picId + index * picDev] = cellId;
		pic_values_in[picId + 0 * picDev] = vec2[0];
		pic_values_in[picId + 1 * picDev] = vec2[1];
		pic_values_in[picId + 2 * picDev] = vec2[2];
		pic_values_in[picId + 3 * picDev] = vec2[3];
		pic_values_in[picId + 4 * picDev] = vec2[4];

	}

}

template<int ratioDt, int gyroNums, picType particle, disType distribution, typename type>
__global__ void GyroAlignedRK4(type* __restrict__ pic1d, type* __restrict__ pic2d, type* __restrict__ pic3d, int* __restrict__ pic_keys_in, type* __restrict__ pic_values_in, type* __restrict__ dP_mid) {

	int illegal;
	int i, j, k;
	int qId, tileId, cellId, picId;

	type flag;
	type li, lj, lk;
	type coes[8] = {};

	type dx, dy, dz, dis, mu;
	type ddt[5] = {};
	type vec0[5] = {};
	type vec1[5] = {};
	type vec2[5] = {};

	type q, q_px, J, J_px, J_py, B, B_px, B_py;
	type gcovxy, gcovyy, gcovyz;
	type gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
	type APhiApt[8] = {};

	type bx, by, bz;
	type rho, bcony;
	type cx, cy, cz;
	type m2e, mu2e;
	type dxy, dxz, dyz;
	type dxB, dyB, dzB;
	type Bstarx, Bstary, Bstarz, Bstar;
	type na, na_px, nb, nb_px, ni, ni_px;
	type ta, ta_px, tb, tb_px, ti, ti_px;
	type V, E, cdwdt;

	type gyroDx, gyroDy;
	type gyroX, gyroY, gyroZ;
	type gconxx, gconxy, gconyy;
	type R0, Z0, R1, Z1, angle, radius;
	type avecxdxA, avecydyA, aveczdzA;
	type avedxPhi, avedyPhi, avedzPhi;
	type avePhipx, avePhipy, avePhipz;
	type aveAptbx, aveAptby, aveAptbz;
	type dvdt1, dxdt1, dydt1;

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

		li = (vec0[0] - xbeg) / gridDx;
		lj = (vec0[1] - ybeg) / gridDy;
		lk = (vec0[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			i = __double2int_rd(li);
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			i = __float2int_rd(li);
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dx = li - i;
		dy = lj - j;
		dz = lk - k;

#pragma unroll
		for (int index = 0; index < 5; index++) {
			vec1[index] = vec0[index];
			vec2[index] = vec0[index];
		}

		qId = i * 28;
		tileId = (j * cellNx + i) * 72;
		cellId *= 64;
		illegal = 0;

		auto dfdt_XVpara = [&]() {

			dxdt1 = 1 / Bstar * (vec1[3] * avecxdxA - avedxPhi);
			dydt1 = 1 / Bstar * (vec1[3] * avecydyA - avedyPhi);
			dvdt1 = -1 / m2e / Bstar * (Bstarx * (avePhipx + aveAptbx) + Bstary * (avePhipy + aveAptby) + Bstarz * (avePhipz + aveAptbz) + mu2e * (avecxdxA * R0 + avecydyA * Z0));

			if constexpr (particle == Ion) {

				if constexpr (distribution == Maxwell) {

					cdwdt = mp * va * va / (ti * kev);

					ddt[4] = (-ni_px / ni + 3 * ti_px / (2 * ti) + (mu * B_px - ti_px / ti * E) * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += IonMass * vec1[3] * cdwdt * dvdt1;

				}
				else {

					cdwdt = 3 * V / (2 * E * V + IonMass * ti * ti * ti);

					ddt[4] = (-ni_px / ni + mu * B_px * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += IonMass * vec1[3] * cdwdt * dvdt1;

					ddt[4] += 3 * ti * ti / (V * V * V + ti * ti * ti) * ti_px * dxdt1;

					if constexpr (distribution == Slowing0) {

						cdwdt = 0;

					}
					else if constexpr (distribution == Slowing1) {

						if constexpr (std::is_same_v<type, double>)
							cdwdt = 2 * exp(-pow((IonVb - V) / IonDeltaV, 2.0)) / (IonMass * V * IonDeltaV * sqrt(pi) * (1 + erf((IonVb - V) / IonDeltaV)));
						else
							cdwdt = 2 * expf(-powf((IonVb - V) / IonDeltaV, 2.0f)) / (IonMass * V * IonDeltaV * sqrtf(pi) * (1 + erff((IonVb - V) / IonDeltaV)));

					}
					else if constexpr (distribution == Slowing2) {

						cdwdt = 2 * mu * (IonLambda0 * E - mu) / (IonDeltaLambda2 * E * E * E);

					}
					else if constexpr (distribution == Slowing3) {

						cdwdt = 2 * mu * (IonLambda0 * E - mu) / (IonDeltaLambda2 * E * E * E);

						if constexpr (std::is_same_v<type, double>)
							cdwdt += 2 * exp(-pow((IonVb - V) / IonDeltaV, 2.0)) / (IonMass * V * IonDeltaV * sqrt(pi) * (1 + erf((IonVb - V) / IonDeltaV)));
						else
							cdwdt += 2 * expf(-powf((IonVb - V) / IonDeltaV, 2.0f)) / (IonMass * V * IonDeltaV * sqrtf(pi) * (1 + erff((IonVb - V) / IonDeltaV)));

					}

					ddt[4] += mu * B_px * cdwdt * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += IonMass * vec1[3] * cdwdt * dvdt1;

				}

			}
			else if constexpr (particle == Alpha) {

				if constexpr (distribution == Maxwell) {

					cdwdt = mp * va * va / (ta * kev);

					ddt[4] = (-na_px / na + 3 * ta_px / (2 * ta) + (mu * B_px - ta_px / ta * E) * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += AlphaMass * vec1[3] * cdwdt * dvdt1;

				}
				else {

					cdwdt = 3 * V / (2 * E * V + AlphaMass * ta * ta * ta);

					ddt[4] = (-na_px / na + mu * B_px * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += AlphaMass * vec1[3] * cdwdt * dvdt1;

					ddt[4] += 3 * ta * ta / (V * V * V + ta * ta * ta) * ta_px * dxdt1;

					if constexpr (distribution == Slowing0) {

						cdwdt = 0;

					}
					else if constexpr (distribution == Slowing1) {

						if constexpr (std::is_same_v<type, double>)
							cdwdt = 2 * exp(-pow((AlphaVb - V) / AlphaDeltaV, 2.0)) / (AlphaMass * V * AlphaDeltaV * sqrt(pi) * (1 + erf((AlphaVb - V) / AlphaDeltaV)));
						else
							cdwdt = 2 * expf(-powf((AlphaVb - V) / AlphaDeltaV, 2.0f)) / (AlphaMass * V * AlphaDeltaV * sqrtf(pi) * (1 + erff((AlphaVb - V) / AlphaDeltaV)));

					}
					else if constexpr (distribution == Slowing2) {

						cdwdt = 2 * mu * (AlphaLambda0 * E - mu) / (AlphaDeltaLambda2 * E * E * E);

					}
					else if constexpr (distribution == Slowing3) {

						cdwdt = 2 * mu * (AlphaLambda0 * E - mu) / (AlphaDeltaLambda2 * E * E * E);

						if constexpr (std::is_same_v<type, double>)
							cdwdt += 2 * exp(-pow((AlphaVb - V) / AlphaDeltaV, 2.0)) / (AlphaMass * V * AlphaDeltaV * sqrt(pi) * (1 + erf((AlphaVb - V) / AlphaDeltaV)));
						else
							cdwdt += 2 * expf(-powf((AlphaVb - V) / AlphaDeltaV, 2.0f)) / (AlphaMass * V * AlphaDeltaV * sqrtf(pi) * (1 + erff((AlphaVb - V) / AlphaDeltaV)));

					}

					ddt[4] += mu * B_px * cdwdt * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += AlphaMass * vec1[3] * cdwdt * dvdt1;

				}

			}
			else if constexpr (particle == Beam) {

				if constexpr (distribution == Maxwell) {

					cdwdt = mp * va * va / (tb * kev);

					ddt[4] = (-nb_px / nb + 3 * tb_px / (2 * tb) + (mu * B_px - tb_px / tb * E) * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += BeamMass * vec1[3] * cdwdt * dvdt1;

				}
				else {

					cdwdt = 3 * V / (2 * E * V + BeamMass * tb * tb * tb);

					ddt[4] = (-nb_px / nb + mu * B_px * cdwdt) * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += BeamMass * vec1[3] * cdwdt * dvdt1;

					ddt[4] += 3 * tb * tb / (V * V * V + tb * tb * tb) * tb_px * dxdt1;

					if constexpr (distribution == Slowing0) {

						cdwdt = 0;

					}
					else if constexpr (distribution == Slowing1) {

						if constexpr (std::is_same_v<type, double>)
							cdwdt = 2 * exp(-pow((BeamVb - V) / BeamDeltaV, 2.0)) / (BeamMass * V * BeamDeltaV * sqrt(pi) * (1 + erf((BeamVb - V) / BeamDeltaV)));
						else
							cdwdt = 2 * expf(-powf((BeamVb - V) / BeamDeltaV, 2.0f)) / (BeamMass * V * BeamDeltaV * sqrtf(pi) * (1 + erff((BeamVb - V) / BeamDeltaV)));

					}
					else if constexpr (distribution == Slowing2) {

						cdwdt = 2 * mu * (BeamLambda0 * E - mu) / (BeamDeltaLambda2 * E * E * E);

					}
					else if constexpr (distribution == Slowing3) {

						cdwdt = 2 * mu * (BeamLambda0 * E - mu) / (BeamDeltaLambda2 * E * E * E);

						if constexpr (std::is_same_v<type, double>)
							cdwdt += 2 * exp(-pow((BeamVb - V) / BeamDeltaV, 2.0)) / (BeamMass * V * BeamDeltaV * sqrt(pi) * (1 + erf((BeamVb - V) / BeamDeltaV)));
						else
							cdwdt += 2 * expf(-powf((BeamVb - V) / BeamDeltaV, 2.0f)) / (BeamMass * V * BeamDeltaV * sqrtf(pi) * (1 + erff((BeamVb - V) / BeamDeltaV)));

					}

					ddt[4] += mu * B_px * cdwdt * dxdt1;
					ddt[4] += mu * B_py * cdwdt * dydt1;
					ddt[4] += BeamMass * vec1[3] * cdwdt * dvdt1;

				}

			}

			ddt[4] *= (dis - vec1[4]);

		};

		auto interpRK4 = [&]() {

#pragma unroll
			for (int index = 0; index < 2; index++)
				coes[index] = (hx[index] + sx[index] * dx);
			FieldGather1d2d<2>(qId, coes, pic1d, q, q_px, na, na_px, nb, nb_px, ni, ni_px, ta, ta_px, tb, tb_px, ti, ti_px);

#pragma unroll
			for (int index = 0; index < 2; index++)
				coes[index + 2] = coes[index];
#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[index] *= (hy[index] + sy[index] * dy);
			FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py,
				gcovxy, gcovyy, gcovyz, gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py,
				gconxx, gconxy, gconyy, R0, Z0);

			if constexpr (particle == Ion) {

				m2e = cm * IonMass / IonChar;
				mu2e = cm * mu / IonChar;
				E = IonMass * vec1[3] * vec1[3] / 2 + mu * B;
				if constexpr (std::is_same_v<type, double>) {
					radius = cm / IonChar * sqrt(2.0 * mu * IonMass / B);
					V = sqrt(2.0 * E / IonMass);
				}
				else {
					radius = cm / IonChar * sqrtf(2.0f * mu * IonMass / B);
					V = sqrtf(2.0f * E / IonMass);
				}

			}
			else if constexpr (particle == Alpha) {

				m2e = cm * AlphaMass / AlphaChar;
				mu2e = cm * mu / AlphaChar;
				E = AlphaMass * vec1[3] * vec1[3] / 2 + mu * B;
				if constexpr (std::is_same_v<type, double>) {
					radius = cm / AlphaChar * sqrt(2.0 * mu * AlphaMass / B);
					V = sqrt(2.0 * E / AlphaMass);
				}
				else {
					radius = cm / AlphaChar * sqrtf(2.0f * mu * AlphaMass / B);
					V = sqrtf(2.0f * E / AlphaMass);
				}

			}
			else if constexpr (particle == Beam) {

				m2e = cm * BeamMass / BeamChar;
				mu2e = cm * mu / BeamChar;
				E = BeamMass * vec1[3] * vec1[3] / 2 + mu * B;
				if constexpr (std::is_same_v<type, double>) {
					radius = cm / BeamChar * sqrt(2.0 * mu * BeamMass / B);
					V = sqrt(2.0 * E / BeamMass);
				}
				else {
					radius = cm / BeamChar * sqrtf(2.0f * mu * BeamMass / B);
					V = sqrtf(2.0f * E / BeamMass);
				}

			}

			rho = rho0 + vec1[0] * drho;
			bcony = 2 * psitmax * drho * rho / (q * J * B);

			bx = bcony * gcovxy;
			by = bcony * gcovyy;
			bz = bcony * gcovyz;

			cx = bcony / J * (gcovyz_py - gcovyz * (J_py / J + B_py / B));
			cy = -bcony / J * (gcovyz_px - gcovyz * (J_px / J + B_px / B) + gcovyz * (drho / rho - q_px / q));
			cz = bcony / J * (gcovyy_px - gcovyy * (J_px / J + B_px / B) - gcovxy_py + gcovxy * (J_py / J + B_py / B) + gcovyy * (drho / rho - q_px / q));

			dxy = bcony / J * gcovyz;
			dxz = bcony / J * gcovyy;
			dyz = bcony / J * gcovxy;

			dxB = dxy * B_py;
			dyB = -dxy * B_px;
			dzB = dxz * B_px - dyz * B_py;

			Bstarx = cx * m2e * vec1[3];
			Bstary = cy * m2e * vec1[3] + B * bcony;
			Bstarz = cz * m2e * vec1[3];

			Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;

			if constexpr (std::is_same_v<type, double>)
				angle = acos(gconxy / sqrt(gconxx * gconyy));
			else
				angle = acosf(gconxy / sqrtf(gconxx * gconyy));

			if (i == gridNx - 2) {
				tileId = (j * cellNx + i - 1) * 72 + 64;
				FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
				if constexpr (std::is_same_v<type, double>)
					gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
				else
					gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
			}
			else {
				tileId = (j * cellNx + i + 1) * 72 + 64;
				FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
				if constexpr (std::is_same_v<type, double>)
					gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
				else
					gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
			}

			tileId = ((j + 1) * cellNx + i) * 72 + 64;
			FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
			if constexpr (std::is_same_v<type, double>)
				gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDy;
			else
				gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDy;

			R0 = B_px; Z0 = B_py;
			avecxdxA = 0; avecydyA = 0; aveczdzA = 0;
			avedxPhi = 0;  avedyPhi = 0;  avedzPhi = 0;
			avePhipx = 0;  avePhipy = 0;  avePhipz = 0;
			aveAptbx = 0;  aveAptby = 0; aveAptbz = 0;

#pragma unroll
			for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

				if constexpr (std::is_same_v<type, double>) {
					gyroX = vec1[0] + sin(2.0 * pi * gyroId / gyroNums + angle / 2.0) / sin(angle) * gyroDx;
					gyroY = vec1[1] + sin(2.0 * pi * gyroId / gyroNums - angle / 2.0) / sin(angle) * gyroDy;
				}
				else {
					gyroX = vec1[0] + sinf(2.0f * pi * gyroId / gyroNums + angle / 2.0f) / sinf(angle) * gyroDx;
					gyroY = vec1[1] + sinf(2.0f * pi * gyroId / gyroNums - angle / 2.0f) / sinf(angle) * gyroDy;
				}

				if (gyroX < 0 || gyroX >= 1)
					continue;

				li = (gyroX - xbeg) / gridDx;
				if constexpr (std::is_same_v<type, double>)
					i = __double2int_rd(li);
				else
					i = __float2int_rd(li);
				dx = li - i;
				qId = i * 28;
				coes[0] = hx[0] + sx[0] * dx;
				coes[1] = hx[1] + sx[1] * dx;
				FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

				gyroZ = vec1[2] - q * (gyroY - vec1[1]) - vec1[1] * (gyroX - q) - (gyroY - vec1[1]) * (gyroX - q);

				if constexpr (std::is_same_v<type, double>) {
					gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
					gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
					gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
				}
				else {
					gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
					gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
					gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
				}

				lj = (gyroY - ybeg) / gridDy;
				lk = (gyroZ - zbeg) / gridDz;

				if constexpr (std::is_same_v<type, double>) {
					j = __double2int_rd(lj);
					k = __double2int_rd(lk);
				}
				else {
					j = __float2int_rd(lj);
					k = __float2int_rd(lk);
				}

				dy = lj - j;
				dz = lk - k;
				tileId = (j * cellNx + i) * 72;
				cellId = (j * cellNxz + i * cellNz + k) * 64;

#pragma unroll
				for (int index = 0; index < 2; index++)
					coes[index + 2] = coes[index];
#pragma unroll
				for (int index = 0; index < 4; index++)
					coes[index] *= (hy[index] + sy[index] * dy);
#pragma unroll
				for (int index = 0; index < 4; index++)
					coes[index + 4] = coes[index];
#pragma unroll
				for (int index = 0; index < 8; index++)
					coes[index] *= (hz[index] + sz[index] * dz);
#pragma unroll
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

			Bstarx += avecxdxA;
			Bstary += avecydyA;
			Bstarz += aveczdzA;
			Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;

			ddt[0] = 1 / Bstar * (vec1[3] * Bstarx - avedxPhi - mu2e * dxB);
			ddt[1] = 1 / Bstar * (vec1[3] * Bstary - avedyPhi - mu2e * dyB);
			ddt[2] = 1 / Bstar * (vec1[3] * Bstarz - avedzPhi - mu2e * dzB);
			ddt[3] = -1 / m2e / Bstar * (Bstarx * (avePhipx + aveAptbx + mu2e * R0) + Bstary * (avePhipy + aveAptby + mu2e * Z0) + Bstarz * (avePhipz + aveAptbz));

			dfdt_XVpara();

		};

		/*-----------------------------1st RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 6;

		flag = vec0[0] + ddt[0] * gridDt * ratioDt / 2;
		if (flag < xbeg || flag >= xend)
			illegal = 1;
		if (!illegal)
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index] + ddt[index] * gridDt * ratioDt / 2;
		else
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index];

		li = (vec1[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
		}
		else {
			vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
		}

		lj = (vec1[1] - ybeg) / gridDy;
		lk = (vec1[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		tileId = (j * cellNx + i) * 72;
		cellId = (j * cellNxz + i * cellNz + k) * 64;

		/*-----------------------------2nd RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 3;

		flag = vec0[0] + ddt[0] * gridDt * ratioDt / 2;
		if (flag < xbeg || flag >= xend)
			illegal = 1;
		if (!illegal)
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index] + ddt[index] * gridDt * ratioDt / 2;
		else
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index];

		li = (vec1[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
		}
		else {
			vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
		}

		lj = (vec1[1] - ybeg) / gridDy;
		lk = (vec1[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		tileId = (j * cellNx + i) * 72;
		cellId = (j * cellNxz + i * cellNz + k) * 64;

		/*-----------------------------3rd RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 3;

		flag = vec0[0] + ddt[0] * gridDt * ratioDt;
		if (flag < xbeg || flag >= xend)
			illegal = 1;
		if (!illegal)
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index] + ddt[index] * gridDt * ratioDt;
		else
			for (int index = 0; index < 5; index++)
				vec1[index] = vec0[index];

		li = (vec1[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec1[2] = vec1[2] + q * floor((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floor((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floor((vec1[2] - zori) / zrange) * zrange;
		}
		else {
			vec1[2] = vec1[2] + q * floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[1] = vec1[1] - floorf((vec1[1] - yori) / yrange) * yrange;
			vec1[2] = vec1[2] - floorf((vec1[2] - zori) / zrange) * zrange;
		}

		lj = (vec1[1] - ybeg) / gridDy;
		lk = (vec1[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		tileId = (j * cellNx + i) * 72;
		cellId = (j * cellNxz + i * cellNz + k) * 64;

		/*-----------------------------4th RK4----------------------------*/

		interpRK4();

		for (int index = 0; index < 5; index++)
			vec2[index] += ddt[index] * gridDt * ratioDt / 6;

		if (vec2[0] < xbeg || vec2[0] >= xend)
			illegal = 1;
		if (illegal) {
			for (int index = 0; index < 5; index++)
				vec2[index] = vec0[index];
			vec2[1] = -vec0[1];
			vec2[4] = 0;
		}

		li = (vec2[0] - xbeg) / gridDx;
		if constexpr (std::is_same_v<type, double>)
			i = __double2int_rd(li);
		else
			i = __float2int_rd(li);
		dx = li - i;
		qId = i * 28;
		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;
		FieldGather1d2d<2>(qId, coes, pic1d, q);

		if constexpr (std::is_same_v<type, double>) {
			vec2[2] = vec2[2] + q * floor((vec2[1] - yori) / yrange) * yrange;
			vec2[1] = vec2[1] - floor((vec2[1] - yori) / yrange) * yrange;
			vec2[2] = vec2[2] - floor((vec2[2] - zori) / zrange) * zrange;
		}
		else {
			vec2[2] = vec2[2] + q * floorf((vec2[1] - yori) / yrange) * yrange;
			vec2[1] = vec2[1] - floorf((vec2[1] - yori) / yrange) * yrange;
			vec2[2] = vec2[2] - floorf((vec2[2] - zori) / zrange) * zrange;
		}

		lj = (vec2[1] - ybeg) / gridDy;
		lk = (vec2[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dy = lj - j;
		dz = lk - k;
		cellId = j * cellNxz + i * cellNz + k;

#pragma unroll
		for (int index = 0; index < 2; index++)
			coes[index + 2] = coes[index];
#pragma unroll
		for (int index = 0; index < 4; index++)
			coes[index] *= (hy[index] + sy[index] * dy);

		tileId = (j * cellNx + i) * 72 + 4;
		FieldGather1d2d<4>(tileId, coes, pic2d, B);
		tileId = (j * cellNx + i) * 72 + 52;
		FieldGather1d2d<4>(tileId, coes, pic2d, gconxx, gconxy, gconyy, R0, Z0);

		if constexpr (particle == Ion) {
			if constexpr (std::is_same_v<type, double>)
				radius = cm / IonChar * sqrt(2.0 * mu * IonMass / B);
			else
				radius = cm / IonChar * sqrtf(2.0f * mu * IonMass / B);
		}
		else if constexpr (particle == Alpha) {
			if constexpr (std::is_same_v<type, double>)
				radius = cm / AlphaChar * sqrt(2.0 * mu * AlphaMass / B);
			else
				radius = cm / AlphaChar * sqrtf(2.0f * mu * AlphaMass / B);
		}
		else if constexpr (particle == Beam) {
			if constexpr (std::is_same_v<type, double>)
				radius = cm / BeamChar * sqrt(2.0 * mu * BeamMass / B);
			else
				radius = cm / BeamChar * sqrtf(2.0f * mu * BeamMass / B);
		}

		if constexpr (std::is_same_v<type, double>)
			angle = acos(gconxy / sqrt(gconxx * gconyy));
		else
			angle = acosf(gconxy / sqrtf(gconxx * gconyy));

		if (i == gridNx - 2) {
			tileId = (j * cellNx + i - 1) * 72 + 64;
			FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
			if constexpr (std::is_same_v<type, double>)
				gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
			else
				gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
		}
		else {
			tileId = (j * cellNx + i + 1) * 72 + 64;
			FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
			if constexpr (std::is_same_v<type, double>)
				gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
			else
				gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
		}

		tileId = ((j + 1) * cellNx + i) * 72 + 64;
		FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
		if constexpr (std::is_same_v<type, double>)
			gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDy;
		else
			gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDy;

#pragma unroll
		for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

			if constexpr (std::is_same_v<type, double>) {
				gyroX = vec2[0] + sin(2.0 * pi * gyroId / gyroNums + angle / 2.0) / sin(angle) * gyroDx;
				gyroY = vec2[1] + sin(2.0 * pi * gyroId / gyroNums - angle / 2.0) / sin(angle) * gyroDy;
			}
			else {
				gyroX = vec2[0] + sinf(2.0f * pi * gyroId / gyroNums + angle / 2.0f) / sinf(angle) * gyroDx;
				gyroY = vec2[1] + sinf(2.0f * pi * gyroId / gyroNums - angle / 2.0f) / sinf(angle) * gyroDy;
			}

			if (gyroX < 0 || gyroX >= 1)
				continue;

			li = (gyroX - xbeg) / gridDx;
			if constexpr (std::is_same_v<type, double>)
				i = __double2int_rd(li);
			else
				i = __float2int_rd(li);
			dx = li - i;
			qId = i * 28;
			coes[0] = hx[0] + sx[0] * dx;
			coes[1] = hx[1] + sx[1] * dx;
			FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

			gyroZ = vec2[2] - q * (gyroY - vec2[1]) - vec2[1] * (gyroX - q);

			if constexpr (std::is_same_v<type, double>) {
				gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
				gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
				gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
			}
			else {
				gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
				gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
				gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
			}

			lj = (gyroY - ybeg) / gridDy;
			lk = (gyroZ - zbeg) / gridDz;

			if constexpr (std::is_same_v<type, double>) {
				j = __double2int_rd(lj);
				k = __double2int_rd(lk);
			}
			else {
				j = __float2int_rd(lj);
				k = __float2int_rd(lk);
			}

			dy = lj - j;
			dz = lk - k;
			tileId = (j * cellNx + i) * 72;

#pragma unroll
			for (int index = 0; index < 2; index++)
				coes[index + 2] = coes[index];
#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[index] *= (hy[index] + sy[index] * dy);
			FieldGather1d2d<4>(tileId, coes, pic2d, J, B);

			if constexpr (particle == Ion)
				dis = (IonMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * IonConst / 2 / gyroNums;
			else if constexpr (particle == Alpha)
				dis = (AlphaMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * AlphaConst / 2 / gyroNums;
			else if constexpr (particle == Beam)
				dis = (BeamMass * vec2[3] * vec2[3] + mu * B) * vec2[4] / J * BeamConst / 2 / gyroNums;

#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[index + 4] = coes[index];
#pragma unroll
			for (int index = 0; index < 8; index++)
				coes[index] *= (hz[index] + sz[index] * dz);

			if (i == 0) {
#pragma unroll
				for (int index = 0; index < 4; index++)
					coes[2 * index] = 0;
			}
			else if (i == gridNx - 2) {
#pragma unroll
				for (int index = 0; index < 4; index++)
					coes[2 * index + 1] = 0;
			}

			k = (k - gridGhost + gridNz) % gridNz;
			qId = j * gridNxz + i * gridNz + k;

			i = gridNz;
			j = gridNxz;
			k = (k + 1 + gridNz) % gridNz - k;

			atomicAdd(&dP_mid[qId], coes[0] * dis);
			atomicAdd(&dP_mid[qId + i], coes[1] * dis);
			atomicAdd(&dP_mid[qId + j], coes[2] * dis);
			atomicAdd(&dP_mid[qId + i + j], coes[3] * dis);
			atomicAdd(&dP_mid[qId + k], coes[4] * dis);
			atomicAdd(&dP_mid[qId + i + k], coes[5] * dis);
			atomicAdd(&dP_mid[qId + j + k], coes[6] * dis);
			atomicAdd(&dP_mid[qId + i + j + k], coes[7] * dis);

		}

#pragma unroll
		for (int index = 0; index < 7; index++)
			pic_keys_in[picId + index * picDev] = cellId;
		pic_values_in[picId + 0 * picDev] = vec2[0];
		pic_values_in[picId + 1 * picDev] = vec2[1];
		pic_values_in[picId + 2 * picDev] = vec2[2];
		pic_values_in[picId + 3 * picDev] = vec2[3];
		pic_values_in[picId + 4 * picDev] = vec2[4];

	}

}

template<typename randState>
__global__ void PICSetupState(randState* __restrict__ randStates) {

	int stateId = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(stateId, 0, 0, &randStates[stateId]);

}

/*-----------------------------------------------------PIC Diagnose-----------------------------------------------------*/

template<picType particle, typename type>
__global__ void PICDiagDensity(type* __restrict__ pic2d, int* __restrict__ pic_keys_in, type* __restrict__ pic_values_in, type* __restrict__ pic_density) {

	int i, j;
	int tileId, picId;

	type J, li, lj, dx, dy, dis;
	type vec0[3] = {};
	type coes[4] = {};

	for (int id = 0; id < pptNums; id++) {

		picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
		vec0[0] = pic_values_in[picId + 0 * picDev];
		vec0[1] = pic_values_in[picId + 1 * picDev];
		vec0[2] = pic_values_in[picId + 4 * picDev];

		li = (vec0[0] - xbeg) / gridDx;
		lj = (vec0[1] - ybeg) / gridDy;

		if constexpr (std::is_same_v<type, double>) {
			i = __double2int_rd(li);
			j = __double2int_rd(lj);
		}
		else {
			i = __float2int_rd(li);
			j = __float2int_rd(lj);
		}

		dx = li - i;
		dy = lj - j;

		tileId = (j * cellNx + i) * 72;

		/*---------------------------Diag Density--------------------------*/

#pragma unroll
		for (int index = 0; index < 4; index++)
			coes[index] = (hx[index] + sx[index] * dx) * (hy[index] + sy[index] * dy);

		FieldGather1d2d<4>(tileId, coes, pic2d, J);

		if constexpr (particle == Ion)
			dis = vec0[2] / J * IonConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz);
		else if constexpr (particle == Alpha)
			dis = vec0[2] / J * AlphaConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz);
		else if constexpr (particle == Beam)
			dis = vec0[2] / J * BeamConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz);

		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;

		if (i == 0)
			coes[0] *= 2;
		else if (i == gridNx - 2)
			coes[1] *= 2;

		atomicAdd(&pic_density[i], coes[0] * dis);
		atomicAdd(&pic_density[i + 1], coes[1] * dis);

	}

}

template<int gyroNums, picType particle, typename type>
__global__ void PICDiagDiffusivity(type* __restrict__ pic1d, type* __restrict__ pic2d, type* __restrict__ pic3d, int* __restrict__ pic_keys_in, type* __restrict__ pic_values_in, type* __restrict__ pic_diffusivity) {

	int i, j, k;
	int qId, tileId, cellId, picId;

	type li, lj, lk;
	type coes[8] = {};

	type dx, dy, dz, dis, mu, dxdt;
	type vec0[5] = {};

	type q, q_px, J, J_px, J_py, B, B_px, B_py;
	type gcovxy, gcovyy, gcovyz;
	type gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py;
	type APhiApt[8] = {};

	type bx, by, bz;
	type rho, bcony;
	type cx, cy, cz;
	type m2e, mu2e;
	type dxy, dxz, dyz;
	type dxB, dyB, dzB;
	type Bstarx, Bstary, Bstarz, Bstar;
	type na, na_px, nb, nb_px, ni, ni_px;
	type ta, ta_px, tb, tb_px, ti, ti_px;

	type gyroDx, gyroDy;
	type gyroX, gyroY, gyroZ;
	type gconxx, gconxy, gconyy;
	type R0, Z0, R1, Z1, angle, radius;
	type avecxdxA, avecydyA, aveczdzA;
	type avedxPhi, avedyPhi, avedzPhi;
	type avePhipx, avePhipy, avePhipz;
	type aveAptbx, aveAptby, aveAptbz;

	for (int id = 0; id < pptNums; id++) {

		picId = blockIdx.x * blockDim.x * pptNums + id * blockDim.x + threadIdx.x;
		cellId = pic_keys_in[picId];
		vec0[0] = pic_values_in[picId + 0 * picDev];
		vec0[1] = pic_values_in[picId + 1 * picDev];
		vec0[2] = pic_values_in[picId + 2 * picDev];
		vec0[3] = pic_values_in[picId + 3 * picDev];
		vec0[4] = pic_values_in[picId + 4 * picDev];
		mu = pic_values_in[picId + 6 * picDev];

		li = (vec0[0] - xbeg) / gridDx;
		lj = (vec0[1] - ybeg) / gridDy;
		lk = (vec0[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			i = __double2int_rd(li);
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			i = __float2int_rd(li);
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dx = li - i;
		dy = lj - j;
		dz = lk - k;

		qId = i * 28;
		tileId = (j * cellNx + i) * 72;
		cellId *= 64;

		auto interpRK4 = [&]() {

#pragma unroll
			for (int index = 0; index < 2; index++)
				coes[index] = (hx[index] + sx[index] * dx);
			FieldGather1d2d<2>(qId, coes, pic1d, q, q_px, na, na_px, nb, nb_px, ni, ni_px, ta, ta_px, tb, tb_px, ti, ti_px);

#pragma unroll
			for (int index = 0; index < 2; index++)
				coes[index + 2] = coes[index];
#pragma unroll
			for (int index = 0; index < 4; index++)
				coes[index] *= (hy[index] + sy[index] * dy);
			FieldGather1d2d<4>(tileId, coes, pic2d, J, B, J_px, J_py, B_px, B_py,
				gcovxy, gcovyy, gcovyz, gcovxy_py, gcovyy_px, gcovyz_px, gcovyz_py,
				gconxx, gconxy, gconyy, R0, Z0);

			if constexpr (particle == Ion) {

				m2e = cm * IonMass / IonChar;
				mu2e = cm * mu / IonChar;
				if constexpr (std::is_same_v<type, double>)
					radius = cm / IonChar * sqrt(2.0 * mu * IonMass / B);
				else
					radius = cm / IonChar * sqrtf(2.0f * mu * IonMass / B);

			}
			else if constexpr (particle == Alpha) {

				m2e = cm * AlphaMass / AlphaChar;
				mu2e = cm * mu / AlphaChar;
				if constexpr (std::is_same_v<type, double>)
					radius = cm / AlphaChar * sqrt(2.0 * mu * AlphaMass / B);
				else
					radius = cm / AlphaChar * sqrtf(2.0f * mu * AlphaMass / B);

			}
			else if constexpr (particle == Beam) {

				m2e = cm * BeamMass / BeamChar;
				mu2e = cm * mu / BeamChar;
				if constexpr (std::is_same_v<type, double>)
					radius = cm / BeamChar * sqrt(2.0 * mu * BeamMass / B);
				else
					radius = cm / BeamChar * sqrtf(2.0f * mu * BeamMass / B);

			}

			rho = rho0 + vec0[0] * drho;
			bcony = 2 * psitmax * drho * rho / (q * J * B);

			bx = bcony * gcovxy;
			by = bcony * gcovyy;
			bz = bcony * gcovyz;

			cx = bcony / J * (gcovyz_py - gcovyz * (J_py / J + B_py / B));
			cy = -bcony / J * (gcovyz_px - gcovyz * (J_px / J + B_px / B) + gcovyz * (drho / rho - q_px / q));
			cz = bcony / J * (gcovyy_px - gcovyy * (J_px / J + B_px / B) - gcovxy_py + gcovxy * (J_py / J + B_py / B) + gcovyy * (drho / rho - q_px / q));

			dxy = bcony / J * gcovyz;
			dxz = bcony / J * gcovyy;
			dyz = bcony / J * gcovxy;

			dxB = dxy * B_py;
			dyB = -dxy * B_px;
			dzB = dxz * B_px - dyz * B_py;

			Bstarx = cx * m2e * vec0[3];
			Bstary = cy * m2e * vec0[3] + B * bcony;
			Bstarz = cz * m2e * vec0[3];

			Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;

			if constexpr (std::is_same_v<type, double>)
				angle = acos(gconxy / sqrt(gconxx * gconyy));
			else
				angle = acosf(gconxy / sqrtf(gconxx * gconyy));

			if (i == gridNx - 2) {
				tileId = (j * cellNx + i - 1) * 72 + 64;
				FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
				if constexpr (std::is_same_v<type, double>)
					gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
				else
					gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
			}
			else {
				tileId = (j * cellNx + i + 1) * 72 + 64;
				FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
				if constexpr (std::is_same_v<type, double>)
					gyroDx = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
				else
					gyroDx = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDx;
			}

			tileId = ((j + 1) * cellNx + i) * 72 + 64;
			FieldGather1d2d<4>(tileId, coes, pic2d, R1, Z1);
			if constexpr (std::is_same_v<type, double>)
				gyroDy = radius / sqrt((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDy;
			else
				gyroDy = radius / sqrtf((R1 - R0) * (R1 - R0) + (Z1 - Z0) * (Z1 - Z0)) * gridDy;

			R0 = B_px; Z0 = B_py;
			avecxdxA = 0; avecydyA = 0; aveczdzA = 0;
			avedxPhi = 0;  avedyPhi = 0;  avedzPhi = 0;
			avePhipx = 0;  avePhipy = 0;  avePhipz = 0;
			aveAptbx = 0;  aveAptby = 0; aveAptbz = 0;

#pragma unroll
			for (int gyroId = 0; gyroId < gyroNums; gyroId++) {

				if constexpr (std::is_same_v<type, double>) {
					gyroX = vec0[0] + sin(2.0 * pi * gyroId / gyroNums + angle / 2.0) / sin(angle) * gyroDx;
					gyroY = vec0[1] + sin(2.0 * pi * gyroId / gyroNums - angle / 2.0) / sin(angle) * gyroDy;
				}
				else {
					gyroX = vec0[0] + sinf(2.0f * pi * gyroId / gyroNums + angle / 2.0f) / sinf(angle) * gyroDx;
					gyroY = vec0[1] + sinf(2.0f * pi * gyroId / gyroNums - angle / 2.0f) / sinf(angle) * gyroDy;
				}

				if (gyroX < 0 || gyroX >= 1)
					continue;

				li = (gyroX - xbeg) / gridDx;
				if constexpr (std::is_same_v<type, double>)
					i = __double2int_rd(li);
				else
					i = __float2int_rd(li);
				dx = li - i;
				qId = i * 28;
				coes[0] = hx[0] + sx[0] * dx;
				coes[1] = hx[1] + sx[1] * dx;
				FieldGather1d2d<2>(qId, coes, pic1d, gyroX);

				gyroZ = vec0[2] - q * (gyroY - vec0[1]) - vec0[1] * (gyroX - q) - (gyroY - vec0[1]) * (gyroX - q);

				if constexpr (std::is_same_v<type, double>) {
					gyroZ = gyroZ + gyroX * floor((gyroY - yori) / yrange) * yrange;
					gyroY = gyroY - floor((gyroY - yori) / yrange) * yrange;
					gyroZ = gyroZ - floor((gyroZ - zori) / zrange) * zrange;
				}
				else {
					gyroZ = gyroZ + gyroX * floorf((gyroY - yori) / yrange) * yrange;
					gyroY = gyroY - floorf((gyroY - yori) / yrange) * yrange;
					gyroZ = gyroZ - floorf((gyroZ - zori) / zrange) * zrange;
				}

				lj = (gyroY - ybeg) / gridDy;
				lk = (gyroZ - zbeg) / gridDz;

				if constexpr (std::is_same_v<type, double>) {
					j = __double2int_rd(lj);
					k = __double2int_rd(lk);
				}
				else {
					j = __float2int_rd(lj);
					k = __float2int_rd(lk);
				}

				dy = lj - j;
				dz = lk - k;
				tileId = (j * cellNx + i) * 72;
				cellId = (j * cellNxz + i * cellNz + k) * 64;

#pragma unroll
				for (int index = 0; index < 2; index++)
					coes[index + 2] = coes[index];
#pragma unroll
				for (int index = 0; index < 4; index++)
					coes[index] *= (hy[index] + sy[index] * dy);
#pragma unroll
				for (int index = 0; index < 4; index++)
					coes[index + 4] = coes[index];
#pragma unroll
				for (int index = 0; index < 8; index++)
					coes[index] *= (hz[index] + sz[index] * dz);
#pragma unroll
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

			Bstarx += avecxdxA;
			Bstary += avecydyA;
			Bstarz += aveczdzA;
			Bstar = Bstarx * bx + Bstary * by + Bstarz * bz;

			dxdt = 1 / Bstar * (vec0[3] * Bstarx - avedxPhi - mu2e * dxB);

		};

		interpRK4();

		li = (vec0[0] - xbeg) / gridDx;
		lj = (vec0[1] - ybeg) / gridDy;
		lk = (vec0[2] - zbeg) / gridDz;

		if constexpr (std::is_same_v<type, double>) {
			i = __double2int_rd(li);
			j = __double2int_rd(lj);
			k = __double2int_rd(lk);
		}
		else {
			i = __float2int_rd(li);
			j = __float2int_rd(lj);
			k = __float2int_rd(lk);
		}

		dx = li - i;
		dy = lj - j;
		dz = lk - k;

		/*---------------------------Diag Diffusivity--------------------------*/

#pragma unroll
		for (int index = 0; index < 4; index++)
			coes[index] = (hx[index] + sx[index] * dx) * (hy[index] + sy[index] * dy);

		tileId = (j * cellNx + i) * 72;
		FieldGather1d2d<4>(tileId, coes, pic2d, J);
		tileId = (j * cellNx + i) * 72 + 52;
		FieldGather1d2d<4>(tileId, coes, pic2d, gconxx);

		if constexpr (std::is_same_v<type, double>)
			dxdt /= sqrt(gconxx);
		else
			dxdt /= sqrtf(gconxx);

		if constexpr (particle == Ion)
			dis = -vec0[4] / J * IonConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz) * dxdt * va / (ni_px * l3);
		else if constexpr (particle == Alpha)
			dis = -vec0[4] / J * AlphaConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz) * dxdt * va / (na_px * l3);
		else if constexpr (particle == Beam)
			dis = -vec0[4] / J * BeamConst * pitchB0 * pitchB0 / 2 / mu0 / (mp * va * va) / (gridNy * gridNz) * dxdt * va / (nb_px * l3);

		coes[0] = hx[0] + sx[0] * dx;
		coes[1] = hx[1] + sx[1] * dx;

		if (i == 0)
			coes[0] *= 2;
		else if (i == gridNx - 2)
			coes[1] *= 2;

		atomicAdd(&pic_diffusivity[i], coes[0] * dis);
		atomicAdd(&pic_diffusivity[i + 1], coes[1] * dis);

	}

}

/*----------------------------------------------------------------------MAIN----------------------------------------------------------------------*/

int main(int argc, char* argv[]) {

	std::string configPath = "../config.json";
	if (argc > 1) {
		configPath = argv[1];
	}
	try {
		loadConfigFromJson(configPath);
	}
	catch (const std::exception& e) {
		std::cout << "Config error: " << e.what() << std::endl;
		return 1;
	}
	h_PICGridDimx = h_picDev / pptNums / PICBlockDimx;
	const dataType dt = h_dt;
	const int continueSteps = h_continueSteps;
	const int totalSteps = h_totalSteps;
	const int diagSteps = h_diagSteps;
	const int outputSteps = h_outputSteps;
	const int outerLoopMax = h_outerLoopMax;
	const int innerLoopMax = h_innerLoopMax;
	const int MHDLoopMax = h_MHDLoopMax;
	const int ppcNums = h_ppcNums;
	const int picDev = h_picDev;
	const int PICGridDimx = h_PICGridDimx;
	const int leftN = h_leftN;
	const int rightN = h_rightN;
	const int diagLeftX = h_diagLeftX;
	const int diagRightX = h_diagRightX;
	const int diagY = h_diagY;

	int myRank, nRanks, localRank = 0;

	//initializing MPI

	MPICHECK(MPI_Init(&argc, &argv));
	MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
	MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));
	if (nRanks != hostNums) {
		std::cout << "Error: nRanks != hostNums." << std::endl;
		return 0;
	}

	//calculating localRank which is used in selecting a GPU

	uint64_t hostHashs[nRanks];
	char hostname[1024];
	getHostName(hostname, 1024);
	hostHashs[myRank] = getHostHash(hostname);
	MPICHECK(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hostHashs, sizeof(uint64_t), MPI_BYTE, MPI_COMM_WORLD));
	for (int p = 0; p < nRanks; p++) {
		if (p == myRank) break;
		if (hostHashs[p] == hostHashs[myRank]) localRank++;
	}

	//generating NCCL unique ID at one process and broadcasting it to all

	ncclUniqueId id;
	ncclComm_t comms[devNums];
	if (myRank == 0) ncclGetUniqueId(&id);
	MPICHECK(MPI_Bcast((void*)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));

	//initializing NCCL, group API is required around ncclCommInitRank as it is
	//called across multiple GPUs in each process

	NCCLCHECK(ncclGroupStart());
	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		syncConfigToDeviceConstants();
		NCCLCHECK(ncclCommInitRank(comms + i, nRanks * devNums, id, myRank * devNums + i));
	}
	NCCLCHECK(ncclGroupEnd());

	//initializing ncclSendRecv offset

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

	//initializing ncclType

	ncclDataType_t ncclType;
	if constexpr (std::is_same_v<dataType, double>)
		ncclType = ncclDouble;
	else
		ncclType = ncclFloat;

	//initializing cuGMEC

	HybridModel<dataType> cuGMEC(
		std::vector<int> { devNums, nRanks, myRank, localRank, gridNx, gridNy, gridNz, gridGhost, ppcNums, h_tubes },
		std::vector<double> { h_B0, h_L0, h_VA0, h_RHO0, h_RHO1, h_PSITMAX, dt },
		std::tuple<bool, double> { ifNablaPerp2A().value, h_perp2A },
		std::tuple<bool, double> { ifNablaPerp2Phi().value, h_perp2Phi },
		std::tuple<bool, double> { ifNablaPerp2dNe().value, h_perp2dNe },
		std::tuple<bool, double> { ifNablaPerp2dTe().value, h_perp2dTe },
		std::tuple<bool, double> { ifNablaPerp2dPi().value, h_perp2dPi },
		std::tuple<bool, double> { ifNablaPerp2dPa().value, h_perp2dPa },
		std::tuple<bool, double> { ifNablaPerp2dPb().value, h_perp2dPb },
		std::tuple<bool, unsigned int> { (ifIon().value&& ifIonSlowing().value) || (ifAlpha().value && ifAlphaSlowing().value) || (ifBeam().value && ifBeamSlowing().value), randMax },
		std::tuple<bool, std::vector<double>> { ifIon().value, std::vector<double>{h_IonMass, h_IonChar, h_IonBeta, h_IonVmin, h_IonVmax, h_IonVb, h_IonDeltaV, h_IonLambda0, h_IonDeltaLambda2} },
		std::tuple<bool, std::vector<double>> { ifAlpha().value, std::vector<double>{h_AlphaMass, h_AlphaChar, h_AlphaBeta, h_AlphaVmin, h_AlphaVmax, h_AlphaVb, h_AlphaDeltaV, h_AlphaLambda0, h_AlphaDeltaLambda2} },
		std::tuple<bool, std::vector<double>> { ifBeam().value, std::vector<double>{h_BeamMass, h_BeamChar, h_BeamBeta, h_BeamVmin, h_BeamVmax, h_BeamVb, h_BeamDeltaV, h_BeamLambda0, h_BeamDeltaLambda2} });

	cuGMEC.allocateHostMemory();
	cuGMEC.allocateDeviceMemory();

	cuGMEC.loadMHDEquilibrium(h_MHDCollocated);
	cuGMEC.loadMHDPerturbation(h_MHDPerturbation);
	cuGMEC.compressCollocatedCoefficient();

	cuGMEC.computeSparseMatrix<Laplacian, trueType, trueType>();

	if constexpr (std::is_same_v<ifNablaPerp2A, trueType>) {
		if constexpr (std::is_same_v<ifStaggered, falseType>)
			cuGMEC.computeSparseMatrix<Resistive, trueType, trueType>();
	}

	if constexpr (std::is_same_v<ifNablaPerp2Phi, trueType>)
		cuGMEC.computeSparseMatrix<Perp2Phi, trueType, trueType>();
	if constexpr (std::is_same_v<ifNablaPerp2dNe, trueType>)
		cuGMEC.computeSparseMatrix<Perp2dNe, trueType, trueType>();
	if constexpr (std::is_same_v<ifNablaPerp2dTe, trueType>)
		cuGMEC.computeSparseMatrix<Perp2dTe, trueType, trueType>();
	if constexpr (std::is_same_v<ifNablaPerp2dPi, trueType>)
		cuGMEC.computeSparseMatrix<Perp2dPi, trueType, trueType>();
	if constexpr (std::is_same_v<ifNablaPerp2dPa, trueType>)
		cuGMEC.computeSparseMatrix<Perp2dPa, trueType, trueType>();
	if constexpr (std::is_same_v<ifNablaPerp2dPb, trueType>)
		cuGMEC.computeSparseMatrix<Perp2dPb, trueType, trueType>();

	if constexpr (std::is_same_v<ifIon, trueType>) {
		if constexpr (std::is_same_v<ifContinue, trueType>) {
			std::string file1, file2, file3, file4;
			file1 = "../IonConst_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file2 = "../IonOffsets_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file3 = "../IonKeys_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file4 = "../IonValues_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			cuGMEC.loadParticles<Ion>(file1, file2, file3, file4);
			if (continueSteps == 0)
				cuGMEC.computeEquilibriumPressure<Ion>();
		}
		else {
			cuGMEC.loadParticles<Ion, IonType, IonMarker>();
			cuGMEC.computeEquilibriumPressure<Ion>();
		}
	}
	if constexpr (std::is_same_v<ifAlpha, trueType>) {
		if constexpr (std::is_same_v<ifContinue, trueType>) {
			std::string file1, file2, file3, file4;
			file1 = "../AlphaConst_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file2 = "../AlphaOffsets_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file3 = "../AlphaKeys_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file4 = "../AlphaValues_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			cuGMEC.loadParticles<Alpha>(file1, file2, file3, file4);
			if (continueSteps == 0)
				cuGMEC.computeEquilibriumPressure<Alpha>();
		}
		else {
			cuGMEC.loadParticles<Alpha, AlphaType, AlphaMarker>();
			cuGMEC.computeEquilibriumPressure<Alpha>();
		}
	}
	if constexpr (std::is_same_v<ifBeam, trueType>) {
		if constexpr (std::is_same_v<ifContinue, trueType>) {
			std::string file1, file2, file3, file4;
			file1 = "../BeamConst_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file2 = "../BeamOffsets_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file3 = "../BeamKeys_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			file4 = "../BeamValues_" + std::to_string(myRank) + "_" + std::to_string(continueSteps) + ".bin";
			cuGMEC.loadParticles<Beam>(file1, file2, file3, file4);
			if (continueSteps == 0)
				cuGMEC.computeEquilibriumPressure<Beam>();
		}
		else {
			cuGMEC.loadParticles<Beam, BeamType, BeamMarker>();
			cuGMEC.computeEquilibriumPressure<Beam>();
		}
	}
	if constexpr ((std::is_same_v<ifIon, trueType> && std::is_same_v<ifIonSlowing, trueType>)
		|| (std::is_same_v<ifAlpha, trueType> && std::is_same_v<ifAlphaSlowing, trueType>)
		|| (std::is_same_v<ifBeam, trueType> && std::is_same_v<ifBeamSlowing, trueType>)) {
		cuGMEC.loadRandom();
	}

	cuGMEC.memcpyHostToDevice();

	if constexpr (std::is_same_v<ifStaggered, trueType>) {
		cuGMEC.loadMHDEquilibrium(h_MHDStaggered);
		cuGMEC.compressStaggeredCoefficient();
		if constexpr (std::is_same_v<ifNablaPerp2A, trueType>)
			cuGMEC.computeSparseMatrix<Resistive, trueType, trueType>();
	}

	if (myRank == 0)
		if ((access(std::string("./result").c_str(), 0)) == -1)
			mkdir(std::string("./result").c_str(), S_IRWXU);

	//initializing cuFFT for filtering toroidal mode

	std::vector<cufftHandle> nPlanR2Cs(devNums);
	std::vector<cufftHandle> nPlanC2Rs(devNums);
	cufftDoubleComplex* nFreqd[devNums];
	cufftComplex* nFreqf[devNums];

	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		if constexpr (std::is_same_v<dataType, double>) {
			CUFFTCHECK(cufftPlan1d(&nPlanR2Cs[i], nFFTTimeSize, CUFFT_D2Z, nFFTBatchSize));
			CUFFTCHECK(cufftPlan1d(&nPlanC2Rs[i], nFFTTimeSize, CUFFT_Z2D, nFFTBatchSize));
			CUDACHECK(cudaMalloc((void**)&nFreqd[i], sizeof(cufftDoubleComplex) * nFFTBatchSize * nFFTFreqSize));
		}
		else {
			CUFFTCHECK(cufftPlan1d(&nPlanR2Cs[i], nFFTTimeSize, CUFFT_R2C, nFFTBatchSize));
			CUFFTCHECK(cufftPlan1d(&nPlanC2Rs[i], nFFTTimeSize, CUFFT_C2R, nFFTBatchSize));
			CUDACHECK(cudaMalloc((void**)&nFreqf[i], sizeof(cufftComplex) * nFFTBatchSize * nFFTFreqSize));
		}
		CUFFTCHECK(cufftSetStream(nPlanR2Cs[i], 0));
		CUFFTCHECK(cufftSetStream(nPlanC2Rs[i], 0));
	}

	//initializing cuFFT for filtering poloidal mode

	std::vector<cufftHandle> mPlanR2Cs(devNums);
	std::vector<cufftHandle> mPlanC2Rs(devNums);
	cufftDoubleComplex* mFreqd[devNums];
	cufftComplex* mFreqf[devNums];

	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		if constexpr (std::is_same_v<dataType, double>) {
			CUFFTCHECK(cufftPlan1d(&mPlanR2Cs[i], mFFTTimeSize, CUFFT_D2Z, mFFTBatchSize));
			CUFFTCHECK(cufftPlan1d(&mPlanC2Rs[i], mFFTTimeSize, CUFFT_Z2D, mFFTBatchSize));
			CUDACHECK(cudaMalloc((void**)&mFreqd[i], sizeof(cufftDoubleComplex) * mFFTBatchSize * mFFTFreqSize));
		}
		else {
			CUFFTCHECK(cufftPlan1d(&mPlanR2Cs[i], mFFTTimeSize, CUFFT_R2C, mFFTBatchSize));
			CUFFTCHECK(cufftPlan1d(&mPlanC2Rs[i], mFFTTimeSize, CUFFT_C2R, mFFTBatchSize));
			CUDACHECK(cudaMalloc((void**)&mFreqf[i], sizeof(cufftComplex) * mFFTBatchSize * mFFTFreqSize));
		}
		CUFFTCHECK(cufftSetStream(mPlanR2Cs[i], 0));
		CUFFTCHECK(cufftSetStream(mPlanC2Rs[i], 0));
	}

	//initializing cuRAND

	randomState* randStates[devNums];

	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		CUDACHECK(cudaMalloc((void**)&randStates[i], sizeof(randomState) * PICGridDimx * PICBlockDimx));
		PICSetupState << <PICGridDimx, PICBlockDimx >> > (randStates[i]);
	}

	dim3 MRK4GridSize(MRK4GridDimx, MRK4GridDimy, MRK4GridDimz);
	dim3 MRK4BlockSize(MRK4BlockDimx, MRK4BlockDimy, MRK4BlockDimz);

	dim3 GhostGridSize(GhostGridDimx, GhostGridDimy, GhostGridDimz);
	dim3 GhostBlockSize(GhostBlockDimx, GhostBlockDimy, GhostBlockDimz);

	dim3 M2PGridSize(M2PGridDimx, M2PGridDimy, M2PGridDimz);
	dim3 M2PBlockSize(M2PBlockDimx, M2PBlockDimy, M2PBlockDimz);

	dim3 PICGridSize(PICGridDimx);
	dim3 PICBlockSize(PICBlockDimx);

	//referencing MHD

	dataType**& qtheta = cuGMEC.d_qtheta;

	dataType**& w_beg = cuGMEC.d_w_beg;
	dataType**& w_midl = cuGMEC.d_w_midl;
	dataType**& w_midr = cuGMEC.d_w_midr;
	dataType**& w_end = cuGMEC.d_w_end;

	dataType**& A_beg = cuGMEC.d_A_beg;
	dataType**& A_midl = cuGMEC.d_A_midl;
	dataType**& A_midr = cuGMEC.d_A_midr;
	dataType**& A_end = cuGMEC.d_A_end;

	dataType**& dNe_beg = cuGMEC.d_dNe_beg;
	dataType**& dNe_midl = cuGMEC.d_dNe_midl;
	dataType**& dNe_midr = cuGMEC.d_dNe_midr;
	dataType**& dNe_end = cuGMEC.d_dNe_end;

	dataType**& dTe_beg = cuGMEC.d_dTe_beg;
	dataType**& dTe_midl = cuGMEC.d_dTe_midl;
	dataType**& dTe_midr = cuGMEC.d_dTe_midr;
	dataType**& dTe_end = cuGMEC.d_dTe_end;

	dataType**& Phi_midl = cuGMEC.d_Phi_midl;
	dataType**& Phi_midr = cuGMEC.d_Phi_midr;
	dataType**& dJpB_midl = cuGMEC.d_dJpB_midl;
	dataType**& dJpB_midr = cuGMEC.d_dJpB_midr;
	dataType**& dPe_midl = cuGMEC.d_dPe_midl;
	dataType**& dPe_midr = cuGMEC.d_dPe_midr;
	dataType**& Apt_midl = cuGMEC.d_Apt_midl;
	dataType**& Apt_midr = cuGMEC.d_Apt_midr;

	dataType**& Ne0 = cuGMEC.d_Ne0;
	dataType**& Te0 = cuGMEC.d_Te0;
	dataType**& Ne0_px = cuGMEC.d_Ne0_px;
	dataType**& Te0_px = cuGMEC.d_Te0_px;

	dataType**& w2Phi = cuGMEC.d_w_Phi;
	dataType**& A2dJpB = cuGMEC.d_A2dJpB;
	dataType**& Phi2w = cuGMEC.d_Phi2w;
	dataType**& wdPAdJpB2w = cuGMEC.d_wdPAdJpB2w;
	dataType**& APhidNe2A = cuGMEC.d_APhidNe2A;
	dataType**& dPePhiAdJpB2dNe = cuGMEC.d_dPePhiAdJpB2dNe;
	dataType**& PhidTedNe2dTe = cuGMEC.d_PhidTedNe2dTe;

	dataType**& wPhi_w = cuGMEC.d_wPhi_w;
	dataType**& AdJpB_w = cuGMEC.d_AdJpB_w;
	dataType**& PhiA_A = cuGMEC.d_PhiA_A;
	dataType**& NeA_A = cuGMEC.d_NeA_A;
	dataType**& dNePhi_dNe = cuGMEC.d_dNePhi_dNe;
	dataType**& PhiTe_dTe = cuGMEC.d_PhiTe_dTe;
	dataType**& PhiTeA_dTe = cuGMEC.d_PhiTeA_dTe;

	//referencing PIC

	dataType**& pic_APhidNe2A = cuGMEC.dpic_APhidNe2A;
	dataType**& pic_PhiA_A = cuGMEC.dpic_PhiA_A;
	dataType**& pic_NeA_A = cuGMEC.dpic_NeA_A;

	dataType**& pic1d = cuGMEC.d_pic1d;
	dataType**& pic2d = cuGMEC.d_pic2d;
	dataType**& pic3d = cuGMEC.d_pic3d;
	dataType**& globalA = cuGMEC.d_globalA;
	dataType**& globalPhi = cuGMEC.d_globalPhi;
	dataType**& globalApt = cuGMEC.d_globalApt;
	dataType**& globalPa = cuGMEC.d_globalPa;
	dataType**& globalPi = cuGMEC.d_globalPi;
	dataType**& globalPb = cuGMEC.d_globalPb;
	dataType**& dPa_midl = cuGMEC.d_dPa_midl;
	dataType**& dPa_midr = cuGMEC.d_dPa_midr;
	dataType**& dPi_midl = cuGMEC.d_dPi_midl;
	dataType**& dPi_midr = cuGMEC.d_dPi_midr;
	dataType**& dPb_midl = cuGMEC.d_dPb_midl;
	dataType**& dPb_midr = cuGMEC.d_dPb_midr;

	int**& Alpha_offsets = cuGMEC.d_Alpha_offsets;
	int**& Alpha_keys_in = cuGMEC.d_Alpha_keys_in;
	int**& Alpha_keys_out = cuGMEC.d_Alpha_keys_out;
	dataType**& Alpha_values_in = cuGMEC.d_Alpha_values_in;
	dataType**& Alpha_values_out = cuGMEC.d_Alpha_values_out;

	int**& Ion_offsets = cuGMEC.d_Ion_offsets;
	int**& Ion_keys_in = cuGMEC.d_Ion_keys_in;
	int**& Ion_keys_out = cuGMEC.d_Ion_keys_out;
	dataType**& Ion_values_in = cuGMEC.d_Ion_values_in;
	dataType**& Ion_values_out = cuGMEC.d_Ion_values_out;

	int**& Beam_offsets = cuGMEC.d_Beam_offsets;
	int**& Beam_keys_in = cuGMEC.d_Beam_keys_in;
	int**& Beam_keys_out = cuGMEC.d_Beam_keys_out;
	dataType**& Beam_values_in = cuGMEC.d_Beam_values_in;
	dataType**& Beam_values_out = cuGMEC.d_Beam_values_out;

	int**& rand_keys_in = cuGMEC.d_rand_keys;
	dataType**& rand_values_in = cuGMEC.d_rand_values;

	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		if constexpr (std::is_same_v<ifIon, trueType>)
			CUDACHECK(cudaMemcpyToSymbol(IonConst, &cuGMEC.IonConst, sizeof(dataType)));
		if constexpr (std::is_same_v<ifAlpha, trueType>)
			CUDACHECK(cudaMemcpyToSymbol(AlphaConst, &cuGMEC.AlphaConst, sizeof(dataType)));
		if constexpr (std::is_same_v<ifBeam, trueType>)
			CUDACHECK(cudaMemcpyToSymbol(BeamConst, &cuGMEC.BeamConst, sizeof(dataType)));
	}

	std::vector<dataType> h_amplitude((totalSteps / diagSteps + 1) * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1));
	std::vector<dataType> h_frequency((totalSteps / diagSteps + 1) * (diagRightX - diagLeftX + 1));
	std::vector<dataType> h_modeReal((totalSteps / diagSteps + 1) * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1));
	std::vector<dataType> h_modeImag((totalSteps / diagSteps + 1) * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1));
	std::vector<dataType> h_Epara((totalSteps / diagSteps + 1) * (diagRightX - diagLeftX + 1));
	std::vector<dataType> h_EparaES((totalSteps / diagSteps + 1) * (diagRightX - diagLeftX + 1));
	std::vector<dataType> h_AlphaDensity((totalSteps / diagSteps + 1) * gridNx);
	std::vector<dataType> h_IonDensity((totalSteps / diagSteps + 1) * gridNx);
	std::vector<dataType> h_BeamDensity((totalSteps / diagSteps + 1) * gridNx);
	std::vector<dataType> h_AlphaDiffusivity((totalSteps / diagSteps + 1) * gridNx);
	std::vector<dataType> h_IonDiffusivity((totalSteps / diagSteps + 1) * gridNx);
	std::vector<dataType> h_BeamDiffusivity((totalSteps / diagSteps + 1) * gridNx);

	dataType* d_amplitude[devNums];
	dataType* d_frequency[devNums];
	dataType* d_modeReal[devNums];
	dataType* d_modeImag[devNums];
	dataType* d_Epara[devNums];
	dataType* d_EparaES[devNums];
	dataType* d_AlphaDensity[devNums];
	dataType* d_IonDensity[devNums];
	dataType* d_BeamDensity[devNums];
	dataType* d_AlphaDiffusivity[devNums];
	dataType* d_IonDiffusivity[devNums];
	dataType* d_BeamDiffusivity[devNums];

	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		CUDACHECK(cudaMalloc((void**)&d_amplitude[i], sizeof(dataType) * h_amplitude.size()));
		CUDACHECK(cudaMalloc((void**)&d_frequency[i], sizeof(dataType) * h_frequency.size()));
		CUDACHECK(cudaMalloc((void**)&d_modeReal[i], sizeof(dataType) * h_modeReal.size()));
		CUDACHECK(cudaMalloc((void**)&d_modeImag[i], sizeof(dataType) * h_modeImag.size()));
		CUDACHECK(cudaMalloc((void**)&d_Epara[i], sizeof(dataType) * h_Epara.size()));
		CUDACHECK(cudaMalloc((void**)&d_EparaES[i], sizeof(dataType) * h_EparaES.size()));
		CUDACHECK(cudaMalloc((void**)&d_AlphaDensity[i], sizeof(dataType) * h_AlphaDensity.size()));
		CUDACHECK(cudaMalloc((void**)&d_IonDensity[i], sizeof(dataType) * h_IonDensity.size()));
		CUDACHECK(cudaMalloc((void**)&d_BeamDensity[i], sizeof(dataType) * h_BeamDensity.size()));
		CUDACHECK(cudaMalloc((void**)&d_AlphaDiffusivity[i], sizeof(dataType) * h_AlphaDiffusivity.size()));
		CUDACHECK(cudaMalloc((void**)&d_IonDiffusivity[i], sizeof(dataType) * h_IonDiffusivity.size()));
		CUDACHECK(cudaMalloc((void**)&d_BeamDiffusivity[i], sizeof(dataType) * h_BeamDiffusivity.size()));
		CUDACHECK(cudaMemcpy(d_amplitude[i], h_amplitude.data(), sizeof(dataType) * h_amplitude.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_frequency[i], h_frequency.data(), sizeof(dataType) * h_frequency.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_modeReal[i], h_modeReal.data(), sizeof(dataType) * h_modeReal.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_modeImag[i], h_modeImag.data(), sizeof(dataType) * h_modeImag.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_Epara[i], h_Epara.data(), sizeof(dataType) * h_Epara.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_EparaES[i], h_EparaES.data(), sizeof(dataType) * h_EparaES.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_AlphaDensity[i], h_AlphaDensity.data(), sizeof(dataType) * h_AlphaDensity.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_IonDensity[i], h_IonDensity.data(), sizeof(dataType) * h_IonDensity.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_BeamDensity[i], h_BeamDensity.data(), sizeof(dataType) * h_BeamDensity.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_AlphaDiffusivity[i], h_AlphaDiffusivity.data(), sizeof(dataType) * h_AlphaDiffusivity.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_IonDiffusivity[i], h_IonDiffusivity.data(), sizeof(dataType) * h_IonDiffusivity.size(), cudaMemcpyHostToDevice));
		CUDACHECK(cudaMemcpy(d_BeamDiffusivity[i], h_BeamDiffusivity.data(), sizeof(dataType) * h_BeamDiffusivity.size(), cudaMemcpyHostToDevice));
	}

	std::vector<dataType> h_totalw((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	std::vector<dataType> h_totalA((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	std::vector<dataType> h_totaldNe((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	std::vector<dataType> h_totaldTe((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	std::vector<dataType> h_totalPhi((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	std::vector<dataType> h_totaldPi((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	std::vector<dataType> h_totaldPa((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	std::vector<dataType> h_totaldPb((size_t)(totalSteps / outputSteps + 1) * (rightN - leftN + 1) * gridNy * gridNxz);
	dataType* d_eachw[devNums];
	dataType* d_eachA[devNums];
	dataType* d_eachdNe[devNums];
	dataType* d_eachdTe[devNums];
	dataType* d_eachPhi[devNums];
	dataType* d_eachdPi[devNums];
	dataType* d_eachdPa[devNums];
	dataType* d_eachdPb[devNums];
	dataType* d_totalw[devNums];
	dataType* d_totalA[devNums];
	dataType* d_totaldNe[devNums];
	dataType* d_totaldTe[devNums];
	dataType* d_totalPhi[devNums];
	dataType* d_totaldPi[devNums];
	dataType* d_totaldPa[devNums];
	dataType* d_totaldPb[devNums];

	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		if constexpr (std::is_same_v<ifOutputw, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachw[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachw[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totalw[i], sizeof(dataType) * h_totalw.size()));
			CUDACHECK(cudaMemsetAsync(d_totalw[i], 0, sizeof(dataType) * h_totalw.size(), 0));
		}
		if constexpr (std::is_same_v<ifOutputA, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachA[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachA[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totalA[i], sizeof(dataType) * h_totalA.size()));
			CUDACHECK(cudaMemsetAsync(d_totalA[i], 0, sizeof(dataType) * h_totalA.size(), 0));
		}
		if constexpr (std::is_same_v<ifOutputdNe, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachdNe[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachdNe[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totaldNe[i], sizeof(dataType) * h_totaldNe.size()));
			CUDACHECK(cudaMemsetAsync(d_totaldNe[i], 0, sizeof(dataType) * h_totaldNe.size(), 0));
		}
		if constexpr (std::is_same_v<ifOutputdTe, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachdTe[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachdTe[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totaldTe[i], sizeof(dataType) * h_totaldTe.size()));
			CUDACHECK(cudaMemsetAsync(d_totaldTe[i], 0, sizeof(dataType) * h_totaldTe.size(), 0));
		}
		if constexpr (std::is_same_v<ifOutputPhi, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachPhi[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachPhi[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totalPhi[i], sizeof(dataType) * h_totalPhi.size()));
			CUDACHECK(cudaMemsetAsync(d_totalPhi[i], 0, sizeof(dataType) * h_totalPhi.size(), 0));
		}
		if constexpr (std::is_same_v<ifOutputdPi, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachdPi[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachdPi[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totaldPi[i], sizeof(dataType) * h_totaldPi.size()));
			CUDACHECK(cudaMemsetAsync(d_totaldPi[i], 0, sizeof(dataType) * h_totaldPi.size(), 0));
		}
		if constexpr (std::is_same_v<ifOutputdPa, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachdPa[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachdPa[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totaldPa[i], sizeof(dataType) * h_totaldPa.size()));
			CUDACHECK(cudaMemsetAsync(d_totaldPa[i], 0, sizeof(dataType) * h_totaldPa.size(), 0));
		}
		if constexpr (std::is_same_v<ifOutputdPb, trueType>) {
			CUDACHECK(cudaMalloc((void**)&d_eachdPb[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz));
			CUDACHECK(cudaMemsetAsync(d_eachdPb[i], 0, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, 0));
			CUDACHECK(cudaMalloc((void**)&d_totaldPb[i], sizeof(dataType) * h_totaldPb.size()));
			CUDACHECK(cudaMemsetAsync(d_totaldPb[i], 0, sizeof(dataType) * h_totaldPb.size(), 0));
		}
	}

	dataType* Phi_yxz[devNums];
	dataType* Phi_xzy[devNums];
	dataType* A_yxz[devNums];
	dataType* A_xzy[devNums];
	dataType* dNe_yxz[devNums];
	dataType* dNe_xzy[devNums];
	dataType* dTe_yxz[devNums];
	dataType* dTe_xzy[devNums];
	dataType* dPi_yxz[devNums];
	dataType* dPi_xzy[devNums];
	dataType* dPa_yxz[devNums];
	dataType* dPa_xzy[devNums];
	dataType* dPb_yxz[devNums];
	dataType* dPb_xzy[devNums];

	for (int i = 0; i < devNums; i++) {
		CUDACHECK(cudaSetDevice(localRank * devNums + i));
		CUDACHECK(cudaMalloc((void**)&Phi_yxz[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(Phi_yxz[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&Phi_xzy[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(Phi_xzy[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&A_yxz[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(A_yxz[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&A_xzy[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(A_xzy[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dNe_yxz[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dNe_yxz[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dNe_xzy[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dNe_xzy[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dTe_yxz[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dTe_yxz[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dTe_xzy[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dTe_xzy[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dPi_yxz[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dPi_yxz[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dPi_xzy[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dPi_xzy[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dPa_yxz[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dPa_yxz[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dPa_xzy[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dPa_xzy[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dPb_yxz[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dPb_yxz[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
		CUDACHECK(cudaMalloc((void**)&dPb_xzy[i], sizeof(dataType) * gridNy * gridNxz));
		CUDACHECK(cudaMemsetAsync(dPb_xzy[i], 0, sizeof(dataType) * gridNy * gridNxz, 0));
	}

	void* Ion_storage[devNums] = {};
	size_t  Ion_storage_bytes[devNums];
	void* Alpha_storage[devNums] = {};
	size_t  Alpha_storage_bytes[devNums];
	void* Beam_storage[devNums] = {};
	size_t  Beam_storage_bytes[devNums];

	if (myRank == 0)
		std::cout << BOLDYELLOW << "Start: Allocate memory for sorting particles." << RESET << std::endl;

	if constexpr (std::is_same_v<ifIon, trueType>) {
		for (int i = 0; i < devNums; i++) {
			CUDACHECK(cudaSetDevice(localRank * devNums + i));
			cub::DeviceSegmentedRadixSort::SortPairs(Ion_storage[i], Ion_storage_bytes[i], Ion_keys_out[i], Ion_keys_in[i], Ion_values_out[i], Ion_values_in[i], picDev * 7, 7, Ion_offsets[i], Ion_offsets[i] + 1);
			CUDACHECK(cudaMalloc(&Ion_storage[i], Ion_storage_bytes[i]));
			cub::DeviceSegmentedRadixSort::SortPairs(Ion_storage[i], Ion_storage_bytes[i], Ion_keys_out[i], Ion_keys_in[i], Ion_values_out[i], Ion_values_in[i], picDev * 7, 7, Ion_offsets[i], Ion_offsets[i] + 1);
		}
	}
	if constexpr (std::is_same_v<ifAlpha, trueType>) {
		for (int i = 0; i < devNums; i++) {
			CUDACHECK(cudaSetDevice(localRank * devNums + i));
			cub::DeviceSegmentedRadixSort::SortPairs(Alpha_storage[i], Alpha_storage_bytes[i], Alpha_keys_out[i], Alpha_keys_in[i], Alpha_values_out[i], Alpha_values_in[i], picDev * 7, 7, Alpha_offsets[i], Alpha_offsets[i] + 1);
			CUDACHECK(cudaMalloc(&Alpha_storage[i], Alpha_storage_bytes[i]));
			cub::DeviceSegmentedRadixSort::SortPairs(Alpha_storage[i], Alpha_storage_bytes[i], Alpha_keys_out[i], Alpha_keys_in[i], Alpha_values_out[i], Alpha_values_in[i], picDev * 7, 7, Alpha_offsets[i], Alpha_offsets[i] + 1);
		}
	}
	if constexpr (std::is_same_v<ifBeam, trueType>) {
		for (int i = 0; i < devNums; i++) {
			CUDACHECK(cudaSetDevice(localRank * devNums + i));
			cub::DeviceSegmentedRadixSort::SortPairs(Beam_storage[i], Beam_storage_bytes[i], Beam_keys_out[i], Beam_keys_in[i], Beam_values_out[i], Beam_values_in[i], picDev * 7, 7, Beam_offsets[i], Beam_offsets[i] + 1);
			CUDACHECK(cudaMalloc(&Beam_storage[i], Beam_storage_bytes[i]));
			cub::DeviceSegmentedRadixSort::SortPairs(Beam_storage[i], Beam_storage_bytes[i], Beam_keys_out[i], Beam_keys_in[i], Beam_values_out[i], Beam_values_in[i], picDev * 7, 7, Beam_offsets[i], Beam_offsets[i] + 1);
		}
	}

	if (myRank == 0) {
		size_t avail, total, used;
		cudaSetDevice(localRank * devNums);
		cudaMemGetInfo(&avail, &total);
		used = total - avail;
		std::cout << BOLDYELLOW << "Device memory used: " << (double)used / 1024 / 1024 / 1024 << " GB." << RESET << std::endl;
		std::cout << BOLDGREEN << "Done." << RESET << std::endl;
		std::cout << std::endl;
	}

	std::vector<std::vector<cudssHandle_t>>& cudssHandles = cuGMEC.cudssHandles;

	std::vector<std::vector<cudssConfig_t>>& laplacianConfigs = cuGMEC.laplacianConfigs;
	std::vector<std::vector<cudssData_t>>& laplacianDatas = cuGMEC.laplacianDatas;
	std::vector<std::vector<cudssMatrix_t>>& laplacianAs = cuGMEC.laplacianAs;
	std::vector<std::vector<cudssMatrix_t>>& laplacianXs = cuGMEC.laplacianXs;
	std::vector<std::vector<cudssMatrix_t>>& laplacianBs = cuGMEC.laplacianBs;

	std::vector<std::vector<cudssConfig_t>>& resistiveConfigs = cuGMEC.resistiveConfigs;
	std::vector<std::vector<cudssData_t>>& resistiveDatas = cuGMEC.resistiveDatas;
	std::vector<std::vector<cudssMatrix_t>>& resistiveAs = cuGMEC.resistiveAs;
	std::vector<std::vector<cudssMatrix_t>>& resistiveXs = cuGMEC.resistiveXs;
	std::vector<std::vector<cudssMatrix_t>>& resistiveBs = cuGMEC.resistiveBs;

	std::vector<std::vector<cudssConfig_t>>& wConfigs = cuGMEC.wConfigs;
	std::vector<std::vector<cudssData_t>>& wDatas = cuGMEC.wDatas;
	std::vector<std::vector<cudssMatrix_t>>& wAs = cuGMEC.wAs;
	std::vector<std::vector<cudssMatrix_t>>& wXs = cuGMEC.wXs;
	std::vector<std::vector<cudssMatrix_t>>& wBs = cuGMEC.wBs;

	std::vector<std::vector<cudssConfig_t>>& dNeConfigs = cuGMEC.dNeConfigs;
	std::vector<std::vector<cudssData_t>>& dNeDatas = cuGMEC.dNeDatas;
	std::vector<std::vector<cudssMatrix_t>>& dNeAs = cuGMEC.dNeAs;
	std::vector<std::vector<cudssMatrix_t>>& dNeXs = cuGMEC.dNeXs;
	std::vector<std::vector<cudssMatrix_t>>& dNeBs = cuGMEC.dNeBs;

	std::vector<std::vector<cudssConfig_t>>& dTeConfigs = cuGMEC.dTeConfigs;
	std::vector<std::vector<cudssData_t>>& dTeDatas = cuGMEC.dTeDatas;
	std::vector<std::vector<cudssMatrix_t>>& dTeAs = cuGMEC.dTeAs;
	std::vector<std::vector<cudssMatrix_t>>& dTeXs = cuGMEC.dTeXs;
	std::vector<std::vector<cudssMatrix_t>>& dTeBs = cuGMEC.dTeBs;

	std::vector<std::vector<cudssConfig_t>>& dPiConfigs = cuGMEC.dPiConfigs;
	std::vector<std::vector<cudssData_t>>& dPiDatas = cuGMEC.dPiDatas;
	std::vector<std::vector<cudssMatrix_t>>& dPiAs = cuGMEC.dPiAs;
	std::vector<std::vector<cudssMatrix_t>>& dPiXs = cuGMEC.dPiXs;
	std::vector<std::vector<cudssMatrix_t>>& dPiBs = cuGMEC.dPiBs;

	std::vector<std::vector<cudssConfig_t>>& dPaConfigs = cuGMEC.dPaConfigs;
	std::vector<std::vector<cudssData_t>>& dPaDatas = cuGMEC.dPaDatas;
	std::vector<std::vector<cudssMatrix_t>>& dPaAs = cuGMEC.dPaAs;
	std::vector<std::vector<cudssMatrix_t>>& dPaXs = cuGMEC.dPaXs;
	std::vector<std::vector<cudssMatrix_t>>& dPaBs = cuGMEC.dPaBs;

	std::vector<std::vector<cudssConfig_t>>& dPbConfigs = cuGMEC.dPbConfigs;
	std::vector<std::vector<cudssData_t>>& dPbDatas = cuGMEC.dPbDatas;
	std::vector<std::vector<cudssMatrix_t>>& dPbAs = cuGMEC.dPbAs;
	std::vector<std::vector<cudssMatrix_t>>& dPbXs = cuGMEC.dPbXs;
	std::vector<std::vector<cudssMatrix_t>>& dPbBs = cuGMEC.dPbBs;

	if (std::is_same_v<ifContinue, trueType> && continueSteps != 0) {

		for (int i = 0; i < devNums; i++) {
			cudaSetDevice(localRank * devNums + i);
			for (int j = 0; j < devNy; j++) {
				cudssMatrixSetValues(laplacianBs[i][j], w_midl[i] + (j + gridGhost) * gridNxz);
				cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, laplacianConfigs[i][j], laplacianDatas[i][j], laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]);
			}
			MHD2dJpBdPePhi<ifNonlinear, ifLocal, ifFLRMHD> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (
				A_midl[i], dJpB_midl[i], A2dJpB[i], w_midl[i], Phi_midl[i], w2Phi[i], dNe_midl[i], dTe_midl[i], dPe_midl[i], Ne0[i], Te0[i]);
		}

	}
	else {

		for (int i = 0; i < devNums; i++) {
			cudaSetDevice(localRank * devNums + i);
			MHD2w<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midl[i], w_midl[i], Phi2w[i]);
			MHD2dJpBdPePhi<ifNonlinear, ifLocal, ifFLRMHD> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (
				A_midl[i], dJpB_midl[i], A2dJpB[i], w_midl[i], Phi_midl[i], w2Phi[i], dNe_midl[i], dTe_midl[i], dPe_midl[i], Ne0[i], Te0[i]);
			cudaMemcpyAsync(w_beg[i], w_midl[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
		}

	}

	ncclGroupStart();
	for (int i = 0; i < devNums; i++) {

		ncclSend(w_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		ncclRecv(w_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclSend(w_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclRecv(w_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

		ncclSend(A_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		ncclRecv(A_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclSend(A_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclRecv(A_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

		ncclSend(dNe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		ncclRecv(dNe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclSend(dNe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclRecv(dNe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

		ncclSend(dTe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		ncclRecv(dTe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclSend(dTe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclRecv(dTe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

		ncclSend(Phi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		ncclRecv(Phi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclSend(Phi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclRecv(Phi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

		ncclSend(dJpB_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		ncclRecv(dJpB_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclSend(dJpB_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclRecv(dJpB_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

		ncclSend(dPe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		ncclRecv(dPe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclSend(dPe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
		ncclRecv(dPe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

		if constexpr (std::is_same_v<ifIon, trueType>) {
			ncclSend(dPi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
			ncclRecv(dPi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
			ncclSend(dPi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
			ncclRecv(dPi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		}

		if constexpr (std::is_same_v<ifAlpha, trueType>) {
			ncclSend(dPa_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
			ncclRecv(dPa_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
			ncclSend(dPa_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
			ncclRecv(dPa_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		}

		if constexpr (std::is_same_v<ifBeam, trueType>) {
			ncclSend(dPb_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
			ncclRecv(dPb_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
			ncclSend(dPb_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
			ncclRecv(dPb_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
		}

	}
	ncclGroupEnd();

	std::vector<cudaEvent_t> start(devNums);
	std::vector<cudaEvent_t> end(devNums);
	std::vector<float> time(devNums);

	int diagIndex = 0;
	int outputIndex = 0;

	if (myRank == 0) {
		std::cout << BOLDGREEN << "Gyrokinetic-MHD hybrid simulation is running." << RESET << std::endl;
		std::cout << std::endl;
	}

	for (int i = 0; i < devNums; i++) {
		cudaSetDevice(localRank * devNums + i);
		cudaEventCreate(&start[i]);
		cudaEventCreate(&end[i]);
	}

	for (int i = 0; i < devNums; i++) {
		cudaSetDevice(localRank * devNums + i);
		cudaEventRecord(start[i]);
	}

	if (diagIndex % diagSteps == 0) {

		if constexpr (std::is_same_v<ifDiagAmplitude, trueType>) {
			for (int i = 0; i < devNums; i++) {
				cudaSetDevice(localRank * devNums + i);
				if constexpr (std::is_same_v<dataType, double>) {
					cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
					MHDDiagAmplitude << <1, diagRightX - diagLeftX + 1, 0, 0 >> > (nFreqd[i],
						d_amplitude[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
						d_modeReal[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
						d_modeImag[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1));
				}
				else {
					cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
					MHDDiagAmplitude << <1, diagRightX - diagLeftX + 1, 0, 0 >> > (nFreqf[i],
						d_amplitude[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
						d_modeReal[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
						d_modeImag[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1));
				}
			}
		}

		if constexpr (std::is_same_v<ifDiagFrequency, trueType>) {
			for (int i = 0; i < devNums; i++) {
				cudaSetDevice(localRank * devNums + i);
				MHDDiagFrequency << <1, diagRightX - diagLeftX + 1, 0, 0 >> > (Phi_midl[i], d_frequency[i] + diagIndex / diagSteps * (diagRightX - diagLeftX + 1));
			}
		}

		if constexpr (std::is_same_v<ifDiagEparallel, trueType>) {
			for (int i = 0; i < devNums; i++) {
				cudaSetDevice(localRank * devNums + i);
				MHDDiagEparallel<ifNonlinear, ifStaggered, ifEparallel> << <8, (diagRightX - diagLeftX + 1) / 8, 0, 0 >> > (qtheta[i],
					A_midl[i], dNe_midl[i], dTe_midl[i], Phi_midl[i], Ne0[i], Te0[i], Ne0_px[i], APhidNe2A[i], PhiA_A[i], NeA_A[i],
					d_Epara[i] + diagIndex / diagSteps * (diagRightX - diagLeftX + 1), d_EparaES[i] + diagIndex / diagSteps * (diagRightX - diagLeftX + 1));
			}
		}

		if constexpr (std::is_same_v<ifDiagDensity, trueType>) {

			for (int i = 0; i < devNums; i++) {
				cudaSetDevice(localRank * devNums + i);
				if constexpr (std::is_same_v<ifIon, trueType>)
					PICDiagDensity<Ion> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic2d[i], Ion_keys_in[i], Ion_values_in[i], d_IonDensity[i] + diagIndex / diagSteps * gridNx);
				if constexpr (std::is_same_v<ifAlpha, trueType>)
					PICDiagDensity<Alpha> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic2d[i], Alpha_keys_in[i], Alpha_values_in[i], d_AlphaDensity[i] + diagIndex / diagSteps * gridNx);
				if constexpr (std::is_same_v<ifBeam, trueType>)
					PICDiagDensity<Beam> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic2d[i], Beam_keys_in[i], Beam_values_in[i], d_BeamDensity[i] + diagIndex / diagSteps * gridNx);
			}

		}

	}

	if (outputIndex % outputSteps == 0) {

		for (int mode = leftN; mode <= rightN; mode++) {

			for (int i = 0; i < devNums; i++) {
				cudaSetDevice(localRank * devNums + i);
				if constexpr (std::is_same_v<dataType, double>) {
					if constexpr (std::is_same_v<ifOutputw, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)w_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachw[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachw[i]);
					}
					if constexpr (std::is_same_v<ifOutputA, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)A_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachA[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachA[i]);
					}
					if constexpr (std::is_same_v<ifOutputdNe, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)dNe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdNe[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdNe[i]);
					}
					if constexpr (std::is_same_v<ifOutputdTe, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)dTe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdTe[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdTe[i]);
					}
					if constexpr (std::is_same_v<ifOutputPhi, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachPhi[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachPhi[i]);
					}
					if constexpr (std::is_same_v<ifOutputdPi, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)dPi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdPi[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPi[i]);
					}
					if constexpr (std::is_same_v<ifOutputdPa, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)dPa_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdPa[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPa[i]);
					}
					if constexpr (std::is_same_v<ifOutputdPb, trueType>) {
						cufftExecD2Z(nPlanR2Cs[i], (double*)dPb_midl[i] + gridGhost * gridNxz, nFreqd[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
						cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdPb[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPb[i]);
					}
				}
				else {
					if constexpr (std::is_same_v<ifOutputw, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)w_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachw[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachw[i]);
					}
					if constexpr (std::is_same_v<ifOutputA, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)A_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachA[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachA[i]);
					}
					if constexpr (std::is_same_v<ifOutputdNe, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)dNe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdNe[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdNe[i]);
					}
					if constexpr (std::is_same_v<ifOutputdTe, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)dTe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdTe[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdTe[i]);
					}
					if constexpr (std::is_same_v<ifOutputPhi, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachPhi[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachPhi[i]);
					}
					if constexpr (std::is_same_v<ifOutputdPi, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)dPi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdPi[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPi[i]);
					}
					if constexpr (std::is_same_v<ifOutputdPa, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)dPa_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdPa[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPa[i]);
					}
					if constexpr (std::is_same_v<ifOutputdPb, trueType>) {
						cufftExecR2C(nPlanR2Cs[i], (float*)dPb_midl[i] + gridGhost * gridNxz, nFreqf[i]);
						MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
						cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdPb[i] + gridGhost * gridNxz);
						MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPb[i]);
					}
				}
			}

			ncclGroupStart();
			for (int i = 0; i < devNums; i++) {
				if constexpr (std::is_same_v<ifOutputw, trueType>)
					ncclAllGather(d_eachw[i] + gridGhost * gridNxz, d_totalw[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				if constexpr (std::is_same_v<ifOutputA, trueType>)
					ncclAllGather(d_eachA[i] + gridGhost * gridNxz, d_totalA[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				if constexpr (std::is_same_v<ifOutputdNe, trueType>)
					ncclAllGather(d_eachdNe[i] + gridGhost * gridNxz, d_totaldNe[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				if constexpr (std::is_same_v<ifOutputdTe, trueType>)
					ncclAllGather(d_eachdTe[i] + gridGhost * gridNxz, d_totaldTe[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				if constexpr (std::is_same_v<ifOutputPhi, trueType>)
					ncclAllGather(d_eachPhi[i] + gridGhost * gridNxz, d_totalPhi[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				if constexpr (std::is_same_v<ifOutputdPi, trueType>)
					ncclAllGather(d_eachdPi[i] + gridGhost * gridNxz, d_totaldPi[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				if constexpr (std::is_same_v<ifOutputdPa, trueType>)
					ncclAllGather(d_eachdPa[i] + gridGhost * gridNxz, d_totaldPa[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				if constexpr (std::is_same_v<ifOutputdPb, trueType>)
					ncclAllGather(d_eachdPb[i] + gridGhost * gridNxz, d_totaldPb[i]
						+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
			}
			ncclGroupEnd();

		}

	}

	for (int outerLoop = 0; outerLoop < outerLoopMax; outerLoop++) {

		if (outerLoopMax / 10 != 0) {

			if (outerLoop % (outerLoopMax / 10) == 0) {
				if (myRank == 0) {
					std::cout << BOLDGREEN << 100 * outerLoop / outerLoopMax << "%" << RESET << std::endl;
				}
			}

		}
		else {

			if (myRank == 0) {
				std::cout << BOLDGREEN << 100 * outerLoop / outerLoopMax << "%" << RESET << std::endl;
			}

		}

		for (int innerLoop = 0; innerLoop < innerLoopMax; innerLoop++) {

			for (int MHDLoop = 0; MHDLoop < MHDLoopMax; MHDLoop++) {

				/*--------------------------------------MHD RK4--------------------------------------*/

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHDLinearRK4<1, ifNonlinear, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
						w_beg[i], w_midl[i], w_midr[i], w_end[i], A_beg[i], A_midl[i], A_midr[i], A_end[i],
						dNe_beg[i], dNe_midl[i], dNe_midr[i], dNe_end[i], dTe_beg[i], dTe_midl[i], dTe_midr[i], dTe_end[i],
						Phi_midl[i], dJpB_midl[i], dPe_midl[i], dPi_midl[i], dPa_midl[i], dPb_midl[i],
						wdPAdJpB2w[i], APhidNe2A[i], dPePhiAdJpB2dNe[i], PhidTedNe2dTe[i]);
					if constexpr (std::is_same_v<ifNonlinear, trueType>) {
						MHDNonlinearRK4<1, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
							w_midl[i], w_midr[i], w_end[i], A_midl[i], A_midr[i], A_end[i], dNe_midl[i], dNe_midr[i], dNe_end[i], dTe_midl[i], dTe_midr[i], dTe_end[i],
							Phi_midl[i], dJpB_midl[i], dPe_midl[i], Ne0[i], Te0[i], Ne0_px[i], Te0_px[i],
							APhidNe2A[i], wPhi_w[i], AdJpB_w[i], PhiA_A[i], NeA_A[i], dNePhi_dNe[i], PhiTe_dTe[i], PhiTeA_dTe[i]);
					}
					for (int j = 0; j < devNy; j++) {
						cudssMatrixSetValues(laplacianBs[i][j], w_midr[i] + (j + gridGhost) * gridNxz);
						cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, laplacianConfigs[i][j], laplacianDatas[i][j], laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]);
					}
					MHD2dJpBdPePhi<ifNonlinear, ifLocal, ifFLRMHD> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (
						A_midr[i], dJpB_midl[i], A2dJpB[i], w_midr[i], Phi_midl[i], w2Phi[i], dNe_midr[i], dTe_midr[i], dPe_midl[i], Ne0[i], Te0[i]);
				}

				ncclGroupStart();
				for (int i = 0; i < devNums; i++) {

					ncclSend(w_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(w_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(w_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(w_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(A_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(A_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(A_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(A_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dNe_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dNe_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dNe_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dNe_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dTe_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dTe_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dTe_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dTe_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(Phi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(Phi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dJpB_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dJpB_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dPe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dPe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

				}
				ncclGroupEnd();

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHDLinearRK4<2, ifNonlinear, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
						w_beg[i], w_midr[i], w_midl[i], w_end[i], A_beg[i], A_midr[i], A_midl[i], A_end[i],
						dNe_beg[i], dNe_midr[i], dNe_midl[i], dNe_end[i], dTe_beg[i], dTe_midr[i], dTe_midl[i], dTe_end[i],
						Phi_midl[i], dJpB_midl[i], dPe_midl[i], dPi_midl[i], dPa_midl[i], dPb_midl[i],
						wdPAdJpB2w[i], APhidNe2A[i], dPePhiAdJpB2dNe[i], PhidTedNe2dTe[i]);
					if constexpr (std::is_same_v<ifNonlinear, trueType>) {
						MHDNonlinearRK4<2, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
							w_midr[i], w_midl[i], w_end[i], A_midr[i], A_midl[i], A_end[i], dNe_midr[i], dNe_midl[i], dNe_end[i], dTe_midr[i], dTe_midl[i], dTe_end[i],
							Phi_midl[i], dJpB_midl[i], dPe_midl[i], Ne0[i], Te0[i], Ne0_px[i], Te0_px[i],
							APhidNe2A[i], wPhi_w[i], AdJpB_w[i], PhiA_A[i], NeA_A[i], dNePhi_dNe[i], PhiTe_dTe[i], PhiTeA_dTe[i]);
					}
					for (int j = 0; j < devNy; j++) {
						cudssMatrixSetValues(laplacianBs[i][j], w_midl[i] + (j + gridGhost) * gridNxz);
						cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, laplacianConfigs[i][j], laplacianDatas[i][j], laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]);
					}
					MHD2dJpBdPePhi<ifNonlinear, ifLocal, ifFLRMHD> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (
						A_midl[i], dJpB_midl[i], A2dJpB[i], w_midl[i], Phi_midl[i], w2Phi[i], dNe_midl[i], dTe_midl[i], dPe_midl[i], Ne0[i], Te0[i]);
				}

				ncclGroupStart();
				for (int i = 0; i < devNums; i++) {

					ncclSend(w_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(w_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(w_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(w_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(A_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(A_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(A_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(A_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dNe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dNe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dNe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dNe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dTe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dTe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dTe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dTe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(Phi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(Phi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dJpB_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dJpB_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dPe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dPe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

				}
				ncclGroupEnd();

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHDLinearRK4<3, ifNonlinear, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
						w_beg[i], w_midl[i], w_midr[i], w_end[i], A_beg[i], A_midl[i], A_midr[i], A_end[i],
						dNe_beg[i], dNe_midl[i], dNe_midr[i], dNe_end[i], dTe_beg[i], dTe_midl[i], dTe_midr[i], dTe_end[i],
						Phi_midl[i], dJpB_midl[i], dPe_midl[i], dPi_midl[i], dPa_midl[i], dPb_midl[i],
						wdPAdJpB2w[i], APhidNe2A[i], dPePhiAdJpB2dNe[i], PhidTedNe2dTe[i]);
					if constexpr (std::is_same_v<ifNonlinear, trueType>) {
						MHDNonlinearRK4<3, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
							w_midl[i], w_midr[i], w_end[i], A_midl[i], A_midr[i], A_end[i], dNe_midl[i], dNe_midr[i], dNe_end[i], dTe_midl[i], dTe_midr[i], dTe_end[i],
							Phi_midl[i], dJpB_midl[i], dPe_midl[i], Ne0[i], Te0[i], Ne0_px[i], Te0_px[i],
							APhidNe2A[i], wPhi_w[i], AdJpB_w[i], PhiA_A[i], NeA_A[i], dNePhi_dNe[i], PhiTe_dTe[i], PhiTeA_dTe[i]);
					}
					for (int j = 0; j < devNy; j++) {
						cudssMatrixSetValues(laplacianBs[i][j], w_midr[i] + (j + gridGhost) * gridNxz);
						cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, laplacianConfigs[i][j], laplacianDatas[i][j], laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]);
					}
					MHD2dJpBdPePhi<ifNonlinear, ifLocal, ifFLRMHD> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (
						A_midr[i], dJpB_midl[i], A2dJpB[i], w_midr[i], Phi_midl[i], w2Phi[i], dNe_midr[i], dTe_midr[i], dPe_midl[i], Ne0[i], Te0[i]);
				}

				ncclGroupStart();
				for (int i = 0; i < devNums; i++) {

					ncclSend(w_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(w_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(w_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(w_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(A_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(A_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(A_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(A_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dNe_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dNe_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dNe_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dNe_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dTe_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dTe_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dTe_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dTe_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(Phi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(Phi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dJpB_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dJpB_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dPe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dPe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

				}
				ncclGroupEnd();

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHDLinearRK4<4, ifNonlinear, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
						w_beg[i], w_midr[i], w_midl[i], w_end[i], A_beg[i], A_midr[i], A_midl[i], A_end[i],
						dNe_beg[i], dNe_midr[i], dNe_midl[i], dNe_end[i], dTe_beg[i], dTe_midr[i], dTe_midl[i], dTe_end[i],
						Phi_midl[i], dJpB_midl[i], dPe_midl[i], dPi_midl[i], dPa_midl[i], dPb_midl[i],
						wdPAdJpB2w[i], APhidNe2A[i], dPePhiAdJpB2dNe[i], PhidTedNe2dTe[i]);
					if constexpr (std::is_same_v<ifNonlinear, trueType>) {
						MHDNonlinearRK4<4, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
							w_midr[i], w_midl[i], w_end[i], A_midr[i], A_midl[i], A_end[i], dNe_midr[i], dNe_midl[i], dNe_end[i], dTe_midr[i], dTe_midl[i], dTe_end[i],
							Phi_midl[i], dJpB_midl[i], dPe_midl[i], Ne0[i], Te0[i], Ne0_px[i], Te0_px[i],
							APhidNe2A[i], wPhi_w[i], AdJpB_w[i], PhiA_A[i], NeA_A[i], dNePhi_dNe[i], PhiTe_dTe[i], PhiTeA_dTe[i]);
					}
					for (int j = 0; j < devNy; j++) {
						cudssMatrixSetValues(laplacianBs[i][j], w_midl[i] + (j + gridGhost) * gridNxz);
						cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, laplacianConfigs[i][j], laplacianDatas[i][j], laplacianAs[i][j], laplacianXs[i][j], laplacianBs[i][j]);
					}
				}

				if constexpr (std::is_same_v<ifNablaPerp2Phi, trueType> || std::is_same_v<ifNablaPerp2A, trueType>
					|| std::is_same_v<ifNablaPerp2dNe, trueType> || std::is_same_v<ifNablaPerp2dTe, trueType>) {

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<ifNablaPerp2Phi, trueType>) {
							for (int j = 0; j < devNy; j++)
								cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, wConfigs[i][j], wDatas[i][j], wAs[i][j], wXs[i][j], wBs[i][j]);
							cudaMemcpyAsync(Phi_midl[i], Phi_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						}
						if constexpr (std::is_same_v<ifNablaPerp2A, trueType>) {
							for (int j = 0; j < devNy; j++)
								cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, resistiveConfigs[i][j], resistiveDatas[i][j], resistiveAs[i][j], resistiveXs[i][j], resistiveBs[i][j]);
							cudaMemcpyAsync(A_midl[i], A_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						}
						if constexpr (std::is_same_v<ifNablaPerp2dNe, trueType>) {
							for (int j = 0; j < devNy; j++)
								cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, dNeConfigs[i][j], dNeDatas[i][j], dNeAs[i][j], dNeXs[i][j], dNeBs[i][j]);
							cudaMemcpyAsync(dNe_midl[i], dNe_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						}
						if constexpr (std::is_same_v<ifNablaPerp2dTe, trueType>) {
							for (int j = 0; j < devNy; j++)
								cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, dTeConfigs[i][j], dTeDatas[i][j], dTeAs[i][j], dTeXs[i][j], dTeBs[i][j]);
							cudaMemcpyAsync(dTe_midl[i], dTe_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						}
					}

				}

				if constexpr (std::is_same_v<ifNablaPara4Phi, trueType> || std::is_same_v<ifNablaPara4A, trueType>
					|| std::is_same_v<ifNablaPara4dNe, trueType> || std::is_same_v<ifNablaPara4dTe, trueType>) {

					ncclGroupStart();
					for (int i = 0; i < devNums; i++) {

						if constexpr (std::is_same_v<ifNablaPara4Phi, trueType>) {
							ncclSend(Phi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(Phi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(Phi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(Phi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4A, trueType>) {
							ncclSend(A_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(A_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(A_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(A_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dNe, trueType>) {
							ncclSend(dNe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dNe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dNe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dNe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dTe, trueType>) {
							ncclSend(dTe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dTe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dTe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dTe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

					}
					ncclGroupEnd();

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<ifNablaPara4Phi, trueType>)
							MHDNablaPara2<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], Phi_midl[i], Phi_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4A, trueType>)
							MHDNablaPara2<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], A_midl[i], A_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dNe, trueType>)
							MHDNablaPara2<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dNe_midl[i], dNe_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dTe, trueType>)
							MHDNablaPara2<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dTe_midl[i], dTe_midr[i]);
					}

					ncclGroupStart();
					for (int i = 0; i < devNums; i++) {

						if constexpr (std::is_same_v<ifNablaPara4Phi, trueType>) {
							ncclSend(Phi_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(Phi_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(Phi_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(Phi_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4A, trueType>) {
							ncclSend(A_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(A_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(A_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(A_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dNe, trueType>) {
							ncclSend(dNe_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dNe_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dNe_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dNe_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dTe, trueType>) {
							ncclSend(dTe_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dTe_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dTe_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dTe_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

					}
					ncclGroupEnd();

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<ifNablaPara4Phi, trueType>)
							MHDNablaPara4<0, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], Phi_midl[i], Phi_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4A, trueType>)
							MHDNablaPara4<1, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], A_midl[i], A_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dNe, trueType>)
							MHDNablaPara4<2, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dNe_midl[i], dNe_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dTe, trueType>)
							MHDNablaPara4<3, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dTe_midl[i], dTe_midr[i]);
					}

				}

				if constexpr (std::is_same_v<ifFilterN_Phi, trueType> || std::is_same_v<ifFilterN_A, trueType>
					|| std::is_same_v<ifFilterN_dNe, trueType> || std::is_same_v<ifFilterN_dTe, trueType>) {

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<dataType, double>) {

							if constexpr (std::is_same_v<ifFilterN_Phi, trueType>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], leftN, rightN);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)Phi_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifFilterN_A, trueType>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)A_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], leftN, rightN);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)A_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifFilterN_dNe, trueType>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)dNe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], leftN, rightN);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dNe_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifFilterN_dTe, trueType>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)dTe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], leftN, rightN);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dTe_midl[i] + gridGhost * gridNxz);
							}

						}
						else {

							if constexpr (std::is_same_v<ifFilterN_Phi, trueType>) {
								cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], leftN, rightN);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)Phi_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifFilterN_A, trueType>) {
								cufftExecR2C(nPlanR2Cs[i], (float*)A_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], leftN, rightN);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)A_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifFilterN_dNe, trueType>) {
								cufftExecR2C(nPlanR2Cs[i], (float*)dNe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], leftN, rightN);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dNe_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifFilterN_dTe, trueType>) {
								cufftExecR2C(nPlanR2Cs[i], (float*)dTe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], leftN, rightN);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dTe_midl[i] + gridGhost * gridNxz);
							}

						}
						if constexpr (std::is_same_v<ifFilterN_Phi, trueType>)
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midl[i]);
						if constexpr (std::is_same_v<ifFilterN_A, trueType>)
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (A_midl[i]);
						if constexpr (std::is_same_v<ifFilterN_dNe, trueType>)
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dNe_midl[i]);
						if constexpr (std::is_same_v<ifFilterN_dTe, trueType>)
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dTe_midl[i]);
					}

				}

				if constexpr (removeN_Phi.size() > 0 || removeN_A.size() > 0 || removeN_dNe.size() > 0 || removeN_dTe.size() > 0) {

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						for (int toroidal : removeN_Phi) {

							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)Phi_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)Phi_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midr[i], Phi_midl[i]);

						}
						for (int toroidal : removeN_A) {

							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)A_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)A_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)A_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)A_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (A_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (A_midr[i], A_midl[i]);

						}
						for (int toroidal : removeN_dNe) {

							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)dNe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dNe_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)dNe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dNe_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dNe_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dNe_midr[i], dNe_midl[i]);

						}
						for (int toroidal : removeN_dTe) {

							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)dTe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dTe_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)dTe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dTe_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dTe_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dTe_midr[i], dTe_midl[i]);

						}
					}

				}

				if constexpr (selectNM_Phi.size() > 0 || selectNM_A.size() > 0 || selectNM_dNe.size() > 0 || selectNM_dTe.size() > 0) {

					for (const auto& [toroidal, poloidalLeft, poloidalRight] : selectNM_Phi) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)Phi_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)Phi_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midr[i], Phi_midl[i]);
						}

						ncclGroupStart();
						for (int i = 0; i < devNums; i++) {
							ncclAllGather(Phi_midr[i] + gridGhost * gridNxz, Phi_yxz[i], devNy * gridNxz, ncclType, comms[i], 0);
						}
						ncclGroupEnd();

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							MHDTransposeLeft << <gridNy, gridNx, 0, 0 >> > (Phi_yxz[i], Phi_xzy[i]);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(mPlanR2Cs[i], (double*)Phi_xzy[i], mFreqd[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqd[i], poloidalLeft, poloidalRight);
								cufftExecZ2D(mPlanC2Rs[i], mFreqd[i], (double*)Phi_xzy[i]);

							}
							else {

								cufftExecR2C(mPlanR2Cs[i], (float*)Phi_xzy[i], mFreqf[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqf[i], poloidalLeft, poloidalRight);
								cufftExecC2R(mPlanC2Rs[i], mFreqf[i], (float*)Phi_xzy[i]);

							}
							MHDTransposeRight << <gridNx, gridNz, 0, 0 >> > (Phi_xzy[i], Phi_yxz[i]);
							cudaMemcpyAsync(Phi_midr[i] + gridGhost * gridNxz, Phi_yxz[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * devNy * gridNxz, cudaMemcpyDeviceToDevice, 0);
							MHDFilterResizeM << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midr[i]);
							MHDAddMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midr[i], Phi_midl[i]);

						}

					}
					for (const auto& [toroidal, poloidalLeft, poloidalRight] : selectNM_A) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)A_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)A_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)A_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)A_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (A_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (A_midr[i], A_midl[i]);
						}

						ncclGroupStart();
						for (int i = 0; i < devNums; i++) {
							ncclAllGather(A_midr[i] + gridGhost * gridNxz, A_yxz[i], devNy * gridNxz, ncclType, comms[i], 0);
						}
						ncclGroupEnd();

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							MHDTransposeLeft << <gridNy, gridNx, 0, 0 >> > (A_yxz[i], A_xzy[i]);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(mPlanR2Cs[i], (double*)A_xzy[i], mFreqd[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqd[i], poloidalLeft, poloidalRight);
								cufftExecZ2D(mPlanC2Rs[i], mFreqd[i], (double*)A_xzy[i]);

							}
							else {

								cufftExecR2C(mPlanR2Cs[i], (float*)A_xzy[i], mFreqf[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqf[i], poloidalLeft, poloidalRight);
								cufftExecC2R(mPlanC2Rs[i], mFreqf[i], (float*)A_xzy[i]);

							}
							MHDTransposeRight << <gridNx, gridNz, 0, 0 >> > (A_xzy[i], A_yxz[i]);
							cudaMemcpyAsync(A_midr[i] + gridGhost * gridNxz, A_yxz[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * devNy * gridNxz, cudaMemcpyDeviceToDevice, 0);
							MHDFilterResizeM << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (A_midr[i]);
							MHDAddMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (A_midr[i], A_midl[i]);

						}

					}
					for (const auto& [toroidal, poloidalLeft, poloidalRight] : selectNM_dNe) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)dNe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dNe_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)dNe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dNe_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dNe_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dNe_midr[i], dNe_midl[i]);
						}

						ncclGroupStart();
						for (int i = 0; i < devNums; i++) {
							ncclAllGather(dNe_midr[i] + gridGhost * gridNxz, dNe_yxz[i], devNy * gridNxz, ncclType, comms[i], 0);
						}
						ncclGroupEnd();

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							MHDTransposeLeft << <gridNy, gridNx, 0, 0 >> > (dNe_yxz[i], dNe_xzy[i]);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(mPlanR2Cs[i], (double*)dNe_xzy[i], mFreqd[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqd[i], poloidalLeft, poloidalRight);
								cufftExecZ2D(mPlanC2Rs[i], mFreqd[i], (double*)dNe_xzy[i]);

							}
							else {

								cufftExecR2C(mPlanR2Cs[i], (float*)dNe_xzy[i], mFreqf[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqf[i], poloidalLeft, poloidalRight);
								cufftExecC2R(mPlanC2Rs[i], mFreqf[i], (float*)dNe_xzy[i]);

							}
							MHDTransposeRight << <gridNx, gridNz, 0, 0 >> > (dNe_xzy[i], dNe_yxz[i]);
							cudaMemcpyAsync(dNe_midr[i] + gridGhost * gridNxz, dNe_yxz[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * devNy * gridNxz, cudaMemcpyDeviceToDevice, 0);
							MHDFilterResizeM << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dNe_midr[i]);
							MHDAddMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dNe_midr[i], dNe_midl[i]);

						}

					}
					for (const auto& [toroidal, poloidalLeft, poloidalRight] : selectNM_dTe) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(nPlanR2Cs[i], (double*)dTe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dTe_midr[i] + gridGhost * gridNxz);

							}
							else {

								cufftExecR2C(nPlanR2Cs[i], (float*)dTe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dTe_midr[i] + gridGhost * gridNxz);

							}
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dTe_midr[i]);
							MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dTe_midr[i], dTe_midl[i]);
						}

						ncclGroupStart();
						for (int i = 0; i < devNums; i++) {
							ncclAllGather(dTe_midr[i] + gridGhost * gridNxz, dTe_yxz[i], devNy * gridNxz, ncclType, comms[i], 0);
						}
						ncclGroupEnd();

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							MHDTransposeLeft << <gridNy, gridNx, 0, 0 >> > (dTe_yxz[i], dTe_xzy[i]);
							if constexpr (std::is_same_v<dataType, double>) {

								cufftExecD2Z(mPlanR2Cs[i], (double*)dTe_xzy[i], mFreqd[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqd[i], poloidalLeft, poloidalRight);
								cufftExecZ2D(mPlanC2Rs[i], mFreqd[i], (double*)dTe_xzy[i]);

							}
							else {

								cufftExecR2C(mPlanR2Cs[i], (float*)dTe_xzy[i], mFreqf[i]);
								MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqf[i], poloidalLeft, poloidalRight);
								cufftExecC2R(mPlanC2Rs[i], mFreqf[i], (float*)dTe_xzy[i]);

							}
							MHDTransposeRight << <gridNx, gridNz, 0, 0 >> > (dTe_xzy[i], dTe_yxz[i]);
							cudaMemcpyAsync(dTe_midr[i] + gridGhost * gridNxz, dTe_yxz[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * devNy * gridNxz, cudaMemcpyDeviceToDevice, 0);
							MHDFilterResizeM << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dTe_midr[i]);
							MHDAddMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dTe_midr[i], dTe_midl[i]);

						}

					}

				}

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHD2w<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (Phi_midl[i], w_midl[i], Phi2w[i]);
					MHD2dJpBdPePhi<ifNonlinear, ifLocal, ifFLRMHD> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (
						A_midl[i], dJpB_midl[i], A2dJpB[i], w_midl[i], Phi_midl[i], w2Phi[i], dNe_midl[i], dTe_midl[i], dPe_midl[i], Ne0[i], Te0[i]);
				}

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					cudaMemcpyAsync(w_beg[i], w_midl[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
					cudaMemcpyAsync(A_beg[i], A_midl[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
					cudaMemcpyAsync(dNe_beg[i], dNe_midl[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
					cudaMemcpyAsync(dTe_beg[i], dTe_midl[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
				}

				ncclGroupStart();
				for (int i = 0; i < devNums; i++) {

					ncclSend(w_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(w_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(w_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(w_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(A_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(A_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(A_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(A_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dNe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dNe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dNe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dNe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dTe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dTe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dTe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dTe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(Phi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(Phi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(Phi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dJpB_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dJpB_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dJpB_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

					ncclSend(dPe_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclSend(dPe_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
					ncclRecv(dPe_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);

				}
				ncclGroupEnd();

				/*----------------------------------------diagnose----------------------------------------*/

				diagIndex++;

				outputIndex++;

				if (diagIndex % diagSteps == 0) {

					if constexpr (std::is_same_v<ifDiagAmplitude, trueType>) {
						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<dataType, double>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDDiagAmplitude << <1, diagRightX - diagLeftX + 1, 0, 0 >> > (nFreqd[i],
									d_amplitude[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
									d_modeReal[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
									d_modeImag[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1));
							}
							else {
								cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDDiagAmplitude << <1, diagRightX - diagLeftX + 1, 0, 0 >> > (nFreqf[i],
									d_amplitude[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
									d_modeReal[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1),
									d_modeImag[i] + diagIndex / diagSteps * (rightN - leftN + 1) * (diagRightX - diagLeftX + 1));
							}
						}
					}

					if constexpr (std::is_same_v<ifDiagFrequency, trueType>) {
						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							MHDDiagFrequency << <1, diagRightX - diagLeftX + 1, 0, 0 >> > (Phi_midl[i], d_frequency[i] + diagIndex / diagSteps * (diagRightX - diagLeftX + 1));
						}
					}

					if constexpr (std::is_same_v<ifDiagEparallel, trueType>) {
						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							MHDDiagEparallel<ifNonlinear, ifStaggered, ifEparallel> << <8, (diagRightX - diagLeftX + 1) / 8, 0, 0 >> > (qtheta[i],
								A_midl[i], dNe_midl[i], dTe_midl[i], Phi_midl[i], Ne0[i], Te0[i], Ne0_px[i], APhidNe2A[i], PhiA_A[i], NeA_A[i],
								d_Epara[i] + diagIndex / diagSteps * (diagRightX - diagLeftX + 1), d_EparaES[i] + diagIndex / diagSteps * (diagRightX - diagLeftX + 1));
						}
					}

					if constexpr (std::is_same_v<ifDiagDensity, trueType>) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<ifIon, trueType>)
								PICDiagDensity<Ion> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic2d[i], Ion_keys_in[i], Ion_values_in[i], d_IonDensity[i] + diagIndex / diagSteps * gridNx);
							if constexpr (std::is_same_v<ifAlpha, trueType>)
								PICDiagDensity<Alpha> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic2d[i], Alpha_keys_in[i], Alpha_values_in[i], d_AlphaDensity[i] + diagIndex / diagSteps * gridNx);
							if constexpr (std::is_same_v<ifBeam, trueType>)
								PICDiagDensity<Beam> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic2d[i], Beam_keys_in[i], Beam_values_in[i], d_BeamDensity[i] + diagIndex / diagSteps * gridNx);
						}

					}

					if constexpr (std::is_same_v<ifDiagDiffusivity, trueType>) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<ifIon, trueType>)
								PICDiagDiffusivity<gyroNums, Ion> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Ion_keys_in[i], Ion_values_in[i], d_IonDiffusivity[i] + diagIndex / diagSteps * gridNx);
							if constexpr (std::is_same_v<ifAlpha, trueType>)
								PICDiagDiffusivity<gyroNums, Alpha> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Alpha_keys_in[i], Alpha_values_in[i], d_AlphaDiffusivity[i] + diagIndex / diagSteps * gridNx);
							if constexpr (std::is_same_v<ifBeam, trueType>)
								PICDiagDiffusivity<gyroNums, Beam> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Beam_keys_in[i], Beam_values_in[i], d_BeamDiffusivity[i] + diagIndex / diagSteps * gridNx);
						}

					}

				}

				if (outputIndex % outputSteps == 0) {

					for (int mode = leftN; mode <= rightN; mode++) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<dataType, double>) {
								if constexpr (std::is_same_v<ifOutputw, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)w_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachw[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachw[i]);
								}
								if constexpr (std::is_same_v<ifOutputA, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)A_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachA[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachA[i]);
								}
								if constexpr (std::is_same_v<ifOutputdNe, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dNe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdNe[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdNe[i]);
								}
								if constexpr (std::is_same_v<ifOutputdTe, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dTe_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdTe[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdTe[i]);
								}
								if constexpr (std::is_same_v<ifOutputPhi, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)Phi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachPhi[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachPhi[i]);
								}
								if constexpr (std::is_same_v<ifOutputdPi, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdPi[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPi[i]);
								}
								if constexpr (std::is_same_v<ifOutputdPa, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPa_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdPa[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPa[i]);
								}
								if constexpr (std::is_same_v<ifOutputdPb, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPb_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], mode);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)d_eachdPb[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPb[i]);
								}
							}
							else {
								if constexpr (std::is_same_v<ifOutputw, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)w_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachw[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachw[i]);
								}
								if constexpr (std::is_same_v<ifOutputA, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)A_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachA[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachA[i]);
								}
								if constexpr (std::is_same_v<ifOutputdNe, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dNe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdNe[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdNe[i]);
								}
								if constexpr (std::is_same_v<ifOutputdTe, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dTe_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdTe[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdTe[i]);
								}
								if constexpr (std::is_same_v<ifOutputPhi, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)Phi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachPhi[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachPhi[i]);
								}
								if constexpr (std::is_same_v<ifOutputdPi, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdPi[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPi[i]);
								}
								if constexpr (std::is_same_v<ifOutputdPa, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPa_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdPa[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPa[i]);
								}
								if constexpr (std::is_same_v<ifOutputdPb, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPb_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], mode);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)d_eachdPb[i] + gridGhost * gridNxz);
									MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (d_eachdPb[i]);
								}
							}
						}

						ncclGroupStart();
						for (int i = 0; i < devNums; i++) {
							if constexpr (std::is_same_v<ifOutputw, trueType>)
								ncclAllGather(d_eachw[i] + gridGhost * gridNxz, d_totalw[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifOutputA, trueType>)
								ncclAllGather(d_eachA[i] + gridGhost * gridNxz, d_totalA[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifOutputdNe, trueType>)
								ncclAllGather(d_eachdNe[i] + gridGhost * gridNxz, d_totaldNe[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifOutputdTe, trueType>)
								ncclAllGather(d_eachdTe[i] + gridGhost * gridNxz, d_totaldTe[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifOutputPhi, trueType>)
								ncclAllGather(d_eachPhi[i] + gridGhost * gridNxz, d_totalPhi[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifOutputdPi, trueType>)
								ncclAllGather(d_eachdPi[i] + gridGhost * gridNxz, d_totaldPi[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifOutputdPa, trueType>)
								ncclAllGather(d_eachdPa[i] + gridGhost * gridNxz, d_totaldPa[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifOutputdPb, trueType>)
								ncclAllGather(d_eachdPb[i] + gridGhost * gridNxz, d_totaldPb[i]
									+ (size_t)outputIndex / outputSteps * (rightN - leftN + 1) * gridNy * gridNxz + (mode - leftN) * gridNy * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
						}
						ncclGroupEnd();

					}

				}

			}

			/*---------------------------------------PIC RK4---------------------------------------*/

			if constexpr (std::is_same_v<ifIon, trueType> || std::is_same_v<ifAlpha, trueType> || std::is_same_v<ifBeam, trueType>) {

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHD2Apt<ifNonlinear, ifLocal, ifStaggered, ifEparallel> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i],
						A_midl[i], dNe_midl[i], dTe_midl[i], Phi_midl[i], Ne0[i], Te0[i], Ne0_px[i], pic_APhidNe2A[i], pic_PhiA_A[i], pic_NeA_A[i], Apt_midl[i]);
					if constexpr (std::is_same_v<ifStaggered, trueType>) {
						MHDStaggered2C<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], A_midl[i], A_midr[i]);
						MHDShifted2A<0, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], A_midr[i], dJpB_midr[i], Phi_midl[i], Phi_midr[i], Apt_midl[i], Apt_midr[i]);
					}
					else {
						MHDShifted2A<0, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], A_midl[i], A_midr[i], Phi_midl[i], Phi_midr[i], Apt_midl[i], Apt_midr[i]);
					}
				}

				ncclGroupStart();
				for (int i = 0; i < devNums; i++) {
					if constexpr (std::is_same_v<ifStaggered, trueType>)
						ncclAllGather(dJpB_midr[i] + gridGhost * gridNxz, globalA[i] + gridGhost * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
					else
						ncclAllGather(A_midr[i] + gridGhost * gridNxz, globalA[i] + gridGhost * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
					ncclAllGather(Phi_midr[i] + gridGhost * gridNxz, globalPhi[i] + gridGhost * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
					ncclAllGather(Apt_midr[i] + gridGhost * gridNxz, globalApt[i] + gridGhost * gridNxz, devNy * gridNxz, ncclType, comms[i], 0);
				}
				ncclGroupEnd();

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHDAlignedGhost<ifLocal> << <GhostGridSize, GhostBlockSize, 0, 0 >> > (qtheta[i], globalA[i], globalPhi[i], globalApt[i]);
					MHD2PIC<ifLocal> << <M2PGridSize, M2PBlockSize, 0, 0 >> > (pic3d[i], globalA[i], globalPhi[i], globalApt[i]);
					if constexpr (std::is_same_v<ifFLRPIC, trueType>) {
						if constexpr (std::is_same_v<ifIon, trueType>)
							GyroAlignedRK4<ratioDt, gyroNums, Ion, IonType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Ion_keys_in[i], Ion_values_in[i], globalPi[i]);
						if constexpr (std::is_same_v<ifAlpha, trueType>)
							GyroAlignedRK4<ratioDt, gyroNums, Alpha, AlphaType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Alpha_keys_in[i], Alpha_values_in[i], globalPa[i]);
						if constexpr (std::is_same_v<ifBeam, trueType>)
							GyroAlignedRK4<ratioDt, gyroNums, Beam, BeamType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Beam_keys_in[i], Beam_values_in[i], globalPb[i]);
					}
					else {
						if constexpr (std::is_same_v<ifIon, trueType>)
							DriftAlignedRK4<ratioDt, Ion, IonType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Ion_keys_in[i], Ion_values_in[i], globalPi[i]);
						if constexpr (std::is_same_v<ifAlpha, trueType>)
							DriftAlignedRK4<ratioDt, Alpha, AlphaType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Alpha_keys_in[i], Alpha_values_in[i], globalPa[i]);
						if constexpr (std::is_same_v<ifBeam, trueType>)
							DriftAlignedRK4<ratioDt, Beam, BeamType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], pic3d[i], Beam_keys_in[i], Beam_values_in[i], globalPb[i]);
					}
				}

				ncclGroupStart();
				for (int i = 0; i < devNums; i++) {
					if constexpr (std::is_same_v<ifIon, trueType>)
						ncclAllReduce(globalPi[i], globalPi[i], (gridNy + 2 * gridGhost) * gridNxz, ncclType, ncclSum, comms[i], 0);
					if constexpr (std::is_same_v<ifAlpha, trueType>)
						ncclAllReduce(globalPa[i], globalPa[i], (gridNy + 2 * gridGhost) * gridNxz, ncclType, ncclSum, comms[i], 0);
					if constexpr (std::is_same_v<ifBeam, trueType>)
						ncclAllReduce(globalPb[i], globalPb[i], (gridNy + 2 * gridGhost) * gridNxz, ncclType, ncclSum, comms[i], 0);
				}
				ncclGroupEnd();

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					if constexpr (std::is_same_v<ifIon, trueType>) {
						PICAlignedGhost<ifLocal> << <GhostGridSize, GhostBlockSize, 0, 0 >> > (qtheta[i], globalPi[i]);
						cudaMemcpyAsync(dPi_midl[i], globalPi[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						cudaMemsetAsync(globalPi[i], 0, sizeof(dataType) * (gridNy + 2 * gridGhost) * gridNxz, 0);
					}
					if constexpr (std::is_same_v<ifAlpha, trueType>) {
						PICAlignedGhost<ifLocal> << <GhostGridSize, GhostBlockSize, 0, 0 >> > (qtheta[i], globalPa[i]);
						cudaMemcpyAsync(dPa_midl[i], globalPa[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						cudaMemsetAsync(globalPa[i], 0, sizeof(dataType) * (gridNy + 2 * gridGhost) * gridNxz, 0);
					}
					if constexpr (std::is_same_v<ifBeam, trueType>) {
						PICAlignedGhost<ifLocal> << <GhostGridSize, GhostBlockSize, 0, 0 >> > (qtheta[i], globalPb[i]);
						cudaMemcpyAsync(dPb_midl[i], globalPb[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						cudaMemsetAsync(globalPb[i], 0, sizeof(dataType) * (gridNy + 2 * gridGhost) * gridNxz, 0);
					}
				}

				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					MHDShifted2A<1, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dPi_midl[i], dPi_midr[i], dPa_midl[i], dPa_midr[i], dPb_midl[i], dPb_midr[i]);
					if constexpr (std::is_same_v<ifIon, trueType>)
						cudaMemcpyAsync(dPi_midl[i], dPi_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
					if constexpr (std::is_same_v<ifAlpha, trueType>)
						cudaMemcpyAsync(dPa_midl[i], dPa_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
					if constexpr (std::is_same_v<ifBeam, trueType>)
						cudaMemcpyAsync(dPb_midl[i], dPb_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
				}

				if constexpr (std::is_same_v<ifNablaPerp2dPi, trueType> || std::is_same_v<ifNablaPerp2dPa, trueType> || std::is_same_v<ifNablaPerp2dPb, trueType>) {

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<ifNablaPerp2dPi, trueType>) {
							for (int j = 0; j < devNy; j++)
								cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, dPiConfigs[i][j], dPiDatas[i][j], dPiAs[i][j], dPiXs[i][j], dPiBs[i][j]);
							cudaMemcpyAsync(dPi_midl[i], dPi_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						}
						if constexpr (std::is_same_v<ifNablaPerp2dPa, trueType>) {
							for (int j = 0; j < devNy; j++)
								cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, dPaConfigs[i][j], dPaDatas[i][j], dPaAs[i][j], dPaXs[i][j], dPaBs[i][j]);
							cudaMemcpyAsync(dPa_midl[i], dPa_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						}
						if constexpr (std::is_same_v<ifNablaPerp2dPb, trueType>) {
							for (int j = 0; j < devNy; j++)
								cudssExecute(cudssHandles[i][j], CUDSS_PHASE_SOLVE, dPbConfigs[i][j], dPbDatas[i][j], dPbAs[i][j], dPbXs[i][j], dPbBs[i][j]);
							cudaMemcpyAsync(dPb_midl[i], dPb_midr[i], sizeof(dataType) * (devNy + 2 * gridGhost) * gridNxz, cudaMemcpyDeviceToDevice, 0);
						}
					}

				}

				if constexpr (std::is_same_v<ifNablaPara4dPi, trueType> || std::is_same_v<ifNablaPara4dPa, trueType> || std::is_same_v<ifNablaPara4dPb, trueType>) {

					ncclGroupStart();
					for (int i = 0; i < devNums; i++) {

						if constexpr (std::is_same_v<ifNablaPara4dPi, trueType>) {
							ncclSend(dPi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dPi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dPi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dPi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dPa, trueType>) {
							ncclSend(dPa_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dPa_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dPa_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dPa_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dPb, trueType>) {
							ncclSend(dPb_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dPb_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dPb_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dPb_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

					}
					ncclGroupEnd();

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<ifNablaPara4dPi, trueType>)
							MHDNablaPara2<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dPi_midl[i], dPi_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dPa, trueType>)
							MHDNablaPara2<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dPa_midl[i], dPa_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dPb, trueType>)
							MHDNablaPara2<ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dPb_midl[i], dPb_midr[i]);
					}

					ncclGroupStart();
					for (int i = 0; i < devNums; i++) {

						if constexpr (std::is_same_v<ifNablaPara4dPi, trueType>) {
							ncclSend(dPi_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dPi_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dPi_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dPi_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dPa, trueType>) {
							ncclSend(dPa_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dPa_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dPa_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dPa_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

						if constexpr (std::is_same_v<ifNablaPara4dPb, trueType>) {
							ncclSend(dPb_midr[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
							ncclRecv(dPb_midr[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclSend(dPb_midr[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
							ncclRecv(dPb_midr[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						}

					}
					ncclGroupEnd();

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<ifNablaPara4dPi, trueType>)
							MHDNablaPara4<4, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dPi_midl[i], dPi_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dPa, trueType>)
							MHDNablaPara4<5, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dPa_midl[i], dPa_midr[i]);
						if constexpr (std::is_same_v<ifNablaPara4dPb, trueType>)
							MHDNablaPara4<6, ifLocal> << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (qtheta[i], dPb_midl[i], dPb_midr[i]);
					}

				}

				if constexpr (std::is_same_v<ifFilterN_dP, trueType>) {

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						if constexpr (std::is_same_v<dataType, double>) {

							if constexpr (std::is_same_v<ifIon, trueType>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)dPi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], leftN, rightN);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPi_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifAlpha, trueType>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)dPa_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], leftN, rightN);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPa_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifBeam, trueType>) {
								cufftExecD2Z(nPlanR2Cs[i], (double*)dPb_midl[i] + gridGhost * gridNxz, nFreqd[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], leftN, rightN);
								cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPb_midl[i] + gridGhost * gridNxz);
							}

						}
						else {

							if constexpr (std::is_same_v<ifIon, trueType>) {
								cufftExecR2C(nPlanR2Cs[i], (float*)dPi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], leftN, rightN);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPi_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifAlpha, trueType>) {
								cufftExecR2C(nPlanR2Cs[i], (float*)dPa_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], leftN, rightN);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPa_midl[i] + gridGhost * gridNxz);
							}

							if constexpr (std::is_same_v<ifBeam, trueType>) {
								cufftExecR2C(nPlanR2Cs[i], (float*)dPb_midl[i] + gridGhost * gridNxz, nFreqf[i]);
								MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], leftN, rightN);
								cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPb_midl[i] + gridGhost * gridNxz);
							}

						}
						if constexpr (std::is_same_v<ifIon, trueType>)
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPi_midl[i]);
						if constexpr (std::is_same_v<ifAlpha, trueType>)
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPa_midl[i]);
						if constexpr (std::is_same_v<ifBeam, trueType>)
							MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPb_midl[i]);
					}

				}

				if constexpr (removeN_dP.size() > 0) {

					for (int i = 0; i < devNums; i++) {
						cudaSetDevice(localRank * devNums + i);
						for (int toroidal : removeN_dP) {

							if constexpr (std::is_same_v<dataType, double>) {

								if constexpr (std::is_same_v<ifIon, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPi_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifAlpha, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPa_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPa_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifBeam, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPb_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPb_midr[i] + gridGhost * gridNxz);
								}

							}
							else {

								if constexpr (std::is_same_v<ifIon, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPi_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifAlpha, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPa_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPa_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifBeam, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPb_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPb_midr[i] + gridGhost * gridNxz);
								}

							}
							if constexpr (std::is_same_v<ifIon, trueType>) {
								MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPi_midr[i]);
								MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPi_midr[i], dPi_midl[i]);
							}
							if constexpr (std::is_same_v<ifAlpha, trueType>) {
								MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPa_midr[i]);
								MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPa_midr[i], dPa_midl[i]);
							}
							if constexpr (std::is_same_v<ifBeam, trueType>) {
								MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPb_midr[i]);
								MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPb_midr[i], dPb_midl[i]);
							}

						}
					}

				}

				if constexpr (selectNM_dP.size() > 0) {

					for (const auto& [toroidal, poloidalLeft, poloidalRight] : selectNM_dP) {

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<dataType, double>) {

								if constexpr (std::is_same_v<ifIon, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPi_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPi_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifAlpha, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPa_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPa_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifBeam, trueType>) {
									cufftExecD2Z(nPlanR2Cs[i], (double*)dPb_midl[i] + gridGhost * gridNxz, nFreqd[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqd[i], toroidal);
									cufftExecZ2D(nPlanC2Rs[i], nFreqd[i], (double*)dPb_midr[i] + gridGhost * gridNxz);
								}

							}
							else {

								if constexpr (std::is_same_v<ifIon, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPi_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPi_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifAlpha, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPa_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPa_midr[i] + gridGhost * gridNxz);
								}

								if constexpr (std::is_same_v<ifBeam, trueType>) {
									cufftExecR2C(nPlanR2Cs[i], (float*)dPb_midl[i] + gridGhost * gridNxz, nFreqf[i]);
									MHDFilterModeN << <devNy, gridNx, 0, 0 >> > (nFreqf[i], toroidal);
									cufftExecC2R(nPlanC2Rs[i], nFreqf[i], (float*)dPb_midr[i] + gridGhost * gridNxz);
								}

							}
							if constexpr (std::is_same_v<ifIon, trueType>) {
								MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPi_midr[i]);
								MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPi_midr[i], dPi_midl[i]);
							}
							if constexpr (std::is_same_v<ifAlpha, trueType>) {
								MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPa_midr[i]);
								MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPa_midr[i], dPa_midl[i]);
							}
							if constexpr (std::is_same_v<ifBeam, trueType>) {
								MHDFilterResizeN << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPb_midr[i]);
								MHDSubtractMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPb_midr[i], dPb_midl[i]);
							}
						}

						ncclGroupStart();
						for (int i = 0; i < devNums; i++) {
							if constexpr (std::is_same_v<ifIon, trueType>)
								ncclAllGather(dPi_midr[i] + gridGhost * gridNxz, dPi_yxz[i], devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifAlpha, trueType>)
								ncclAllGather(dPa_midr[i] + gridGhost * gridNxz, dPa_yxz[i], devNy * gridNxz, ncclType, comms[i], 0);
							if constexpr (std::is_same_v<ifBeam, trueType>)
								ncclAllGather(dPb_midr[i] + gridGhost * gridNxz, dPb_yxz[i], devNy * gridNxz, ncclType, comms[i], 0);
						}
						ncclGroupEnd();

						for (int i = 0; i < devNums; i++) {
							cudaSetDevice(localRank * devNums + i);
							if constexpr (std::is_same_v<ifIon, trueType>)
								MHDTransposeLeft << <gridNy, gridNx, 0, 0 >> > (dPi_yxz[i], dPi_xzy[i]);
							if constexpr (std::is_same_v<ifAlpha, trueType>)
								MHDTransposeLeft << <gridNy, gridNx, 0, 0 >> > (dPa_yxz[i], dPa_xzy[i]);
							if constexpr (std::is_same_v<ifBeam, trueType>)
								MHDTransposeLeft << <gridNy, gridNx, 0, 0 >> > (dPb_yxz[i], dPb_xzy[i]);
							if constexpr (std::is_same_v<dataType, double>) {

								if constexpr (std::is_same_v<ifIon, trueType>) {
									cufftExecD2Z(mPlanR2Cs[i], (double*)dPi_xzy[i], mFreqd[i]);
									MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqd[i], poloidalLeft, poloidalRight);
									cufftExecZ2D(mPlanC2Rs[i], mFreqd[i], (double*)dPi_xzy[i]);
								}

								if constexpr (std::is_same_v<ifAlpha, trueType>) {
									cufftExecD2Z(mPlanR2Cs[i], (double*)dPa_xzy[i], mFreqd[i]);
									MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqd[i], poloidalLeft, poloidalRight);
									cufftExecZ2D(mPlanC2Rs[i], mFreqd[i], (double*)dPa_xzy[i]);
								}

								if constexpr (std::is_same_v<ifBeam, trueType>) {
									cufftExecD2Z(mPlanR2Cs[i], (double*)dPb_xzy[i], mFreqd[i]);
									MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqd[i], poloidalLeft, poloidalRight);
									cufftExecZ2D(mPlanC2Rs[i], mFreqd[i], (double*)dPb_xzy[i]);
								}

							}
							else {

								if constexpr (std::is_same_v<ifIon, trueType>) {
									cufftExecR2C(mPlanR2Cs[i], (float*)dPi_xzy[i], mFreqf[i]);
									MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqf[i], poloidalLeft, poloidalRight);
									cufftExecC2R(mPlanC2Rs[i], mFreqf[i], (float*)dPi_xzy[i]);
								}

								if constexpr (std::is_same_v<ifAlpha, trueType>) {
									cufftExecR2C(mPlanR2Cs[i], (float*)dPa_xzy[i], mFreqf[i]);
									MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqf[i], poloidalLeft, poloidalRight);
									cufftExecC2R(mPlanC2Rs[i], mFreqf[i], (float*)dPa_xzy[i]);
								}

								if constexpr (std::is_same_v<ifBeam, trueType>) {
									cufftExecR2C(mPlanR2Cs[i], (float*)dPb_xzy[i], mFreqf[i]);
									MHDFilterModeM << <gridNx, gridNz, 0, 0 >> > (mFreqf[i], poloidalLeft, poloidalRight);
									cufftExecC2R(mPlanC2Rs[i], mFreqf[i], (float*)dPb_xzy[i]);
								}

							}
							if constexpr (std::is_same_v<ifIon, trueType>) {
								MHDTransposeRight << <gridNx, gridNz, 0, 0 >> > (dPi_xzy[i], dPi_yxz[i]);
								cudaMemcpyAsync(dPi_midr[i] + gridGhost * gridNxz, dPi_yxz[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * devNy * gridNxz, cudaMemcpyDeviceToDevice, 0);
								MHDFilterResizeM << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPi_midr[i]);
								MHDAddMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPi_midr[i], dPi_midl[i]);
							}
							if constexpr (std::is_same_v<ifAlpha, trueType>) {
								MHDTransposeRight << <gridNx, gridNz, 0, 0 >> > (dPa_xzy[i], dPa_yxz[i]);
								cudaMemcpyAsync(dPa_midr[i] + gridGhost * gridNxz, dPa_yxz[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * devNy * gridNxz, cudaMemcpyDeviceToDevice, 0);
								MHDFilterResizeM << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPa_midr[i]);
								MHDAddMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPa_midr[i], dPa_midl[i]);
							}
							if constexpr (std::is_same_v<ifBeam, trueType>) {
								MHDTransposeRight << <gridNx, gridNz, 0, 0 >> > (dPb_xzy[i], dPb_yxz[i]);
								cudaMemcpyAsync(dPb_midr[i] + gridGhost * gridNxz, dPb_yxz[i] + (myRank * hostNy + i * devNy) * gridNxz, sizeof(dataType) * devNy * gridNxz, cudaMemcpyDeviceToDevice, 0);
								MHDFilterResizeM << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPb_midr[i]);
								MHDAddMode << <MRK4GridSize, MRK4BlockSize, 0, 0 >> > (dPb_midr[i], dPb_midl[i]);
							}
						}

					}

				}

				//for (int i = 0; i < devNums; i++) {
				//	cudaSetDevice(localRank * devNums + i);
				//	if constexpr (std::is_same_v<ifIon, trueType> && std::is_same_v<ifIonSlowing, trueType>)
				//		PICSlowingDown<ratioDt, Ion, IonType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], Ion_keys_in[i], Ion_values_in[i], rand_keys_in[i], rand_values_in[i], randStates[i]);
				//	if constexpr (std::is_same_v<ifAlpha, trueType> && std::is_same_v<ifAlphaSlowing, trueType>)
				//		PICSlowingDown<ratioDt, Alpha, AlphaType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], Alpha_keys_in[i], Alpha_values_in[i], rand_keys_in[i], rand_values_in[i], randStates[i]);
				//	if constexpr (std::is_same_v<ifBeam, trueType> && std::is_same_v<ifBeamSlowing, trueType>)
				//		PICSlowingDown<ratioDt, Beam, BeamType> << <PICGridSize, PICBlockSize, 0, 0 >> > (pic1d[i], pic2d[i], Beam_keys_in[i], Beam_values_in[i], rand_keys_in[i], rand_values_in[i], randStates[i]);
				//}

				ncclGroupStart();
				for (int i = 0; i < devNums; i++) {

					if constexpr (std::is_same_v<ifIon, trueType>) {
						ncclSend(dPi_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						ncclRecv(dPi_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
						ncclSend(dPi_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
						ncclRecv(dPi_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					}

					if constexpr (std::is_same_v<ifAlpha, trueType>) {
						ncclSend(dPa_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						ncclRecv(dPa_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
						ncclSend(dPa_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
						ncclRecv(dPa_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					}

					if constexpr (std::is_same_v<ifBeam, trueType>) {
						ncclSend(dPb_midl[i] + ncclLeftSend, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
						ncclRecv(dPb_midl[i] + ncclRightRecv, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
						ncclSend(dPb_midl[i] + ncclRightSend, gridGhost * gridNxz, ncclType, ncclRightNei[i], comms[i], 0);
						ncclRecv(dPb_midl[i] + ncclLeftRecv, gridGhost * gridNxz, ncclType, ncclLeftNei[i], comms[i], 0);
					}

				}
				ncclGroupEnd();

			}

		}

		if constexpr (std::is_same_v<ifIon, trueType> || std::is_same_v<ifAlpha, trueType> || std::is_same_v<ifBeam, trueType>) {

			if constexpr (std::is_same_v<ifIon, trueType>) {
				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					cub::DeviceSegmentedRadixSort::SortPairs(Ion_storage[i], Ion_storage_bytes[i], Ion_keys_in[i], Ion_keys_out[i], Ion_values_in[i], Ion_values_out[i], picDev * 7, 7, Ion_offsets[i], Ion_offsets[i] + 1);
					cudaMemcpyAsync(Ion_keys_in[i], Ion_keys_out[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToDevice);
					cudaMemcpyAsync(Ion_values_in[i], Ion_values_out[i], sizeof(dataType) * picDev * 7, cudaMemcpyDeviceToDevice);
				}
			}
			if constexpr (std::is_same_v<ifAlpha, trueType>) {
				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					cub::DeviceSegmentedRadixSort::SortPairs(Alpha_storage[i], Alpha_storage_bytes[i], Alpha_keys_in[i], Alpha_keys_out[i], Alpha_values_in[i], Alpha_values_out[i], picDev * 7, 7, Alpha_offsets[i], Alpha_offsets[i] + 1);
					cudaMemcpyAsync(Alpha_keys_in[i], Alpha_keys_out[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToDevice);
					cudaMemcpyAsync(Alpha_values_in[i], Alpha_values_out[i], sizeof(dataType) * picDev * 7, cudaMemcpyDeviceToDevice);
				}
			}
			if constexpr (std::is_same_v<ifBeam, trueType>) {
				for (int i = 0; i < devNums; i++) {
					cudaSetDevice(localRank * devNums + i);
					cub::DeviceSegmentedRadixSort::SortPairs(Beam_storage[i], Beam_storage_bytes[i], Beam_keys_in[i], Beam_keys_out[i], Beam_values_in[i], Beam_values_out[i], picDev * 7, 7, Beam_offsets[i], Beam_offsets[i] + 1);
					cudaMemcpyAsync(Beam_keys_in[i], Beam_keys_out[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToDevice);
					cudaMemcpyAsync(Beam_values_in[i], Beam_values_out[i], sizeof(dataType) * picDev * 7, cudaMemcpyDeviceToDevice);
				}
			}

		}

	}

	if (myRank == 0) {
		std::cout << BOLDGREEN << 100 << "%" << RESET << std::endl;
		std::cout << std::endl;
	}

	for (int i = 0; i < devNums; i++) {
		cudaSetDevice(localRank * devNums + i);
		cudaEventRecord(end[i]);
		cudaEventSynchronize(end[i]);
		cudaEventElapsedTime(&time[i], start[i], end[i]);
		CUDACHECK(cudaGetLastError());
	}

	if (myRank == 0) {
		for (int i = 1; i < devNums; i++)
			time[0] += time[i];
		time[0] /= devNums;
		if (time[0] > 1000)
			std::cout << BOLDGREEN << "Time used: " << std::setprecision(10) << time[0] / 1000 << "s." << RESET << std::endl;
		else
			std::cout << BOLDGREEN << "Time used: " << std::setprecision(10) << time[0] << "ms." << RESET << std::endl;
		std::cout << std::endl;
	}

	if (myRank == 0)
		std::cout << BOLDYELLOW << "Start: Output results of this simulation." << RESET << std::endl;

	NCCLCHECK(ncclGroupStart());
	for (int i = 0; i < devNums; i++) {
		NCCLCHECK(ncclAllReduce(d_IonDensity[i], d_IonDensity[i], h_IonDensity.size(), ncclType, ncclSum, comms[i], 0));
		NCCLCHECK(ncclAllReduce(d_AlphaDensity[i], d_AlphaDensity[i], h_AlphaDensity.size(), ncclType, ncclSum, comms[i], 0));
		NCCLCHECK(ncclAllReduce(d_BeamDensity[i], d_BeamDensity[i], h_BeamDensity.size(), ncclType, ncclSum, comms[i], 0));
		NCCLCHECK(ncclAllReduce(d_IonDiffusivity[i], d_IonDiffusivity[i], h_IonDiffusivity.size(), ncclType, ncclSum, comms[i], 0));
		NCCLCHECK(ncclAllReduce(d_AlphaDiffusivity[i], d_AlphaDiffusivity[i], h_AlphaDiffusivity.size(), ncclType, ncclSum, comms[i], 0));
		NCCLCHECK(ncclAllReduce(d_BeamDiffusivity[i], d_BeamDiffusivity[i], h_BeamDiffusivity.size(), ncclType, ncclSum, comms[i], 0));
		NCCLCHECK(ncclAllGather(w_midl[i] + gridGhost * gridNxz, cuGMEC.d_w[i], devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(A_midl[i] + gridGhost * gridNxz, cuGMEC.d_A[i], devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(dNe_midl[i] + gridGhost * gridNxz, cuGMEC.d_dNe[i], devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(dTe_midl[i] + gridGhost * gridNxz, cuGMEC.d_dTe[i], devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(Phi_midl[i] + gridGhost * gridNxz, cuGMEC.d_Phi[i], devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(dJpB_midl[i] + gridGhost * gridNxz, cuGMEC.d_dJpB[i], devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(dPe_midl[i] + gridGhost * gridNxz, cuGMEC.d_dPe[i], devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(dPi_midl[i] + gridGhost * gridNxz, cuGMEC.d_globalPi[i] + gridGhost * gridNxz, devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(dPa_midl[i] + gridGhost * gridNxz, cuGMEC.d_globalPa[i] + gridGhost * gridNxz, devNy * gridNxz, ncclType, comms[i], 0));
		NCCLCHECK(ncclAllGather(dPb_midl[i] + gridGhost * gridNxz, cuGMEC.d_globalPb[i] + gridGhost * gridNxz, devNy * gridNxz, ncclType, comms[i], 0));
	}
	NCCLCHECK(ncclGroupEnd());

	CUDACHECK(cudaSetDevice(localRank * devNums));
	CUDACHECK(cudaMemcpy(h_IonDensity.data(), d_IonDensity[0], sizeof(dataType) * h_IonDensity.size(), cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(h_AlphaDensity.data(), d_AlphaDensity[0], sizeof(dataType) * h_AlphaDensity.size(), cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(h_BeamDensity.data(), d_BeamDensity[0], sizeof(dataType) * h_BeamDensity.size(), cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(h_IonDiffusivity.data(), d_IonDiffusivity[0], sizeof(dataType) * h_IonDiffusivity.size(), cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(h_AlphaDiffusivity.data(), d_AlphaDiffusivity[0], sizeof(dataType) * h_AlphaDiffusivity.size(), cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(h_BeamDiffusivity.data(), d_BeamDiffusivity[0], sizeof(dataType) * h_BeamDiffusivity.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputw, trueType>)
		CUDACHECK(cudaMemcpy(h_totalw.data(), d_totalw[0], sizeof(dataType) * h_totalw.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputA, trueType>)
		CUDACHECK(cudaMemcpy(h_totalA.data(), d_totalA[0], sizeof(dataType) * h_totalA.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputdNe, trueType>)
		CUDACHECK(cudaMemcpy(h_totaldNe.data(), d_totaldNe[0], sizeof(dataType) * h_totaldNe.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputdTe, trueType>)
		CUDACHECK(cudaMemcpy(h_totaldTe.data(), d_totaldTe[0], sizeof(dataType) * h_totaldTe.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputPhi, trueType>)
		CUDACHECK(cudaMemcpy(h_totalPhi.data(), d_totalPhi[0], sizeof(dataType) * h_totalPhi.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputdPi, trueType>)
		CUDACHECK(cudaMemcpy(h_totaldPi.data(), d_totaldPi[0], sizeof(dataType) * h_totaldPi.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputdPa, trueType>)
		CUDACHECK(cudaMemcpy(h_totaldPa.data(), d_totaldPa[0], sizeof(dataType) * h_totaldPa.size(), cudaMemcpyDeviceToHost));
	if constexpr (std::is_same_v<ifOutputdPb, trueType>)
		CUDACHECK(cudaMemcpy(h_totaldPb.data(), d_totaldPb[0], sizeof(dataType) * h_totaldPb.size(), cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_w[0][0] + gridGhost * gridNxz, cuGMEC.d_w[0], sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_A[0][0] + gridGhost * gridNxz, cuGMEC.d_A[0], sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_dNe[0][0] + gridGhost * gridNxz, cuGMEC.d_dNe[0], sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_dTe[0][0] + gridGhost * gridNxz, cuGMEC.d_dTe[0], sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_Phi[0][0] + gridGhost * gridNxz, cuGMEC.d_Phi[0], sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_dJpB[0][0] + gridGhost * gridNxz, cuGMEC.d_dJpB[0], sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_dPe[0][0] + gridGhost * gridNxz, cuGMEC.d_dPe[0], sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_globalPi[0][0] + gridGhost * gridNxz, cuGMEC.d_globalPi[0] + gridGhost * gridNxz, sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_globalPa[0][0] + gridGhost * gridNxz, cuGMEC.d_globalPa[0] + gridGhost * gridNxz, sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));
	CUDACHECK(cudaMemcpy(cuGMEC.h_globalPb[0][0] + gridGhost * gridNxz, cuGMEC.d_globalPb[0] + gridGhost * gridNxz, sizeof(dataType) * gridNy * gridNxz, cudaMemcpyDeviceToHost));

	if constexpr (std::is_same_v<ifIon, trueType>) {

		for (int i = 0; i < devNums; i++) {
			cudaSetDevice(localRank * devNums + i);
			if constexpr (std::is_same_v<ifIon, trueType>) {
				CUDACHECK(cudaMemcpy(cuGMEC.h_Ion_keys[i], Ion_keys_in[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToHost));
				CUDACHECK(cudaMemcpy(cuGMEC.h_Ion_values[i], Ion_values_in[i], sizeof(dataType) * picDev * 7, cudaMemcpyDeviceToHost));
			}
		}

		std::ofstream output;
		std::string fileName;

		fileName = "./result/IonConst_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(&cuGMEC.IonConst), sizeof(dataType));
		output.close();

		fileName = "./result/IonOffsets_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Ion_offsets[0]), sizeof(int) * devNums * 8);
		output.close();

		fileName = "./result/IonKeys_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Ion_keys[0]), sizeof(int) * devNums * picDev * 7);
		output.close();

		fileName = "./result/IonValues_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Ion_values[0]), sizeof(dataType) * devNums * picDev * 7);
		output.close();

	}
	if constexpr (std::is_same_v<ifAlpha, trueType>) {

		for (int i = 0; i < devNums; i++) {
			cudaSetDevice(localRank * devNums + i);
			if constexpr (std::is_same_v<ifAlpha, trueType>) {
				CUDACHECK(cudaMemcpy(cuGMEC.h_Alpha_keys[i], Alpha_keys_in[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToHost));
				CUDACHECK(cudaMemcpy(cuGMEC.h_Alpha_values[i], Alpha_values_in[i], sizeof(dataType) * picDev * 7, cudaMemcpyDeviceToHost));
			}
		}

		std::ofstream output;
		std::string fileName;

		fileName = "./result/AlphaConst_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(&cuGMEC.AlphaConst), sizeof(dataType));
		output.close();

		fileName = "./result/AlphaOffsets_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Alpha_offsets[0]), sizeof(int) * devNums * 8);
		output.close();

		fileName = "./result/AlphaKeys_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Alpha_keys[0]), sizeof(int) * devNums * picDev * 7);
		output.close();

		fileName = "./result/AlphaValues_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Alpha_values[0]), sizeof(dataType) * devNums * picDev * 7);
		output.close();

	}
	if constexpr (std::is_same_v<ifBeam, trueType>) {

		for (int i = 0; i < devNums; i++) {
			cudaSetDevice(localRank * devNums + i);
			if constexpr (std::is_same_v<ifBeam, trueType>) {
				CUDACHECK(cudaMemcpy(cuGMEC.h_Beam_keys[i], Beam_keys_in[i], sizeof(int) * picDev * 7, cudaMemcpyDeviceToHost));
				CUDACHECK(cudaMemcpy(cuGMEC.h_Beam_values[i], Beam_values_in[i], sizeof(dataType) * picDev * 7, cudaMemcpyDeviceToHost));
			}
		}

		std::ofstream output;
		std::string fileName;

		fileName = "./result/BeamConst_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(&cuGMEC.BeamConst), sizeof(dataType));
		output.close();

		fileName = "./result/BeamOffsets_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Beam_offsets[0]), sizeof(int) * devNums * 8);
		output.close();

		fileName = "./result/BeamKeys_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Beam_keys[0]), sizeof(int) * devNums * picDev * 7);
		output.close();

		fileName = "./result/BeamValues_" + std::to_string(myRank) + "_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Beam_values[0]), sizeof(dataType) * devNums * picDev * 7);
		output.close();

	}

	if (hostNums == 1) {
		if (devNums == 1) {
			CUDACHECK(cudaSetDevice(0));
			CUDACHECK(cudaMemcpy(h_amplitude.data(), d_amplitude[0], sizeof(dataType) * h_amplitude.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_frequency.data(), d_frequency[0], sizeof(dataType) * h_frequency.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_modeReal.data(), d_modeReal[0], sizeof(dataType) * h_modeReal.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_modeImag.data(), d_modeImag[0], sizeof(dataType) * h_modeImag.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_Epara.data(), d_Epara[0], sizeof(dataType) * h_Epara.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_EparaES.data(), d_EparaES[0], sizeof(dataType) * h_EparaES.size(), cudaMemcpyDeviceToHost));
		}
		else {
			CUDACHECK(cudaSetDevice(devNums / 2));
			CUDACHECK(cudaMemcpy(h_amplitude.data(), d_amplitude[devNums / 2], sizeof(dataType) * h_amplitude.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_frequency.data(), d_frequency[devNums / 2], sizeof(dataType) * h_frequency.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_modeReal.data(), d_modeReal[devNums / 2], sizeof(dataType) * h_modeReal.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_modeImag.data(), d_modeImag[devNums / 2], sizeof(dataType) * h_modeImag.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_Epara.data(), d_Epara[devNums / 2], sizeof(dataType) * h_Epara.size(), cudaMemcpyDeviceToHost));
			CUDACHECK(cudaMemcpy(h_EparaES.data(), d_EparaES[devNums / 2], sizeof(dataType) * h_EparaES.size(), cudaMemcpyDeviceToHost));
		}
	}
	else {
		CUDACHECK(cudaSetDevice(localRank * devNums));
		CUDACHECK(cudaMemcpy(h_amplitude.data(), d_amplitude[0], sizeof(dataType) * h_amplitude.size(), cudaMemcpyDeviceToHost));
		CUDACHECK(cudaMemcpy(h_frequency.data(), d_frequency[0], sizeof(dataType) * h_frequency.size(), cudaMemcpyDeviceToHost));
		CUDACHECK(cudaMemcpy(h_modeReal.data(), d_modeReal[0], sizeof(dataType) * h_modeReal.size(), cudaMemcpyDeviceToHost));
		CUDACHECK(cudaMemcpy(h_modeImag.data(), d_modeImag[0], sizeof(dataType) * h_modeImag.size(), cudaMemcpyDeviceToHost));
		CUDACHECK(cudaMemcpy(h_Epara.data(), d_Epara[0], sizeof(dataType) * h_Epara.size(), cudaMemcpyDeviceToHost));
		CUDACHECK(cudaMemcpy(h_EparaES.data(), d_EparaES[0], sizeof(dataType) * h_EparaES.size(), cudaMemcpyDeviceToHost));
	}

	if (myRank == hostNums / 2) {

		std::ofstream output;

		output.open("./result/amplitude.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_amplitude.data()), sizeof(dataType) * h_amplitude.size());
		output.close();

		output.open("./result/frequency.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_frequency.data()), sizeof(dataType) * h_frequency.size());
		output.close();

		output.open("./result/RealMode.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_modeReal.data()), sizeof(dataType) * h_modeReal.size());
		output.close();

		output.open("./result/ImagMode.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_modeImag.data()), sizeof(dataType) * h_modeImag.size());
		output.close();

		output.open("./result/Epara.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_Epara.data()), sizeof(dataType) * h_Epara.size());
		output.close();

		output.open("./result/EparaES.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_EparaES.data()), sizeof(dataType) * h_EparaES.size());
		output.close();

	}

	if (myRank == 0) {

		std::ofstream output;

		output.open("./result/IonDensity.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_IonDensity.data()), sizeof(dataType) * h_IonDensity.size());
		output.close();

		output.open("./result/AlphaDensity.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_AlphaDensity.data()), sizeof(dataType) * h_AlphaDensity.size());
		output.close();

		output.open("./result/BeamDensity.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_BeamDensity.data()), sizeof(dataType) * h_BeamDensity.size());
		output.close();

		output.open("./result/IonDiffusivity.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_IonDiffusivity.data()), sizeof(dataType) * h_IonDiffusivity.size());
		output.close();

		output.open("./result/AlphaDiffusivity.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_AlphaDiffusivity.data()), sizeof(dataType) * h_AlphaDiffusivity.size());
		output.close();

		output.open("./result/BeamDiffusivity.bin", std::ios::out | std::ios::binary);
		output.write((char*)(h_BeamDiffusivity.data()), sizeof(dataType) * h_BeamDiffusivity.size());
		output.close();

		if constexpr (std::is_same_v<ifOutputw, trueType>) {
			output.open("./result/totalw.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totalw.data()), sizeof(dataType) * h_totalw.size());
			output.close();
		}
		if constexpr (std::is_same_v<ifOutputA, trueType>) {
			output.open("./result/totalA.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totalA.data()), sizeof(dataType) * h_totalA.size());
			output.close();
		}
		if constexpr (std::is_same_v<ifOutputdNe, trueType>) {
			output.open("./result/totaldNe.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totaldNe.data()), sizeof(dataType) * h_totaldNe.size());
			output.close();
		}
		if constexpr (std::is_same_v<ifOutputdTe, trueType>) {
			output.open("./result/totaldTe.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totaldTe.data()), sizeof(dataType) * h_totaldTe.size());
			output.close();
		}
		if constexpr (std::is_same_v<ifOutputPhi, trueType>) {
			output.open("./result/totalPhi.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totalPhi.data()), sizeof(dataType) * h_totalPhi.size());
			output.close();
		}
		if constexpr (std::is_same_v<ifOutputdPi, trueType>) {
			output.open("./result/totaldPi.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totaldPi.data()), sizeof(dataType) * h_totaldPi.size());
			output.close();
		}
		if constexpr (std::is_same_v<ifOutputdPa, trueType>) {
			output.open("./result/totaldPa.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totaldPa.data()), sizeof(dataType) * h_totaldPa.size());
			output.close();
		}
		if constexpr (std::is_same_v<ifOutputdPb, trueType>) {
			output.open("./result/totaldPb.bin", std::ios::out | std::ios::binary);
			output.write((char*)(h_totaldPb.data()), sizeof(dataType) * h_totaldPb.size());
			output.close();
		}

		output.open("./result/w.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_w[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/A.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_A[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/dNe.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_dNe[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/dTe.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_dTe[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/Phi.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_Phi[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/dJpB.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_dJpB[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/dPe.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_dPe[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/dPi.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_globalPi[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/dPa.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_globalPa[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		output.open("./result/dPb.bin", std::ios::out | std::ios::binary);
		output.write((char*)(cuGMEC.h_globalPb[0][0] + gridGhost * gridNxz), sizeof(dataType) * gridNy * gridNxz);
		output.close();

		std::vector<dataType> MHDFinalPerturbation(gridNy * gridNxz * 7);

		for (int j = 0; j < gridNy; j++) {
			for (int i = 0; i < gridNx; i++) {
				for (int k = 0; k < gridNz; k++) {
					MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 0 * gridNy * gridNxz] = cuGMEC.h_w[j + gridGhost][i][k];
					MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 1 * gridNy * gridNxz] = cuGMEC.h_A[j + gridGhost][i][k];
					MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 2 * gridNy * gridNxz] = cuGMEC.h_dNe[j + gridGhost][i][k];
					MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 3 * gridNy * gridNxz] = cuGMEC.h_dTe[j + gridGhost][i][k];
					MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 4 * gridNy * gridNxz] = cuGMEC.h_globalPi[j + gridGhost][i][k];
					MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 5 * gridNy * gridNxz] = cuGMEC.h_globalPa[j + gridGhost][i][k];
					MHDFinalPerturbation[j * gridNxz + i * gridNz + k + 6 * gridNy * gridNxz] = cuGMEC.h_globalPb[j + gridGhost][i][k];
				}
			}
		}

		std::string fileName;
		fileName = "./result/MHDPerturbation_" + std::to_string(continueSteps + totalSteps) + ".bin";
		output.open(fileName.c_str(), std::ios::out | std::ios::binary);
		output.write((char*)(MHDFinalPerturbation.data()), sizeof(dataType) * MHDFinalPerturbation.size());
		output.close();

	}

	if (myRank == 0) {
		std::cout << BOLDGREEN << "Done." << RESET << std::endl;
		std::cout << std::endl;
	}

	for (int i = 0; i < devNums; i++) {

		CUDACHECK(cudaSetDevice(localRank * devNums + i));

		if constexpr (std::is_same_v<dataType, double>)
			CUDACHECK(cudaFree(nFreqd[i]));
		else
			CUDACHECK(cudaFree(nFreqf[i]));

		CUFFTCHECK(cufftDestroy(nPlanR2Cs[i]));
		CUFFTCHECK(cufftDestroy(nPlanC2Rs[i]));

		if constexpr (std::is_same_v<dataType, double>)
			CUDACHECK(cudaFree(mFreqd[i]));
		else
			CUDACHECK(cudaFree(mFreqf[i]));

		CUFFTCHECK(cufftDestroy(mPlanR2Cs[i]));
		CUFFTCHECK(cufftDestroy(mPlanC2Rs[i]));

		CUDACHECK(cudaFree(randStates[i]));

		CUDACHECK(cudaFree(d_amplitude[i]));
		CUDACHECK(cudaFree(d_frequency[i]));
		CUDACHECK(cudaFree(d_modeReal[i]));
		CUDACHECK(cudaFree(d_modeImag[i]));
		CUDACHECK(cudaFree(d_Epara[i]));
		CUDACHECK(cudaFree(d_EparaES[i]));
		CUDACHECK(cudaFree(d_AlphaDensity[i]));
		CUDACHECK(cudaFree(d_IonDensity[i]));
		CUDACHECK(cudaFree(d_BeamDensity[i]));
		CUDACHECK(cudaFree(d_AlphaDiffusivity[i]));
		CUDACHECK(cudaFree(d_IonDiffusivity[i]));
		CUDACHECK(cudaFree(d_BeamDiffusivity[i]));

		if constexpr (std::is_same_v<ifOutputw, trueType>) {
			CUDACHECK(cudaFree(d_eachw[i]));
			CUDACHECK(cudaFree(d_totalw[i]));
		}
		if constexpr (std::is_same_v<ifOutputA, trueType>) {
			CUDACHECK(cudaFree(d_eachA[i]));
			CUDACHECK(cudaFree(d_totalA[i]));
		}
		if constexpr (std::is_same_v<ifOutputdNe, trueType>) {
			CUDACHECK(cudaFree(d_eachdNe[i]));
			CUDACHECK(cudaFree(d_totaldNe[i]));
		}
		if constexpr (std::is_same_v<ifOutputdTe, trueType>) {
			CUDACHECK(cudaFree(d_eachdTe[i]));
			CUDACHECK(cudaFree(d_totaldTe[i]));
		}
		if constexpr (std::is_same_v<ifOutputPhi, trueType>) {
			CUDACHECK(cudaFree(d_eachPhi[i]));
			CUDACHECK(cudaFree(d_totalPhi[i]));
		}
		if constexpr (std::is_same_v<ifOutputdPi, trueType>) {
			CUDACHECK(cudaFree(d_eachdPi[i]));
			CUDACHECK(cudaFree(d_totaldPi[i]));
		}
		if constexpr (std::is_same_v<ifOutputdPa, trueType>) {
			CUDACHECK(cudaFree(d_eachdPa[i]));
			CUDACHECK(cudaFree(d_totaldPa[i]));
		}
		if constexpr (std::is_same_v<ifOutputdPb, trueType>) {
			CUDACHECK(cudaFree(d_eachdPb[i]));
			CUDACHECK(cudaFree(d_totaldPb[i]));
		}

		CUDACHECK(cudaFree(Phi_yxz[i]));
		CUDACHECK(cudaFree(Phi_xzy[i]));
		CUDACHECK(cudaFree(A_yxz[i]));
		CUDACHECK(cudaFree(A_xzy[i]));
		CUDACHECK(cudaFree(dNe_yxz[i]));
		CUDACHECK(cudaFree(dNe_xzy[i]));
		CUDACHECK(cudaFree(dTe_yxz[i]));
		CUDACHECK(cudaFree(dTe_xzy[i]));
		CUDACHECK(cudaFree(dPi_yxz[i]));
		CUDACHECK(cudaFree(dPi_xzy[i]));
		CUDACHECK(cudaFree(dPa_yxz[i]));
		CUDACHECK(cudaFree(dPa_xzy[i]));
		CUDACHECK(cudaFree(dPb_yxz[i]));
		CUDACHECK(cudaFree(dPb_xzy[i]));

		if constexpr (std::is_same_v<ifIon, trueType>)
			CUDACHECK(cudaFree(Ion_storage[i]));
		if constexpr (std::is_same_v<ifAlpha, trueType>)
			CUDACHECK(cudaFree(Alpha_storage[i]));
		if constexpr (std::is_same_v<ifBeam, trueType>)
			CUDACHECK(cudaFree(Beam_storage[i]));

	}

	cuGMEC.releaseDeviceMemory();
	cuGMEC.releaseHostMemory();

	//finalizing NCCL
	for (int i = 0; i < devNums; i++) {
		ncclCommDestroy(comms[i]);
	}

	if (myRank == 0) {
		std::cout << BOLDYELLOW << "Start: Exit the program." << RESET << std::endl;
		std::cout << BOLDGREEN << "Done." << RESET << std::endl;
		std::cout << std::endl;
	}

	//finalizing MPI
	MPICHECK(MPI_Finalize());

	return 0;
}