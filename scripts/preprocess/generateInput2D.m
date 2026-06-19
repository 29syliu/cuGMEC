%% cuGMEC 二维输入文件生成脚本

%{

1. 脚本功能
将 DESC 的二维平衡数据和给定的径向粒子剖面转换为 cuGMEC 的平衡输入文件。

2. 代码块 (1) 的参数含义
addpath：MATLAB 预处理脚本目录和 BSI 目录的位置。
inputPath：DESC 输出数据和 NTP.mat 所在目录。当前脚本会从该目录读取 standard2D.mat，NTP.mat。
outputPath：MATLAB 生成 cuGMEC 输入文件的输出目录。
equilibriumName：最终生成的平衡输入二进制文件名，例如 MHDEquilibrium_256_32.bin。

3. 代码块 (2) 的参数含义
NTP.mat 需要包含：
rhoSample，neSample，TeSample
niSample，TiSample，PiSample
naSample，TaSample，PaSample
nbSample，TbSample，PbSample

rhoSample 表示离散的 rho，必须严格递增，所有采样数组都应为 1 x N 向量。
n 为密度，单位 10^19 / m^3。
P 为压强，单位 Pa。
T 在麦克斯韦分布下表示温度，单位 keV。
T 在慢化分布下表示决定速度，单位 m / s。

IonMass：背景热离子质量，以质子质量归一化。
IonType，AlphaType，BeamType：三类离子的分布类型。
Type = 1：麦克斯韦分布，使用 n 和 T。
Type = 2：慢化分布，使用 n，T 和 P。
Type = 3：该粒子组分不存在。

*Order 对应各剖面拟合时使用的多项式最高次数。
*Window 对应各剖面拟合时使用的窗口长度。通常先取 Window = size(rhoSample, 2)，拟合不理想时可以减小窗口长度。

4. 代码块 (3) 的输出
运行后，outputPath 下会生成：
equilibriumName，即 cuGMEC 的平衡输入文件，对应 cuGMEC 参数中的 MHDEquilibrium。
normalization2D.mat，即后续设置 cuGMEC 参数所需的归一化文件。
normalization2D.mat 中包含 B0，L0，VA0，RHO0，RHO1，PSITMAX，IonBeta，AlphaBeta，BeamBeta。

%}

%% (1)

clear all
close all
addpath('C:\Users\Desktop\preprocess')
addpath('C:\Users\Desktop\preprocess\BSI')

inputPath = 'C:\Users\Desktop\';
outputPath = 'C:\Users\Desktop\';
equilibriumName = 'MHDEquilibrium_256_32.bin';

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

load(strcat(inputPath,'standard2D.mat'));

%% (2)

load(strcat(inputPath,'NTP.mat'))

IonMass = 2.5;
IonType = 1;
AlphaType = 2;
BeamType = 2;

neOrder = 12;
neWindow = size(rhoSample,2);
TeOrder = 12;
TeWindow = size(rhoSample,2);

niOrder = 12;
niWindow = size(rhoSample,2);
TiOrder = 12;
TiWindow = size(rhoSample,2);
PiOrder = 0;
PiWindow = 0;

naOrder = 12;
naWindow = size(rhoSample,2);
TaOrder = 12;
TaWindow = size(rhoSample,2);
PaOrder = 12;
PaWindow = size(rhoSample,2);

nbOrder = 6;
nbWindow = 22;
TbOrder = 0;
TbWindow = size(rhoSample,2);
PbOrder = 6;
PbWindow = 22;

run writeProfile1D

%% (3)

run DESC2PEST2SFA2D
run writeNormalization2D
run writeEquilibrium2D
