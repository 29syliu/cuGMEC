
%%

%{

进行任何操作前，请仔细阅读这一段注释。
进行任何操作前，请仔细阅读这一段注释。
进行任何操作前，请仔细阅读这一段注释。

下方需要你自定义参数的代码块都已经用数字编号，分别为(1)，(2)，(3)，(4)。
请按照以下步骤进行操作:

1. 进入代码块(1)，在inputPath指定前一步DESC输出数据所在的文件夹，
   假设是D:\descFiles\，则inputPath = 'D:\descFiles\'。
   同时你应当查看这个文件夹下是否已经存在DESC输出的数据。
   例如在D:\descFiles\下应当有:

   standard2D.mat
   refined2D.mat
   plot2D.mat

   指定outputPath，即在之后由MATLAB生成的cuGMEC输入文件的目标文件夹。
   例如outputPath = 'D:\descFiles\'。
   指定equilibriumName，即在之后由MATLAB生成的cuGMEC输入文件的文件名，代表平衡文件名。
   例如equilibriumName = 'MHDCollocated_256_64.bin'。
   然后运行代码块(1)。
   
2. 进入代码块(2)，你需要给定粒子数据。

   在运行之前你需要完成以下步骤，首先将离散数据保存在名为NTP.mat的文件内，
   并且该文件已置于inputPath下，NTP.mat含有以下变量：
   rhoSample，neSample，TeSample，
   niSample，TiSample，PiSample，
   naSample，TaSample，PaSample，
   nbSample，TbSample，PbSample。
   其中rhoSample的值代表离散的rho(归一化环向磁通开根号)，其他Sample的值则代表在对应rho上的其他物理量。
   其中rhoSample的值应当严格递增，以上所有变量的维度都为(1,N)且N > 1。
   
   n代表密度，用10^19/m^3归一化。
   P代表压强，用1Pa进行归一化。
   T代表温度，用1keV归一化(麦克斯韦分布时)。
   T代表决定速度，用1m/s归一化(慢化分布时)。
   
   对于背景电子，需要给定neSample(密度)和TeSample(温度)。
   对于离子(背景热离子，第一种快离子a，第二种快离子b)，需要给定其分布类型，并按此给定T与P。
   Type = 1 -> 麦克斯韦分布：需要给定nSample(密度)和TSample(温度)。PSample(压强)不会被读取，无需给定。
   Type = 2 -> 慢化分布：需要给定nSample(密度)，TSample(决定速度)，PSample(压强)。
   Type = 3 -> 该离子不存在：无需给定n，T，PSample。

   IonMass为背景热离子质量，用质子质量归一化。
   
   给定拟合每一个剖面的多项式最高次数和窗口长度。
   例如neOrder = 12;neWindow = 50;
   代表拟合ne时多项式最高次数为12，且以50为窗口长度，每次向右移动1个点，得到多个窗口，进行分段多项式拟合。
   当剖面比较简单时，可以令Window = size(rhoSample,2)，即只采用一个窗口的全局多项式拟合。
   建议优先采用Window = size(rhoSample,2)，当拟合结果不理想时再调整窗口长度。

   然后运行代码块(2)，会输出所有剖面及其一阶导数图，可以由此调整拟合参数，反复运行代码块(2)。
   然后运行代码块(3)，会输出坐标转换时各个量的验证图。

3. 进入代码块(4)。
   指定perturbationName，即在之后由MATLAB生成的cuGMEC输入文件的文件名，代表扰动文件名。
   给定tube，即模拟1/tube个环向区域。
   给定leftN和rightN，
   环向模数n的范围为：leftN*tube ≤ n ≤ rightN*tube，且间隔为tube。
   例如tube = 6，leftN = 1，rightN = 1，代表在1/6环向区域施加n=6的初始扰动。
   例如tube = 6，leftN = 3，rightN = 3，代表在1/6环向区域施加n=18的初始扰动。
   例如tube = 6，leftN = 0，rightN = 4，代表在1/6环向区域施加n=0，6，12，18，24的初始扰动。
   给定环向网格数N_phi。
   建议N_phi ≥ rightN*8，例如N_phi = rightN*16。
   扰动为高斯分布，给定最大位置在radialIndex，宽度为width，幅度为amplitude。
   然后运行代码块(4)。会输出初始扰动在三个截面上的图，你可以由此调整，反复运行代码块(4)。

   

当完成1.2.3.后，你会在outputPath下看到生成的两个文件，文件名分别为你定义的equilibriumName和perturbationName。
它们是cuGMEC的两个输入文件，两个文件名也分别是cuGMEC参数中MHDCollocated和MHDPerturbation的值。
在outputPath下还会生成Normalization2D.mat，里面包含B0，L0，VA0，RHO0，RHO1，PSITMAX，IonBeta，AlphaBeta，BeamBeta，
在cuGMEC中存在同名参数，它们的值由Normalization2D.mat决定。接下来请参考下一步文档设置cuGMEC参数并编译运行。

   
%}


%% (1)

clear all
close all
addpath('../lib/BSI')

inputPath = 'C:\Users\ALFVEN\Desktop\ITER\0.08-0.98(256 32)\';
outputPath = 'C:\Users\ALFVEN\Desktop\ITER\0.08-0.98(256 32)\test\';
equilibriumName = 'MHDCollocated_256_16.bin';

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

load(strcat(inputPath,'standard2D.mat'));
load(strcat(inputPath,'refined2D.mat'));

%% (2)

load(strcat(inputPath,'NTP.mat'))

naSample = naSample_pdf;
nbSample = nbSample_pdf;

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

run ShowProfile1D

%% (3)

run DESC2PEST2SFA2D
run OutputNormalization2D
run OutputEquilibrium2D

%% (4)

perturbationName = 'MHDPerturbation_256_16_16_30_1.bin';

tube = 30;
leftN = 1;
rightN = 7;
N_phi = 16;

radialIndex = 80;
width = 0.04;
amplitude = 1.0e-9;

run OutputMHDPerturbation3D

%%

for tube = 30:2:40
    % 动态生成文件名，将 tube 值填入对应位置
    perturbationName = sprintf('MHDPerturbation_512_64_16_%d_1_160.bin', tube);
    
    % 其他参数固定（若已在外部定义，可省略赋值）
    leftN = 1;
    rightN = 1;
    N_phi = 16;
    radialIndex = 160;
    width = 0.03;
    amplitude = 2.5e-10;
    
    % 运行主脚本
    run OutputMHDPerturbation3D
end