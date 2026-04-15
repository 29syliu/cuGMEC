
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

   normal_equilibrium.mat
   more_1d.mat
   /more(子文件夹)和其内部的more_3d_i.mat
   plot.mat

   指定outputPath，即在之后由MATLAB生成的cuGMEC输入文件的目标文件夹。
   例如outputPath = 'D:\descFiles\'。
   指定equilibriumName和perturbationName，即在之后由MATLAB生成的cuGMEC输入文件的文件名，分别代表平衡文件名和扰动文件名。
   例如equilibriumName = 'MHDCollocated.bin'; perturbationName = 'MHDPerturbation.bin';
   然后运行代码块(1)。
   
2. 你需要给定电子的密度和温度剖面，背景热离子以及两种快粒子的密度，温度和压强剖面。
   e，i，a，b分别代表电子，热离子，第一种快粒子(通常是阿尔法粒子)和第二种快粒子(通常是束粒子)。
   n，T，P分别代表密度，温度，压强。将密度用10^19/m^3归一化，温度用1keV归一化，压强用1Pa归一化，再执行以下操作。
   你需要确定以解析表达式还是离散数据给定剖面。
   
   如果以解析表达式给定剖面:
   进入代码块(2)，给定背景热离子质量，用质子质量归一化。
   给定多项式每个次数的系数，例如若nepoly = [0.7, -0.3, 0, 0.1];
   则ne(rho) = 0.7*rho^0 + (-0.3)*rho^1 + (0)*rho^2 + (0.1)*rho^3; 其中rho = r/a。
   给定离子的分布类型。
   Type = 1 -> 麦克斯韦分布。此时必须且仅需给定npoly和Tpoly。Ppoly的值无效。
   Type = 2 -> 慢化分布。此时必须且仅需给定npoly和Ppoly。Tpoly的值无效。
   Type = 3 -> 该离子不存在。npoly，Tpoly，Ppoly的值都无效。
   如果是其他解析形式，请计算离散数据，然后选择下面的分支。
   然后运行代码块(2)。

   如果以离散数据给定剖面:
   进入代码块(3)，给定背景热离子质量，用质子质量归一化。
   我们要求将离散数据保存在名为NTP.mat的文件内，并且该文件已置于inputPath下，NTP.mat含有以下变量，
   rhoSample，neSample，TeSample，
   niSample，TiSample，PiSample，
   naSample，TaSample，PaSample，
   nbSample，TbSample，PbSample。
   其中rhoSample的值代表离散的rho，其他Sample的值则代表在对应rho上的其他物理量。
   其中rhoSample的值应当严格递增，以上所有变量的维度都为(1,N)且N > 1。
   给定离子的分布类型。
   Type = 1 -> 麦克斯韦分布。此时必须且仅需给定nSample和TSample。PSample的值无效。
   Type = 2 -> 慢化分布。此时必须且仅需给定nSample和PSample。TSample的值无效。
   Type = 3 -> 该离子不存在。nSample，TSample，PSample的值都无效。
   给定拟合多项式的最高次数。例如neorder = 10;
   然后运行代码块(3)。
   然后运行代码块(4)。

   请注意:
   代码块(2)与(3)仅执行一个，取决于剖面为解析还是离散。
   值无效则代表取任何值都不影响MATLAB生成的文件，且该变量可以不存在。
   执行完代码块(2)或(3)后，会输出所有剖面及其一阶导数的图，你可以由此调整，反复运行代码块(2)或(3), 然后运行代码块(4)。

3. 进入代码块(5)。
   给定环向网格数N_phi。
   给定tube，即模拟1/tube个环向区域。
   给定初始扰动类型以及环向模数:
   PerturType = 0 -> 环向模数只有一个，为leftN=rightN=tube。
   PerturType = 1 -> 环向模数有多个，为leftN:tube:rightN，且leftN和rightN为tube的正整数倍。
   扰动为高斯分布，给定最大位置在radialIndex，宽度为width，幅度为amplitude。
   然后运行代码块(5)。
   
   请注意:
   执行完代码块(5)后，会输出初始扰动在三个截面上的图，你可以由此调整，反复运行代码块(5)。


当完成1.2.3.后，你会在outputPath下看到生成的两个文件，文件名分别为你定义的equilibriumName和perturbationName。
它们是cuGMEC的两个输入文件，两个文件名也分别是cuGMEC参数中MHDCollocated和MHDPerturbation的值。
在outputPath下还会生成Normalizatio.mat，里面包含B0，L0，VA0，CHI0，CHI1，IonBeta，AlphaBeta，BeamBeta，
在cuGMEC中存在同名参数，它们的值由Normalization.mat决定。接下来请参考下一步文档设置cuGMEC参数并编译运行。

   
%}


%% (1)

clear all
close all
addpath('../lib/BSI')
order = uint64(5);

inputPath = 'D:\CPC\R835\0.05-0.90(256_64_16)\KBM\Beta1.8%\TEST\';
outputPath = inputPath;
equilibriumName = 'MHDCollocated.bin';
perturbationName = 'MHDPerturbation.bin';

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

if ~exist(strcat(inputPath,'collocated.mat'), 'file')

    if ~exist(strcat(inputPath,'more_equilibrium.mat'), 'file')
        morePath = strcat(inputPath,'more\');
        matFiles = dir([morePath, '*.mat']);
        for i = 1:length(matFiles)
            filePath = fullfile(morePath, matFiles(i).name);
            load(filePath);
        end
        run MergeMore2D
        
        clearvars -except order inputPath outputPath equilibriumName perturbationName
        addpath('../lib/BSI')
    end

    load(strcat(inputPath,'more_1d.mat'));
    load(strcat(inputPath,'more_equilibrium.mat'));
    load(strcat(inputPath,'normal_equilibrium.mat'));
    allVars = who;
    varsToExclude = {'order', 'inputPath', 'outputPath','equilibriumName','perturbationName'};
    varsToSave = setdiff(allVars, varsToExclude);
    save(strcat(inputPath,'collocated.mat'), varsToSave{:});

    clearvars -except order inputPath outputPath equilibriumName perturbationName
    addpath('../lib/BSI')

end

MHDstaggered = 0;
load(strcat(inputPath,'collocated.mat'));
N_rho = size(R,1);
N_theta = size(R,2);
ghost = 2;


%% (2)

IonMass = 1.0;
IonType = 1;
AlphaType = 1;
BeamType = 1;

nepoly = [];
Tepoly = [];

nipoly = [];
Tipoly = [];
Pipoly = [];

napoly = [];
Tapoly = [];
Papoly = [];

nbpoly = [];
Tbpoly = [];
Pbpoly = [];

run ShowAnalyticProfile

ni_c = polyval(fliplr(nipoly), 0)*1e19;

run Desc2Boozer2ShiftAlign
run OutputNormalization
run OutputMHDEquilibrium

%% (3)

IonMass = 1.0;
IonType = 1;
AlphaType = 3;
BeamType = 3;

load(strcat(inputPath,'NTP.mat'))

neorder = 10;
Teorder = 10;

niorder = 10;
Tiorder = 10;
Piorder = 10;

naorder = 10;
Taorder = 10;
Paorder = 10;

nborder = 10;
Tborder = 10;
Pborder = 10;

p = polyfit(rhoSample, neSample, neorder);
nepoly = fliplr(p);
p = polyfit(rhoSample, TeSample, Teorder);
Tepoly = fliplr(p);

if IonType~=3
    p = polyfit(rhoSample, niSample, niorder);
    nipoly = fliplr(p);
    if IonType==1
        p = polyfit(rhoSample, TiSample, Tiorder);
        Tipoly = fliplr(p);
    elseif IonType==2
        p = polyfit(rhoSample, PiSample, Piorder);
        Pipoly = fliplr(p);
    end
end

if AlphaType~=3
    p = polyfit(rhoSample, naSample, naorder);
    napoly = fliplr(p);
    if AlphaType==1
        p = polyfit(rhoSample, TaSample, Taorder);
        Tapoly = fliplr(p);
    elseif AlphaType==2
        p = polyfit(rhoSample, PaSample, Paorder);
        Papoly = fliplr(p);
    end
end

if BeamType~=3
    p = polyfit(rhoSample, nbSample, nborder);
    nbpoly = fliplr(p);
    if BeamType==1
        p = polyfit(rhoSample, TbSample, Tborder);
        Tbpoly = fliplr(p);
    elseif BeamType==2
        p = polyfit(rhoSample, PbSample, Pborder);
        Pbpoly = fliplr(p);
    end
end

run ShowNumericalProfile

ni_c = polyval(fliplr(nipoly), 0)*1e19;

%% (4)

run Desc2Boozer2ShiftAlign
run OutputNormalization
run OutputMHDEquilibrium

%% (5)

N_phi = 16;
tube = 10;

PerturType = 1;
leftN = 10;
rightN = 10;

radialIndex = 115;
width = 0.05;
amplitude = 1.0e-5;

%amplitude*exp(-(x-x(radialIndex))^2/(2*width^2))*cos(M*theta - N*zeta);

run OutputMHDPerturbation


