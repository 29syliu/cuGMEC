%% cuGMEC MHD 诊断统一可视化脚本

%{
功能：
统一可视化 MHD 诊断和输出场量。

输入目录：
将相关输入、诊断和输出文件放在同一个 inputDir 中。

必须包含：
cuGMEC_param.h, NTP.mat，以及以下二者之一：
standard2D.mat, plot2D.mat, normalization2D.mat
standard3D.mat, plot3D.mat, normalization3D.mat

按开关读取：
ifDiagAmplitude, ifDiagFrequency, ifDiagEparallel, ifDiagZFDrive, ifDiagShearing,
ifOutputPhi, ifOutputA, ifOutputdNe, ifOutputdTe,
ifOutputdPi, ifOutputdPa, ifOutputdPb

可能包含：
amplitude.bin, RealMode.bin, ImagMode.bin, frequency.bin, Epara.bin,
EparaES.bin, MaxwellDrive.bin, ReynoldsDrive.bin, ZonalDrive.bin, shearing.bin,
Phi.bin, A.bin, dNe.bin, dTe.bin, dPi.bin, dPa.bin, dPb.bin,
totalPhi.bin, totalA.bin, totaldNe.bin, totaldTe.bin,
totaldPi.bin, totaldPa.bin, totaldPb.bin
%}

%% 用户设置

clear; close all;

inputDir = 'C:\Users\Desktop\test';

paramFile = fullfile(inputDir, 'cuGMEC_param.h');
[standardFile, plotFile, normalizationFile, equilibriumDim] = resolveMHDInputFiles(inputDir);
NTPFile = fullfile(inputDir, 'NTP.mat');
scriptFile = which('visualizeMHD');
scriptDir = fileparts(scriptFile);
repoBsiDir = fullfile(scriptDir, '..', 'preprocess', 'BSI');

%% 读取所有输入

assert(isfolder(inputDir), '缺少输入目录：%s', inputDir);
assert(isfile(paramFile), '缺少参数文件：%s', paramFile);
assert(isfile(standardFile), '缺少 MAT 文件：%s', standardFile);
assert(isfile(plotFile), '缺少 MAT 文件：%s', plotFile);
assert(isfile(normalizationFile), '缺少 MAT 文件：%s', normalizationFile);
assert(isfile(NTPFile), '缺少 MAT 文件：%s', NTPFile);
if isfolder(repoBsiDir)
    addpath(repoBsiDir);
else
    fprintf('[skip] BSI: 目录不存在\n');
end

paramText = fileread(paramFile);
mhdInput = loadMHDInputData(standardFile, plotFile, normalizationFile, NTPFile);
meta = readMHDMetadata(paramText);
meta.equilibriumDim = equilibriumDim;
meta.NFP = readNormalizationNFP(mhdInput);
meta.isTokamak = (meta.NFP == 1);
meta.isStellarator = (meta.NFP ~= 1);
plotPhi = [];
if isfield(mhdInput, 'phi')
    plotPhi = mhdInput.phi;
end
mhdFieldGeom = buildMHDFieldPlotGeometry(mhdInput.q, mhdInput.theta_pest, mhdInput.rho, ...
    mhdInput.qplot, mhdInput.rhoplot, mhdInput.Rplot, mhdInput.Zplot, ...
    meta.yGrid, meta.zGrid, meta.tubes, meta.gridGhost, meta.NFP, plotPhi);
mhdData = readAllMHDDiagnostics(inputDir, meta);
mhdDiagnosticContext = buildMHDDiagnosticContext(meta, mhdInput);

mhdWorkspace = struct();
mhdWorkspace.input = mhdInput;
mhdWorkspace.meta = meta;
mhdWorkspace.data = mhdData;
shifted = [];
aligned = [];

%% 可视化单时间 MHD 场量

mhdFieldOpt = struct( ...
    'enabled', true, ...
    'names', "Phi", ...
    'modeN', meta.physicalNAll(1), ...
    'toroidalAngle', 0.0, ...
    'toroidalIndex', 1, ...
    'colormapIndex', 1);
%{
enabled      : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
toroidalAngle : 托卡马克绘图所在环向角 phi。
toroidalIndex : 仿星器绘图所在环向网格索引。
colormapIndex : 色表序号，可选 1-5。
%}
mhdWorkspace.field = runMHDFieldDiagnostic( ...
    inputDir, meta, mhdInput, mhdFieldGeom, mhdFieldOpt, 'snapshot', 'surface');
[shifted, aligned] = updateFieldAliases(mhdWorkspace.field, shifted, aligned);

%% 可视化单时间 MHD 极向模数

mhdPoloidalModeOpt = struct( ...
    'enabled', true, ...
    'names', "Phi", ...
    'modeN', meta.physicalNAll(1), ...
    'mRange', [1, 20], ...
    'minFFT', 1e-8, ...
    'normalize', true, ...
    'xLim', [0, 1], ...
    'xTicks', 0:0.2:1, ...
    'yLim', [0, 1], ...
    'yTicks', 0:0.25:1);
%{
enabled       : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
mRange        : 绘制的极向模数 m 范围；示例 [1 20]。
minFFT        : 只绘制峰值与全局峰值之比超过该值的 m。
normalize     : 是否用所选 m 范围内的最大 FFT 幅值归一化。
xLim/xTicks   : 横坐标范围和刻度；[] 表示自动。
yLim/yTicks   : 纵坐标范围和刻度；[] 表示自动。
%}
mhdWorkspace.poloidalMode = runMHDFieldDiagnostic( ...
    inputDir, meta, mhdInput, mhdFieldGeom, mhdPoloidalModeOpt, 'snapshot', 'poloidal');

%% 可视化 MHD 幅值（amplitude）

mhdAmplitudeOpt = struct( ...
    'enabled', meta.switch.ifDiagAmplitude, ...
    'timeAxis', 'ta', ...
    'radialIndex', 90, ...
    'modeN', [], ...
    'growthRange', [], ...
    'growthUnit', '1/wa', ...
    'yLim', [-10 0], ...
    'yTicks', -8:2:0, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
timeAxis      : 'ta', 'ms' 或 'steps'。
radialIndex   : 固定径向索引（1 到 gridNx）。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
growthRange   : [] 表示自动范围；仅交互模式使用。
growthUnit    : '1/wa' 或 '1/s'。
yLim          : [] 表示自动范围；示例 [-8 -2]。
yTicks        : [] 表示自动刻度；示例 -8:1:-2。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
mhdWorkspace.amplitude = runMHDAmplitudePlot(mhdData, mhdDiagnosticContext, mhdAmplitudeOpt);

%% 可视化 MHD 短时傅里叶频率

mhdMultipleFrequencyOpt = struct( ...
    'enabled', meta.switch.ifDiagAmplitude, ...
    'timeAxis', 'ms', ...
    'radialIndex', 85, ...
    'modeN', meta.physicalNAll(1), ...
    'windowLength', 1024, ...
    'windowStep', 64, ...
    'nFFT', 15000, ...
    'frequencyRangeHz', [], ...
    'windowType', 'hann', ...
    'removeMean', true, ...
    'interactive', 2);
%{
enabled         : 是否绘图。
timeAxis        : 'ta', 'ms' 或 'steps'。
radialIndex     : 固定径向索引（1 到 gridNx）。
modeN           : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
windowLength    : 短时 FFT 窗口长度，单位为诊断点数。
windowStep      : 短时 FFT 窗口滑动步长，单位为诊断点数。
nFFT            : 时间 FFT 点数；[] 表示 2^nextpow2(windowLength)。
frequencyRangeHz: 选峰频率范围，单位 Hz；示例 [-2e5 2e5]。
windowType      : 'hann' 或 'rect'。
removeMean      : 每个窗口内去除复信号均值。
interactive     : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
mhdWorkspace.multipleFrequency = runMHDMultipleFrequencyPlot( ...
    mhdData, mhdDiagnosticContext, mhdMultipleFrequencyOpt);

%% 可视化 MHD 相位频率

mhdPhaseFrequencyOpt = struct( ...
    'enabled', meta.switch.ifDiagAmplitude, ...
    'timeAxis', 'ms', ...
    'radialIndex', 20, ...
    'modeN', meta.physicalNAll(1), ...
    'phaseStep', 2000, ...
    'smoothWindow', 1, ...
    'amplitudeFloor', 0, ...
    'yLim', [], ...
    'yTicks', [], ...
    'interactive', 0);
%{
enabled       : 是否绘图。
timeAxis      : 'ta', 'ms' 或 'steps'。
radialIndex   : 固定径向索引（1 到 gridNx）。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
phaseStep     : 相位差分间隔，单位为诊断点数。
smoothWindow  : 频率曲线平滑窗口，单位为差分后的点数；1 表示不平滑。
amplitudeFloor: 信号幅度低于该值时对应频率置为 NaN。
yLim          : [] 表示自动范围；示例 [0 2e5]。
yTicks        : [] 表示自动刻度；示例 0:5e4:2e5。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
mhdWorkspace.phaseFrequency = runMHDPhaseFrequencyPlot(mhdData, mhdDiagnosticContext, mhdPhaseFrequencyOpt);

%% 可视化 MHD 等高线频率诊断

mhdContourFrequencyOpt = struct( ...
    'enabled', meta.switch.ifDiagAmplitude, ...
    'modeN', meta.physicalNAll(1), ...
    'timeIndexRange', [], ...
    'frequencyUnit', 'omegaA', ...
    'contourLevels', linspace(0.2, 1, 10), ...
    'showContourLine', true, ...
    'yLim', [0 1], ...
    'yTicks', 0:0.2:0.8);
%{
enabled       : 是否绘图。
modeN         : 真实物理环向模数 n；一次绘制一个 n；示例 6。
timeIndexRange: 参与 FFT 的诊断时间索引范围；[] 表示全部时间；示例 [1000 5000]。
frequencyUnit : 'omegaA' 表示 omega/omega_A；'kHz' 表示 f/kHz。
contourLevels : 归一化 FFT 强度等值级别；示例 linspace(0.2, 1, 10)。
showContourLine: 是否显示等高线边界；false 时只显示填色色块。
yLim          : [] 表示自动范围；单位跟随 frequencyUnit；示例 [0 1]。
yTicks        : [] 表示自动刻度；单位跟随 frequencyUnit；示例 0:0.2:0.8。
%}
mhdWorkspace.contourFrequency = runMHDContourFrequencyPlot( ...
    mhdData, mhdDiagnosticContext, mhdContourFrequencyOpt);

%% 可视化单点 MHD Phi 信号

mhdSingleSignalOpt = struct( ...
    'enabled', meta.switch.ifDiagFrequency, ...
    'timeAxis', 'ta', ...
    'radialIndex', 128, ...
    'logFloor', realmin, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
timeAxis      : 'ta', 'ms' 或 'steps'。
radialIndex   : 固定径向索引（1 到 gridNx）。
logFloor      : 避免 log(0)。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
mhdWorkspace.singleSignal = runMHDSingleSignalPlot(mhdData, mhdDiagnosticContext, mhdSingleSignalOpt);

%% 可视化 MHD Epara / EparaES

mhdEparaOpt = struct( ...
    'enabled', meta.switch.ifDiagEparallel, ...
    'timeIndex', 12000, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
timeIndex     : 固定诊断时间索引（1 到 nDiagTime）。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
mhdWorkspace.epara = runMHDEparaPlot(mhdData, mhdDiagnosticContext, mhdEparaOpt);

%% 可视化 MHD 剪切率（shearing）

mhdShearingOpt = struct( ...
    'enabled', meta.switch.ifDiagShearing, ...
    'plotType', 1, ...
    'timeIndex', 12000, ...
    'radialIndex', 90, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
plotType      : 1 固定 t 的 x 剖面；2 固定 x 的 t 剖面；3 x-t 二维图。
timeIndex     : plotType = 1 时的固定诊断时间索引（1 到 nDiagTime）。
radialIndex   : plotType = 2 时的固定径向索引（1 到 gridNx）。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
mhdWorkspace.shearing = runMHDShearingPlot(mhdData, mhdDiagnosticContext, mhdShearingOpt);

%% 可视化 MHD MaxwellDrive / ReynoldsDrive / ZonalDrive

mhdZonalDriveOpt = struct( ...
    'enabled', meta.switch.ifDiagZFDrive, ...
    'timeIndex', 15000, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
timeIndex     : 固定诊断时间索引（1 到 nDiagTime）。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
mhdWorkspace.zonalDrive = runMHDZonalDrivePlot(mhdData, mhdDiagnosticContext, mhdZonalDriveOpt);

%% 可视化单时间 total MHD 场量

mhdTotalFieldOpt = struct( ...
    'enabled', true, ...
    'names', "Phi", ...
    'timeIndex', meta.nOutputTime, ...
    'modeN', meta.physicalNAll, ...
    'toroidalAngle', 0.0, ...
    'toroidalIndex', 1, ...
    'colormapIndex', 3);
%{
enabled      : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"。
timeIndex     : total 场量输出时间索引（1 到 nOutputTime）。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
toroidalAngle : 托卡马克绘图所在环向角 phi。
toroidalIndex : 仿星器绘图所在环向网格索引。
colormapIndex : 色表序号，可选 1-5。
%}
mhdWorkspace.totalField = runMHDFieldDiagnostic( ...
    mhdData.total, meta, mhdInput, mhdFieldGeom, mhdTotalFieldOpt, 'total', 'surface');
[shifted, aligned] = updateFieldAliases(mhdWorkspace.totalField, shifted, aligned);

%% 可视化单时间 total MHD 极向模数

mhdTotalPoloidalModeOpt = struct( ...
    'enabled', true, ...
    'names', "Phi", ...
    'timeIndex', meta.nOutputTime, ...
    'modeN', meta.physicalNAll, ...
    'mRange', [1, 20], ...
    'minFFT', 1e-8, ...
    'normalize', true, ...
    'xLim', [0, 1], ...
    'xTicks', 0:0.2:1, ...
    'yLim', [0, 1], ...
    'yTicks', 0:0.25:1);
%{
enabled       : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"；对应 totalPhi 等 total 场。
timeIndex     : total 场量输出时间索引（1 到 nOutputTime）。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
mRange        : 绘制的极向模数 m 范围；示例 [1 20]。
minFFT        : 只绘制峰值与全局峰值之比超过该值的 m。
normalize     : 是否用所选 m 范围内的最大 FFT 幅值归一化。
xLim/xTicks   : 横坐标范围和刻度；[] 表示自动。
yLim/yTicks   : 纵坐标范围和刻度；[] 表示自动。
%}
mhdWorkspace.totalPoloidalMode = runMHDFieldDiagnostic( ...
    mhdData.total, meta, mhdInput, mhdFieldGeom, mhdTotalPoloidalModeOpt, 'total', 'poloidal');

%% 局部函数

function [standardFile, plotFile, normalizationFile, equilibriumDim] = resolveMHDInputFiles(inputDir)

    [standardFile, plotFile, normalizationFile] = buildMHDInputFileSet(inputDir, 3);
    if isfile(standardFile) && isfile(plotFile) && isfile(normalizationFile)
        equilibriumDim = 3;
        return;
    end

    [standardFile, plotFile, normalizationFile] = buildMHDInputFileSet(inputDir, 2);
    equilibriumDim = 2;
end

function [standardFile, plotFile, normalizationFile] = buildMHDInputFileSet(inputDir, equilibriumDim)
    dimensionSuffix = sprintf('%dD.mat', equilibriumDim);
    standardFile = fullfile(inputDir, ['standard' dimensionSuffix]);
    plotFile = fullfile(inputDir, ['plot' dimensionSuffix]);
    normalizationFile = fullfile(inputDir, ['normalization' dimensionSuffix]);
end

function NFP = readNormalizationNFP(mhdInput)

    assert(isfield(mhdInput, 'NFP'), '归一化 MAT 文件缺少字段 "NFP"。');
    NFP = double(mhdInput.NFP);
    assert(isscalar(NFP) && isfinite(NFP) && NFP == floor(NFP) && NFP >= 1, ...
        'NFP 必须为正整数。');
end

function [shifted, aligned] = updateFieldAliases(workspace, shifted, aligned)

    fieldNames = fieldnames(workspace);
    for iField = 1:numel(fieldNames)
        fieldResult = workspace.(fieldNames{iField});
        if isstruct(fieldResult) && isfield(fieldResult, 'shifted') && isfield(fieldResult, 'aligned')
            shifted = fieldResult.shifted;
            aligned = fieldResult.aligned;
            return;
        end
    end
end

function data = loadMHDInputData(standardFile, plotFile, normalizationFile, NTPFile)

    data = struct();
    files = {standardFile, plotFile, normalizationFile, NTPFile};
    for iFile = 1:numel(files)
        raw = load(files{iFile});
        fields = fieldnames(raw);
        for iField = 1:numel(fields)
            data.(fields{iField}) = raw.(fields{iField});
        end
    end

    requiredFields = {'q', 'theta_pest', 'rho', 'qplot', 'rhoplot', 'Rplot', 'Zplot', ...
        'B0', 'L0', 'VA0', 'TeSample'};
    for iField = 1:numel(requiredFields)
        assert(isfield(data, requiredFields{iField}), ...
            '输入 MAT 文件缺少字段 "%s"。', requiredFields{iField});
        if isnumeric(data.(requiredFields{iField}))
            data.(requiredFields{iField}) = double(data.(requiredFields{iField}));
        end
    end
end

function meta = readMHDMetadata(paramText)

    integerNames = {'gridNx', 'gridNy', 'gridNz', 'leftN', 'rightN', 'tubes'};
    meta = readNamedMHDParameters(struct(), paramText, integerNames, @readIntParam);
    meta.gridGhost = 2;
    stepNames = {'totalSteps', 'diagSteps', 'outputSteps'};
    meta = readNamedMHDParameters(meta, paramText, stepNames, @readIntParam);
    meta.dt = readFloatParam(paramText, 'dt');
    meta.mhdPrecision = readPrecisionParam(paramText);

    assert(meta.rightN >= meta.leftN, 'rightN 必须大于或等于 leftN。');
    assert(meta.gridNx > 1 && meta.gridNy > 0 && meta.gridNz > 0, ...
        'gridNx/gridNy/gridNz 取值不合法。');
    assert(meta.tubes > 0, 'tubes 必须为正整数。');
    assert(meta.totalSteps >= 0 && meta.diagSteps > 0 && meta.outputSteps > 0, ...
        'totalSteps/diagSteps/outputSteps 取值不合法。');

    switchNames = { ...
        'ifDiagAmplitude', 'ifDiagFrequency', 'ifDiagEparallel', 'ifDiagZFDrive', 'ifDiagShearing', ...
        'ifOutputPhi', 'ifOutputA', 'ifOutputdNe', 'ifOutputdTe', 'ifOutputdPi', 'ifOutputdPa', 'ifOutputdPb'};
    meta.switch = readNamedMHDParameters(struct(), paramText, switchNames, @readSwitchParam);

    meta.nDiagTime = floor(meta.totalSteps / meta.diagSteps) + 1;
    meta.nOutputTime = floor(meta.totalSteps / meta.outputSteps) + 1;
    meta.nMode = meta.rightN - meta.leftN + 1;
    meta.modeIndexAll = meta.leftN:meta.rightN;
    meta.physicalNAll = meta.modeIndexAll * meta.tubes;
    meta.tDiag = (0:meta.nDiagTime - 1) * meta.diagSteps * meta.dt;
    meta.tOutput = (0:meta.nOutputTime - 1) * meta.outputSteps * meta.dt;
    meta.xGrid = linspace(0.0, 1.0, meta.gridNx);
    meta.yGrid = ((0:meta.gridNy - 1) + 0.5) / meta.gridNy * 2 * pi - pi;
    meta.zGrid = ((0:meta.gridNz - 1) + 0.5) / meta.gridNz * 2 * pi / meta.tubes - pi / meta.tubes;
end

function values = readNamedMHDParameters(values, paramText, names, readerFcn)
    for iName = 1:numel(names)
        name = names{iName};
        values.(name) = readerFcn(paramText, name);
    end
end

function data = readAllMHDDiagnostics(inputDir, meta)

    diagnosticSpecs = { ...
        'amplitude', 'ifDiagAmplitude', 'amplitude.bin', 'mode'; ...
        'RealMode', 'ifDiagAmplitude', 'RealMode.bin', 'mode'; ...
        'ImagMode', 'ifDiagAmplitude', 'ImagMode.bin', 'mode'; ...
        'frequency', 'ifDiagFrequency', 'frequency.bin', 'radial'; ...
        'Epara', 'ifDiagEparallel', 'Epara.bin', 'radial'; ...
        'EparaES', 'ifDiagEparallel', 'EparaES.bin', 'radial'; ...
        'Shearing', 'ifDiagShearing', 'shearing.bin', 'radial'; ...
        'MaxwellDrive', 'ifDiagZFDrive', 'MaxwellDrive.bin', 'radial'; ...
        'ReynoldsDrive', 'ifDiagZFDrive', 'ReynoldsDrive.bin', 'radial'; ...
        'ZonalDrive', 'ifDiagZFDrive', 'ZonalDrive.bin', 'radial'; ...
        'Phi', 'ifOutputPhi', 'totalPhi.bin', 'output'; ...
        'A', 'ifOutputA', 'totalA.bin', 'output'; ...
        'dNe', 'ifOutputdNe', 'totaldNe.bin', 'output'; ...
        'dTe', 'ifOutputdTe', 'totaldTe.bin', 'output'; ...
        'dPi', 'ifOutputdPi', 'totaldPi.bin', 'output'; ...
        'dPa', 'ifOutputdPa', 'totaldPa.bin', 'output'; ...
        'dPb', 'ifOutputdPb', 'totaldPb.bin', 'output'};

    ordinaryNames = diagnosticSpecs(1:10, 1);
    data = cell2struct(repmat({[]}, size(ordinaryNames)), ordinaryNames, 1);
    data.total = struct();
    readerFcns.mode = @(filePath) readModeDiagnosticAsTXN( ...
        filePath, meta.mhdPrecision, meta.nDiagTime, meta.gridNx, meta.nMode);
    readerFcns.radial = @(filePath) readRadialDiagnosticAsTX( ...
        filePath, meta.mhdPrecision, meta.nDiagTime, meta.gridNx);
    readerFcns.output = @(filePath) readOutputAsTYXZ( ...
        filePath, meta.mhdPrecision, meta.nOutputTime, meta.gridNy, meta.gridNx, meta.gridNz);

    for iSpec = 1:size(diagnosticSpecs, 1)
        fieldName = diagnosticSpecs{iSpec, 1};
        switchName = diagnosticSpecs{iSpec, 2};
        fileName = diagnosticSpecs{iSpec, 3};
        dataKind = diagnosticSpecs{iSpec, 4};
        isOutput = strcmp(dataKind, 'output');
        logName = fieldName;
        if isOutput
            logName = ['total' fieldName];
        end

        if ~meta.switch.(switchName)
            logSkipped(logName, ['开关 ' switchName ' 为 false']);
            continue;
        end

        filePath = fullfile(inputDir, fileName);
        if ~isfile(filePath)
            logSkipped(fileName, '文件不存在');
            continue;
        end

        readerFcn = readerFcns.(dataKind);
        loadedData = readerFcn(filePath);

        if isOutput
            data.total.(fieldName) = loadedData;
        else
            data.(fieldName) = loadedData;
        end
        logLoaded(logName, loadedData);
    end
end

function runInteractivePlot(interactiveMode, staticPlotFcn, interactivePlotFcn)

    interactiveMode = requireIntegerInRange(interactiveMode, 0, 2, 'interactive');
    switch interactiveMode
        case 0
            staticPlotFcn();
        case {1, 2}
            if isempty(interactivePlotFcn)
                staticPlotFcn();
            else
                interactivePlotFcn(interactiveMode == 2);
            end
    end
end

function [workspace, ready] = beginMHDDiagnostic(opt, data, requiredFields, plotName, missingReason)

    workspace = struct('options', opt);
    ready = false;
    if ~requireLogicalScalar(opt.enabled, 'enabled')
        logSkipped(plotName, '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(data, requiredFields)
        logSkipped(plotName, missingReason);
        return;
    end
    ready = true;
end

function value = getOptionValue(opt, fieldName, defaultValue)

    value = defaultValue;
    if isfield(opt, fieldName) && ~isempty(opt.(fieldName))
        value = opt.(fieldName);
    end
end

function value = requireLogicalScalar(value, label)

    isLogicalLike = (islogical(value) || isnumeric(value)) && isreal(value) && isscalar(value);
    assert(isLogicalLike && isfinite(double(value)) && (double(value) == 0 || double(value) == 1), ...
        '%s 必须是逻辑标量 true/false（或 0/1）。', label);
    value = logical(value);
end

function value = requireFiniteScalarInRange(value, minValue, maxValue, label)

    assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value), ...
        '%s 必须是有限数值标量。', label);
    value = double(value);
    assert(value >= minValue && value <= maxValue, ...
        '%s 必须位于 [%g, %g]。', label, minValue, maxValue);
end

function value = requireIntegerInRange(value, minValue, maxValue, label)

    value = requireFiniteScalarInRange(value, minValue, maxValue, label);
    assert(value == floor(value), '%s 必须是整数。', label);
end

function workspace = runMHDFieldDiagnostic( ...
    inputSource, meta, mhdInput, mhdFieldGeom, opt, sourceKind, diagnosticKind)

    workspace = struct('options', opt);
    switch sourceKind
        case 'snapshot'
            useTotalField = false;
        case 'total'
            useTotalField = true;
        otherwise
            error('未知 MHD 场量来源：%s。', sourceKind);
    end
    switch diagnosticKind
        case 'surface'
            isSurface = true;
        case 'poloidal'
            isSurface = false;
        otherwise
            error('未知 MHD 场量诊断类型：%s。', diagnosticKind);
    end
    if isSurface && useTotalField
        plotName = 'total field plot';
        modeContext = 'total 场量';
    elseif isSurface
        plotName = 'MHD field plot';
        modeContext = 'MHD 场量';
    else
        plotName = 'MHD poloidal mode plot';
        modeContext = 'MHD 极向模数';
    end
    if ~requireLogicalScalar(opt.enabled, 'enabled')
        logSkipped(plotName, '绘图开关为 false');
        return;
    end
    if isSurface && ~hasBSpline()
        logSkipped(plotName, '缺少 bspline');
        return;
    end

    fieldSource = struct();
    fieldSource.kind = sourceKind;
    fieldSource.input = inputSource;
    fieldSource.meta = meta;
    fieldSource.options = opt;
    fieldNames = parseMHDFieldNames(opt.names);
    for iField = 1:numel(fieldNames)
        [fieldZXY, sourceInfo, ready] = readMHDFieldSource(fieldSource, fieldNames(iField));
        if ~ready
            continue;
        end

        if isSurface && ~useTotalField
            fieldYXZ = mhdFieldZXYToYXZ(fieldZXY);
            logLoaded(sourceInfo.fieldName, fieldYXZ);
        end
        [fieldModeIndex, fieldPhysicalN] = parseModeN( ...
            opt.modeN, meta.modeIndexAll, meta.physicalNAll, modeContext);
        fieldZXY = filterMHDToroidalModes(fieldZXY, fieldModeIndex);

        if isSurface
            plotContext = '';
            if useTotalField
                plotContext = sprintf('timeIndex=%d, t_a=%.6g', ...
                    sourceInfo.timeIndex, meta.tOutput(sourceInfo.timeIndex));
            end
            [shifted, aligned] = plotMHDFieldOnPoloidalPlane(fieldZXY, sourceInfo.plotName, ...
                mhdFieldGeom, opt, fieldPhysicalN, plotContext);
            fieldResult = struct();
            if useTotalField
                fieldResult.timeIndex = sourceInfo.timeIndex;
            else
                % 保留历史 workspace 字段名；其中数组的实际维序为 YXZ。
                fieldResult.fieldNYXZ = fieldYXZ;
            end
            fieldResult.shifted = shifted;
            fieldResult.aligned = aligned;
            fieldResult.physicalN = fieldPhysicalN;
            workspace.(sourceInfo.fieldName) = fieldResult;
            continue;
        end

        [plotField, fftField, xVec, xLabelText, selectedM, yData, rawAmplitude] = ...
            calculateMHDPoloidalModeFFT(fieldZXY, mhdInput, mhdFieldGeom, opt);
        if isempty(selectedM)
            logSkipped([sourceInfo.plotName ' poloidal mode plot'], '没有可绘制的极向模数');
            continue;
        end

        if useNormalizedPoloidalMode(opt)
            yLabelText = '$|U_m|/\max |U_m|$';
        else
            yLabelText = '$|\delta\phi_m|$';
        end
        if useTotalField
            statusText = sprintf('%s poloidal FFT: timeIndex=%d, t_a=%.6g, n=%s, m=%s', ...
                sourceInfo.plotName, sourceInfo.timeIndex, meta.tOutput(sourceInfo.timeIndex), ...
                formatNumberList(fieldPhysicalN), formatNumberList(selectedM));
        else
            statusText = sprintf('%s poloidal FFT: n=%s, m=%s', ...
                sourceInfo.plotName, formatNumberList(fieldPhysicalN), formatNumberList(selectedM));
        end
        plotData = buildLinePlotData(xVec, yData, xLabelText, yLabelText, ...
            sprintf('%s poloidal FFT', sourceInfo.plotName), compose('$m=%d$', selectedM), statusText);
        drawPlot(plotData, opt, @renderLinePlot);
        fprintf('[plot] %s\n', statusText);

        workspace.(sourceInfo.fieldName) = struct( ...
            'poloidalPlane', plotField, ...
            'fft', fftField, ...
            'rawAmplitude', rawAmplitude, ...
            'amplitude', yData, ...
            'x', xVec, ...
            'm', selectedM, ...
            'physicalN', fieldPhysicalN);
        if useTotalField
            workspace.(sourceInfo.fieldName).timeIndex = sourceInfo.timeIndex;
        end
    end
end

function [fieldZXY, sourceInfo, ready] = readMHDFieldSource(fieldSource, fieldName)

    fieldNameText = char(fieldName);
    fieldZXY = [];
    sourceInfo = struct();
    sourceInfo.fieldName = fieldNameText;
    sourceInfo.plotName = fieldNameText;
    sourceInfo.timeIndex = [];
    ready = false;
    switch fieldSource.kind
        case 'total'
            sourceInfo.plotName = ['total' fieldNameText];
            if ~isfield(fieldSource.input, fieldNameText) || isempty(fieldSource.input.(fieldNameText))
                logSkipped(sourceInfo.plotName, '未读取该场量');
                return;
            end
            totalData = fieldSource.input.(fieldNameText);
            sourceInfo.timeIndex = parseTimeIndex(fieldSource.options.timeIndex, size(totalData, 1));
            fieldZXY = totalMHDFieldTimeSliceAsZXY(totalData, sourceInfo.timeIndex);
        case 'snapshot'
            fieldFile = fullfile(fieldSource.input, [fieldNameText '.bin']);
            if ~isfile(fieldFile)
                logSkipped(fieldNameText, '文件不存在');
                return;
            end
            fieldZXY = readMHDFieldAsZXY(fieldSource.input, fieldName, fieldSource.meta.mhdPrecision, ...
                fieldSource.meta.gridNy, fieldSource.meta.gridNx, fieldSource.meta.gridNz);
        otherwise
            error('未知 MHD 场量来源：%s。', fieldSource.kind);
    end
    ready = true;
end

function fieldNames = parseMHDFieldNames(names)

    allowedNames = ["Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"];
    fieldNames = string(names);
    assert(~isempty(fieldNames) && isvector(fieldNames) && all(ismember(fieldNames(:), allowedNames)), ...
        'names 必须至少包含一个支持的场量：Phi、A、dNe、dTe、dPi、dPa 或 dPb。');
    fieldNames = fieldNames(:).';
end

function [plotField, fftField, xVec, xLabelText, selectedM, yData, rawAmplitude] = ...
    calculateMHDPoloidalModeFFT(fieldZXY, mhdInput, geom, opt)

    validatePoloidalModeOptions(opt);
    pestRefined = prepareMHDPESTField(fieldZXY, geom, opt, geom.nPlotTheta);

    plotField = squeeze(pestRefined(1, :, :));
    fftField = fft(pestRefined, [], 3);
    nTheta = size(pestRefined, 3);
    amplitude = squeeze(max(abs(fftField), [], 1)) * (2 / nTheta);
    amplitude(:, 1) = amplitude(:, 1) / 2;
    if mod(nTheta, 2) == 0
        amplitude(:, nTheta / 2 + 1) = amplitude(:, nTheta / 2 + 1) / 2;
    end
    amplitude(~isfinite(amplitude)) = 0;

    candidateM = poloidalModeCandidates(opt, size(amplitude, 2));
    selectedM = selectPoloidalModeNumbers(amplitude, candidateM, opt);
    if isempty(selectedM)
        xVec = [];
        xLabelText = '';
        yData = [];
        rawAmplitude = [];
        return;
    end

    rawAmplitude = amplitude(:, selectedM + 1);
    rangeAmplitude = amplitude(:, candidateM + 1);
    maxAmplitude = max(rangeAmplitude(:));
    if useNormalizedPoloidalMode(opt) && maxAmplitude > 0
        yData = rawAmplitude ./ maxAmplitude;
    else
        yData = rawAmplitude;
    end

    xVec = radialVectorWithLength(mhdInput, geom, size(amplitude, 1));
    xLabelText = '$\sqrt{s}$';
end

function validatePoloidalModeOptions(opt)

    mRange = opt.mRange;
    isValidMRange = isnumeric(mRange) && isreal(mRange) && numel(mRange) == 2 && ...
        all(isfinite(mRange(:))) && all(mRange(:) == floor(mRange(:))) && all(mRange(:) >= 0);
    assert(isValidMRange, 'mRange 必须包含两个有限非负整数。');
    requireFiniteScalarInRange(opt.minFFT, 0, 1, 'minFFT');
    requireLogicalScalar(getOptionValue(opt, 'normalize', true), 'normalize');
end

function [pestRefined, alignedField, refinedYGrid] = prepareMHDPESTField(fieldZXY, geom, opt, targetThetaCount)

    alignedField = undoMHDFieldAlignedShift(fieldZXY, geom);
    alignedGhost = buildFieldAlignedGhost(alignedField, geom);
    [alignedRefined, refinedYGrid, refinedIndex] = refineFieldAlignedY(alignedGhost, geom, opt, targetThetaCount);
    pestRefined = shiftAlignedToRefinedPest(alignedRefined, refinedIndex, geom);
end

function alignedGhost = buildFieldAlignedGhost(alignedField, geom)

    [nZ, nX, nY] = size(alignedField);
    gridGhost = geom.gridGhost;
    assert(nY > gridGhost, 'gridNy 必须大于 gridGhost。');

    alignedGhost = zeros(nZ, nX, nY + 2 * gridGhost);
    alignedGhost(:, :, gridGhost + 1:gridGhost + nY) = alignedField;

    deltaQtheta = (geom.qtheta(:, 2) - geom.qtheta(:, 1)) * nY;
    for iX = 1:nX
        for iGhost = 1:gridGhost
            leftSource = nY - gridGhost + iGhost;
            rightSource = iGhost;
            alignedGhost(:, iX, iGhost) = shiftMHDFieldZ(alignedField(:, iX, leftSource), ...
                -deltaQtheta(iX), geom);
            alignedGhost(:, iX, gridGhost + nY + iGhost) = shiftMHDFieldZ(alignedField(:, iX, rightSource), ...
                deltaQtheta(iX), geom);
        end
    end
end

function [alignedRefined, refinedYGrid, refinedIndex] = refineFieldAlignedY(alignedGhost, geom, opt, targetThetaCount)

    [nZ, nX, ~] = size(alignedGhost);
    nY = numel(geom.yGrid);
    gridGhost = geom.gridGhost;
    refinedTimes = poloidalFFTRefinedTimes(opt, nY, targetThetaCount);
    refinedNy = nY * refinedTimes;
    dtheta = 2.0 * pi / nY;
    theta0 = geom.yGrid(1);

    refinedIndex = ((0:refinedNy - 1) + 0.5) / refinedTimes - 0.5;
    refinedYGrid = theta0 + refinedIndex * dtheta;
    alignedRefined = zeros(nZ, nX, refinedNy);
    for iY = 1:refinedNy
        shiftJ = floor(refinedIndex(iY));
        shiftDj = refinedIndex(iY) - shiftJ;
        alignedRefined(:, :, iY) = lagrangeFourPoint( ...
            alignedGhost(:, :, shiftJ + gridGhost), ...
            alignedGhost(:, :, shiftJ + gridGhost + 1), ...
            alignedGhost(:, :, shiftJ + gridGhost + 2), ...
            alignedGhost(:, :, shiftJ + gridGhost + 3), shiftDj);
    end
end

function refinedTimes = poloidalFFTRefinedTimes(opt, gridNy, targetThetaCount)

    refinedTimes = 2;
    if isfield(opt, 'mRange') && numel(opt.mRange) >= 2
        maxM = max(0, max(round(double(opt.mRange(1:2)))));
    else
        maxM = floor(gridNy / 2);
    end

    minRefinedNy = 16 * maxM;
    if ~isempty(targetThetaCount)
        minRefinedNy = max(minRefinedNy, double(targetThetaCount));
    end

    while gridNy * refinedTimes < minRefinedNy
        refinedTimes = refinedTimes * 2;
    end
end

function pestRefined = shiftAlignedToRefinedPest(alignedRefined, refinedIndex, geom)

    [nZ, nX, refinedNy] = size(alignedRefined);

    qtheta0 = geom.qtheta(:, 1);
    dqtheta = geom.qtheta(:, 2) - geom.qtheta(:, 1);
    pestRefined = zeros(nZ, nX, refinedNy);
    for iX = 1:nX
        qthetaRefined = qtheta0(iX) + refinedIndex * dqtheta(iX);
        for iY = 1:refinedNy
            pestRefined(:, iX, iY) = shiftMHDFieldZ(alignedRefined(:, iX, iY), ...
                -qthetaRefined(iY), geom);
        end
    end
end

function fieldShifted = shiftMHDFieldZ(fieldZ, shiftZ, geom)

    nZ = numel(geom.zGrid);
    gridDz = 2.0 * pi / geom.tubes / nZ;
    shiftLk = shiftZ / gridDz;
    shiftK = floor(shiftLk);
    shiftDk = shiftLk - shiftK;

    fieldZ = fieldZ(:);
    kIndex = (0:nZ - 1).';
    fieldShifted = lagrangeFourPoint( ...
        fieldZ(mod(kIndex + shiftK - 1, nZ) + 1), ...
        fieldZ(mod(kIndex + shiftK + 0, nZ) + 1), ...
        fieldZ(mod(kIndex + shiftK + 1, nZ) + 1), ...
        fieldZ(mod(kIndex + shiftK + 2, nZ) + 1), shiftDk);
end

function field = lagrangeFourPoint(field0, field1, field2, field3, shiftD)

    field = -shiftD .* (-1 + shiftD) .* (-2 + shiftD) / 6 .* field0 + ...
        (1 + shiftD) .* (-1 + shiftD) .* (-2 + shiftD) / 2 .* field1 - ...
        shiftD .* (1 + shiftD) .* (-2 + shiftD) / 2 .* field2 + ...
        shiftD .* (1 + shiftD) .* (-1 + shiftD) / 6 .* field3;
end

function candidateM = poloidalModeCandidates(opt, nTheta)

    maxM = max(0, floor(nTheta / 2));
    mRange = sort(double(opt.mRange(:)'));
    mRange(1) = max(0, mRange(1));
    mRange(2) = min(maxM, mRange(2));

    if mRange(2) < mRange(1)
        candidateM = [];
    else
        candidateM = mRange(1):mRange(2);
    end
end

function selectedM = selectPoloidalModeNumbers(amplitude, candidateM, opt)

    if isempty(candidateM)
        selectedM = [];
        return;
    end

    candidateAmplitude = amplitude(:, candidateM + 1);
    radialMax = max(candidateAmplitude, [], 1);
    radialMax(~isfinite(radialMax)) = 0;
    globalMax = max(radialMax);
    if globalMax <= 0
        selectedM = candidateM;
        return;
    end

    selectedMask = radialMax > globalMax * opt.minFFT;
    if ~any(selectedMask)
        [~, maxIndex] = max(radialMax);
        selectedMask(maxIndex) = true;
    end

    selectedM = candidateM(selectedMask);
end

function normalizeFFT = useNormalizedPoloidalMode(opt)

    normalizeFFT = logical(getOptionValue(opt, 'normalize', true));
end

function xVec = radialVectorWithLength(mhdInput, geom, nRadial)

    if isfield(geom, 'rhoGrid') && ~isempty(geom.rhoGrid) && numel(geom.rhoGrid) == nRadial
        radialValues = geom.rhoGrid;
    elseif isfield(geom, 'rhoplot') && ~isempty(geom.rhoplot)
        radialValues = geom.rhoplot;
    elseif isfield(mhdInput, 'rhoplot') && ~isempty(mhdInput.rhoplot)
        radialValues = mhdInput.rhoplot;
    elseif isfield(geom, 'rhoGrid') && ~isempty(geom.rhoGrid)
        radialValues = geom.rhoGrid;
    else
        radialValues = linspace(0, 1, nRadial);
    end

    if ~isvector(radialValues)
        radialValues = radialValues(:, 1);
    end
    radialValues = radialValues(:);

    if numel(radialValues) == nRadial
        xVec = radialValues;
    else
        xVec = linspace(0, 1, nRadial).';
    end
end

function ctx = buildMHDDiagnosticContext(meta, mhdInput)

    ctx = struct();
    ctx.xGrid = meta.xGrid;
    ctx.tDiag = meta.tDiag;
    ctx.diagSteps = meta.diagSteps;
    ctx.modeIndexAll = meta.modeIndexAll;
    ctx.physicalNAll = meta.physicalNAll;
    ctx.B0 = mhdInput.B0;
    ctx.L0 = mhdInput.L0;
    ctx.VA0 = mhdInput.VA0;
    ctx.TeSample = mhdInput.TeSample;
    ctx.rho = mhdInput.rho;
end

function validateDiagnosticGrid(ctx, minTime, contextText, data, fieldNames)

    dataSize = size(data.(fieldNames{1}));
    for iField = 2:numel(fieldNames)
        assert(isequal(size(data.(fieldNames{iField})), dataSize), ...
            '%s 诊断数组尺寸必须一致。', contextText);
    end
    assert(dataSize(1) >= minTime, '%s 至少需要 %d 个诊断时间点。', contextText, minTime);
    assert(numel(ctx.tDiag) == dataSize(1), '%s 时间维度必须匹配 tDiag。', contextText);
    assert(numel(ctx.xGrid) == dataSize(2), '%s 径向维度必须匹配 xGrid。', contextText);
end

function workspace = runMHDTimeDiagnostic(mhdData, ctx, opt, fieldNames, plotName, ...
    missingReason, minTime, staticPlotFcn, interactivePlotFcn)

    [workspace, ready] = beginMHDDiagnostic(opt, mhdData, fieldNames, plotName, missingReason);
    if ~ready
        return;
    end
    validateDiagnosticGrid(ctx, minTime, plotName, mhdData, fieldNames);
    runInteractivePlot(getOptionValue(opt, 'interactive', 0), staticPlotFcn, interactivePlotFcn);
end

function workspace = runMHDIndexedRadialDiagnostic( ...
    mhdData, ctx, opt, fieldNames, plotName, missingReason, figName, buildPlotDataFcn)

    [workspace, ready] = beginMHDDiagnostic(opt, mhdData, fieldNames, plotName, missingReason);
    if ~ready
        return;
    end

    dataSize = size(mhdData.(fieldNames{1}));
    for iField = 2:numel(fieldNames)
        assert(isequal(size(mhdData.(fieldNames{iField})), dataSize), ...
            '%s 诊断数组尺寸必须一致。', plotName);
    end
    assert(size(ctx.rho, 1) == dataSize(2), '%s 径向维度必须匹配 rho(:, 1)。', plotName);

    runInteractivePlot(opt.interactive, ...
        @() drawStaticDiagnostic(buildPlotDataFcn(opt), opt, @renderLinePlot), ...
        @(dynamicUpdate) plotIndexOptionInteractive( ...
        figName, 'timeIndex', dataSize(1), dynamicUpdate, opt, buildPlotDataFcn));
end

function workspace = runMHDAmplitudePlot(mhdData, ctx, opt)

    workspace = runMHDTimeDiagnostic( ...
        mhdData, ctx, opt, {'amplitude'}, 'amplitude plot', '未读取 amplitude 数据', 1, ...
        @() drawStaticDiagnostic(buildAmplitudePlotData(mhdData.amplitude, ctx, opt, false), opt, @renderLinePlot), ...
        @(dynamicUpdate) plotAmplitudeInteractive(mhdData.amplitude, ctx, opt, dynamicUpdate));
end

function workspace = runMHDMultipleFrequencyPlot(mhdData, ctx, opt)

    workspace = runMHDTimeDiagnostic( ...
        mhdData, ctx, opt, {'RealMode', 'ImagMode'}, 'mode frequency plot', ...
        '未读取 RealMode 或 ImagMode 数据', 2, ...
        @() drawStaticDiagnostic( ...
        buildMultipleFrequencyPlotData(mhdData.RealMode, mhdData.ImagMode, ctx, opt), opt, @renderLinePlot), ...
        @(dynamicUpdate) plotMultipleFrequencyInteractive( ...
        mhdData.RealMode, mhdData.ImagMode, ctx, opt, dynamicUpdate));
end

function workspace = runMHDPhaseFrequencyPlot(mhdData, ctx, opt)

    workspace = runMHDTimeDiagnostic( ...
        mhdData, ctx, opt, {'RealMode', 'ImagMode'}, 'phase frequency plot', ...
        '未读取 RealMode 或 ImagMode 数据', 2, ...
        @() drawStaticDiagnostic( ...
        buildPhaseFrequencyPlotData(mhdData.RealMode, mhdData.ImagMode, ctx, opt), opt, @renderLinePlot), ...
        @(dynamicUpdate) plotPhaseFrequencyInteractive( ...
        mhdData.RealMode, mhdData.ImagMode, ctx, opt, dynamicUpdate));
end

function workspace = runMHDContourFrequencyPlot(mhdData, ctx, opt)

    workspace = runMHDTimeDiagnostic( ...
        mhdData, ctx, opt, {'amplitude', 'RealMode'}, 'contour frequency plot', ...
        '未读取 amplitude 或 RealMode 数据', 2, ...
        @() drawStaticDiagnostic(buildContourFrequencyPlotData( ...
        mhdData.amplitude, mhdData.RealMode, ctx, opt), opt, @renderMapPlot), []);
end

function workspace = runMHDSingleSignalPlot(mhdData, ctx, opt)

    workspace = runMHDTimeDiagnostic( ...
        mhdData, ctx, opt, {'frequency'}, 'single signal plot', '未读取 frequency.bin 数据', ...
        1, ...
        @() drawStaticDiagnostic(buildSingleSignalPlotData(mhdData.frequency, ctx, opt), opt, @renderLinePlot), ...
        @(dynamicUpdate) plotSingleSignalInteractive(mhdData.frequency, ctx, opt, dynamicUpdate));
end

function workspace = runMHDEparaPlot(mhdData, ctx, opt)

    workspace = runMHDIndexedRadialDiagnostic( ...
        mhdData, ctx, opt, {'Epara', 'EparaES'}, 'Epara plot', '未读取 Epara 或 EparaES 数据', ...
        'Epara / EparaES', @(tempOpt) buildEparaPlotData(mhdData.Epara, mhdData.EparaES, ctx, tempOpt));
end

function workspace = runMHDShearingPlot(mhdData, ctx, opt)

    [workspace, ready] = beginMHDDiagnostic( ...
        opt, mhdData, {'Shearing'}, 'shearing plot', '未读取 Shearing 数据');
    if ~ready
        return;
    end

    validateDiagnosticGrid(ctx, 1, 'shearing', mhdData, {'Shearing'});

    plotType = requireIntegerInRange(opt.plotType, 1, 3, 'plotType');
    switch plotType
        case 1
            runInteractivePlot(opt.interactive, ...
                @() drawStaticDiagnostic( ...
                buildShearingTimeSlicePlotData(mhdData.Shearing, ctx, opt), opt, @renderLinePlot), ...
                @(dynamicUpdate) plotIndexOptionInteractive( ...
                'Shearing', 'timeIndex', size(mhdData.Shearing, 1), dynamicUpdate, opt, ...
                @(tempOpt) buildShearingTimeSlicePlotData(mhdData.Shearing, ctx, tempOpt)));
        case 2
            runInteractivePlot(opt.interactive, ...
                @() drawStaticDiagnostic( ...
                buildShearingRadialTracePlotData(mhdData.Shearing, ctx, opt), opt, @renderLinePlot), ...
                @(dynamicUpdate) plotIndexOptionInteractive( ...
                'Shearing trace', 'radialIndex', size(mhdData.Shearing, 2), dynamicUpdate, opt, ...
                @(tempOpt) buildShearingRadialTracePlotData(mhdData.Shearing, ctx, tempOpt)));
        case 3
            runInteractivePlot(opt.interactive, ...
                @() drawStaticDiagnostic(buildShearingMapPlotData(mhdData.Shearing, ctx), opt, @renderMapPlot), []);
    end
end

function workspace = runMHDZonalDrivePlot(mhdData, ctx, opt)

    workspace = runMHDIndexedRadialDiagnostic( ...
        mhdData, ctx, opt, {'MaxwellDrive', 'ReynoldsDrive', 'ZonalDrive'}, ...
        'ZF drive plot', '未读取 MaxwellDrive/ReynoldsDrive/ZonalDrive 数据', ...
        'MaxwellDrive / ReynoldsDrive / ZonalDrive', ...
        @(tempOpt) buildZonalDrivePlotData(mhdData.MaxwellDrive, mhdData.ReynoldsDrive, ...
        mhdData.ZonalDrive, ctx, tempOpt));
end

function [modeIndex, physicalN, modePosition] = parseModeN(modeN, modeIndexAll, physicalNAll, contextText)

    assert(~isempty(modeIndexAll) && numel(modeIndexAll) == numel(physicalNAll), ...
        '%s modeN 列表为空或尺寸不一致。', contextText);
    if isempty(modeN)
        modeIndex = modeIndexAll;
        physicalN = physicalNAll;
        modePosition = 1:numel(modeIndexAll);
        return;
    end

    if ischar(modeN) || isstring(modeN)
        key = lower(strtrim(char(modeN)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            idx = max(1, round(numel(physicalNAll) / 2));
            modeIndex = modeIndexAll(idx);
            physicalN = physicalNAll(idx);
            modePosition = idx;
            return;
        elseif ismember(key, {'end', 'last'})
            modeIndex = modeIndexAll(end);
            physicalN = physicalNAll(end);
            modePosition = numel(modeIndexAll);
            return;
        elseif ismember(key, {'first', 'begin'})
            modeIndex = modeIndexAll(1);
            physicalN = physicalNAll(1);
            modePosition = 1;
            return;
        else
            physicalN = str2double(key);
        end
    else
        physicalN = double(modeN(:)');
    end

    assert(all(isfinite(physicalN)) && all(physicalN == floor(physicalN)) && all(physicalN >= 0), ...
        '%s modeN 必须为非负整数物理环向模数。', contextText);
    [isMember, memberIdx] = ismember(physicalN, physicalNAll);
    assert(all(isMember), '%s modeN 必须位于有效物理 n 集合 [%s]。当前值：[%s]。', ...
        contextText, formatNumberList(physicalNAll), formatNumberList(physicalN));
    modeIndex = modeIndexAll(memberIdx);
    modePosition = memberIdx;
end

function idx = parseIndex(indexText, nIndex, label)

    if isnumeric(indexText)
        idx = double(indexText);
    else
        key = lower(strtrim(char(indexText)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            idx = max(1, round(nIndex / 2));
        elseif ismember(key, {'end', 'last'})
            idx = nIndex;
        else
            idx = str2double(key);
        end
    end

    idx = requireIndex(idx, nIndex, label);
end

function idx = parseTimeIndex(timeIndex, nTime)

    idx = parseIndex(timeIndex, nTime, 'timeIndex');
end

function idx = requireIndex(idx, maxIndex, label)

    idx = requireIntegerInRange(idx, 1, maxIndex, [label ' 下标']);
end

function logLoaded(name, data)
    fprintf('[load] %s: size=[%s]\n', char(name), formatSize(size(data)));
end

function logSkipped(name, reason)
    fprintf('[skip] %s: %s\n', char(name), reason);
end

function ok = hasBSpline()
    fileStatus = exist('bspline', 'file');
    builtinStatus = exist('bspline', 'builtin');
    ok = any(fileStatus == [2, 3, 6]) || builtinStatus == 5;
end

function ok = hasMHDDataFields(dataStruct, fieldNames)

    ok = all(cellfun(@(fieldName) ...
        isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName)), fieldNames));
end

function text = formatSize(dataSize)
    text = strjoin(cellstr(compose('%d', dataSize(:))), ' ');
end

function text = formatNumberList(values)
    if isempty(values)
        text = '';
    else
        text = strjoin(cellstr(compose('%g', values(:))), ', ');
    end
end

function value = readIntParam(paramText, name)
    token = regexp(paramText, ['const\s+int\s+' name '\s*=\s*([0-9]+)\s*;'], 'tokens', 'once');
    assert(~isempty(token), '找不到非负整数参数：%s', name);
    value = str2double(token{1});
end

function value = readFloatParam(paramText, name)
    pattern = ['const\s+(?:double|float|mhdReal|picReal)\s+' name ...
        '\s*=\s*([-+]?[0-9eE+\-\.]+)\s*;'];
    token = regexp(paramText, pattern, ...
        'tokens', 'once');
    assert(~isempty(token), '找不到浮点参数：%s', name);
    value = str2double(token{1});
end

function value = readSwitchParam(paramText, name)
    token = regexp(paramText, ['using\s+' name '\s*=\s*(trueType|falseType)\s*;'], 'tokens', 'once');
    assert(~isempty(token), '找不到开关参数：%s', name);
    value = strcmp(token{1}, 'trueType');
end

function precision = readPrecisionParam(paramText)
    token = regexp(paramText, 'using\s+mhdReal\s*=\s*(double|float)\s*;', 'tokens', 'once');
    assert(~isempty(token), '找不到 mhdReal 精度定义。');
    precision = token{1};
end

function raw = readBinaryVector(filePath, precision)
    assert(isfile(filePath), '缺少文件：%s', filePath);

    fid = fopen(filePath, 'rb');
    assert(fid >= 0, '无法打开文件：%s', filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    raw = fread(fid, inf, freadPrecision(precision));
end

function precisionText = freadPrecision(precision)
    switch lower(strtrim(char(precision)))
        case 'double'
            precisionText = 'double=>double';
        case {'float', 'single'}
            precisionText = 'single=>double';
        otherwise
            error('不支持的精度类型：%s。', precision);
    end
end

% 以下读取函数返回的数组维序为：TXN=[time,x,mode]，TX=[time,x]，
% TYXZ=[time,y,x,z]，ZXY=[z,x,y]，YXZ=[y,x,z]。
function data = readModeDiagnosticAsTXN(filePath, precision, nTime, gridNx, nMode)
    dimensionText = sprintf('nTime=%d, gridNx=%d, nMode=%d', nTime, gridNx, nMode);
    data = readBinaryArray(filePath, precision, [nMode, gridNx, nTime], [3, 2, 1], dimensionText);
end

function data = readRadialDiagnosticAsTX(filePath, precision, nTime, gridNx)
    dimensionText = sprintf('nTime=%d, gridNx=%d', nTime, gridNx);
    data = readBinaryArray(filePath, precision, [gridNx, nTime], [2, 1], dimensionText);
end

function data = readOutputAsTYXZ(filePath, precision, nTime, gridNy, gridNx, gridNz)
    dimensionText = sprintf('nTime=%d, gridNy=%d, gridNx=%d, gridNz=%d', ...
        nTime, gridNy, gridNx, gridNz);
    data = readBinaryArray(filePath, precision, ...
        [gridNz, gridNx, gridNy, nTime], [4, 3, 2, 1], dimensionText);
end

function fieldZXY = readMHDFieldAsZXY(inputDir, fieldName, precision, gridNy, gridNx, gridNz)
    fieldNameText = char(fieldName);
    filePath = fullfile(inputDir, [fieldNameText '.bin']);
    dimensionText = sprintf('gridNy=%d, gridNx=%d, gridNz=%d', gridNy, gridNx, gridNz);
    fieldZXY = readBinaryArray(filePath, precision, ...
        [gridNz, gridNx, gridNy], [1, 2, 3], dimensionText);
end

function data = readBinaryArray(filePath, precision, storageShape, permutation, dimensionText)
    raw = readBinaryVector(filePath, precision);
    expectedCount = prod(storageShape);
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（%s）。', ...
        filePath, numel(raw), expectedCount, dimensionText);

    data = permute(reshape(raw, storageShape), permutation);
end

function fieldZXY = totalMHDFieldTimeSliceAsZXY(totalDataTYXZ, timeIndex)
    fieldYXZ = reshape(totalDataTYXZ(timeIndex, :, :, :), ...
        [size(totalDataTYXZ, 2), size(totalDataTYXZ, 3), size(totalDataTYXZ, 4)]);
    fieldZXY = permute(fieldYXZ, [3, 2, 1]);
end

function fieldYXZ = mhdFieldZXYToYXZ(fieldZXY)
    fieldYXZ = permute(fieldZXY, [3, 2, 1]);
end

function fieldFilteredZXY = filterMHDToroidalModes(fieldZXY, modeIndexKeep)
    modeIndexKeep = unique(modeIndexKeep(:)');
    nZ = size(fieldZXY, 1);
    assert(all(abs(modeIndexKeep) <= floor(nZ / 2)), 'modeIndexKeep 超出 z 向 FFT 可解析范围。');

    modeMask = false(nZ, 1);
    for iMode = 1:numel(modeIndexKeep)
        modeAbs = abs(modeIndexKeep(iMode));
        modeMask(mod(modeAbs, nZ) + 1) = true;
        if modeAbs ~= 0
            modeMask(mod(-modeAbs, nZ) + 1) = true;
        end
    end

    fieldSpectrum = fft(fieldZXY, [], 1);
    fieldSpectrum(~modeMask, :, :) = 0;
    fieldFilteredZXY = real(ifft(fieldSpectrum, [], 1));
end

function geom = buildMHDFieldPlotGeometry( ...
    q, theta_pest, rho, qplot, rhoplot, Rplot, Zplot, yGrid, zGrid, tubes, gridGhost, NFP, phi)
    [nPlotRho, nPlotTheta, nPlotPhi] = size(qplot);
    thetaPlotGrid = (0.5:nPlotTheta - 0.5) / nPlotTheta * 2 * pi - pi;
    q2D = mhdPlotToroidalSlice(q, 1);
    theta2D = mhdPlotToroidalSlice(theta_pest, 1);
    rho2D = mhdPlotToroidalSlice(rho, 1);
    has3DPlot = (ndims(qplot) >= 3 && nPlotPhi > 1);
    qplot2D = mhdPlotToroidalSlice(qplot, 1);
    thetaPlot2D = repmat(thetaPlotGrid, nPlotRho, 1);
    geom = struct();
    geom.qtheta = q2D .* theta2D;
    geom.rhoGrid = rho2D(:, 1);
    geom.qplotAll = qplot;
    geom.rhoplotAll = rhoplot;
    geom.RplotAll = Rplot;
    geom.ZplotAll = Zplot;
    geom.qplot = qplot2D;
    geom.Rplot = mhdPlotToroidalSlice(Rplot, 1);
    geom.Zplot = mhdPlotToroidalSlice(Zplot, 1);
    geom.rhoplot = mhdPlotToroidalSlice(rhoplot, 1);
    geom.thplot = thetaPlot2D;
    geom.qthetaplot = qplot2D .* thetaPlot2D;
    geom.nPlotTheta = nPlotTheta;
    geom.nPlotPhi = nPlotPhi;
    geom.plotPhiGrid = mhdPlotPhiGrid(phi, nPlotPhi, zGrid);
    geom.has3DPlot = has3DPlot;
    geom.NFP = NFP;
    geom.isTokamak = (NFP == 1);
    geom.isStellarator = (NFP ~= 1);
    geom.yGrid = yGrid;
    geom.zGrid = zGrid;
    geom.tubes = tubes;
    geom.gridGhost = gridGhost;
end

function phiGrid = mhdPlotPhiGrid(phi, nPlotPhi, zGrid)

    if isempty(phi)
        if nPlotPhi == numel(zGrid)
            phiGrid = zGrid(:);
        else
            phiGrid = zeros(nPlotPhi, 1);
        end
        return;
    end

    if ndims(phi) >= 3
        phiGrid = squeeze(phi(1, 1, :));
    elseif isvector(phi)
        phiGrid = phi(:);
    else
        phiGrid = squeeze(phi(1, :));
    end

    assert(numel(phiGrid) == nPlotPhi, ...
        'phi 网格长度必须匹配 plot3D.mat 的环向维度。');
    phiGrid = double(phiGrid(:));
end

function [planeGeom, phi0, toroidalIndex] = mhdPlotPlaneGeometry(geom, opt)

    if geom.isStellarator
        assert(geom.has3DPlot, '仿星器绘图需要 plot3D.mat。');
        requestedIndex = getOptionValue(opt, 'toroidalIndex', 1);
        toroidalIndex = parseToroidalIndex(requestedIndex, geom.nPlotPhi);
        phi0 = geom.plotPhiGrid(toroidalIndex);
    else
        toroidalIndex = 1;
        phi0 = requireFiniteScalarInRange( ...
            getOptionValue(opt, 'toroidalAngle', 0.0), -Inf, Inf, 'toroidalAngle');
    end

    planeGeom = geom;
    planeGeom.qplot = mhdPlotToroidalSlice(geom.qplotAll, toroidalIndex);
    planeGeom.rhoplot = mhdPlotToroidalSlice(geom.rhoplotAll, toroidalIndex);
    planeGeom.Rplot = mhdPlotToroidalSlice(geom.RplotAll, toroidalIndex);
    planeGeom.Zplot = mhdPlotToroidalSlice(geom.ZplotAll, toroidalIndex);
    planeGeom.qthetaplot = planeGeom.qplot .* geom.thplot;
end

function value2D = mhdPlotToroidalSlice(value, toroidalIndex)

    if ndims(value) >= 3 && size(value, 3) > 1
        value2D = squeeze(value(:, :, toroidalIndex));
    else
        value2D = squeeze(value);
    end
end

function toroidalIndex = parseToroidalIndex(toroidalIndex, nPlotPhi)

    toroidalIndex = double(toroidalIndex);
    toroidalIndex = requireIndex(toroidalIndex, nPlotPhi, 'toroidalIndex');
end

function text = mhdToroidalLocationText(geom, phi0, toroidalIndex)

    if geom.isStellarator
        text = sprintf('phiIndex=%d, phi=%.6g', toroidalIndex, phi0);
    else
        text = sprintf('phi=%.6g', phi0);
    end
end

function [shifted, aligned] = plotMHDFieldOnPoloidalPlane(fieldZXY, fieldName, geom, opt, physicalN, contextText)
    [pestRefined, fieldNonShiftZXY, refinedYGrid] = prepareMHDPESTField( ...
        fieldZXY, geom, opt, geom.nPlotTheta);
    shifted = mhdFieldZXYToYXZ(fieldZXY);
    aligned = mhdFieldZXYToYXZ(fieldNonShiftZXY);
    fieldPlotZXY = normalizeMHDFieldForPlot(pestRefined);
    [result, planeGeom, phi0, toroidalIndex] = interpolateMHDFieldToPlotGrid( ...
        fieldPlotZXY, refinedYGrid, geom, opt);
    drawMHDFieldSurface(result, fieldName, planeGeom, opt);

    if strlength(string(contextText)) > 0
        contextText = [char(contextText), ', '];
    end
    fprintf('[plot] %s: %s%s, n=[%s]\n', ...
        char(fieldName), char(contextText), ...
        mhdToroidalLocationText(geom, phi0, toroidalIndex), formatNumberList(physicalN));
end

function fieldPlotZXY = normalizeMHDFieldForPlot(fieldZXY)
    scale = max(abs(fieldZXY(:)));
    fieldPlotZXY = fieldZXY;
    if isfinite(scale) && scale > 0
        fieldPlotZXY = fieldZXY / scale;
    end
end

function fieldNonShiftZXY = undoMHDFieldAlignedShift(fieldZXY, geom)
    [nZ, nX, nY] = size(fieldZXY);
    assert(numel(geom.zGrid) == nZ, 'zGrid 长度必须匹配场量 z 维度。');
    assert(isequal(size(geom.qtheta), [nX, nY]), 'qtheta 尺寸必须匹配场量 x/y 维度。');

    fieldNonShiftZXY = zeros(size(fieldZXY));

    for iX = 1:nX
        for iY = 1:nY
            fieldNonShiftZXY(:, iX, iY) = shiftMHDFieldZ(fieldZXY(:, iX, iY), ...
                geom.qtheta(iX, iY), geom);
        end
    end
end

function [result, planeGeom, phi0, toroidalIndex] = interpolateMHDFieldToPlotGrid( ...
    fieldPestZXY, refinedYGrid, geom, opt)
    [planeGeom, phi0, toroidalIndex] = mhdPlotPlaneGeometry(geom, opt);
    plotZ = phi0 + zeros(size(planeGeom.rhoplot));
    coordinate3D = [plotZ(:), planeGeom.rhoplot(:), planeGeom.thplot(:)];
    derivative = uint64([0, 0, 0])';
    isPeriodic = [true, false, false]';
    rangeZRhoTheta = { ...
        [geom.zGrid(1), geom.zGrid(1) + 2 * pi / geom.tubes]; ...
        [geom.rhoGrid(1), geom.rhoGrid(end)]; ...
        [refinedYGrid(1), refinedYGrid(1) + 2 * pi]};

    result = bspline(uint64(4), isPeriodic, rangeZRhoTheta, ...
        fieldPestZXY, coordinate3D, derivative);
    result = reshape(result, size(planeGeom.Rplot));
end

function drawMHDFieldSurface(result, fieldName, geom, opt)
    figure;
    axHandle = gca;
    pcolor(axHandle, geom.Rplot, geom.Zplot, result);
    shading(axHandle, 'interp');
    axis(axHandle, 'equal');
    hold(axHandle, 'on');

    finiteResult = result(isfinite(result));
    if isempty(finiteResult)
        colorMin = -1;
        colorMax = 1;
    else
        colorMin = min(finiteResult);
        colorMax = max(finiteResult);
        if colorMin == colorMax
            colorMin = colorMin - 1;
            colorMax = colorMax + 1;
        end
    end

    plot(axHandle, geom.Rplot(1, :), geom.Zplot(1, :), 'b', 'LineWidth', 1.5);
    plot(axHandle, geom.Rplot(end, :), geom.Zplot(end, :), 'b', 'LineWidth', 1.5);

    xlim(axHandle, [min(geom.Rplot(:)) - 0.5, max(geom.Rplot(:)) + 0.5]);
    ylim(axHandle, [min(geom.Zplot(:)) - 0.5, max(geom.Zplot(:)) + 0.5]);

    clim(axHandle, [colorMin, colorMax]);
    colormap(axHandle, mhdFieldColormap(opt.colormapIndex, colorMin, colorMax, 256));

    applyPlotAxesStyle(axHandle);
    set(axHandle, 'LineWidth', 1.2);

    xlabel(axHandle, '$R/\mathrm{m}$', 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(axHandle, '$Z/\mathrm{m}$', 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', 14);
    title(axHandle, char(fieldName), 'Interpreter', 'none', ...
        'FontName', 'Times New Roman', 'FontSize', 14);
    box(axHandle, 'on');
end

function cmap = mhdFieldColormap(colormapIndex, valueMin, valueMax, nColor)
    colormapIndex = requireIntegerInRange(colormapIndex, 1, 5, 'colormapIndex');
    switch colormapIndex
        case 1
            x = [valueMin, 2 / 3 * valueMin, 1 / 3 * valueMin, 0, 1 / 3 * valueMax, 2 / 3 * valueMax, valueMax];
            colors = [0, 1, 1; 0.3, 0.3, 0.8; 0, 0, 1; 1, 1, 1; 1, 0.2, 0.2; 1, 0.3, 0.6; 1, 1, 0];
        case 2
            x = [valueMin, 1 / 2 * valueMin, 0, 1 / 2 * valueMax, valueMax];
            colors = [0.3, 0.3, 0.8; 0, 0, 1; 1, 1, 1; 1, 0.2, 0.2; 1, 0.3, 0.6];
        case 3
            x = [valueMin, 1 / 2 * valueMin, 0, 1 / 2 * valueMax, valueMax];
            colors = [0, 1, 1; 0, 0, 1; 1, 1, 1; 1, 0.2, 0.2; 1, 1, 0];
        case 4
            x = [valueMin, 0, valueMax];
            colors = [0, 0, 1; 1, 1, 1; 1, 0, 0];
        case 5
            cmap = jet(nColor);
            return;
    end

    xq = linspace(valueMin, valueMax, nColor);
    if valueMin < 0 && valueMax > 0
        cmap = interp1(x, colors, xq, 'linear');
    else
        cmap = interp1(linspace(valueMin, valueMax, size(colors, 1)), colors, xq, 'linear');
    end
end

function [xData, xLabelText] = diagnosticTimeAxis(tAValues, stepValues, L0, VA0, timeAxis, errorText)

    switch char(lower(timeAxis))
        case {'ta', 'alfven', 'alfven_time'}
            xData = tAValues;
            xLabelText = '$t_a$';
        case {'ms', 'millisecond', 'milliseconds'}
            xData = tAValues * L0 / VA0 * 1000;
            xLabelText = '$t/\mathrm{ms}$';
        case {'steps', 'step'}
            xData = stepValues;
            xLabelText = '$\mathrm{steps}$';
        otherwise
            error(errorText, timeAxis);
    end
end

function plotAmplitudeInteractive(amplitude, ctx, opt, dynamicUpdate)
    nRadial = size(amplitude, 2);
    radialIndex0 = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    opt.radialIndex = radialIndex0;
    [xData0, ~, ~, ~, ~] = ...
        calculateAmplitudeData(amplitude, ctx, opt);
    [growthStart0, growthEnd0] = amplitudeGrowthInitialRange(xData0, opt);
    controls = [ ...
        buildIntegerSliderControl('radialIndex', 'radialIndex', radialIndex0, 1, nRadial), ...
        buildNumericSliderControl( ...
            'growthStart', 'growthStart', growthStart0, xData0(1), xData0(end), numel(xData0)), ...
        buildNumericSliderControl('growthEnd', 'growthEnd', growthEnd0, xData0(1), xData0(end), numel(xData0))];
    plotInteractiveLineWithSliders('Amplitude', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.radialIndex = values.radialIndex;
        tempOpt.growthRange = sort([values.growthStart, values.growthEnd]);
        data = buildAmplitudePlotData(amplitude, ctx, tempOpt, true);
    end
end

function plotData = buildAmplitudePlotData(amplitude, ctx, opt, includeGrowth)
    [xData, logAmplitude, xLabelText, toroidalModeN, radialX, lnAmplitude] = ...
        calculateAmplitudeData(amplitude, ctx, opt);
    plotData = buildLinePlotData(xData, logAmplitude, xLabelText, '', '$e\delta\phi/T_e$', ...
        compose('$n=%d$', toroidalModeN), '');
    plotData.logText = sprintf('amplitude: x=%.6g, radialIndex=%d, n=[%s]', ...
        radialX, opt.radialIndex, formatNumberList(toroidalModeN));
    if ~includeGrowth
        return;
    end

    [growthRate, growthUnit, growthRange] = calculateAmplitudeGrowthRate(lnAmplitude, xData, ctx, opt);
    rangeMask = xData(:) >= growthRange(1) & xData(:) <= growthRange(2);
    maxItems = compose('n=%d: %.6g', toroidalModeN(:), max(logAmplitude(rangeMask, :), [], 1).');
    plotData.status = amplitudeGrowthStatus( ...
        radialX, opt.radialIndex, growthRange, growthUnit, toroidalModeN, growthRate);
    plotData.status = sprintf('%s, max log10(ephi/Te): %s', ...
        plotData.status, strjoin(cellstr(maxItems), ', '));
    plotData.xLines = growthRange;
end

function [xData, logAmplitude, xLabelText, toroidalModeN, radialX, lnAmplitude] = ...
    calculateAmplitudeData(amplitude, ctx, opt)
    nTime = size(amplitude, 1);

    [~, toroidalModeN, modeIdx] = parseModeN( ...
        opt.modeN, ctx.modeIndexAll, ctx.physicalNAll, 'Amplitude');

    nRadial = size(amplitude, 2);
    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = ctx.xGrid(radialIdx);
    [xData, xLabelText] = diagnosticTimeAxis(ctx.tDiag, (0:nTime - 1) * ctx.diagSteps, ...
        ctx.L0, ctx.VA0, opt.timeAxis, '未知 amplitude 时间轴：%s。可选 ''ta''、''ms'' 或 ''steps''。');

    amplitudeScale = ctx.B0 * ctx.L0 * ctx.VA0 / (ctx.TeSample(1, 1) * 1000);
    amplitudeToPlot = reshape(amplitude(:, radialIdx, modeIdx), nTime, []);
    scaledAmplitude = amplitudeScale * amplitudeToPlot;
    logAmplitude = log10(scaledAmplitude);
    lnAmplitude = log(scaledAmplitude);
end

function [growthStart, growthEnd] = amplitudeGrowthInitialRange(xData, opt)
    xData = xData(:);
    assert(numel(xData) >= 2, 'Amplitude 增长率至少需要两个时间点。');

    axisMin = xData(1);
    axisMax = xData(end);
    if isfield(opt, 'growthRange') && numel(opt.growthRange) == 2 && all(isfinite(opt.growthRange))
        growthRange = sort(opt.growthRange);
    else
        growthRange = axisMin + [0.4, 0.6] * (axisMax - axisMin);
    end

    growthRange = min(max(growthRange, axisMin), axisMax);
    if growthRange(1) == growthRange(2)
        growthRange = [axisMin, axisMax];
    end

    growthStart = growthRange(1);
    growthEnd = growthRange(2);
end

function [growthRate, growthUnit, growthRange] = ...
    calculateAmplitudeGrowthRate(lnAmplitude, xData, ctx, opt)
    xData = xData(:);
    tDiag = ctx.tDiag(:);
    growthRange = sort(opt.growthRange);
    growthRange = min(max(growthRange, xData(1)), xData(end));
    assert(growthRange(2) > growthRange(1), 'Amplitude growthRange 必须包含两个不同的 x 位置。');

    tRange = interp1(xData, tDiag, growthRange, 'linear');
    logStart = interp1(xData, lnAmplitude, growthRange(1), 'linear');
    logEnd = interp1(xData, lnAmplitude, growthRange(2), 'linear');

    % 增长率使用自然对数，绘图曲线仍使用 log10。
    growthRate = (logEnd - logStart) / (tRange(2) - tRange(1));
    growthUnitOpt = getOptionValue(opt, 'growthUnit', '1/wa');
    [unitScale, growthUnit] = amplitudeGrowthUnitScale(growthUnitOpt, ctx.L0, ctx.VA0);
    growthRate = growthRate(:)' * unitScale;
end

function [unitScale, growthUnit] = amplitudeGrowthUnitScale(growthUnit, L0, VA0)
    growthUnitText = char(lower(growthUnit));
    switch growthUnitText
        case {'1/wa', '1/w_a', 'wa', 'omegaa', 'omega_a'}
            unitScale = 1.0;
            growthUnit = '1/wa';
        case {'1/s', 's^-1', 's^{-1}', 'per_s'}
            unitScale = VA0 / L0;
            growthUnit = '1/s';
        otherwise
            error('未知 amplitude 增长率单位：%s。可选 ''1/wa'' 或 ''1/s''。', char(growthUnit));
    end
end

function statusText = amplitudeGrowthStatus(radialX, radialIndex, growthRange, growthUnit, toroidalModeN, growthRate)
    growthItems = compose('n=%d: %.6g', toroidalModeN(:), growthRate(:));
    statusText = sprintf('x = %.6g, radialIndex = %d, growthRange = [%.6g, %.6g], growth(%s): %s', ...
        radialX, radialIndex, growthRange(1), growthRange(2), growthUnit, strjoin(cellstr(growthItems), ', '));
end

function plotData = buildMultipleFrequencyPlotData(RealMode, ImagMode, ctx, opt)
    [xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
        calculateMultipleFrequencyData(RealMode, ImagMode, ctx, opt);
    plotData = buildLinePlotData(xData, peakFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
        '$\mathrm{Short\mbox{-}Time\;Fourier\;Transform}$', compose('$n=%d$', toroidalModeN), ...
        sprintf('x = %.6g, radialIndex = %d, windowLength = %d, windowStep = %d, nFFT = %d', ...
        radialX, opt.radialIndex, opt.windowLength, opt.windowStep, multipleFrequencyNFFT(opt, opt.windowLength)));
    plotData.logText = sprintf('mode frequency: x=%.6g, radialIndex=%d, n=[%s]', ...
        radialX, opt.radialIndex, formatNumberList(toroidalModeN));
end

function plotMultipleFrequencyInteractive(RealMode, ImagMode, ctx, opt, dynamicUpdate)
    nTime = size(RealMode, 1);
    [windowLength0, windowStep0, nFFT0] = parseMultipleFrequencyOptions(opt, nTime);
    nFFTMax = max(nFFT0, 2 ^ nextpow2(nTime));
    windowStepMax = max(windowStep0, max(1, nTime - 1));
    controls = [ ...
        buildIntegerSliderControl('windowLength', 'windowLength', windowLength0, 2, nTime), ...
        buildIntegerSliderControl('windowStep', 'windowStep', windowStep0, 1, windowStepMax), ...
        buildIntegerSliderControl('nFFT', 'nFFT', nFFT0, 2, nFFTMax)];
    controls(3).minField = 'windowLength';

    plotInteractiveOptionsWithSliders('Short-Time Fourier Transform', controls, dynamicUpdate, opt, ...
        @(tempOpt) buildMultipleFrequencyPlotData(RealMode, ImagMode, ctx, tempOpt));
end

function nFFT = multipleFrequencyNFFT(opt, windowLength)

    nFFT = getOptionValue(opt, 'nFFT', 2 ^ nextpow2(windowLength));
end

function [windowLength, windowStep, nFFT, frequencyRangeHz, windowType, removeMean] = ...
    parseMultipleFrequencyOptions(opt, nTime)

    windowLength = requireIntegerInRange(opt.windowLength, 2, nTime, 'windowLength');
    windowStep = requireIntegerInRange(opt.windowStep, 1, Inf, 'windowStep');
    nFFT = requireIntegerInRange(multipleFrequencyNFFT(opt, windowLength), windowLength, Inf, 'nFFT');
    removeMean = requireLogicalScalar(opt.removeMean, 'removeMean');

    frequencyRangeHz = opt.frequencyRangeHz;
    if ~isempty(frequencyRangeHz)
        isValidRange = isnumeric(frequencyRangeHz) && isreal(frequencyRangeHz) && ...
            numel(frequencyRangeHz) == 2 && all(isfinite(frequencyRangeHz(:)));
        assert(isValidRange && frequencyRangeHz(1) <= frequencyRangeHz(2), ...
            'frequencyRangeHz 必须是按升序排列的两个有限数，或使用 []。');
        frequencyRangeHz = double(frequencyRangeHz(:)');
    end

    windowTypeInput = opt.windowType;
    isTextScalar = (ischar(windowTypeInput) && isrow(windowTypeInput)) || ...
        (isstring(windowTypeInput) && isscalar(windowTypeInput));
    assert(isTextScalar, 'windowType 必须为字符向量或字符串标量。');
    windowType = lower(string(windowTypeInput));
    assert(any(windowType == ["hann", "rect"]), ...
        'windowType 必须为 ''hann'' 或 ''rect''。');
end

function [xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
    calculateMultipleFrequencyData(RealMode, ImagMode, ctx, opt)
    nTime = size(RealMode, 1);
    [windowLength, windowStep, nFFT, frequencyRangeHz, windowType, removeMean] = ...
        parseMultipleFrequencyOptions(opt, nTime);

    [~, toroidalModeN, modeIdx] = parseModeN( ...
        opt.modeN, ctx.modeIndexAll, ctx.physicalNAll, 'Mode frequency');

    nRadial = size(RealMode, 2);
    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = ctx.xGrid(radialIdx);

    dtPhysical = (ctx.tDiag(2) - ctx.tDiag(1)) * ctx.L0 / ctx.VA0;
    sampleRateHz = 1 / dtPhysical;
    frequencyHz = (-floor(nFFT / 2):ceil(nFFT / 2) - 1)' * sampleRateHz / nFFT;
    frequencyMask = true(size(frequencyHz));
    if ~isempty(frequencyRangeHz)
        frequencyMask = frequencyHz >= frequencyRangeHz(1) & frequencyHz <= frequencyRangeHz(2);
    end
    assert(any(frequencyMask), 'frequencyRangeHz 未包含任何 FFT 频率点。');

    switch char(windowType)
        case 'hann'
            windowData = 0.5 - 0.5 * cos(2 * pi * (0:windowLength - 1)' / (windowLength - 1));
        case 'rect'
            windowData = ones(windowLength, 1);
    end

    windowStart = 1:windowStep:(nTime - windowLength + 1);
    nWindow = numel(windowStart);
    peakFrequencyHz = nan(nWindow, numel(modeIdx));
    centerIndex = windowStart + floor((windowLength - 1) / 2);

    for iMode = 1:numel(modeIdx)
        signal = RealMode(:, radialIdx, modeIdx(iMode)) + 1i * ImagMode(:, radialIdx, modeIdx(iMode));

        for iWindow = 1:nWindow
            idx = windowStart(iWindow):(windowStart(iWindow) + windowLength - 1);
            segment = signal(idx);
            if removeMean
                segment = segment - mean(segment);
            end

            spectrum = fftshift(fft(segment(:) .* windowData, nFFT));
            powerSpectrum = abs(spectrum) .^ 2;
            powerSpectrum(~frequencyMask) = -inf;
            [~, peakIdx] = max(powerSpectrum);
            peakFrequencyHz(iWindow, iMode) = -frequencyHz(peakIdx);
        end
    end

    [xData, xLabelText] = diagnosticTimeAxis( ...
        ctx.tDiag(centerIndex), (centerIndex - 1) * ctx.diagSteps, ctx.L0, ctx.VA0, opt.timeAxis, ...
        '未知 mode frequency 时间轴：%s。可选 ''ta''、''ms'' 或 ''steps''。');
end

function plotData = buildPhaseFrequencyPlotData(RealMode, ImagMode, ctx, opt)
    [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
        calculatePhaseFrequencyData(RealMode, ImagMode, ctx, opt);
    plotData = buildLinePlotData(xData, phaseFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
        '', compose('$n=%d$', toroidalModeN), ...
        sprintf('x = %.6g, radialIndex = %d, phaseStep = %d, smoothWindow = %d', ...
        radialX, opt.radialIndex, opt.phaseStep, getOptionValue(opt, 'smoothWindow', 1)));
    plotData.logText = sprintf('phase frequency: x=%.6g, radialIndex=%d, n=[%s], phaseStep=%d', ...
        radialX, opt.radialIndex, formatNumberList(toroidalModeN), opt.phaseStep);
end

function plotPhaseFrequencyInteractive(RealMode, ImagMode, ctx, opt, dynamicUpdate)
    nTime = size(RealMode, 1);
    [phaseStep0, smoothWindow0] = parsePhaseFrequencyOptions(opt, nTime);
    smoothWindowMax = max(smoothWindow0, max(1, min(1001, nTime - 1)));
    controls = [ ...
        buildIntegerSliderControl('phaseStep', 'phaseStep', phaseStep0, 1, nTime - 1), ...
        buildIntegerSliderControl('smoothWindow', 'smoothWindow', smoothWindow0, 1, smoothWindowMax)];
    plotInteractiveOptionsWithSliders('Phase Frequency', controls, dynamicUpdate, opt, ...
        @(tempOpt) buildPhaseFrequencyPlotData(RealMode, ImagMode, ctx, tempOpt));
end

function [phaseStep, smoothWindow, amplitudeFloor] = parsePhaseFrequencyOptions(opt, nTime)

    phaseStep = requireIntegerInRange(opt.phaseStep, 1, nTime - 1, 'phaseStep');
    smoothWindow = requireIntegerInRange(getOptionValue(opt, 'smoothWindow', 1), 1, Inf, 'smoothWindow');
    amplitudeFloor = requireFiniteScalarInRange( ...
        getOptionValue(opt, 'amplitudeFloor', 0), 0, Inf, 'amplitudeFloor');
end

function [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
    calculatePhaseFrequencyData(RealMode, ImagMode, ctx, opt)
    nTime = size(RealMode, 1);
    [phaseStep, smoothWindow, amplitudeFloor] = parsePhaseFrequencyOptions(opt, nTime);

    [~, toroidalModeN, modeIdx] = parseModeN( ...
        opt.modeN, ctx.modeIndexAll, ctx.physicalNAll, 'Phase frequency');

    nRadial = size(RealMode, 2);
    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = ctx.xGrid(radialIdx);

    dtPhysical = (ctx.tDiag(2) - ctx.tDiag(1)) * ctx.L0 / ctx.VA0;
    startIndex = (1:(nTime - phaseStep))';
    endIndex = startIndex + phaseStep;
    phaseFrequencyHz = nan(numel(startIndex), numel(modeIdx));

    for iMode = 1:numel(modeIdx)
        signal = RealMode(:, radialIdx, modeIdx(iMode)) + 1i * ImagMode(:, radialIdx, modeIdx(iMode));
        signal = signal(:);
        phase = unwrap(angle(signal));

        phaseFrequencyHz(:, iMode) = -(phase(endIndex) - phase(startIndex)) / (2 * pi * dtPhysical * phaseStep);

        if amplitudeFloor > 0
            signalAmplitude = abs(signal);
            lowAmplitude = signalAmplitude(startIndex) < amplitudeFloor | signalAmplitude(endIndex) < amplitudeFloor;
            phaseFrequencyHz(lowAmplitude, iMode) = NaN;
        end

        effectiveSmoothWindow = min(smoothWindow, numel(startIndex));
        if effectiveSmoothWindow > 1
            phaseFrequencyHz(:, iMode) = smoothdata(phaseFrequencyHz(:, iMode), 'movmean', effectiveSmoothWindow);
        end
    end

    tACenter = 0.5 * (ctx.tDiag(startIndex) + ctx.tDiag(endIndex));
    stepCenter = ((startIndex - 1) + phaseStep / 2) * ctx.diagSteps;
    [xData, xLabelText] = diagnosticTimeAxis(tACenter, stepCenter, ctx.L0, ctx.VA0, opt.timeAxis, ...
        '未知 phase frequency 时间轴：%s。可选 ''ta''、''ms'' 或 ''steps''。');
end

function plotData = buildContourFrequencyPlotData(amplitude, RealMode, ctx, opt)

    nTime = size(RealMode, 1);

    [~, toroidalModeN, modeIdx] = parseModeN( ...
        opt.modeN, ctx.modeIndexAll, ctx.physicalNAll, 'Contour frequency');
    assert(isscalar(modeIdx), 'Contour frequency modeN 一次只能选择一个物理环向模数。');

    timeIndexRange = opt.timeIndexRange;
    if isempty(timeIndexRange)
        timeIndex = 1:nTime;
    else
        timeIndexRange = sort(reshape(double(timeIndexRange), 1, []));
        isValidRange = numel(timeIndexRange) == 2 && all(isfinite(timeIndexRange)) && ...
            all(timeIndexRange == floor(timeIndexRange));
        assert(isValidRange, ...
            'Contour frequency timeIndexRange 必须为两个有限整数。');
        assert(timeIndexRange(1) >= 1 && timeIndexRange(2) <= nTime && timeIndexRange(2) > timeIndexRange(1), ...
            'Contour frequency timeIndexRange 必须位于 [1, %d] 且至少包含两个时间点。', nTime);
        timeIndex = timeIndexRange(1):timeIndexRange(2);
    end

    amplitude2D = amplitude(timeIndex, :, modeIdx);
    field2D = RealMode(timeIndex, :, modeIdx);
    tDiag = ctx.tDiag(timeIndex);
    nTime = numel(timeIndex);

    fieldRMax = max(amplitude2D, [], 2);
    fieldRMax(fieldRMax == 0) = NaN;
    normalizedField = bsxfun(@rdivide, field2D, fieldRMax);
    normalizedField(~isfinite(normalizedField)) = 0;

    spectrum = abs(fft(normalizedField / nTime, [], 1));
    nFrequency = floor(nTime / 2) + 1;
    frequencyIntensity = spectrum(1:nFrequency, :);
    if mod(nTime, 2) == 0 && nFrequency > 2
        frequencyIntensity(2:end - 1, :) = 2 * frequencyIntensity(2:end - 1, :);
    elseif mod(nTime, 2) == 1 && nFrequency > 1
        frequencyIntensity(2:end, :) = 2 * frequencyIntensity(2:end, :);
    end

    globalMax = max(frequencyIntensity(:));
    if globalMax > 0
        frequencyIntensity = frequencyIntensity ./ globalMax;
    end

    omegaVec = (0:nFrequency - 1)' / ((tDiag(2) - tDiag(1)) * nTime) * 2 * pi;
    switch char(lower(opt.frequencyUnit))
        case {'omegaa', 'omega_a', 'normalized'}
            frequencyVec = omegaVec;
            yLabelText = '$\omega/\omega_A$';
            frequencyUnit = 'omegaA';
        case {'khz'}
            frequencyVec = omegaVec * ctx.VA0 / ctx.L0 / (2 * pi) / 1000;
            yLabelText = '$f/\mathrm{kHz}$';
            frequencyUnit = 'kHz';
        otherwise
            error('未知 frequencyUnit：%s。可选 ''omegaA'' 或 ''kHz''。', opt.frequencyUnit);
    end

    plotData = buildMapPlotData(ctx.xGrid(:).', frequencyVec, frequencyIntensity, '$r/a$', ...
        yLabelText, '$I/I_{\mathrm{max}}$', '$\mathrm{Frequency\;contour}$');
    plotData.logText = sprintf('contour frequency: n=%d, timeIndexRange=[%d %d], frequencyUnit=%s', ...
        toroidalModeN, timeIndex(1), timeIndex(end), frequencyUnit);
end

function plotData = buildSingleSignalPlotData(singlePointPhi, ctx, opt)
    [xData, yData, xLabelText, radialX] = calculateSingleSignalData(singlePointPhi, ctx, opt);
    plotData = buildLinePlotData(xData, yData, xLabelText, '', '$\log(|\delta\phi|)$', [], ...
        sprintf('x = %.6g, radialIndex = %d', radialX, opt.radialIndex));
    plotData.logText = sprintf('single signal: x=%.6g, radialIndex=%d', radialX, opt.radialIndex);
end

function plotSingleSignalInteractive(singlePointPhi, ctx, opt, dynamicUpdate)
    nRadial = size(singlePointPhi, 2);
    radialIndex0 = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    controls = buildIntegerSliderControl('radialIndex', 'radialIndex', radialIndex0, 1, nRadial);
    plotInteractiveOptionsWithSliders('singlePointPhi', controls, dynamicUpdate, opt, ...
        @(tempOpt) buildSingleSignalPlotData(singlePointPhi, ctx, tempOpt));
end

function [xData, yData, xLabelText, radialX] = calculateSingleSignalData(singlePointPhi, ctx, opt)
    nTime = size(singlePointPhi, 1);

    nRadial = size(singlePointPhi, 2);

    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = ctx.xGrid(radialIdx);
    assert(isscalar(opt.logFloor) && opt.logFloor > 0, 'single signal logFloor 必须为正数。');

    [xData, xLabelText] = diagnosticTimeAxis(ctx.tDiag, (0:nTime - 1) * ctx.diagSteps, ...
        ctx.L0, ctx.VA0, opt.timeAxis, 'single signal 的 timeAxis 必须为 "ta"、"ms" 或 "steps"。当前值：%s。');

    yData = log(max(abs(singlePointPhi(:, radialIdx)), opt.logFloor));
end

function plotIndexOptionInteractive(figName, optionField, nIndex, dynamicUpdate, opt, buildPlotDataFcn)
    initialIndex = parseIndex(opt.(optionField), nIndex, optionField);
    controls = buildIntegerSliderControl(optionField, optionField, initialIndex, 1, nIndex);
    plotInteractiveOptionsWithSliders(figName, controls, dynamicUpdate, opt, buildPlotDataFcn);
end

function plotData = buildEparaPlotData(Epara, EparaES, ctx, opt)
    [xData, yData, timeIdx] = calculateEparaData(Epara, EparaES, ctx, opt);
    plotData = buildLinePlotData(xData, yData, '$r/a$', '', '', ...
        {'$E_{\parallel}$', '$E_{\parallel}^{\mathrm{ES}}$', '$\partial \delta A_{\parallel}/\partial t$'}, ...
        sprintf('timeIndex = %d', timeIdx));
    plotData.logText = sprintf('Epara/EparaES: timeIndex=%d', timeIdx);
end

function [xData, yData, timeIdx] = calculateEparaData(Epara, EparaES, ctx, opt)
    timeIdx = parseTimeIndex(opt.timeIndex, size(Epara, 1));
    xData = ctx.rho(:, 1);
    yData = [Epara(timeIdx, :)', EparaES(timeIdx, :)', (EparaES(timeIdx, :) - Epara(timeIdx, :))'];
end

function plotData = buildShearingTimeSlicePlotData(Shearing, ctx, opt)
    [xData, yData, timeIdx] = calculateShearingTimeSliceData(Shearing, ctx, opt);
    plotData = buildLinePlotData(xData, yData, '$x$', '$\gamma_E\;[\mathrm{s}^{-1}]$', ...
        '$\mathrm{Hahm\mbox{-}Burrell\;shearing}$', [], sprintf('timeIndex = %d', timeIdx));
    plotData.logText = sprintf('Shearing: timeIndex=%d', timeIdx);
end

function [xData, yData, timeIdx] = calculateShearingTimeSliceData(Shearing, ctx, opt)
    timeIdx = parseTimeIndex(opt.timeIndex, size(Shearing, 1));
    xData = ctx.xGrid(:);
    yData = Shearing(timeIdx, :)';
end

function plotData = buildShearingRadialTracePlotData(Shearing, ctx, opt)
    [xData, yData, radialIdx, radialX] = calculateShearingRadialTraceData(Shearing, ctx, opt);
    plotData = buildLinePlotData(xData, yData, '$t_a$', '$\gamma_E\;[\mathrm{s}^{-1}]$', ...
        '$\mathrm{Hahm\mbox{-}Burrell\;shearing}$', [], ...
        sprintf('x = %.6g, radialIndex = %d', radialX, radialIdx));
    plotData.logText = sprintf('Shearing trace: x=%.6g, radialIndex=%d', radialX, radialIdx);
end

function [xData, yData, radialIdx, radialX] = calculateShearingRadialTraceData(Shearing, ctx, opt)
    radialIdx = parseIndex(opt.radialIndex, size(Shearing, 2), 'radialIndex');
    radialX = ctx.xGrid(radialIdx);
    xData = ctx.tDiag(:);
    yData = Shearing(:, radialIdx);
end

function plotData = buildShearingMapPlotData(Shearing, ctx)
    nTime = size(Shearing, 1);
    nRadial = size(Shearing, 2);
    plotData = buildMapPlotData(ctx.xGrid(:).', ctx.tDiag(:), Shearing, '$x$', '$t_a$', ...
        '$\gamma_E\;[\mathrm{s}^{-1}]$', '$\mathrm{Hahm\mbox{-}Burrell\;shearing}$');
    plotData.logText = sprintf('Shearing map: nTime=%d, gridNx=%d', nTime, nRadial);
end

function plotData = buildZonalDrivePlotData(MaxwellDrive, ReynoldsDrive, ZonalDrive, ctx, opt)
    [xData, yData, timeIdx] = calculateZonalDriveData( ...
        MaxwellDrive, ReynoldsDrive, ZonalDrive, ctx, opt);
    plotData = buildLinePlotData(xData, yData, '$r/a$', '', '$\mathrm{zonal\;flow\;drive}$', ...
        {'$\mathrm{MaxwellDrive}$', '$\mathrm{ReynoldsDrive}$', '$\mathrm{ZonalDrive}$'}, ...
        sprintf('timeIndex = %d', timeIdx));
    plotData.logText = sprintf('ZF drive: timeIndex=%d', timeIdx);
end

function [xData, yData, timeIdx] = calculateZonalDriveData( ...
    MaxwellDrive, ReynoldsDrive, ZonalDrive, ctx, opt)
    timeIdx = parseTimeIndex(opt.timeIndex, size(MaxwellDrive, 1));
    xData = ctx.rho(:, 1);
    yData = [MaxwellDrive(timeIdx, :)', ReynoldsDrive(timeIdx, :)', ZonalDrive(timeIdx, :)'];
end

function applyLinePlotColorOrder(axHandle)
    set(axHandle, 'ColorOrder', getPreferredLineColors(), 'NextPlot', 'replacechildren');
end

function value = clampInteger(value, minValue, maxValue)
    value = round(value);
    value = min(max(value, minValue), maxValue);
end

function value = clampSliderValue(value, control)
    if control.isInteger
        value = clampInteger(value, control.min, control.max);
    else
        value = min(max(value, control.min), control.max);
    end
end

function control = buildIntegerSliderControl(fieldName, labelText, value, minValue, maxValue)
    control = struct();
    control.field = fieldName;
    control.label = labelText;
    control.value = clampInteger(value, minValue, maxValue);
    control.min = minValue;
    control.max = maxValue;
    control.isInteger = true;
    control.nStep = [];
    control.minField = '';
end

function control = buildNumericSliderControl(fieldName, labelText, value, minValue, maxValue, nStep)
    control = struct();
    control.field = fieldName;
    control.label = labelText;
    control.value = min(max(value, minValue), maxValue);
    control.min = minValue;
    control.max = maxValue;
    control.isInteger = false;
    control.nStep = nStep;
    control.minField = '';
end

function plotData = buildLinePlotData(xData, yData, xLabelText, yLabelText, titleText, legendText, statusText)
    plotData = struct();
    plotData.xVec = xData;
    plotData.yVec = yData;
    plotData.xlabelText = xLabelText;
    plotData.ylabelText = yLabelText;
    plotData.titleText = titleText;
    plotData.legendText = legendText;
    plotData.status = statusText;
end

function plotData = buildMapPlotData(xData, yData, zData, xLabelText, yLabelText, colorbarLabel, titleText)
    plotData = struct();
    plotData.xVec = xData;
    plotData.yVec = yData;
    plotData.Z = zData;
    plotData.xlabelText = xLabelText;
    plotData.ylabelText = yLabelText;
    plotData.colorbarLabel = colorbarLabel;
    plotData.titleText = titleText;
end

function drawPlot(plotData, opt, renderFcn)
    figure;
    axHandle = gca;
    renderFcn(axHandle, plotData, opt);
end

function drawStaticDiagnostic(plotData, opt, renderFcn)
    drawPlot(plotData, opt, renderFcn);
    fprintf('[plot] %s\n', plotData.logText);
end

function renderMapPlot(axHandle, plotData, opt)
    figHandle = ancestor(axHandle, 'figure');
    delete(findall(figHandle, 'Type', 'ColorBar'));
    delete(allchild(axHandle));

    contourLevels = getOptionValue(opt, 'contourLevels', []);
    if ~isempty(contourLevels)
        showContourLine = requireLogicalScalar(getOptionValue(opt, 'showContourLine', true), 'showContourLine');
        if ~showContourLine
            contourf(axHandle, plotData.xVec, plotData.yVec, plotData.Z, contourLevels, 'LineStyle', 'none');
        else
            contourf(axHandle, plotData.xVec, plotData.yVec, plotData.Z, contourLevels);
        end
    else
        imagesc(axHandle, plotData.xVec, plotData.yVec, plotData.Z);
    end
    set(axHandle, 'YDir', 'normal');
    colormap(axHandle, jet(256));
    colorbarHandle = colorbar(axHandle);
    colorbarHandle.FontName = 'Times New Roman';
    colorbarHandle.FontSize = 12;
    ylabel(colorbarHandle, plotData.colorbarLabel, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', 14);

    finiteZ = plotData.Z(isfinite(plotData.Z));
    if ~isempty(contourLevels)
        clim(axHandle, [min(contourLevels), max(contourLevels)]);
    elseif isempty(finiteZ)
        clim(axHandle, [0, 1]);
    elseif min(finiteZ) == max(finiteZ)
        centerValue = finiteZ(1);
        halfWidth = max(1, 0.05 * abs(centerValue));
        clim(axHandle, centerValue + [-halfWidth, halfWidth]);
    else
        clim(axHandle, [min(finiteZ), max(finiteZ)]);
    end

    applyOptionalAxesSetting(axHandle, opt, 'yLim', @ylim);
    applyOptionalAxesSetting(axHandle, opt, 'yTicks', @yticks);

    box(axHandle, 'on');
    applyPlotAxesStyle(axHandle);
    applyPlotText(axHandle, plotData);
end

function renderLinePlot(axHandle, plotData, opt)
    delete(allchild(axHandle));
    applyLinePlotColorOrder(axHandle);
    lineHandles = plot(axHandle, plotData.xVec, plotData.yVec, 'LineWidth', 1.5);
    box(axHandle, 'on');
    applyPlotAxesStyle(axHandle);

    applyOptionalAxesSetting(axHandle, opt, 'yLim', @ylim);
    applyOptionalAxesSetting(axHandle, opt, 'yTicks', @yticks);
    applyOptionalAxesSetting(axHandle, opt, 'xLim', @xlim);
    applyOptionalAxesSetting(axHandle, opt, 'xTicks', @xticks);

    if isfield(plotData, 'xLines') && ~isempty(plotData.xLines)
        drawVerticalReferenceLines(axHandle, plotData.xLines);
    end

    applyPlotText(axHandle, plotData);
    if ~isempty(plotData.legendText)
        legend(lineHandles, plotData.legendText, 'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', 'FontSize', 12, 'Location', 'best');
    end
end

function applyPlotText(axHandle, plotData)
    xlabel(axHandle, plotData.xlabelText, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', 14);
    if ~isempty(plotData.ylabelText)
        ylabel(axHandle, plotData.ylabelText, 'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', 'FontSize', 14);
    end
    if ~isempty(plotData.titleText)
        title(axHandle, plotData.titleText, 'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', 'FontSize', 14);
    end
end

function drawVerticalReferenceLines(axHandle, xLines)
    yLimits = ylim(axHandle);
    holdState = ishold(axHandle);
    hold(axHandle, 'on');
    for iLine = 1:numel(xLines)
        plot(axHandle, [xLines(iLine), xLines(iLine)], yLimits, 'k--', ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    if ~holdState
        hold(axHandle, 'off');
    end
end

function applyOptionalAxesSetting(axHandle, opt, fieldName, setterFcn)
    value = getOptionValue(opt, fieldName, []);
    if ~isempty(value)
        setterFcn(axHandle, value);
    end
end

function plotInteractiveLineWithSliders(figName, controls, dynamicUpdate, opt, buildPlotDataFcn)
    figHandle = figure('Name', figName, 'Color', 'w', 'Position', [120, 120, 920, 620]);
    nControl = numel(controls);
    controlBottom = 0.035;
    controlSpacing = 0.045;
    statusBottom = controlBottom + controlSpacing * max(nControl, 1) + 0.010;
    axesBottom = statusBottom + 0.105;
    axHandle = axes('Parent', figHandle, 'Units', 'normalized', ...
        'Position', [0.10, axesBottom, 0.86, 0.93 - axesBottom]);
    statusText = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.10, statusBottom, 0.86, 0.035], 'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', 'FontName', 'Times New Roman', 'FontSize', 12);
    sliderLabels = gobjects(nControl, 1);
    sliders = gobjects(nControl, 1);

    for iControl = 1:nControl
        yPos = controlBottom + controlSpacing * (nControl - iControl);
        sliderLabels(iControl) = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.10, yPos - 0.01, 0.19, 0.035], ...
            'String', sliderLabelText(controls(iControl), controls(iControl).value), ...
            'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', 12);
        sliders(iControl) = uicontrol(figHandle, 'Style', 'slider', 'Units', 'normalized', ...
            'Position', [0.30, yPos, 0.61, 0.025], 'Min', controls(iControl).min, ...
            'Max', controls(iControl).max, 'Value', controls(iControl).value, ...
            'SliderStep', sliderStepForControl(controls(iControl)), ...
            'Callback', @refreshPlot);
    end

    if dynamicUpdate
        dynamicListeners = {};
        try
            for iControl = 1:nControl
                dynamicListeners{end + 1} = addlistener( ...
                    sliders(iControl), 'ContinuousValueChange', @refreshPlot); %#ok<AGROW>
            end
            setappdata(figHandle, 'dynamicSliderListeners', dynamicListeners);
        catch
            warning('visualizeMHD:DynamicSliderUnavailable', ...
                '动态滑块更新不可用，将使用回调更新。');
            dynamicUpdate = false;
        end
    end

    refreshPlot();

    function refreshPlot(varargin)
        if ~isgraphics(axHandle)
            return;
        end

        values = struct();
        for iSlider = 1:nControl
            control = controls(iSlider);
            sliderValue = get(sliders(iSlider), 'Value');
            if ~isempty(control.minField)
                control.min = max(control.min, values.(control.minField));
                if sliderValue < control.min
                    sliderValue = control.min;
                    set(sliders(iSlider), 'Value', sliderValue);
                end
                set(sliders(iSlider), 'Min', control.min);
            end
            sliderValue = clampSliderValue(sliderValue, control);
            values.(controls(iSlider).field) = sliderValue;
            if ~dynamicUpdate
                set(sliders(iSlider), 'Value', sliderValue);
            end
            set(sliderLabels(iSlider), 'String', sliderLabelText(controls(iSlider), sliderValue));
        end

        data = buildPlotDataFcn(values);
        renderLinePlot(axHandle, data, opt);
        set(statusText, 'String', data.status);
    end
end

function plotInteractiveOptionsWithSliders(figName, controls, dynamicUpdate, opt, buildPlotDataFcn)
    plotInteractiveLineWithSliders(figName, controls, dynamicUpdate, opt, @applySliderValues);

    function data = applySliderValues(values)
        tempOpt = opt;
        valueNames = fieldnames(values);
        for iValue = 1:numel(valueNames)
            tempOpt.(valueNames{iValue}) = values.(valueNames{iValue});
        end
        data = buildPlotDataFcn(tempOpt);
    end
end

function text = sliderLabelText(control, value)
    if control.isInteger
        text = sprintf('%s = %d', control.label, value);
    else
        text = sprintf('%s = %.6g', control.label, value);
    end
end

function step = sliderStepForControl(control)
    if control.isInteger
        step = calculateSliderStep(control.max - control.min + 1);
    elseif isfield(control, 'nStep') && ~isempty(control.nStep) && control.nStep > 1
        step = calculateSliderStep(control.nStep);
    else
        step = [0.005, 0.05];
    end
end

function applyPlotAxesStyle(axHandle)
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 12, ...
        'TickLabelInterpreter', 'latex', ...
        'XGrid', 'on', 'YGrid', 'on', 'GridLineWidth', 1.0, 'GridAlpha', 0.3, ...
        'MinorGridAlpha', 0.3, 'GridColor', 'k', 'MinorGridColor', 'k', ...
        'GridLineStyle', '-', 'MinorGridLineStyle', '-');
end

function step = calculateSliderStep(maxValue)
    if maxValue <= 1
        step = [1, 1];
    else
        smallStep = 1 / (maxValue - 1);
        largeStep = min(1, max(1, round(maxValue / 20)) / (maxValue - 1));
        step = [smallStep, largeStep];
    end
end

function colors = getPreferredLineColors()
    colors = [ ...
        1.0000, 0.0000, 0.0000;  % #FF0000
        0.0000, 0.0000, 1.0000;  % #0000FF
        0.0000, 0.5020, 0.0000;  % #008000
        0.0000, 1.0000, 0.0000;  % #00FF00
        1.0000, 0.6471, 0.0000;  % #FFA500
        0.5020, 0.0000, 0.5020;  % #800080
        0.0000, 1.0000, 1.0000;  % #00FFFF
        1.0000, 0.7529, 0.7961]; % #FFC0CB
end
