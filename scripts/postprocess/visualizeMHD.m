%% cuGMEC MHD 诊断统一可视化脚本

%{
功能：
统一可视化 MHD 诊断和输出场量。

输入目录：
将相关输入、诊断和输出文件放在同一个 inputDir 中。

必须包含：
cuGMEC_param.h, standard2D.mat, plot2D.mat, normalization2D.mat, NTP.mat

按开关读取：
ifDiagAmplitude, ifDiagFrequency, ifDiagEparallel, ifDiagZFDrive,
ifOutputPhi, ifOutputA, ifOutputdNe, ifOutputdTe,
ifOutputdPi, ifOutputdPa, ifOutputdPb

可能包含：
amplitude.bin, RealMode.bin, ImagMode.bin, frequency.bin, Epara.bin,
EparaES.bin, MaxwellDrive.bin, ReynoldsDrive.bin, ZonalDrive.bin,
Phi.bin, A.bin, dNe.bin, dTe.bin, dPi.bin, dPa.bin, dPb.bin,
totalPhi.bin, totalA.bin, totaldNe.bin, totaldTe.bin,
totaldPi.bin, totaldPa.bin, totaldPb.bin
%}

%% 用户设置

clear; close all;

inputDir = 'C:\Users\Desktop\test';

paramFile = fullfile(inputDir, 'cuGMEC_param.h');
standardFile = fullfile(inputDir, 'standard2D.mat');
plotFile = fullfile(inputDir, 'plot2D.mat');
normalizationFile = fullfile(inputDir, 'normalization2D.mat');
NTPFile = fullfile(inputDir, 'NTP.mat');
bsiDir = fullfile(inputDir, 'BSI');
scriptDir = fileparts(mfilename('fullpath'));
repoBsiDir = fullfile(scriptDir, '..', 'preprocess', 'BSI');

%% 读取所有输入

assert(isfolder(inputDir), '缺少输入目录：%s', inputDir);
assert(isfile(paramFile), '缺少参数文件：%s', paramFile);
assert(isfile(standardFile), '缺少 MAT 文件：%s', standardFile);
assert(isfile(plotFile), '缺少 MAT 文件：%s', plotFile);
assert(isfile(normalizationFile), '缺少 MAT 文件：%s', normalizationFile);
assert(isfile(NTPFile), '缺少 MAT 文件：%s', NTPFile);
if isfolder(bsiDir)
    addpath(bsiDir);
elseif isfolder(repoBsiDir)
    addpath(repoBsiDir);
else
    fprintf('[skip] BSI: 目录不存在\n');
end

paramText = fileread(paramFile);
mhdInput = loadMHDInputData(standardFile, plotFile, normalizationFile, NTPFile);
meta = readMHDMetadata(paramText);
mhdFieldGeom = mhdFieldPlotGeometry(mhdInput.q, mhdInput.theta_pest, mhdInput.rho, ...
    mhdInput.qplot, mhdInput.rhoplot, mhdInput.Rplot, mhdInput.Zplot, ...
    meta.yGrid, meta.zGrid, meta.tubes);
mhdData = readAllMHDDiagnostics(inputDir, meta);

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
    'colormapIndex', 1, ...
    'bsplineOrder', 4, ...
    'titleFontSize', 14, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled      : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
toroidalAngle : 绘图所在环向角 phi。
colormapIndex : 色表序号，可选 1-5。
bsplineOrder  : B 样条阶数。
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.field = runMHDFieldPlot(inputDir, meta, mhdFieldGeom, mhdFieldOpt);
[tmpShifted, tmpAligned, hasFieldResult] = extractShiftedAligned(mhdWorkspace.field);
if hasFieldResult
    shifted = tmpShifted;
    aligned = tmpAligned;
end

%% 可视化单时间MHD极向模数

mhdPoloidalModeOpt = struct( ...
    'enabled', true, ...
    'names', "Phi", ...
    'modeN', mhdFieldOpt.modeN, ...
    'toroidalAngle', mhdFieldOpt.toroidalAngle, ...
    'mRange', [1, 20], ...
    'plotThreshold', 1e-8, ...
    'minFFT', 0, ...
    'normalize', true, ...
    'radialAxis', 'rho', ...
    'bsplineOrder', mhdFieldOpt.bsplineOrder, ...
    'xLim', [0, 1], ...
    'xTicks', 0:0.2:1, ...
    'yLim', [0, 1], ...
    'yTicks', 0:0.25:1, ...
    'titleFontSize', 14, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled       : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
toroidalAngle : 做极向截面的环向角 phi。
mRange        : 绘制的极向模数 m 范围；示例 [1 20]。
plotThreshold : 只绘制峰值超过全局峰值该比例的 m。
minFFT        : 只绘制归一化峰值超过该值的 m。
normalize     : 是否用所选 m 范围内的最大 FFT 幅值归一化。
radialAxis    : 'rho' 表示横坐标为 sqrt(s)；'psi' 表示近似使用 rho^2。
bsplineOrder  : B 样条阶数。
xLim/xTicks   : 横坐标范围和刻度；[] 表示自动。
yLim/yTicks   : 纵坐标范围和刻度；[] 表示自动。
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.poloidalMode = runMHDPoloidalModePlot(inputDir, meta, mhdInput, mhdFieldGeom, mhdPoloidalModeOpt);

%% 可视化 amplitude

mhdAmplitudeOpt = struct( ...
    'enabled', meta.switch.ifDiagAmplitude, ...
    'timeAxis', 'ta', ...
    'radialIndex', 90, ...
    'modeN', [], ...
    'growthRange', [], ...
    'growthUnit', '1/wa', ...
    'yLim', [-10 0], ...
    'yTicks', -8:2:0, ...
    'interactive', 2, ...
    'titleFontSize', 16, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
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
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.amplitude = runMHDAmplitudePlot(mhdData, meta, mhdInput, mhdAmplitudeOpt);


%% 可视化短时傅里叶频率

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
    'interactive', 2, ...
    'titleFontSize', 14, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
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
titleFontSize   : 标题字号。
labelFontSize   : 坐标轴标签字号。
axisFontSize    : 坐标轴字号。
%}
mhdWorkspace.multipleFrequency = runMHDMultipleFrequencyPlot(mhdData, meta, mhdInput, mhdMultipleFrequencyOpt);


%% 可视化相位频率

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
    'interactive', 0, ...
    'titleFontSize', 14, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
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
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.phaseFrequency = runMHDPhaseFrequencyPlot(mhdData, meta, mhdInput, mhdPhaseFrequencyOpt);

%% 可视化 phase 方法二维频率强度图

mhdPhaseFrequencyMapOpt = struct( ...
    'enabled', meta.switch.ifDiagAmplitude, ...
    'timeAxis', 'ms', ...
    'radialIndex', 20, ...
    'modeN', meta.physicalNAll(1), ...
    'phaseStep', 2000, ...
    'smoothWindow', 1, ...
    'amplitudeFloor', 0, ...
    'frequencyRangeHz', [], ...
    'frequencyBins', 256, ...
    'frequencySmoothingBins', 3, ...
    'intensityPower', 2, ...
    'normalizeIntensity', true, ...
    'logIntensity', false, ...
    'interactive', 2, ...
    'titleFontSize', 14, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled               : 是否绘图。
timeAxis              : 'ta', 'ms' 或 'steps'。
radialIndex           : 固定径向索引（1 到 gridNx）。
modeN                 : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
phaseStep             : 相位差分间隔，单位为诊断点数。
smoothWindow          : phase 频率曲线平滑窗口；1 表示不平滑。
amplitudeFloor        : 信号幅度低于该值时不沉积强度。
frequencyRangeHz      : 频率轴范围，单位 Hz；[] 表示自动范围。
frequencyBins         : 频率方向 bin 数量。
frequencySmoothingBins: 频率方向强度平滑 bin 数；1 表示不平滑。
intensityPower        : 频率 bin 权重使用 abs(signal)^intensityPower。
normalizeIntensity    : true 时每个时间列归一化，颜色表示该频率占比。
logIntensity          : true 时颜色绘制 log10 强度。
interactive           : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
titleFontSize         : 标题字号。
labelFontSize         : 坐标轴标签字号。
axisFontSize          : 坐标轴字号。
%}
mhdWorkspace.phaseFrequencyMap = runMHDPhaseFrequencyMapPlot(mhdData, meta, mhdInput, mhdPhaseFrequencyMapOpt);

%% 可视化单点 Phi 信号

mhdSingleSignalOpt = struct( ...
    'enabled', meta.switch.ifDiagFrequency, ...
    'timeAxis', 'ta', ...
    'radialIndex', 128, ...
    'logFloor', realmin, ...
    'interactive', 2, ...
    'titleFontSize', 16, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled       : 是否绘图。
timeAxis      : 'ta', 'ms' 或 'steps'。
radialIndex   : 固定径向索引（1 到 gridNx）。
logFloor      : 避免 log(0)。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.singleSignal = runMHDSingleSignalPlot(mhdData, meta, mhdInput, mhdSingleSignalOpt);

%% 可视化 Epara / EparaES

mhdEparaOpt = struct( ...
    'enabled', meta.switch.ifDiagEparallel, ...
    'timeIndex', 12000, ...
    'interactive', 2, ...
    'titleFontSize', 16, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled       : 是否绘图。
timeIndex     : 固定诊断时间索引（1 到 nDiagTime）。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.epara = runMHDEparaPlot(mhdData, mhdInput, mhdEparaOpt);

%% 可视化 MaxwellDrive / ReynoldsDrive / ZonalDrive

mhdZonalDriveOpt = struct( ...
    'enabled', meta.switch.ifDiagZFDrive, ...
    'timeIndex', 15000, ...
    'interactive', 2, ...
    'titleFontSize', 16, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled       : 是否绘图。
timeIndex     : 固定诊断时间索引（1 到 nDiagTime）。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.zonalDrive = runMHDZonalDrivePlot(mhdData, mhdInput, mhdZonalDriveOpt);

%% 可视化 totalPhi / totalA / totaldNe / totaldTe / totaldPi / totaldPa / totaldPb

mhdTotalFieldOpt = struct( ...
    'enabled', true, ...
    'names', "Phi", ...
    'timeIndex', meta.nOutputTime, ...
    'modeN', meta.physicalNAll, ...
    'toroidalAngle', 0.0, ...
    'colormapIndex', 3, ...
    'bsplineOrder', 4, ...
    'titleFontSize', 14, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled      : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"。
timeIndex     : total 场量输出时间索引（1 到 nOutputTime）。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
toroidalAngle : 绘图所在环向角 phi。
colormapIndex : 色表序号，可选 1-5。
bsplineOrder  : B 样条阶数。
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.totalField = runMHDTotalFieldPlot(mhdData.total, meta, mhdFieldGeom, mhdTotalFieldOpt);
[tmpShifted, tmpAligned, hasTotalFieldResult] = extractShiftedAligned(mhdWorkspace.totalField);
if hasTotalFieldResult
    shifted = tmpShifted;
    aligned = tmpAligned;
end
clear tmpShifted tmpAligned hasFieldResult hasTotalFieldResult

%% 可视化单时间MHD极向模数

mhdTotalPoloidalModeOpt = struct( ...
    'enabled', true, ...
    'names', mhdTotalFieldOpt.names, ...
    'timeIndex', meta.nOutputTime, ...
    'modeN', mhdTotalFieldOpt.modeN, ...
    'toroidalAngle', mhdTotalFieldOpt.toroidalAngle, ...
    'mRange', [1, 20], ...
    'plotThreshold', 1e-8, ...
    'minFFT', 0, ...
    'normalize', true, ...
    'radialAxis', 'rho', ...
    'bsplineOrder', mhdTotalFieldOpt.bsplineOrder, ...
    'xLim', [0, 1], ...
    'xTicks', 0:0.2:1, ...
    'yLim', [0, 1], ...
    'yTicks', 0:0.25:1, ...
    'titleFontSize', 14, ...
    'labelFontSize', 14, ...
    'axisFontSize', 12);
%{
enabled       : 是否绘图。
names         : 可选 "Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"；对应 totalPhi 等 total 场。
timeIndex     : total 场量输出时间索引（1 到 nOutputTime）。
modeN         : 真实物理环向模数 n；[] 表示全部有效 n；示例 [6 18 36]。
toroidalAngle : 做极向截面的环向角 phi。
mRange        : 绘制的极向模数 m 范围；示例 [1 20]。
plotThreshold : 只绘制峰值超过全局峰值该比例的 m。
minFFT        : 只绘制归一化峰值超过该值的 m。
normalize     : 是否用所选 m 范围内的最大 FFT 幅值归一化。
radialAxis    : 'rho' 表示横坐标为 sqrt(s)，'psi' 表示近似使用 rho^2。
bsplineOrder  : B 样条阶数。
xLim/xTicks   : 横坐标范围和刻度；[] 表示自动。
yLim/yTicks   : 纵坐标范围和刻度；[] 表示自动。
titleFontSize : 标题字号。
labelFontSize : 坐标轴标签字号。
axisFontSize  : 坐标轴字号。
%}
mhdWorkspace.totalPoloidalMode = runMHDPoloidalModePlot(mhdData.total, meta, mhdInput, mhdFieldGeom, mhdTotalPoloidalModeOpt);

%% 局部函数

function [shifted, aligned, hasResult] = extractShiftedAligned(workspace)

    shifted = [];
    aligned = [];
    hasResult = false;

    if ~isstruct(workspace)
        return;
    end

    if isfield(workspace, 'shifted') && isfield(workspace, 'aligned')
        shifted = workspace.shifted;
        aligned = workspace.aligned;
        hasResult = true;
        return;
    end

    fields = fieldnames(workspace);
    for iField = 1:numel(fields)
        fieldName = fields{iField};
        fieldValue = workspace.(fieldName);
        if isstruct(fieldValue) && isfield(fieldValue, 'shifted') && isfield(fieldValue, 'aligned')
            shifted = fieldValue.shifted;
            aligned = fieldValue.aligned;
            hasResult = true;
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

    meta.gridNx = readIntParam(paramText, 'gridNx');
    meta.gridNy = readIntParam(paramText, 'gridNy');
    meta.gridNz = readIntParam(paramText, 'gridNz');
    meta.leftN = readIntParam(paramText, 'leftN');
    meta.rightN = readIntParam(paramText, 'rightN');
    meta.tubes = readIntParam(paramText, 'tubes');
    meta.totalSteps = readIntParam(paramText, 'totalSteps');
    meta.diagSteps = readIntParam(paramText, 'diagSteps');
    meta.outputSteps = readIntParam(paramText, 'outputSteps');
    meta.dt = readFloatParam(paramText, 'dt');
    meta.mhdPrecision = readPrecisionParam(paramText);

    assert(meta.rightN >= meta.leftN, 'rightN 必须大于或等于 leftN。');
    assert(meta.gridNx > 1 && meta.gridNy > 0 && meta.gridNz > 0, ...
        'gridNx/gridNy/gridNz 取值不合法。');
    assert(meta.tubes > 0, 'tubes 必须为正整数。');
    assert(meta.totalSteps >= 0 && meta.diagSteps > 0 && meta.outputSteps > 0, ...
        'totalSteps/diagSteps/outputSteps 取值不合法。');

    meta.switch.ifDiagAmplitude = readSwitchParam(paramText, 'ifDiagAmplitude');
    meta.switch.ifDiagFrequency = readSwitchParam(paramText, 'ifDiagFrequency');
    meta.switch.ifDiagEparallel = readSwitchParam(paramText, 'ifDiagEparallel');
    meta.switch.ifDiagZFDrive = readSwitchParam(paramText, 'ifDiagZFDrive');
    meta.switch.ifOutputPhi = readSwitchParam(paramText, 'ifOutputPhi');
    meta.switch.ifOutputA = readSwitchParam(paramText, 'ifOutputA');
    meta.switch.ifOutputdNe = readSwitchParam(paramText, 'ifOutputdNe');
    meta.switch.ifOutputdTe = readSwitchParam(paramText, 'ifOutputdTe');
    meta.switch.ifOutputdPi = readSwitchParam(paramText, 'ifOutputdPi');
    meta.switch.ifOutputdPa = readSwitchParam(paramText, 'ifOutputdPa');
    meta.switch.ifOutputdPb = readSwitchParam(paramText, 'ifOutputdPb');

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

function data = readAllMHDDiagnostics(inputDir, meta)

    data = struct( ...
        'amplitude', [], ...
        'RealMode', [], ...
        'ImagMode', [], ...
        'frequency', [], ...
        'Epara', [], ...
        'EparaES', [], ...
        'MaxwellDrive', [], ...
        'ReynoldsDrive', [], ...
        'ZonalDrive', [], ...
        'total', struct());

    diagnosticSpecs = { ...
        'amplitude', 'ifDiagAmplitude', 'amplitude.bin', 'mode'; ...
        'RealMode', 'ifDiagAmplitude', 'RealMode.bin', 'mode'; ...
        'ImagMode', 'ifDiagAmplitude', 'ImagMode.bin', 'mode'; ...
        'frequency', 'ifDiagFrequency', 'frequency.bin', 'radial'; ...
        'Epara', 'ifDiagEparallel', 'Epara.bin', 'radial'; ...
        'EparaES', 'ifDiagEparallel', 'EparaES.bin', 'radial'; ...
        'MaxwellDrive', 'ifDiagZFDrive', 'MaxwellDrive.bin', 'radial'; ...
        'ReynoldsDrive', 'ifDiagZFDrive', 'ReynoldsDrive.bin', 'radial'; ...
        'ZonalDrive', 'ifDiagZFDrive', 'ZonalDrive.bin', 'radial'};

    for iSpec = 1:size(diagnosticSpecs, 1)
        fieldName = diagnosticSpecs{iSpec, 1};
        switchName = diagnosticSpecs{iSpec, 2};
        fileName = diagnosticSpecs{iSpec, 3};
        dataKind = diagnosticSpecs{iSpec, 4};

        if ~meta.switch.(switchName)
            logSkipped(fieldName, ['开关 ' switchName ' 为 false']);
            continue;
        end

        filePath = fullfile(inputDir, fileName);
        if ~isfile(filePath)
            logSkipped(fileName, '文件不存在');
            continue;
        end

        switch dataKind
            case 'mode'
                data.(fieldName) = readModeDiagnosticAsTXN(filePath, ...
                    meta.mhdPrecision, meta.nDiagTime, meta.gridNx, meta.nMode);
            case 'radial'
                data.(fieldName) = readRadialDiagnosticAsTX(filePath, ...
                    meta.mhdPrecision, meta.nDiagTime, meta.gridNx);
            otherwise
                error('未知 MHD 诊断类型：%s。', dataKind);
        end
        logLoaded(fieldName, data.(fieldName));
    end

    totalSpecs = { ...
        'Phi', 'ifOutputPhi', 'totalPhi.bin'; ...
        'A', 'ifOutputA', 'totalA.bin'; ...
        'dNe', 'ifOutputdNe', 'totaldNe.bin'; ...
        'dTe', 'ifOutputdTe', 'totaldTe.bin'; ...
        'dPi', 'ifOutputdPi', 'totaldPi.bin'; ...
        'dPa', 'ifOutputdPa', 'totaldPa.bin'; ...
        'dPb', 'ifOutputdPb', 'totaldPb.bin'};
    for iSpec = 1:size(totalSpecs, 1)
        fieldName = totalSpecs{iSpec, 1};
        switchName = totalSpecs{iSpec, 2};
        fileName = totalSpecs{iSpec, 3};
        if ~meta.switch.(switchName)
            logSkipped(['total' fieldName], ['开关 ' switchName ' 为 false']);
            continue;
        end

        filePath = fullfile(inputDir, fileName);
        if ~isfile(filePath)
            logSkipped(fileName, '文件不存在');
            continue;
        end

        data.total.(fieldName) = readOutputAsTNYXZ(filePath, ...
            meta.mhdPrecision, meta.nOutputTime, meta.gridNy, meta.gridNx, meta.gridNz);
        logLoaded(['total' fieldName], data.total.(fieldName));
    end
end

function runInteractivePlot(interactiveMode, staticPlotFcn, interactivePlotFcn)

    switch interactiveMode
        case 0
            staticPlotFcn();
        case {1, 2}
            if isempty(interactivePlotFcn)
                staticPlotFcn();
            else
                interactivePlotFcn(interactiveMode == 2);
            end
        otherwise
            error('interactive 必须为 0、1 或 2。');
    end
end

function workspace = runMHDFieldPlot(inputDir, meta, mhdFieldGeom, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('MHD field plot', '绘图开关为 false');
        return;
    end
    if ~hasBSpline()
        logSkipped('MHD field plot', '缺少 bspline');
        return;
    end

    fieldNames = string(opt.names);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames(iField);
        fieldNameText = char(fieldName);
        fieldFile = fullfile(inputDir, [fieldNameText '.bin']);
        if ~isfile(fieldFile)
            logSkipped(fieldNameText, '文件不存在');
            continue;
        end

        fieldZYX = readMHDFieldAsZYX(inputDir, fieldName, meta.mhdPrecision, meta.gridNy, meta.gridNx, meta.gridNz);
        fieldNYXZ = mhdFieldZYXToNYXZ(fieldZYX);
        logLoaded(fieldNameText, fieldNYXZ);
        [fieldModeIndex, fieldPhysicalN] = parseModeN(opt.modeN, meta.modeIndexAll, meta.physicalNAll, 'MHD 场量');
        fieldZYX = filterMHDToroidalModes(fieldZYX, meta.modeIndexAll, fieldModeIndex);
        [shifted, aligned] = plotMHDFieldOnPoloidalPlane(fieldZYX, fieldName, mhdFieldGeom, opt, fieldPhysicalN);

        workspace.(fieldNameText) = struct( ...
            'fieldNYXZ', fieldNYXZ, ...
            'shifted', shifted, ...
            'aligned', aligned, ...
            'physicalN', fieldPhysicalN);
    end
end

function workspace = runMHDPoloidalModePlot(inputSource, meta, mhdInput, mhdFieldGeom, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('MHD poloidal mode plot', '绘图开关为 false');
        return;
    end
    if ~hasBSpline()
        logSkipped('MHD poloidal mode plot', '缺少 bspline');
        return;
    end

    useTotalField = isstruct(inputSource) && isfield(opt, 'timeIndex') && ~isempty(opt.timeIndex);
    fieldNames = string(opt.names);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames(iField);
        fieldNameText = char(fieldName);
        if useTotalField
            if ~isfield(inputSource, fieldNameText) || isempty(inputSource.(fieldNameText))
                logSkipped(['total' fieldNameText], '未读取该场量');
                continue;
            end
            totalData = inputSource.(fieldNameText);
            timeIdx = parseTimeIndex(opt.timeIndex, size(totalData, 1));
            fieldZYX = totalMHDFieldTimeSliceAsZYX(totalData, timeIdx);
            fieldNameTextForPlot = ['total' fieldNameText];
        else
            fieldFile = fullfile(inputSource, [fieldNameText '.bin']);
            if ~isfile(fieldFile)
                logSkipped(fieldNameText, '文件不存在');
                continue;
            end

            fieldZYX = readMHDFieldAsZYX(inputSource, fieldName, meta.mhdPrecision, meta.gridNy, meta.gridNx, meta.gridNz);
            fieldNameTextForPlot = fieldNameText;
        end
        [fieldModeIndex, fieldPhysicalN] = parseModeN(opt.modeN, meta.modeIndexAll, meta.physicalNAll, 'MHD 极向模数');
        fieldZYX = filterMHDToroidalModes(fieldZYX, meta.modeIndexAll, fieldModeIndex);

        [plotField, fftField, xVec, xLabelText, selectedM, yData, rawAmplitude] = ...
            calculateMHDPoloidalModeFFT(fieldZYX, mhdInput, mhdFieldGeom, opt);
        if isempty(selectedM)
            logSkipped([fieldNameTextForPlot ' poloidal mode plot'], '没有可绘制的极向模数');
            continue;
        end

        if useNormalizedPoloidalMode(opt)
            yLabelText = '$|U_m|/\max |U_m|$';
        else
            yLabelText = '$|U_m|$';
        end
        titleText = sprintf('%s poloidal FFT', fieldNameTextForPlot);
        legendText = compose('$m=%d$', selectedM);
        if useTotalField
            statusText = sprintf('%s poloidal FFT: timeIndex=%d, t_a=%.6g, n=%s, m=%s', ...
                fieldNameTextForPlot, timeIdx, meta.tOutput(timeIdx), ...
                formatNumberList(fieldPhysicalN), formatNumberList(selectedM));
        else
            statusText = sprintf('%s poloidal FFT: n=%s, m=%s', ...
                fieldNameTextForPlot, formatNumberList(fieldPhysicalN), formatNumberList(selectedM));
        end

        plotData = linePlotData(xVec, yData, xLabelText, yLabelText, titleText, legendText, statusText);
        drawLinePlot(plotData, opt);
        fprintf('[plot] %s\n', statusText);

        workspace.(fieldNameText) = struct( ...
            'poloidalPlane', plotField, ...
            'fft', fftField, ...
            'rawAmplitude', rawAmplitude, ...
            'amplitude', yData, ...
            'x', xVec, ...
            'm', selectedM, ...
            'physicalN', fieldPhysicalN);
        if useTotalField
            workspace.(fieldNameText).timeIndex = timeIdx;
        end
    end
end

function [plotField, fftField, xVec, xLabelText, selectedM, yData, rawAmplitude] = ...
    calculateMHDPoloidalModeFFT(fieldZYX, mhdInput, geom, opt)

    fieldNonShiftZYX = undoMHDFieldAlignedShift(fieldZYX, geom, opt);
    plotField = interpolateMHDFieldToPlotGrid(fieldNonShiftZYX, geom, opt);
    fftField = fft(plotField, [], 2);
    amplitude = abs(fftField);
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

    [xVec, xLabelText] = poloidalModeRadialAxis(mhdInput, geom, opt, size(plotField, 1));
end

function candidateM = poloidalModeCandidates(opt, nTheta)

    maxM = max(0, floor(nTheta / 2));
    if isfield(opt, 'mRange') && numel(opt.mRange) >= 2
        mRange = sort(round(double(opt.mRange(1:2))));
    else
        mRange = [0, maxM];
    end
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

    plotThreshold = 0;
    if isfield(opt, 'plotThreshold') && ~isempty(opt.plotThreshold)
        plotThreshold = opt.plotThreshold;
    end
    minFFT = 0;
    if isfield(opt, 'minFFT') && ~isempty(opt.minFFT)
        minFFT = opt.minFFT;
    end

    selectedMask = radialMax > globalMax * plotThreshold;
    selectedMask = selectedMask & (radialMax ./ globalMax > minFFT);
    if ~any(selectedMask)
        [~, maxIndex] = max(radialMax);
        selectedMask(maxIndex) = true;
    end

    selectedM = candidateM(selectedMask);
end

function normalizeFFT = useNormalizedPoloidalMode(opt)

    normalizeFFT = true;
    if isfield(opt, 'normalize') && ~isempty(opt.normalize)
        normalizeFFT = opt.normalize;
    end
end

function [xVec, xLabelText] = poloidalModeRadialAxis(mhdInput, geom, opt, nRadial)

    radialAxis = 'rho';
    if isfield(opt, 'radialAxis') && ~isempty(opt.radialAxis)
        radialAxis = char(opt.radialAxis);
    end

    switch lower(radialAxis)
        case {'psi', 'psip', 'psi_p'}
            xVec = radialVectorWithLength(mhdInput, geom, nRadial);
            xVec = xVec .^ 2;
            xLabelText = '$\psi_p$';
        otherwise
            xVec = radialVectorWithLength(mhdInput, geom, nRadial);
            xLabelText = '$\sqrt{s}$';
    end
end

function xVec = radialVectorWithLength(mhdInput, geom, nRadial)

    if isfield(geom, 'rhoplot') && ~isempty(geom.rhoplot)
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

function workspace = runMHDAmplitudePlot(mhdData, meta, mhdInput, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('amplitude plot', '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(mhdData, {'amplitude'})
        logSkipped('amplitude plot', '未读取 amplitude 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotAmplitude(mhdData.amplitude, meta.xGrid, meta.tDiag, meta.diagSteps, ...
        meta.modeIndexAll, meta.physicalNAll, mhdInput.B0, mhdInput.L0, mhdInput.VA0, mhdInput.TeSample, opt), ...
        @(dynamicUpdate) plotAmplitudeInteractive(mhdData.amplitude, meta.xGrid, meta.tDiag, meta.diagSteps, ...
        meta.modeIndexAll, meta.physicalNAll, mhdInput.B0, mhdInput.L0, mhdInput.VA0, mhdInput.TeSample, opt, dynamicUpdate));
end

function workspace = runMHDMultipleFrequencyPlot(mhdData, meta, mhdInput, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('mode frequency plot', '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(mhdData, {'RealMode', 'ImagMode'})
        logSkipped('mode frequency plot', '未读取 RealMode 或 ImagMode 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotMultipleFrequency(mhdData.RealMode, mhdData.ImagMode, meta.xGrid, meta.tDiag, meta.diagSteps, ...
        meta.modeIndexAll, meta.physicalNAll, mhdInput.L0, mhdInput.VA0, opt), ...
        @(dynamicUpdate) plotMultipleFrequencyInteractive(mhdData.RealMode, mhdData.ImagMode, meta.xGrid, ...
        meta.tDiag, meta.diagSteps, meta.modeIndexAll, meta.physicalNAll, mhdInput.L0, mhdInput.VA0, opt, dynamicUpdate));
end

function workspace = runMHDPhaseFrequencyPlot(mhdData, meta, mhdInput, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('phase frequency plot', '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(mhdData, {'RealMode', 'ImagMode'})
        logSkipped('phase frequency plot', '未读取 RealMode 或 ImagMode 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotPhaseFrequency(mhdData.RealMode, mhdData.ImagMode, meta.xGrid, meta.tDiag, meta.diagSteps, ...
        meta.modeIndexAll, meta.physicalNAll, mhdInput.L0, mhdInput.VA0, opt), ...
        @(dynamicUpdate) plotPhaseFrequencyInteractive(mhdData.RealMode, mhdData.ImagMode, meta.xGrid, ...
        meta.tDiag, meta.diagSteps, meta.modeIndexAll, meta.physicalNAll, mhdInput.L0, mhdInput.VA0, opt, dynamicUpdate));
end

function workspace = runMHDPhaseFrequencyMapPlot(mhdData, meta, mhdInput, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('phase frequency map', '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(mhdData, {'RealMode', 'ImagMode'})
        logSkipped('phase frequency map', '未读取 RealMode 或 ImagMode 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotPhaseFrequencyMap(mhdData.RealMode, mhdData.ImagMode, meta.xGrid, meta.tDiag, meta.diagSteps, ...
        meta.modeIndexAll, meta.physicalNAll, mhdInput.L0, mhdInput.VA0, opt), ...
        @(dynamicUpdate) plotPhaseFrequencyMapInteractive(mhdData.RealMode, mhdData.ImagMode, meta.xGrid, ...
        meta.tDiag, meta.diagSteps, meta.modeIndexAll, meta.physicalNAll, mhdInput.L0, mhdInput.VA0, opt, dynamicUpdate));
end

function workspace = runMHDSingleSignalPlot(mhdData, meta, mhdInput, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('single signal plot', '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(mhdData, {'frequency'})
        logSkipped('single signal plot', '未读取 frequency.bin 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotSingleSignal(mhdData.frequency, meta.xGrid, meta.tDiag, meta.diagSteps, mhdInput.L0, mhdInput.VA0, opt), ...
        @(dynamicUpdate) plotSingleSignalInteractive(mhdData.frequency, meta.xGrid, meta.tDiag, ...
        meta.diagSteps, mhdInput.L0, mhdInput.VA0, opt, dynamicUpdate));
end

function workspace = runMHDEparaPlot(mhdData, mhdInput, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('Epara plot', '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(mhdData, {'Epara', 'EparaES'})
        logSkipped('Epara plot', '未读取 Epara 或 EparaES 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotEpara(mhdData.Epara, mhdData.EparaES, mhdInput.rho, opt), ...
        @(dynamicUpdate) plotEparaInteractive(mhdData.Epara, mhdData.EparaES, mhdInput.rho, opt, dynamicUpdate));
end

function workspace = runMHDZonalDrivePlot(mhdData, mhdInput, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('ZF drive plot', '绘图开关为 false');
        return;
    end
    if ~hasMHDDataFields(mhdData, {'MaxwellDrive', 'ReynoldsDrive', 'ZonalDrive'})
        logSkipped('ZF drive plot', '未读取 MaxwellDrive/ReynoldsDrive/ZonalDrive 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotZonalDrive(mhdData.MaxwellDrive, mhdData.ReynoldsDrive, mhdData.ZonalDrive, mhdInput.rho, opt), ...
        @(dynamicUpdate) plotZonalDriveInteractive(mhdData.MaxwellDrive, mhdData.ReynoldsDrive, ...
        mhdData.ZonalDrive, mhdInput.rho, opt, dynamicUpdate));
end

function workspace = runMHDTotalFieldPlot(totalFields, meta, mhdFieldGeom, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('total field plot', '绘图开关为 false');
        return;
    end
    if ~hasBSpline()
        logSkipped('total field plot', '缺少 bspline');
        return;
    end

    totalFieldNames = string(opt.names);
    for iField = 1:numel(totalFieldNames)
        fieldName = totalFieldNames(iField);
        fieldNameText = char(fieldName);
        if ~isfield(totalFields, fieldNameText) || isempty(totalFields.(fieldNameText))
            logSkipped(['total' fieldNameText], '未读取该场量');
            continue;
        end

        totalData = totalFields.(fieldNameText);
        timeIdx = parseTimeIndex(opt.timeIndex, size(totalData, 1));
        fieldZYX = totalMHDFieldTimeSliceAsZYX(totalData, timeIdx);
        [fieldModeIndex, fieldPhysicalN] = parseModeN(opt.modeN, meta.modeIndexAll, meta.physicalNAll, 'total 场量');
        fieldZYX = filterMHDToroidalModes(fieldZYX, meta.modeIndexAll, fieldModeIndex);
        plotContext = sprintf('timeIndex=%d, t_a=%.6g', timeIdx, meta.tOutput(timeIdx));
        [shifted, aligned] = plotMHDFieldOnPoloidalPlane( ...
            fieldZYX, "total" + fieldName, mhdFieldGeom, opt, fieldPhysicalN, plotContext);

        workspace.(fieldNameText) = struct( ...
            'timeIndex', timeIdx, ...
            'shifted', shifted, ...
            'aligned', aligned, ...
            'physicalN', fieldPhysicalN);
    end
end

function [modeIndex, physicalN] = parseModeN(modeN, modeIndexAll, physicalNAll, contextText)

    assert(~isempty(modeIndexAll) && numel(modeIndexAll) == numel(physicalNAll), ...
        '%s modeN 列表为空或尺寸不一致。', contextText);
    if isempty(modeN)
        modeIndex = modeIndexAll;
        physicalN = physicalNAll;
        return;
    end

    if ischar(modeN) || isstring(modeN)
        key = lower(strtrim(char(modeN)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            idx = max(1, round(numel(physicalNAll) / 2));
            modeIndex = modeIndexAll(idx);
            physicalN = physicalNAll(idx);
            return;
        elseif ismember(key, {'end', 'last'})
            modeIndex = modeIndexAll(end);
            physicalN = physicalNAll(end);
            return;
        elseif ismember(key, {'first', 'begin'})
            modeIndex = modeIndexAll(1);
            physicalN = physicalNAll(1);
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

    idx = clampIndex(idx, nIndex, label);
end

function idx = parseTimeIndex(timeIndex, nTime)

    idx = parseIndex(timeIndex, nTime, 'timeIndex');
end

function idx = clampIndex(idx, maxIndex, label)

    assert(isscalar(idx) && isfinite(idx) && idx == floor(idx), ...
        '%s 下标必须是有限整数。', label);
    assert(idx >= 1 && idx <= maxIndex, ...
        '%s 下标 %d 超出有效范围 [1, %d]。', label, idx, maxIndex);
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

    ok = true;
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        ok = ok && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName));
    end
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
    token = regexp(paramText, ['const\s+(?:double|float|mhdReal|picReal)\s+' name '\s*=\s*([-+]?[0-9eE+\-\.]+)\s*;'], ...
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

function data = readModeDiagnosticAsTXN(filePath, precision, nTime, gridNx, nMode)
    raw = readBinaryVector(filePath, precision);
    expectedCount = nTime * gridNx * nMode;
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（nTime=%d, gridNx=%d, nMode=%d）。', ...
        filePath, numel(raw), expectedCount, nTime, gridNx, nMode);

    data = reshape(raw, [nMode, gridNx, nTime]);
    data = permute(data, [3, 2, 1]);
end

function data = readRadialDiagnosticAsTX(filePath, precision, nTime, gridNx)
    raw = readBinaryVector(filePath, precision);
    expectedCount = nTime * gridNx;
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（nTime=%d, gridNx=%d）。', ...
        filePath, numel(raw), expectedCount, nTime, gridNx);

    data = reshape(raw, [gridNx, nTime]);
    data = permute(data, [2, 1]);
end

function data = readOutputAsTNYXZ(filePath, precision, nTime, gridNy, gridNx, gridNz)
    raw = readBinaryVector(filePath, precision);
    expectedCount = nTime * gridNy * gridNx * gridNz;
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（nTime=%d, gridNy=%d, gridNx=%d, gridNz=%d）。', ...
        filePath, numel(raw), expectedCount, nTime, gridNy, gridNx, gridNz);

    data = reshape(raw, [gridNz, gridNx, gridNy, nTime]);
    data = permute(data, [4, 3, 2, 1]);
end

function fieldZYX = readMHDFieldAsZYX(inputDir, fieldName, precision, gridNy, gridNx, gridNz)
    fieldNameText = char(fieldName);
    filePath = fullfile(inputDir, [fieldNameText '.bin']);
    raw = readBinaryVector(filePath, precision);
    expectedCount = gridNy * gridNx * gridNz;
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（gridNy=%d, gridNx=%d, gridNz=%d）。', ...
        filePath, numel(raw), expectedCount, gridNy, gridNx, gridNz);

    fieldZYX = reshape(raw, [gridNz, gridNx, gridNy]);
end

function fieldZYX = totalMHDFieldTimeSliceAsZYX(totalDataTNYXZ, timeIndex)
    nTime = size(totalDataTNYXZ, 1);
    assert(isscalar(timeIndex) && timeIndex == floor(timeIndex) && timeIndex >= 1 && timeIndex <= nTime, ...
        'total timeIndex 必须是 [1, %d] 内的整数。', nTime);

    sliceNYXZ = reshape(totalDataTNYXZ(timeIndex, :, :, :), ...
        [size(totalDataTNYXZ, 2), size(totalDataTNYXZ, 3), size(totalDataTNYXZ, 4)]);
    fieldZYX = permute(sliceNYXZ, [3, 2, 1]);
end

function fieldNYXZ = mhdFieldZYXToNYXZ(fieldZYX)
    fieldNYXZ = permute(fieldZYX, [3, 2, 1]);
end

function fieldFilteredZYX = filterMHDToroidalModes(fieldZYX, modeIndexAll, modeIndexKeep)
    if isempty(modeIndexKeep)
        modeIndexKeep = modeIndexAll;
    end

    modeIndexKeep = unique(modeIndexKeep(:)');
    assert(all(modeIndexKeep == floor(modeIndexKeep)), 'modeIndexKeep 必须包含整数模数。');
    assert(all(ismember(modeIndexKeep, modeIndexAll)), 'modeIndexKeep 必须位于 modeIndexAll 内。');

    nZ = size(fieldZYX, 1);
    assert(all(abs(modeIndexKeep) <= floor(nZ / 2)), 'modeIndexKeep 超出 z 向 FFT 可解析范围。');

    modeMask = false(nZ, 1);
    for iMode = 1:numel(modeIndexKeep)
        modeAbs = abs(modeIndexKeep(iMode));
        modeMask(mod(modeAbs, nZ) + 1) = true;
        if modeAbs ~= 0
            modeMask(mod(-modeAbs, nZ) + 1) = true;
        end
    end

    fieldSpectrum = fft(fieldZYX, [], 1);
    fieldSpectrum(~modeMask, :, :) = 0;
    fieldFilteredZYX = real(ifft(fieldSpectrum, [], 1));
end

function geom = mhdFieldPlotGeometry(q, theta_pest, rho, qplot, rhoplot, Rplot, Zplot, yGrid, zGrid, tubes)
    [nPlotRho, nPlotTheta] = size(qplot);
    thetaPlotGrid = (0.5:nPlotTheta - 0.5) / nPlotTheta * 2 * pi - pi;

    geom.qtheta = q .* theta_pest;
    geom.rhoGrid = rho(:, 1);
    geom.Rplot = Rplot;
    geom.Zplot = Zplot;
    geom.rhoplot = rhoplot;
    geom.thplot = repmat(thetaPlotGrid, nPlotRho, 1);
    geom.qthetaplot = qplot .* geom.thplot;
    geom.yGrid = yGrid;
    geom.zGrid = zGrid;
    geom.tubes = tubes;
end

function [shifted, aligned] = plotMHDFieldOnPoloidalPlane(fieldZYX, fieldName, geom, opt, physicalN, contextText)
    if nargin < 5
        physicalN = opt.modeN;
    end
    if nargin < 6
        contextText = '';
    end

    fieldNonShiftZYX = undoMHDFieldAlignedShift(fieldZYX, geom, opt);
    shifted = mhdFieldZYXToNYXZ(fieldZYX);
    aligned = mhdFieldZYXToNYXZ(fieldNonShiftZYX);
    fieldPlotZYX = normalizeMHDFieldForPlot(fieldNonShiftZYX);
    result = interpolateMHDFieldToPlotGrid(fieldPlotZYX, geom, opt);
    drawMHDFieldSurface(result, fieldName, geom, opt);

    if strlength(string(contextText)) > 0
        contextText = [char(contextText), ', '];
    end
    fprintf('[plot] %s: %sphi=%.6g, n=[%s]\n', ...
        char(fieldName), char(contextText), opt.toroidalAngle, formatNumberList(physicalN));
end

function fieldPlotZYX = normalizeMHDFieldForPlot(fieldZYX)
    scale = max(abs(fieldZYX(:)));
    fieldPlotZYX = fieldZYX;
    if isfinite(scale) && scale > 0
        fieldPlotZYX = fieldZYX / scale;
    end
end

function fieldNonShiftZYX = undoMHDFieldAlignedShift(fieldZYX, geom, opt)
    [nZ, nX, nY] = size(fieldZYX);
    assert(numel(geom.zGrid) == nZ, 'zGrid 长度必须匹配场量 z 维度。');
    assert(isequal(size(geom.qtheta), [nX, nY]), 'qtheta 尺寸必须匹配场量 x/y 维度。');

    fieldNonShiftZYX = zeros(size(fieldZYX));
    rangeZ = {[geom.zGrid(1), geom.zGrid(1) + 2 * pi / geom.tubes]};
    derivative = uint64(0);

    for iX = 1:nX
        for iY = 1:nY
            shiftedZ = geom.zGrid + geom.qtheta(iX, iY);
            fieldNonShiftZYX(:, iX, iY) = bspline(uint64(opt.bsplineOrder), true, rangeZ, ...
                fieldZYX(:, iX, iY), shiftedZ', derivative);
        end
    end
end

function result = interpolateMHDFieldToPlotGrid(fieldNonShiftZYX, geom, opt)
    plotZ = opt.toroidalAngle - geom.qthetaplot;
    coordinate3D = [plotZ(:), geom.rhoplot(:), geom.thplot(:)];
    derivative = uint64([0, 0, 0])';
    isPeriodic = [true, false, false]';
    rangeZRhoTheta = { ...
        [geom.zGrid(1), geom.zGrid(1) + 2 * pi / geom.tubes]; ...
        [geom.rhoGrid(1), geom.rhoGrid(end)]; ...
        [geom.yGrid(1), geom.yGrid(1) + 2 * pi]};

    result = bspline(uint64(opt.bsplineOrder), isPeriodic, rangeZRhoTheta, ...
        fieldNonShiftZYX, coordinate3D, derivative);
    result = reshape(result, size(geom.Rplot));
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

    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', opt.axisFontSize, ...
        'LineWidth', 1.2, 'TickLabelInterpreter', 'latex', ...
        'XGrid', 'on', 'YGrid', 'on', 'GridLineWidth', 1.0, 'GridAlpha', 0.3, ...
        'MinorGridAlpha', 0.3, 'GridColor', 'k', 'MinorGridColor', 'k', ...
        'GridLineStyle', '-', 'MinorGridLineStyle', '-');

    xlabel(axHandle, '$R/\mathrm{m}$', 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
    ylabel(axHandle, '$Z/\mathrm{m}$', 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
    title(axHandle, char(fieldName), 'Interpreter', 'none', ...
        'FontName', 'Times New Roman', 'FontSize', opt.titleFontSize);
    box(axHandle, 'on');
end

function cmap = mhdFieldColormap(colormapIndex, valueMin, valueMax, nColor)
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
        otherwise
            error('colormapIndex 必须是 1 到 5 的整数。');
    end

    xq = linspace(valueMin, valueMax, nColor);
    if valueMin < 0 && valueMax > 0
        cmap = interp1(x, colors, xq, 'linear');
    else
        cmap = interp1(linspace(valueMin, valueMax, size(colors, 1)), colors, xq, 'linear');
    end
end

function plotAmplitude(amplitude, xGrid, tDiag, diagSteps, modeIndexAll, physicalNAll, ...
    B0, L0, VA0, TeSample, opt)
    [xData, logAmplitude, xLabelText, toroidalModeN, radialX] = ...
        calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, modeIndexAll, physicalNAll, B0, L0, VA0, TeSample, opt);

    drawLinePlot(linePlotData(xData, logAmplitude, xLabelText, '', '$e\delta\phi/T_e$', ...
        compose('$n=%d$', toroidalModeN), ''), opt);

    fprintf('[plot] amplitude: x=%.6g, radialIndex=%d, n=[%s]\n', ...
        radialX, opt.radialIndex, formatNumberList(toroidalModeN));
end

function plotAmplitudeInteractive(amplitude, xGrid, tDiag, diagSteps, modeIndexAll, physicalNAll, ...
    B0, L0, VA0, TeSample, opt, dynamicUpdate)
    nRadial = size(amplitude, 2);
    radialIndex0 = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    opt.radialIndex = radialIndex0;
    [xData0, ~, ~, ~, ~] = ...
        calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, modeIndexAll, physicalNAll, B0, L0, VA0, TeSample, opt);
    [growthStart0, growthEnd0] = amplitudeGrowthInitialRange(xData0, opt);
    controls = [ ...
        integerSliderControl('radialIndex', 'radialIndex', radialIndex0, 1, nRadial), ...
        numericSliderControl('growthStart', 'growthStart', growthStart0, xData0(1), xData0(end), numel(xData0)), ...
        numericSliderControl('growthEnd', 'growthEnd', growthEnd0, xData0(1), xData0(end), numel(xData0))];
    plotInteractiveLineWithSliders('Amplitude', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.radialIndex = values.radialIndex;
        tempOpt.growthRange = sort([values.growthStart, values.growthEnd]);
        [xData, logAmplitude, xLabelText, toroidalModeN, radialX, lnAmplitude] = ...
            calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, ...
            modeIndexAll, physicalNAll, B0, L0, VA0, TeSample, tempOpt);
        [growthRate, growthUnit, growthRange] = ...
            calculateAmplitudeGrowthRate(lnAmplitude, xData, tDiag, L0, VA0, tempOpt);
        data = linePlotData(xData, logAmplitude, xLabelText, '', '$e\delta\phi/T_e$', ...
            compose('$n=%d$', toroidalModeN), ...
            amplitudeGrowthStatus(radialX, tempOpt.radialIndex, growthRange, growthUnit, toroidalModeN, growthRate));
        data.xLines = growthRange;
    end
end

function [xData, logAmplitude, xLabelText, toroidalModeN, radialX, lnAmplitude] = ...
    calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, modeIndexAll, physicalNAll, ...
        B0, L0, VA0, TeSample, opt)
    nTime = size(amplitude, 1);
    assert(numel(tDiag) == nTime, 'tDiag 长度必须匹配 amplitude 时间维度。');

    [modeIndex, toroidalModeN] = parseModeN(opt.modeN, modeIndexAll, physicalNAll, 'Amplitude');

    nRadial = size(amplitude, 2);
    assert(numel(xGrid) == nRadial, 'xGrid 长度必须匹配 amplitude 径向维度。');
    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = xGrid(radialIdx);
    [~, modeIdx] = ismember(modeIndex, modeIndexAll);

    switch char(lower(opt.timeAxis))
        case {'ta', 'alfven', 'alfven_time'}
            xData = tDiag;
            xLabelText = '$t_a$';
        case {'ms', 'millisecond', 'milliseconds'}
            xData = tDiag * L0 / VA0 * 1000;
            xLabelText = '$t/\mathrm{ms}$';
        case {'steps', 'step'}
            xData = (0:nTime - 1) * diagSteps;
            xLabelText = '$\mathrm{steps}$';
        otherwise
            error('未知 amplitude 时间轴：%s。可选 ''ta''、''ms'' 或 ''steps''。', opt.timeAxis);
    end

    amplitudeScale = B0 * L0 * VA0 / (TeSample(1, 1) * 1000);
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
    calculateAmplitudeGrowthRate(lnAmplitude, xData, tDiag, L0, VA0, opt)
    xData = xData(:);
    tDiag = tDiag(:);
    growthRange = sort(opt.growthRange);
    growthRange = min(max(growthRange, xData(1)), xData(end));
    assert(growthRange(2) > growthRange(1), 'Amplitude growthRange 必须包含两个不同的 x 位置。');

    tRange = interp1(xData, tDiag, growthRange, 'linear');
    logStart = interp1(xData, lnAmplitude, growthRange(1), 'linear');
    logEnd = interp1(xData, lnAmplitude, growthRange(2), 'linear');

    % Growth rate uses natural log; the plotted curve still uses log10.
    growthRate = (logEnd - logStart) / (tRange(2) - tRange(1));
    growthUnitOpt = '1/wa';
    if isfield(opt, 'growthUnit') && ~isempty(opt.growthUnit)
        growthUnitOpt = opt.growthUnit;
    end
    [unitScale, growthUnit] = amplitudeGrowthUnitScale(growthUnitOpt, L0, VA0);
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

function plotMultipleFrequency(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    modeIndexAll, physicalNAll, L0, VA0, opt)
    [xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
        calculateMultipleFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        modeIndexAll, physicalNAll, L0, VA0, opt);

    drawLinePlot(linePlotData(xData, peakFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
        '$\mathrm{Short\mbox{-}Time\;Fourier\;Transform}$', compose('$n=%d$', toroidalModeN), ''), opt);

    fprintf('[plot] mode frequency: x=%.6g, radialIndex=%d, n=[%s]\n', ...
        radialX, opt.radialIndex, formatNumberList(toroidalModeN));
end

function plotMultipleFrequencyInteractive(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    modeIndexAll, physicalNAll, L0, VA0, opt, dynamicUpdate)
    nTime = size(RealMode, 1);
    if isempty(opt.nFFT)
        nFFT0 = 2 ^ nextpow2(opt.windowLength);
    else
        nFFT0 = opt.nFFT;
    end
    nFFTMax = max(nFFT0, 2 ^ nextpow2(nTime));
    controls = [ ...
        integerSliderControl('windowLength', 'windowLength', opt.windowLength, 2, nTime), ...
        integerSliderControl('windowStep', 'windowStep', opt.windowStep, 1, max(1, nTime - 1)), ...
        integerSliderControl('nFFT', 'nFFT', nFFT0, 2, nFFTMax)];

    plotInteractiveLineWithSliders('Short-Time Fourier Transform', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.windowLength = values.windowLength;
        tempOpt.windowStep = values.windowStep;
        tempOpt.nFFT = max(values.nFFT, tempOpt.windowLength);
        [xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
            calculateMultipleFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
            modeIndexAll, physicalNAll, L0, VA0, tempOpt);
        data = linePlotData(xData, peakFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
            '$\mathrm{Short\mbox{-}Time\;Fourier\;Transform}$', compose('$n=%d$', toroidalModeN), ...
            sprintf('x = %.6g, radialIndex = %d, windowLength = %d, windowStep = %d, nFFT = %d', ...
            radialX, tempOpt.radialIndex, tempOpt.windowLength, tempOpt.windowStep, tempOpt.nFFT));
    end
end

function nFFT = multipleFrequencyNFFT(opt, windowLength)

    if isempty(opt.nFFT)
        nFFT = 2 ^ nextpow2(windowLength);
    else
        nFFT = opt.nFFT;
    end
end

function [xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
    calculateMultipleFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        modeIndexAll, physicalNAll, L0, VA0, opt)
    nTime = size(RealMode, 1);
    assert(isequal(size(RealMode), size(ImagMode)), 'RealMode 和 ImagMode 尺寸必须一致。');
    assert(numel(tDiag) == nTime, 'tDiag 长度必须匹配 RealMode 时间维度。');
    assert(nTime >= 2, '频率分析至少需要两个诊断时间点。');

    [modeIndex, toroidalModeN] = parseModeN(opt.modeN, modeIndexAll, physicalNAll, 'Mode frequency');

    nRadial = size(RealMode, 2);
    assert(numel(xGrid) == nRadial, 'xGrid 长度必须匹配 RealMode 径向维度。');
    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = xGrid(radialIdx);

    windowLength = opt.windowLength;
    windowStep = opt.windowStep;
    assert(windowLength >= 2 && windowLength == floor(windowLength) && windowLength <= nTime, ...
        'windowLength 必须是 [2, %d] 内的整数。', nTime);
    assert(windowStep >= 1 && windowStep == floor(windowStep), 'windowStep 必须是正整数。');

    nFFT = multipleFrequencyNFFT(opt, windowLength);
    assert(nFFT >= windowLength && nFFT == floor(nFFT), 'nFFT 必须是不小于 windowLength 的整数。');

    dtPhysical = (tDiag(2) - tDiag(1)) * L0 / VA0;
    sampleRateHz = 1 / dtPhysical;
    frequencyHz = (-floor(nFFT / 2):ceil(nFFT / 2) - 1)' * sampleRateHz / nFFT;
    frequencyMask = true(size(frequencyHz));
    if ~isempty(opt.frequencyRangeHz)
        frequencyMask = frequencyHz >= opt.frequencyRangeHz(1) & frequencyHz <= opt.frequencyRangeHz(2);
    end
    assert(any(frequencyMask), 'frequencyRangeHz 未包含任何 FFT 频率点。');

    switch char(lower(opt.windowType))
        case 'hann'
            windowData = 0.5 - 0.5 * cos(2 * pi * (0:windowLength - 1)' / (windowLength - 1));
        case 'rect'
            windowData = ones(windowLength, 1);
        otherwise
            error('未知 windowType：%s。可选 ''hann'' 或 ''rect''。', opt.windowType);
    end

    windowStart = 1:windowStep:(nTime - windowLength + 1);
    nWindow = numel(windowStart);
    [~, modeIdx] = ismember(modeIndex, modeIndexAll);
    peakFrequencyHz = nan(nWindow, numel(modeIndex));
    centerIndex = windowStart + floor((windowLength - 1) / 2);

    for iMode = 1:numel(modeIdx)
        signal = RealMode(:, radialIdx, modeIdx(iMode)) + 1i * ImagMode(:, radialIdx, modeIdx(iMode));

        for iWindow = 1:nWindow
            idx = windowStart(iWindow):(windowStart(iWindow) + windowLength - 1);
            segment = signal(idx);
            if opt.removeMean
                segment = segment - mean(segment);
            end

            spectrum = fftshift(fft(segment(:) .* windowData, nFFT));
            powerSpectrum = abs(spectrum) .^ 2;
            powerSpectrum(~frequencyMask) = -inf;
            [~, peakIdx] = max(powerSpectrum);
            peakFrequencyHz(iWindow, iMode) = -frequencyHz(peakIdx);
        end
    end

    switch char(lower(opt.timeAxis))
        case {'ta', 'alfven', 'alfven_time'}
            xData = tDiag(centerIndex);
            xLabelText = '$t_a$';
        case {'ms', 'millisecond', 'milliseconds'}
            xData = tDiag(centerIndex) * L0 / VA0 * 1000;
            xLabelText = '$t/\mathrm{ms}$';
        case {'steps', 'step'}
            xData = (centerIndex - 1) * diagSteps;
            xLabelText = '$\mathrm{steps}$';
        otherwise
            error('未知 mode frequency 时间轴：%s。可选 ''ta''、''ms'' 或 ''steps''。', opt.timeAxis);
    end
end

function plotPhaseFrequency(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    modeIndexAll, physicalNAll, L0, VA0, opt)
    [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
        calculatePhaseFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        modeIndexAll, physicalNAll, L0, VA0, opt);

    drawLinePlot(linePlotData(xData, phaseFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
        '', compose('$n=%d$', toroidalModeN), ''), opt);

    fprintf('[plot] phase frequency: x=%.6g, radialIndex=%d, n=[%s], phaseStep=%d\n', ...
        radialX, opt.radialIndex, formatNumberList(toroidalModeN), opt.phaseStep);
end

function plotPhaseFrequencyInteractive(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    modeIndexAll, physicalNAll, L0, VA0, opt, dynamicUpdate)
    nTime = size(RealMode, 1);
    controls = [ ...
        integerSliderControl('phaseStep', 'phaseStep', opt.phaseStep, 1, nTime - 1), ...
        integerSliderControl('smoothWindow', 'smoothWindow', opt.smoothWindow, 1, max(1, min(1001, nTime - 1)))];
    plotInteractiveLineWithSliders('Phase Frequency', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.phaseStep = values.phaseStep;
        tempOpt.smoothWindow = values.smoothWindow;
        [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
            calculatePhaseFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
            modeIndexAll, physicalNAll, L0, VA0, tempOpt);
        data = linePlotData(xData, phaseFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
            '', compose('$n=%d$', toroidalModeN), ...
            sprintf('x = %.6g, radialIndex = %d, phaseStep = %d, smoothWindow = %d', ...
            radialX, tempOpt.radialIndex, tempOpt.phaseStep, tempOpt.smoothWindow));
    end
end

function [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
    calculatePhaseFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        modeIndexAll, physicalNAll, L0, VA0, opt)
    nTime = size(RealMode, 1);
    assert(isequal(size(RealMode), size(ImagMode)), 'RealMode 和 ImagMode 尺寸必须一致。');
    assert(numel(tDiag) == nTime, 'tDiag 长度必须匹配 RealMode 时间维度。');
    assert(nTime >= 2, '相位频率分析至少需要两个诊断时间点。');

    [modeIndex, toroidalModeN] = parseModeN(opt.modeN, modeIndexAll, physicalNAll, 'Phase frequency');

    nRadial = size(RealMode, 2);
    assert(numel(xGrid) == nRadial, 'xGrid 长度必须匹配 RealMode 径向维度。');
    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = xGrid(radialIdx);

    phaseStep = opt.phaseStep;
    assert(isscalar(phaseStep) && phaseStep == floor(phaseStep) && phaseStep >= 1 && phaseStep <= nTime - 1, ...
        'phaseStep 必须是 [1, %d] 内的整数。', nTime - 1);

    smoothWindow = opt.smoothWindow;
    if isempty(smoothWindow)
        smoothWindow = 1;
    end
    assert(isscalar(smoothWindow) && smoothWindow == floor(smoothWindow) && smoothWindow >= 1, ...
        'smoothWindow 必须是正整数。');

    amplitudeFloor = opt.amplitudeFloor;
    if isempty(amplitudeFloor)
        amplitudeFloor = 0;
    end
    assert(isscalar(amplitudeFloor) && amplitudeFloor >= 0, 'amplitudeFloor 必须非负。');

    dtPhysical = (tDiag(2) - tDiag(1)) * L0 / VA0;
    startIndex = (1:(nTime - phaseStep))';
    endIndex = startIndex + phaseStep;
    [~, modeIdx] = ismember(modeIndex, modeIndexAll);
    phaseFrequencyHz = nan(numel(startIndex), numel(modeIndex));

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

    switch char(lower(opt.timeAxis))
        case {'ta', 'alfven', 'alfven_time'}
            xData = 0.5 * (tDiag(startIndex) + tDiag(endIndex));
            xLabelText = '$t_a$';
        case {'ms', 'millisecond', 'milliseconds'}
            xData = 0.5 * (tDiag(startIndex) + tDiag(endIndex)) * L0 / VA0 * 1000;
            xLabelText = '$t/\mathrm{ms}$';
        case {'steps', 'step'}
            xData = ((startIndex - 1) + phaseStep / 2) * diagSteps;
            xLabelText = '$\mathrm{steps}$';
        otherwise
            error('未知 phase frequency 时间轴：%s。可选 ''ta''、''ms'' 或 ''steps''。', opt.timeAxis);
    end
end

function plotPhaseFrequencyMap(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    modeIndexAll, physicalNAll, L0, VA0, opt)

    plotData = calculatePhaseFrequencyMapData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        modeIndexAll, physicalNAll, L0, VA0, opt);
    drawFrequencyMap(plotData, opt);

    fprintf('[plot] phase frequency map: x=%.6g, radialIndex=%d, n=[%s]\n', ...
        plotData.radialX, plotData.radialIndex, formatNumberList(plotData.toroidalModeN));
end

function plotPhaseFrequencyMapInteractive(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    modeIndexAll, physicalNAll, L0, VA0, opt, dynamicUpdate)

    nTime = size(RealMode, 1);
    nRadial = size(RealMode, 2);
    radialIndex0 = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    controls = [ ...
        integerSliderControl('radialIndex', 'radialIndex', radialIndex0, 1, nRadial), ...
        integerSliderControl('phaseStep', 'phaseStep', opt.phaseStep, 1, nTime - 1), ...
        integerSliderControl('smoothWindow', 'smoothWindow', opt.smoothWindow, 1, max(1, min(1001, nTime - 1))), ...
        integerSliderControl('frequencyBins', 'frequencyBins', opt.frequencyBins, 16, 2048)];
    plotInteractiveFrequencyMapWithSliders('Phase Frequency Intensity Map', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.radialIndex = values.radialIndex;
        tempOpt.phaseStep = values.phaseStep;
        tempOpt.smoothWindow = values.smoothWindow;
        tempOpt.frequencyBins = values.frequencyBins;
        data = calculatePhaseFrequencyMapData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
            modeIndexAll, physicalNAll, L0, VA0, tempOpt);
    end
end

function plotData = calculatePhaseFrequencyMapData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    modeIndexAll, physicalNAll, L0, VA0, opt)

    [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX, intensity] = ...
        calculatePhaseFrequencyAndIntensity(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        modeIndexAll, physicalNAll, L0, VA0, opt);
    [frequencyVec, intensityMap] = depositPhaseFrequencyIntensity(phaseFrequencyHz, intensity, opt);
    [intensityMap, colorbarLabel] = scalePhaseFrequencyIntensity(intensityMap, opt);
    titleText = '$\mathrm{Phase\;frequency\;fraction}$';
    if isfield(opt, 'normalizeIntensity') && ~isempty(opt.normalizeIntensity) && ~opt.normalizeIntensity
        titleText = '$\mathrm{Phase\;frequency\;intensity}$';
    end

    plotData = frequencyMapPlotData(xData, frequencyVec, intensityMap, xLabelText, ...
        '$f/\mathrm{Hz}$', colorbarLabel, ...
        titleText, ...
        sprintf('x = %.6g, radialIndex = %d, n = [%s], phaseStep = %d, smoothWindow = %d', ...
        radialX, parseIndex(opt.radialIndex, size(RealMode, 2), 'radialIndex'), ...
        formatNumberList(toroidalModeN), opt.phaseStep, opt.smoothWindow));
    plotData.radialX = radialX;
    plotData.radialIndex = parseIndex(opt.radialIndex, size(RealMode, 2), 'radialIndex');
    plotData.toroidalModeN = toroidalModeN;
end

function [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX, intensity] = ...
    calculatePhaseFrequencyAndIntensity(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        modeIndexAll, physicalNAll, L0, VA0, opt)

    nTime = size(RealMode, 1);
    assert(isequal(size(RealMode), size(ImagMode)), 'RealMode 和 ImagMode 尺寸必须一致。');
    assert(numel(tDiag) == nTime, 'tDiag 长度必须匹配 RealMode 时间维度。');
    assert(nTime >= 2, '相位频率图至少需要两个诊断时间点。');

    [modeIndex, toroidalModeN] = parseModeN(opt.modeN, modeIndexAll, physicalNAll, 'Phase frequency map');
    nRadial = size(RealMode, 2);
    assert(numel(xGrid) == nRadial, 'xGrid 长度必须匹配 RealMode 径向维度。');
    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = xGrid(radialIdx);

    phaseStep = opt.phaseStep;
    assert(isscalar(phaseStep) && phaseStep == floor(phaseStep) && phaseStep >= 1 && phaseStep <= nTime - 1, ...
        'phaseStep 必须是 [1, %d] 内的整数。', nTime - 1);

    smoothWindow = opt.smoothWindow;
    if isempty(smoothWindow)
        smoothWindow = 1;
    end
    assert(isscalar(smoothWindow) && smoothWindow == floor(smoothWindow) && smoothWindow >= 1, ...
        'smoothWindow 必须是正整数。');

    amplitudeFloor = opt.amplitudeFloor;
    if isempty(amplitudeFloor)
        amplitudeFloor = 0;
    end
    assert(isscalar(amplitudeFloor) && amplitudeFloor >= 0, 'amplitudeFloor 必须非负。');

    intensityPower = opt.intensityPower;
    if isempty(intensityPower)
        intensityPower = 2;
    end
    assert(isscalar(intensityPower) && isfinite(intensityPower) && intensityPower > 0, ...
        'intensityPower 必须为正数。');

    dtPhysical = (tDiag(2) - tDiag(1)) * L0 / VA0;
    startIndex = (1:(nTime - phaseStep))';
    endIndex = startIndex + phaseStep;
    [~, modeIdx] = ismember(modeIndex, modeIndexAll);
    phaseFrequencyHz = nan(numel(startIndex), numel(modeIndex));
    intensity = nan(numel(startIndex), numel(modeIndex));

    for iMode = 1:numel(modeIdx)
        signal = RealMode(:, radialIdx, modeIdx(iMode)) + 1i * ImagMode(:, radialIdx, modeIdx(iMode));
        signal = signal(:);
        phase = unwrap(angle(signal));

        phaseFrequencyHz(:, iMode) = -(phase(endIndex) - phase(startIndex)) / (2 * pi * dtPhysical * phaseStep);
        intensity(:, iMode) = 0.5 * (abs(signal(startIndex)) .^ intensityPower + abs(signal(endIndex)) .^ intensityPower);

        if amplitudeFloor > 0
            lowAmplitude = abs(signal(startIndex)) < amplitudeFloor | abs(signal(endIndex)) < amplitudeFloor;
            phaseFrequencyHz(lowAmplitude, iMode) = NaN;
            intensity(lowAmplitude, iMode) = NaN;
        end

        effectiveSmoothWindow = min(smoothWindow, numel(startIndex));
        if effectiveSmoothWindow > 1
            phaseFrequencyHz(:, iMode) = smoothdata(phaseFrequencyHz(:, iMode), 'movmean', effectiveSmoothWindow);
        end
    end

    switch char(lower(opt.timeAxis))
        case {'ta', 'alfven', 'alfven_time'}
            xData = 0.5 * (tDiag(startIndex) + tDiag(endIndex));
            xLabelText = '$t_a$';
        case {'ms', 'millisecond', 'milliseconds'}
            xData = 0.5 * (tDiag(startIndex) + tDiag(endIndex)) * L0 / VA0 * 1000;
            xLabelText = '$t/\mathrm{ms}$';
        case {'steps', 'step'}
            xData = ((startIndex - 1) + phaseStep / 2) * diagSteps;
            xLabelText = '$\mathrm{steps}$';
        otherwise
            error('未知 phase frequency map 时间轴：%s。可选 ''ta''、''ms'' 或 ''steps''。', opt.timeAxis);
    end
end

function [frequencyVec, intensityMap] = depositPhaseFrequencyIntensity(phaseFrequencyHz, intensity, opt)

    nFrequency = opt.frequencyBins;
    assert(isscalar(nFrequency) && nFrequency == floor(nFrequency) && nFrequency >= 2, ...
        'frequencyBins 必须是不小于 2 的整数。');

    frequencyRange = phaseFrequencyMapRange(phaseFrequencyHz, opt);
    frequencyEdges = linspace(frequencyRange(1), frequencyRange(2), nFrequency + 1);
    frequencyVec = 0.5 * (frequencyEdges(1:end - 1) + frequencyEdges(2:end));
    intensityMap = zeros(nFrequency, size(phaseFrequencyHz, 1));

    for iMode = 1:size(phaseFrequencyHz, 2)
        freq = phaseFrequencyHz(:, iMode);
        modeIntensity = intensity(:, iMode);
        valid = isfinite(freq) & isfinite(modeIntensity) & modeIntensity > 0 & ...
            freq >= frequencyRange(1) & freq <= frequencyRange(2);
        timeIndex = find(valid);
        if isempty(timeIndex)
            continue;
        end

        binIndex = floor((freq(valid) - frequencyRange(1)) / diff(frequencyRange) * nFrequency) + 1;
        binIndex = min(max(binIndex, 1), nFrequency);
        for iPoint = 1:numel(timeIndex)
            intensityMap(binIndex(iPoint), timeIndex(iPoint)) = ...
                intensityMap(binIndex(iPoint), timeIndex(iPoint)) + modeIntensity(timeIndex(iPoint));
        end
    end

    smoothBins = opt.frequencySmoothingBins;
    if isempty(smoothBins)
        smoothBins = 1;
    end
    assert(isscalar(smoothBins) && smoothBins == floor(smoothBins) && smoothBins >= 1, ...
        'frequencySmoothingBins 必须是正整数。');
    if smoothBins > 1
        intensityMap = smoothdata(intensityMap, 1, 'movmean', smoothBins);
    end
    intensityMap(intensityMap <= 0) = NaN;
end

function frequencyRange = phaseFrequencyMapRange(phaseFrequencyHz, opt)

    if isfield(opt, 'frequencyRangeHz') && ~isempty(opt.frequencyRangeHz)
        frequencyRange = sort(reshape(double(opt.frequencyRangeHz), 1, []));
        assert(numel(frequencyRange) == 2 && all(isfinite(frequencyRange)) && frequencyRange(1) < frequencyRange(2), ...
            'frequencyRangeHz 必须为有限的 [min max]。');
        return;
    end

    finiteFrequency = phaseFrequencyHz(isfinite(phaseFrequencyHz));
    if isempty(finiteFrequency)
        frequencyRange = [-1, 1];
        return;
    end

    frequencyRange = [min(finiteFrequency), max(finiteFrequency)];
    if frequencyRange(1) == frequencyRange(2)
        padding = max(1, 0.05 * abs(frequencyRange(1)));
    else
        padding = 0.05 * diff(frequencyRange);
    end
    frequencyRange = frequencyRange + [-padding, padding];
end

function [intensityMap, colorbarLabel] = scalePhaseFrequencyIntensity(intensityMap, opt)

    normalizeIntensity = true;
    if isfield(opt, 'normalizeIntensity') && ~isempty(opt.normalizeIntensity)
        normalizeIntensity = logical(opt.normalizeIntensity);
    end

    if normalizeIntensity
        finiteIntensity = intensityMap;
        finiteIntensity(~isfinite(finiteIntensity)) = 0;
        columnSum = sum(finiteIntensity, 1);
        validColumn = columnSum > 0;
        intensityMap = nan(size(intensityMap));
        intensityMap(:, validColumn) = bsxfun(@rdivide, finiteIntensity(:, validColumn), columnSum(validColumn));
        intensityMap(intensityMap <= 0) = NaN;
    end

    if isfield(opt, 'logIntensity') && opt.logIntensity
        intensityMap = log10(intensityMap);
        if normalizeIntensity
            colorbarLabel = '$\log_{10}\mathrm{fraction}$';
        else
            colorbarLabel = '$\log_{10} I$';
        end
    elseif normalizeIntensity
        colorbarLabel = '$\mathrm{fraction}$';
    else
        colorbarLabel = '$I$';
    end
end

function plotSingleSignal(singlePointPhi, xGrid, tDiag, diagSteps, L0, VA0, opt)
    [xData, yData, xLabelText, radialX] = calculateSingleSignalData(singlePointPhi, xGrid, tDiag, diagSteps, L0, VA0, opt);
    drawLinePlot(linePlotData(xData, yData, xLabelText, '', '$\log(|\delta\phi|)$', [], ''), opt);

    fprintf('[plot] single signal: x=%.6g, radialIndex=%d\n', radialX, opt.radialIndex);
end

function plotSingleSignalInteractive(singlePointPhi, xGrid, tDiag, diagSteps, L0, VA0, opt, dynamicUpdate)
    nRadial = size(singlePointPhi, 2);
    radialIndex0 = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    controls = integerSliderControl('radialIndex', 'radialIndex', radialIndex0, 1, nRadial);
    plotInteractiveLineWithSliders('singlePointPhi', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.radialIndex = values.radialIndex;
        [xData, yData, xLabelText, radialX] = calculateSingleSignalData(singlePointPhi, xGrid, tDiag, diagSteps, L0, VA0, tempOpt);
        data = linePlotData(xData, yData, xLabelText, '', '$\log(|\delta\phi|)$', [], ...
            sprintf('x = %.6g, radialIndex = %d', radialX, tempOpt.radialIndex));
    end
end

function [xData, yData, xLabelText, radialX] = calculateSingleSignalData(singlePointPhi, xGrid, tDiag, diagSteps, L0, VA0, opt)
    nTime = size(singlePointPhi, 1);
    assert(numel(tDiag) == nTime, 'tDiag 长度必须匹配 singlePointPhi 时间维度。');

    nRadial = size(singlePointPhi, 2);
    assert(numel(xGrid) == nRadial, 'xGrid 长度必须匹配 singlePointPhi 径向维度。');

    radialIdx = parseIndex(opt.radialIndex, nRadial, 'radialIndex');
    radialX = xGrid(radialIdx);
    assert(isscalar(opt.logFloor) && opt.logFloor > 0, 'single signal logFloor 必须为正数。');

    switch char(lower(opt.timeAxis))
        case {'ta', 'alfven', 'alfven_time'}
            xData = tDiag;
            xLabelText = '$t_a$';
        case {'ms', 'millisecond', 'milliseconds'}
            xData = tDiag * L0 / VA0 * 1000;
            xLabelText = '$t/\mathrm{ms}$';
        case {'steps', 'step'}
            xData = (0:nTime - 1) * diagSteps;
            xLabelText = '$\mathrm{steps}$';
        otherwise
            error('single signal 的 timeAxis 必须为 "ta"、"ms" 或 "steps"。当前值：%s。', opt.timeAxis);
    end

    yData = log(max(abs(singlePointPhi(:, radialIdx)), opt.logFloor));
end

function plotEpara(Epara, EparaES, rho, opt)
    assert(isequal(size(Epara), size(EparaES)), 'Epara 和 EparaES 尺寸必须一致。');

    nTime = size(Epara, 1);
    nRadial = size(Epara, 2);
    timeIdx = parseTimeIndex(opt.timeIndex, nTime);
    opt.timeIndex = timeIdx;
    assert(size(rho, 1) == nRadial, 'rho(:, 1) 长度必须匹配 Epara 径向维度。');

    [xData, yData] = calculateEparaData(Epara, EparaES, rho, opt);
    drawLinePlot(linePlotData(xData, yData, '$r/a$', '', '', ...
        {'$E_{\parallel}$', '$E_{\parallel}^{\mathrm{ES}}$', '$\partial \delta A_{\parallel}/\partial t$'}, ''), opt);

    fprintf('[plot] Epara/EparaES: timeIndex=%d\n', timeIdx);
end

function plotEparaInteractive(Epara, EparaES, rho, opt, dynamicUpdate)
    nTime = size(Epara, 1);
    timeIndex0 = parseTimeIndex(opt.timeIndex, nTime);
    controls = integerSliderControl('timeIndex', 'timeIndex', timeIndex0, 1, nTime);
    plotInteractiveLineWithSliders('Epara / EparaES', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.timeIndex = values.timeIndex;
        [xData, yData] = calculateEparaData(Epara, EparaES, rho, tempOpt);
        data = linePlotData(xData, yData, '$r/a$', '', '', ...
            {'$E_{\parallel}$', '$E_{\parallel}^{\mathrm{ES}}$', '$\partial \delta A_{\parallel}/\partial t$'}, ...
            sprintf('timeIndex = %d', tempOpt.timeIndex));
    end
end

function [xData, yData] = calculateEparaData(Epara, EparaES, rho, opt)
    nTime = size(Epara, 1);
    timeIdx = parseTimeIndex(opt.timeIndex, nTime);
    xData = rho(:, 1);
    yData = [Epara(timeIdx, :)', EparaES(timeIdx, :)', (EparaES(timeIdx, :) - Epara(timeIdx, :))'];
end

function plotZonalDrive(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, opt)
    assert(isequal(size(MaxwellDrive), size(ReynoldsDrive), size(ZonalDrive)), ...
        'MaxwellDrive、ReynoldsDrive 和 ZonalDrive 尺寸必须一致。');

    nTime = size(MaxwellDrive, 1);
    nRadial = size(MaxwellDrive, 2);
    timeIdx = parseTimeIndex(opt.timeIndex, nTime);
    opt.timeIndex = timeIdx;
    assert(size(rho, 1) == nRadial, 'rho(:, 1) 长度必须匹配 zonal drive 径向维度。');

    [xData, yData] = calculateZonalDriveData(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, opt);
    drawLinePlot(linePlotData(xData, yData, '$r/a$', '', '$\mathrm{zonal\;flow\;drive}$', ...
        {'$\mathrm{MaxwellDrive}$', '$\mathrm{ReynoldsDrive}$', '$\mathrm{ZonalDrive}$'}, ''), opt);

    fprintf('[plot] ZF drive: timeIndex=%d\n', timeIdx);
end

function plotZonalDriveInteractive(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, opt, dynamicUpdate)
    nTime = size(MaxwellDrive, 1);
    timeIndex0 = parseTimeIndex(opt.timeIndex, nTime);
    controls = integerSliderControl('timeIndex', 'timeIndex', timeIndex0, 1, nTime);
    plotInteractiveLineWithSliders('MaxwellDrive / ReynoldsDrive / ZonalDrive', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.timeIndex = values.timeIndex;
        [xData, yData] = calculateZonalDriveData(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, tempOpt);
        data = linePlotData(xData, yData, '$r/a$', '', '$\mathrm{zonal\;flow\;drive}$', ...
            {'$\mathrm{MaxwellDrive}$', '$\mathrm{ReynoldsDrive}$', '$\mathrm{ZonalDrive}$'}, ...
            sprintf('timeIndex = %d', tempOpt.timeIndex));
    end
end

function [xData, yData] = calculateZonalDriveData(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, opt)
    nTime = size(MaxwellDrive, 1);
    timeIdx = parseTimeIndex(opt.timeIndex, nTime);
    xData = rho(:, 1);
    yData = [MaxwellDrive(timeIdx, :)', ReynoldsDrive(timeIdx, :)', ZonalDrive(timeIdx, :)'];
end

function applyLinePlotColorOrder(axHandle)
    set(axHandle, 'ColorOrder', preferredLineColors(), 'NextPlot', 'replacechildren');
end

function value = clampInteger(value, minValue, maxValue)
    value = round(value);
    value = min(max(value, minValue), maxValue);
end

function value = clampSliderValue(value, control)
    if isfield(control, 'allowedValues') && ~isempty(control.allowedValues)
        allowedValues = double(control.allowedValues(:));
        [~, nearestIndex] = min(abs(allowedValues - double(value)));
        value = allowedValues(nearestIndex);
    elseif control.isInteger
        value = clampInteger(value, control.min, control.max);
    else
        value = min(max(value, control.min), control.max);
    end
end

function control = integerSliderControl(fieldName, labelText, value, minValue, maxValue)
    control.field = fieldName;
    control.label = labelText;
    control.value = clampInteger(value, minValue, maxValue);
    control.min = minValue;
    control.max = maxValue;
    control.isInteger = true;
    control.nStep = [];
    control.allowedValues = [];
end

function control = numericSliderControl(fieldName, labelText, value, minValue, maxValue, nStep)
    if nargin < 6
        nStep = [];
    end

    control.field = fieldName;
    control.label = labelText;
    control.value = min(max(value, minValue), maxValue);
    control.min = minValue;
    control.max = maxValue;
    control.isInteger = false;
    control.nStep = nStep;
    control.allowedValues = [];
end

function plotData = linePlotData(xData, yData, xLabelText, yLabelText, titleText, legendText, statusText)
    plotData.xVec = xData;
    plotData.yVec = yData;
    plotData.xlabelText = xLabelText;
    plotData.ylabelText = yLabelText;
    plotData.titleText = titleText;
    plotData.legendText = legendText;
    plotData.status = statusText;
end

function plotData = frequencyMapPlotData(xData, yData, zData, xLabelText, yLabelText, colorbarLabel, titleText, statusText)
    plotData.xVec = xData;
    plotData.yVec = yData;
    plotData.Z = zData;
    plotData.xlabelText = xLabelText;
    plotData.ylabelText = yLabelText;
    plotData.colorbarLabel = colorbarLabel;
    plotData.titleText = titleText;
    plotData.status = statusText;
end

function drawLinePlot(plotData, opt)
    figure;
    axHandle = gca;
    renderLinePlot(axHandle, plotData, opt);
end

function drawFrequencyMap(plotData, opt)
    figure;
    axHandle = gca;
    renderFrequencyMap(axHandle, plotData, opt);
end

function renderFrequencyMap(axHandle, plotData, opt)
    figHandle = ancestor(axHandle, 'figure');
    delete(findall(figHandle, 'Type', 'ColorBar'));
    delete(allchild(axHandle));

    imagesc(axHandle, plotData.xVec, plotData.yVec, plotData.Z);
    set(axHandle, 'YDir', 'normal');
    colormap(axHandle, jet(256));
    cb = colorbar(axHandle);
    cb.FontName = 'Times New Roman';
    cb.FontSize = opt.axisFontSize;
    ylabel(cb, plotData.colorbarLabel, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);

    finiteZ = plotData.Z(isfinite(plotData.Z));
    if isempty(finiteZ)
        clim(axHandle, [0, 1]);
    elseif min(finiteZ) == max(finiteZ)
        centerValue = finiteZ(1);
        halfWidth = max(1, 0.05 * abs(centerValue));
        clim(axHandle, centerValue + [-halfWidth, halfWidth]);
    else
        clim(axHandle, [min(finiteZ), max(finiteZ)]);
    end

    if isfield(opt, 'frequencyRangeHz') && ~isempty(opt.frequencyRangeHz)
        ylim(axHandle, sort(opt.frequencyRangeHz));
    end

    box(axHandle, 'on');
    applyLinePlotAxesStyle(axHandle, opt.axisFontSize);
    xlabel(axHandle, plotData.xlabelText, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
    ylabel(axHandle, plotData.ylabelText, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
    title(axHandle, plotData.titleText, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.titleFontSize);
end

function renderLinePlot(axHandle, plotData, opt)
    delete(allchild(axHandle));
    applyLinePlotColorOrder(axHandle);
    plot(axHandle, plotData.xVec, plotData.yVec, 'LineWidth', 1.5);
    box(axHandle, 'on');
    applyLinePlotAxesStyle(axHandle, opt.axisFontSize);

    if isfield(opt, 'yLim') && ~isempty(opt.yLim)
        ylim(axHandle, opt.yLim);
    end
    if isfield(opt, 'yTicks') && ~isempty(opt.yTicks)
        yticks(axHandle, opt.yTicks);
    end
    if isfield(opt, 'xLim') && ~isempty(opt.xLim)
        xlim(axHandle, opt.xLim);
    end
    if isfield(opt, 'xTicks') && ~isempty(opt.xTicks)
        xticks(axHandle, opt.xTicks);
    end

    if isfield(plotData, 'xLines') && ~isempty(plotData.xLines)
        drawVerticalReferenceLines(axHandle, plotData.xLines);
    end

    xlabel(axHandle, plotData.xlabelText, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
    if ~isempty(plotData.ylabelText)
        ylabel(axHandle, plotData.ylabelText, 'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
    end
    if ~isempty(plotData.titleText)
        title(axHandle, plotData.titleText, 'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', 'FontSize', opt.titleFontSize);
    end
    if ~isempty(plotData.legendText)
        legend(axHandle, plotData.legendText, 'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', 'FontSize', opt.axisFontSize, 'Location', 'best');
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

function plotInteractiveLineWithSliders(figName, controls, dynamicUpdate, opt, computePlotData)
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
        'HorizontalAlignment', 'left', 'FontName', 'Times New Roman', 'FontSize', opt.axisFontSize);
    sliderLabels = gobjects(nControl, 1);
    sliders = gobjects(nControl, 1);

    for iControl = 1:nControl
        yPos = controlBottom + controlSpacing * (nControl - iControl);
        sliderLabels(iControl) = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.10, yPos - 0.01, 0.19, 0.035], ...
            'String', sliderLabelText(controls(iControl), controls(iControl).value), ...
            'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', opt.axisFontSize);
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
                dynamicListeners{end + 1} = addlistener(sliders(iControl), 'ContinuousValueChange', @refreshPlot); %#ok<AGROW>
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
            sliderValue = clampSliderValue(get(sliders(iSlider), 'Value'), controls(iSlider));
            values.(controls(iSlider).field) = sliderValue;
            if ~dynamicUpdate
                set(sliders(iSlider), 'Value', sliderValue);
            end
            set(sliderLabels(iSlider), 'String', sliderLabelText(controls(iSlider), sliderValue));
        end

        data = computePlotData(values);
        renderLinePlot(axHandle, data, opt);
        set(statusText, 'String', data.status);
    end
end

function plotInteractiveFrequencyMapWithSliders(figName, controls, dynamicUpdate, opt, computePlotData)
    figHandle = figure('Name', figName, 'Color', 'w', 'Position', [120, 80, 980, 700]);
    nControl = numel(controls);
    controlBottom = 0.035;
    controlSpacing = 0.045;
    statusBottom = controlBottom + controlSpacing * max(nControl, 1) + 0.010;
    axesBottom = statusBottom + 0.105;
    axHandle = axes('Parent', figHandle, 'Units', 'normalized', ...
        'Position', [0.10, axesBottom, 0.78, 0.93 - axesBottom]);
    statusText = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.10, statusBottom, 0.86, 0.035], 'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', 'FontName', 'Times New Roman', 'FontSize', opt.axisFontSize);
    sliderLabels = gobjects(nControl, 1);
    sliders = gobjects(nControl, 1);

    for iControl = 1:nControl
        yPos = controlBottom + controlSpacing * (nControl - iControl);
        sliderLabels(iControl) = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.10, yPos - 0.01, 0.19, 0.035], ...
            'String', sliderLabelText(controls(iControl), controls(iControl).value), ...
            'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', opt.axisFontSize);
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
                dynamicListeners{end + 1} = addlistener(sliders(iControl), 'ContinuousValueChange', @refreshPlot); %#ok<AGROW>
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
            sliderValue = clampSliderValue(get(sliders(iSlider), 'Value'), controls(iSlider));
            values.(controls(iSlider).field) = sliderValue;
            if ~dynamicUpdate
                set(sliders(iSlider), 'Value', sliderValue);
            end
            set(sliderLabels(iSlider), 'String', sliderLabelText(controls(iSlider), sliderValue));
        end

        data = computePlotData(values);
        renderFrequencyMap(axHandle, data, opt);
        set(statusText, 'String', data.status);
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
    if isfield(control, 'allowedValues') && ~isempty(control.allowedValues)
        step = sliderStep(numel(control.allowedValues));
    elseif control.isInteger
        step = sliderStep(control.max - control.min + 1);
    elseif isfield(control, 'nStep') && ~isempty(control.nStep) && control.nStep > 1
        step = sliderStep(control.nStep);
    else
        step = [0.005, 0.05];
    end
end

function applyLinePlotAxesStyle(axHandle, axisFontSize)
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', axisFontSize, ...
        'TickLabelInterpreter', 'latex', ...
        'XGrid', 'on', 'YGrid', 'on', 'GridLineWidth', 1.0, 'GridAlpha', 0.3, ...
        'MinorGridAlpha', 0.3, 'GridColor', 'k', 'MinorGridColor', 'k', ...
        'GridLineStyle', '-', 'MinorGridLineStyle', '-');
end

function step = sliderStep(maxValue)
    if maxValue <= 1
        step = [1, 1];
    else
        smallStep = 1 / (maxValue - 1);
        largeStep = min(1, max(1, round(maxValue / 20)) / (maxValue - 1));
        step = [smallStep, largeStep];
    end
end

function colors = preferredLineColors()
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
