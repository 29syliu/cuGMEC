%% cuGMEC MHD 诊断可视化脚本

%{

本脚本用于可视化 MHD 部分的诊断。
它会根据 cuGMEC_param.h 中的开关，读取已启用的 MHD 相关诊断和输出。

请确保 inputDir 下有：
standard2D.mat, plot2D.mat, normalization2D.mat, NTP.mat, cuGMEC_param.h, BSI

请确保 outputDir 下有当前已启用诊断和输出对应的 .bin 文件，例如：
amplitude.bin, RealMode.bin, ImagMode.bin, frequency.bin, Epara.bin,
EparaES.bin, MaxwellDrive.bin, ReynoldsDrive.bin, ZonalDrive.bin,
Phi.bin, A.bin, dNe.bin, dTe.bin, dPi.bin, dPa.bin, dPb.bin,
totalPhi.bin, totalA.bin, totaldNe.bin, totaldTe.bin,
totaldPi.bin, totaldPa.bin, totaldPb.bin

%}

%% 用户设置


inputDir = 'C:\Users\Desktop\test';
outputDir = 'C:\Users\Desktop\test';

paramFile = fullfile(inputDir, 'cuGMEC_param.h');
standardFile = fullfile(inputDir, 'standard2D.mat');
plotFile = fullfile(inputDir, 'plot2D.mat');
normalizationFile = fullfile(inputDir, 'normalization2D.mat');
NTPFile = fullfile(inputDir, 'NTP.mat');
bsiDir = fullfile(inputDir, 'BSI');
assert(isfile(paramFile), 'Missing parameter file: %s', paramFile);
assert(isfile(standardFile), 'Missing MAT file: %s', standardFile);
assert(isfile(plotFile), 'Missing MAT file: %s', plotFile);
assert(isfile(normalizationFile), 'Missing MAT file: %s', normalizationFile);
assert(isfile(NTPFile), 'Missing MAT file: %s', NTPFile);
assert(isfolder(bsiDir), 'Missing BSI directory: %s', bsiDir);
addpath(bsiDir);


%% 读取参数和开关


load(standardFile);
load(plotFile);
load(normalizationFile);
load(NTPFile);

paramText = fileread(paramFile);

gridNx = readIntParam(paramText, 'gridNx');
gridNy = readIntParam(paramText, 'gridNy');
gridNz = readIntParam(paramText, 'gridNz');
leftN = readIntParam(paramText, 'leftN');
rightN = readIntParam(paramText, 'rightN');
tubes = readIntParam(paramText, 'tubes');

totalSteps = readIntParam(paramText, 'totalSteps');
diagSteps = readIntParam(paramText, 'diagSteps');
outputSteps = readIntParam(paramText, 'outputSteps');
dt = readFloatParam(paramText, 'dt');

ifDiagAmplitude = readSwitchParam(paramText, 'ifDiagAmplitude');
ifDiagFrequency = readSwitchParam(paramText, 'ifDiagFrequency');
ifDiagEparallel = readSwitchParam(paramText, 'ifDiagEparallel');
ifDiagZFDrive = readSwitchParam(paramText, 'ifDiagZFDrive');

ifOutputPhi = readSwitchParam(paramText, 'ifOutputPhi');
ifOutputA = readSwitchParam(paramText, 'ifOutputA');
ifOutputdNe = readSwitchParam(paramText, 'ifOutputdNe');
ifOutputdTe = readSwitchParam(paramText, 'ifOutputdTe');
ifOutputdPi = readSwitchParam(paramText, 'ifOutputdPi');
ifOutputdPa = readSwitchParam(paramText, 'ifOutputdPa');
ifOutputdPb = readSwitchParam(paramText, 'ifOutputdPb');

mhdPrecision = readPrecisionParam(paramText);


%% 派生尺寸与坐标轴


nDiagTime = floor(totalSteps / diagSteps) + 1;
nOutputTime = floor(totalSteps / outputSteps) + 1;
nMode = rightN - leftN + 1;

modeIndexAll = leftN:rightN;
tDiag = (0:nDiagTime - 1) * diagSteps * dt;
tOutput = (0:nOutputTime - 1) * outputSteps * dt;

xGrid = linspace(0.0, 1.0, gridNx);
yGrid = ((0:gridNy-1) + 0.5) / gridNy * 2 * pi - pi;
zGrid = ((0:gridNz-1) + 0.5) / gridNz * 2 * pi / tubes - pi / tubes;


%% 读取并可视化单时间 MHD 场量


fieldPlot.names = "Phi";
% 可选："Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"

fieldPlot.modeN = 1;
% modeN 为 FFT 模索引；[] 表示 leftN:rightN；示例：[1 3 6] 或 0:2:6

fieldPlot.toroidalAngle = 0.0;
% 绘图所在环向角 phi

fieldPlot.colormapIndex = 1;
% 色表序号，可选 1-5

fieldPlot.bsplineOrder = 4;
fieldPlot.titleFontSize = 14;
fieldPlot.labelFontSize = 14;
fieldPlot.axisFontSize = 12;

mhdFieldGeom = mhdFieldPlotGeometry(q, theta_pest, rho, qplot, rhoplot, Rplot, Zplot, yGrid, zGrid, tubes);

fieldNames = string(fieldPlot.names);
for iField = 1:numel(fieldNames)
    fieldName = fieldNames(iField);
    fieldNameText = char(fieldName);
    fieldFile = fullfile(outputDir, [fieldNameText '.bin']);
    if ~isfile(fieldFile)
        logSkipped(fieldNameText, 'file does not exist');
        continue;
    end

    fieldZYX = readMHDFieldAsZYX(outputDir, fieldName, mhdPrecision, gridNy, gridNx, gridNz);
    fieldNYXZ = mhdFieldZYXToNYXZ(fieldZYX);
    assignin('base', fieldNameText, fieldNYXZ);
    logLoaded(fieldNameText, fieldNYXZ);
    fieldModeIndex = resolveModeIndex(fieldPlot.modeN, leftN, rightN, 'Field plot');
    fieldZYX = filterMHDToroidalModes(fieldZYX, modeIndexAll, fieldModeIndex);
    [shifted, aligned] = plotMHDFieldOnPoloidalPlane(fieldZYX, fieldName, mhdFieldGeom, fieldPlot, fieldModeIndex);
    assignin('base', 'shifted', shifted);
    assignin('base', 'aligned', aligned);
end


%% 读取 amplitude / RealMode / ImagMode


if ifDiagAmplitude
    amplitude = readModeDiagnosticAsTXN(fullfile(outputDir, 'amplitude.bin'), mhdPrecision, nDiagTime, gridNx, nMode);
    RealMode = readModeDiagnosticAsTXN(fullfile(outputDir, 'RealMode.bin'), mhdPrecision, nDiagTime, gridNx, nMode);
    ImagMode = readModeDiagnosticAsTXN(fullfile(outputDir, 'ImagMode.bin'), mhdPrecision, nDiagTime, gridNx, nMode);

    logLoaded('amplitude', amplitude);
    logLoaded('RealMode', RealMode);
    logLoaded('ImagMode', ImagMode);
end


%% 可视化 amplitude


amplitudePlot.timeAxis = 'ta';
% 'ta', 'ms', 或 'steps'

amplitudePlot.radialIndex = 90;
% 固定径向网格点

amplitudePlot.modeN = [];
% modeN 为 FFT 模索引；[] 表示 leftN:rightN；示例：[1 3 6] 或 0:2:6

amplitudePlot.growthRange = [];
% [] 表示自动范围；仅交互模式使用

amplitudePlot.growthUnit = '1/wa';
% '1/wa' 或 '1/s'

amplitudePlot.yLim = [-10 0];
% [] 表示自动范围；示例：[-8 -2]

amplitudePlot.yTicks = -8:2:0;
% [] 表示自动刻度；示例：-8:1:-2

amplitudePlot.interactive = 2;
% 0 不交互；1 静态交互；2 动态交互

amplitudePlot.titleFontSize = 16;
amplitudePlot.labelFontSize = 14;
amplitudePlot.axisFontSize = 12;

if ifDiagAmplitude
    runInteractivePlot(amplitudePlot.interactive, ...
        @() plotAmplitude(amplitude, xGrid, tDiag, diagSteps, leftN, rightN, tubes, B0, L0, VA0, TeSample, amplitudePlot), ...
        @(dynamicUpdate) plotAmplitudeInteractive(amplitude, xGrid, tDiag, diagSteps, ...
        leftN, rightN, tubes, B0, L0, VA0, TeSample, amplitudePlot, dynamicUpdate));
end


%% 可视化短时傅里叶频率


multipleFrequencyPlot.timeAxis = 'ms';
% 'ta', 'ms', 或 'steps'

multipleFrequencyPlot.radialIndex = 85;
% 固定径向网格点

multipleFrequencyPlot.modeN = [1];
% modeN 为 FFT 模索引；[] 表示 leftN:rightN；示例：[1 3 6] 或 0:2:6

multipleFrequencyPlot.windowLength = 1024;
% 短时 FFT 窗口长度，单位为诊断点数

multipleFrequencyPlot.windowStep = 64;
% 短时 FFT 窗口滑动步长，单位为诊断点数

multipleFrequencyPlot.nFFT = 15000;
% 时间 FFT 点数；[] 表示 2^nextpow2(windowLength)

multipleFrequencyPlot.frequencyRangeHz = [];
% 选峰频率范围，单位 Hz；示例：[-2e5 2e5]

multipleFrequencyPlot.windowType = 'hann';
% 'hann' 或 'rect'

multipleFrequencyPlot.removeMean = true;
% 每个窗口内去除复信号均值

multipleFrequencyPlot.interactive = 2;
% 0 不交互；1 静态交互；2 动态交互

multipleFrequencyPlot.titleFontSize = 14;
multipleFrequencyPlot.labelFontSize = 14;
multipleFrequencyPlot.axisFontSize = 12;

if ifDiagAmplitude
    runInteractivePlot(multipleFrequencyPlot.interactive, ...
        @() plotMultipleFrequency(RealMode, ImagMode, xGrid, tDiag, diagSteps, leftN, rightN, tubes, L0, VA0, multipleFrequencyPlot), ...
        @(dynamicUpdate) plotMultipleFrequencyInteractive(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        leftN, rightN, tubes, L0, VA0, multipleFrequencyPlot, dynamicUpdate));
end


%% 可视化相位频率


phaseFrequencyPlot.timeAxis = 'ms';
% 'ta', 'ms', 或 'steps'

phaseFrequencyPlot.radialIndex = 20;
% 固定径向网格点

phaseFrequencyPlot.modeN = [1];
% modeN 为 FFT 模索引；[] 表示 leftN:rightN；示例：[1 3 6] 或 0:2:6

phaseFrequencyPlot.phaseStep = 2000;
% 相位差分间隔，单位为诊断点数

phaseFrequencyPlot.smoothWindow = 1;
% 频率曲线平滑窗口，单位为差分后的点数；1 表示不平滑

phaseFrequencyPlot.amplitudeFloor = 0;
% 信号幅度低于该值时对应频率置为 NaN

phaseFrequencyPlot.yLim = [];
% [] 表示自动范围；示例：[0 2e5]

phaseFrequencyPlot.yTicks = [];
% [] 表示自动刻度；示例：0:5e4:2e5

phaseFrequencyPlot.interactive = 0;
% 0 不交互；1 静态交互；2 动态交互

phaseFrequencyPlot.titleFontSize = 14;
phaseFrequencyPlot.labelFontSize = 14;
phaseFrequencyPlot.axisFontSize = 12;

if ifDiagAmplitude
    runInteractivePlot(phaseFrequencyPlot.interactive, ...
        @() plotPhaseFrequency(RealMode, ImagMode, xGrid, tDiag, diagSteps, leftN, rightN, tubes, L0, VA0, phaseFrequencyPlot), ...
        @(dynamicUpdate) plotPhaseFrequencyInteractive(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
        leftN, rightN, tubes, L0, VA0, phaseFrequencyPlot, dynamicUpdate));
end


%% 读取单点信号


if ifDiagFrequency
    frequency = readRadialDiagnosticAsTX(fullfile(outputDir, 'frequency.bin'), mhdPrecision, nDiagTime, gridNx);

    logLoaded('frequency', frequency);
end


%% 可视化单点信号


singleFrequencyPlot.timeAxis = 'ta';
% 'ta', 'ms', 或 'steps'

singleFrequencyPlot.radialIndex = 128;
% 固定径向网格点

singleFrequencyPlot.logFloor = realmin;
% 避免 log(0)

singleFrequencyPlot.interactive = 2;
% 0 不交互；1 静态交互；2 动态交互

singleFrequencyPlot.titleFontSize = 16;
singleFrequencyPlot.labelFontSize = 14;
singleFrequencyPlot.axisFontSize = 12;

if ifDiagFrequency
    runInteractivePlot(singleFrequencyPlot.interactive, ...
        @() plotSingleFrequency(frequency, xGrid, tDiag, diagSteps, L0, VA0, singleFrequencyPlot), ...
        @(dynamicUpdate) plotSingleFrequencyInteractive(frequency, xGrid, tDiag, diagSteps, L0, VA0, singleFrequencyPlot, dynamicUpdate));
end


%% 读取 Epara / EparaES


if ifDiagEparallel
    Epara = readRadialDiagnosticAsTX(fullfile(outputDir, 'Epara.bin'), mhdPrecision, nDiagTime, gridNx);
    EparaES = readRadialDiagnosticAsTX(fullfile(outputDir, 'EparaES.bin'), mhdPrecision, nDiagTime, gridNx);

    logLoaded('Epara', Epara);
    logLoaded('EparaES', EparaES);
end


%% 可视化 Epara / EparaES


eparaPlot.timeIndex = 12000;
% 固定诊断时间点

eparaPlot.interactive = 2;
% 0 不交互；1 静态交互；2 动态交互

eparaPlot.titleFontSize = 16;
eparaPlot.labelFontSize = 14;
eparaPlot.axisFontSize = 12;

if ifDiagEparallel
    runInteractivePlot(eparaPlot.interactive, ...
        @() plotEpara(Epara, EparaES, rho, eparaPlot), ...
        @(dynamicUpdate) plotEparaInteractive(Epara, EparaES, rho, eparaPlot, dynamicUpdate));
end


%% 读取 MaxwellDrive / ReynoldsDrive / ZonalDrive


if ifDiagZFDrive
    MaxwellDrive = readRadialDiagnosticAsTX(fullfile(outputDir, 'MaxwellDrive.bin'), mhdPrecision, nDiagTime, gridNx);
    ReynoldsDrive = readRadialDiagnosticAsTX(fullfile(outputDir, 'ReynoldsDrive.bin'), mhdPrecision, nDiagTime, gridNx);
    ZonalDrive = readRadialDiagnosticAsTX(fullfile(outputDir, 'ZonalDrive.bin'), mhdPrecision, nDiagTime, gridNx);

    logLoaded('MaxwellDrive', MaxwellDrive);
    logLoaded('ReynoldsDrive', ReynoldsDrive);
    logLoaded('ZonalDrive', ZonalDrive);
end


%% 可视化 MaxwellDrive / ReynoldsDrive / ZonalDrive


zonalDrivePlot.timeIndex = 15000;
% 固定诊断时间点

zonalDrivePlot.interactive = 2;
% 0 不交互；1 静态交互；2 动态交互

zonalDrivePlot.titleFontSize = 16;
zonalDrivePlot.labelFontSize = 14;
zonalDrivePlot.axisFontSize = 12;

if ifDiagZFDrive
    runInteractivePlot(zonalDrivePlot.interactive, ...
        @() plotZonalDrive(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, zonalDrivePlot), ...
        @(dynamicUpdate) plotZonalDriveInteractive(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, zonalDrivePlot, dynamicUpdate));
end


%% 读取 totalPhi / totalA / totaldNe / totaldTe / totaldPi / totaldPa / totaldPb


if ifOutputPhi
    totalPhi = readOutputAsTNYXZ(fullfile(outputDir, 'totalPhi.bin'), mhdPrecision, nOutputTime, gridNy, gridNx, gridNz);
    assignin('base', 'totalPhi', totalPhi);
    logLoaded('totalPhi', totalPhi);
end

if ifOutputA
    totalA = readOutputAsTNYXZ(fullfile(outputDir, 'totalA.bin'), mhdPrecision, nOutputTime, gridNy, gridNx, gridNz);
    assignin('base', 'totalA', totalA);
    logLoaded('totalA', totalA);
end

if ifOutputdNe
    totaldNe = readOutputAsTNYXZ(fullfile(outputDir, 'totaldNe.bin'), mhdPrecision, nOutputTime, gridNy, gridNx, gridNz);
    assignin('base', 'totaldNe', totaldNe);
    logLoaded('totaldNe', totaldNe);
end

if ifOutputdTe
    totaldTe = readOutputAsTNYXZ(fullfile(outputDir, 'totaldTe.bin'), mhdPrecision, nOutputTime, gridNy, gridNx, gridNz);
    assignin('base', 'totaldTe', totaldTe);
    logLoaded('totaldTe', totaldTe);
end

if ifOutputdPi
    totaldPi = readOutputAsTNYXZ(fullfile(outputDir, 'totaldPi.bin'), mhdPrecision, nOutputTime, gridNy, gridNx, gridNz);
    assignin('base', 'totaldPi', totaldPi);
    logLoaded('totaldPi', totaldPi);
end

if ifOutputdPa
    totaldPa = readOutputAsTNYXZ(fullfile(outputDir, 'totaldPa.bin'), mhdPrecision, nOutputTime, gridNy, gridNx, gridNz);
    assignin('base', 'totaldPa', totaldPa);
    logLoaded('totaldPa', totaldPa);
end

if ifOutputdPb
    totaldPb = readOutputAsTNYXZ(fullfile(outputDir, 'totaldPb.bin'), mhdPrecision, nOutputTime, gridNy, gridNx, gridNz);
    assignin('base', 'totaldPb', totaldPb);
    logLoaded('totaldPb', totaldPb);
end


%% 可视化 totalPhi / totalA / totaldNe / totaldTe / totaldPi / totaldPa / totaldPb


totalFieldPlot.names = "Phi";
% 可选："Phi", "A", "dNe", "dTe", "dPi", "dPa", "dPb"

totalFieldPlot.timeIndex = nOutputTime;
% total 场量的输出时间点；1 为初始输出，nOutputTime 为最后一次输出

totalFieldPlot.modeN = leftN:rightN;
% modeN 为 FFT 模索引；[] 表示 leftN:rightN；示例：[1 3 6]

totalFieldPlot.toroidalAngle = 0.0;
% 绘图所在环向角 phi

totalFieldPlot.colormapIndex = 3;
% 色表序号，可选 1-5

totalFieldPlot.bsplineOrder = 4;
totalFieldPlot.titleFontSize = 14;
totalFieldPlot.labelFontSize = 14;
totalFieldPlot.axisFontSize = 12;

totalFields = struct();
if ifOutputPhi, totalFields.Phi = totalPhi; end
if ifOutputA, totalFields.A = totalA; end
if ifOutputdNe, totalFields.dNe = totaldNe; end
if ifOutputdTe, totalFields.dTe = totaldTe; end
if ifOutputdPi, totalFields.dPi = totaldPi; end
if ifOutputdPa, totalFields.dPa = totaldPa; end
if ifOutputdPb, totalFields.dPb = totaldPb; end

totalFieldNames = string(totalFieldPlot.names);
for iField = 1:numel(totalFieldNames)
    fieldName = totalFieldNames(iField);
    fieldNameText = char(fieldName);
    if ~isfield(totalFields, fieldNameText)
        logSkipped(['total' fieldNameText], 'field was not loaded');
        continue;
    end

    totalData = totalFields.(fieldNameText);
    fieldZYX = totalMHDFieldTimeSliceAsZYX(totalData, totalFieldPlot.timeIndex);
    fieldModeIndex = resolveModeIndex(totalFieldPlot.modeN, leftN, rightN, 'Total field plot');
    fieldZYX = filterMHDToroidalModes(fieldZYX, modeIndexAll, fieldModeIndex);
    plotContext = sprintf('timeIndex=%d, t_a=%.6g', totalFieldPlot.timeIndex, tOutput(totalFieldPlot.timeIndex));
    [shifted, aligned] = plotMHDFieldOnPoloidalPlane( ...
        fieldZYX, "total" + fieldName, mhdFieldGeom, totalFieldPlot, fieldModeIndex, plotContext);
    assignin('base', 'shifted', shifted);
    assignin('base', 'aligned', aligned);
end


%% 局部函数


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
        error('interactive must be 0, 1, or 2.');
end
end

function modeIndex = resolveModeIndex(modeN, leftN, rightN, contextText)
if isempty(modeN)
    modeIndex = leftN:rightN;
else
    modeIndex = modeN(:)';
end

assert(all(modeIndex == floor(modeIndex)), '%s modeN must contain integer values.', contextText);
assert(all(modeIndex >= leftN & modeIndex <= rightN), ...
    '%s modeN must be within [%d, %d].', contextText, leftN, rightN);
end

function logLoaded(name, data)
fprintf('[load] %s: size=[%s]\n', char(name), formatSize(size(data)));
end

function logSkipped(name, reason)
fprintf('[skip] %s: %s\n', char(name), reason);
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
assert(~isempty(token), 'Cannot find integer parameter: %s', name);
value = str2double(token{1});
end

function value = readFloatParam(paramText, name)
token = regexp(paramText, ['const\s+(?:double|float|mhdReal|picReal)\s+' name '\s*=\s*([0-9eE+\-\.]+)\s*;'], ...
    'tokens', 'once');
assert(~isempty(token), 'Cannot find float parameter: %s', name);
value = str2double(token{1});
end

function value = readSwitchParam(paramText, name)
token = regexp(paramText, ['using\s+' name '\s*=\s*(trueType|falseType)\s*;'], 'tokens', 'once');
assert(~isempty(token), 'Cannot find switch parameter: %s', name);
value = strcmp(token{1}, 'trueType');
end

function precision = readPrecisionParam(paramText)
token = regexp(paramText, 'using\s+mhdReal\s*=\s*(double|float)\s*;', 'tokens', 'once');
assert(~isempty(token), 'Cannot find mhdReal precision.');
precision = token{1};
end

function raw = readBinaryVector(filePath, precision)
assert(isfile(filePath), 'Missing file: %s', filePath);

fid = fopen(filePath, 'rb');
assert(fid >= 0, 'Cannot open file: %s', filePath);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

raw = fread(fid, inf, ['*' precision]);
end

function data = readModeDiagnosticAsTXN(filePath, precision, nTime, gridNx, nMode)
raw = readBinaryVector(filePath, precision);
expectedCount = nTime * gridNx * nMode;
assert(numel(raw) == expectedCount, 'Size mismatch in %s', filePath);

data = reshape(raw, [nMode, gridNx, nTime]);
data = permute(data, [3, 2, 1]);
end

function data = readRadialDiagnosticAsTX(filePath, precision, nTime, gridNx)
raw = readBinaryVector(filePath, precision);
expectedCount = nTime * gridNx;
assert(numel(raw) == expectedCount, 'Size mismatch in %s', filePath);

data = reshape(raw, [gridNx, nTime]);
data = permute(data, [2, 1]);
end

function data = readOutputAsTNYXZ(filePath, precision, nTime, gridNy, gridNx, gridNz)
raw = readBinaryVector(filePath, precision);
expectedCount = nTime * gridNy * gridNx * gridNz;
assert(numel(raw) == expectedCount, 'Size mismatch in %s', filePath);

data = reshape(raw, [gridNz, gridNx, gridNy, nTime]);
data = permute(data, [4, 3, 2, 1]);
end

function fieldZYX = readMHDFieldAsZYX(outputDir, fieldName, precision, gridNy, gridNx, gridNz)
fieldNameText = char(fieldName);
filePath = fullfile(outputDir, [fieldNameText '.bin']);
raw = readBinaryVector(filePath, precision);
expectedCount = gridNy * gridNx * gridNz;
assert(numel(raw) == expectedCount, 'Size mismatch in %s', filePath);

fieldZYX = reshape(raw, [gridNz, gridNx, gridNy]);
end

function fieldZYX = totalMHDFieldTimeSliceAsZYX(totalDataTNYXZ, timeIndex)
nTime = size(totalDataTNYXZ, 1);
assert(isscalar(timeIndex) && timeIndex == floor(timeIndex) && timeIndex >= 1 && timeIndex <= nTime, ...
    'total timeIndex must be an integer within [1, %d].', nTime);

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
assert(all(modeIndexKeep == floor(modeIndexKeep)), 'modeIndexKeep must contain integer mode numbers.');
assert(all(ismember(modeIndexKeep, modeIndexAll)), 'modeIndexKeep must be within modeIndexAll.');

nZ = size(fieldZYX, 1);
assert(all(abs(modeIndexKeep) <= floor(nZ / 2)), 'modeIndexKeep contains modes larger than the nz FFT range.');

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

function [shifted, aligned] = plotMHDFieldOnPoloidalPlane(fieldZYX, fieldName, geom, opt, modeIndex, contextText)
if nargin < 5
    modeIndex = opt.modeN;
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
fprintf('[plot] %s: %sphi=%.6g, modeIndex=[%s], toroidal n=[%s]\n', ...
    char(fieldName), char(contextText), opt.toroidalAngle, ...
    formatNumberList(modeIndex), formatNumberList(modeIndex * geom.tubes));
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
assert(numel(geom.zGrid) == nZ, 'zGrid length must match field nz dimension.');
assert(isequal(size(geom.qtheta), [nX, nY]), 'qtheta size must match field nx/ny dimensions.');

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
surf(axHandle, geom.Rplot, geom.Zplot, result, 'EdgeColor', 'none');
shading(axHandle, 'interp');
view(axHandle, 0, 90);
axis(axHandle, 'equal');
hold(axHandle, 'on');

finiteResult = result(isfinite(result));
if isempty(finiteResult)
    surfaceTop = 0;
    colorMin = -1;
    colorMax = 1;
else
    surfaceTop = max(finiteResult);
    colorMin = min(finiteResult);
    colorMax = max(finiteResult);
    if colorMin == colorMax
        colorMin = colorMin - 1;
        colorMax = colorMax + 1;
    end
end

boundaryZ = ones(size(geom.Rplot(1, :))) * (surfaceTop + max(1, abs(surfaceTop)) * 1e-6);
plot3(axHandle, geom.Rplot(1, :), geom.Zplot(1, :), boundaryZ, 'b', 'LineWidth', 1.5);
plot3(axHandle, geom.Rplot(end, :), geom.Zplot(end, :), boundaryZ, 'b', 'LineWidth', 1.5);

xlim(axHandle, [min(geom.Rplot(:)) - 0.5, max(geom.Rplot(:)) + 0.5]);
ylim(axHandle, [min(geom.Zplot(:)) - 0.5, max(geom.Zplot(:)) + 0.5]);

caxis(axHandle, [colorMin, colorMax]);
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
        error('colormapIndex must be an integer from 1 to 5.');
end

xq = linspace(valueMin, valueMax, nColor);
if valueMin < 0 && valueMax > 0
    cmap = interp1(x, colors, xq, 'linear');
else
    cmap = interp1(linspace(valueMin, valueMax, size(colors, 1)), colors, xq, 'linear');
end
end

function plotAmplitude(amplitude, xGrid, tDiag, diagSteps, leftN, rightN, tubes, ...
    B0, L0, VA0, TeSample, opt)
[xData, logAmplitude, xLabelText, toroidalModeN, radialX] = ...
    calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, leftN, rightN, tubes, B0, L0, VA0, TeSample, opt);

drawLinePlot(linePlotData(xData, logAmplitude, xLabelText, '', '$e\delta\phi/T_e$', ...
    compose('$n=%d$', toroidalModeN), ''), opt);

fprintf('[plot] amplitude: x=%.6g, radialIndex=%d, toroidal n=[%s]\n', ...
    radialX, opt.radialIndex, formatNumberList(toroidalModeN));
end

function plotAmplitudeInteractive(amplitude, xGrid, tDiag, diagSteps, leftN, rightN, tubes, ...
    B0, L0, VA0, TeSample, opt, dynamicUpdate)
nRadial = size(amplitude, 2);
[xData0, ~, ~, ~, ~] = ...
    calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, leftN, rightN, tubes, B0, L0, VA0, TeSample, opt);
[growthStart0, growthEnd0] = amplitudeGrowthInitialRange(xData0, opt);
controls = [ ...
    integerSliderControl('radialIndex', 'radialIndex', opt.radialIndex, 1, nRadial), ...
    numericSliderControl('growthStart', 'growthStart', growthStart0, xData0(1), xData0(end), numel(xData0)), ...
    numericSliderControl('growthEnd', 'growthEnd', growthEnd0, xData0(1), xData0(end), numel(xData0))];
plotInteractiveLineWithSliders('Amplitude', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.radialIndex = values.radialIndex;
        tempOpt.growthRange = sort([values.growthStart, values.growthEnd]);
        [xData, logAmplitude, xLabelText, toroidalModeN, radialX, lnAmplitude] = ...
            calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, ...
            leftN, rightN, tubes, B0, L0, VA0, TeSample, tempOpt);
        [growthRate, growthUnit, growthRange] = ...
            calculateAmplitudeGrowthRate(lnAmplitude, xData, tDiag, L0, VA0, tempOpt);
        assignin('base', 'amplitudeGrowthModeN', toroidalModeN);
        assignin('base', 'amplitudeGrowthRate', growthRate);
        assignin('base', 'amplitudeGrowthRange', growthRange);
        assignin('base', 'amplitudeGrowthUnit', growthUnit);
        data = linePlotData(xData, logAmplitude, xLabelText, '', '$e\delta\phi/T_e$', ...
            compose('$n=%d$', toroidalModeN), ...
            amplitudeGrowthStatus(radialX, tempOpt.radialIndex, growthRange, growthUnit, toroidalModeN, growthRate));
        data.xLines = growthRange;
    end
end

function [xData, logAmplitude, xLabelText, toroidalModeN, radialX, lnAmplitude] = ...
    calculateAmplitudeData(amplitude, xGrid, tDiag, diagSteps, leftN, rightN, tubes, ...
    B0, L0, VA0, TeSample, opt)
nTime = size(amplitude, 1);
assert(numel(tDiag) == nTime, 'tDiag length must match amplitude time dimension.');

modeIndex = resolveModeIndex(opt.modeN, leftN, rightN, 'Amplitude');

nRadial = size(amplitude, 2);
assert(numel(xGrid) == nRadial, 'xGrid length must match amplitude radial dimension.');
radialIdx = opt.radialIndex;
assert(isscalar(radialIdx) && radialIdx == floor(radialIdx) && radialIdx >= 1 && radialIdx <= nRadial, ...
    'Amplitude radial index must be an integer within [1, %d].', nRadial);
radialX = xGrid(radialIdx);
modeIdx = modeIndex - leftN + 1;
toroidalModeN = modeIndex * tubes;

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
        error('Unknown amplitude time axis: %s. Use ''ta'', ''ms'', or ''steps''.', opt.timeAxis);
end

amplitudeScale = B0 * L0 * VA0 / (TeSample(1, 1) * 1000);
amplitudeToPlot = reshape(amplitude(:, radialIdx, modeIdx), nTime, []);
scaledAmplitude = amplitudeScale * amplitudeToPlot;
logAmplitude = log10(scaledAmplitude);
lnAmplitude = log(scaledAmplitude);
end

function [growthStart, growthEnd] = amplitudeGrowthInitialRange(xData, opt)
xData = xData(:);
assert(numel(xData) >= 2, 'Amplitude growth rate needs at least two time points.');

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
assert(growthRange(2) > growthRange(1), 'Amplitude growth range must contain two different x positions.');

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
        error('Unknown amplitude growth unit: %s. Use ''1/wa'' or ''1/s''.', char(growthUnit));
end
end

function statusText = amplitudeGrowthStatus(radialX, radialIndex, growthRange, growthUnit, toroidalModeN, growthRate)
growthItems = compose('n=%d: %.6g', toroidalModeN(:), growthRate(:));
statusText = sprintf('x = %.6g, radialIndex = %d, growthRange = [%.6g, %.6g], growth(%s): %s', ...
    radialX, radialIndex, growthRange(1), growthRange(2), growthUnit, strjoin(cellstr(growthItems), ', '));
end

function plotMultipleFrequency(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt)
[xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX, sampleRateHz] = ...
    calculateMultipleFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt);

drawLinePlot(linePlotData(xData, peakFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
    '$\mathrm{Short\mbox{-}Time\;Fourier\;Transform}$', compose('$n=%d$', toroidalModeN), ''), opt);

fprintf('[plot] mode frequency: x=%.6g, radialIndex=%d, toroidal n=[%s], fs=%.6g Hz\n', ...
    radialX, opt.radialIndex, formatNumberList(toroidalModeN), sampleRateHz);
end

function plotMultipleFrequencyInteractive(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt, dynamicUpdate)
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
        [xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX, sampleRateHz] = ...
            calculateMultipleFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
            leftN, rightN, tubes, L0, VA0, tempOpt);
        data = linePlotData(xData, peakFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
            '$\mathrm{Short\mbox{-}Time\;Fourier\;Transform}$', compose('$n=%d$', toroidalModeN), ...
            sprintf('x = %.6g, radialIndex = %d, windowLength = %d, windowStep = %d, nFFT = %d, fs = %.6g Hz', ...
            radialX, tempOpt.radialIndex, tempOpt.windowLength, tempOpt.windowStep, tempOpt.nFFT, sampleRateHz));
    end
end

function [xData, peakFrequencyHz, xLabelText, toroidalModeN, radialX, sampleRateHz] = ...
    calculateMultipleFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt)
nTime = size(RealMode, 1);
assert(isequal(size(RealMode), size(ImagMode)), 'RealMode and ImagMode must have the same size.');
assert(numel(tDiag) == nTime, 'tDiag length must match RealMode time dimension.');
assert(nTime >= 2, 'At least two diagnostic time points are required for frequency analysis.');

modeIndex = resolveModeIndex(opt.modeN, leftN, rightN, 'Mode frequency');

nRadial = size(RealMode, 2);
assert(numel(xGrid) == nRadial, 'xGrid length must match RealMode radial dimension.');
radialIdx = opt.radialIndex;
assert(isscalar(radialIdx) && radialIdx == floor(radialIdx) && radialIdx >= 1 && radialIdx <= nRadial, ...
    'Mode frequency radial index must be an integer within [1, %d].', nRadial);
radialX = xGrid(radialIdx);

windowLength = opt.windowLength;
windowStep = opt.windowStep;
assert(windowLength >= 2 && windowLength == floor(windowLength) && windowLength <= nTime, ...
    'windowLength must be an integer within [2, %d].', nTime);
assert(windowStep >= 1 && windowStep == floor(windowStep), 'windowStep must be a positive integer.');

if isempty(opt.nFFT)
    nFFT = 2 ^ nextpow2(windowLength);
else
    nFFT = opt.nFFT;
end
assert(nFFT >= windowLength && nFFT == floor(nFFT), 'nFFT must be an integer no smaller than windowLength.');

dtPhysical = (tDiag(2) - tDiag(1)) * L0 / VA0;
sampleRateHz = 1 / dtPhysical;
frequencyHz = (-floor(nFFT / 2):ceil(nFFT / 2) - 1)' * sampleRateHz / nFFT;
frequencyMask = true(size(frequencyHz));
if ~isempty(opt.frequencyRangeHz)
    frequencyMask = frequencyHz >= opt.frequencyRangeHz(1) & frequencyHz <= opt.frequencyRangeHz(2);
end
assert(any(frequencyMask), 'frequencyRangeHz does not include any FFT frequency bin.');

switch char(lower(opt.windowType))
    case 'hann'
        windowData = 0.5 - 0.5 * cos(2 * pi * (0:windowLength - 1)' / (windowLength - 1));
    case 'rect'
        windowData = ones(windowLength, 1);
    otherwise
        error('Unknown windowType: %s. Use ''hann'' or ''rect''.', opt.windowType);
end

windowStart = 1:windowStep:(nTime - windowLength + 1);
nWindow = numel(windowStart);
modeIdx = modeIndex - leftN + 1;
toroidalModeN = modeIndex * tubes;
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
        error('Unknown mode frequency time axis: %s. Use ''ta'', ''ms'', or ''steps''.', opt.timeAxis);
end
end

function plotPhaseFrequency(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt)
[xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
    calculatePhaseFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt);

drawLinePlot(linePlotData(xData, phaseFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
    '', compose('$n=%d$', toroidalModeN), ''), opt);

fprintf('[plot] phase frequency: x=%.6g, radialIndex=%d, toroidal n=[%s], phaseStep=%d\n', ...
    radialX, opt.radialIndex, formatNumberList(toroidalModeN), opt.phaseStep);
end

function plotPhaseFrequencyInteractive(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt, dynamicUpdate)
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
            leftN, rightN, tubes, L0, VA0, tempOpt);
        data = linePlotData(xData, phaseFrequencyHz, xLabelText, '$f/\mathrm{Hz}$', ...
            '', compose('$n=%d$', toroidalModeN), ...
            sprintf('x = %.6g, radialIndex = %d, phaseStep = %d, smoothWindow = %d', ...
            radialX, tempOpt.radialIndex, tempOpt.phaseStep, tempOpt.smoothWindow));
    end
end

function [xData, phaseFrequencyHz, xLabelText, toroidalModeN, radialX] = ...
    calculatePhaseFrequencyData(RealMode, ImagMode, xGrid, tDiag, diagSteps, ...
    leftN, rightN, tubes, L0, VA0, opt)
nTime = size(RealMode, 1);
assert(isequal(size(RealMode), size(ImagMode)), 'RealMode and ImagMode must have the same size.');
assert(numel(tDiag) == nTime, 'tDiag length must match RealMode time dimension.');
assert(nTime >= 2, 'At least two diagnostic time points are required for phase frequency analysis.');

modeIndex = resolveModeIndex(opt.modeN, leftN, rightN, 'Phase frequency');

nRadial = size(RealMode, 2);
assert(numel(xGrid) == nRadial, 'xGrid length must match RealMode radial dimension.');
radialIdx = opt.radialIndex;
assert(isscalar(radialIdx) && radialIdx == floor(radialIdx) && radialIdx >= 1 && radialIdx <= nRadial, ...
    'Phase frequency radial index must be an integer within [1, %d].', nRadial);
radialX = xGrid(radialIdx);

phaseStep = opt.phaseStep;
assert(isscalar(phaseStep) && phaseStep == floor(phaseStep) && phaseStep >= 1 && phaseStep <= nTime - 1, ...
    'phaseStep must be an integer within [1, %d].', nTime - 1);

smoothWindow = opt.smoothWindow;
if isempty(smoothWindow)
    smoothWindow = 1;
end
assert(isscalar(smoothWindow) && smoothWindow == floor(smoothWindow) && smoothWindow >= 1, ...
    'smoothWindow must be a positive integer.');

amplitudeFloor = opt.amplitudeFloor;
if isempty(amplitudeFloor)
    amplitudeFloor = 0;
end
assert(isscalar(amplitudeFloor) && amplitudeFloor >= 0, 'amplitudeFloor must be non-negative.');

dtPhysical = (tDiag(2) - tDiag(1)) * L0 / VA0;
startIndex = (1:(nTime - phaseStep))';
endIndex = startIndex + phaseStep;
modeIdx = modeIndex - leftN + 1;
toroidalModeN = modeIndex * tubes;
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
        error('Unknown phase frequency time axis: %s. Use ''ta'', ''ms'', or ''steps''.', opt.timeAxis);
end
end

function plotSingleFrequency(frequency, xGrid, tDiag, diagSteps, L0, VA0, opt)
[xData, yData, xLabelText, radialX] = calculateSingleFrequencyData(frequency, xGrid, tDiag, diagSteps, L0, VA0, opt);
drawLinePlot(linePlotData(xData, yData, xLabelText, '', '$\log(|\delta\phi|)$', [], ''), opt);

fprintf('[plot] single frequency: x=%.6g, radialIndex=%d\n', radialX, opt.radialIndex);
end

function plotSingleFrequencyInteractive(frequency, xGrid, tDiag, diagSteps, L0, VA0, opt, dynamicUpdate)
nRadial = size(frequency, 2);
controls = integerSliderControl('radialIndex', 'radialIndex', opt.radialIndex, 1, nRadial);
plotInteractiveLineWithSliders('log(delta phi)', controls, dynamicUpdate, opt, @computePlotData);

    function data = computePlotData(values)
        tempOpt = opt;
        tempOpt.radialIndex = values.radialIndex;
        [xData, yData, xLabelText, radialX] = calculateSingleFrequencyData(frequency, xGrid, tDiag, diagSteps, L0, VA0, tempOpt);
        data = linePlotData(xData, yData, xLabelText, '', '$\log(|\delta\phi|)$', [], ...
            sprintf('x = %.6g, radialIndex = %d', radialX, tempOpt.radialIndex));
    end
end

function [xData, yData, xLabelText, radialX] = calculateSingleFrequencyData(frequency, xGrid, tDiag, diagSteps, L0, VA0, opt)
nTime = size(frequency, 1);
assert(numel(tDiag) == nTime, 'tDiag length must match frequency time dimension.');

nRadial = size(frequency, 2);
assert(numel(xGrid) == nRadial, 'xGrid length must match frequency radial dimension.');

radialIdx = opt.radialIndex;
assert(isscalar(radialIdx) && radialIdx == floor(radialIdx) && radialIdx >= 1 && radialIdx <= nRadial, ...
    'Frequency radial index must be an integer within [1, %d].', nRadial);
radialX = xGrid(radialIdx);
assert(isscalar(opt.logFloor) && opt.logFloor > 0, 'Frequency logFloor must be positive.');

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
        error('Unknown single frequency time axis: %s. Use ''ta'', ''ms'', or ''steps''.', opt.timeAxis);
end

yData = log(max(abs(frequency(:, radialIdx)), opt.logFloor));
end

function plotEpara(Epara, EparaES, rho, opt)
assert(isequal(size(Epara), size(EparaES)), 'Epara and EparaES must have the same size.');

nTime = size(Epara, 1);
nRadial = size(Epara, 2);
timeIdx = opt.timeIndex;
assert(isscalar(timeIdx) && timeIdx == floor(timeIdx) && timeIdx >= 1 && timeIdx <= nTime, ...
    'Epara time index must be an integer within [1, %d].', nTime);
assert(size(rho, 1) == nRadial, 'rho(:, 1) length must match Epara radial dimension.');

[xData, yData] = calculateEparaData(Epara, EparaES, rho, opt);
drawLinePlot(linePlotData(xData, yData, '$r/a$', '', '', ...
    {'$E_{\parallel}$', '$E_{\parallel}^{\mathrm{ES}}$', '$\partial \delta A_{\parallel}/\partial t$'}, ''), opt);

fprintf('[plot] Epara/EparaES: timeIndex=%d\n', timeIdx);
end

function plotEparaInteractive(Epara, EparaES, rho, opt, dynamicUpdate)
nTime = size(Epara, 1);
controls = integerSliderControl('timeIndex', 'timeIndex', opt.timeIndex, 1, nTime);
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
timeIdx = opt.timeIndex;
assert(timeIdx >= 1 && timeIdx <= nTime, 'Epara time index must be within [1, %d].', nTime);
xData = rho(:, 1);
yData = [Epara(timeIdx, :)', EparaES(timeIdx, :)', (EparaES(timeIdx, :) - Epara(timeIdx, :))'];
end

function plotZonalDrive(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, opt)
assert(isequal(size(MaxwellDrive), size(ReynoldsDrive), size(ZonalDrive)), ...
    'MaxwellDrive, ReynoldsDrive, and ZonalDrive must have the same size.');

nTime = size(MaxwellDrive, 1);
nRadial = size(MaxwellDrive, 2);
timeIdx = opt.timeIndex;
assert(isscalar(timeIdx) && timeIdx == floor(timeIdx) && timeIdx >= 1 && timeIdx <= nTime, ...
    'Zonal drive time index must be an integer within [1, %d].', nTime);
assert(size(rho, 1) == nRadial, 'rho(:, 1) length must match zonal drive radial dimension.');

[xData, yData] = calculateZonalDriveData(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, opt);
drawLinePlot(linePlotData(xData, yData, '$r/a$', '', '$\mathrm{zonal\;flow\;drive}$', ...
    {'$\mathrm{MaxwellDrive}$', '$\mathrm{ReynoldsDrive}$', '$\mathrm{ZonalDrive}$'}, ''), opt);

fprintf('[plot] ZF drive: timeIndex=%d\n', timeIdx);
end

function plotZonalDriveInteractive(MaxwellDrive, ReynoldsDrive, ZonalDrive, rho, opt, dynamicUpdate)
nTime = size(MaxwellDrive, 1);
controls = integerSliderControl('timeIndex', 'timeIndex', opt.timeIndex, 1, nTime);
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
timeIdx = opt.timeIndex;
assert(timeIdx >= 1 && timeIdx <= nTime, 'Zonal drive time index must be within [1, %d].', nTime);
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
if control.isInteger
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
end

function plotData = linePlotData(xData, yData, xLabelText, yLabelText, titleText, legendText, statusText)
plotData.x = xData;
plotData.y = yData;
plotData.xLabel = xLabelText;
plotData.yLabel = yLabelText;
plotData.title = titleText;
plotData.legend = legendText;
plotData.status = statusText;
end

function drawLinePlot(plotData, opt)
figure;
axHandle = gca;
renderLinePlot(axHandle, plotData, opt);
end

function renderLinePlot(axHandle, plotData, opt)
delete(allchild(axHandle));
applyLinePlotColorOrder(axHandle);
plot(axHandle, plotData.x, plotData.y, 'LineWidth', 1.5);
box(axHandle, 'on');
applyLinePlotAxesStyle(axHandle, opt.axisFontSize);

if isfield(opt, 'yLim') && ~isempty(opt.yLim)
    ylim(axHandle, opt.yLim);
end
if isfield(opt, 'yTicks') && ~isempty(opt.yTicks)
    yticks(axHandle, opt.yTicks);
end

if isfield(plotData, 'xLines') && ~isempty(plotData.xLines)
    drawVerticalReferenceLines(axHandle, plotData.xLines);
end

xlabel(axHandle, plotData.xLabel, 'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
if ~isempty(plotData.yLabel)
    ylabel(axHandle, plotData.yLabel, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.labelFontSize);
end
if ~isempty(plotData.title)
    title(axHandle, plotData.title, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', 'FontSize', opt.titleFontSize);
end
if ~isempty(plotData.legend)
    legend(axHandle, plotData.legend, 'Interpreter', 'latex', ...
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
            'Dynamic slider update unavailable; using callback updates.');
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

function text = sliderLabelText(control, value)
if control.isInteger
    text = sprintf('%s = %d', control.label, value);
else
    text = sprintf('%s = %.6g', control.label, value);
end
end

function step = sliderStepForControl(control)
if control.isInteger
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
