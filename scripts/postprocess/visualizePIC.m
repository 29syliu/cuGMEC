%% cuGMEC PIC 诊断统一可视化脚本

%{
功能：
统一可视化 PIC phase / pitch / orbit 诊断。

输入目录：
将相关输入、诊断和输出文件放在同一个 inputDir 中。

必须包含：
cuGMEC_param.h, normalization2D.mat

按开关读取：
ifIon, ifAlpha, ifBeam,
ifDiagDiffusivity,
ifOutputPhaseSpaceOrbit, ifOutputPhaseSpaceJacobian, ifOutputPhaseSpaceF0,
ifOutputPhaseSpaceDeltaF, ifOutputPhaseSpacePower,
ifOutputPitchSpaceJacobian, ifOutputPitchSpaceF0,
ifOutputPitchSpaceDeltaF, ifOutputPitchSpacePower

可能包含：
Ion/Alpha/BeamPhaseSpaceOrbit.bin,
Ion/Alpha/BeamPhaseSpaceJacobian.bin,
Ion/Alpha/BeamPhaseSpaceF0.bin,
Ion/Alpha/BeamPhaseDeltaF.bin,
Ion/Alpha/BeamPhasePower.bin,
Ion/Alpha/BeamPitchSpaceJacobian.bin,
Ion/Alpha/BeamPitchSpaceF0.bin,
Ion/Alpha/BeamPitchDeltaF.bin,
Ion/Alpha/BeamPitchPower.bin,
Ion/Alpha/BeamDiffusivity.bin
%}

%% 用户设置

clear; close all;

inputDir = 'C:\Users\Desktop\test';

paramFile = fullfile(inputDir, 'cuGMEC_param.h');
normalizationFile = fullfile(inputDir, 'normalization2D.mat');

speciesList = {'Ion', 'Alpha', 'Beam'};

%% 读取所有输入

assert(isfile(paramFile), '缺少参数文件：%s', paramFile);
assert(isfile(normalizationFile), '缺少 MAT 文件：%s', normalizationFile);
assert(isfolder(inputDir), '缺少输入目录：%s', inputDir);

paramText = fileread(paramFile);
normData = load(normalizationFile);
meta = readPICMetadata(paramText, normData);
speciesList = enabledPICSpecies(speciesList, meta);
searchDirs = {inputDir};

picPhaseData = initializePhaseSpaceData(speciesList, normData, meta);
picPitchData = initializePitchSpaceData(speciesList, paramText, meta);
picDiffusivityData = initializeDiffusivityData(speciesList, meta);

orbitRaw = readAllOrbitRaw(speciesList, searchDirs, meta);
picPhaseData = readAllPhaseDiagnostics(picPhaseData, speciesList, searchDirs, meta);
picPitchData = readAllPitchDiagnostics(picPitchData, speciesList, searchDirs, meta);
picDiffusivityData = readAllDiffusivityDiagnostics(picDiffusivityData, speciesList, searchDirs, meta);

picWorkspace = struct();
picWorkspace.meta = meta;
picWorkspace.phase = picPhaseData;
picWorkspace.pitch = picPitchData;
picWorkspace.diffusivity = picDiffusivityData;

%% 处理 orbit.bin

picOrbitProcessOpt = struct( ...
    'enabled', true, ...
    'plotConservationDiagnostics', false);
%{
enabled                    : 是否处理 orbit 文件。
plotConservationDiagnostics: 是否绘制守恒量误差。
%}

phaseSpaceOrbit = struct();
PhaseSpaceOrbitSummary = struct();

if picOrbitProcessOpt.enabled
    [phaseSpaceOrbit, PhaseSpaceOrbitSummary] = processAllOrbitRaw( ...
        orbitRaw, speciesList, picPhaseData, picOrbitProcessOpt);
    printOrbitSummaryIndex(PhaseSpaceOrbitSummary);
else
    logSkipped('PhaseSpaceOrbit processing', '处理开关为 false');
end

picWorkspace.orbit = phaseSpaceOrbit;
picWorkspace.orbitSummary = PhaseSpaceOrbitSummary;

%% 可视化轨道频率

picOrbitFrequencyOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'fixedCoordinate', 'E', ...
    'slice', max(1, round(meta.gridE / 2)), ...
    'unit', 'w', ...
    'branch', 'para', ...
    'colormapIndex', 1, ...
    'contourCount', 20, ...
    'interactive', 2);
%{
enabled         : 是否绘图。
species         : 'Ion', 'Alpha', 'Beam'。
fixedCoordinate : 'E', 'Pphi', 'Lambda'。
slice           : 固定坐标方向的整数索引；默认固定 E。
unit            : 'Hz' 或 'w'。
branch          : 'para' 或 'anti'。
colormapIndex   : 非负色表，可选 1-2；有正负固定红蓝。
contourCount    : 等值线数量；0 表示不绘制。
interactive     : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
picWorkspace.orbitFrequency = runPICOrbitFrequencyPlot(phaseSpaceOrbit, meta, picOrbitFrequencyOpt);

%% 可视化共振线

picResonanceLineOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'fixedCoordinate', 'E', ...
    'slice', max(1, round(meta.gridE / 2)), ...
    'branch', 'para', ...
    'frequencyHz', 65e3, ...
    'frequencyHzRange', [0, 150e3], ...
    'toroidalMode', 30, ...
    'toroidalModeRange', [1, 80], ...
    'poloidalMode', 33, ...
    'poloidalModeRange', [1, 80], ...
    'harmonic', 0, ...
    'harmonicRange', [-20, 20], ...
    'colormapIndex', 1, ...
    'interactive', 2);
%{
enabled         : 是否绘图。
species         : 'Ion', 'Alpha', 'Beam'。
fixedCoordinate : 'E', 'Pphi', 'Lambda'。
slice           : 固定坐标方向的整数索引；默认固定 E。
branch          : 'para', 'anti' 或 'trapped'。
frequencyHz     : 有符号模频率，单位 Hz。
frequencyHzRange: 频率滑块范围，例如 [-150e3 150e3]。
toroidalMode    : 真实物理环向模数 n。
toroidalModeRange: n 滑块范围。
poloidalMode    : 正极向模数 m；trapped 分支不使用。
poloidalModeRange: m 滑块范围。
harmonic         : 轨道谐波 l。
harmonicRange    : l 滑块范围。
colormapIndex    : 共振残差固定红蓝。
interactive      : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
picWorkspace.resonanceLine = runPICResonanceLinePlot(phaseSpaceOrbit, meta, picResonanceLineOpt);

%% 可视化 phase: J / F0 / DF

picPhaseQuantityOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'quantity', 'J', ...
    'fixedCoordinate', 'Lambda', ...
    'slice', max(1, round(meta.gridLambda / 2)), ...
    'timeIndex', meta.nOutputTime, ...
    'colormapIndex', 1, ...
    'contourCount', 0, ...
    'resonance', picResonanceOptions(false, 'para', 65e3, [0, 150e3], 30, [1, 80], 33, [1, 80], 0, [-20, 20]), ...
    'interactive', 2);
%{
enabled         : 是否绘图。
species         : 'Ion', 'Alpha', 'Beam'。
quantity        : 'J', 'F0', 'DF', 'f0', 'df', 'df/f0'。
fixedCoordinate : 'E', 'Pphi', 'Lambda'。
slice           : 固定坐标方向的整数索引；默认固定 Lambda。
timeIndex       : 输出时间索引（1 到 nOutputTime）；仅用于 DF / df / df/f0。
colormapIndex   : 非负色表，可选 1-2；有正负固定红蓝。
contourCount    : 等值线数量；0 表示不绘制。
resonance       : 共振线选项，由 picResonanceOptions 构造。
resonance.enabled: 是否叠加共振线。
resonance.branch : 'para', 'anti' 或 'trapped'。
resonance.frequencyHz     : 共振线有符号模频率，单位 Hz。
resonance.frequencyHzRange: 共振线频率滑块范围。
resonance.toroidalMode    : 共振线真实物理环向模数 n。
resonance.toroidalModeRange: 共振线 n 滑块范围。
resonance.poloidalMode    : 共振线正极向模数 m；trapped 分支不使用。
resonance.poloidalModeRange: 共振线 m 滑块范围。
resonance.harmonic         : 共振线轨道谐波 l；叠加到 phase 图时也可设为 [lMin lMax]。
resonance.harmonicRange    : l min / l max 滑块允许范围。
interactive      : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
% For phase quantity interactive=1/2, resonance.branch may also be a list,
% for example {'para','anti'}, {'para','trapped'}, {'anti','trapped'},
% {'para','anti','trapped'}, ["para","anti"], or 'all'. The remaining
% resonance parameters and sliders are shared by all selected branches.
picWorkspace.phaseQuantity = runPICPhaseQuantityPlot(picPhaseData, phaseSpaceOrbit, meta, picPhaseQuantityOpt);

%% 可视化 phase power

picPhasePowerOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'fixedCoordinate', 'Lambda', ...
    'slice', max(1, round(meta.gridLambda / 2)), ...
    'modeN', meta.physicalNAll(max(1, round(numel(meta.physicalNAll) / 2))), ...
    'timeIndex', meta.nOutputTime, ...
    'colormapIndex', 1, ...
    'contourCount', 0, ...
    'resonance', picResonanceOptions(false, 'para', 65e3, [0, 150e3], 30, [1, 80], 33, [1, 80], 0, [-20, 20]), ...
    'interactive', 2);
%{
enabled         : 是否绘图。
species         : 'Ion', 'Alpha', 'Beam'。
fixedCoordinate : 'E', 'Pphi', 'Lambda'。
slice           : 固定坐标方向的整数索引；默认固定 Lambda。
modeN           : 真实物理环向模数 n，不是 mode 数组下标。
timeIndex       : 输出时间索引（1 到 nOutputTime）。
colormapIndex   : 非负色表，可选 1-2；有正负固定红蓝。
contourCount    : 等值线数量；0 表示不绘制。
resonance       : 共振线选项，由 picResonanceOptions 构造。
resonance.enabled: 是否叠加共振线。
resonance.branch : 'para', 'anti' 或 'trapped'。
resonance.frequencyHz     : 共振线有符号模频率，单位 Hz。
resonance.frequencyHzRange: 共振线频率滑块范围。
resonance.toroidalMode    : 共振线真实物理环向模数 n。
resonance.toroidalModeRange: 共振线 n 滑块范围。
resonance.poloidalMode    : 共振线正极向模数 m；trapped 分支不使用。
resonance.poloidalModeRange: 共振线 m 滑块范围。
resonance.harmonic         : 共振线轨道谐波 l；叠加到 phase 图时也可设为 [lMin lMax]。
resonance.harmonicRange    : l min / l max 滑块允许范围。
interactive      : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
picWorkspace.phasePower = runPICPhasePowerPlot(picPhaseData, phaseSpaceOrbit, meta, picPhasePowerOpt);

%% 可视化 pitch: J / F0 / DF

picPitchQuantityOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'quantity', 'J', ...
    'timeIndex', meta.nOutputTime, ...
    'colormapIndex', 1, ...
    'contourCount', 0, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
species       : 'Ion', 'Alpha', 'Beam'。
quantity      : 'J', 'F0', 'DF', 'f0', 'df', 'df/f0'。
timeIndex     : 输出时间索引（1 到 nOutputTime）；仅用于 DF / df / df/f0。
colormapIndex : 非负色表，可选 1-2；有正负固定红蓝。
contourCount  : 等值线数量；0 表示不绘制。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
picWorkspace.pitchQuantity = runPICPitchQuantityPlot(picPitchData, picPitchQuantityOpt);

%% 可视化 pitch power

picPitchPowerOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'modeN', meta.physicalNAll(max(1, round(numel(meta.physicalNAll) / 2))), ...
    'timeIndex', meta.nOutputTime, ...
    'colormapIndex', 1, ...
    'contourCount', 0, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
species       : 'Ion', 'Alpha', 'Beam'。
modeN         : 真实物理环向模数 n，不是 mode 数组下标。
timeIndex     : 输出时间索引（1 到 nOutputTime）。
colormapIndex : 非负色表，可选 1-2；有正负固定红蓝。
contourCount  : 等值线数量；0 表示不绘制。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
picWorkspace.pitchPower = runPICPitchPowerPlot(picPitchData, picPitchPowerOpt);

%% 可视化径向扩散系数

picDiffusivityOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'plotType', 1, ...
    'nRange', [min(meta.physicalNAll), max(meta.physicalNAll)], ...
    'timeIndex', max(1, round(meta.nDiagTime / 2)), ...
    'radialIndex', max(1, round(meta.gridNx / 2)), ...
    'timeAxis', 'ms', ...
    'radialAxis', 'rho', ...
    'colormapIndex', 1, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
species       : 'Ion', 'Alpha', 'Beam'。
plotType      : 1 径向剖面；2 时间曲线；3 二维图。
nRange        : 真实物理 n 范围；plotType = 1/2 求和，plotType = 3 作为 n 滑块范围。
timeIndex     : plotType = 1 使用的诊断时间索引（1 到 nDiagTime）。
radialIndex   : plotType = 2 使用的径向索引（1 到 gridNx）。
timeAxis      : 'ta', 'ms', 's' 或 'steps'。
radialAxis    : 'rho' 或 'x'。
colormapIndex : 二维图非负色表，可选 1-2；有正负固定红蓝。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
picWorkspace.diffusivityPlot = runPICDiffusivityPlot(picDiffusivityData, meta, picDiffusivityOpt);

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
            error('interactive 必须为 0、1 或 2。');
    end
end

function workspace = runPICOrbitFrequencyPlot(phaseSpaceOrbit, meta, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('orbit frequency plot', '绘图开关为 false');
        return;
    end
    if ~hasSpeciesData(phaseSpaceOrbit, opt.species)
        logSkipped('orbit frequency plot', '未读取对应物种的 orbit 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotOrbitFrequency(phaseSpaceOrbit, meta, opt), ...
        @(dynamicUpdate) plotOrbitFrequencyInteractive(phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICResonanceLinePlot(phaseSpaceOrbit, meta, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('resonance line plot', '绘图开关为 false');
        return;
    end
    if ~hasSpeciesData(phaseSpaceOrbit, opt.species)
        logSkipped('resonance line plot', '未读取对应物种的 orbit 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotResonanceLine(phaseSpaceOrbit, meta, opt), ...
        @(dynamicUpdate) plotResonanceLineInteractive(phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICPhaseQuantityPlot(picPhaseData, phaseSpaceOrbit, meta, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('phase quantity plot', '绘图开关为 false');
        return;
    end
    if ~hasRequestedPhaseQuantity(picPhaseData, opt.species, opt.quantity)
        logSkipped('phase quantity plot', '未读取所需 phase 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotPhaseQuantity(picPhaseData, phaseSpaceOrbit, meta, opt), ...
        @(dynamicUpdate) plotPhaseQuantityInteractive(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICPhasePowerPlot(picPhaseData, phaseSpaceOrbit, meta, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('PhasePower plot', '绘图开关为 false');
        return;
    end
    if ~hasSpeciesFieldData(picPhaseData, opt.species, 'Power')
        logSkipped('PhasePower plot', '未读取 PhasePower 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotPhasePower(picPhaseData, phaseSpaceOrbit, meta, opt), ...
        @(dynamicUpdate) plotPhasePowerInteractive(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICPitchQuantityPlot(picPitchData, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('pitch quantity plot', '绘图开关为 false');
        return;
    end
    if ~hasRequestedPhaseQuantity(picPitchData, opt.species, opt.quantity)
        logSkipped('pitch quantity plot', '未读取所需 pitch 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotPitchQuantity(picPitchData, opt), ...
        @(dynamicUpdate) plotPitchQuantityInteractive(picPitchData, opt, dynamicUpdate));
end

function workspace = runPICPitchPowerPlot(picPitchData, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('PitchPower plot', '绘图开关为 false');
        return;
    end
    if ~hasSpeciesFieldData(picPitchData, opt.species, 'Power')
        logSkipped('PitchPower plot', '未读取 PitchPower 数据');
        return;
    end

    runInteractivePlot(opt.interactive, ...
        @() plotPitchPower(picPitchData, opt), ...
        @(dynamicUpdate) plotPitchPowerInteractive(picPitchData, opt, dynamicUpdate));
end

function workspace = runPICDiffusivityPlot(picDiffusivityData, meta, opt)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped('Diffusivity plot', '绘图开关为 false');
        return;
    end
    if ~hasSpeciesFieldData(picDiffusivityData, opt.species, 'Diffusivity')
        logSkipped('Diffusivity plot', '未读取 Diffusivity 数据');
        return;
    end

    switch opt.plotType
        case 1
            runInteractivePlot(opt.interactive, ...
                @() plotDiffusivityRadial(picDiffusivityData, opt), ...
                @(dynamicUpdate) plotDiffusivityRadialInteractive(picDiffusivityData, opt, dynamicUpdate));
        case 2
            runInteractivePlot(opt.interactive, ...
                @() plotDiffusivityTime(picDiffusivityData, meta, opt), ...
                @(dynamicUpdate) plotDiffusivityTimeInteractive(picDiffusivityData, meta, opt, dynamicUpdate));
        case 3
            runInteractivePlot(opt.interactive, ...
                @() plotDiffusivityMap(picDiffusivityData, meta, opt), ...
                @(dynamicUpdate) plotDiffusivityMapInteractive(picDiffusivityData, meta, opt, dynamicUpdate));
        otherwise
            error('plotType 必须为 1、2 或 3。');
    end
end

function resonance = picResonanceOptions(enabled, branch, frequencyHz, frequencyHzRange, toroidalMode, ...
    toroidalModeRange, poloidalMode, poloidalModeRange, harmonic, harmonicRange)

    harmonicBounds = sort(reshape(double(harmonic), 1, []));
    if isscalar(harmonicBounds)
        harmonicBounds = [harmonicBounds, harmonicBounds];
    end

    resonance.enabled = enabled;
    resonance.branch = branch;
    resonance.frequencyHz = frequencyHz;
    resonance.frequencyHzRange = frequencyHzRange;
    resonance.toroidalMode = toroidalMode;
    resonance.toroidalModeRange = toroidalModeRange;
    resonance.poloidalMode = poloidalMode;
    resonance.poloidalModeRange = poloidalModeRange;
    resonance.harmonic = harmonicBounds(1);
    resonance.harmonicMin = harmonicBounds(1);
    resonance.harmonicMax = harmonicBounds(end);
    resonance.harmonicRange = harmonicRange;
end

function meta = readPICMetadata(paramText, normData)

    meta.speciesEnabled.Ion = readSwitchParam(paramText, 'ifIon');
    meta.speciesEnabled.Alpha = readSwitchParam(paramText, 'ifAlpha');
    meta.speciesEnabled.Beam = readSwitchParam(paramText, 'ifBeam');
    meta.switch.ifDiagDiffusivity = readSwitchParam(paramText, 'ifDiagDiffusivity');
    meta.switch.ifOutputPhaseSpaceOrbit = readSwitchParam(paramText, 'ifOutputPhaseSpaceOrbit');
    meta.switch.ifOutputPhaseSpaceJacobian = readSwitchParam(paramText, 'ifOutputPhaseSpaceJacobian');
    meta.switch.ifOutputPhaseSpaceF0 = readSwitchParam(paramText, 'ifOutputPhaseSpaceF0');
    meta.switch.ifOutputPhaseSpaceDeltaF = readSwitchParam(paramText, 'ifOutputPhaseSpaceDeltaF');
    meta.switch.ifOutputPhaseSpacePower = readSwitchParam(paramText, 'ifOutputPhaseSpacePower');
    meta.switch.ifOutputPitchSpaceJacobian = readSwitchParam(paramText, 'ifOutputPitchSpaceJacobian');
    meta.switch.ifOutputPitchSpaceF0 = readSwitchParam(paramText, 'ifOutputPitchSpaceF0');
    meta.switch.ifOutputPitchSpaceDeltaF = readSwitchParam(paramText, 'ifOutputPitchSpaceDeltaF');
    meta.switch.ifOutputPitchSpacePower = readSwitchParam(paramText, 'ifOutputPitchSpacePower');

    meta.paramGridE = readIntParam(paramText, 'gridE');
    meta.paramGridPphi = readIntParam(paramText, 'gridPphi');
    meta.paramGridLambda = readIntParam(paramText, 'gridLambda');
    meta.gridVpara = readIntParam(paramText, 'gridVpara');
    meta.gridVperp = readIntParam(paramText, 'gridVperp');
    meta.gridNx = readIntParam(paramText, 'gridNx');
    meta.leftN = readIntParam(paramText, 'leftN');
    meta.rightN = readIntParam(paramText, 'rightN');
    meta.tubes = readIntParam(paramText, 'tubes');
    assert(meta.leftN <= meta.rightN, 'leftN 必须小于或等于 rightN。');
    assert(meta.tubes > 0, 'tubes 必须为正整数。');
    meta.modeIndexAll = meta.leftN:meta.rightN;
    meta.physicalNAll = meta.modeIndexAll * meta.tubes;

    [meta.gridE, meta.gridPphi, meta.gridLambda] = readPhaseGrid(normData);
    assert(meta.paramGridE == meta.gridE && meta.paramGridPphi == meta.gridPphi && ...
        meta.paramGridLambda == meta.gridLambda, ...
        'cuGMEC_param.h 与 normalization2D.mat 中的 gridE/gridPphi/gridLambda 不一致。');

    meta.totalSteps = readIntParam(paramText, 'totalSteps');
    meta.outputSteps = readIntParam(paramText, 'outputSteps');
    meta.diagSteps = readIntParam(paramText, 'diagSteps');
    meta.dt = readFloatParam(paramText, 'dt');
    assert(meta.gridNx > 1, 'gridNx 必须大于 1。');
    assert(meta.totalSteps >= 0 && meta.outputSteps > 0 && meta.diagSteps > 0, ...
        'totalSteps/outputSteps/diagSteps 取值不合法。');
    meta.nOutputTime = floor(meta.totalSteps / meta.outputSteps) + 1;
    meta.nDiagTime = floor(meta.totalSteps / meta.diagSteps) + 1;
    meta.mhdPrecision = readMHDPrecisionParam(paramText);
    meta.L0 = readPositiveScalar(normData, 'L0');
    meta.VA0 = readPositiveScalar(normData, 'VA0');
    assert(isfield(normData, 'RHO0') && isfield(normData, 'RHO1'), ...
        'normalization2D.mat 缺少 RHO0 或 RHO1。');
    meta.RHO0 = validateFiniteScalar(normData.RHO0, 'RHO0');
    meta.RHO1 = validateFiniteScalar(normData.RHO1, 'RHO1');
    assert(meta.RHO0 < meta.RHO1, 'RHO0 必须小于 RHO1。');
    meta.tDiag = (0:meta.nDiagTime - 1) * meta.diagSteps * meta.dt;
    meta.timeSeconds = meta.tDiag * meta.L0 / meta.VA0;
    meta.xGrid = linspace(0, 1, meta.gridNx);
    meta.rhoGrid = linspace(meta.RHO0, meta.RHO1, meta.gridNx);

end

function speciesList = enabledPICSpecies(speciesList, meta)

    keep = false(size(speciesList));
    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        keep(speciesIndex) = isfield(meta.speciesEnabled, speciesName) && meta.speciesEnabled.(speciesName);
        if ~keep(speciesIndex)
            logSkipped(speciesName, '物种开关为 false');
        end
    end
    speciesList = speciesList(keep);
end

function picPhaseData = initializePhaseSpaceData(speciesList, normData, meta)

    picPhaseData = struct();
    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        phaseRange = readSpeciesRange(normData, [speciesName 'EPphiLambda']);

        picPhaseData.(speciesName) = struct( ...
            'species', speciesName, ...
            'gridE', meta.gridE, ...
            'gridPphi', meta.gridPphi, ...
            'gridLambda', meta.gridLambda, ...
            'modeIndexAll', meta.modeIndexAll, ...
            'physicalNAll', meta.physicalNAll, ...
            'tubes', meta.tubes, ...
            'EPphiLambda', phaseRange, ...
            'E1d', linspace(phaseRange(1), phaseRange(2), meta.gridE), ...
            'Pphi1d', linspace(phaseRange(3), phaseRange(4), meta.gridPphi), ...
            'Lambda1d', linspace(phaseRange(5), phaseRange(6), meta.gridLambda), ...
            'J', [], ...
            'F0', [], ...
            'DF', [], ...
            'Power', []);
    end
end

function picPitchData = initializePitchSpaceData(speciesList, paramText, meta)

    picPitchData = struct();
    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        vmax = readFloatParam(paramText, [speciesName 'Vmax']);
        if strcmp(speciesName, 'Beam')
            minVpara = 0;
        else
            minVpara = -vmax;
        end

        picPitchData.(speciesName) = struct( ...
            'species', speciesName, ...
            'gridVpara', meta.gridVpara, ...
            'gridVperp', meta.gridVperp, ...
            'modeIndexAll', meta.modeIndexAll, ...
            'physicalNAll', meta.physicalNAll, ...
            'tubes', meta.tubes, ...
            'Vpara1d', linspace(minVpara, vmax, meta.gridVpara), ...
            'Vperp1d', linspace(0, vmax, meta.gridVperp), ...
            'J', [], ...
            'F0', [], ...
            'DF', [], ...
            'Power', []);
    end
end

function picDiffusivityData = initializeDiffusivityData(speciesList, meta)

    picDiffusivityData = struct();
    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});

        picDiffusivityData.(speciesName) = struct( ...
            'species', speciesName, ...
            'gridNx', meta.gridNx, ...
            'modeIndexAll', meta.modeIndexAll, ...
            'physicalNAll', meta.physicalNAll, ...
            'tubes', meta.tubes, ...
            'diagSteps', meta.diagSteps, ...
            'xGrid', meta.xGrid, ...
            'rhoGrid', meta.rhoGrid, ...
            'tDiag', meta.tDiag, ...
            'timeSeconds', meta.timeSeconds, ...
            'Diffusivity', []);
    end
end

function orbitRaw = readAllOrbitRaw(speciesList, searchDirs, meta)

    orbitRaw = struct();
    if ~meta.switch.ifOutputPhaseSpaceOrbit
        logSkipped('PhaseSpaceOrbit', '开关 ifOutputPhaseSpaceOrbit 为 false');
        return;
    end

    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        fileName = [speciesName 'PhaseSpaceOrbit.bin'];
        filePath = findExistingFile(searchDirs, fileName);
        if isempty(filePath)
            orbitRaw.(speciesName) = [];
            logSkipped(fileName, '文件不存在');
            continue;
        end

        orbitRaw.(speciesName) = struct('file', filePath, 'data', readOrbitBinary(filePath));
        logLoaded(fileName, orbitRaw.(speciesName).data);
    end
end

function picPhaseData = readAllPhaseDiagnostics(picPhaseData, speciesList, searchDirs, meta)

    phaseSpecs = { ...
        'J', 'ifOutputPhaseSpaceJacobian', 'PhaseSpaceJacobian', 'double', 'phase3d'; ...
        'F0', 'ifOutputPhaseSpaceF0', 'PhaseSpaceF0', 'double', 'phase3d'; ...
        'DF', 'ifOutputPhaseSpaceDeltaF', 'PhaseDeltaF', meta.mhdPrecision, 'phase4d'; ...
        'Power', 'ifOutputPhaseSpacePower', 'PhasePower', meta.mhdPrecision, 'phasePower'};

    for iSpec = 1:size(phaseSpecs, 1)
        fieldName = phaseSpecs{iSpec, 1};
        switchName = phaseSpecs{iSpec, 2};
        fileSuffix = phaseSpecs{iSpec, 3};
        precision = phaseSpecs{iSpec, 4};
        dataKind = phaseSpecs{iSpec, 5};

        if ~meta.switch.(switchName)
            logSkipped(fileSuffix, ['开关 ' switchName ' 为 false']);
            continue;
        end

        for speciesIndex = 1:numel(speciesList)
            speciesName = char(speciesList{speciesIndex});
            fileName = [speciesName fileSuffix '.bin'];
            filePath = findExistingFile(searchDirs, fileName);
            if isempty(filePath)
                logSkipped(fileName, '文件不存在');
                continue;
            end

            switch dataKind
                case 'phase3d'
                    picPhaseData.(speciesName).(fieldName) = readPhase3D(filePath, ...
                        precision, meta.gridE, meta.gridPphi, meta.gridLambda);
                case 'phase4d'
                    picPhaseData.(speciesName).(fieldName) = readPhase4D(filePath, ...
                        precision, meta.gridE, meta.gridPphi, meta.gridLambda, meta.nOutputTime);
                case 'phasePower'
                    picPhaseData.(speciesName).(fieldName) = readPhasePower5D(filePath, ...
                        precision, meta.gridE, meta.gridPphi, meta.gridLambda, numel(meta.modeIndexAll), meta.nOutputTime);
                otherwise
                    error('未知 phase 诊断类型：%s。', dataKind);
            end
            logLoaded(fileName, picPhaseData.(speciesName).(fieldName));
        end
    end
end

function picPitchData = readAllPitchDiagnostics(picPitchData, speciesList, searchDirs, meta)

    pitchSpecs = { ...
        'J', 'ifOutputPitchSpaceJacobian', 'PitchSpaceJacobian', 'double', 'pitch2d'; ...
        'F0', 'ifOutputPitchSpaceF0', 'PitchSpaceF0', 'double', 'pitch2d'; ...
        'DF', 'ifOutputPitchSpaceDeltaF', 'PitchDeltaF', meta.mhdPrecision, 'pitch3d'; ...
        'Power', 'ifOutputPitchSpacePower', 'PitchPower', meta.mhdPrecision, 'pitchPower'};

    for iSpec = 1:size(pitchSpecs, 1)
        fieldName = pitchSpecs{iSpec, 1};
        switchName = pitchSpecs{iSpec, 2};
        fileSuffix = pitchSpecs{iSpec, 3};
        precision = pitchSpecs{iSpec, 4};
        dataKind = pitchSpecs{iSpec, 5};

        if ~meta.switch.(switchName)
            logSkipped(fileSuffix, ['开关 ' switchName ' 为 false']);
            continue;
        end

        for speciesIndex = 1:numel(speciesList)
            speciesName = char(speciesList{speciesIndex});
            fileName = [speciesName fileSuffix '.bin'];
            filePath = findExistingFile(searchDirs, fileName);
            if isempty(filePath)
                logSkipped(fileName, '文件不存在');
                continue;
            end

            switch dataKind
                case 'pitch2d'
                    picPitchData.(speciesName).(fieldName) = readPitch2D(filePath, ...
                        precision, meta.gridVpara, meta.gridVperp);
                case 'pitch3d'
                    picPitchData.(speciesName).(fieldName) = readPitch3D(filePath, ...
                        precision, meta.gridVpara, meta.gridVperp, meta.nOutputTime);
                case 'pitchPower'
                    picPitchData.(speciesName).(fieldName) = readPitchPower4D(filePath, ...
                        precision, meta.gridVpara, meta.gridVperp, numel(meta.modeIndexAll), meta.nOutputTime);
                otherwise
                    error('未知 pitch 诊断类型：%s。', dataKind);
            end
            logLoaded(fileName, picPitchData.(speciesName).(fieldName));
        end
    end
end

function picDiffusivityData = readAllDiffusivityDiagnostics(picDiffusivityData, speciesList, searchDirs, meta)

    if ~meta.switch.ifDiagDiffusivity
        logSkipped('Diffusivity', '开关 ifDiagDiffusivity 为 false');
        return;
    end

    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        fileName = [speciesName 'Diffusivity.bin'];
        filePath = findExistingFile(searchDirs, fileName);
        if isempty(filePath)
            logSkipped(fileName, '文件不存在');
            continue;
        end

        picDiffusivityData.(speciesName).Diffusivity = readDiffusivity3D(filePath, ...
            meta.mhdPrecision, meta.nDiagTime, numel(meta.modeIndexAll), meta.gridNx);
        logLoaded(fileName, picDiffusivityData.(speciesName).Diffusivity);
    end
end

function [phaseSpaceOrbit, PhaseSpaceOrbitSummary] = processAllOrbitRaw(orbitRaw, speciesList, picPhaseData, opt)

    phaseSpaceOrbit = struct();
    PhaseSpaceOrbitSummary = struct();
    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        if ~isfield(orbitRaw, speciesName) || isempty(orbitRaw.(speciesName))
            continue;
        end

        fprintf('\n============================================================\n');
        fprintf('%s phase-space orbit analysis\n', speciesName);
        fprintf('============================================================\n');

        result = analyzeSpeciesOrbitRecords(speciesName, orbitRaw.(speciesName), picPhaseData.(speciesName));
        phaseSpaceOrbit.(speciesName) = result.phaseSpaceData;
        PhaseSpaceOrbitSummary.(speciesName) = result.summary;
        printSpeciesSummary(result.summary);

        if optionOrDefault(opt, 'plotConservationDiagnostics', false)
            plotConservationDiagnostics(speciesName, result.diagnostics, ...
                result.phaseSpaceData.E1d, result.phaseSpaceData.Pphi1d, result.phaseSpaceData.Lambda1d);
        end
    end
end

function result = analyzeSpeciesOrbitRecords(speciesName, rawOrbit, phaseGrid)

    records = orbitRecordColumns(rawOrbit.data);
    assertNoNaN(speciesName, records);

    nPhase = phaseGrid.gridE * phaseGrid.gridPphi * phaseGrid.gridLambda;
    expectedRecords = 2 * nPhase;
    numRecords = numel(records.Ids);
    assert(numRecords == expectedRecords, ...
        '%s 尺寸不匹配：读到 %d 条记录，期望 %d 条（2*gridE*gridPphi*gridLambda，gridE=%d, gridPphi=%d, gridLambda=%d）。', ...
        rawOrbit.file, numRecords, expectedRecords, phaseGrid.gridE, phaseGrid.gridPphi, phaseGrid.gridLambda);

    phaseSpaceData = emptyOrbitPhaseData(phaseGrid);
    recordBranch = ones(numRecords, 1);
    recordBranch(nPhase + 1:end) = -1;

    validRecord = ~isPadRecord(records.Ids);
    localId = abs(double(records.Ids));
    validLocalId = validRecord & localId >= 0 & localId < nPhase & localId == floor(localId);
    if any(validRecord & ~validLocalId)
        badIndex = find(validRecord & ~validLocalId, 1);
        error('%s 中存在非法相空间 id：record=%d, id=%d。', rawOrbit.file, badIndex, records.Ids(badIndex));
    end

    warnIfSignedIdOrderLooksWrong(speciesName, records.Ids, recordBranch, validRecord);
    [plusIndex, minusIndex, duplicateCounts] = buildBranchIndex(localId, recordBranch, validLocalId, nPhase);

    diagnostics = struct();
    summary = initializeOrbitSummary(speciesName, rawOrbit.file, numRecords, nPhase, validRecord, records.orbits);
    summary.duplicatePlusIds = duplicateCounts.plus;
    summary.duplicateMinusIds = duplicateCounts.minus;

    [phaseSpaceData.trapped, diagnostics.trapped, summary.trapped] = extractTrappedOrbits( ...
        phaseSpaceData.trapped, plusIndex, minusIndex, records, phaseGrid);
    [phaseSpaceData.para, diagnostics.para, summary.para] = extractPassingOrbits( ...
        2.5, phaseSpaceData.para, validLocalId, localId, recordBranch, plusIndex, minusIndex, records, phaseGrid);
    [phaseSpaceData.anti, diagnostics.anti, summary.anti] = extractPassingOrbits( ...
        3.5, phaseSpaceData.anti, validLocalId, localId, recordBranch, plusIndex, minusIndex, records, phaseGrid);

    summary.classifiedBranchCount = 2 * summary.trapped.count + summary.para.count + summary.anti.count;
    summary.classifiedRatioTotal = safeDivide(summary.classifiedBranchCount, expectedRecords);
    summary.classifiedRatioInitialized = safeDivide(summary.classifiedBranchCount, summary.initializedRecords);
    summary.unclassifiedInitializedRecords = summary.initializedRecords - summary.classifiedBranchCount;
    summary.unclassifiedInitializedRatio = safeDivide(summary.unclassifiedInitializedRecords, summary.initializedRecords);

    result = struct('phaseSpaceData', phaseSpaceData, 'diagnostics', diagnostics, 'summary', summary);
end

function phaseSpaceData = emptyOrbitPhaseData(phaseGrid)

    phaseSpaceData = struct( ...
        'species', phaseGrid.species, ...
        'gridE', phaseGrid.gridE, ...
        'gridPphi', phaseGrid.gridPphi, ...
        'gridLambda', phaseGrid.gridLambda, ...
        'EPphiLambda', phaseGrid.EPphiLambda, ...
        'E1d', phaseGrid.E1d, ...
        'Pphi1d', phaseGrid.Pphi1d, ...
        'Lambda1d', phaseGrid.Lambda1d, ...
        'trapped', emptyOrbitClass(phaseGrid), ...
        'para', emptyOrbitClass(phaseGrid), ...
        'anti', emptyOrbitClass(phaseGrid));
end

function classData = emptyOrbitClass(phaseGrid)

    zeroArray = zeros(phaseGrid.gridE, phaseGrid.gridPphi, phaseGrid.gridLambda);
    classData = struct('dtheta', zeroArray, 'dphiTotal', zeroArray, 'dphiVpara', zeroArray, 'dT', zeroArray);
end

function [plusIndex, minusIndex, duplicateCounts] = buildBranchIndex(localId, recordBranch, validRecord, nPhase)

    plusIndex = zeros(nPhase, 1);
    minusIndex = zeros(nPhase, 1);
    plusRecords = find(validRecord & recordBranch > 0);
    minusRecords = find(validRecord & recordBranch < 0);
    plusLocalIndex = localId(plusRecords) + 1;
    minusLocalIndex = localId(minusRecords) + 1;

    duplicateCounts = struct( ...
        'plus', numel(plusLocalIndex) - numel(unique(plusLocalIndex)), ...
        'minus', numel(minusLocalIndex) - numel(unique(minusLocalIndex)));
    plusIndex(plusLocalIndex) = plusRecords;
    minusIndex(minusLocalIndex) = minusRecords;
end

function [classData, diagnostic, summary] = extractTrappedOrbits(classData, plusIndex, minusIndex, records, phaseGrid)

    pairedLocalIndex = find(plusIndex > 0 & minusIndex > 0);
    plusRecord = plusIndex(pairedLocalIndex);
    minusRecord = minusIndex(pairedLocalIndex);
    trappedOrbit = isOrbit(records.orbits(plusRecord), 4.5) & isOrbit(records.orbits(minusRecord), 4.5);
    dTRelDiff = relativePairDifference(records.dTs(plusRecord), records.dTs(minusRecord));
    accepted = trappedOrbit & dTRelDiff < trappedOrbitRelativeTolerance();

    acceptedLocalId = pairedLocalIndex(accepted) - 1;
    acceptedPlusRecord = plusRecord(accepted);
    acceptedMinusRecord = minusRecord(accepted);
    classData = fillOrbitClass(classData, acceptedLocalId, ...
        averagePair(records.dtheta(acceptedPlusRecord), records.dtheta(acceptedMinusRecord)), ...
        averagePair(records.dphiTotal(acceptedPlusRecord), records.dphiTotal(acceptedMinusRecord)), ...
        averagePair(records.dphiVpara(acceptedPlusRecord), records.dphiVpara(acceptedMinusRecord)), ...
        averagePair(records.dTs(acceptedPlusRecord), records.dTs(acceptedMinusRecord)), phaseGrid);

    diagnostic = struct( ...
        'localId', acceptedLocalId(:), ...
        'E', averagePair(records.Es(acceptedPlusRecord), records.Es(acceptedMinusRecord)), ...
        'Pphi', averagePair(records.Pphis(acceptedPlusRecord), records.Pphis(acceptedMinusRecord)), ...
        'Lambda', averagePair(records.Lambdas(acceptedPlusRecord), records.Lambdas(acceptedMinusRecord)), ...
        'dTRelDiff', dTRelDiff(accepted), ...
        'plusRecord', acceptedPlusRecord(:), ...
        'minusRecord', acceptedMinusRecord(:));

    summary = struct( ...
        'count', numel(acceptedLocalId), ...
        'branchCount', 2 * numel(acceptedLocalId), ...
        'pairedCandidateCount', numel(pairedLocalIndex), ...
        'bothTrappedCount', sum(trappedOrbit), ...
        'rejectedByDTCount', sum(trappedOrbit & ~accepted), ...
        'maxAcceptedDTRelDiff', maxOrNaN(dTRelDiff(accepted)), ...
        'meanAcceptedDTRelDiff', meanOrNaN(dTRelDiff(accepted)));
end

function [classData, diagnostic, summary] = extractPassingOrbits(targetOrbit, classData, validRecord, localId, ...
    recordBranch, plusIndex, minusIndex, records, phaseGrid)

    targetRecord = find(validRecord & isOrbit(records.orbits, targetOrbit));
    targetLocalId = localId(targetRecord);
    targetBranch = recordBranch(targetRecord);
    counterpartRecord = zeros(size(targetRecord));
    plusTarget = targetBranch > 0;
    minusTarget = targetBranch < 0;
    counterpartRecord(plusTarget) = minusIndex(targetLocalId(plusTarget) + 1);
    counterpartRecord(minusTarget) = plusIndex(targetLocalId(minusTarget) + 1);

    hasCounterpart = counterpartRecord > 0;
    counterpartIsTrapped = false(size(targetRecord));
    counterpartIsTrapped(hasCounterpart) = isOrbit(records.orbits(counterpartRecord(hasCounterpart)), 4.5);
    keep = ~counterpartIsTrapped;
    keptRecord = targetRecord(keep);
    keptLocalId = targetLocalId(keep);
    keptHasCounterpart = hasCounterpart(keep);
    [uniqueLocalId, uniquePosition] = unique(keptLocalId, 'stable');
    selectedRecord = keptRecord(uniquePosition);
    selectedHasCounterpart = keptHasCounterpart(uniquePosition);

    classData = fillOrbitClass(classData, uniqueLocalId, ...
        records.dtheta(selectedRecord), records.dphiTotal(selectedRecord), ...
        records.dphiVpara(selectedRecord), records.dTs(selectedRecord), phaseGrid);

    diagnostic = struct( ...
        'localId', uniqueLocalId(:), ...
        'E', records.Es(selectedRecord), ...
        'Pphi', records.Pphis(selectedRecord), ...
        'Lambda', records.Lambdas(selectedRecord), ...
        'record', selectedRecord(:));

    summary = struct( ...
        'count', numel(selectedRecord), ...
        'singleCount', sum(~selectedHasCounterpart), ...
        'pairedCount', sum(selectedHasCounterpart), ...
        'candidateCount', numel(targetRecord), ...
        'rejectedByTrappedCounterpartCount', sum(counterpartIsTrapped), ...
        'duplicateTargetIdCount', numel(keptLocalId) - numel(uniqueLocalId));
end

function classData = fillOrbitClass(classData, localIds, dthetaValues, dphiTotalValues, dphiVparaValues, dTValues, phaseGrid)

    if isempty(localIds)
        return;
    end

    linearIndex = localIdsToLinear(localIds, phaseGrid.gridE, phaseGrid.gridPphi, phaseGrid.gridLambda);
    classData.dtheta(linearIndex) = dthetaValues(:);
    classData.dphiTotal(linearIndex) = dphiTotalValues(:);
    classData.dphiVpara(linearIndex) = dphiVparaValues(:);
    classData.dT(linearIndex) = dTValues(:);
end

function plotOrbitFrequency(phaseSpaceOrbit, meta, opt)

    [orbitData, speciesLabel] = resolveSpeciesData(phaseSpaceOrbit, opt.species, 'orbit');
    branchName = normalizeFrequencyBranch(opt.branch);
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx = parsePhaseSlice(opt.slice, dim, orbitData);
    [unitScale, colorbarLabel] = frequencyUnitScale(opt.unit);
    frequency = calculateOrbitFrequencyField(orbitData, branchName, unitScale, meta);
    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField(dim, idx, frequency, orbitData);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s %s orbit frequency 在该切片中没有有限值数据。\n', speciesLabel, branchName);
        return;
    end

    Z = fillEnclosedBlankRegions(Z);
    titleText = simpleTitleText(speciesLabel, orbitFrequencyLatex(branchName), sliceTitle, '');
    plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, colorbarLabel, ...
        opt.colormapIndex, opt.contourCount, sprintf('%s %s orbit frequency', speciesLabel, branchName));
    drawMapFigure(sprintf('%s %s orbit frequency', speciesLabel, branchName), plotData);
end

function plotOrbitFrequencyInteractive(phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [orbitData, speciesLabel] = resolveSpeciesData(phaseSpaceOrbit, opt.species, 'orbit');
    branchName = normalizeFrequencyBranch(opt.branch);
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx0 = parsePhaseSlice(opt.slice, dim, orbitData);
    [unitScale, colorbarLabel] = frequencyUnitScale(opt.unit);
    frequency = calculateOrbitFrequencyField(orbitData, branchName, unitScale, meta);
    nSlice = phaseDimensionSize(dim, orbitData);
    contourCount = optionOrDefault(opt, 'contourCount', 0);
    controls = [ ...
        integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], idx0, 1, nSlice), ...
        integerSliderControl('contourCount', 'contours', contourCount, 0, max(40, 2 * contourCount))];

    plotInteractiveMap(sprintf('%s %s orbit frequency', speciesLabel, branchName), controls, dynamicUpdate, ...
        @computePlotData, @renderMap);

    function plotData = computePlotData(values)
        [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
            dim, values.sliceIndex, frequency, orbitData);
        if any(isfinite(Z(:)))
            Z = fillEnclosedBlankRegions(Z);
        end
        titleText = simpleTitleText(speciesLabel, orbitFrequencyLatex(branchName), sliceTitle, '');
        statusText = sprintf('%s %s orbit frequency, %s index = %d, contours = %d', ...
            speciesLabel, branchName, char(opt.fixedCoordinate), values.sliceIndex, values.contourCount);
        plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, colorbarLabel, ...
            opt.colormapIndex, values.contourCount, statusText);
    end
end

function plotResonanceLine(phaseSpaceOrbit, meta, opt)

    [orbitData, speciesLabel] = resolveSpeciesData(phaseSpaceOrbit, opt.species, 'orbit');
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx = parsePhaseSlice(opt.slice, dim, orbitData);
    resOpt = validatedResonanceOptions(opt);
    [residualHz, physicalN] = calculateResonanceResidualField(orbitData, resOpt, meta);
    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField(dim, idx, residualHz, orbitData);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s %s resonance 在该切片中没有有限值数据。\n', speciesLabel, resOpt.branch);
        return;
    end

    titleText = resonanceTitleText(speciesLabel, resOpt, physicalN, sliceTitle);
    plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        '$\Delta f_{\mathrm{res}}/\mathrm{Hz}$', opt.colormapIndex, 0, ...
        resonanceStatusText(speciesLabel, resOpt, physicalN, residualHasZeroContour(Z)));
    plotData.forceSigned = true;
    plotData.resonanceZ = Z;
    drawMapFigure(sprintf('%s %s resonance residual', speciesLabel, resOpt.branch), plotData);
end

function plotResonanceLineInteractive(phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [orbitData, speciesLabel] = resolveSpeciesData(phaseSpaceOrbit, opt.species, 'orbit');
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx0 = parsePhaseSlice(opt.slice, dim, orbitData);
    resOpt0 = validatedResonanceOptions(opt);
    controls = [integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], ...
        idx0, 1, phaseDimensionSize(dim, orbitData))];
    controls = appendResonanceControls(controls, resOpt0);

    plotInteractiveMap(sprintf('%s %s resonance residual', speciesLabel, resOpt0.branch), controls, dynamicUpdate, ...
        @computePlotData, @renderMap);

    function plotData = computePlotData(values)
        resOpt = resonanceOptionsFromValues(resOpt0, values);
        [residualHz, physicalN] = calculateResonanceResidualField(orbitData, resOpt, meta);
        [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
            dim, values.sliceIndex, residualHz, orbitData);

        titleText = resonanceTitleText(speciesLabel, resOpt, physicalN, sliceTitle);
        hasZero = residualHasZeroContour(Z);
        plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
            '$\Delta f_{\mathrm{res}}/\mathrm{Hz}$', opt.colormapIndex, 0, ...
            resonanceStatusText(speciesLabel, resOpt, physicalN, hasZero));
        plotData.forceSigned = true;
        plotData.resonanceZ = Z;
    end
end

function plotPhaseQuantity(picPhaseData, phaseSpaceOrbit, meta, opt)

    [speciesData, speciesLabel] = resolveSpeciesData(picPhaseData, opt.species, 'phase');
    [quantityData, quantityField, quantityLatex, isTimeDependent] = resolvePhaseQuantity(speciesData, opt.quantity);
    [fieldData, timeIndexText] = selectQuantityFrame(quantityData, isTimeDependent, opt.timeIndex);
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx = parsePhaseSlice(opt.slice, dim, speciesData);
    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField(dim, idx, fieldData, speciesData);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s %s 在该切片中没有有限值数据。\n', speciesLabel, quantityField);
        return;
    end

    titleText = phaseQuantityTitleText(speciesLabel, quantityLatex, sliceTitle, timeIndexText, opt.resonance);
    statusText = sprintf('%s %s, %s index = %d%s', speciesLabel, quantityField, char(opt.fixedCoordinate), idx, timeIndexText);
    plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantityLatex '$'], opt.colormapIndex, opt.contourCount, statusText);
    plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, idx, opt.resonance, meta);
    drawMapFigure(sprintf('%s %s phase-space slice', speciesLabel, quantityField), plotData);
end

function plotPhaseQuantityInteractive(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPhaseData, opt.species, 'phase');
    [quantityData, quantityField, quantityLatex, isTimeDependent] = resolvePhaseQuantity(speciesData, opt.quantity);
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx0 = parsePhaseSlice(opt.slice, dim, speciesData);
    contourCount = optionOrDefault(opt, 'contourCount', 0);
    controls = integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], ...
        idx0, 1, phaseDimensionSize(dim, speciesData));

    if isTimeDependent
        controls = [controls, integerSliderControl('timeIndex', 'timeIndex', ...
            parseTimeIndex(opt.timeIndex, timeDimensionSize(quantityData)), 1, timeDimensionSize(quantityData))];
    end
    controls = [controls, integerSliderControl('contourCount', 'contours', contourCount, 0, max(20, 2 * contourCount))];
    controls = appendResonanceControls(controls, opt.resonance, true, true);

    plotInteractiveMap(sprintf('%s %s phase-space slice', speciesLabel, quantityField), controls, dynamicUpdate, ...
        @computePlotData, @renderMap);

    function plotData = computePlotData(values)
        tempTimeIndex = optionOrDefault(opt, 'timeIndex', 1);
        if isTimeDependent
            tempTimeIndex = values.timeIndex;
        end
        [fieldData, timeIndexText] = selectQuantityFrame(quantityData, isTimeDependent, tempTimeIndex);
        [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
            dim, values.sliceIndex, fieldData, speciesData);
        resOpt = resonanceOptionsFromValues(opt.resonance, values, true);
        titleText = phaseQuantityTitleText(speciesLabel, quantityLatex, sliceTitle, timeIndexText, resOpt);
        statusText = sprintf('%s %s, %s index = %d%s, contours = %d', ...
            speciesLabel, quantityField, char(opt.fixedCoordinate), values.sliceIndex, timeIndexText, values.contourCount);
        plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
            ['$' quantityLatex '$'], opt.colormapIndex, values.contourCount, statusText);
        plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, values.sliceIndex, resOpt, meta);
    end
end

function plotPhasePower(picPhaseData, phaseSpaceOrbit, meta, opt)

    [speciesData, speciesLabel] = resolveSpeciesData(picPhaseData, opt.species, 'phase');
    powerData = requireQuantityData(speciesData, 'Power');
    [modeIdx, modeN] = parseModeN(opt.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    timeIndex = parseTimeIndex(opt.timeIndex, size(powerData, 5));
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx = parsePhaseSlice(opt.slice, dim, speciesData);
    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField(dim, idx, powerData(:, :, :, modeIdx, timeIndex), speciesData);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s PhasePower 在该切片中没有有限值数据。\n', speciesLabel);
        return;
    end

    quantityLatex = 'P_{\mathrm{PIC}}';
    timeIndexText = sprintf(', timeIndex = %d', timeIndex);
    titleText = phasePowerTitleText(speciesLabel, quantityLatex, sliceTitle, modeN, timeIndexText, opt.resonance);
    statusText = sprintf('%s PhasePower, %s index = %d%s%s', ...
        speciesLabel, char(opt.fixedCoordinate), idx, modeStatusText(modeN), timeIndexText);
    plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantityLatex '$'], opt.colormapIndex, opt.contourCount, statusText);
    plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, idx, opt.resonance, meta);
    drawMapFigure(sprintf('%s PhasePower phase-space slice', speciesLabel), plotData);
end

function plotPhasePowerInteractive(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPhaseData, opt.species, 'phase');
    powerData = requireQuantityData(speciesData, 'Power');
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    idx0 = parsePhaseSlice(opt.slice, dim, speciesData);
    [~, modeN0] = parseModeN(opt.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    timeIndex0 = parseTimeIndex(opt.timeIndex, size(powerData, 5));
    contourCount = optionOrDefault(opt, 'contourCount', 0);
    controls = integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], ...
        idx0, 1, phaseDimensionSize(dim, speciesData));
    if numel(speciesData.modeIndexAll) > 1
        modeControl = integerSliderControl('modeN', 'n', modeN0, ...
            min(speciesData.physicalNAll), max(speciesData.physicalNAll));
        modeControl.allowedValues = speciesData.physicalNAll;
        controls = [controls, modeControl];
    end
    if size(powerData, 5) > 1
        controls = [controls, integerSliderControl('timeIndex', 'timeIndex', timeIndex0, 1, size(powerData, 5))];
    end
    controls = [controls, integerSliderControl('contourCount', 'contours', contourCount, 0, max(20, 2 * contourCount))];
    controls = appendResonanceControls(controls, opt.resonance, true);

    plotInteractiveMap(sprintf('%s PhasePower phase-space slice', speciesLabel), controls, dynamicUpdate, ...
        @computePlotData, @renderMap);

    function plotData = computePlotData(values)
        tempModeN = modeN0;
        if isfield(values, 'modeN')
            tempModeN = values.modeN;
        end
        [tempModeIdx, tempModeN] = parseModeN(tempModeN, speciesData.modeIndexAll, speciesData.physicalNAll);
        tempTimeIndex = timeIndex0;
        if isfield(values, 'timeIndex')
            tempTimeIndex = values.timeIndex;
        end
        [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
            dim, values.sliceIndex, powerData(:, :, :, tempModeIdx, tempTimeIndex), speciesData);
        resOpt = resonanceOptionsFromValues(opt.resonance, values);
        quantityLatex = 'P_{\mathrm{PIC}}';
        timeIndexText = sprintf(', timeIndex = %d', tempTimeIndex);
        titleText = phasePowerTitleText(speciesLabel, quantityLatex, sliceTitle, tempModeN, timeIndexText, resOpt);
        statusText = sprintf('%s PhasePower, %s index = %d%s%s, contours = %d', ...
            speciesLabel, char(opt.fixedCoordinate), values.sliceIndex, ...
            modeStatusText(tempModeN), timeIndexText, values.contourCount);
        plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
            ['$' quantityLatex '$'], opt.colormapIndex, values.contourCount, statusText);
        plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, values.sliceIndex, resOpt, meta);
    end
end

function plotPitchQuantity(picPitchData, opt)

    [speciesData, speciesLabel] = resolveSpeciesData(picPitchData, opt.species, 'pitch');
    [quantityData, quantityField, quantityLatex, isTimeDependent] = resolvePhaseQuantity(speciesData, opt.quantity);
    [fieldData, timeIndexText] = selectQuantityFrame(quantityData, isTimeDependent, opt.timeIndex);
    [Z, xVec, yVec, xlabelText, ylabelText] = pitchMapField(fieldData, speciesData);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s pitch %s 没有有限值数据。\n', speciesLabel, quantityField);
        return;
    end

    titleText = pitchTitleText(speciesLabel, quantityLatex, timeIndexText);
    plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantityLatex '$'], opt.colormapIndex, opt.contourCount, ...
        sprintf('%s pitch %s%s', speciesLabel, quantityField, timeIndexText));
    drawMapFigure(sprintf('%s %s pitch-space map', speciesLabel, quantityField), plotData);
end

function plotPitchQuantityInteractive(picPitchData, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPitchData, opt.species, 'pitch');
    [quantityData, quantityField, quantityLatex, isTimeDependent] = resolvePhaseQuantity(speciesData, opt.quantity);
    contourCount = optionOrDefault(opt, 'contourCount', 0);
    controls = integerSliderControl('contourCount', 'contours', contourCount, 0, max(20, 2 * contourCount));
    if isTimeDependent
        controls = [integerSliderControl('timeIndex', 'timeIndex', ...
            parseTimeIndex(opt.timeIndex, timeDimensionSize(quantityData)), 1, timeDimensionSize(quantityData)), controls];
    end

    plotInteractiveMap(sprintf('%s %s pitch-space map', speciesLabel, quantityField), controls, dynamicUpdate, ...
        @computePlotData, @renderMap);

    function plotData = computePlotData(values)
        tempTimeIndex = optionOrDefault(opt, 'timeIndex', 1);
        if isTimeDependent
            tempTimeIndex = values.timeIndex;
        end
        [fieldData, timeIndexText] = selectQuantityFrame(quantityData, isTimeDependent, tempTimeIndex);
        [Z, xVec, yVec, xlabelText, ylabelText] = pitchMapField(fieldData, speciesData);
        titleText = pitchTitleText(speciesLabel, quantityLatex, timeIndexText);
        plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
            ['$' quantityLatex '$'], opt.colormapIndex, values.contourCount, ...
            sprintf('%s pitch %s%s, contours = %d', speciesLabel, quantityField, timeIndexText, values.contourCount));
    end
end

function plotPitchPower(picPitchData, opt)

    [speciesData, speciesLabel] = resolveSpeciesData(picPitchData, opt.species, 'pitch');
    powerData = requireQuantityData(speciesData, 'Power');
    [modeIdx, modeN] = parseModeN(opt.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    timeIndex = parseTimeIndex(opt.timeIndex, size(powerData, 4));
    [Z, xVec, yVec, xlabelText, ylabelText] = pitchMapField(powerData(:, :, modeIdx, timeIndex), speciesData);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s PitchPower 没有有限值数据。\n', speciesLabel);
        return;
    end

    quantityLatex = 'P_{\mathrm{PIC}}';
    timeIndexText = sprintf(', timeIndex = %d', timeIndex);
    titleText = pitchPowerTitleText(speciesLabel, quantityLatex, modeN, timeIndexText);
    plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantityLatex '$'], opt.colormapIndex, opt.contourCount, ...
        sprintf('%s PitchPower%s%s', speciesLabel, modeStatusText(modeN), timeIndexText));
    drawMapFigure(sprintf('%s PitchPower pitch-space map', speciesLabel), plotData);
end

function plotPitchPowerInteractive(picPitchData, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPitchData, opt.species, 'pitch');
    powerData = requireQuantityData(speciesData, 'Power');
    [~, modeN0] = parseModeN(opt.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    timeIndex0 = parseTimeIndex(opt.timeIndex, size(powerData, 4));
    contourCount = optionOrDefault(opt, 'contourCount', 0);
    controls = struct([]);
    if numel(speciesData.modeIndexAll) > 1
        modeControl = integerSliderControl('modeN', 'n', modeN0, ...
            min(speciesData.physicalNAll), max(speciesData.physicalNAll));
        modeControl.allowedValues = speciesData.physicalNAll;
        controls = [controls, modeControl];
    end
    if size(powerData, 4) > 1
        controls = [controls, integerSliderControl('timeIndex', 'timeIndex', timeIndex0, 1, size(powerData, 4))];
    end
    controls = [controls, integerSliderControl('contourCount', 'contours', contourCount, 0, max(20, 2 * contourCount))];

    plotInteractiveMap(sprintf('%s PitchPower pitch-space map', speciesLabel), controls, dynamicUpdate, ...
        @computePlotData, @renderMap);

    function plotData = computePlotData(values)
        tempModeN = modeN0;
        if isfield(values, 'modeN')
            tempModeN = values.modeN;
        end
        [tempModeIdx, tempModeN] = parseModeN(tempModeN, speciesData.modeIndexAll, speciesData.physicalNAll);
        tempTimeIndex = timeIndex0;
        if isfield(values, 'timeIndex')
            tempTimeIndex = values.timeIndex;
        end
        [Z, xVec, yVec, xlabelText, ylabelText] = pitchMapField(powerData(:, :, tempModeIdx, tempTimeIndex), speciesData);
        quantityLatex = 'P_{\mathrm{PIC}}';
        timeIndexText = sprintf(', timeIndex = %d', tempTimeIndex);
        titleText = pitchPowerTitleText(speciesLabel, quantityLatex, tempModeN, timeIndexText);
        plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
            ['$' quantityLatex '$'], opt.colormapIndex, values.contourCount, ...
            sprintf('%s PitchPower%s%s, contours = %d', speciesLabel, ...
            modeStatusText(tempModeN), timeIndexText, values.contourCount));
    end
end

function plotDiffusivityRadial(picDiffusivityData, opt)

    [speciesData, speciesLabel] = resolveSpeciesData(picDiffusivityData, opt.species, 'diffusivity');
    [Dsum, selectedN] = diffusivityModeSum(speciesData, opt.nRange);
    timeIndex = parseTimeIndex(opt.timeIndex, size(Dsum, 1));
    [xVec, xlabelText] = diffusivityRadialAxis(speciesData, opt.radialAxis);
    yVec = Dsum(timeIndex, :);

    if ~any(isfinite(yVec(:)))
        fprintf('[plot] %s Diffusivity 径向剖面没有有限值数据。\n', speciesLabel);
        return;
    end

    titleText = diffusivityTitleText(speciesLabel, selectedN, sprintf(',\\quad \\mathrm{timeIndex} = %d', timeIndex));
    statusText = sprintf('%s Diffusivity radial, n = [%d, %d], timeIndex = %d', ...
        speciesLabel, min(selectedN), max(selectedN), timeIndex);
    plotData = linePlotData(xVec, yVec, xlabelText, '$D/(\mathrm{m}^{2}/\mathrm{s})$', titleText, statusText);
    drawLineFigure(sprintf('%s Diffusivity radial', speciesLabel), plotData);
end

function plotDiffusivityRadialInteractive(picDiffusivityData, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picDiffusivityData, opt.species, 'diffusivity');
    D = requireQuantityData(speciesData, 'Diffusivity');
    nRange0 = initialDiffusivityNRange(speciesData, opt.nRange);
    timeIndex0 = parseTimeIndex(opt.timeIndex, size(D, 1));
    controls = diffusivityRangeControls(speciesData, nRange0);
    if size(D, 1) > 1
        controls = [controls, integerSliderControl('timeIndex', 'timeIndex', timeIndex0, 1, size(D, 1))];
    end

    plotInteractiveMap(sprintf('%s Diffusivity radial', speciesLabel), controls, dynamicUpdate, ...
        @computePlotData, @renderLine);

    function plotData = computePlotData(values)
        tempOpt = diffusivityOptionsFromValues(opt, values);
        [Dsum, selectedN] = diffusivityModeSum(speciesData, tempOpt.nRange);
        tempTimeIndex = parseTimeIndex(tempOpt.timeIndex, size(Dsum, 1));
        [xVec, xlabelText] = diffusivityRadialAxis(speciesData, tempOpt.radialAxis);
        titleText = diffusivityTitleText(speciesLabel, selectedN, ...
            sprintf(',\\quad \\mathrm{timeIndex} = %d', tempTimeIndex));
        statusText = sprintf('%s Diffusivity radial, n = [%d, %d], timeIndex = %d', ...
            speciesLabel, min(selectedN), max(selectedN), tempTimeIndex);
        plotData = linePlotData(xVec, Dsum(tempTimeIndex, :), xlabelText, ...
            '$D/(\mathrm{m}^{2}/\mathrm{s})$', titleText, statusText);
    end
end

function plotDiffusivityTime(picDiffusivityData, meta, opt)

    [speciesData, speciesLabel] = resolveSpeciesData(picDiffusivityData, opt.species, 'diffusivity');
    [Dsum, selectedN] = diffusivityModeSum(speciesData, opt.nRange);
    radialIndex = parseDiffusivityIndex(opt.radialIndex, size(Dsum, 2), 'radial');
    [xVec, xlabelText] = diffusivityTimeAxis(speciesData, meta, opt.timeAxis);
    [radialVec, ~, radialName] = diffusivityRadialAxis(speciesData, opt.radialAxis);
    radialValue = radialVec(radialIndex);
    yVec = Dsum(:, radialIndex);

    if ~any(isfinite(yVec(:)))
        fprintf('[plot] %s Diffusivity 时间曲线没有有限值数据。\n', speciesLabel);
        return;
    end

    titleText = diffusivityTitleText(speciesLabel, selectedN, ...
        sprintf(',\\quad %s = %.6g', radialName, radialValue));
    statusText = sprintf('%s Diffusivity time, n = [%d, %d], radialIndex = %d', ...
        speciesLabel, min(selectedN), max(selectedN), radialIndex);
    plotData = linePlotData(xVec, yVec, xlabelText, '$D/(\mathrm{m}^{2}/\mathrm{s})$', titleText, statusText);
    drawLineFigure(sprintf('%s Diffusivity time', speciesLabel), plotData);
end

function plotDiffusivityTimeInteractive(picDiffusivityData, meta, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picDiffusivityData, opt.species, 'diffusivity');
    D = requireQuantityData(speciesData, 'Diffusivity');
    nRange0 = initialDiffusivityNRange(speciesData, opt.nRange);
    radialIndex0 = parseDiffusivityIndex(opt.radialIndex, size(D, 3), 'radial');
    controls = diffusivityRangeControls(speciesData, nRange0);
    if size(D, 3) > 1
        controls = [controls, integerSliderControl('radialIndex', 'radialIndex', radialIndex0, 1, size(D, 3))];
    end

    plotInteractiveMap(sprintf('%s Diffusivity time', speciesLabel), controls, dynamicUpdate, ...
        @computePlotData, @renderLine);

    function plotData = computePlotData(values)
        tempOpt = diffusivityOptionsFromValues(opt, values);
        [Dsum, selectedN] = diffusivityModeSum(speciesData, tempOpt.nRange);
        tempRadialIndex = parseDiffusivityIndex(tempOpt.radialIndex, size(Dsum, 2), 'radial');
        [xVec, xlabelText] = diffusivityTimeAxis(speciesData, meta, tempOpt.timeAxis);
        [radialVec, ~, radialName] = diffusivityRadialAxis(speciesData, tempOpt.radialAxis);
        titleText = diffusivityTitleText(speciesLabel, selectedN, ...
            sprintf(',\\quad %s = %.6g', radialName, radialVec(tempRadialIndex)));
        statusText = sprintf('%s Diffusivity time, n = [%d, %d], radialIndex = %d', ...
            speciesLabel, min(selectedN), max(selectedN), tempRadialIndex);
        plotData = linePlotData(xVec, Dsum(:, tempRadialIndex), xlabelText, ...
            '$D/(\mathrm{m}^{2}/\mathrm{s})$', titleText, statusText);
    end
end

function plotDiffusivityMap(picDiffusivityData, meta, opt)

    [speciesData, speciesLabel] = resolveSpeciesData(picDiffusivityData, opt.species, 'diffusivity');
    [Dn, selectedN] = diffusivityModeSlice(speciesData, opt.nRange);
    [xVec, xlabelText] = diffusivityTimeAxis(speciesData, meta, opt.timeAxis);
    [yVec, ylabelText] = diffusivityRadialAxis(speciesData, opt.radialAxis);
    Z = Dn.';

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s Diffusivity 二维图没有有限值数据。\n', speciesLabel);
        return;
    end

    titleText = diffusivityTitleText(speciesLabel, selectedN, '');
    statusText = sprintf('%s Diffusivity map, n = %d', speciesLabel, selectedN);
    plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        '$D/(\mathrm{m}^{2}/\mathrm{s})$', opt.colormapIndex, 0, statusText);
    drawMapFigure(sprintf('%s Diffusivity map', speciesLabel), plotData);
end

function plotDiffusivityMapInteractive(picDiffusivityData, meta, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picDiffusivityData, opt.species, 'diffusivity');
    [nSliderMin, nSliderMax, n0, nAllowed] = diffusivitySingleNControlValues(speciesData, opt.nRange);
    controls = struct([]);
    if nSliderMin < nSliderMax
        controls = integerSliderControl('nSingle', 'n', n0, nSliderMin, nSliderMax);
        controls.allowedValues = nAllowed;
    end

    plotInteractiveMap(sprintf('%s Diffusivity map', speciesLabel), controls, dynamicUpdate, ...
        @computePlotData, @renderMap);

    function plotData = computePlotData(values)
        tempOpt = opt;
        tempOpt.nRange = n0;
        if isfield(values, 'nSingle')
            tempOpt.nRange = values.nSingle;
        end
        [Dn, selectedN] = diffusivityModeSlice(speciesData, tempOpt.nRange);
        [xVec, xlabelText] = diffusivityTimeAxis(speciesData, meta, tempOpt.timeAxis);
        [yVec, ylabelText] = diffusivityRadialAxis(speciesData, tempOpt.radialAxis);
        titleText = diffusivityTitleText(speciesLabel, selectedN, '');
        statusText = sprintf('%s Diffusivity map, n = %d', speciesLabel, selectedN);
        plotData = mapPlotData(Dn.', xVec, yVec, xlabelText, ylabelText, titleText, ...
            '$D/(\mathrm{m}^{2}/\mathrm{s})$', opt.colormapIndex, 0, statusText);
    end
end

function plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, idx, resOpt, meta)

    plotData.resonanceZ = [];
    plotData.resonanceOverlays = struct('Z', {}, 'label', {}, 'harmonic', {});
    if ~isfield(resOpt, 'enabled') || ~resOpt.enabled
        return;
    end

    if ~hasSpeciesData(phaseSpaceOrbit, speciesLabel)
        plotData.status = [plotData.status ', resonance = no orbit data'];
        return;
    end

    orbitData = phaseSpaceOrbit.(speciesLabel);
    resOpt = validatedResonanceOverlayOptions(resOpt);
    branchNames = resonanceBranchNames(resOpt);
    harmonicValues = resonanceHarmonicValues(resOpt);
    hasZero = false(1, numel(branchNames) * numel(harmonicValues));
    physicalN = resOpt.toroidalMode;
    iOverlay = 0;
    for iBranch = 1:numel(branchNames)
        for iHarmonic = 1:numel(harmonicValues)
            iOverlay = iOverlay + 1;
            lineOpt = resOpt;
            lineOpt.branch = branchNames{iBranch};
            lineOpt.harmonic = harmonicValues(iHarmonic);
            [residualHz, physicalN] = calculateResonanceResidualField(orbitData, lineOpt, meta);
            [resZ, ~, ~, ~, ~, ~] = slicePhaseField(dim, idx, residualHz, orbitData);
            plotData.resonanceOverlays(iOverlay).Z = resZ;
            plotData.resonanceOverlays(iOverlay).label = resonanceOverlayLabel( ...
                branchNames{iBranch}, harmonicValues(iHarmonic), numel(branchNames));
            plotData.resonanceOverlays(iOverlay).harmonic = harmonicValues(iHarmonic);
            hasZero(iOverlay) = residualHasZeroContour(resZ);
        end
    end
    plotData.status = [plotData.status ', ' resonanceOverlayStatusText('', resOpt, physicalN, branchNames, hasZero)];
end

function [quantityData, quantityField, quantityLatex, isTimeDependent] = resolvePhaseQuantity(speciesData, quantityName)

    quantityText = strtrim(char(quantityName));
    key = lower(strrep(strrep(quantityText, '_', ''), ' ', ''));

    if strcmp(quantityText, 'J') || ismember(key, {'j', 'jacobian', 'phasespacejacobian', 'pitchspacejacobian'})
        quantityField = 'J';
        quantityLatex = '\mathcal{J}';
        isTimeDependent = false;
        quantityData = requireQuantityData(speciesData, 'J');
    elseif strcmp(quantityText, 'F0') || ismember(key, {'phasespacef0', 'pitchspacef0'})
        quantityField = 'F0';
        quantityLatex = 'F_0';
        isTimeDependent = false;
        quantityData = requireQuantityData(speciesData, 'F0');
    elseif strcmp(quantityText, 'DF') || ismember(key, {'deltaf', 'phasedeltaf', 'pitchdeltaf'})
        quantityField = 'DF';
        quantityLatex = '\delta F';
        isTimeDependent = true;
        quantityData = requireQuantityData(speciesData, 'DF');
    elseif strcmp(quantityText, 'f0')
        quantityField = 'f0';
        quantityLatex = 'f_0';
        isTimeDependent = false;
        quantityData = divideQuantity(requireQuantityData(speciesData, 'F0'), ...
            requireQuantityData(speciesData, 'J'), 'J', 1e-6, 'absolute');
    elseif strcmp(quantityText, 'df')
        quantityField = 'df';
        quantityLatex = '\delta f';
        isTimeDependent = true;
        quantityData = divideQuantity(requireQuantityData(speciesData, 'DF'), ...
            requireQuantityData(speciesData, 'J'), 'J', 1e-6, 'absolute');
    elseif ismember(key, {'df/f0', 'dff0', 'df2f0', 'dfoverf0'})
        quantityField = 'df/f0';
        quantityLatex = '\delta f / f_0';
        isTimeDependent = true;
        quantityData = divideQuantity(requireQuantityData(speciesData, 'DF'), ...
            requireQuantityData(speciesData, 'F0'), 'F0', 1e-6, 'relative');
    else
        error('quantity 必须为 "J"、"F0"、"DF"、"f0"、"df" 或 "df/f0"。');
    end
end

function data = requireQuantityData(speciesData, quantityField)

    data = speciesData.(quantityField);
    if isempty(data)
        error('%s %s 尚未读取。请检查文件路径。', speciesData.species, quantityField);
    end
end

function data = divideQuantity(numerator, denominator, denominatorName, floorValue, floorMode)

    denominatorSize = size(denominator);
    numeratorSize = size(numerator);
    assert(numel(numeratorSize) >= numel(denominatorSize) && ...
        all(numeratorSize(1:numel(denominatorSize)) == denominatorSize), ...
        '分子和分母的网格尺寸不一致。');

    finiteDenominator = denominator(isfinite(denominator));
    if strcmpi(floorMode, 'absolute')
        threshold = floorValue;
    else
        threshold = 0;
        if ~isempty(finiteDenominator)
            threshold = floorValue * max(abs(finiteDenominator));
        end
    end

    invalidDenominator = ~isfinite(denominator) | abs(denominator) < threshold;
    data = bsxfun(@rdivide, numerator, denominator);
    if ndims(data) > ndims(invalidDenominator)
        repeatSize = ones(1, ndims(data));
        for dimIndex = ndims(invalidDenominator) + 1:ndims(data)
            repeatSize(dimIndex) = size(data, dimIndex);
        end
        invalidDenominator = repmat(invalidDenominator, repeatSize);
    end
    data(invalidDenominator) = NaN;

    if any(invalidDenominator(:))
        fprintf('[note] %s 中接近零或非有限的分母已置为 NaN。\n', denominatorName);
    end
end

function [fieldData, timeIndexText] = selectQuantityFrame(quantityData, isTimeDependent, timeIndex)

    if ~isTimeDependent
        fieldData = quantityData;
        timeIndexText = '';
        return;
    end

    nTime = timeDimensionSize(quantityData);
    idx = parseTimeIndex(timeIndex, nTime);
    if ndims(quantityData) == 4
        fieldData = quantityData(:, :, :, idx);
    elseif ndims(quantityData) == 3
        fieldData = quantityData(:, :, idx);
    else
        error('含时间的诊断量必须为 3D 或 4D 数组。');
    end
    timeIndexText = sprintf(', timeIndex = %d', idx);
end

function nTime = timeDimensionSize(quantityData)

    if ndims(quantityData) == 4
        nTime = size(quantityData, 4);
    elseif ndims(quantityData) == 3
        nTime = size(quantityData, 3);
    else
        error('含时间的诊断量必须为 3D 或 4D 数组。');
    end
end

function frequency = calculateOrbitFrequencyField(phaseSpaceData, branchName, unitScale, meta)

    dT = phaseSpaceData.(branchName).dT;
    trappedDT = phaseSpaceData.trapped.dT;
    fillFromTrapped = (dT == 0 | ~isfinite(dT)) & isfinite(trappedDT) & trappedDT > 0;
    dT(fillFromTrapped) = trappedDT(fillFromTrapped);
    dT(~isfinite(dT) | dT <= 0) = NaN;
    frequency = unitScale ./ (dT .* (meta.L0 / meta.VA0));
    frequency(~isfinite(frequency)) = NaN;
end

function [residualHz, physicalToroidalMode] = calculateResonanceResidualField(phaseSpaceData, resOpt, meta)

    branchName = normalizeResonanceBranch(resOpt.branch);
    physicalToroidalMode = resOpt.toroidalMode;
    classData = phaseSpaceData.(branchName);
    orbitTime = classData.dT .* (meta.L0 / meta.VA0);
    validOrbit = isfinite(orbitTime) & orbitTime > 0;
    phaseRate = nan(size(orbitTime));

    switch branchName
        case {'para', 'anti'}
            phaseAdvance = physicalToroidalMode .* classData.dphiTotal - ...
                resOpt.poloidalMode .* classData.dtheta + 2 * pi * resOpt.harmonic;
            phaseRate(validOrbit) = phaseAdvance(validOrbit) ./ orbitTime(validOrbit);
        case 'trapped'
            precessionPhase = classData.dphiTotal - classData.dphiVpara;
            precessionFrequency = nan(size(orbitTime));
            bounceFrequency = nan(size(orbitTime));
            precessionFrequency(validOrbit) = precessionPhase(validOrbit) ./ orbitTime(validOrbit);
            bounceFrequency(validOrbit) = 2 * pi ./ orbitTime(validOrbit);
            phaseRate(validOrbit) = physicalToroidalMode .* precessionFrequency(validOrbit) + ...
                resOpt.harmonic .* bounceFrequency(validOrbit);
    end

    residualHz = (2 * pi * resOpt.frequencyHz - phaseRate) ./ (2 * pi);
    residualHz(~isfinite(residualHz)) = NaN;
end

function [Z, xVec, yVec, xlabelText, ylabelText, titleText] = slicePhaseField(dim, idx, fieldData, speciesData)

    switch dim
        case 1
            idx = clampIndex(idx, numel(speciesData.E1d), 'E');
            Z = reshape(fieldData(idx, :, :), numel(speciesData.Pphi1d), numel(speciesData.Lambda1d));
            xVec = speciesData.Lambda1d;
            yVec = -speciesData.Pphi1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$P_{\varphi}$';
            titleText = sprintf('$E = %.6g$', speciesData.E1d(idx));
        case 2
            idx = clampIndex(idx, numel(speciesData.Pphi1d), 'Pphi');
            Z = reshape(fieldData(:, idx, :), numel(speciesData.E1d), numel(speciesData.Lambda1d));
            xVec = speciesData.Lambda1d;
            yVec = speciesData.E1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$E$';
            titleText = sprintf('$P_{\\varphi} = %.6g$', -speciesData.Pphi1d(idx));
        case 3
            idx = clampIndex(idx, numel(speciesData.Lambda1d), 'Lambda');
            Z = reshape(fieldData(:, :, idx), numel(speciesData.E1d), numel(speciesData.Pphi1d));
            xVec = -speciesData.Pphi1d;
            yVec = speciesData.E1d;
            xlabelText = '$P_{\varphi}$';
            ylabelText = '$E$';
            titleText = sprintf('$\\Lambda = %.6g$', speciesData.Lambda1d(idx));
        otherwise
            error('切片维度必须为 1、2 或 3。');
    end

    Z(~isfinite(Z)) = NaN;
end

function [Z, xVec, yVec, xlabelText, ylabelText] = pitchMapField(fieldData, speciesData)

    assert(isequal(size(fieldData), [numel(speciesData.Vpara1d), numel(speciesData.Vperp1d)]), ...
        'pitch 空间数据尺寸必须为 [gridVpara, gridVperp]。');
    Z = fieldData.';
    Z(~isfinite(Z)) = NaN;
    xVec = speciesData.Vpara1d;
    yVec = speciesData.Vperp1d;
    xlabelText = '$v_{\parallel}$';
    ylabelText = '$v_{\perp}$';
end

function [Dsum, selectedN] = diffusivityModeSum(speciesData, nRange)

    data = requireQuantityData(speciesData, 'Diffusivity');
    assert(ndims(data) == 3 && size(data, 1) == numel(speciesData.tDiag) && ...
        size(data, 2) == numel(speciesData.modeIndexAll) && size(data, 3) == speciesData.gridNx, ...
        'Diffusivity 数据尺寸必须为 [time, modeN, radial]。');

    [modeMask, selectedN] = diffusivityModeMask(speciesData, nRange);
    selectedData = data(:, modeMask, :);
    validData = isfinite(selectedData);
    selectedData(~validData) = 0;
    Dsum = sum(selectedData, 2);
    validCount = sum(validData, 2);
    Dsum(validCount == 0) = NaN;
    Dsum = reshape(Dsum, size(data, 1), size(data, 3));
end

function [Dn, selectedN] = diffusivityModeSlice(speciesData, nRange)

    data = requireQuantityData(speciesData, 'Diffusivity');
    physicalN = diffusivityPhysicalN(speciesData);
    nTarget = diffusivitySingleN(speciesData, nRange);
    [~, modeIdx] = min(abs(physicalN - nTarget));
    selectedN = physicalN(modeIdx);
    Dn = reshape(data(:, modeIdx, :), size(data, 1), size(data, 3));
    Dn(~isfinite(Dn)) = NaN;
end

function [modeMask, selectedN] = diffusivityModeMask(speciesData, nRange)

    physicalN = diffusivityPhysicalN(speciesData);
    range = normalizeDiffusivityNRange(nRange);
    modeMask = physicalN >= range(1) & physicalN <= range(2);
    if ~any(modeMask)
        [~, nearestIdx] = min(abs(physicalN - mean(range)));
        modeMask(nearestIdx) = true;
    end
    selectedN = physicalN(modeMask);
end

function physicalN = diffusivityPhysicalN(speciesData)

    physicalN = speciesData.physicalNAll;
end

function range = normalizeDiffusivityNRange(nRange)

    range = reshape(double(nRange), 1, []);
    assert(numel(range) == 1 || numel(range) == 2, 'nRange 必须为标量或 [min max]。');
    assert(all(isfinite(range)), 'nRange 必须为有限数。');
    if isscalar(range)
        range = [range, range];
    else
        range = sort(range);
    end
end

function range = initialDiffusivityNRange(speciesData, nRange)

    [~, selectedN] = diffusivityModeMask(speciesData, nRange);
    range = [min(selectedN), max(selectedN)];
end

function nValue = diffusivitySingleN(speciesData, nRange)

    physicalN = diffusivityPhysicalN(speciesData);
    rawRange = normalizeDiffusivityNRange(nRange);
    inRangeN = physicalN(physicalN >= rawRange(1) & physicalN <= rawRange(2));
    if isempty(inRangeN)
        [~, nearestIdx] = min(abs(physicalN - mean(rawRange)));
        nValue = physicalN(nearestIdx);
    else
        nValue = inRangeN(max(1, round(numel(inRangeN) / 2)));
    end
end

function [nSliderMin, nSliderMax, n0, nAllowed] = diffusivitySingleNControlValues(speciesData, nRange)

    [~, selectedN] = diffusivityModeMask(speciesData, nRange);
    nSliderMin = min(selectedN);
    nSliderMax = max(selectedN);
    n0 = diffusivitySingleN(speciesData, nRange);
    nAllowed = selectedN;
end

function controls = diffusivityRangeControls(speciesData, nRange)

    range = initialDiffusivityNRange(speciesData, nRange);
    if range(1) == range(2)
        controls = struct([]);
        return;
    end

    controls = [ ...
        integerSliderControl('nMin', 'n min', range(1), range(1), range(2)), ...
        integerSliderControl('nMax', 'n max', range(2), range(1), range(2))];
    nAllowed = speciesData.physicalNAll(speciesData.physicalNAll >= range(1) & speciesData.physicalNAll <= range(2));
    controls(1).allowedValues = nAllowed;
    controls(2).allowedValues = nAllowed;
end

function opt = diffusivityOptionsFromValues(opt, values)

    if isfield(values, 'nMin') && isfield(values, 'nMax')
        opt.nRange = sort([values.nMin, values.nMax]);
    end
    if isfield(values, 'timeIndex')
        opt.timeIndex = values.timeIndex;
    end
    if isfield(values, 'radialIndex')
        opt.radialIndex = values.radialIndex;
    end
end

function [xVec, xlabelText] = diffusivityTimeAxis(speciesData, meta, axisType)

    key = lower(strtrim(char(axisType)));
    switch key
        case 'ta'
            xVec = speciesData.tDiag;
            xlabelText = '$t_a$';
        case 'ms'
            xVec = speciesData.timeSeconds * 1e3;
            xlabelText = '$t/\mathrm{ms}$';
        case 's'
            xVec = speciesData.timeSeconds;
            xlabelText = '$t/\mathrm{s}$';
        case 'steps'
            xVec = (0:numel(speciesData.tDiag) - 1) * meta.diagSteps;
            xlabelText = '$\mathrm{step}$';
        otherwise
            error('timeAxis 必须为 "ta"、"ms"、"s" 或 "steps"。');
    end
end

function [xVec, xlabelText, axisName] = diffusivityRadialAxis(speciesData, axisType)

    key = lower(strtrim(char(axisType)));
    switch key
        case 'rho'
            xVec = speciesData.rhoGrid;
            xlabelText = '$\rho$';
            axisName = '\rho';
        case 'x'
            xVec = speciesData.xGrid;
            xlabelText = '$x$';
            axisName = 'x';
        otherwise
            error('radialAxis 必须为 "rho" 或 "x"。');
    end
end

function idx = parseDiffusivityIndex(indexText, nIndex, label)

    idx = parseIndex(indexText, nIndex, label);
end

function titleText = diffusivityTitleText(speciesLabel, selectedN, extraLatex)

    if numel(selectedN) == 1
        nLatex = sprintf('n = %d', selectedN);
    else
        nLatex = sprintf('n \\in [%d,%d]', min(selectedN), max(selectedN));
    end
    titleText = sprintf('$\\mathrm{%s}\\quad D,\\quad %s%s$', speciesLabel, nLatex, extraLatex);
end

function plotData = linePlotData(xVec, yVec, xlabelText, ylabelText, titleText, statusText)

    xVec = reshape(xVec, 1, []);
    yVec = reshape(yVec, 1, []);
    assert(numel(xVec) == numel(yVec), '线图横纵坐标长度不一致。');
    plotData = struct( ...
        'xVec', xVec, ...
        'yVec', yVec, ...
        'xlabelText', xlabelText, ...
        'ylabelText', ylabelText, ...
        'titleText', titleText, ...
        'status', statusText);
end

function drawLineFigure(figureName, plotData)

    figHandle = figure('Name', figureName, 'Color', 'w', 'Position', [120, 120, 900, 560]);
    axHandle = axes('Parent', figHandle, 'Units', 'normalized', 'Position', [0.12, 0.14, 0.82, 0.76]);
    renderLine(axHandle, plotData);
end

function renderLine(axHandle, plotData)

    resetMapAxes(axHandle);
    validData = isfinite(plotData.xVec) & isfinite(plotData.yVec);
    if any(validData)
        plot(axHandle, plotData.xVec, plotData.yVec, 'LineWidth', 1.8, 'Color', [0.20, 0.40, 0.80]);
    else
        text(axHandle, 0.5, 0.5, 'no finite data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 14);
    end

    xlabel(axHandle, plotData.xlabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(axHandle, plotData.ylabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(axHandle, plotData.titleText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 13);
    grid(axHandle, 'on');
    axis(axHandle, 'tight');
    box(axHandle, 'on');
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14, ...
        'Color', 'white', 'Layer', 'top', 'TickDir', 'out', 'LineWidth', 1);
end

function drawMapFigure(figureName, plotData)

    figHandle = figure('Name', figureName, 'Color', 'w', 'Position', [100, 80, 980, 760]);
    axHandle = axes('Parent', figHandle, 'Units', 'normalized', 'Position', [0.11, 0.12, 0.74, 0.78]);
    renderMap(axHandle, plotData);
end

function renderMap(axHandle, plotData)

    resetMapAxes(axHandle);
    [X, Y] = meshgrid(plotData.xVec, plotData.yVec);
    pcolor(axHandle, X, Y, plotData.Z);
    shading(axHandle, 'interp');
    hold(axHandle, 'on');

    finiteZ = plotData.Z(isfinite(plotData.Z));
    if plotData.contourCount > 0 && ~isempty(finiteZ) && min(finiteZ) < max(finiteZ)
        contour(axHandle, X, Y, plotData.Z, plotData.contourCount, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end

    forceSigned = isfield(plotData, 'forceSigned') && plotData.forceSigned;
    isSigned = forceSigned || (~isempty(finiteZ) && any(finiteZ < 0));
    if isSigned
        [colorMin, colorMax] = symmetricFiniteColorLimits(plotData.Z);
    else
        [colorMin, colorMax] = finiteColorLimits(plotData.Z);
    end
    clim(axHandle, [colorMin, colorMax]);
    colormap(axHandle, phaseSpaceColormap(isSigned, plotData.colormapIndex, 256));

    cb = colorbar(axHandle);
    cb.FontName = 'Times New Roman';
    cb.FontSize = 13;
    ylabel(cb, plotData.colorbarLabel, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 13);

    if isfield(plotData, 'resonanceOverlays') && ~isempty(plotData.resonanceOverlays)
        drawResonanceOverlays(axHandle, X, Y, plotData.resonanceOverlays);
    elseif isfield(plotData, 'resonanceZ') && ~isempty(plotData.resonanceZ) && residualHasZeroContour(plotData.resonanceZ)
        contour(axHandle, X, Y, plotData.resonanceZ, [0, 0], 'EdgeColor', 'w', 'LineWidth', 2.6);
        contour(axHandle, X, Y, plotData.resonanceZ, [0, 0], 'EdgeColor', 'k', 'LineWidth', 1.2);
    end

    xlabel(axHandle, plotData.xlabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(axHandle, plotData.ylabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(axHandle, plotData.titleText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 13);
    axis(axHandle, 'tight');
    box(axHandle, 'on');
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14, ...
        'Color', 'white', 'Layer', 'top', 'TickDir', 'out', 'LineWidth', 1);
    hold(axHandle, 'off');
end

function drawResonanceOverlays(axHandle, X, Y, overlays)

    legendHandles = gobjects(0);
    legendLabels = {};
    visibleIndex = 0;
    for iOverlay = 1:numel(overlays)
        resZ = overlays(iOverlay).Z;
        if ~residualHasZeroContour(resZ)
            continue;
        end

        visibleIndex = visibleIndex + 1;
        lineColor = resonanceOverlayColor(visibleIndex, numel(overlays));
        contour(axHandle, X, Y, resZ, [0, 0], ...
            'EdgeColor', 'w', 'LineWidth', 2.8, 'HandleVisibility', 'off');
        [~, lineHandle] = contour(axHandle, X, Y, resZ, [0, 0], ...
            'EdgeColor', lineColor, 'LineWidth', 1.35);
        legendHandles(end + 1) = lineHandle;
        legendLabels{end + 1} = overlays(iOverlay).label;
    end

    if ~isempty(legendHandles)
        legend(axHandle, legendHandles, legendLabels, 'Interpreter', 'none', ...
            'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 10, 'Box', 'on');
    end
end

function lineColor = resonanceOverlayColor(index, totalCount)

    if totalCount <= 1
        lineColor = [0, 0, 0];
        return;
    end
    colorTable = lines(max(totalCount, 7));
    lineColor = colorTable(index, :);
end

function resetMapAxes(axHandle)

    figHandle = ancestor(axHandle, 'figure');
    delete(findall(figHandle, 'Type', 'ColorBar'));
    delete(findall(figHandle, 'Type', 'Legend'));
    cla(axHandle);
end

function plotData = mapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
    colorbarLabel, colormapIndex, contourCount, statusText)

    plotData = struct( ...
        'Z', Z, ...
        'xVec', xVec, ...
        'yVec', yVec, ...
        'xlabelText', xlabelText, ...
        'ylabelText', ylabelText, ...
        'titleText', titleText, ...
        'colorbarLabel', colorbarLabel, ...
        'colormapIndex', colormapIndex, ...
        'contourCount', contourCount, ...
        'status', statusText, ...
        'forceSigned', false, ...
        'resonanceZ', [], ...
        'resonanceOverlays', struct('Z', {}, 'label', {}, 'harmonic', {}));
end

function plotInteractiveMap(figName, controls, dynamicUpdate, computePlotData, renderPlotData)

    figHandle = figure('Name', figName, 'Color', 'w', 'Position', [80, 40, 1080, 860]);
    nControl = numel(controls);
    controlBottom = 0.030;
    controlSpacing = 0.035;
    statusBottom = controlBottom + controlSpacing * max(nControl, 1) + 0.012;
    axesBottom = statusBottom + 0.065;
    axHandle = axes('Parent', figHandle, 'Units', 'normalized', ...
        'Position', [0.10, axesBottom, 0.74, 0.94 - axesBottom]);
    statusText = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.10, statusBottom, 0.82, 0.035], 'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', 'FontName', 'Times New Roman', 'FontSize', 11);

    sliderLabels = gobjects(nControl, 1);
    sliders = gobjects(nControl, 1);
    for iControl = 1:nControl
        yPos = controlBottom + controlSpacing * (nControl - iControl);
        sliderLabels(iControl) = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.10, yPos - 0.007, 0.19, 0.030], ...
            'String', sliderLabelText(controls(iControl), controls(iControl).value), ...
            'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', 11);
        sliders(iControl) = uicontrol(figHandle, 'Style', 'slider', 'Units', 'normalized', ...
            'Position', [0.30, yPos, 0.62, 0.022], ...
            'Min', controls(iControl).min, 'Max', controls(iControl).max, ...
            'Value', controls(iControl).value, ...
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
            warning('visualizePIC:DynamicSliderUnavailable', ...
                '当前 MATLAB 环境不支持滑块连续更新，将在释放滑块后更新。');
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

        plotData = computePlotData(values);
        renderPlotData(axHandle, plotData);
        set(statusText, 'String', plotData.status);
    end
end

function controls = appendResonanceControls(controls, resOpt, useHarmonicBounds, allowMultipleBranches)

    if ~isfield(resOpt, 'enabled') || ~resOpt.enabled
        return;
    end
    if nargin < 3
        useHarmonicBounds = false;
    end
    if nargin < 4
        allowMultipleBranches = false;
    end

    if allowMultipleBranches
        resOpt = validatedResonanceOverlayOptions(resOpt);
        branchNames = resonanceBranchNames(resOpt);
    else
        resOpt = validatedResonanceOptions(resOpt);
        branchNames = {resOpt.branch};
    end
    freqRange = normalizeSliderRange(resOpt.frequencyHzRange, resOpt.frequencyHz, 'frequencyHzRange', false);
    nRange = normalizeSliderRange(resOpt.toroidalModeRange, resOpt.toroidalMode, 'toroidalModeRange', true);
    lBounds = [resOpt.harmonicMin, resOpt.harmonicMax];
    lRange = normalizeSliderRangeForValues(resOpt.harmonicRange, lBounds, 'harmonicRange', true);
    nRange(1) = max(0, nRange(1));
    if nRange(1) == nRange(2)
        nRange(2) = nRange(1) + 1;
    end

    controls = [controls, ...
        numericSliderControl('resFrequencyHz', 'f/Hz', resOpt.frequencyHz, freqRange(1), freqRange(2), 501), ...
        integerSliderControl('resToroidalMode', 'n', resOpt.toroidalMode, nRange(1), nRange(2))];
    if any(~strcmp(branchNames, 'trapped'))
        mRange = normalizeSliderRange(resOpt.poloidalModeRange, resOpt.poloidalMode, 'poloidalModeRange', true);
        mRange(1) = max(0, mRange(1));
        if mRange(1) == mRange(2)
            mRange(2) = mRange(1) + 1;
        end
        controls = [controls, integerSliderControl('resPoloidalMode', 'm', resOpt.poloidalMode, mRange(1), mRange(2))];
    end
    if useHarmonicBounds
        controls = [controls, ...
            integerSliderControl('resHarmonicMin', 'l min', lBounds(1), lRange(1), lRange(2)), ...
            integerSliderControl('resHarmonicMax', 'l max', lBounds(2), lRange(1), lRange(2))];
    else
        controls = [controls, integerSliderControl('resHarmonic', 'l', resOpt.harmonic, lRange(1), lRange(2))];
    end
end

function resOpt = resonanceOptionsFromValues(resOpt, values, allowMultipleBranches)

    if ~isfield(resOpt, 'enabled') || ~resOpt.enabled
        return;
    end
    if nargin < 3
        allowMultipleBranches = false;
    end
    if isfield(values, 'resFrequencyHz')
        resOpt.frequencyHz = values.resFrequencyHz;
    end
    if isfield(values, 'resToroidalMode')
        resOpt.toroidalMode = values.resToroidalMode;
    end
    if isfield(values, 'resPoloidalMode')
        resOpt.poloidalMode = values.resPoloidalMode;
    end
    if isfield(values, 'resHarmonic')
        resOpt.harmonic = values.resHarmonic;
        resOpt.harmonicMin = values.resHarmonic;
        resOpt.harmonicMax = values.resHarmonic;
    end
    if isfield(values, 'resHarmonicMin')
        resOpt.harmonicMin = values.resHarmonicMin;
    end
    if isfield(values, 'resHarmonicMax')
        resOpt.harmonicMax = values.resHarmonicMax;
    end
    if allowMultipleBranches
        resOpt = validatedResonanceOverlayOptions(resOpt);
    else
        resOpt = validatedResonanceOptions(resOpt);
    end
end

function resOpt = validatedResonanceOptions(resOpt)

    if ~isfield(resOpt, 'enabled')
        resOpt.enabled = true;
    end
    resOpt.branch = normalizeResonanceBranch(resOpt.branch);
    resOpt.frequencyHz = validateFiniteScalar(resOpt.frequencyHz, 'frequencyHz');
    resOpt.toroidalMode = readNonnegativeIntegerScalar(resOpt.toroidalMode, 'toroidalMode');
    if ~strcmp(resOpt.branch, 'trapped')
        resOpt.poloidalMode = readNonnegativeIntegerScalar(resOpt.poloidalMode, 'poloidalMode');
    end
    resOpt = normalizeResonanceHarmonicOptions(resOpt);
end

function resOpt = validatedResonanceOverlayOptions(resOpt)

    if ~isfield(resOpt, 'enabled')
        resOpt.enabled = true;
    end
    branchNames = normalizeResonanceBranches(resOpt.branch);
    resOpt.branch = branchNames;
    resOpt.frequencyHz = validateFiniteScalar(resOpt.frequencyHz, 'frequencyHz');
    resOpt.toroidalMode = readNonnegativeIntegerScalar(resOpt.toroidalMode, 'toroidalMode');
    if any(~strcmp(branchNames, 'trapped'))
        resOpt.poloidalMode = readNonnegativeIntegerScalar(resOpt.poloidalMode, 'poloidalMode');
    elseif ~isfield(resOpt, 'poloidalMode') || isempty(resOpt.poloidalMode)
        resOpt.poloidalMode = 0;
    end
    resOpt = normalizeResonanceHarmonicOptions(resOpt);
end

function resOpt = normalizeResonanceHarmonicOptions(resOpt)

    if ~isfield(resOpt, 'harmonic') || isempty(resOpt.harmonic)
        resOpt.harmonic = 0;
    end
    rawHarmonic = reshape(double(resOpt.harmonic), 1, []);
    assert(~isempty(rawHarmonic) && all(isfinite(rawHarmonic)) && all(rawHarmonic == floor(rawHarmonic)), ...
        'harmonic 必须是整数标量或整数范围。');

    if isfield(resOpt, 'harmonicMin') && ~isempty(resOpt.harmonicMin)
        harmonicMin = validateIntegerScalar(resOpt.harmonicMin, 'harmonicMin');
    else
        harmonicMin = min(rawHarmonic);
    end
    if isfield(resOpt, 'harmonicMax') && ~isempty(resOpt.harmonicMax)
        harmonicMax = validateIntegerScalar(resOpt.harmonicMax, 'harmonicMax');
    else
        harmonicMax = max(rawHarmonic);
    end
    if harmonicMin > harmonicMax
        tmp = harmonicMin;
        harmonicMin = harmonicMax;
        harmonicMax = tmp;
    end

    resOpt.harmonicMin = harmonicMin;
    resOpt.harmonicMax = harmonicMax;
    resOpt.harmonic = validateIntegerScalar(rawHarmonic(1), 'harmonic');
end

function harmonicValues = resonanceHarmonicValues(resOpt)

    harmonicValues = resOpt.harmonicMin:resOpt.harmonicMax;
end

function branchNames = resonanceBranchNames(resOpt)

    branchNames = normalizeResonanceBranches(resOpt.branch);
end

function label = resonanceOverlayLabel(branchName, harmonic, nBranch)

    if nBranch > 1
        label = sprintf('%s, l = %d', branchName, harmonic);
    else
        label = sprintf('l = %d', harmonic);
    end
end

function titleText = simpleTitleText(speciesLabel, quantityLatex, sliceTitle, extraLatex)

    sliceLatex = stripMathDelimiters(sliceTitle);
    titleText = sprintf('$\\mathrm{%s}\\quad %s,\\quad %s%s$', ...
        speciesLabel, quantityLatex, sliceLatex, extraLatex);
end

function text = orbitFrequencyLatex(branchName)

    text = sprintf('\\mathrm{%s}\\ \\mathrm{orbit}\\ \\mathrm{frequency}', branchName);
end

function titleText = phaseQuantityTitleText(speciesLabel, quantityLatex, sliceTitle, timeIndexText, resOpt)

    extraLatex = timeTextToLatex(timeIndexText);
    if isfield(resOpt, 'enabled') && resOpt.enabled
        extraLatex = [extraLatex ',\quad \mathrm{res}'];
    end
    titleText = simpleTitleText(speciesLabel, quantityLatex, sliceTitle, extraLatex);
end

function titleText = phasePowerTitleText(speciesLabel, quantityLatex, sliceTitle, modeN, timeIndexText, resOpt)

    extraLatex = sprintf(',\\quad n = %d%s', modeN, timeTextToLatex(timeIndexText));
    if isfield(resOpt, 'enabled') && resOpt.enabled
        extraLatex = [extraLatex ',\quad \mathrm{res}'];
    end
    titleText = simpleTitleText(speciesLabel, quantityLatex, sliceTitle, extraLatex);
end

function titleText = pitchTitleText(speciesLabel, quantityLatex, timeIndexText)

    titleText = sprintf('$\\mathrm{%s}\\quad %s\\quad \\mathrm{pitch}%s$', ...
        speciesLabel, quantityLatex, timeTextToLatex(timeIndexText));
end

function titleText = pitchPowerTitleText(speciesLabel, quantityLatex, modeN, timeIndexText)

    titleText = sprintf('$\\mathrm{%s}\\quad %s\\quad \\mathrm{pitch},\\quad n = %d%s$', ...
        speciesLabel, quantityLatex, modeN, timeTextToLatex(timeIndexText));
end

function titleText = resonanceTitleText(speciesLabel, resOpt, physicalN, sliceTitle)

    sliceLatex = stripMathDelimiters(sliceTitle);
    resLatex = sprintf('\\mathrm{%s}\\ \\mathrm{resonance}', resOpt.branch);
    if strcmp(resOpt.branch, 'trapped')
        titleText = sprintf('$\\mathrm{%s}\\quad %s,\\quad %s,\\quad n = %d,\\quad l = %d,\\quad f = %.6g\\,\\mathrm{Hz}$', ...
            speciesLabel, resLatex, sliceLatex, physicalN, resOpt.harmonic, resOpt.frequencyHz);
    else
        titleText = sprintf('$\\mathrm{%s}\\quad %s,\\quad %s,\\quad n = %d,\\quad m = %d,\\quad l = %d,\\quad f = %.6g\\,\\mathrm{Hz}$', ...
            speciesLabel, resLatex, sliceLatex, physicalN, resOpt.poloidalMode, resOpt.harmonic, resOpt.frequencyHz);
    end
end

function text = resonanceStatusText(speciesLabel, resOpt, physicalN, hasZero)

    zeroText = resonanceZeroText(hasZero);
    harmonicText = resonanceHarmonicStatusText(resOpt);
    prefix = '';
    if ~isempty(speciesLabel)
        prefix = [speciesLabel ' '];
    end
    if strcmp(resOpt.branch, 'trapped')
        text = sprintf('%sres %s, f = %.6g Hz, n = %d, %s, zero = %s', ...
            prefix, resOpt.branch, resOpt.frequencyHz, physicalN, harmonicText, zeroText);
    else
        text = sprintf('%sres %s, f = %.6g Hz, n = %d, m = %d, %s, zero = %s', ...
            prefix, resOpt.branch, resOpt.frequencyHz, physicalN, resOpt.poloidalMode, harmonicText, zeroText);
    end
end

function text = resonanceOverlayStatusText(speciesLabel, resOpt, physicalN, branchNames, hasZero)

    zeroText = resonanceZeroText(hasZero);
    harmonicText = resonanceHarmonicStatusText(resOpt);
    branchText = strjoin(branchNames, '/');
    prefix = '';
    if ~isempty(speciesLabel)
        prefix = [speciesLabel ' '];
    end
    if any(~strcmp(branchNames, 'trapped'))
        text = sprintf('%sres %s, f = %.6g Hz, n = %d, m = %d, %s, zero = %s', ...
            prefix, branchText, resOpt.frequencyHz, physicalN, resOpt.poloidalMode, harmonicText, zeroText);
    else
        text = sprintf('%sres %s, f = %.6g Hz, n = %d, %s, zero = %s', ...
            prefix, branchText, resOpt.frequencyHz, physicalN, harmonicText, zeroText);
    end
end

function text = resonanceZeroText(hasZero)

    if numel(hasZero) > 1
        text = sprintf('%d/%d', nnz(hasZero), numel(hasZero));
    elseif hasZero
        text = 'yes';
    else
        text = 'no';
    end
end

function text = resonanceHarmonicStatusText(resOpt)

    if isfield(resOpt, 'harmonicMin') && isfield(resOpt, 'harmonicMax') && ...
            resOpt.harmonicMin ~= resOpt.harmonicMax
        text = sprintf('l = [%d, %d]', resOpt.harmonicMin, resOpt.harmonicMax);
    elseif isfield(resOpt, 'harmonicMin')
        text = sprintf('l = %d', resOpt.harmonicMin);
    else
        text = sprintf('l = %d', resOpt.harmonic);
    end
end

function text = timeTextToLatex(timeIndexText)

    text = '';
    if ~isempty(timeIndexText)
        text = strrep(timeIndexText, ', timeIndex = ', ',\quad \mathrm{timeIndex} = ');
    end
end

function text = stripMathDelimiters(text)

    if startsWith(text, '$') && endsWith(text, '$') && strlength(text) >= 2
        text = extractBetween(string(text), 2, strlength(text) - 1);
        text = char(text);
    end
end

function [speciesData, speciesLabel] = resolveSpeciesData(dataStruct, speciesName, dataLabel)

    speciesLabel = strtrim(char(speciesName));
    if isfield(dataStruct, speciesLabel) && ~isempty(dataStruct.(speciesLabel))
        speciesData = dataStruct.(speciesLabel);
        return;
    end

    fields = fieldnames(dataStruct);
    for fieldIndex = 1:numel(fields)
        if strcmpi(fields{fieldIndex}, speciesLabel) && ~isempty(dataStruct.(fields{fieldIndex}))
            speciesLabel = fields{fieldIndex};
            speciesData = dataStruct.(speciesLabel);
            return;
        end
    end

    error('没有找到物种 "%s" 的 %s 数据。', speciesLabel, dataLabel);
end

function ok = hasSpeciesData(dataStruct, speciesName)

    speciesName = strtrim(char(speciesName));
    ok = isfield(dataStruct, speciesName) && ~isempty(dataStruct.(speciesName));
end

function ok = hasSpeciesFieldData(dataStruct, speciesName, fieldName)

    ok = false;
    speciesName = strtrim(char(speciesName));
    if isfield(dataStruct, speciesName) && isfield(dataStruct.(speciesName), fieldName)
        ok = ~isempty(dataStruct.(speciesName).(fieldName));
    end
end

function ok = hasRequestedPhaseQuantity(dataStruct, speciesName, quantityName)

    speciesName = strtrim(char(speciesName));
    if ~isfield(dataStruct, speciesName)
        ok = false;
        return;
    end

    quantityText = strtrim(char(quantityName));
    key = lower(strrep(strrep(quantityText, '_', ''), ' ', ''));
    if strcmp(quantityText, 'J') || ismember(key, {'j', 'jacobian', 'phasespacejacobian', 'pitchspacejacobian'})
        requiredFields = {'J'};
    elseif strcmp(quantityText, 'F0') || ismember(key, {'phasespacef0', 'pitchspacef0'})
        requiredFields = {'F0'};
    elseif strcmp(quantityText, 'DF') || ismember(key, {'deltaf', 'phasedeltaf', 'pitchdeltaf'})
        requiredFields = {'DF'};
    elseif strcmp(quantityText, 'f0')
        requiredFields = {'F0', 'J'};
    elseif strcmp(quantityText, 'df')
        requiredFields = {'DF', 'J'};
    elseif ismember(key, {'df/f0', 'dff0', 'df2f0', 'dfoverf0'})
        requiredFields = {'DF', 'F0'};
    else
        requiredFields = {quantityText};
    end

    ok = true;
    for iField = 1:numel(requiredFields)
        ok = ok && isfield(dataStruct.(speciesName), requiredFields{iField}) && ...
            ~isempty(dataStruct.(speciesName).(requiredFields{iField}));
    end
end

function dim = phaseCoordinateToDimension(coordinateText)

    key = lower(strtrim(char(coordinateText)));
    key = strrep(key, '_', '');
    key = strrep(key, '{', '');
    key = strrep(key, '}', '');
    key = strrep(key, '\', '');

    switch key
        case {'e', 'energy'}
            dim = 1;
        case {'pphi', 'pvarphi'}
            dim = 2;
        case 'lambda'
            dim = 3;
        otherwise
            error('coordinate 必须为 "E"、"Pphi" 或 "Lambda"。');
    end
end

function idx = parsePhaseSlice(sliceText, dim, speciesData)

    if isnumeric(sliceText)
        idx = double(sliceText);
    else
        key = lower(strtrim(char(sliceText)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            idx = defaultSliceIndex(dim, speciesData.gridE, speciesData.gridPphi, speciesData.gridLambda);
        elseif ismember(key, {'end', 'last'})
            idx = phaseDimensionSize(dim, speciesData);
        else
            idx = str2double(key);
        end
    end

    idx = clampIndex(idx, phaseDimensionSize(dim, speciesData), char(phaseDimensionLabel(dim)));
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

function [modeIdx, modeN] = parseModeN(modeNText, modeIndexAll, physicalNAll)

    assert(~isempty(modeIndexAll) && numel(modeIndexAll) == numel(physicalNAll), ...
        'modeN 列表为空或尺寸不一致。');
    if isnumeric(modeNText)
        modeN = double(modeNText);
    else
        key = lower(strtrim(char(modeNText)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            modeIdx = max(1, round(numel(modeIndexAll) / 2));
            modeN = physicalNAll(modeIdx);
            return;
        elseif ismember(key, {'end', 'last'})
            modeIdx = numel(modeIndexAll);
            modeN = physicalNAll(modeIdx);
            return;
        elseif ismember(key, {'first', 'begin'})
            modeIdx = 1;
            modeN = physicalNAll(modeIdx);
            return;
        else
            modeN = str2double(key);
        end
    end

    assert(isscalar(modeN) && isfinite(modeN) && modeN == floor(modeN) && modeN >= 0, ...
        'modeN 必须为非负整数物理环向模数。');
    modeIdx = find(physicalNAll == modeN, 1);
    assert(~isempty(modeIdx), 'modeN 必须位于有效物理 n 集合 [%s]。当前值：%d。', ...
        formatNumberList(physicalNAll), modeN);
end

function text = modeStatusText(modeN)

    text = sprintf(', n = %d', modeN);
end

function nSlice = phaseDimensionSize(dim, speciesData)

    switch dim
        case 1
            nSlice = speciesData.gridE;
        case 2
            nSlice = speciesData.gridPphi;
        case 3
            nSlice = speciesData.gridLambda;
        otherwise
            error('切片维度必须为 1、2 或 3。');
    end
end

function label = phaseDimensionLabel(dim)

    labels = {'E', 'Pphi', 'Lambda'};
    label = labels{dim};
end

function idx = defaultSliceIndex(dim, gridE, gridPphi, gridLambda)

    switch dim
        case 1
            idx = max(1, round(gridE / 2));
        case 2
            idx = max(1, round(gridPphi / 2));
        case 3
            idx = max(1, round(gridLambda / 2));
        otherwise
            error('切片维度必须为 1、2 或 3。');
    end
end

function idx = clampIndex(idx, maxIndex, label)

    assert(isscalar(idx) && isfinite(idx) && idx == floor(idx), ...
        '%s 下标必须是有限整数。', label);
    assert(idx >= 1 && idx <= maxIndex, ...
        '%s 下标 %d 超出有效范围 [1, %d]。', label, idx, maxIndex);
end

function branchName = normalizeFrequencyBranch(branchText)

    branchName = lower(strtrim(char(branchText)));
    if ~ismember(branchName, {'para', 'anti'})
        error('frequency branch 必须为 "para" 或 "anti"。');
    end
end

function branchName = normalizeResonanceBranch(branchText)

    branchName = lower(strtrim(char(branchText)));
    if ~ismember(branchName, {'para', 'anti', 'trapped'})
        error('resonance branch 必须为 "para"、"anti" 或 "trapped"。');
    end
end

function branchNames = normalizeResonanceBranches(branchText)

    if iscell(branchText)
        rawNames = branchText(:).';
    elseif isstring(branchText)
        rawNames = cellstr(branchText(:).');
    else
        branchText = char(branchText);
        rawNames = regexp(strtrim(branchText), '[,;|\s]+', 'split');
    end

    branchNames = {};
    for iName = 1:numel(rawNames)
        nameText = lower(strtrim(char(rawNames{iName})));
        if isempty(nameText)
            continue;
        elseif strcmp(nameText, 'all')
            branchNames = [branchNames, {'para', 'anti', 'trapped'}]; %#ok<AGROW>
        else
            branchNames{end + 1} = normalizeResonanceBranch(nameText); %#ok<AGROW>
        end
    end

    assert(~isempty(branchNames), ...
        'resonance branch must include at least one of "para", "anti", or "trapped".');
    branchNames = unique(branchNames, 'stable');
end

function [unitScale, colorbarLabel] = frequencyUnitScale(unitText)

    key = lower(strtrim(char(unitText)));
    switch key
        case 'hz'
            unitScale = 1;
            colorbarLabel = '$f/\mathrm{Hz}$';
        case {'w', 'omega'}
            unitScale = 2 * pi;
            colorbarLabel = '$\omega/(\mathrm{rad/s})$';
        otherwise
            error('frequency unit 必须为 "Hz" 或 "w"。');
    end
end

function filePath = findExistingFile(searchDirs, fileName)

    filePath = '';
    for iDir = 1:numel(searchDirs)
        candidate = fullfile(searchDirs{iDir}, fileName);
        if isfile(candidate)
            filePath = candidate;
            return;
        end
    end
end

function data = readOrbitBinary(orbitFile)

    fid = fopen(orbitFile, 'rb');
    assert(fid >= 0, '无法打开文件：%s', orbitFile);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fseek(fid, 0, 'eof');
    fileSize = ftell(fid);
    fseek(fid, 0, 'bof');

    bytesPerRecord = 9 * 8;
    assert(mod(fileSize, bytesPerRecord) == 0, ...
        '%s 尺寸不匹配：文件大小 %d 字节不能被单条记录 %d 字节整除。', orbitFile, fileSize, bytesPerRecord);
    numRecords = fileSize / bytesPerRecord;
    rawData = fread(fid, numRecords * 9, 'double=>double');
    assert(numel(rawData) == numRecords * 9, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（numRecords=%d, columns=9）。', ...
        orbitFile, numel(rawData), numRecords * 9, numRecords);
    data = reshape(rawData, 9, numRecords).';
end

function records = orbitRecordColumns(data)

    records = struct( ...
        'Ids', int32(data(:, 1)), ...
        'orbits', data(:, 2), ...
        'dtheta', data(:, 3), ...
        'dphiTotal', data(:, 4), ...
        'dphiVpara', data(:, 5), ...
        'dTs', data(:, 6), ...
        'Es', data(:, 7), ...
        'Pphis', data(:, 8), ...
        'Lambdas', data(:, 9));
end

function data = readPhase3D(filePath, precision, gridE, gridPphi, gridLambda)

    raw = readBinaryVector(filePath, precision);
    expectedCount = gridE * gridPphi * gridLambda;
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（gridE=%d, gridPphi=%d, gridLambda=%d）。', ...
        filePath, numel(raw), expectedCount, gridE, gridPphi, gridLambda);
    data = reshape(raw, [gridLambda, gridPphi, gridE]);
    data = permute(data, [3, 2, 1]);
end

function data = readPhase4D(filePath, precision, gridE, gridPphi, gridLambda, expectedTime)

    raw = readBinaryVector(filePath, precision);
    perFrame = gridE * gridPphi * gridLambda;
    assert(mod(numel(raw), perFrame) == 0, ...
        '%s 尺寸不匹配：读到 %d 个数，不能被单帧长度 %d 整除（gridE=%d, gridPphi=%d, gridLambda=%d）。', ...
        filePath, numel(raw), perFrame, gridE, gridPphi, gridLambda);
    nTime = numel(raw) / perFrame;
    assert(nTime == expectedTime, ...
        '%s 尺寸不匹配：读到 %d 帧，期望 %d 帧（每帧 %d 个数，总数 %d）。', ...
        filePath, nTime, expectedTime, perFrame, numel(raw));
    data = reshape(raw, [gridLambda, gridPphi, gridE, nTime]);
    data = permute(data, [3, 2, 1, 4]);
end

function data = readPhasePower5D(filePath, precision, gridE, gridPphi, gridLambda, modeCount, expectedTime)

    raw = readBinaryVector(filePath, precision);
    perTime = gridE * gridPphi * gridLambda * modeCount;
    assert(mod(numel(raw), perTime) == 0, ...
        '%s 尺寸不匹配：读到 %d 个数，不能被单个输出时刻长度 %d 整除（gridE=%d, gridPphi=%d, gridLambda=%d, modeCount=%d）。', ...
        filePath, numel(raw), perTime, gridE, gridPphi, gridLambda, modeCount);
    nTime = numel(raw) / perTime;
    assert(nTime == expectedTime, ...
        '%s 尺寸不匹配：读到 %d 个输出时刻，期望 %d 个（每时刻 %d 个数，总数 %d）。', ...
        filePath, nTime, expectedTime, perTime, numel(raw));
    data = reshape(raw, [gridLambda, gridPphi, gridE, modeCount, nTime]);
    data = permute(data, [3, 2, 1, 4, 5]);
end

function data = readPitch2D(filePath, precision, gridVpara, gridVperp)

    raw = readBinaryVector(filePath, precision);
    expectedCount = gridVpara * gridVperp;
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（gridVpara=%d, gridVperp=%d）。', ...
        filePath, numel(raw), expectedCount, gridVpara, gridVperp);
    data = reshape(raw, [gridVperp, gridVpara]);
    data = permute(data, [2, 1]);
end

function data = readPitch3D(filePath, precision, gridVpara, gridVperp, expectedTime)

    raw = readBinaryVector(filePath, precision);
    perFrame = gridVpara * gridVperp;
    assert(mod(numel(raw), perFrame) == 0, ...
        '%s 尺寸不匹配：读到 %d 个数，不能被单帧长度 %d 整除（gridVpara=%d, gridVperp=%d）。', ...
        filePath, numel(raw), perFrame, gridVpara, gridVperp);
    nTime = numel(raw) / perFrame;
    assert(nTime == expectedTime, ...
        '%s 尺寸不匹配：读到 %d 帧，期望 %d 帧（每帧 %d 个数，总数 %d）。', ...
        filePath, nTime, expectedTime, perFrame, numel(raw));
    data = reshape(raw, [gridVperp, gridVpara, nTime]);
    data = permute(data, [2, 1, 3]);
end

function data = readPitchPower4D(filePath, precision, gridVpara, gridVperp, modeCount, expectedTime)

    raw = readBinaryVector(filePath, precision);
    perTime = gridVpara * gridVperp * modeCount;
    assert(mod(numel(raw), perTime) == 0, ...
        '%s 尺寸不匹配：读到 %d 个数，不能被单个输出时刻长度 %d 整除（gridVpara=%d, gridVperp=%d, modeCount=%d）。', ...
        filePath, numel(raw), perTime, gridVpara, gridVperp, modeCount);
    nTime = numel(raw) / perTime;
    assert(nTime == expectedTime, ...
        '%s 尺寸不匹配：读到 %d 个输出时刻，期望 %d 个（每时刻 %d 个数，总数 %d）。', ...
        filePath, nTime, expectedTime, perTime, numel(raw));
    data = reshape(raw, [gridVperp, gridVpara, modeCount, nTime]);
    data = permute(data, [2, 1, 3, 4]);
end

function data = readDiffusivity3D(filePath, precision, expectedTime, modeCount, gridNx)

    raw = readBinaryVector(filePath, precision);
    expectedCount = expectedTime * modeCount * gridNx;
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（expectedTime=%d, modeCount=%d, gridNx=%d）。', ...
        filePath, numel(raw), expectedCount, expectedTime, modeCount, gridNx);
    data = reshape(raw, [gridNx, modeCount, expectedTime]);
    data = permute(data, [3, 2, 1]);
end

function raw = readBinaryVector(filePath, precision)

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

function [gridE, gridPphi, gridLambda] = readPhaseGrid(normData)

    requiredFields = {'gridE', 'gridPphi', 'gridLambda'};
    for fieldIndex = 1:numel(requiredFields)
        assert(isfield(normData, requiredFields{fieldIndex}), ...
            'normalization2D.mat 缺少 "%s"。', requiredFields{fieldIndex});
    end

    gridE = readPositiveIntegerScalar(normData.gridE, 'gridE');
    gridPphi = readPositiveIntegerScalar(normData.gridPphi, 'gridPphi');
    gridLambda = readPositiveIntegerScalar(normData.gridLambda, 'gridLambda');
end

function phaseRange = readSpeciesRange(normData, rangeField)

    assert(isfield(normData, rangeField), 'normalization2D.mat 缺少 "%s"。', rangeField);
    phaseRange = reshape(double(normData.(rangeField)), 1, []);
    assert(numel(phaseRange) == 6 && all(isfinite(phaseRange)), ...
        '%s 必须包含 [minE maxE minPphi maxPphi minLambda maxLambda]。', rangeField);
end

function value = readIntParam(paramText, name)

    token = regexp(paramText, ['const\s+int\s+' name '\s*=\s*(\d+)\s*;'], 'tokens', 'once');
    assert(~isempty(token), '找不到整数参数：%s', name);
    value = str2double(token{1});
end

function value = readSwitchParam(paramText, name)

    token = regexp(paramText, ['using\s+' name '\s*=\s*(trueType|falseType)\s*;'], 'tokens', 'once');
    assert(~isempty(token), '找不到开关参数：%s', name);
    value = strcmp(token{1}, 'trueType');
end

function value = readFloatParam(paramText, name)

    token = regexp(paramText, ['const\s+(?:double|float|mhdReal|picReal)\s+' name ...
        '\s*=\s*([-+]?[0-9eE+\-\.]+)\s*;'], 'tokens', 'once');
    assert(~isempty(token), '找不到浮点参数：%s', name);
    value = str2double(token{1});
end

function precision = readMHDPrecisionParam(paramText)

    token = regexp(paramText, 'using\s+mhdReal\s*=\s*(double|float)\s*;', 'tokens', 'once');
    assert(~isempty(token), '找不到 mhdReal 精度定义。');
    precision = token{1};
end

function value = readPositiveScalar(normData, fieldName)

    assert(isfield(normData, fieldName), 'normalization2D.mat 缺少 "%s"。', fieldName);
    value = double(normData.(fieldName));
    assert(isscalar(value) && isfinite(value) && value > 0, '%s 必须是正标量。', fieldName);
end

function value = readPositiveIntegerScalar(rawValue, fieldName)

    value = double(rawValue);
    assert(isscalar(value) && isfinite(value) && value == floor(value) && value > 0, ...
        '%s 必须是正整数标量。', fieldName);
end

function value = readNonnegativeIntegerScalar(rawValue, fieldName)

    value = double(rawValue);
    assert(isscalar(value) && isfinite(value) && value == floor(value) && value >= 0, ...
        '%s 必须是非负整数标量。', fieldName);
end

function value = validateFiniteScalar(rawValue, fieldName)

    value = double(rawValue);
    assert(isscalar(value) && isfinite(value), '%s 必须是有限标量。', fieldName);
end

function value = validateIntegerScalar(rawValue, fieldName)

    value = double(rawValue);
    assert(isscalar(value) && isfinite(value) && value == floor(value), '%s 必须是整数标量。', fieldName);
end

function range = normalizeSliderRange(rawRange, initialValue, fieldName, isInteger)

    range = reshape(double(rawRange), 1, []);
    assert(numel(range) == 2 && all(isfinite(range)), '%s 必须为 [min max]。', fieldName);
    range = sort(range);
    initialValue = double(initialValue);
    range(1) = min(range(1), initialValue);
    range(2) = max(range(2), initialValue);

    if isInteger
        range = [floor(range(1)), ceil(range(2))];
        if range(1) == range(2)
            range = [range(1) - 1, range(2) + 1];
        end
    elseif range(1) == range(2)
        deltaValue = max(1, 0.1 * abs(initialValue));
        range = [range(1) - deltaValue, range(2) + deltaValue];
    end
end

function range = normalizeSliderRangeForValues(rawRange, initialValues, fieldName, isInteger)

    range = reshape(double(rawRange), 1, []);
    initialValues = reshape(double(initialValues), 1, []);
    assert(numel(range) == 2 && all(isfinite(range)), '%s 必须为 [min max]。', fieldName);
    assert(~isempty(initialValues) && all(isfinite(initialValues)), '%s 初始值必须有限。', fieldName);
    range = sort(range);
    range(1) = min([range(1), initialValues]);
    range(2) = max([range(2), initialValues]);

    if isInteger
        range = [floor(range(1)), ceil(range(2))];
        if range(1) == range(2)
            range = [range(1) - 1, range(2) + 1];
        end
    elseif range(1) == range(2)
        deltaValue = max(1, 0.1 * max(abs(initialValues)));
        range = [range(1) - deltaValue, range(2) + deltaValue];
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

function value = clampInteger(value, minValue, maxValue)

    value = round(double(value));
    value = min(max(value, minValue), maxValue);
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

function step = sliderStep(maxValue)

    if maxValue <= 1
        step = [1, 1];
    else
        smallStep = 1 / (maxValue - 1);
        largeStep = min(1, max(1, round(maxValue / 20)) / (maxValue - 1));
        step = [smallStep, largeStep];
    end
end

function data = localIdsToLinear(localIds, gridE, gridPphi, gridLambda)

    [iE, iPphi, iLambda] = localIdsToSubscripts(localIds, gridPphi, gridLambda);
    data = sub2ind([gridE, gridPphi, gridLambda], iE, iPphi, iLambda);
end

function [iE, iPphi, iLambda] = localIdsToSubscripts(localIds, gridPphi, gridLambda)

    localIds = double(localIds(:));
    strideE = gridPphi * gridLambda;
    iE = floor(localIds / strideE) + 1;
    remainder = mod(localIds, strideE);
    iPphi = floor(remainder / gridLambda) + 1;
    iLambda = mod(remainder, gridLambda) + 1;
end

function [baseE, basePphi, baseLambda] = initialCoordinatesFromLocalId(localIds, E1d, Pphi1d, Lambda1d)

    [iE, iPphi, iLambda] = localIdsToSubscripts(localIds, numel(Pphi1d), numel(Lambda1d));
    baseE = E1d(iE(:));
    basePphi = Pphi1d(iPphi(:));
    baseLambda = Lambda1d(iLambda(:));
end

function plotConservationDiagnostics(speciesName, diagnostics, E1d, Pphi1d, Lambda1d)

    classNames = {'trapped', 'para', 'anti'};
    for classIndex = 1:numel(classNames)
        className = classNames{classIndex};
        if ~isfield(diagnostics, className) || isempty(diagnostics.(className).localId)
            fprintf('[plot] %s %s 无守恒量诊断点。\n', speciesName, className);
            continue;
        end
        plotInvariantErrors(speciesName, className, diagnostics.(className), E1d, Pphi1d, Lambda1d);
    end
end

function plotInvariantErrors(speciesName, classLabel, diagnostic, E1d, Pphi1d, Lambda1d)

    [baseE, basePphi, baseLambda] = initialCoordinatesFromLocalId(diagnostic.localId, E1d, Pphi1d, Lambda1d);
    invariantLabels = {'E', 'Pphi', 'Lambda'};
    invariantErrors = {mixedConservationError(diagnostic.E, baseE), ...
        mixedConservationError(diagnostic.Pphi, basePphi), ...
        mixedConservationError(diagnostic.Lambda, baseLambda)};

    figure('Name', [speciesName ' ' classLabel ' invariant error'], 'Color', 'w', 'Position', [100, 100, 900, 760]);
    for iPlot = 1:numel(invariantLabels)
        axHandle = subplot(3, 1, iPlot);
        plot(axHandle, invariantErrors{iPlot}, '.');
        grid(axHandle, 'on');
        ylabel(axHandle, invariantLabels{iPlot}, 'FontName', 'Times New Roman', 'FontSize', 14);
        set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14);
        if iPlot == 1
            title(axHandle, sprintf('%s %s relative error', speciesName, classLabel), ...
                'Interpreter', 'none', 'FontName', 'Times New Roman', 'FontSize', 14);
        elseif iPlot == numel(invariantLabels)
            xlabel(axHandle, 'particle index', 'FontName', 'Times New Roman', 'FontSize', 14);
        end
    end
end

function errorValue = mixedConservationError(finalValue, baselineValue)

    finalValue = finalValue(:);
    baselineValue = baselineValue(:);
    errorValue = finalValue - baselineValue;
    relativeMask = abs(baselineValue) > 1e-12;
    errorValue(relativeMask) = errorValue(relativeMask) ./ baselineValue(relativeMask);
end

function Z = fillEnclosedBlankRegions(Z)

    [ny, nx] = size(Z);
    if ny < 3 || nx < 3
        return;
    end

    sourceZ = Z;
    blankMask = ~isfinite(sourceZ);
    if ~any(blankMask(:))
        return;
    end

    visited = false(ny, nx);
    offsets = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 1; 1, -1; 1, 0; 1, 1];
    blankIndex = find(blankMask);

    for startIndex = blankIndex(:)'
        if visited(startIndex)
            continue;
        end

        queue = zeros(numel(blankIndex), 1);
        component = zeros(numel(blankIndex), 1);
        boundaryIndex = zeros(numel(blankIndex) * 8, 1);
        head = 1;
        tail = 1;
        nComponent = 0;
        nBoundary = 0;
        touchesBoundary = false;
        queue(tail) = startIndex;
        visited(startIndex) = true;

        while head <= tail
            currentIndex = queue(head);
            head = head + 1;
            nComponent = nComponent + 1;
            component(nComponent) = currentIndex;
            [i, j] = ind2sub([ny, nx], currentIndex);
            touchesBoundary = touchesBoundary || i == 1 || i == ny || j == 1 || j == nx;

            for k = 1:8
                ni = i + offsets(k, 1);
                nj = j + offsets(k, 2);
                if ni < 1 || ni > ny || nj < 1 || nj > nx
                    touchesBoundary = true;
                    continue;
                end
                neighborIndex = sub2ind([ny, nx], ni, nj);
                if blankMask(neighborIndex)
                    if ~visited(neighborIndex)
                        tail = tail + 1;
                        queue(tail) = neighborIndex;
                        visited(neighborIndex) = true;
                    end
                else
                    nBoundary = nBoundary + 1;
                    boundaryIndex(nBoundary) = neighborIndex;
                end
            end
        end

        if ~touchesBoundary && nBoundary > 0
            boundaryIndex = unique(boundaryIndex(1:nBoundary));
            Z(component(1:nComponent)) = mean(sourceZ(boundaryIndex));
        end
    end
end

function tf = residualHasZeroContour(Z)

    if ~any(isfinite(Z(:)))
        tf = false;
        return;
    end
    if any(Z(:) == 0)
        tf = true;
        return;
    end
    if size(Z, 1) < 2 || size(Z, 2) < 2
        tf = false;
        return;
    end

    z00 = Z(1:end-1, 1:end-1);
    z10 = Z(2:end, 1:end-1);
    z01 = Z(1:end-1, 2:end);
    z11 = Z(2:end, 2:end);
    finiteCell = isfinite(z00) & isfinite(z10) & isfinite(z01) & isfinite(z11);
    cellMin = min(min(z00, z10), min(z01, z11));
    cellMax = max(max(z00, z10), max(z01, z11));
    tf = any(finiteCell(:) & cellMin(:) < 0 & cellMax(:) > 0);
end

function [colorMin, colorMax] = finiteColorLimits(data)

    finiteData = data(isfinite(data));
    if isempty(finiteData)
        colorMin = -1;
        colorMax = 1;
        return;
    end

    colorMin = min(finiteData);
    colorMax = max(finiteData);
    robustMax = finitePercentile(finiteData, 99.5);
    if isfinite(robustMax) && robustMax > colorMin && colorMax > 5 * robustMax
        colorMax = robustMax;
    end
    if colorMin == colorMax
        colorMin = colorMin - 1;
        colorMax = colorMax + 1;
    end
end

function [colorMin, colorMax] = symmetricFiniteColorLimits(data)

    finiteData = data(isfinite(data));
    if isempty(finiteData)
        colorMin = -1;
        colorMax = 1;
        return;
    end

    absData = abs(finiteData);
    colorLimit = max(absData);
    robustLimit = finitePercentile(absData, 99.5);
    if isfinite(robustLimit) && robustLimit > 0 && colorLimit > 5 * robustLimit
        colorLimit = robustLimit;
    end
    if colorLimit <= 0
        colorLimit = 1;
    end
    colorMin = -colorLimit;
    colorMax = colorLimit;
end

function value = finitePercentile(data, percent)

    data = sort(data(isfinite(data)));
    if isempty(data)
        value = NaN;
        return;
    end
    if isscalar(data)
        value = data(1);
        return;
    end

    position = 1 + (numel(data) - 1) * percent / 100;
    lowerIndex = floor(position);
    upperIndex = ceil(position);
    weight = position - lowerIndex;
    value = (1 - weight) * data(lowerIndex) + weight * data(upperIndex);
end

function cmap = phaseSpaceColormap(isSigned, colormapIndex, nColor)

    if isSigned
        cmap = redblue(nColor);
        return;
    end

    switch colormapIndex
        case 1
            if exist('turbo', 'file') == 2 || exist('turbo', 'builtin') == 5
                cmap = turbo(nColor);
            else
                cmap = parula(nColor);
            end
        case 2
            cmap = jet(nColor);
        otherwise
            error('colormapIndex 必须为 1 或 2。');
    end
end

function cmap = redblue(nColor)

    if nargin < 1
        nColor = 256;
    end
    anchorX = [-1.0; -0.5; 0.0; 0.5; 1.0];
    anchorC = [ ...
        0.05, 0.10, 0.40; ...
        0.00, 0.00, 1.00; ...
        0.97, 0.97, 0.97; ...
        1.00, 0.00, 0.00; ...
        0.40, 0.00, 0.05];
    x = linspace(-1, 1, nColor).';
    cmap = interp1(anchorX, anchorC, x);
end

function summary = initializeOrbitSummary(speciesName, orbitFile, numRecords, nPhase, validRecord, orbits)

    expectedRecords = 2 * nPhase;
    initializedRecords = sum(validRecord);
    orbitCounts = struct( ...
        'pad', sum(orbits < 1), ...
        'loss', sum(orbits > 1 & orbits < 2), ...
        'unknown', sum(orbits > 5), ...
        'para', sum(orbits > 2 & orbits < 3), ...
        'anti', sum(orbits > 3 & orbits < 4), ...
        'trapped', sum(orbits > 4 & orbits < 5), ...
        'exactPara', sum(isOrbit(orbits, 2.5)), ...
        'exactAnti', sum(isOrbit(orbits, 3.5)), ...
        'exactTrapped', sum(isOrbit(orbits, 4.5)));

    summary = struct( ...
        'species', speciesName, ...
        'orbitFile', orbitFile, ...
        'records', numRecords, ...
        'expectedRecords', expectedRecords, ...
        'nPhase', nPhase, ...
        'initializedRecords', initializedRecords, ...
        'initializedRatio', safeDivide(initializedRecords, expectedRecords), ...
        'orbitCounts', orbitCounts);
end

function printSpeciesSummary(summary)

    totalRecords = summary.expectedRecords;
    initialized = summary.initializedRecords;
    fprintf('pad      : %d (%.2f%% of total)\n', summary.orbitCounts.pad, 100 * safeDivide(summary.orbitCounts.pad, totalRecords));
    fprintf('effective: %d (%.2f%% of total)\n\n', initialized, 100 * summary.initializedRatio);
    fprintf('loss     : %d (%.2f%% of total, %.2f%% of effective)\n', ...
        summary.orbitCounts.loss, 100 * safeDivide(summary.orbitCounts.loss, totalRecords), ...
        100 * safeDivide(summary.orbitCounts.loss, initialized));
    fprintf('unknown  : %d (%.2f%% of total, %.2f%% of effective)\n', ...
        summary.orbitCounts.unknown, 100 * safeDivide(summary.orbitCounts.unknown, totalRecords), ...
        100 * safeDivide(summary.orbitCounts.unknown, initialized));
    fprintf('para     : %d (%.2f%% of total, %.2f%% of effective)\n', ...
        summary.para.count, 100 * safeDivide(summary.para.count, totalRecords), ...
        100 * safeDivide(summary.para.count, initialized));
    fprintf('anti     : %d (%.2f%% of total, %.2f%% of effective)\n', ...
        summary.anti.count, 100 * safeDivide(summary.anti.count, totalRecords), ...
        100 * safeDivide(summary.anti.count, initialized));
    fprintf('trapped  : %d (%.2f%% of total, %.2f%% of effective)\n', ...
        summary.trapped.branchCount, 100 * safeDivide(summary.trapped.branchCount, totalRecords), ...
        100 * safeDivide(summary.trapped.branchCount, initialized));
end

function assertNoNaN(speciesName, records)

    nanCount = sum(isnan(double(records.Ids))) + sum(isnan(records.orbits)) + ...
        sum(isnan(records.dtheta)) + sum(isnan(records.dphiTotal)) + ...
        sum(isnan(records.dphiVpara)) + sum(isnan(records.dTs)) + ...
        sum(isnan(records.Es)) + sum(isnan(records.Pphis)) + sum(isnan(records.Lambdas));
    assert(nanCount == 0, '%s PhaseSpaceOrbit.bin 中存在 NaN，数量为 %d。', speciesName, nanCount);
end

function warnIfSignedIdOrderLooksWrong(speciesName, Ids, recordBranch, validRecord)

    positiveBranchBad = any(Ids(validRecord & recordBranch > 0) < 0);
    negativeBranchBad = any(Ids(validRecord & recordBranch < 0 & Ids ~= 0) > 0);
    if positiveBranchBad || negativeBranchBad
        warning('%s 的 ID 符号与 [正向分支; 反向分支] 文件顺序不完全一致，将继续按文件顺序判断分支。', speciesName);
    end
end

function tf = isPadRecord(Ids)

    tf = Ids == int32(20251106);
end

function tf = isOrbit(orbitValues, targetOrbit)

    tf = abs(orbitValues - targetOrbit) <= orbitTolerance();
end

function value = orbitTolerance()

    value = 1e-12;
end

function value = trappedOrbitRelativeTolerance()

    value = 0.05;
end

function value = averagePair(a, b)

    value = 0.5 * (a(:) + b(:));
end

function relDiff = relativePairDifference(a, b)

    a = a(:);
    b = b(:);
    denominator = max(abs(a), abs(b));
    denominator = max(denominator, 1e-12);
    relDiff = abs(a - b) ./ denominator;
end

function value = safeDivide(numerator, denominator)

    if denominator == 0
        value = NaN;
    else
        value = numerator / denominator;
    end
end

function value = maxOrNaN(x)

    if isempty(x)
        value = NaN;
    else
        value = max(x);
    end
end

function value = meanOrNaN(x)

    if isempty(x)
        value = NaN;
    else
        value = mean(x);
    end
end

function value = optionOrDefault(opt, fieldName, defaultValue)

    if isfield(opt, fieldName) && ~isempty(opt.(fieldName))
        value = opt.(fieldName);
    else
        value = defaultValue;
    end
end

function logLoaded(name, data)

    fprintf('[load] %s: size=[%s]\n', char(name), formatSize(size(data)));
end

function logSkipped(name, reason)

    fprintf('[skip] %s: %s\n', char(name), reason);
end

function printOrbitSummaryIndex(summaryStruct)

    fields = fieldnames(summaryStruct);
    if isempty(fields)
        fprintf('[orbit] 未处理任何 PhaseSpaceOrbit.bin。\n');
    else
        fprintf('[orbit] 已处理物种：%s\n', strjoin(fields, ', '));
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
