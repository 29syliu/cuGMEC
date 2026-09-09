%% cuGMEC PIC 诊断统一可视化脚本

%{
功能：
统一可视化 PIC phase / pitch / orbit 诊断。

输入目录：
将相关输入、诊断和输出文件放在同一个 inputDir 中。

必须包含：
cuGMEC_param.h, normalization2D.mat, NTP.mat,
plot2D.mat 或 plot3D.mat

按开关读取：
ifIon, ifAlpha, ifBeam,
ifDiagDensity, ifDiagDiffusivity,
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
Ion/Alpha/BeamDiffusivity.bin,
Ion/Alpha/BeamDensity.bin,
Ion/Alpha/BeamPhaseSpaceMapping.bin
%}

%% 用户设置

clear; close all;

inputDir = 'C:\Users\Desktop\test';

paramFile = fullfile(inputDir, 'cuGMEC_param.h');
normalizationFile = fullfile(inputDir, 'normalization2D.mat');
ntpFile = fullfile(inputDir, 'NTP.mat');
resonanceDetuningPlotFile = fullfile(inputDir, 'plot2D.mat');

speciesList = {'Ion', 'Alpha', 'Beam'};

%% 读取所有输入

assert(isfile(paramFile), '缺少参数文件：%s', paramFile);
assert(isfile(normalizationFile), '缺少 MAT 文件：%s', normalizationFile);
assert(isfolder(inputDir), '缺少输入目录：%s', inputDir);

paramText = fileread(paramFile);
normData = load(normalizationFile);
ntpData = struct();
if isfile(ntpFile)
    ntpData = load(ntpFile);
end
meta = readPICMetadata(paramText, normData);
speciesList = enabledPICSpecies(speciesList, meta);
searchDirs = {inputDir};

picPhaseData = initializePICSpeciesData(speciesList, ...
    @(speciesName) initializePhaseSpecies(speciesName, normData, meta));
picPitchData = initializePICSpeciesData(speciesList, ...
    @(speciesName) initializePitchSpecies(speciesName, paramText, meta));
picDiffusivityData = initializePICSpeciesData(speciesList, ...
    @(speciesName) initializeDiffusivitySpecies(speciesName, meta));
picDensityData = initializePICSpeciesData(speciesList, ...
    @(speciesName) initializeDensitySpecies(speciesName, meta, ntpData));

orbitRaw = readAllOrbitRaw(speciesList, searchDirs, meta);
picPhaseData = readAllPhaseDiagnostics(picPhaseData, speciesList, searchDirs, meta);
picPitchData = readAllPitchDiagnostics(picPitchData, speciesList, searchDirs, meta);
picDiffusivityData = readAllDiffusivityDiagnostics(picDiffusivityData, speciesList, searchDirs, meta);
picDensityData = readAllDensityDiagnostics(picDensityData, speciesList, searchDirs, meta);
picResonanceDetuningInput = readAllResonanceDetuningInputs( ...
    resonanceDetuningPlotFile, inputDir, speciesList, meta);

picWorkspace = struct();
picWorkspace.meta = meta;
picWorkspace.phase = picPhaseData;
picWorkspace.pitch = picPitchData;
picWorkspace.diffusivity = picDiffusivityData;
picWorkspace.density = picDensityData;

%% 处理 orbit.bin

picOrbitProcessOpt = struct( ...
    'enabled', true, ...
    'plotConservationDiagnostics', false);
%{
enabled                    : 是否处理 orbit 文件。
plotConservationDiagnostics: 是否绘制守恒量误差。
%}

phaseSpaceOrbit = struct();
phaseSpaceOrbitSummary = struct();

if picOrbitProcessOpt.enabled
    [phaseSpaceOrbit, phaseSpaceOrbitSummary] = processAllOrbitRaw( ...
        orbitRaw, speciesList, picPhaseData, picOrbitProcessOpt);
    printOrbitSummarySpecies(phaseSpaceOrbitSummary);
else
    logSkipped('PhaseSpaceOrbit processing', '处理开关为 false');
end

picWorkspace.orbit = phaseSpaceOrbit;
picWorkspace.orbitSummary = phaseSpaceOrbitSummary;

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
    'detuningLine', picPhaseDetuningLineOptions(false, [], [], [], []), ...
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
detuningLine    : 固定 Lambda 图上的黑色失谐路径线段，由 picPhaseDetuningLineOptions 构造。
detuningLine.enabled        : 是否叠加该线段；仅 fixedCoordinate='Lambda' 时生效。
detuningLine.E0             : 线段中心能量坐标。
detuningLine.Pphi0          : 图上显示的 Pphi 中心坐标。
detuningLine.deltaPphiLeft  : 图上 Pphi 向左长度。
detuningLine.deltaPphiRight : 图上 Pphi 向右长度。
interactive      : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
% For phase quantity interactive=1/2, resonance.branch may also be a list,
% for example {'para','anti'}, {'para','trapped'}, {'anti','trapped'},
% {'para','anti','trapped'}, ["para","anti"], or 'all'. Frequency, n, and m
% sliders are shared; each selected branch gets its own l min / l max pair.
picWorkspace.phaseQuantity = runPICPhaseQuantityPlot(picPhaseData, phaseSpaceOrbit, meta, picPhaseQuantityOpt);

%% 可视化共振失谐

picResonanceDetuningOpt = struct( ...
    'enabled', false, ...
    'species', 'Alpha', ...
    'plotFile', resonanceDetuningPlotFile, ...
    'branch', 'trapped', ...
    'frequencyHz', 65e3, ...
    'toroidalMode', 30, ...
    'poloidalMode', 33, ...
    'harmonic', 0, ...
    'E0', [], ...
    'Pphi0', [], ...
    'Lambda0', [], ...
    'deltaPphiLeft', [], ...
    'deltaPphiRight', [], ...
    'sampleCount', 201);
%{
enabled            : 是否绘制共振失谐 R(x)。
species            : 'Ion', 'Alpha', 'Beam'。
plotFile           : plot2D.mat 或 plot3D.mat；当前按 NFP=1 的 Rplot(:,thetaMid) 建立 rho->R 映射。
branch             : 'para', 'anti' 或 'trapped'。
frequencyHz        : 有符号模频率，单位 Hz。
toroidalMode       : 真实物理环向模数 n。
poloidalMode       : 正极向模数 m；trapped 分支不使用。
harmonic           : 轨道谐波 l。
E0                 : 参考能量坐标，使用真实网格坐标。
Pphi0              : 图上显示的 Pphi 坐标。
Lambda0            : 参考 Lambda 坐标，使用真实网格坐标。
deltaPphiLeft      : 图上 Pphi 向左采样长度。
deltaPphiRight     : 图上 Pphi 向右采样长度。
sampleCount        : Pphi 路径采样点数，建议为奇数。
%}
picWorkspace.resonanceDetuning = runPICResonanceDetuningPlot( ...
    picPhaseData, phaseSpaceOrbit, meta, picResonanceDetuningInput, picResonanceDetuningOpt);

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
picWorkspace.diffusivityPlot = runPICDiffusivityPlot(picDiffusivityData, picDiffusivityOpt);

%% 可视化扰动密度

picDensityOpt = struct( ...
    'enabled', true, ...
    'species', 'Alpha', ...
    'plotType', 1, ...
    'timeIndex', max(1, round(meta.nDiagTime / 2)), ...
    'radialIndex', max(1, round(meta.gridNx / 2)), ...
    'timeAxis', 'ms', ...
    'radialAxis', 'rho', ...
    'radialDensityMode', 'delta', ...
    'colormapIndex', 1, ...
    'interactive', 2);
%{
enabled       : 是否绘图。
species       : 'Ion', 'Alpha', 'Beam'。
plotType      : 1 径向剖面；2 时间曲线；3 二维图。
timeIndex     : plotType = 1 使用的诊断时间索引（1 到 nDiagTime）。
radialIndex   : plotType = 2 使用的径向索引（1 到 gridNx）。
timeAxis      : 'ta', 'ms', 's' 或 'steps'。
radialAxis    : 'rho' 或 'x'。
radialDensityMode : 'delta' 画扰动密度；'total' 在径向图同画初始密度和总密度。
colormapIndex : 二维图非负色表，可选 1-2；有正负固定红蓝。
interactive   : 0 不交互；1 滑块释放后更新；2 拖动滑块时连续更新。
%}
picWorkspace.densityPlot = runPICDensityPlot(picDensityData, picDensityOpt);

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

    workspace = runPICPlotDiagnostic(opt, 'orbit frequency plot', ...
        @() hasSpeciesFields(phaseSpaceOrbit, opt.species, {}), '未读取对应物种的 orbit 数据', ...
        @(dynamicUpdate) plotOrbitFrequency(phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICResonanceLinePlot(phaseSpaceOrbit, meta, opt)

    workspace = runPICPlotDiagnostic(opt, 'resonance line plot', ...
        @() hasSpeciesFields(phaseSpaceOrbit, opt.species, {}), '未读取对应物种的 orbit 数据', ...
        @(dynamicUpdate) plotResonanceLine(phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICPhaseQuantityPlot(picPhaseData, phaseSpaceOrbit, meta, opt)

    workspace = runPICPlotDiagnostic(opt, 'phase quantity plot', ...
        @() hasRequestedPICQuantity(picPhaseData, opt.species, opt.quantity), '未读取所需 phase 数据', ...
        @(dynamicUpdate) plotPhaseQuantity(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICResonanceDetuningPlot( ...
    picPhaseData, phaseSpaceOrbit, meta, detuningInput, opt)

    workspace = struct('options', opt, 'residualHz', [], 'mapping', [], 'gridPoint', [], ...
        'center', [], 'path', [], 'fit', []);
    if ~opt.enabled
        logSkipped('resonance detuning plot', '绘图开关为 false');
        return;
    end
    if ~hasSpeciesFields(picPhaseData, opt.species, {})
        logSkipped('resonance detuning plot', '未读取对应物种的 phase 数据');
        return;
    end
    if ~hasSpeciesFields(phaseSpaceOrbit, opt.species, {})
        logSkipped('resonance detuning plot', '未读取对应物种的 orbit 数据');
        return;
    end

    opt = normalizeResonanceDetuningOptions(opt, meta);
    plotInput = buildResonanceDetuningPlotInput(detuningInput.plot, opt.plotFile, meta);
    [speciesData, speciesLabel] = resolveSpeciesData(picPhaseData, opt.species, 'phase');
    [orbitData, ~] = resolveSpeciesData(phaseSpaceOrbit, speciesLabel, 'orbit');

    resOpt = struct( ...
        'branch', opt.branch, ...
        'frequencyHz', opt.frequencyHz, ...
        'toroidalMode', opt.toroidalMode, ...
        'poloidalMode', opt.poloidalMode, ...
        'harmonic', opt.harmonic);
    [residualHz, ~] = calculateResonanceResidualField(orbitData, resOpt, meta);

    mapping = buildResonanceDetuningMapping( ...
        detuningInput.mapping, speciesLabel, speciesData, opt.branch, plotInput, meta);
    gridPoint = nearestResonanceDetuningGridPoint(speciesData, mapping, residualHz, opt);
    center = projectResonanceDetuningCenter(speciesData, mapping, residualHz, gridPoint);
    workspace.options = opt;
    workspace.species = speciesLabel;
    workspace.residualHz = residualHz;
    workspace.mapping = mapping;
    workspace.gridPoint = gridPoint;
    workspace.center = center;
    path = sampleResonanceDetuningPath(speciesData, mapping, residualHz, center, opt);
    workspace.path = path;
    if isfield(path, 'aborted') && path.aborted
        return;
    end
    fit = fitResonanceDetuningPath(path);

    printResonanceDetuningFitSummary(fit);
    plotResonanceDetuningPath(path, fit, center, opt);

    workspace.fit = fit;
end

function workspace = runPICPhasePowerPlot(picPhaseData, phaseSpaceOrbit, meta, opt)

    workspace = runPICPlotDiagnostic(opt, 'PhasePower plot', ...
        @() hasSpeciesFields(picPhaseData, opt.species, {'Power'}), '未读取 PhasePower 数据', ...
        @(dynamicUpdate) plotPhasePower(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate));
end

function workspace = runPICPitchQuantityPlot(picPitchData, opt)

    workspace = runPICPlotDiagnostic(opt, 'pitch quantity plot', ...
        @() hasRequestedPICQuantity(picPitchData, opt.species, opt.quantity), '未读取所需 pitch 数据', ...
        @(dynamicUpdate) plotPitchQuantity(picPitchData, opt, dynamicUpdate));
end

function workspace = runPICPitchPowerPlot(picPitchData, opt)

    workspace = runPICPlotDiagnostic(opt, 'PitchPower plot', ...
        @() hasSpeciesFields(picPitchData, opt.species, {'Power'}), '未读取 PitchPower 数据', ...
        @(dynamicUpdate) plotPitchPower(picPitchData, opt, dynamicUpdate));
end

function workspace = runPICDiffusivityPlot(picDiffusivityData, opt)

    workspace = runPICPlotDiagnostic(opt, 'Diffusivity plot', ...
        @() hasSpeciesFields(picDiffusivityData, opt.species, {'Diffusivity'}), '未读取 Diffusivity 数据', ...
        @(dynamicUpdate) plotDiffusivity(picDiffusivityData, opt, dynamicUpdate));
end

function workspace = runPICDensityPlot(picDensityData, opt)

    workspace = runPICPlotDiagnostic(opt, 'Density plot', ...
        @() hasSpeciesFields(picDensityData, opt.species, {'Density'}), '未读取 Density 数据', ...
        @(dynamicUpdate) plotDensity(picDensityData, opt, dynamicUpdate));
end

function workspace = runPICPlotDiagnostic(opt, plotName, availabilityFcn, missingReason, plotFcn)

    workspace = struct('options', opt);
    if ~opt.enabled
        logSkipped(plotName, '绘图开关为 false');
        return;
    end
    if ~availabilityFcn()
        logSkipped(plotName, missingReason);
        return;
    end
    runInteractivePlot(opt.interactive, @() plotFcn([]), plotFcn);
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

function detuningLine = picPhaseDetuningLineOptions(enabled, E0, Pphi0, deltaPphiLeft, deltaPphiRight)

    detuningLine.enabled = enabled;
    detuningLine.E0 = E0;
    detuningLine.Pphi0 = Pphi0;
    detuningLine.deltaPphiLeft = deltaPphiLeft;
    detuningLine.deltaPphiRight = deltaPphiRight;
end

function meta = readPICMetadata(paramText, normData)

    speciesNames = {'Ion', 'Alpha', 'Beam'};
    speciesSwitchNames = {'ifIon', 'ifAlpha', 'ifBeam'};
    meta.speciesEnabled = readNamedPICParameters( ...
        struct(), paramText, speciesNames, speciesSwitchNames, @readSwitchParam);

    switchNames = { ...
        'ifDiagDensity', 'ifDiagDiffusivity', 'ifOutputPhaseSpaceOrbit', ...
        'ifOutputPhaseSpaceJacobian', 'ifOutputPhaseSpaceF0', ...
        'ifOutputPhaseSpaceDeltaF', 'ifOutputPhaseSpacePower', ...
        'ifOutputPitchSpaceJacobian', 'ifOutputPitchSpaceF0', ...
        'ifOutputPitchSpaceDeltaF', 'ifOutputPitchSpacePower'};
    meta.switch = readNamedPICParameters( ...
        struct(), paramText, switchNames, switchNames, @readSwitchParam);

    integerFields = { ...
        'paramGridE', 'paramGridPphi', 'paramGridLambda', ...
        'gridVpara', 'gridVperp', 'gridNx', 'leftN', 'rightN', 'tubes', ...
        'totalSteps', 'outputSteps', 'diagSteps'};
    integerNames = { ...
        'gridE', 'gridPphi', 'gridLambda', ...
        'gridVpara', 'gridVperp', 'gridNx', 'leftN', 'rightN', 'tubes', ...
        'totalSteps', 'outputSteps', 'diagSteps'};
    meta = readNamedPICParameters( ...
        meta, paramText, integerFields, integerNames, @readIntParam);
    assert(meta.leftN <= meta.rightN, 'leftN 必须小于或等于 rightN。');
    assert(meta.tubes > 0, 'tubes 必须为正整数。');
    meta.modeIndexAll = meta.leftN:meta.rightN;
    meta.physicalNAll = meta.modeIndexAll * meta.tubes;

    meta.hasPhaseGrid = all(isfield(normData, {'gridE', 'gridPphi', 'gridLambda'}));
    if meta.hasPhaseGrid
        [meta.gridE, meta.gridPphi, meta.gridLambda] = readPhaseGrid(normData);
        assert(meta.paramGridE == meta.gridE && meta.paramGridPphi == meta.gridPphi && ...
            meta.paramGridLambda == meta.gridLambda, ...
            'cuGMEC_param.h 与 normalization2D.mat 中的 gridE/gridPphi/gridLambda 不一致。');
    else
        meta.gridE = meta.paramGridE;
        meta.gridPphi = meta.paramGridPphi;
        meta.gridLambda = meta.paramGridLambda;
    end

    meta.dt = readFloatParam(paramText, 'dt');
    assert(meta.gridNx > 1, 'gridNx 必须大于 1。');
    assert(meta.totalSteps >= 0 && meta.outputSteps > 0 && meta.diagSteps > 0, ...
        'totalSteps/outputSteps/diagSteps 取值不合法。');
    meta.nOutputTime = floor(meta.totalSteps / meta.outputSteps) + 1;
    meta.nDiagTime = floor(meta.totalSteps / meta.diagSteps) + 1;
    meta.mhdPrecision = readMHDPrecisionParam(paramText);
    meta.B0 = NaN;
    if isfield(normData, 'B0')
        meta.B0 = readPositiveScalar(normData, 'B0');
    end
    meta.L0 = readPositiveScalar(normData, 'L0');
    meta.VA0 = readPositiveScalar(normData, 'VA0');
    meta.MP = NaN;
    if isfield(normData, 'MP')
        meta.MP = readPositiveScalar(normData, 'MP');
    end
    meta.QE = NaN;
    if isfield(normData, 'QE')
        meta.QE = readPositiveScalar(normData, 'QE');
    end
    if isfield(normData, 'NFP')
        meta.NFP = requireFiniteScalarInRange(normData.NFP, -Inf, Inf, 'NFP');
    else
        meta.NFP = NaN;
    end
    assert(isfield(normData, 'RHO0') && isfield(normData, 'RHO1'), ...
        'normalization2D.mat 缺少 RHO0 或 RHO1。');
    meta.RHO0 = requireFiniteScalarInRange(normData.RHO0, -Inf, Inf, 'RHO0');
    meta.RHO1 = requireFiniteScalarInRange(normData.RHO1, -Inf, Inf, 'RHO1');
    assert(meta.RHO0 < meta.RHO1, 'RHO0 必须小于 RHO1。');
    meta.tDiag = (0:meta.nDiagTime - 1) * meta.diagSteps * meta.dt;
    meta.timeSeconds = meta.tDiag * meta.L0 / meta.VA0;
    meta.xGrid = linspace(0, 1, meta.gridNx);
    meta.rhoGrid = linspace(meta.RHO0, meta.RHO1, meta.gridNx);

end

function values = readNamedPICParameters(values, paramText, fieldNames, parameterNames, readerFcn)

    assert(numel(fieldNames) == numel(parameterNames), ...
        'PIC 参数字段名和源参数名数量必须一致。');
    for parameterIndex = 1:numel(fieldNames)
        values.(fieldNames{parameterIndex}) = readerFcn(paramText, parameterNames{parameterIndex});
    end
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

function speciesData = initializePhaseSpecies(speciesName, normData, meta)

    if ~meta.hasPhaseGrid
        speciesData = struct('species', speciesName, ...
            'J', [], 'F0', [], 'DF', [], 'Power', []);
        return;
    end
    phaseRange = readSpeciesRange(normData, [speciesName 'EPphiLambda']);
    speciesData = struct( ...
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
        'J', [], 'F0', [], 'DF', [], 'Power', []);
end

function speciesData = initializePitchSpecies(speciesName, paramText, meta)

    vmax = readFloatParam(paramText, [speciesName 'Vmax']);
    minVpara = -vmax;
    if strcmp(speciesName, 'Beam')
        minVpara = 0;
    end
    speciesData = struct( ...
        'species', speciesName, ...
        'gridVpara', meta.gridVpara, ...
        'gridVperp', meta.gridVperp, ...
        'modeIndexAll', meta.modeIndexAll, ...
        'physicalNAll', meta.physicalNAll, ...
        'tubes', meta.tubes, ...
        'Vpara1d', linspace(minVpara, vmax, meta.gridVpara), ...
        'Vperp1d', linspace(0, vmax, meta.gridVperp), ...
        'J', [], 'F0', [], 'DF', [], 'Power', []);
end

function speciesData = initializeDiffusivitySpecies(speciesName, meta)

    speciesData = struct( ...
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

function speciesData = initializeDensitySpecies(speciesName, meta, ntpData)

    [initialDensity, initialDensityError] = safeInitialDensityFromNTP(ntpData, speciesName, meta.rhoGrid);
    speciesData = struct( ...
        'species', speciesName, ...
        'gridNx', meta.gridNx, ...
        'diagSteps', meta.diagSteps, ...
        'xGrid', meta.xGrid, ...
        'rhoGrid', meta.rhoGrid, ...
        'tDiag', meta.tDiag, ...
        'timeSeconds', meta.timeSeconds, ...
        'InitialDensity', initialDensity, ...
        'InitialDensityError', initialDensityError, ...
        'Density', []);
end

function speciesData = initializePICSpeciesData(speciesList, initializerFcn)

    speciesData = struct();
    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        speciesData.(speciesName) = initializerFcn(speciesName);
    end
end

function [initialDensity, initialDensityError] = safeInitialDensityFromNTP(ntpData, speciesName, rhoGrid)

    initialDensity = [];
    initialDensityError = '';
    try
        initialDensity = initialDensityFromNTP(ntpData, speciesName, rhoGrid);
    catch err
        initialDensityError = err.message;
    end
end

function initialDensity = initialDensityFromNTP(ntpData, speciesName, rhoGrid)

    initialDensity = [];
    densityField = ntpDensityFieldName(speciesName);
    if ~isstruct(ntpData) || ~isfield(ntpData, 'rhoSample') || ~isfield(ntpData, densityField)
        return;
    end

    rhoSample = double(ntpData.rhoSample(:));
    densitySample = double(ntpData.(densityField)(:)) * 1e19;
    validSample = isfinite(rhoSample) & isfinite(densitySample);
    rhoSample = rhoSample(validSample);
    densitySample = densitySample(validSample);
    assert(numel(rhoSample) >= 2, ...
        'NTP.mat 中 rhoSample/%s 有效点数必须至少为 2。', densityField);

    [rhoSample, order] = sort(rhoSample);
    densitySample = densitySample(order);
    [rhoSample, uniqueIndex] = unique(rhoSample, 'stable');
    densitySample = densitySample(uniqueIndex);
    assert(numel(rhoSample) >= 2, ...
        'NTP.mat 中 rhoSample 去重后有效点数必须至少为 2。');

    rhoGrid = reshape(double(rhoGrid), 1, []);
    initialDensity = interp1(rhoSample, densitySample, rhoGrid, 'linear');
    if any(~isfinite(initialDensity))
        error(['NTP.mat 中 rhoSample 范围 [%.6g, %.6g] 不能覆盖模拟 rhoGrid 范围 ' ...
            '[%.6g, %.6g]，无法插值得到 %s 初始密度。'], ...
            min(rhoSample), max(rhoSample), min(rhoGrid), max(rhoGrid), speciesName);
    end
end

function densityField = ntpDensityFieldName(speciesName)

    switch lower(char(speciesName))
        case 'ion'
            densityField = 'niSample';
        case 'alpha'
            densityField = 'naSample';
        case 'beam'
            densityField = 'nbSample';
        otherwise
            error('不支持的 PIC 物种：%s。', speciesName);
    end
end

function orbitRaw = readAllOrbitRaw(speciesList, searchDirs, meta)

    orbitRaw = struct();
    if ~meta.switch.ifOutputPhaseSpaceOrbit
        logSkipped('PhaseSpaceOrbit', '开关 ifOutputPhaseSpaceOrbit 为 false');
        return;
    end
    if ~meta.hasPhaseGrid
        logSkipped('PhaseSpaceOrbit', 'normalization2D.mat 缺少相空间网格');
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
        'J', 'ifOutputPhaseSpaceJacobian', 'PhaseSpaceJacobian', 'double', ...
            @(filePath, precision) readPhase3D(filePath, precision, meta.gridE, meta.gridPphi, meta.gridLambda); ...
        'F0', 'ifOutputPhaseSpaceF0', 'PhaseSpaceF0', 'double', ...
            @(filePath, precision) readPhase3D(filePath, precision, meta.gridE, meta.gridPphi, meta.gridLambda); ...
        'DF', 'ifOutputPhaseSpaceDeltaF', 'PhaseDeltaF', meta.mhdPrecision, ...
            @(filePath, precision) readPhase4D(filePath, precision, meta.gridE, meta.gridPphi, meta.gridLambda, meta.nOutputTime); ...
        'Power', 'ifOutputPhaseSpacePower', 'PhasePower', meta.mhdPrecision, ...
            @(filePath, precision) readPhasePower5D(filePath, precision, meta.gridE, meta.gridPphi, ...
                meta.gridLambda, numel(meta.modeIndexAll), meta.nOutputTime)};

    if ~meta.hasPhaseGrid
        logSkipped('PhaseSpace diagnostics', 'normalization2D.mat 缺少相空间网格');
        return;
    end
    picPhaseData = readPICDiagnosticGroup( ...
        picPhaseData, speciesList, searchDirs, meta.switch, phaseSpecs);
end

function picPitchData = readAllPitchDiagnostics(picPitchData, speciesList, searchDirs, meta)

    pitchSpecs = { ...
        'J', 'ifOutputPitchSpaceJacobian', 'PitchSpaceJacobian', 'double', ...
            @(filePath, precision) readPitch2D(filePath, precision, meta.gridVpara, meta.gridVperp); ...
        'F0', 'ifOutputPitchSpaceF0', 'PitchSpaceF0', 'double', ...
            @(filePath, precision) readPitch2D(filePath, precision, meta.gridVpara, meta.gridVperp); ...
        'DF', 'ifOutputPitchSpaceDeltaF', 'PitchDeltaF', meta.mhdPrecision, ...
            @(filePath, precision) readPitch3D(filePath, precision, meta.gridVpara, meta.gridVperp, meta.nOutputTime); ...
        'Power', 'ifOutputPitchSpacePower', 'PitchPower', meta.mhdPrecision, ...
            @(filePath, precision) readPitchPower4D(filePath, precision, meta.gridVpara, ...
                meta.gridVperp, numel(meta.modeIndexAll), meta.nOutputTime)};
    picPitchData = readPICDiagnosticGroup( ...
        picPitchData, speciesList, searchDirs, meta.switch, pitchSpecs);
end

function picDiffusivityData = readAllDiffusivityDiagnostics(picDiffusivityData, speciesList, searchDirs, meta)

    diffusivitySpecs = { ...
        'Diffusivity', 'ifDiagDiffusivity', 'Diffusivity', meta.mhdPrecision, ...
            @(filePath, precision) readDiffusivity3D(filePath, precision, ...
                meta.nDiagTime, numel(meta.modeIndexAll), meta.gridNx)};
    picDiffusivityData = readPICDiagnosticGroup( ...
        picDiffusivityData, speciesList, searchDirs, meta.switch, diffusivitySpecs);
end

function picDensityData = readAllDensityDiagnostics(picDensityData, speciesList, searchDirs, meta)

    densitySpecs = { ...
        'Density', 'ifDiagDensity', 'Density', meta.mhdPrecision, ...
            @(filePath, precision) readDensity2D(filePath, precision, meta.nDiagTime, meta.gridNx)};
    picDensityData = readPICDiagnosticGroup( ...
        picDensityData, speciesList, searchDirs, meta.switch, densitySpecs);
end

function diagnosticData = readPICDiagnosticGroup( ...
    diagnosticData, speciesList, searchDirs, diagnosticSwitches, diagnosticSpecs)

    for diagnosticIndex = 1:size(diagnosticSpecs, 1)
        fieldName = diagnosticSpecs{diagnosticIndex, 1};
        switchName = diagnosticSpecs{diagnosticIndex, 2};
        fileSuffix = diagnosticSpecs{diagnosticIndex, 3};
        precision = diagnosticSpecs{diagnosticIndex, 4};
        readerFcn = diagnosticSpecs{diagnosticIndex, 5};
        if ~diagnosticSwitches.(switchName)
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

            loadedData = readerFcn(filePath, precision);
            diagnosticData.(speciesName).(fieldName) = loadedData;
            logLoaded(fileName, loadedData);
        end
    end
end

function [phaseSpaceOrbit, phaseSpaceOrbitSummary] = processAllOrbitRaw(orbitRaw, speciesList, picPhaseData, opt)

    phaseSpaceOrbit = struct();
    phaseSpaceOrbitSummary = struct();
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
        phaseSpaceOrbitSummary.(speciesName) = result.summary;
        printSpeciesSummary(result.summary);

        if getOptionValue(opt, 'plotConservationDiagnostics', false)
            plotConservationDiagnostics(speciesName, result.diagnostics, ...
                result.phaseSpaceData.E1d, result.phaseSpaceData.Pphi1d, result.phaseSpaceData.Lambda1d);
        end
    end
end

function result = analyzeSpeciesOrbitRecords(speciesName, rawOrbit, phaseGrid)

    records = orbitRecordColumns(rawOrbit.data);
    [records, nanRecordCount] = convertNaNOrbitRecordsToPad(records);
    if nanRecordCount > 0
        fprintf('[orbit] %s PhaseSpaceOrbit.bin: %d 个含 NaN 的轨道粒子/记录已按 pad 处理。\n', ...
            speciesName, nanRecordCount);
    end
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
        'maxAcceptedDTRelDiff', statisticOrNaN(dTRelDiff(accepted), @max, false), ...
        'meanAcceptedDTRelDiff', statisticOrNaN(dTRelDiff(accepted), @mean, false));
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

function plotOrbitFrequency(phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [orbitData, speciesLabel] = resolveSpeciesData(phaseSpaceOrbit, opt.species, 'orbit');
    branchName = normalizeFrequencyBranch(opt.branch);
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    initialSliceIndex = parsePhaseSlice(opt.slice, dim, orbitData);
    [unitScale, colorbarLabel] = frequencyUnitScale(opt.unit);
    frequency = calculateOrbitFrequencyField(orbitData, branchName, unitScale, meta);
    contourCount = getOptionValue(opt, 'contourCount', 0);
    baseValues = struct('sliceIndex', initialSliceIndex, 'contourCount', contourCount);
    controls = struct([]);
    if ~isempty(dynamicUpdate)
        controls = [ ...
            integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], ...
            initialSliceIndex, 1, phaseDimensionSize(dim, orbitData)), ...
            integerSliderControl('contourCount', 'contours', contourCount, 0, max(40, 2 * contourCount))];
    end
    figureName = sprintf('%s %s orbit frequency', speciesLabel, branchName);
    buildPlotData = @(values) buildOrbitFrequencyPlotData(orbitData, speciesLabel, ...
        branchName, dim, frequency, colorbarLabel, opt, ...
        applyControlValues(baseValues, values), ~isempty(dynamicUpdate));
    runPICPlotMode(figureName, controls, dynamicUpdate, buildPlotData, ...
        @renderMap, sprintf('%s 在该切片中没有有限值数据。', figureName));
end

function plotData = buildOrbitFrequencyPlotData(orbitData, speciesLabel, branchName, ...
    dim, frequency, colorbarLabel, opt, values, showControlStatus)

    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
        dim, values.sliceIndex, frequency, orbitData);
    if any(isfinite(Z(:)))
        Z = fillEnclosedBlankRegions(Z);
    end
    statusText = sprintf('%s %s orbit frequency', speciesLabel, branchName);
    if showControlStatus
        statusText = sprintf('%s, %s index = %d, contours = %d', ...
            statusText, char(opt.fixedCoordinate), values.sliceIndex, values.contourCount);
    end
    frequencyLatex = sprintf('\\mathrm{%s}\\ \\mathrm{orbit}\\ \\mathrm{frequency}', branchName);
    titleText = simpleTitleText(speciesLabel, frequencyLatex, sliceTitle, '');
    plotData = buildMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        colorbarLabel, opt.colormapIndex, values.contourCount, statusText);
end

function plotResonanceLine(phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [orbitData, speciesLabel] = resolveSpeciesData(phaseSpaceOrbit, opt.species, 'orbit');
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    initialSliceIndex = parsePhaseSlice(opt.slice, dim, orbitData);
    initialResonanceOptions = normalizeResonanceOptions(opt);
    baseValues = struct('sliceIndex', initialSliceIndex);
    controls = struct([]);
    if ~isempty(dynamicUpdate)
        controls = integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], ...
            initialSliceIndex, 1, phaseDimensionSize(dim, orbitData));
        controls = appendResonanceControls(controls, initialResonanceOptions);
    end
    figureName = sprintf('%s %s resonance residual', speciesLabel, initialResonanceOptions.branch);
    buildPlotData = @(values) buildResonancePlotData(orbitData, speciesLabel, dim, ...
        initialResonanceOptions, meta, opt, applyControlValues(baseValues, values));
    runPICPlotMode(figureName, controls, dynamicUpdate, buildPlotData, ...
        @renderMap, ...
        sprintf('%s %s resonance 在该切片中没有有限值数据。', speciesLabel, initialResonanceOptions.branch));
end

function plotData = buildResonancePlotData(orbitData, speciesLabel, dim, resOpt, meta, opt, values)

    resOpt = resonanceOptionsFromValues(resOpt, values);
    [residualHz, physicalN] = calculateResonanceResidualField(orbitData, resOpt, meta);
    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
        dim, values.sliceIndex, residualHz, orbitData);
    titleText = resonanceTitleText(speciesLabel, resOpt, physicalN, sliceTitle);
    plotData = buildMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        '$\Delta f_{\mathrm{res}}/\mathrm{Hz}$', opt.colormapIndex, 0, ...
        resonanceStatusText(speciesLabel, resOpt, physicalN, residualHasZeroContour(Z)));
    plotData.forceSigned = true;
    plotData.resonanceZ = Z;
end

function plotPhaseQuantity(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPhaseData, opt.species, 'phase');
    [quantityData, quantitySpec] = computePICQuantity(speciesData, opt.quantity);
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    initialSliceIndex = parsePhaseSlice(opt.slice, dim, speciesData);
    initialTimeIndex = getOptionValue(opt, 'timeIndex', 1);
    if quantitySpec.timeDependent
        initialTimeIndex = parseTimeIndex(initialTimeIndex, timeDimensionSize(quantityData));
    end
    contourCount = getOptionValue(opt, 'contourCount', 0);
    baseValues = struct('sliceIndex', initialSliceIndex, 'timeIndex', initialTimeIndex, 'contourCount', contourCount);
    controls = struct([]);
    if ~isempty(dynamicUpdate)
        controls = integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], ...
            initialSliceIndex, 1, phaseDimensionSize(dim, speciesData));
        if quantitySpec.timeDependent
            controls = [controls, integerSliderControl('timeIndex', 'timeIndex', ...
                initialTimeIndex, 1, timeDimensionSize(quantityData))];
        end
        controls = [controls, integerSliderControl('contourCount', 'contours', ...
            contourCount, 0, max(20, 2 * contourCount))];
        controls = appendResonanceControls(controls, opt.resonance, true, true);
    end
    figureName = sprintf('%s %s phase-space slice', speciesLabel, quantitySpec.field);
    buildPlotData = @(values) buildPhaseQuantityPlotData(speciesData, speciesLabel, ...
        quantityData, quantitySpec, dim, phaseSpaceOrbit, ...
        meta, opt, applyControlValues(baseValues, values), ~isempty(dynamicUpdate));
    runPICPlotMode(figureName, controls, dynamicUpdate, buildPlotData, ...
        @renderMap, ...
        sprintf('%s %s 在该切片中没有有限值数据。', speciesLabel, quantitySpec.field));
end

function plotData = buildPhaseQuantityPlotData(speciesData, speciesLabel, quantityData, ...
    quantitySpec, dim, phaseSpaceOrbit, meta, opt, values, isInteractive)

    [fieldData, timeIndexText] = selectQuantityFrame( ...
        quantityData, quantitySpec.timeDependent, values.timeIndex);
    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
        dim, values.sliceIndex, fieldData, speciesData);
    resOpt = opt.resonance;
    if isInteractive || any(isfinite(Z(:)))
        resOpt = resonanceOptionsFromValues(resOpt, values, true);
    end
    titleText = phaseQuantityTitleText( ...
        speciesLabel, quantitySpec.latex, sliceTitle, timeIndexText, resOpt);
    statusText = sprintf('%s %s, %s index = %d%s', ...
        speciesLabel, quantitySpec.field, char(opt.fixedCoordinate), values.sliceIndex, timeIndexText);
    if isInteractive
        statusText = sprintf('%s, contours = %d', statusText, values.contourCount);
    end
    plotData = buildMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantitySpec.latex '$'], opt.colormapIndex, values.contourCount, statusText);
    if isInteractive || any(isfinite(Z(:)))
        plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, values.sliceIndex, resOpt, meta);
        plotData = attachPhaseDetuningLine(plotData, dim, opt.detuningLine, resOpt, meta);
    end
end

function plotPhasePower(picPhaseData, phaseSpaceOrbit, meta, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPhaseData, opt.species, 'phase');
    powerData = requireQuantityData(speciesData, 'Power');
    dim = phaseCoordinateToDimension(opt.fixedCoordinate);
    initialSliceIndex = parsePhaseSlice(opt.slice, dim, speciesData);
    [~, initialModeN] = parseModeN(opt.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    initialTimeIndex = parseTimeIndex(opt.timeIndex, size(powerData, 5));
    contourCount = getOptionValue(opt, 'contourCount', 0);
    baseValues = struct('sliceIndex', initialSliceIndex, 'modeN', initialModeN, ...
        'timeIndex', initialTimeIndex, 'contourCount', contourCount);
    controls = struct([]);
    if ~isempty(dynamicUpdate)
        controls = integerSliderControl('sliceIndex', [char(opt.fixedCoordinate) ' index'], ...
            initialSliceIndex, 1, phaseDimensionSize(dim, speciesData));
        if numel(speciesData.modeIndexAll) > 1
            modeControl = integerSliderControl('modeN', 'n', initialModeN, ...
                min(speciesData.physicalNAll), max(speciesData.physicalNAll));
            modeControl.allowedValues = speciesData.physicalNAll;
            controls = [controls, modeControl];
        end
        if size(powerData, 5) > 1
            controls = [controls, integerSliderControl('timeIndex', 'timeIndex', ...
                initialTimeIndex, 1, size(powerData, 5))];
        end
        controls = [controls, integerSliderControl('contourCount', 'contours', ...
            contourCount, 0, max(20, 2 * contourCount))];
        controls = appendResonanceControls(controls, opt.resonance, true);
    end
    figureName = sprintf('%s PhasePower phase-space slice', speciesLabel);
    buildPlotData = @(values) buildPhasePowerPlotData(speciesData, speciesLabel, powerData, ...
        dim, phaseSpaceOrbit, meta, opt, applyControlValues(baseValues, values), ~isempty(dynamicUpdate));
    runPICPlotMode(figureName, controls, dynamicUpdate, buildPlotData, ...
        @renderMap, sprintf('%s PhasePower 在该切片中没有有限值数据。', speciesLabel));
end

function plotData = buildPhasePowerPlotData(speciesData, speciesLabel, powerData, ...
    dim, phaseSpaceOrbit, meta, opt, values, isInteractive)

    [modeIndex, modeN] = parseModeN(values.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    timeIndex = parseTimeIndex(values.timeIndex, size(powerData, 5));
    [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = slicePhaseField( ...
        dim, values.sliceIndex, powerData(:, :, :, modeIndex, timeIndex), speciesData);
    resOpt = opt.resonance;
    if isInteractive || any(isfinite(Z(:)))
        resOpt = resonanceOptionsFromValues(resOpt, values);
    end
    quantityLatex = 'P_{\mathrm{PIC}}';
    timeIndexText = sprintf(', timeIndex = %d', timeIndex);
    titleText = phasePowerTitleText(speciesLabel, quantityLatex, sliceTitle, modeN, timeIndexText, resOpt);
    statusText = sprintf('%s PhasePower, %s index = %d, n = %d%s', ...
        speciesLabel, char(opt.fixedCoordinate), values.sliceIndex, modeN, timeIndexText);
    if isInteractive
        statusText = sprintf('%s, contours = %d', statusText, values.contourCount);
    end
    plotData = buildMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantityLatex '$'], opt.colormapIndex, values.contourCount, statusText);
    if isInteractive || any(isfinite(Z(:)))
        plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, values.sliceIndex, resOpt, meta);
    end
end

function plotPitchQuantity(picPitchData, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPitchData, opt.species, 'pitch');
    [quantityData, quantitySpec] = computePICQuantity(speciesData, opt.quantity);
    contourCount = getOptionValue(opt, 'contourCount', 0);
    timeIndex = getOptionValue(opt, 'timeIndex', 1);
    if quantitySpec.timeDependent
        timeIndex = parseTimeIndex(timeIndex, timeDimensionSize(quantityData));
    end
    baseValues = struct('timeIndex', timeIndex, 'contourCount', contourCount);
    controls = struct([]);
    if ~isempty(dynamicUpdate)
        if quantitySpec.timeDependent
            controls = integerSliderControl('timeIndex', 'timeIndex', ...
                timeIndex, 1, timeDimensionSize(quantityData));
        end
        controls = [controls, integerSliderControl('contourCount', 'contours', ...
            contourCount, 0, max(20, 2 * contourCount))];
    end

    figureName = sprintf('%s %s pitch-space map', speciesLabel, quantitySpec.field);
    buildPlotData = @(values) buildPitchQuantityPlotData(speciesData, speciesLabel, ...
        quantityData, quantitySpec, opt, applyControlValues(baseValues, values), ~isempty(dynamicUpdate));
    runPICPlotMode(figureName, controls, dynamicUpdate, buildPlotData, ...
        @renderMap, ...
        sprintf('%s pitch %s 没有有限值数据。', speciesLabel, quantitySpec.field));
end

function plotData = buildPitchQuantityPlotData(speciesData, speciesLabel, quantityData, ...
    quantitySpec, opt, values, isInteractive)

    [fieldData, timeIndexText] = selectQuantityFrame( ...
        quantityData, quantitySpec.timeDependent, values.timeIndex);
    [Z, xVec, yVec, xlabelText, ylabelText] = pitchMapField(fieldData, speciesData);
    titleText = sprintf('$\\mathrm{%s}\\quad %s\\quad \\mathrm{pitch}%s$', ...
        speciesLabel, quantitySpec.latex, timeTextToLatex(timeIndexText));
    statusText = sprintf('%s pitch %s%s', speciesLabel, quantitySpec.field, timeIndexText);
    if isInteractive
        statusText = sprintf('%s, contours = %d', statusText, values.contourCount);
    end
    plotData = buildMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantitySpec.latex '$'], opt.colormapIndex, values.contourCount, statusText);
end

function plotPitchPower(picPitchData, opt, dynamicUpdate)

    [speciesData, speciesLabel] = resolveSpeciesData(picPitchData, opt.species, 'pitch');
    powerData = requireQuantityData(speciesData, 'Power');
    [~, initialModeN] = parseModeN(opt.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    initialTimeIndex = parseTimeIndex(opt.timeIndex, size(powerData, 4));
    contourCount = getOptionValue(opt, 'contourCount', 0);
    controls = struct([]);
    if ~isempty(dynamicUpdate)
        if numel(speciesData.modeIndexAll) > 1
            modeControl = integerSliderControl('modeN', 'n', initialModeN, ...
                min(speciesData.physicalNAll), max(speciesData.physicalNAll));
            modeControl.allowedValues = speciesData.physicalNAll;
            controls = [controls, modeControl];
        end
        if size(powerData, 4) > 1
            controls = [controls, integerSliderControl('timeIndex', 'timeIndex', ...
                initialTimeIndex, 1, size(powerData, 4))];
        end
        controls = [controls, integerSliderControl('contourCount', 'contours', ...
            contourCount, 0, max(20, 2 * contourCount))];
    end

    baseValues = struct('modeN', initialModeN, 'timeIndex', initialTimeIndex, 'contourCount', contourCount);
    figureName = sprintf('%s PitchPower pitch-space map', speciesLabel);
    buildPlotData = @(values) buildPitchPowerPlotData(speciesData, speciesLabel, ...
        powerData, opt, applyControlValues(baseValues, values), ~isempty(dynamicUpdate));
    runPICPlotMode(figureName, controls, dynamicUpdate, buildPlotData, ...
        @renderMap, sprintf('%s PitchPower 没有有限值数据。', speciesLabel));
end

function plotData = buildPitchPowerPlotData(speciesData, speciesLabel, powerData, opt, values, isInteractive)

    [modeIndex, modeN] = parseModeN(values.modeN, speciesData.modeIndexAll, speciesData.physicalNAll);
    timeIndex = parseTimeIndex(values.timeIndex, size(powerData, 4));
    [Z, xVec, yVec, xlabelText, ylabelText] = pitchMapField(powerData(:, :, modeIndex, timeIndex), speciesData);
    quantityLatex = 'P_{\mathrm{PIC}}';
    timeIndexText = sprintf(', timeIndex = %d', timeIndex);
    titleText = sprintf('$\\mathrm{%s}\\quad %s\\quad \\mathrm{pitch},\\quad n = %d%s$', ...
        speciesLabel, quantityLatex, modeN, timeTextToLatex(timeIndexText));
    statusText = sprintf('%s PitchPower, n = %d%s', speciesLabel, modeN, timeIndexText);
    if isInteractive
        statusText = sprintf('%s, contours = %d', statusText, values.contourCount);
    end
    plotData = buildMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        ['$' quantityLatex '$'], opt.colormapIndex, values.contourCount, statusText);
end

function plotDiffusivity(picDiffusivityData, opt, dynamicUpdate)

    plotType = requireIntegerInRange(opt.plotType, 1, 3, 'plotType');
    [speciesData, speciesLabel] = resolveSpeciesData(picDiffusivityData, opt.species, 'diffusivity');
    D = requireQuantityData(speciesData, 'Diffusivity');
    baseOpt = opt;
    controls = struct([]);
    switch plotType
        case 1
            plotTypeName = 'radial';
            renderPlotData = @renderLine;
            noFiniteMessage = sprintf('%s Diffusivity 径向剖面没有有限值数据。', speciesLabel);
            if ~isempty(dynamicUpdate)
                baseOpt.nRange = initialDiffusivityNRange(speciesData, opt.nRange);
                baseOpt.timeIndex = parseTimeIndex(opt.timeIndex, size(D, 1));
                controls = diffusivityRangeControls(speciesData, baseOpt.nRange);
                controls = appendIndexControl(controls, 'timeIndex', baseOpt.timeIndex, size(D, 1));
            end
            buildPlotData = @(values) buildDiffusivityRadialPlotData(speciesData, speciesLabel, ...
                diffusivityOptionsFromValues(baseOpt, values));
        case 2
            plotTypeName = 'time';
            renderPlotData = @renderLine;
            noFiniteMessage = sprintf('%s Diffusivity 时间曲线没有有限值数据。', speciesLabel);
            if ~isempty(dynamicUpdate)
                baseOpt.nRange = initialDiffusivityNRange(speciesData, opt.nRange);
                baseOpt.radialIndex = parseIndex(opt.radialIndex, size(D, 3), 'radial');
                controls = diffusivityRangeControls(speciesData, baseOpt.nRange);
                controls = appendIndexControl(controls, 'radialIndex', baseOpt.radialIndex, size(D, 3));
            end
            buildPlotData = @(values) buildDiffusivityTimePlotData(speciesData, speciesLabel, ...
                diffusivityOptionsFromValues(baseOpt, values));
        case 3
            plotTypeName = 'map';
            renderPlotData = @renderMap;
            noFiniteMessage = sprintf('%s Diffusivity 二维图没有有限值数据。', speciesLabel);
            if ~isempty(dynamicUpdate)
                [~, nAllowed] = diffusivityModeMask(speciesData, opt.nRange);
                nMin = min(nAllowed);
                nMax = max(nAllowed);
                baseOpt.nRange = diffusivitySingleN(speciesData, opt.nRange);
                if nMin < nMax
                    controls = integerSliderControl('nRange', 'n', baseOpt.nRange, nMin, nMax);
                    controls.allowedValues = nAllowed;
                end
            end
            buildPlotData = @(values) buildDiffusivityMapPlotData(speciesData, speciesLabel, ...
                applyControlValues(baseOpt, values));
    end
    runPICPlotMode(sprintf('%s Diffusivity %s', speciesLabel, plotTypeName), ...
        controls, dynamicUpdate, buildPlotData, renderPlotData, noFiniteMessage);
end

function plotData = buildDiffusivityRadialPlotData(speciesData, speciesLabel, opt)

    [Dsum, selectedN] = diffusivityModeSum(speciesData, opt.nRange);
    timeIndex = parseTimeIndex(opt.timeIndex, size(Dsum, 1));
    [xVec, xlabelText] = diagnosticRadialAxis(speciesData, opt.radialAxis);
    titleText = diffusivityTitleText(speciesLabel, selectedN, ...
        sprintf(',\\quad \\mathrm{timeIndex} = %d', timeIndex));
    statusText = sprintf('%s Diffusivity radial, n = [%d, %d], timeIndex = %d', ...
        speciesLabel, min(selectedN), max(selectedN), timeIndex);
    plotData = buildLinePlotData(xVec, Dsum(timeIndex, :), xlabelText, ...
        '$D/(\mathrm{m}^{2}/\mathrm{s})$', titleText, statusText);
end

function plotData = buildDiffusivityTimePlotData(speciesData, speciesLabel, opt)

    [Dsum, selectedN] = diffusivityModeSum(speciesData, opt.nRange);
    radialIndex = parseIndex(opt.radialIndex, size(Dsum, 2), 'radial');
    [xVec, xlabelText] = diagnosticTimeAxis(speciesData, opt.timeAxis);
    [radialVec, ~, radialName] = diagnosticRadialAxis(speciesData, opt.radialAxis);
    titleText = diffusivityTitleText(speciesLabel, selectedN, ...
        sprintf(',\\quad %s = %.6g', radialName, radialVec(radialIndex)));
    statusText = sprintf('%s Diffusivity time, n = [%d, %d], radialIndex = %d', ...
        speciesLabel, min(selectedN), max(selectedN), radialIndex);
    plotData = buildLinePlotData(xVec, Dsum(:, radialIndex), xlabelText, ...
        '$D/(\mathrm{m}^{2}/\mathrm{s})$', titleText, statusText);
end

function plotData = buildDiffusivityMapPlotData(speciesData, speciesLabel, opt)

    [Dn, selectedN] = diffusivityModeSlice(speciesData, opt.nRange);
    [xVec, xlabelText] = diagnosticTimeAxis(speciesData, opt.timeAxis);
    [yVec, ylabelText] = diagnosticRadialAxis(speciesData, opt.radialAxis);
    titleText = diffusivityTitleText(speciesLabel, selectedN, '');
    statusText = sprintf('%s Diffusivity map, n = %d', speciesLabel, selectedN);
    plotData = buildMapPlotData(Dn.', xVec, yVec, xlabelText, ylabelText, titleText, ...
        '$D/(\mathrm{m}^{2}/\mathrm{s})$', opt.colormapIndex, 0, statusText);
end

function plotDensity(picDensityData, opt, dynamicUpdate)

    plotType = requireIntegerInRange(opt.plotType, 1, 3, 'plotType');
    [speciesData, speciesLabel] = resolveSpeciesData(picDensityData, opt.species, 'density');
    densityData = requireQuantityData(speciesData, 'Density');
    assert(ismatrix(densityData) && size(densityData, 2) == speciesData.gridNx, ...
        'Density 数据尺寸必须为 [time, radial]。');
    controls = struct([]);
    switch plotType
        case 1
            plotTypeName = 'radial';
            timeIndex = parseTimeIndex(opt.timeIndex, size(densityData, 1));
            controls = appendIndexControl(controls, 'timeIndex', timeIndex, ...
                size(densityData, 1), ~isempty(dynamicUpdate));
            buildPlotData = @(values) buildDensityRadialPlotData(speciesData, speciesLabel, ...
                densityData, opt, applyControlValues(struct('timeIndex', timeIndex), values), ...
                ~isempty(dynamicUpdate));
            renderPlotData = @renderLine;
            noFiniteMessage = sprintf('%s Density 径向剖面没有有限值数据。', speciesLabel);
        case 2
            plotTypeName = 'time';
            radialIndex = parseIndex(opt.radialIndex, size(densityData, 2), 'radial');
            controls = appendIndexControl(controls, 'radialIndex', radialIndex, ...
                size(densityData, 2), ~isempty(dynamicUpdate));
            baseOpt = opt;
            baseOpt.radialIndex = radialIndex;
            buildPlotData = @(values) buildDensityTimePlotData(speciesData, speciesLabel, ...
                densityData, applyControlValues(baseOpt, values));
            renderPlotData = @renderLine;
            noFiniteMessage = sprintf('%s Density 时间曲线没有有限值数据。', speciesLabel);
        case 3
            plotTypeName = 'map';
            buildPlotData = @(ignored) buildDensityMapPlotData( ...
                speciesData, speciesLabel, densityData, opt);
            renderPlotData = @renderMap;
            noFiniteMessage = sprintf('%s Density 二维图没有有限值数据。', speciesLabel);
    end
    runPICPlotMode(sprintf('%s Density %s', speciesLabel, plotTypeName), ...
        controls, dynamicUpdate, buildPlotData, renderPlotData, noFiniteMessage);
end

function plotData = buildDensityRadialPlotData(speciesData, speciesLabel, densityData, opt, values, isInteractive)

    timeIndex = values.timeIndex;
    deltaDensity = densityData(timeIndex, :);
    if ~isInteractive && ~any(isfinite(deltaDensity(:)))
        plotData = buildLinePlotData(1:numel(deltaDensity), deltaDensity, '', '', '', '');
        return;
    end

    [xVec, xlabelText] = diagnosticRadialAxis(speciesData, opt.radialAxis);
    timeText = sprintf(',\\quad \\mathrm{timeIndex} = %d', timeIndex);

    if isTotalDensityMode(opt)
        initialDensity = requireInitialDensity(speciesData);
        totalDensity = initialDensity + reshape(deltaDensity, 1, []);
        titleText = sprintf('$\\mathrm{%s}\\quad n%s$', speciesLabel, timeText);
        statusText = sprintf('%s Density radial total, timeIndex = %d', speciesLabel, timeIndex);
        plotData = buildLinePlotData(xVec, initialDensity, xlabelText, ...
            '$n\,[\mathrm{m}^{-3}]$', titleText, statusText);
        plotData.lineLabel = '$n_0$';
        plotData.lineOverlays = struct('xVec', xVec, 'yVec', totalDensity, ...
            'label', '$n_0 + \delta n$');
        plotData.finiteData = deltaDensity;
        return;
    end

    titleText = densityTitleText(speciesLabel, timeText);
    statusText = sprintf('%s Density radial, timeIndex = %d', speciesLabel, timeIndex);
    plotData = buildLinePlotData(xVec, deltaDensity, xlabelText, '$\delta n$', titleText, statusText);
end

function plotData = buildDensityTimePlotData(speciesData, speciesLabel, densityData, opt)

    radialIndex = parseIndex(opt.radialIndex, size(densityData, 2), 'radial');
    [xVec, xlabelText] = diagnosticTimeAxis(speciesData, opt.timeAxis);
    [radialVec, ~, radialName] = diagnosticRadialAxis(speciesData, opt.radialAxis);
    titleText = densityTitleText(speciesLabel, ...
        sprintf(',\\quad %s = %.6g', radialName, radialVec(radialIndex)));
    statusText = sprintf('%s Density time, radialIndex = %d', speciesLabel, radialIndex);
    plotData = buildLinePlotData(xVec, densityData(:, radialIndex), xlabelText, ...
        '$\delta n$', titleText, statusText);
end

function plotData = buildDensityMapPlotData(speciesData, speciesLabel, densityData, opt)

    [xVec, xlabelText] = diagnosticTimeAxis(speciesData, opt.timeAxis);
    [yVec, ylabelText] = diagnosticRadialAxis(speciesData, opt.radialAxis);
    titleText = densityTitleText(speciesLabel, '');
    statusText = sprintf('%s Density map', speciesLabel);
    plotData = buildMapPlotData(densityData.', xVec, yVec, xlabelText, ylabelText, titleText, ...
        '$\delta n$', opt.colormapIndex, 0, statusText);
end

function plotData = attachResonanceOverlay(plotData, phaseSpaceOrbit, speciesLabel, dim, sliceIndex, resOpt, meta)

    plotData.resonanceZ = [];
    plotData.resonanceOverlays = struct('Z', {}, 'label', {}, 'harmonic', {});
    if ~isfield(resOpt, 'enabled') || ~resOpt.enabled
        return;
    end

    if ~hasSpeciesFields(phaseSpaceOrbit, speciesLabel, {})
        plotData.status = [plotData.status ', resonance = no orbit data'];
        return;
    end

    orbitData = phaseSpaceOrbit.(speciesLabel);
    resOpt = normalizeResonanceOverlayOptions(resOpt);
    branchNames = normalizeResonanceBranches(resOpt.branch);
    hasZero = false(1, 0);
    physicalN = resOpt.toroidalMode;
    overlayIndex = 0;
    for branchIndex = 1:numel(branchNames)
        harmonicValues = resOpt.harmonicMinList(branchIndex):resOpt.harmonicMaxList(branchIndex);
        for harmonicIndex = 1:numel(harmonicValues)
            overlayIndex = overlayIndex + 1;
            lineOpt = resOpt;
            lineOpt.branch = branchNames{branchIndex};
            lineOpt.harmonic = harmonicValues(harmonicIndex);
            [residualHz, physicalN] = calculateResonanceResidualField(orbitData, lineOpt, meta);
            [resonanceZ, ~, ~, ~, ~, ~] = slicePhaseField(dim, sliceIndex, residualHz, orbitData);
            plotData.resonanceOverlays(overlayIndex).Z = resonanceZ;
            plotData.resonanceOverlays(overlayIndex).label = resonanceOverlayLabel( ...
                branchNames{branchIndex}, harmonicValues(harmonicIndex), numel(branchNames));
            plotData.resonanceOverlays(overlayIndex).harmonic = harmonicValues(harmonicIndex);
            hasZero(overlayIndex) = residualHasZeroContour(resonanceZ);
        end
    end
    plotData.status = [plotData.status ', ' resonanceOverlayStatusText('', resOpt, physicalN, branchNames, hasZero)];
end

function plotData = attachPhaseDetuningLine(plotData, dim, lineOpt, resOpt, meta)

    plotData.detuningLine = [];
    if ~isfield(lineOpt, 'enabled') || ~lineOpt.enabled
        return;
    end
    if dim ~= 3
        plotData.status = [plotData.status ', detuning line skipped: fixedCoordinate must be Lambda'];
        return;
    end

    lineOpt = normalizePhaseDetuningLineOptions(lineOpt);
    frequencyHz = requireFiniteScalarInRange(resOpt.frequencyHz, -Inf, Inf, 'resonance.frequencyHz');
    toroidalMode = requireIntegerInRange(resOpt.toroidalMode, 1, Inf, 'resonance.toroidalMode');
    kappa = resonanceDetuningKappa(meta, frequencyHz, toroidalMode);

    deltaPphiPlot = [-lineOpt.deltaPphiLeft, lineOpt.deltaPphiRight];
    xPlot = lineOpt.Pphi0 + deltaPphiPlot;
    deltaPphiRaw = -deltaPphiPlot;
    yE = lineOpt.E0 + kappa * deltaPphiRaw;

    plotData.detuningLine = struct( ...
        'x', xPlot, ...
        'y', yE, ...
        'label', 'detuning line');
    plotData.status = [plotData.status ', detuning line = on'];
end

function lineOpt = normalizePhaseDetuningLineOptions(lineOpt)

    requiredFields = {'E0', 'Pphi0', 'deltaPphiLeft', 'deltaPphiRight'};
    for fieldIndex = 1:numel(requiredFields)
        fieldName = requiredFields{fieldIndex};
        assert(isfield(lineOpt, fieldName) && ~isempty(lineOpt.(fieldName)), ...
            'detuningLine.%s must be provided.', fieldName);
    end
    lineOpt.E0 = requireFiniteScalarInRange(lineOpt.E0, -Inf, Inf, 'detuningLine.E0');
    lineOpt.Pphi0 = requireFiniteScalarInRange(lineOpt.Pphi0, -Inf, Inf, 'detuningLine.Pphi0');
    lineOpt.deltaPphiLeft = requireFiniteScalarInRange( ...
        lineOpt.deltaPphiLeft, -Inf, Inf, 'detuningLine.deltaPphiLeft');
    lineOpt.deltaPphiRight = requireFiniteScalarInRange( ...
        lineOpt.deltaPphiRight, -Inf, Inf, 'detuningLine.deltaPphiRight');
    assert(lineOpt.deltaPphiLeft > 0 && lineOpt.deltaPphiRight > 0, ...
        'detuningLine.deltaPphiLeft/deltaPphiRight must be positive.');
end

function [kappa, unitRatio, omegaOverN] = resonanceDetuningKappa(meta, frequencyHz, toroidalMode)

    assert(all(isfinite([meta.QE, meta.B0, meta.L0, meta.MP, meta.VA0])) && ...
        all([meta.QE, meta.B0, meta.L0, meta.MP, meta.VA0] > 0), ...
        'detuning diagnostics need positive QE/B0/L0/MP/VA0 in normalization2D.mat.');
    unitRatio = meta.QE * meta.B0 * meta.L0^2 / (meta.MP * meta.VA0^2);
    omegaOverN = 2 * pi * frequencyHz / toroidalMode;
    kappa = unitRatio * omegaOverN;
end

function [quantityData, quantitySpec] = computePICQuantity(speciesData, quantityName)

    quantitySpec = resolvePICQuantitySpec(quantityName);
    if ~quantitySpec.supported
        error('quantity 必须为 "J"、"F0"、"DF"、"f0"、"df" 或 "df/f0"。');
    end

    numerator = requireQuantityData(speciesData, quantitySpec.requiredFields{1});
    quantityData = numerator;
    if numel(quantitySpec.requiredFields) > 1
        denominatorName = quantitySpec.requiredFields{2};
        denominator = requireQuantityData(speciesData, denominatorName);
        quantityData = divideQuantity( ...
            numerator, denominator, denominatorName, 1e-6, quantitySpec.floorMode);
    end
end

function quantitySpec = resolvePICQuantitySpec(quantityName)

    quantityText = strtrim(char(quantityName));
    normalizedName = lower(strrep(strrep(quantityText, '_', ''), ' ', ''));
    definitions = { ...
        'J',     '\mathcal{J}',      false, {'J'},       '',         {},     ...
            {'j', 'jacobian', 'phasespacejacobian', 'pitchspacejacobian'}; ...
        'F0',    'F_0',              false, {'F0'},      '',         {'F0'}, ...
            {'phasespacef0', 'pitchspacef0'}; ...
        'DF',    '\delta F',         true,  {'DF'},      '',         {'DF'}, ...
            {'deltaf', 'phasedeltaf', 'pitchdeltaf'}; ...
        'f0',    'f_0',              false, {'F0', 'J'}, 'absolute', {'f0'}, {}; ...
        'df',    '\delta f',         true,  {'DF', 'J'}, 'absolute', {'df'}, {}; ...
        'df/f0', '\delta f / f_0',   true,  {'DF', 'F0'}, 'relative', {}, ...
            {'df/f0', 'dff0', 'df2f0', 'dfoverf0'}};

    quantitySpec = struct( ...
        'field', quantityText, ...
        'latex', '', ...
        'timeDependent', false, ...
        'requiredFields', {{quantityText}}, ...
        'floorMode', '', ...
        'supported', false);
    for quantityDefinitionIndex = 1:size(definitions, 1)
        exactNames = definitions{quantityDefinitionIndex, 6};
        normalizedNames = definitions{quantityDefinitionIndex, 7};
        if ismember(quantityText, exactNames) || ismember(normalizedName, normalizedNames)
            quantitySpec.field = definitions{quantityDefinitionIndex, 1};
            quantitySpec.latex = definitions{quantityDefinitionIndex, 2};
            quantitySpec.timeDependent = definitions{quantityDefinitionIndex, 3};
            quantitySpec.requiredFields = definitions{quantityDefinitionIndex, 4};
            quantitySpec.floorMode = definitions{quantityDefinitionIndex, 5};
            quantitySpec.supported = true;
            return;
        end
    end
end

function data = requireQuantityData(speciesData, quantityField)

    if ~isfield(speciesData, quantityField) || isempty(speciesData.(quantityField))
        error('%s %s 尚未读取。请检查文件路径。', speciesData.species, quantityField);
    end
    data = speciesData.(quantityField);
end

function initialDensity = requireInitialDensity(speciesData)

    if ~isfield(speciesData, 'InitialDensity') || isempty(speciesData.InitialDensity)
        if isfield(speciesData, 'InitialDensityError') && ~isempty(speciesData.InitialDensityError)
            error('%s total 密度径向图无法使用 NTP 初始密度：%s', ...
                speciesData.species, speciesData.InitialDensityError);
        end
        densityField = ntpDensityFieldName(speciesData.species);
        error(['%s total 密度径向图需要 NTP.mat 中的 rhoSample 和 %s。' ...
            '请确认 inputDir 下存在 NTP.mat，且该物种初始密度字段可用。'], ...
            speciesData.species, densityField);
    end

    initialDensity = reshape(speciesData.InitialDensity, 1, []);
    assert(numel(initialDensity) == speciesData.gridNx, ...
        '%s InitialDensity 长度必须等于 gridNx。', speciesData.species);
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

function [fieldData, timeIndexText] = selectQuantityFrame(quantityData, isTimeDependent, requestedTimeIndex)

    if ~isTimeDependent
        fieldData = quantityData;
        timeIndexText = '';
        return;
    end

    nTime = timeDimensionSize(quantityData);
    selectedTimeIndex = parseTimeIndex(requestedTimeIndex, nTime);
    if ndims(quantityData) == 4
        fieldData = quantityData(:, :, :, selectedTimeIndex);
    elseif ndims(quantityData) == 3
        fieldData = quantityData(:, :, selectedTimeIndex);
    else
        error('含时间的诊断量必须为 3D 或 4D 数组。');
    end
    timeIndexText = sprintf(', timeIndex = %d', selectedTimeIndex);
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

function opt = normalizeResonanceDetuningOptions(opt, meta)

    requiredFields = {'species', 'plotFile', 'branch', 'frequencyHz', 'toroidalMode', ...
        'poloidalMode', 'harmonic', 'E0', 'Pphi0', 'Lambda0', ...
        'deltaPphiLeft', 'deltaPphiRight', 'sampleCount'};
    for fieldIndex = 1:numel(requiredFields)
        fieldName = requiredFields{fieldIndex};
        assert(isfield(opt, fieldName) && ~isempty(opt.(fieldName)), ...
            'picResonanceDetuningOpt.%s must be provided.', fieldName);
    end

    opt.branch = normalizeResonanceBranch(opt.branch);
    opt.frequencyHz = requireFiniteScalarInRange(opt.frequencyHz, -Inf, Inf, 'frequencyHz');
    opt.toroidalMode = requireIntegerInRange(opt.toroidalMode, 1, Inf, 'toroidalMode');
    opt.harmonic = requireIntegerInRange(opt.harmonic, -Inf, Inf, 'harmonic');
    if strcmp(opt.branch, 'trapped')
        opt.poloidalMode = 0;
    else
        opt.poloidalMode = requireIntegerInRange(opt.poloidalMode, 0, Inf, 'poloidalMode');
    end

    opt.E0 = requireFiniteScalarInRange(opt.E0, -Inf, Inf, 'E0');
    opt.PphiPlot0 = requireFiniteScalarInRange(opt.Pphi0, -Inf, Inf, 'Pphi0');
    opt.PphiRaw0 = -opt.PphiPlot0;
    opt.Lambda0 = requireFiniteScalarInRange(opt.Lambda0, -Inf, Inf, 'Lambda0');
    opt.deltaPphiLeft = requireFiniteScalarInRange(opt.deltaPphiLeft, -Inf, Inf, 'deltaPphiLeft');
    opt.deltaPphiRight = requireFiniteScalarInRange(opt.deltaPphiRight, -Inf, Inf, 'deltaPphiRight');
    assert(opt.deltaPphiLeft > 0 && opt.deltaPphiRight > 0, ...
        'deltaPphiLeft and deltaPphiRight must be positive.');
    opt.sampleCount = requireIntegerInRange(opt.sampleCount, 1, Inf, 'sampleCount');
    assert(opt.sampleCount >= 5, 'sampleCount must be at least 5.');

    [opt.kappa, unitRatio, omegaOverN] = resonanceDetuningKappa(meta, opt.frequencyHz, opt.toroidalMode);
    fprintf(['[resonance detuning] kappa=dE/dPphi(raw)=%.16g, displayed slope dE/dPphi(plot)=%.16g.\n', ...
        '[resonance detuning] kappa = (QE*B0*L0^2)/(MP*VA0^2) * 2*pi*f/n = %.16g * %.16g\n'], ...
        opt.kappa, -opt.kappa, unitRatio, omegaOverN);
end

function detuningInput = readAllResonanceDetuningInputs(plotFile, inputDir, speciesList, meta)

    detuningInput = struct();
    detuningInput.plot = captureResonanceDetuningInput( ...
        plotFile, @() load(plotFile)); %#ok<LOAD>
    detuningInput.mapping = struct();

    nPhase = meta.gridE * meta.gridPphi * meta.gridLambda;
    for speciesIndex = 1:numel(speciesList)
        speciesName = char(speciesList{speciesIndex});
        mappingFile = fullfile(inputDir, [speciesName 'PhaseSpaceMapping.bin']);
        detuningInput.mapping.(speciesName) = captureResonanceDetuningInput( ...
            mappingFile, @() readResonanceDetuningMappingRaw(mappingFile, nPhase));
    end
end

function capturedInput = captureResonanceDetuningInput(filePath, readerFcn)

    capturedInput = struct( ...
        'path', filePath, ...
        'exists', isfile(filePath), ...
        'data', [], ...
        'readException', []);
    if ~capturedInput.exists
        return;
    end

    try
        capturedInput.data = readerFcn();
    catch err
        capturedInput.readException = err;
    end
end

function data = requireCapturedResonanceDetuningInput(capturedInput)

    if ~isempty(capturedInput.readException)
        rethrow(capturedInput.readException);
    end
    data = capturedInput.data;
end

function plotInput = buildResonanceDetuningPlotInput(capturedPlot, requestedPlotFile, meta)

    plotFile = capturedPlot.path;
    assert(strcmp(char(requestedPlotFile), char(plotFile)) || ...
        (ispc && strcmpi(char(requestedPlotFile), char(plotFile))), ...
        'plotFile 已在“读取所有输入”阶段固定；修改后请重新运行该 section。');
    assert(capturedPlot.exists, '缺少 plot 文件：%s', plotFile);
    if isfield(meta, 'NFP') && isfinite(meta.NFP)
        assert(meta.NFP == 1, '共振失谐诊断当前只按 NFP=1 实现；当前 NFP=%g。', meta.NFP);
    end

    raw = requireCapturedResonanceDetuningInput(capturedPlot);
    assert(isfield(raw, 'rhoplot') && isfield(raw, 'Rplot'), ...
        '%s 必须包含 rhoplot 和 Rplot。', plotFile);

    Rplot = firstToroidalSlice(raw.Rplot);
    rhoplot = firstToroidalSlice(raw.rhoplot);
    assert(ismatrix(Rplot) && size(Rplot, 1) >= 2 && size(Rplot, 2) >= 1, ...
        'Rplot 必须是 [rho,theta] 或 [rho,theta,phi] 数组。');

    thetaIndex = max(1, round(size(Rplot, 2) / 2));
    RAxis = double(Rplot(:, thetaIndex));
    if isvector(rhoplot)
        rhoAxis = double(rhoplot(:));
    else
        assert(size(rhoplot, 1) == size(Rplot, 1), 'rhoplot 第一维必须匹配 Rplot 第一维。');
        rhoAxis = double(rhoplot(:, min(thetaIndex, size(rhoplot, 2))));
    end
    assert(numel(rhoAxis) == numel(RAxis), 'rhoplot 与 Rplot(:,thetaIndex) 长度不一致。');

    valid = isfinite(rhoAxis) & isfinite(RAxis);
    rhoAxis = rhoAxis(valid);
    RAxis = RAxis(valid);
    [rhoAxis, order] = sort(rhoAxis);
    RAxis = RAxis(order);
    [rhoAxis, uniqueIndex] = unique(rhoAxis, 'stable');
    RAxis = RAxis(uniqueIndex);
    assert(numel(rhoAxis) >= 2, 'rhoplot/Rplot 有效径向点不足，无法建立 rho->R 映射。');

    plotInput = struct('plotFile', plotFile, 'thetaIndex', thetaIndex, ...
        'rhoAxis', rhoAxis, 'RAxis', RAxis);
    fprintf('[load] %s: resonance detuning uses Rplot(:,%d), rho=[%.6g, %.6g], R=[%.6g, %.6g].\n', ...
        plotFile, thetaIndex, min(rhoAxis), max(rhoAxis), min(RAxis), max(RAxis));
end

function data2D = firstToroidalSlice(data)

    data = double(data);
    if ndims(data) <= 2
        data2D = data;
    else
        data2D = data(:, :, 1);
    end
end

function rawMapping = readResonanceDetuningMappingRaw(filePath, nPhase)

    nRecord = 2 * nPhase;
    fid = fopen(filePath, 'rb');
    assert(fid >= 0, '无法打开文件：%s', filePath);
    cleanupObj = onCleanup(@() fclose(fid));

    ids = fread(fid, nRecord, 'int32=>int32');
    rho = fread(fid, nRecord, 'double=>double');
    vpara = fread(fid, nRecord, 'double=>double');
    mu = fread(fid, nRecord, 'double=>double');
    assert(numel(ids) == nRecord && numel(rho) == nRecord && numel(vpara) == nRecord && numel(mu) == nRecord, ...
        '%s 尺寸不匹配：期望每个数组长度为 2*gridE*gridPphi*gridLambda=%d。', filePath, nRecord);
    clear cleanupObj;

    rawMapping = struct( ...
        'file', filePath, ...
        'ids', ids, ...
        'rho', rho, ...
        'vpara', vpara, ...
        'mu', mu);
end

function mapping = buildResonanceDetuningMapping( ...
    capturedMappings, speciesLabel, speciesData, branchName, plotInput, meta)

    assert(isfield(capturedMappings, speciesLabel), ...
        '物种 %s 的 mapping 输入未在“读取所有输入”阶段登记。', speciesLabel);
    capturedMapping = capturedMappings.(speciesLabel);
    filePath = capturedMapping.path;
    assert(capturedMapping.exists, '缺少 mapping 文件：%s', filePath);
    rawMapping = requireCapturedResonanceDetuningInput(capturedMapping);

    nPhase = speciesData.gridE * speciesData.gridPphi * speciesData.gridLambda;
    ids = rawMapping.ids;
    rho = rawMapping.rho;
    vpara = rawMapping.vpara;
    mu = rawMapping.mu;

    [ids, rho, vpara, mu] = selectResonanceDetuningMappingBranch(ids, rho, vpara, mu, nPhase, branchName);
    valid = ~isPadRecord(ids) & isfinite(rho) & isfinite(vpara) & isfinite(mu);
    assert(all(rho(valid) >= 0 & rho(valid) <= 1), ...
        'mapping 文件中的 rho 应为 (rho-RHO0)/(RHO1-RHO0) 归一化坐标，有效值必须位于 [0,1]。');
    rho(~valid) = NaN;
    vpara(~valid) = NaN;
    mu(~valid) = NaN;

    physicalRho = meta.RHO0 + double(rho) * (meta.RHO1 - meta.RHO0);
    minorR = interp1(plotInput.rhoAxis, plotInput.RAxis, physicalRho, 'linear', NaN) - meta.L0;
    physicalRho(~valid) = NaN;
    minorR(~valid) = NaN;

    mapping = struct( ...
        'file', filePath, ...
        'branch', branchName, ...
        'id', phaseVectorToGrid(double(ids), speciesData), ...
        'rhoFile', phaseVectorToGrid(rho, speciesData), ...
        'rho', phaseVectorToGrid(physicalRho, speciesData), ...
        'vpara', phaseVectorToGrid(vpara, speciesData), ...
        'mu', phaseVectorToGrid(mu, speciesData), ...
        'r', phaseVectorToGrid(minorR, speciesData), ...
        'valid', phaseVectorToGrid(valid, speciesData));
    fprintf('[load] %s: branch=%s, valid=%d/%d.\n', filePath, branchName, nnz(valid), nPhase);
end

function [ids, rho, vpara, mu] = selectResonanceDetuningMappingBranch(idsAll, rhoAll, vparaAll, muAll, nPhase, branchName)

    paraIndex = 1:nPhase;
    antiIndex = nPhase + (1:nPhase);
    switch branchName
        case 'para'
            ids = idsAll(paraIndex);
            rho = rhoAll(paraIndex);
            vpara = vparaAll(paraIndex);
            mu = muAll(paraIndex);
        case 'anti'
            ids = idsAll(antiIndex);
            rho = rhoAll(antiIndex);
            vpara = vparaAll(antiIndex);
            mu = muAll(antiIndex);
        case 'trapped'
            idsPara = idsAll(paraIndex);
            idsAnti = idsAll(antiIndex);
            rhoPara = rhoAll(paraIndex);
            rhoAnti = rhoAll(antiIndex);
            validPara = ~isPadRecord(idsPara) & isfinite(rhoPara) & isfinite(muAll(paraIndex));
            validAnti = ~isPadRecord(idsAnti) & isfinite(rhoAnti) & isfinite(muAll(antiIndex));
            useAnti = (~validPara & validAnti) | (validPara & validAnti & rhoAnti > rhoPara);
            ids = idsPara;
            rho = rhoPara;
            vpara = vparaAll(paraIndex);
            mu = muAll(paraIndex);
            ids(useAnti) = abs(idsAnti(useAnti));
            rho(useAnti) = rhoAnti(useAnti);
            vpara(useAnti) = vparaAll(antiIndex(useAnti));
            mu(useAnti) = muAll(antiIndex(useAnti));
        otherwise
            error('未知 mapping 分支：%s。', branchName);
    end
end

function gridData = phaseVectorToGrid(vectorData, speciesData)

    gridData = reshape(vectorData, [speciesData.gridLambda, speciesData.gridPphi, speciesData.gridE]);
    gridData = permute(gridData, [3, 2, 1]);
end

function gridPoint = nearestResonanceDetuningGridPoint(speciesData, mapping, residualHz, opt)

    [~, eIndex] = min(abs(speciesData.E1d - opt.E0));
    [~, pphiIndex] = min(abs(speciesData.Pphi1d - opt.PphiRaw0));
    [~, lambdaIndex] = min(abs(speciesData.Lambda1d - opt.Lambda0));

    gridPoint = resonanceDetuningPoint( ...
        speciesData.E1d(eIndex), speciesData.Pphi1d(pphiIndex), speciesData.Lambda1d(lambdaIndex), ...
        mapping.mu(eIndex, pphiIndex, lambdaIndex), mapping.rho(eIndex, pphiIndex, lambdaIndex), ...
        mapping.r(eIndex, pphiIndex, lambdaIndex), residualHz(eIndex, pphiIndex, lambdaIndex));
    gridPoint.index = [eIndex, pphiIndex, lambdaIndex];

    fprintf(['[resonance detuning] nearest grid point:\n', ...
        '  input E=%.16g, Pphi(plot)=%.16g, Pphi(raw)=%.16g, Lambda=%.16g\n', ...
        '  grid  E=%.16g, Pphi(plot)=%.16g, Pphi(raw)=%.16g, Lambda=%.16g, mu=%.16g, rho=%.16g, r=%.16g, R=%.6g Hz, index=[%d %d %d]\n', ...
        '  delta dE=%.6g, dPphi(plot)=%.6g, dPphi(raw)=%.6g, dLambda=%.6g\n'], ...
        opt.E0, opt.PphiPlot0, opt.PphiRaw0, opt.Lambda0, ...
        gridPoint.E, -gridPoint.Pphi, gridPoint.Pphi, gridPoint.Lambda, gridPoint.mu, gridPoint.rho, gridPoint.r, gridPoint.residualHz, ...
        eIndex, pphiIndex, lambdaIndex, gridPoint.E - opt.E0, -gridPoint.Pphi - opt.PphiPlot0, ...
        gridPoint.Pphi - opt.PphiRaw0, gridPoint.Lambda - opt.Lambda0);
end

function point = resonanceDetuningPoint(E, Pphi, Lambda, mu, rho, r, residualHz)

    point = struct('E', E, 'Pphi', Pphi, 'Lambda', Lambda, 'mu', mu, ...
        'rho', rho, 'r', r, 'residualHz', residualHz);
end

function center = projectResonanceDetuningCenter(speciesData, mapping, residualHz, gridPoint)

    finiteResidual = residualHz(isfinite(residualHz));
    assert(~isempty(finiteResidual), 'R 场没有有限值，无法投影代表点。');
    assert(min(finiteResidual) <= 0 && max(finiteResidual) >= 0, ...
        '三维相空间中没有 R=0 等值面，无法投影代表点。');

    [PphiGrid, EGrid, LambdaGrid] = meshgrid(speciesData.Pphi1d, speciesData.E1d, speciesData.Lambda1d);
    surfaceData = isosurface(PphiGrid, EGrid, LambdaGrid, residualHz, 0);
    assert(isstruct(surfaceData) && isfield(surfaceData, 'vertices') && isfield(surfaceData, 'faces') && ...
        ~isempty(surfaceData.vertices) && ~isempty(surfaceData.faces), ...
        '三维相空间中没有提取到 R=0 等值面。');

    queryPEL = [gridPoint.Pphi, gridPoint.E, gridPoint.Lambda];
    [projectedPEL, distance, faceIndex] = nearestPointOnTriangulatedSurface( ...
        surfaceData.vertices, surfaceData.faces, queryPEL);

    residualInterpolant = phaseFieldInterpolant(speciesData, residualHz);
    rhoInterpolant = phaseFieldInterpolant(speciesData, mapping.rho);
    radiusInterpolant = phaseFieldInterpolant(speciesData, mapping.r);
    mappingMuInterpolant = phaseFieldInterpolant(speciesData, mapping.mu);

    Pphi = projectedPEL(1);
    E = projectedPEL(2);
    Lambda = projectedPEL(3);
    residualAtCenter = residualInterpolant(E, Pphi, Lambda);
    rho = rhoInterpolant(E, Pphi, Lambda);
    r = radiusInterpolant(E, Pphi, Lambda);
    muFromCoordinates = E * Lambda;
    muFromMapping = mappingMuInterpolant(E, Pphi, Lambda);

    center = resonanceDetuningPoint(E, Pphi, Lambda, muFromCoordinates, rho, r, residualAtCenter);
    center.muMapInterpolated = muFromMapping;
    center.projectionDistance = distance;
    center.projectionFaceIndex = faceIndex;

    fprintf(['[resonance detuning] projected R=0 reference:\n', ...
        '  center E=%.16g, Pphi(plot)=%.16g, Pphi(raw)=%.16g, Lambda=%.16g, mu=%.16g, rho=%.16g, r=%.16g, R=%.6g Hz\n', ...
        '  mapping-interpolated mu=%.16g, projectionDistance=%.6g, face=%d\n', ...
        '  delta from grid: dE=%.6g, dPphi(plot)=%.6g, dPphi(raw)=%.6g, dLambda=%.6g, dmu=%.6g, dr=%.6g, dR=%.6g Hz\n'], ...
        center.E, -center.Pphi, center.Pphi, center.Lambda, center.mu, center.rho, center.r, center.residualHz, ...
        center.muMapInterpolated, center.projectionDistance, center.projectionFaceIndex, ...
        center.E - gridPoint.E, -center.Pphi + gridPoint.Pphi, center.Pphi - gridPoint.Pphi, ...
        center.Lambda - gridPoint.Lambda, center.mu - gridPoint.mu, center.r - gridPoint.r, ...
        center.residualHz - gridPoint.residualHz);
end

function interpolant = phaseFieldInterpolant(speciesData, fieldData)

    interpolant = griddedInterpolant( ...
        {speciesData.E1d, speciesData.Pphi1d, speciesData.Lambda1d}, ...
        fieldData, 'linear', 'none');
end

function [point, distance, faceIndex] = nearestPointOnTriangulatedSurface(vertices, faces, queryPoint)

    bestDistanceSquared = Inf;
    point = vertices(1, :);
    faceIndex = 1;
    for candidateFaceIndex = 1:size(faces, 1)
        tri = vertices(faces(candidateFaceIndex, :), :);
        candidate = nearestPointOnTriangle(queryPoint, tri(1, :), tri(2, :), tri(3, :));
        distanceSquared = sum((candidate - queryPoint) .^ 2);
        if distanceSquared < bestDistanceSquared
            bestDistanceSquared = distanceSquared;
            point = candidate;
            faceIndex = candidateFaceIndex;
        end
    end
    distance = sqrt(bestDistanceSquared);
end

function point = nearestPointOnTriangle(p, a, b, c)

    ab = b - a;
    ac = c - a;
    ap = p - a;
    d1 = dot(ab, ap);
    d2 = dot(ac, ap);
    if d1 <= 0 && d2 <= 0
        point = a;
        return;
    end

    bp = p - b;
    d3 = dot(ab, bp);
    d4 = dot(ac, bp);
    if d3 >= 0 && d4 <= d3
        point = b;
        return;
    end

    vc = d1 * d4 - d3 * d2;
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1 / (d1 - d3);
        point = a + v * ab;
        return;
    end

    cp = p - c;
    d5 = dot(ab, cp);
    d6 = dot(ac, cp);
    if d6 >= 0 && d5 <= d6
        point = c;
        return;
    end

    vb = d5 * d2 - d1 * d6;
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2 / (d2 - d6);
        point = a + w * ac;
        return;
    end

    va = d3 * d6 - d5 * d4;
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        point = b + w * (c - b);
        return;
    end

    denominator = 1 / (va + vb + vc);
    v = vb * denominator;
    w = vc * denominator;
    point = a + ab * v + ac * w;
end

function path = sampleResonanceDetuningPath(speciesData, mapping, residualHz, center, opt)

    deltaPphiPlot = linspace(-opt.deltaPphiLeft, opt.deltaPphiRight, opt.sampleCount).';
    deltaPphiPlot = unique([deltaPphiPlot; 0]);
    PphiPlot = -center.Pphi + deltaPphiPlot;
    Pphi = -PphiPlot;
    deltaPphiRaw = Pphi - center.Pphi;
    E = center.E + opt.kappa * deltaPphiRaw;
    Lambda = center.mu ./ E;

    path = struct( ...
        'deltaPphiPlot', deltaPphiPlot, ...
        'deltaPphiRaw', deltaPphiRaw, ...
        'E', E, ...
        'PphiPlot', PphiPlot, ...
        'Pphi', Pphi, ...
        'Lambda', Lambda, ...
        'rho', NaN(size(E)), ...
        'r', NaN(size(E)), ...
        'x', NaN(size(E)), ...
        'R', NaN(size(E)), ...
        'valid', false(size(E)), ...
        'aborted', false, ...
        'abortReason', '');

    outside = ~isfinite(E) | ~isfinite(Pphi) | ~isfinite(Lambda) | ...
        E < min(speciesData.E1d) | E > max(speciesData.E1d) | ...
        Pphi < min(speciesData.Pphi1d) | Pphi > max(speciesData.Pphi1d) | ...
        Lambda < min(speciesData.Lambda1d) | Lambda > max(speciesData.Lambda1d);
    if any(outside)
        path.aborted = true;
        path.abortReason = sprintf('%d/%d sampled points are outside the phase-space grid.', nnz(outside), numel(outside));
        warning(['[resonance detuning] path exceeds phase-space grid; stop this diagnostic before interpolation/fit/plot. ', ...
            'outside=%d/%d, E=[%.6g, %.6g] grid=[%.6g, %.6g], Pphi(raw)=[%.6g, %.6g] grid=[%.6g, %.6g], ', ...
            'Lambda=[%.6g, %.6g] grid=[%.6g, %.6g].'], ...
            nnz(outside), numel(outside), ...
            statisticOrNaN(E, @min, true), statisticOrNaN(E, @max, true), min(speciesData.E1d), max(speciesData.E1d), ...
            statisticOrNaN(Pphi, @min, true), statisticOrNaN(Pphi, @max, true), min(speciesData.Pphi1d), max(speciesData.Pphi1d), ...
            statisticOrNaN(Lambda, @min, true), statisticOrNaN(Lambda, @max, true), min(speciesData.Lambda1d), max(speciesData.Lambda1d));
        return;
    end

    residualInterpolant = phaseFieldInterpolant(speciesData, residualHz);
    radiusInterpolant = phaseFieldInterpolant(speciesData, mapping.r);
    rhoInterpolant = phaseFieldInterpolant(speciesData, mapping.rho);

    R = residualInterpolant(E, Pphi, Lambda);
    r = radiusInterpolant(E, Pphi, Lambda);
    rho = rhoInterpolant(E, Pphi, Lambda);
    x = r - center.r;
    centerIndex = find(deltaPphiPlot == 0, 1);
    R(centerIndex) = 0;
    r(centerIndex) = center.r;
    rho(centerIndex) = center.rho;
    x(centerIndex) = 0;
    valid = isfinite(E) & isfinite(Pphi) & isfinite(Lambda) & isfinite(R) & isfinite(r) & isfinite(x);

    path.rho = rho;
    path.r = r;
    path.x = x;
    path.R = R;
    path.valid = valid;

    fprintf(['[resonance detuning] sampled path: sampleCount=%d, valid=%d/%d, ', ...
        'Pphi(plot)=[%.6g, %.6g], Pphi(raw)=[%.6g, %.6g], x=[%.6g, %.6g].\n'], ...
        numel(deltaPphiPlot), nnz(valid), numel(valid), min(PphiPlot), max(PphiPlot), ...
        min(Pphi), max(Pphi), statisticOrNaN(x(valid), @min, true), statisticOrNaN(x(valid), @max, true));
end

function fit = fitResonanceDetuningPath(path)

    valid = path.valid(:) & isfinite(path.x(:)) & isfinite(path.R(:));
    x = path.x(valid);
    R = path.R(valid);

    [x, order] = sort(x);
    R = R(order);
    [~, centerIndex] = min(abs(x));
    linearIndex = centerIndex - 2:centerIndex + 2;

    xLinear = x(linearIndex);
    RLinear = R(linearIndex);
    fit1 = noInterceptPolynomialFit(xLinear, RLinear, 1);
    fit2 = noInterceptPolynomialFit(x, R, 2);

    a2 = fit2.coeff(1);
    b2 = fit2.coeff(2);
    ratio2 = b2 * x / max(abs(a2), realmin) * sign(a2);

    fit = struct( ...
        'x', x, ...
        'R', R, ...
        'linearX', xLinear, ...
        'linearR', RLinear, ...
        'fit1', fit1, ...
        'fit2', fit2, ...
        'linearYFit', fit1.coeff(1) * x, ...
        'quadraticRatio', ratio2, ...
        'maxAbsQuadraticRatio', max(abs(ratio2)));
end

function fit = noInterceptPolynomialFit(x, y, order)

    X = zeros(numel(x), order);
    for termIndex = 1:order
        X(:, termIndex) = x .^ termIndex;
    end
    coeff = X \ y;
    yFit = X * coeff;
    rmse = sqrt(mean((y - yFit) .^ 2));
    yRange = max(y) - min(y);
    fit = struct('order', order, 'coeff', coeff(:).', 'yFit', yFit, ...
        'rmse', rmse, 'nrmse', rmse / max(abs(yRange), realmin));
end

function printResonanceDetuningFitSummary(fit)

    c1 = fit.fit1.coeff;
    c2 = fit.fit2.coeff;
    fprintf(['[resonance detuning] fit summary:\n', ...
        '  R=a*x, 5 points near R=0: a=%.6g, nrmse=%.6g, x=[%.6g, %.6g]\n', ...
        '  R=a*x+b*x^2, all points:  a=%.6g, b=%.6g, nrmse=%.6g, max|b*x/a|=%.6g\n'], ...
        c1(1), fit.fit1.nrmse, min(fit.linearX), max(fit.linearX), ...
        c2(1), c2(2), fit.fit2.nrmse, fit.maxAbsQuadraticRatio);
end

function plotResonanceDetuningPath(path, fit, center, opt)

    figure('Name', 'resonance detuning R(x)', 'Color', 'w', 'Position', [100, 80, 980, 760]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');

    nexttile;
    plot(path.x(path.valid), path.R(path.valid), 'ko', 'MarkerSize', 4, 'DisplayName', 'R samples'); hold on;
    plot(fit.x, fit.linearYFit, 'b--', 'LineWidth', 1.2, 'DisplayName', 'a x');
    plot(fit.x, fit.fit2.yFit, 'm-', 'LineWidth', 1.2, 'DisplayName', 'a x + b x^2');
    xline(0, 'k:');
    yline(0, 'k:');
    grid on;
    xlabel('x = r - r_0');
    ylabel('R / Hz');
    title(sprintf('R(x), E0=%.4g, Pphi0(plot)=%.4g, Lambda0=%.4g, f=%.4g Hz, n=%g', ...
        center.E, -center.Pphi, center.Lambda, opt.frequencyHz, opt.toroidalMode), 'Interpreter', 'none');
    legend('Location', 'best');

    nexttile;
    plot(fit.x, fit.quadraticRatio, 'm-', 'LineWidth', 1.2, 'DisplayName', 'b x / a'); hold on;
    xline(0, 'k:');
    yline(0, 'k:');
    grid on;
    xlabel('x = r - r_0');
    ylabel('relative detuning contribution');
    title(sprintf('max |b x/a|=%.3g', fit.maxAbsQuadraticRatio));
    legend('Location', 'best');
end

function [Z, xVec, yVec, xlabelText, ylabelText, titleText] = slicePhaseField(dim, sliceIndex, fieldData, speciesData)

    switch dim
        case 1
            sliceIndex = clampIndex(sliceIndex, numel(speciesData.E1d), 'E');
            Z = reshape(fieldData(sliceIndex, :, :), numel(speciesData.Pphi1d), numel(speciesData.Lambda1d));
            xVec = speciesData.Lambda1d;
            yVec = -speciesData.Pphi1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$P_{\varphi}$';
            titleText = sprintf('$E = %.6g$', speciesData.E1d(sliceIndex));
        case 2
            sliceIndex = clampIndex(sliceIndex, numel(speciesData.Pphi1d), 'Pphi');
            Z = reshape(fieldData(:, sliceIndex, :), numel(speciesData.E1d), numel(speciesData.Lambda1d));
            xVec = speciesData.Lambda1d;
            yVec = speciesData.E1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$E$';
            titleText = sprintf('$P_{\\varphi} = %.6g$', -speciesData.Pphi1d(sliceIndex));
        case 3
            sliceIndex = clampIndex(sliceIndex, numel(speciesData.Lambda1d), 'Lambda');
            Z = reshape(fieldData(:, :, sliceIndex), numel(speciesData.E1d), numel(speciesData.Pphi1d));
            xVec = -speciesData.Pphi1d;
            yVec = speciesData.E1d;
            xlabelText = '$P_{\varphi}$';
            ylabelText = '$E$';
            titleText = sprintf('$\\Lambda = %.6g$', speciesData.Lambda1d(sliceIndex));
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
    physicalN = speciesData.physicalNAll;
    nTarget = diffusivitySingleN(speciesData, nRange);
    [~, modeIndex] = min(abs(physicalN - nTarget));
    selectedN = physicalN(modeIndex);
    Dn = reshape(data(:, modeIndex, :), size(data, 1), size(data, 3));
    Dn(~isfinite(Dn)) = NaN;
end

function [modeMask, selectedN] = diffusivityModeMask(speciesData, nRange)

    physicalN = speciesData.physicalNAll;
    range = normalizeDiffusivityNRange(nRange);
    modeMask = physicalN >= range(1) & physicalN <= range(2);
    if ~any(modeMask)
        [~, nearestIndex] = min(abs(physicalN - mean(range)));
        modeMask(nearestIndex) = true;
    end
    selectedN = physicalN(modeMask);
end

function range = normalizeDiffusivityNRange(nRange)

    range = reshape(double(nRange), 1, []);
    assert(ismember(numel(range), [1, 2]) && all(isfinite(range)), ...
        'nRange 必须为有限标量或 [min max]。');
    if isscalar(range)
        range = [range, range];
    end
    range = sort(range);
end

function range = initialDiffusivityNRange(speciesData, nRange)

    [~, selectedN] = diffusivityModeMask(speciesData, nRange);
    range = [min(selectedN), max(selectedN)];
end

function nValue = diffusivitySingleN(speciesData, nRange)

    physicalN = speciesData.physicalNAll;
    rawRange = normalizeDiffusivityNRange(nRange);
    inRangeN = physicalN(physicalN >= rawRange(1) & physicalN <= rawRange(2));
    if isempty(inRangeN)
        [~, nearestIndex] = min(abs(physicalN - mean(rawRange)));
        nValue = physicalN(nearestIndex);
    else
        nValue = inRangeN(max(1, round(numel(inRangeN) / 2)));
    end
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

function controls = appendIndexControl(controls, fieldName, value, count, enabled)

    if nargin < 5
        enabled = true;
    end
    if enabled && count > 1
        controls = [controls, integerSliderControl(fieldName, fieldName, value, 1, count)];
    end
end

function opt = diffusivityOptionsFromValues(opt, values)

    opt = applyControlValues(opt, values);
    if isfield(values, 'nMin') && isfield(values, 'nMax')
        opt.nRange = sort([values.nMin, values.nMax]);
    end
end

function [xVec, xlabelText] = diagnosticTimeAxis(speciesData, axisType)

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
            xVec = (0:numel(speciesData.tDiag) - 1) * speciesData.diagSteps;
            xlabelText = '$\mathrm{step}$';
        otherwise
            error('timeAxis 必须为 "ta"、"ms"、"s" 或 "steps"。');
    end
end

function [xVec, xlabelText, axisName] = diagnosticRadialAxis(speciesData, axisType)

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

function titleText = diffusivityTitleText(speciesLabel, selectedN, extraLatex)

    if isscalar(selectedN)
        nLatex = sprintf('n = %d', selectedN);
    else
        nLatex = sprintf('n \\in [%d,%d]', min(selectedN), max(selectedN));
    end
    titleText = sprintf('$\\mathrm{%s}\\quad D,\\quad %s%s$', speciesLabel, nLatex, extraLatex);
end

function titleText = densityTitleText(speciesLabel, extraLatex)

    titleText = sprintf('$\\mathrm{%s}\\quad \\delta n%s$', speciesLabel, extraLatex);
end

function flag = isTotalDensityMode(opt)

    modeText = lower(strtrim(char(getOptionValue(opt, 'radialDensityMode', 'delta'))));
    assert(ismember(modeText, {'delta', 'total'}), ...
        'radialDensityMode 必须为 "delta" 或 "total"。');
    flag = strcmp(modeText, 'total');
end

function values = applyControlValues(values, controlValues)

    fieldNames = fieldnames(controlValues);
    for fieldIndex = 1:numel(fieldNames)
        fieldName = fieldNames{fieldIndex};
        values.(fieldName) = controlValues.(fieldName);
    end
end

function runPICPlotMode(figureName, controls, dynamicUpdate, buildPlotData, ...
    renderPlotData, noFiniteMessage)

    if ~isempty(dynamicUpdate)
        plotInteractiveDiagnostic(figureName, controls, dynamicUpdate, buildPlotData, renderPlotData);
        return;
    end

    plotData = buildPlotData(struct());
    if ~any(isfinite(plotData.finiteData(:)))
        fprintf('[plot] %s\n', noFiniteMessage);
        return;
    end
    drawPICFigure(figureName, plotData, renderPlotData);
end

function plotData = buildLinePlotData(xVec, yVec, xlabelText, ylabelText, titleText, statusText)

    xVec = reshape(xVec, 1, []);
    yVec = reshape(yVec, 1, []);
    assert(numel(xVec) == numel(yVec), '线图横纵坐标长度不一致。');
    plotData = struct( ...
        'xVec', xVec, ...
        'yVec', yVec, ...
        'xlabelText', xlabelText, ...
        'ylabelText', ylabelText, ...
        'titleText', titleText, ...
        'status', statusText, ...
        'finiteData', yVec);
end

function drawPICFigure(figureName, plotData, renderPlotData)

    layoutIndex = 1 + isfield(plotData, 'Z');
    figurePositions = {[120, 120, 900, 560], [100, 80, 980, 760]};
    axesPositions = {[0.12, 0.14, 0.82, 0.76], [0.11, 0.12, 0.74, 0.78]};
    figHandle = figure('Name', figureName, 'Color', 'w', ...
        'Position', figurePositions{layoutIndex});
    axHandle = axes('Parent', figHandle, 'Units', 'normalized', ...
        'Position', axesPositions{layoutIndex});
    renderPlotData(axHandle, plotData);
end

function renderLine(axHandle, plotData)

    resetPlotAxes(axHandle);
    hold(axHandle, 'on');
    maxLegendCount = 1;
    if isfield(plotData, 'lineOverlays') && ~isempty(plotData.lineOverlays)
        maxLegendCount = maxLegendCount + numel(plotData.lineOverlays);
    end
    legendHandles = gobjects(1, maxLegendCount);
    legendLabels = cell(1, maxLegendCount);
    legendCount = 0;
    hasFiniteLine = false;

    validData = isfinite(plotData.xVec) & isfinite(plotData.yVec);
    if any(validData)
        lineHandle = plot(axHandle, plotData.xVec, plotData.yVec, ...
            'LineWidth', 1.8, 'Color', [0.20, 0.40, 0.80]);
        hasFiniteLine = true;
        if isfield(plotData, 'lineLabel') && ~isempty(plotData.lineLabel)
            legendCount = legendCount + 1;
            legendHandles(legendCount) = lineHandle;
            legendLabels{legendCount} = plotData.lineLabel;
        end
    end

    if isfield(plotData, 'lineOverlays') && ~isempty(plotData.lineOverlays)
        overlays = plotData.lineOverlays;
        colorTable = lines(max(numel(overlays) + 1, 7));
        for overlayIndex = 1:numel(overlays)
            overlayX = reshape(overlays(overlayIndex).xVec, 1, []);
            overlayY = reshape(overlays(overlayIndex).yVec, 1, []);
            assert(numel(overlayX) == numel(overlayY), ...
                '线图叠加曲线横纵坐标长度不一致。');
            overlayValid = isfinite(overlayX) & isfinite(overlayY);
            if any(overlayValid)
                lineHandle = plot(axHandle, overlayX, overlayY, ...
                    'LineWidth', 1.8, 'Color', colorTable(overlayIndex + 1, :));
                hasFiniteLine = true;
                if isfield(overlays, 'label') && ~isempty(overlays(overlayIndex).label)
                    legendCount = legendCount + 1;
                    legendHandles(legendCount) = lineHandle;
                    legendLabels{legendCount} = overlays(overlayIndex).label;
                end
            end
        end
    end

    if ~hasFiniteLine
        text(axHandle, 0.5, 0.5, 'no finite data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 14);
    end

    if legendCount > 0
        legend(axHandle, legendHandles(1:legendCount), legendLabels(1:legendCount), 'Interpreter', 'latex', ...
            'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
    end

    xlabel(axHandle, plotData.xlabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(axHandle, plotData.ylabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(axHandle, plotData.titleText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 13);
    grid(axHandle, 'on');
    axis(axHandle, 'tight');
    box(axHandle, 'on');
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14, ...
        'Color', 'white', 'Layer', 'top', 'TickDir', 'out', 'LineWidth', 1);
    hold(axHandle, 'off');
end

function renderMap(axHandle, plotData)

    resetPlotAxes(axHandle);
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
    [colorMin, colorMax] = finiteColorLimits(plotData.Z, isSigned);
    clim(axHandle, [colorMin, colorMax]);
    colormap(axHandle, phaseSpaceColormap(isSigned, plotData.colormapIndex, 256));

    colorbarHandle = colorbar(axHandle);
    colorbarHandle.FontName = 'Times New Roman';
    colorbarHandle.FontSize = 13;
    ylabel(colorbarHandle, plotData.colorbarLabel, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 13);

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
    xLimits = xlim(axHandle);
    yLimits = ylim(axHandle);
    if isfield(plotData, 'detuningLine') && ~isempty(plotData.detuningLine)
        plot(axHandle, plotData.detuningLine.x, plotData.detuningLine.y, ...
            'k-', 'LineWidth', 1.8, 'HandleVisibility', 'off');
        xlim(axHandle, xLimits);
        ylim(axHandle, yLimits);
    end
    box(axHandle, 'on');
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14, ...
        'Color', 'white', 'Layer', 'top', 'TickDir', 'out', 'LineWidth', 1);
    hold(axHandle, 'off');
end

function drawResonanceOverlays(axHandle, X, Y, overlays)

    legendHandles = gobjects(0);
    legendLabels = {};
    visibleIndex = 0;
    for overlayIndex = 1:numel(overlays)
        resonanceZ = overlays(overlayIndex).Z;
        if ~residualHasZeroContour(resonanceZ)
            continue;
        end

        visibleIndex = visibleIndex + 1;
        lineColor = resonanceOverlayColor(visibleIndex, numel(overlays));
        contour(axHandle, X, Y, resonanceZ, [0, 0], ...
            'EdgeColor', 'w', 'LineWidth', 2.8, 'HandleVisibility', 'off');
        [~, lineHandle] = contour(axHandle, X, Y, resonanceZ, [0, 0], ...
            'EdgeColor', lineColor, 'LineWidth', 1.35);
        legendHandles(end + 1) = lineHandle;
        legendLabels{end + 1} = overlays(overlayIndex).label;
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

function resetPlotAxes(axHandle)

    figHandle = ancestor(axHandle, 'figure');
    delete(findall(figHandle, 'Type', 'ColorBar'));
    delete(findall(figHandle, 'Type', 'Legend'));
    cla(axHandle);
end

function plotData = buildMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
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
        'finiteData', Z, ...
        'forceSigned', false, ...
        'resonanceZ', [], ...
        'resonanceOverlays', struct('Z', {}, 'label', {}, 'harmonic', {}), ...
        'detuningLine', []);
end

function plotInteractiveDiagnostic(figureName, controls, dynamicUpdate, computePlotData, renderPlotData)

    figHandle = figure('Name', figureName, 'Color', 'w', 'Position', [80, 40, 1080, 860]);
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
    for controlIndex = 1:nControl
        yPos = controlBottom + controlSpacing * (nControl - controlIndex);
        sliderLabels(controlIndex) = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.10, yPos - 0.007, 0.19, 0.030], ...
            'String', sliderLabelText(controls(controlIndex), controls(controlIndex).value), ...
            'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', 11);
        sliders(controlIndex) = uicontrol(figHandle, 'Style', 'slider', 'Units', 'normalized', ...
            'Position', [0.30, yPos, 0.62, 0.022], ...
            'Min', controls(controlIndex).min, 'Max', controls(controlIndex).max, ...
            'Value', controls(controlIndex).value, ...
            'SliderStep', sliderStepForControl(controls(controlIndex)), ...
            'Callback', @refreshPlot);
    end

    if dynamicUpdate
        dynamicListeners = {};
        try
            for controlIndex = 1:nControl
                dynamicListeners{end + 1} = addlistener(sliders(controlIndex), 'ContinuousValueChange', @refreshPlot); %#ok<AGROW>
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
        for sliderIndex = 1:nControl
            sliderValue = clampSliderValue(get(sliders(sliderIndex), 'Value'), controls(sliderIndex));
            values.(controls(sliderIndex).field) = sliderValue;
            if ~dynamicUpdate
                set(sliders(sliderIndex), 'Value', sliderValue);
            end
            set(sliderLabels(sliderIndex), 'String', sliderLabelText(controls(sliderIndex), sliderValue));
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
        resOpt = normalizeResonanceOverlayOptions(resOpt);
        branchNames = normalizeResonanceBranches(resOpt.branch);
    else
        resOpt = normalizeResonanceOptions(resOpt);
        branchNames = {resOpt.branch};
    end
    freqRange = normalizeSliderRange(resOpt.frequencyHzRange, resOpt.frequencyHz, 'frequencyHzRange', false);
    nRange = normalizeSliderRange(resOpt.toroidalModeRange, resOpt.toroidalMode, 'toroidalModeRange', true);
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
        if allowMultipleBranches && numel(branchNames) > 1
            for branchIndex = 1:numel(branchNames)
                lBounds = [resOpt.harmonicMinList(branchIndex), resOpt.harmonicMaxList(branchIndex)];
                lRange = normalizeSliderRange(resOpt.harmonicRange, lBounds, 'harmonicRange', true);
                controls = [controls, ...
                    integerSliderControl(sprintf('resHarmonicMin%d', branchIndex), ...
                    [branchNames{branchIndex} ' l min'], lBounds(1), lRange(1), lRange(2)), ...
                    integerSliderControl(sprintf('resHarmonicMax%d', branchIndex), ...
                    [branchNames{branchIndex} ' l max'], lBounds(2), lRange(1), lRange(2))]; %#ok<AGROW>
            end
        else
            lBounds = [resOpt.harmonicMin, resOpt.harmonicMax];
            lRange = normalizeSliderRange(resOpt.harmonicRange, lBounds, 'harmonicRange', true);
            controls = [controls, ...
                integerSliderControl('resHarmonicMin', 'l min', lBounds(1), lRange(1), lRange(2)), ...
                integerSliderControl('resHarmonicMax', 'l max', lBounds(2), lRange(1), lRange(2))];
        end
    else
        lBounds = [resOpt.harmonicMin, resOpt.harmonicMax];
        lRange = normalizeSliderRange(resOpt.harmonicRange, lBounds, 'harmonicRange', true);
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
    valueMap = { ...
        'resFrequencyHz', 'frequencyHz'; 'resToroidalMode', 'toroidalMode'; ...
        'resPoloidalMode', 'poloidalMode'; 'resHarmonicMin', 'harmonicMin'; ...
        'resHarmonicMax', 'harmonicMax'};
    for fieldIndex = 1:size(valueMap, 1)
        if isfield(values, valueMap{fieldIndex, 1})
            resOpt.(valueMap{fieldIndex, 2}) = values.(valueMap{fieldIndex, 1});
        end
    end
    if isfield(values, 'resHarmonic')
        resOpt.harmonic = values.resHarmonic;
        resOpt.harmonicMin = values.resHarmonic;
        resOpt.harmonicMax = values.resHarmonic;
    end
    if allowMultipleBranches
        branchNames = normalizeResonanceBranches(resOpt.branch);
        if numel(branchNames) > 1
            resOpt = normalizeResonanceHarmonicOptions(resOpt);
            resOpt = normalizeResonanceOverlayHarmonicOptions(resOpt, numel(branchNames));
            for branchIndex = 1:numel(branchNames)
                minField = sprintf('resHarmonicMin%d', branchIndex);
                maxField = sprintf('resHarmonicMax%d', branchIndex);
                if isfield(values, minField)
                    resOpt.harmonicMinList(branchIndex) = values.(minField);
                end
                if isfield(values, maxField)
                    resOpt.harmonicMaxList(branchIndex) = values.(maxField);
                end
            end
        end
        resOpt = normalizeResonanceOverlayOptions(resOpt);
    else
        resOpt = normalizeResonanceOptions(resOpt);
    end
end

function resOpt = normalizeResonanceOptions(resOpt)

    resOpt.branch = normalizeResonanceBranch(resOpt.branch);
    resOpt = normalizeResonanceCoreOptions(resOpt);
    if ~strcmp(resOpt.branch, 'trapped')
        resOpt.poloidalMode = requireIntegerInRange(resOpt.poloidalMode, 0, Inf, 'poloidalMode');
    end
    resOpt = normalizeResonanceHarmonicOptions(resOpt);
end

function resOpt = normalizeResonanceOverlayOptions(resOpt)

    branchNames = normalizeResonanceBranches(resOpt.branch);
    resOpt.branch = branchNames;
    resOpt = normalizeResonanceCoreOptions(resOpt);
    if any(~strcmp(branchNames, 'trapped'))
        resOpt.poloidalMode = requireIntegerInRange(resOpt.poloidalMode, 0, Inf, 'poloidalMode');
    elseif ~isfield(resOpt, 'poloidalMode') || isempty(resOpt.poloidalMode)
        resOpt.poloidalMode = 0;
    end
    resOpt = normalizeResonanceHarmonicOptions(resOpt);
    resOpt = normalizeResonanceOverlayHarmonicOptions(resOpt, numel(branchNames));
end


function resOpt = normalizeResonanceCoreOptions(resOpt)

    if ~isfield(resOpt, 'enabled')
        resOpt.enabled = true;
    end
    resOpt.frequencyHz = requireFiniteScalarInRange(resOpt.frequencyHz, -Inf, Inf, 'frequencyHz');
    resOpt.toroidalMode = requireIntegerInRange(resOpt.toroidalMode, 0, Inf, 'toroidalMode');
end

function resOpt = normalizeResonanceHarmonicOptions(resOpt)

    if ~isfield(resOpt, 'harmonic') || isempty(resOpt.harmonic)
        resOpt.harmonic = 0;
    end
    rawHarmonic = reshape(double(resOpt.harmonic), 1, []);
    assert(~isempty(rawHarmonic) && all(isfinite(rawHarmonic)) && all(rawHarmonic == floor(rawHarmonic)), ...
        'harmonic 必须是整数标量或整数范围。');

    [resOpt.harmonicMin, resOpt.harmonicMax] = normalizeIntegerBounds( ...
        getOptionValue(resOpt, 'harmonicMin', min(rawHarmonic)), ...
        getOptionValue(resOpt, 'harmonicMax', max(rawHarmonic)), 1, 'harmonic bounds');
    resOpt.harmonic = requireIntegerInRange(rawHarmonic(1), -Inf, Inf, 'harmonic');
end

function resOpt = normalizeResonanceOverlayHarmonicOptions(resOpt, nBranch)

    [resOpt.harmonicMinList, resOpt.harmonicMaxList] = normalizeIntegerBounds( ...
        getOptionValue(resOpt, 'harmonicMinList', resOpt.harmonicMin), ...
        getOptionValue(resOpt, 'harmonicMaxList', resOpt.harmonicMax), ...
        nBranch, 'harmonicMinList/harmonicMaxList');
end

function [lowerValues, upperValues] = normalizeIntegerBounds(lowerValues, upperValues, count, fieldName)

    lowerValues = reshape(double(lowerValues), 1, []);
    upperValues = reshape(double(upperValues), 1, []);
    if isscalar(lowerValues)
        lowerValues = repmat(lowerValues, 1, count);
    end
    if isscalar(upperValues)
        upperValues = repmat(upperValues, 1, count);
    end
    assert(numel(lowerValues) == count && numel(upperValues) == count && ...
        all(isfinite(lowerValues)) && all(isfinite(upperValues)) && ...
        all(lowerValues == floor(lowerValues)) && all(upperValues == floor(upperValues)), ...
        '%s 必须包含 %d 组有限整数边界。', fieldName, count);
    bounds = sort([lowerValues; upperValues], 1);
    lowerValues = bounds(1, :);
    upperValues = bounds(2, :);
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

    if isfield(resOpt, 'harmonicMinList') && isfield(resOpt, 'harmonicMaxList') && ...
            numel(resOpt.harmonicMinList) > 1
        text = 'l = branch-specific';
        return;
    end

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

function ok = hasSpeciesFields(dataStruct, speciesName, fieldNames)

    speciesName = strtrim(char(speciesName));
    ok = isfield(dataStruct, speciesName) && ~isempty(dataStruct.(speciesName));
    for fieldIndex = 1:numel(fieldNames)
        fieldName = fieldNames{fieldIndex};
        ok = ok && isfield(dataStruct.(speciesName), fieldName) && ...
            ~isempty(dataStruct.(speciesName).(fieldName));
    end
end

function ok = hasRequestedPICQuantity(dataStruct, speciesName, quantityName)

    ok = hasSpeciesFields(dataStruct, speciesName, {});
    if ~ok, return; end
    quantitySpec = resolvePICQuantitySpec(quantityName);
    ok = hasSpeciesFields(dataStruct, speciesName, quantitySpec.requiredFields);
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

function sliceIndex = parsePhaseSlice(sliceText, dim, speciesData)

    if isnumeric(sliceText)
        sliceIndex = double(sliceText);
    else
        key = lower(strtrim(char(sliceText)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            sliceIndex = defaultSliceIndex(dim, speciesData.gridE, speciesData.gridPphi, speciesData.gridLambda);
        elseif ismember(key, {'end', 'last'})
            sliceIndex = phaseDimensionSize(dim, speciesData);
        else
            sliceIndex = str2double(key);
        end
    end

    coordinateLabels = {'E', 'Pphi', 'Lambda'};
    sliceIndex = clampIndex(sliceIndex, phaseDimensionSize(dim, speciesData), coordinateLabels{dim});
end

function index = parseIndex(indexText, nIndex, label)

    if isnumeric(indexText)
        index = double(indexText);
    else
        key = lower(strtrim(char(indexText)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            index = max(1, round(nIndex / 2));
        elseif ismember(key, {'end', 'last'})
            index = nIndex;
        else
            index = str2double(key);
        end
    end

    index = clampIndex(index, nIndex, label);
end

function timeIndex = parseTimeIndex(indexText, nTime)

    timeIndex = parseIndex(indexText, nTime, 'timeIndex');
end

function [modeIndex, modeN] = parseModeN(modeNText, modeIndexAll, physicalNAll)

    assert(~isempty(modeIndexAll) && numel(modeIndexAll) == numel(physicalNAll), ...
        'modeN 列表为空或尺寸不一致。');
    if isnumeric(modeNText)
        modeN = double(modeNText);
    else
        key = lower(strtrim(char(modeNText)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            modeIndex = max(1, round(numel(modeIndexAll) / 2));
            modeN = physicalNAll(modeIndex);
            return;
        elseif ismember(key, {'end', 'last'})
            modeIndex = numel(modeIndexAll);
            modeN = physicalNAll(modeIndex);
            return;
        elseif ismember(key, {'first', 'begin'})
            modeIndex = 1;
            modeN = physicalNAll(modeIndex);
            return;
        else
            modeN = str2double(key);
        end
    end

    assert(isscalar(modeN) && isfinite(modeN) && modeN == floor(modeN) && modeN >= 0, ...
        'modeN 必须为非负整数物理环向模数。');
    modeIndex = find(physicalNAll == modeN, 1);
    assert(~isempty(modeIndex), 'modeN 必须位于有效物理 n 集合 [%s]。当前值：%d。', ...
        strjoin(cellstr(compose('%g', physicalNAll(:))), ', '), modeN);
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

function sliceIndex = defaultSliceIndex(dim, gridE, gridPphi, gridLambda)

    switch dim
        case 1
            sliceIndex = max(1, round(gridE / 2));
        case 2
            sliceIndex = max(1, round(gridPphi / 2));
        case 3
            sliceIndex = max(1, round(gridLambda / 2));
        otherwise
            error('切片维度必须为 1、2 或 3。');
    end
end

function index = clampIndex(index, maxIndex, label)

    assert(isscalar(index) && isfinite(index) && index == floor(index), ...
        '%s 下标必须是有限整数。', label);
    assert(index >= 1 && index <= maxIndex, ...
        '%s 下标 %d 超出有效范围 [1, %d]。', label, index, maxIndex);
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
    for nameIndex = 1:numel(rawNames)
        nameText = lower(strtrim(char(rawNames{nameIndex})));
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
    for directoryIndex = 1:numel(searchDirs)
        candidate = fullfile(searchDirs{directoryIndex}, fileName);
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

    dimensionText = sprintf('gridE=%d, gridPphi=%d, gridLambda=%d', gridE, gridPphi, gridLambda);
    data = readPICBinaryArray(filePath, precision, ...
        [gridLambda, gridPphi, gridE], [3, 2, 1], dimensionText);
end

function data = readPhase4D(filePath, precision, gridE, gridPphi, gridLambda, expectedTime)

    dimensionText = sprintf('gridE=%d, gridPphi=%d, gridLambda=%d, expectedTime=%d', ...
        gridE, gridPphi, gridLambda, expectedTime);
    data = readPICBinaryArray(filePath, precision, ...
        [gridLambda, gridPphi, gridE, expectedTime], [3, 2, 1, 4], dimensionText);
end

function data = readPhasePower5D(filePath, precision, gridE, gridPphi, gridLambda, modeCount, expectedTime)

    dimensionText = sprintf( ...
        'gridE=%d, gridPphi=%d, gridLambda=%d, modeCount=%d, expectedTime=%d', ...
        gridE, gridPphi, gridLambda, modeCount, expectedTime);
    data = readPICBinaryArray(filePath, precision, ...
        [gridLambda, gridPphi, gridE, modeCount, expectedTime], [3, 2, 1, 4, 5], dimensionText);
end

function data = readPitch2D(filePath, precision, gridVpara, gridVperp)

    dimensionText = sprintf('gridVpara=%d, gridVperp=%d', gridVpara, gridVperp);
    data = readPICBinaryArray(filePath, precision, ...
        [gridVperp, gridVpara], [2, 1], dimensionText);
end

function data = readPitch3D(filePath, precision, gridVpara, gridVperp, expectedTime)

    dimensionText = sprintf('gridVpara=%d, gridVperp=%d, expectedTime=%d', ...
        gridVpara, gridVperp, expectedTime);
    data = readPICBinaryArray(filePath, precision, ...
        [gridVperp, gridVpara, expectedTime], [2, 1, 3], dimensionText);
end

function data = readPitchPower4D(filePath, precision, gridVpara, gridVperp, modeCount, expectedTime)

    dimensionText = sprintf( ...
        'gridVpara=%d, gridVperp=%d, modeCount=%d, expectedTime=%d', ...
        gridVpara, gridVperp, modeCount, expectedTime);
    data = readPICBinaryArray(filePath, precision, ...
        [gridVperp, gridVpara, modeCount, expectedTime], [2, 1, 3, 4], dimensionText);
end

function data = readDiffusivity3D(filePath, precision, expectedTime, modeCount, gridNx)

    dimensionText = sprintf('expectedTime=%d, modeCount=%d, gridNx=%d', ...
        expectedTime, modeCount, gridNx);
    data = readPICBinaryArray(filePath, precision, ...
        [gridNx, modeCount, expectedTime], [3, 2, 1], dimensionText);
end

function data = readDensity2D(filePath, precision, expectedTime, gridNx)

    dimensionText = sprintf('expectedTime=%d, gridNx=%d', expectedTime, gridNx);
    data = readPICBinaryArray(filePath, precision, ...
        [gridNx, expectedTime], [2, 1], dimensionText);
end

function data = readPICBinaryArray(filePath, precision, storageShape, permutation, dimensionText)

    raw = readBinaryVector(filePath, precision);
    expectedCount = prod(storageShape);
    assert(numel(raw) == expectedCount, ...
        '%s 尺寸不匹配：读到 %d 个数，期望 %d 个（%s）。', ...
        filePath, numel(raw), expectedCount, dimensionText);
    data = permute(reshape(raw, storageShape), permutation);
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

    gridE = requireIntegerInRange(normData.gridE, 1, Inf, 'gridE');
    gridPphi = requireIntegerInRange(normData.gridPphi, 1, Inf, 'gridPphi');
    gridLambda = requireIntegerInRange(normData.gridLambda, 1, Inf, 'gridLambda');
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

function value = requireFiniteScalarInRange(value, minValue, maxValue, fieldName)

    value = double(value);
    assert(isscalar(value) && isfinite(value), '%s 必须是有限标量。', fieldName);
    assert((~isfinite(minValue) || value >= minValue) && ...
        (~isfinite(maxValue) || value <= maxValue), ...
        '%s 必须位于 [%g, %g]。', fieldName, minValue, maxValue);
end

function value = requireIntegerInRange(value, minValue, maxValue, fieldName)

    value = requireFiniteScalarInRange(value, minValue, maxValue, fieldName);
    assert(value == floor(value), '%s 必须是整数。', fieldName);
end

function range = normalizeSliderRange(rawRange, initialValues, fieldName, isInteger)

    range = reshape(double(rawRange), 1, []);
    initialValues = reshape(double(initialValues), 1, []);
    assert(numel(range) == 2 && all(isfinite(range)), '%s 必须为 [min max]。', fieldName);
    assert(~isempty(initialValues) && all(isfinite(initialValues)), '%s 初始值必须有限。', fieldName);
    range = sort(range);
    range = [min([range(1), initialValues]), max([range(2), initialValues])];

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

    [eIndex, pphiIndex, lambdaIndex] = localIdsToSubscripts(localIds, gridPphi, gridLambda);
    data = sub2ind([gridE, gridPphi, gridLambda], eIndex, pphiIndex, lambdaIndex);
end

function [eIndex, pphiIndex, lambdaIndex] = localIdsToSubscripts(localIds, gridPphi, gridLambda)

    localIds = double(localIds(:));
    strideE = gridPphi * gridLambda;
    eIndex = floor(localIds / strideE) + 1;
    remainder = mod(localIds, strideE);
    pphiIndex = floor(remainder / gridLambda) + 1;
    lambdaIndex = mod(remainder, gridLambda) + 1;
end

function [baseE, basePphi, baseLambda] = initialCoordinatesFromLocalId(localIds, E1d, Pphi1d, Lambda1d)

    [eIndex, pphiIndex, lambdaIndex] = localIdsToSubscripts(localIds, numel(Pphi1d), numel(Lambda1d));
    baseE = E1d(eIndex(:));
    basePphi = Pphi1d(pphiIndex(:));
    baseLambda = Lambda1d(lambdaIndex(:));
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
    for plotIndex = 1:numel(invariantLabels)
        axHandle = subplot(3, 1, plotIndex);
        plot(axHandle, invariantErrors{plotIndex}, '.');
        grid(axHandle, 'on');
        ylabel(axHandle, invariantLabels{plotIndex}, 'FontName', 'Times New Roman', 'FontSize', 14);
        set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14);
        if plotIndex == 1
            title(axHandle, sprintf('%s %s relative error', speciesName, classLabel), ...
                'Interpreter', 'none', 'FontName', 'Times New Roman', 'FontSize', 14);
        elseif plotIndex == numel(invariantLabels)
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

    [rowCount, columnCount] = size(Z);
    if rowCount < 3 || columnCount < 3
        return;
    end

    sourceZ = Z;
    blankMask = ~isfinite(sourceZ);
    if ~any(blankMask(:))
        return;
    end

    visited = false(rowCount, columnCount);
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
            [rowIndex, columnIndex] = ind2sub([rowCount, columnCount], currentIndex);
            touchesBoundary = touchesBoundary || rowIndex == 1 || rowIndex == rowCount || columnIndex == 1 || columnIndex == columnCount;

            for neighborOffsetIndex = 1:8
                neighborRowIndex = rowIndex + offsets(neighborOffsetIndex, 1);
                neighborColumnIndex = columnIndex + offsets(neighborOffsetIndex, 2);
                if neighborRowIndex < 1 || neighborRowIndex > rowCount || neighborColumnIndex < 1 || neighborColumnIndex > columnCount
                    touchesBoundary = true;
                    continue;
                end
                neighborIndex = sub2ind([rowCount, columnCount], neighborRowIndex, neighborColumnIndex);
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

function [colorMin, colorMax] = finiteColorLimits(data, isSymmetric)

    if nargin < 2
        isSymmetric = false;
    end
    finiteData = data(isfinite(data));
    if isempty(finiteData)
        colorMin = -1;
        colorMax = 1;
        return;
    end

    if isSymmetric
        finiteData = abs(finiteData);
        colorMin = 0;
    else
        colorMin = min(finiteData);
    end
    colorMax = max(finiteData);
    robustMax = finitePercentile(finiteData, 99.5);
    if isfinite(robustMax) && robustMax > colorMin && colorMax > 5 * robustMax
        colorMax = robustMax;
    end
    if isSymmetric
        if colorMax <= 0
            colorMax = 1;
        end
        colorMin = -colorMax;
    elseif colorMin == colorMax
        colorMin = colorMin - 1;
        colorMax = colorMax + 1;
    end
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

function [records, nanRecordCount] = convertNaNOrbitRecordsToPad(records)

    nanRecord = isnan(double(records.Ids)) | isnan(records.orbits) | ...
        isnan(records.dtheta) | isnan(records.dphiTotal) | ...
        isnan(records.dphiVpara) | isnan(records.dTs) | ...
        isnan(records.Es) | isnan(records.Pphis) | isnan(records.Lambdas);
    nanRecordCount = sum(nanRecord);
    if nanRecordCount == 0
        return;
    end

    records.Ids(nanRecord) = int32(20251106);
    records.orbits(nanRecord) = 0.5;
    records.dtheta(nanRecord) = 0;
    records.dphiTotal(nanRecord) = 0;
    records.dphiVpara(nanRecord) = 0;
    records.dTs(nanRecord) = 0;
    records.Es(nanRecord) = 0;
    records.Pphis(nanRecord) = 0;
    records.Lambdas(nanRecord) = 0;
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

function value = statisticOrNaN(x, statisticFcn, finiteOnly)

    if finiteOnly
        x = x(isfinite(x));
    end
    if isempty(x)
        value = NaN;
    else
        value = statisticFcn(x);
    end
end

function value = getOptionValue(opt, fieldName, defaultValue)

    if isfield(opt, fieldName) && ~isempty(opt.(fieldName))
        value = opt.(fieldName);
    else
        value = defaultValue;
    end
end

function logLoaded(name, data)

    sizeText = strjoin(cellstr(compose('%d', size(data).')), ' ');
    fprintf('[load] %s: size=[%s]\n', char(name), sizeText);
end

function logSkipped(name, reason)

    fprintf('[skip] %s: %s\n', char(name), reason);
end

function printOrbitSummarySpecies(summaryStruct)

    fields = fieldnames(summaryStruct);
    if isempty(fields)
        fprintf('[orbit] 未处理任何 PhaseSpaceOrbit.bin。\n');
    else
        fprintf('[orbit] 已处理物种：%s\n', strjoin(fields, ', '));
    end
end
