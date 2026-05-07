%% cuGMEC 相空间轨道频率可视化脚本

%{

脚本流程：

(1) 用户设置
设置输入目录, normalization2D.mat 路径, Orbit.bin 所在目录, 以及需要处理的物种。

(2) 读取 normalization2D.mat
读取相空间网格。各物种 EPphiLambda 范围会在处理对应物种时读取。

(3) 逐物种处理 Orbit.bin
读取 Ion/Alpha/BeamPhaseSpaceOrbit.bin, 剔除 pad 记录, 按 trapped, para,
anti 的轨道判据分类，打印分类占比，并绘制守恒量诊断。

(4) 绘制相空间轨道频率
在本节选择绘图物种, 固定坐标切片, 频率单位, 分支方向和色表。

(5) 绘制相空间共振残差线
直接使用用户填写的真实物理环向模数 n、极向模数 m 和轨道谐波 l,
按累计相位公式绘制 residual = 0 等值线。

%}

%% (1) 用户设置

inputPath = 'C:\Users\Desktop\test';
normalizationFile = fullfile(inputPath, 'normalization2D.mat');

orbitPath = inputPath;

speciesList = {'Ion', 'Alpha', 'Beam'};

%% (2) 读取 normalization2D.mat

normData = load(normalizationFile);
[gridE, gridPphi, gridLambda] = readPhaseGrid(normData);

%% (3) 逐物种处理 Orbit.bin

phaseSpaceResults = struct();
PhaseSpaceOrbitSummary = struct();

for speciesIndex = 1:numel(speciesList)
    speciesName = speciesList{speciesIndex};
    orbitFile = fullfile(orbitPath, [speciesName 'PhaseSpaceOrbit.bin']);

    if exist(orbitFile, 'file') ~= 2
        fprintf('[skip] %s 不存在，跳过 %s。\n', orbitFile, speciesName);
        continue;
    end

    rangeField = [speciesName 'EPphiLambda'];
    phaseRange = readSpeciesRange(normData, rangeField);

    fprintf('\n============================================================\n');
    fprintf('%s phase-space orbit analysis\n', speciesName);
    fprintf('============================================================\n');

    result = analyzeSpeciesOrbit(speciesName, orbitFile, phaseRange, gridE, gridPphi, gridLambda);

    dataField = [speciesName 'PhaseSpaceData'];
    phaseSpaceResults.(dataField) = result.phaseSpaceData;
    PhaseSpaceOrbitSummary.(speciesName) = result.summary;

    printSpeciesSummary(result.summary);

    plotConservationDiagnostics(speciesName, result.diagnostics, ...
        result.E1d, result.Pphi1d, result.Lambda1d);
end

%% (4) 绘制相空间轨道频率


orbitFrequencySpecies = 'Alpha';
% 'Ion', 'Alpha', 'Beam'

orbitFrequencyFixedCoordinate = 'E';
% 'E', 'Pphi', 'Lambda'

orbitFrequencySlice = '32';
% 'middle' 或者索引字符串，例如 '32'

orbitFrequencyUnit = 'w';
% 'Hz' 或 'w'

orbitFrequencyBranch = 'para';
% 'para' 表示 para+trapped，'anti' 表示 anti+trapped

orbitFrequencyColormapIndex = 1;
% 色表序号，可选 1-2；1 为自定义色表，2 为 jet

orbitFrequencyContourCount = 20;
% 等高线数目；设为 0 时不绘制等高线

orbitFrequencyInteractive = 2;
% 0 为静态图；1 为滑块释放后更新；2 为拖动滑块时连续更新

runInteractivePlot(orbitFrequencyInteractive, ...
    @() plotPhaseSpaceOrbitFrequency(phaseSpaceResults, normData, ...
    orbitFrequencySpecies, orbitFrequencyFixedCoordinate, orbitFrequencySlice, ...
    orbitFrequencyUnit, orbitFrequencyBranch, orbitFrequencyColormapIndex, orbitFrequencyContourCount), ...
    @(dynamicUpdate) plotPhaseSpaceOrbitFrequencyInteractive(phaseSpaceResults, normData, ...
    orbitFrequencySpecies, orbitFrequencyFixedCoordinate, orbitFrequencySlice, ...
    orbitFrequencyUnit, orbitFrequencyBranch, orbitFrequencyColormapIndex, orbitFrequencyContourCount, dynamicUpdate));

%% (5) 绘制相空间共振残差线

resonanceLineSpecies = 'Alpha';
% 'Ion', 'Alpha', 'Beam'

resonanceLineFixedCoordinate = 'E';
% 'E', 'Pphi', 'Lambda'

resonanceLineSlice = '60';
% 'middle' 或者索引字符串，例如 '32'

resonanceLineBranch = 'para';
% 'para', 'anti' 或 'trapped'

resonanceLineFrequencyHz = 65e3;
% visualizeMHD.m 诊断得到的有符号模频率，单位 Hz。

resonanceLineFrequencyHzRange = [0, 150e3];
% 频率滑块范围；若需要负频率可设为例如 [-150e3, 150e3]

resonanceLineToroidalMode = 30;
% 真实物理环向模数 n；这里不再用 tubes 做换算。

resonanceLineToroidalModeRange = [1, 80];
% n 滑块范围

resonanceLinePoloidalMode = 33;
% 相位 n*zeta - m*theta - omega*t 中的正极向模数 m；仅 para/anti 使用。

resonanceLinePoloidalModeRange = [1, 80];
% m 滑块范围；trapped 分支中该滑块不参与公式

resonanceLineHarmonic = 0;
% 单个轨道谐波 l，可为负整数、0 或正整数。

resonanceLineHarmonicRange = [-20, 20];
% l 滑块范围

resonanceLineColormapIndex = 1;
% 色表序号，可选 1-2；1 为自定义色表，2 为 jet。

resonanceLineInteractive = 2;
% 0 为静态图；1 为滑块释放后更新；2 为拖动滑块时连续更新

runInteractivePlot(resonanceLineInteractive, ...
    @() plotPhaseSpaceResonanceLine(phaseSpaceResults, normData, ...
    resonanceLineSpecies, resonanceLineFixedCoordinate, resonanceLineSlice, ...
    resonanceLineBranch, resonanceLineFrequencyHz, resonanceLineToroidalMode, ...
    resonanceLinePoloidalMode, resonanceLineHarmonic, resonanceLineColormapIndex), ...
    @(dynamicUpdate) plotPhaseSpaceResonanceLineInteractive(phaseSpaceResults, normData, ...
    resonanceLineSpecies, resonanceLineFixedCoordinate, resonanceLineSlice, resonanceLineBranch, ...
    resonanceLineFrequencyHz, resonanceLineFrequencyHzRange, ...
    resonanceLineToroidalMode, resonanceLineToroidalModeRange, ...
    resonanceLinePoloidalMode, resonanceLinePoloidalModeRange, ...
    resonanceLineHarmonic, resonanceLineHarmonicRange, ...
    resonanceLineColormapIndex, dynamicUpdate));

%% local functions

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

function result = analyzeSpeciesOrbit(speciesName, orbitFile, phaseRange, gridE, gridPphi, gridLambda)

    records = orbitRecordColumns(readOrbitBinary(orbitFile));
    assertNoNaN(speciesName, records);

    nPhase = gridE * gridPphi * gridLambda;
    expectedRecords = 2 * nPhase;
    numRecords = numel(records.Ids);
    assert(numRecords == expectedRecords, ...
        '%s 记录数为 %d，但 gridE*gridPphi*gridLambda*2 = %d。', ...
        orbitFile, numRecords, expectedRecords);

    E1d = linspace(phaseRange(1), phaseRange(2), gridE);
    Pphi1d = linspace(phaseRange(3), phaseRange(4), gridPphi);
    Lambda1d = linspace(phaseRange(5), phaseRange(6), gridLambda);

    phaseSpaceData = emptyPhaseSpaceData(gridE, gridPphi, gridLambda, E1d, Pphi1d, Lambda1d, phaseRange);

    recordBranch = ones(numRecords, 1);
    recordBranch(nPhase + 1:end) = -1;

    validRecord = ~isPadRecord(records.Ids);
    localId = abs(double(records.Ids));
    validLocalId = validRecord & localId >= 0 & localId < nPhase & localId == floor(localId);

    if any(validRecord & ~validLocalId)
        badIndex = find(validRecord & ~validLocalId, 1);
        error('%s 中存在非法相空间 id：record=%d, id=%d。', orbitFile, badIndex, records.Ids(badIndex));
    end

    warnIfSignedIdOrderLooksWrong(speciesName, records.Ids, recordBranch, validRecord);

    [plusIndex, minusIndex, duplicateCounts] = buildBranchIndex(localId, recordBranch, validLocalId, nPhase);

    diagnostics = struct();
    summary = initializeSummary(speciesName, orbitFile, numRecords, nPhase, validRecord, records.orbits);
    summary.duplicatePlusIds = duplicateCounts.plus;
    summary.duplicateMinusIds = duplicateCounts.minus;

    [phaseSpaceData.trapped, diagnostics.trapped, trappedSummary] = extractTrappedOrbits( ...
        phaseSpaceData.trapped, plusIndex, minusIndex, records, gridE, gridPphi, gridLambda);
    summary.trapped = trappedSummary;

    [phaseSpaceData.para, diagnostics.para, paraSummary] = extractPassingOrbits( ...
        2.5, phaseSpaceData.para, validLocalId, localId, recordBranch, plusIndex, minusIndex, ...
        records, gridE, gridPphi, gridLambda);
    summary.para = paraSummary;

    [phaseSpaceData.anti, diagnostics.anti, antiSummary] = extractPassingOrbits( ...
        3.5, phaseSpaceData.anti, validLocalId, localId, recordBranch, plusIndex, minusIndex, ...
        records, gridE, gridPphi, gridLambda);
    summary.anti = antiSummary;

    summary.classifiedBranchCount = 2 * summary.trapped.count + summary.para.count + summary.anti.count;
    summary.classifiedRatioTotal = safeDivide(summary.classifiedBranchCount, expectedRecords);
    summary.classifiedRatioInitialized = safeDivide(summary.classifiedBranchCount, summary.initializedRecords);
    summary.unclassifiedInitializedRecords = summary.initializedRecords - summary.classifiedBranchCount;
    summary.unclassifiedInitializedRatio = safeDivide(summary.unclassifiedInitializedRecords, summary.initializedRecords);

    result = struct( ...
        'phaseSpaceData', phaseSpaceData, ...
        'diagnostics', diagnostics, ...
        'summary', summary, ...
        'E1d', E1d, ...
        'Pphi1d', Pphi1d, ...
        'Lambda1d', Lambda1d);
end

function data = readOrbitBinary(orbitFile)

    fid = fopen(orbitFile, 'rb');
    if fid == -1
        error('无法打开文件：%s', orbitFile);
    end
    closeFile = onCleanup(@() fclose(fid));

    fseek(fid, 0, 'eof');
    fileSize = ftell(fid);
    fseek(fid, 0, 'bof');

    bytesPerRecord = 9 * 8;
    if mod(fileSize, bytesPerRecord) ~= 0
        error('文件大小不是 %d 的整数倍，数据可能不完整：%s', bytesPerRecord, orbitFile);
    end

    numRecords = fileSize / bytesPerRecord;
    rawData = fread(fid, numRecords * 9, 'double=>double');
    assert(numel(rawData) == numRecords * 9, '读取文件失败或数据长度不完整：%s', orbitFile);
    data = reshape(rawData, 9, numRecords)';

    clear closeFile;
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

function phaseSpaceData = emptyPhaseSpaceData(gridE, gridPphi, gridLambda, E1d, Pphi1d, Lambda1d, phaseRange)

    phaseSpaceData = struct();
    classNames = {'trapped', 'para', 'anti'};
    for classIndex = 1:numel(classNames)
        phaseSpaceData.(classNames{classIndex}) = emptyOrbitClass(gridE, gridPphi, gridLambda);
    end

    phaseSpaceData.gridE = gridE;
    phaseSpaceData.gridPphi = gridPphi;
    phaseSpaceData.gridLambda = gridLambda;
    phaseSpaceData.EPphiLambda = phaseRange;
    phaseSpaceData.E1d = E1d;
    phaseSpaceData.Pphi1d = Pphi1d;
    phaseSpaceData.Lambda1d = Lambda1d;
end

function classData = emptyOrbitClass(gridE, gridPphi, gridLambda)

    zeroArray = zeros(gridE, gridPphi, gridLambda);
    classData = struct( ...
        'dtheta', zeroArray, ...
        'dphiTotal', zeroArray, ...
        'dphiVpara', zeroArray, ...
        'dT', zeroArray);
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

function [classData, diagnostic, summary] = extractTrappedOrbits(classData, plusIndex, minusIndex, records, ...
    gridE, gridPphi, gridLambda)

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
        averagePair(records.dTs(acceptedPlusRecord), records.dTs(acceptedMinusRecord)), ...
        gridE, gridPphi, gridLambda);

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
    recordBranch, plusIndex, minusIndex, records, gridE, gridPphi, gridLambda)

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
        records.dphiVpara(selectedRecord), records.dTs(selectedRecord), ...
        gridE, gridPphi, gridLambda);

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

function classData = fillOrbitClass(classData, localIds, dthetaValues, dphiTotalValues, dphiVparaValues, dTValues, ...
    gridE, gridPphi, gridLambda)

    if isempty(localIds)
        return;
    end

    linearIndex = localIdsToLinear(localIds, gridE, gridPphi, gridLambda);
    classData.dtheta(linearIndex) = dthetaValues(:);
    classData.dphiTotal(linearIndex) = dphiTotalValues(:);
    classData.dphiVpara(linearIndex) = dphiVparaValues(:);
    classData.dT(linearIndex) = dTValues(:);
end

function linearIndex = localIdsToLinear(localIds, gridE, gridPphi, gridLambda)

    [iE, iPphi, iLambda] = localIdsToSubscripts(localIds, gridPphi, gridLambda);
    linearIndex = sub2ind([gridE, gridPphi, gridLambda], iE, iPphi, iLambda);
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

function plotPhaseSpaceOrbitFrequency(phaseSpaceResults, normData, speciesName, fixedCoordinate, sliceText, unitText, branchText, colormapIndex, contourCount)

    if nargin < 8 || isempty(colormapIndex)
        colormapIndex = 1;
    end
    if nargin < 9 || isempty(contourCount)
        contourCount = 10;
    end
    colormapIndex = validateColormapIndex(colormapIndex);
    contourCount = validateContourCount(contourCount);

    [phaseSpaceData, speciesLabel] = resolveSpeciesPhaseSpaceData(phaseSpaceResults, speciesName);
    dim = phaseCoordinateToDimension(fixedCoordinate);
    idx = parseOrbitFrequencySlice(sliceText, dim, phaseSpaceData);
    branchName = normalizeFrequencyBranch(branchText);
    [unitScale, colorbarLabel] = frequencyUnitScale(unitText);

    L0 = readPositiveScalar(normData, 'L0');
    VA0 = readPositiveScalar(normData, 'VA0');
    frequency = calculateOrbitFrequencyField(phaseSpaceData, branchName, unitScale, L0, VA0);

    [Z, xVec, yVec, xlabelText, ylabelText, titleText] = sliceField( ...
        dim, idx, frequency, phaseSpaceData.E1d, phaseSpaceData.Pphi1d, phaseSpaceData.Lambda1d);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s %s has no orbit-frequency data in this slice.\n', speciesLabel, branchName);
        return;
    end

    Z = fillEnclosedBlankRegions(Z);
    figureName = sprintf('%s %s orbit frequency', speciesLabel, branchName);
    drawPhaseSpaceColorMap(Z, xVec, yVec, xlabelText, ylabelText, titleText, ...
        figureName, colorbarLabel, colormapIndex, contourCount);
end

function plotPhaseSpaceOrbitFrequencyInteractive(phaseSpaceResults, normData, speciesName, fixedCoordinate, sliceText, ...
    unitText, branchText, colormapIndex, contourCount, dynamicUpdate)

    if nargin < 8 || isempty(colormapIndex)
        colormapIndex = 1;
    end
    if nargin < 9 || isempty(contourCount)
        contourCount = 10;
    end

    colormapIndex = validateColormapIndex(colormapIndex);
    contourCount = validateContourCount(contourCount);

    [phaseSpaceData, speciesLabel] = resolveSpeciesPhaseSpaceData(phaseSpaceResults, speciesName);
    dim = phaseCoordinateToDimension(fixedCoordinate);
    idx0 = parseOrbitFrequencySlice(sliceText, dim, phaseSpaceData);
    branchName = normalizeFrequencyBranch(branchText);
    [unitScale, colorbarLabel] = frequencyUnitScale(unitText);

    L0 = readPositiveScalar(normData, 'L0');
    VA0 = readPositiveScalar(normData, 'VA0');
    frequency = calculateOrbitFrequencyField(phaseSpaceData, branchName, unitScale, L0, VA0);

    nSlice = phaseDimensionSize(dim, phaseSpaceData);
    maxContourCount = max(40, 2 * contourCount);
    controls = [ ...
        integerSliderControl('sliceIndex', [char(fixedCoordinate) ' index'], idx0, 1, nSlice), ...
        integerSliderControl('contourCount', 'contours', contourCount, 0, maxContourCount)];

    figureName = sprintf('%s %s orbit frequency', speciesLabel, branchName);
    plotInteractivePhaseSpaceMapWithSliders(figureName, controls, dynamicUpdate, @computePlotData, @renderPlotData);

    function plotData = computePlotData(values)
        [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = sliceField( ...
            dim, values.sliceIndex, frequency, phaseSpaceData.E1d, phaseSpaceData.Pphi1d, phaseSpaceData.Lambda1d);

        hasData = any(isfinite(Z(:)));
        if hasData
            Z = fillEnclosedBlankRegions(Z);
        end

        plotTitle = sprintf('%s %s, %s', speciesLabel, branchName, sliceTitle);
        statusText = sprintf('%s %s | slice=%d/%d | contours=%d | finite=%d', ...
            speciesLabel, branchName, values.sliceIndex, nSlice, values.contourCount, sum(isfinite(Z(:))));

        plotData = phaseSpaceMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, plotTitle, ...
            colorbarLabel, colormapIndex, statusText);
        plotData.contourCount = values.contourCount;
    end

    function renderPlotData(axHandle, plotData)
        renderPhaseSpaceColorMap(axHandle, plotData.Z, plotData.xVec, plotData.yVec, ...
            plotData.xlabelText, plotData.ylabelText, plotData.titleText, ...
            plotData.colorbarLabel, plotData.colormapIndex, plotData.contourCount);
    end
end

function frequency = calculateOrbitFrequencyField(phaseSpaceData, branchName, unitScale, L0, VA0)

    dT = phaseSpaceData.(branchName).dT;
    trappedDT = phaseSpaceData.trapped.dT;
    fillFromTrapped = (dT == 0 | ~isfinite(dT)) & isfinite(trappedDT) & trappedDT > 0;
    dT(fillFromTrapped) = trappedDT(fillFromTrapped);
    dT(~isfinite(dT) | dT <= 0) = NaN;

    frequency = unitScale ./ (dT .* (L0 / VA0));
    frequency(~isfinite(frequency)) = NaN;
end

function plotPhaseSpaceResonanceLine(phaseSpaceResults, normData, speciesName, fixedCoordinate, sliceText, branchText, ...
    modeFrequencyHz, toroidalMode, poloidalMode, harmonic, colormapIndex)

    if nargin < 11 || isempty(colormapIndex)
        colormapIndex = 1;
    end

    colormapIndex = validateColormapIndex(colormapIndex);
    modeFrequencyHz = validateFiniteScalar(modeFrequencyHz, 'resonanceLineFrequencyHz');
    toroidalMode = readPositiveIntegerScalar(toroidalMode, 'resonanceLineToroidalMode');
    poloidalMode = readPositiveIntegerScalar(poloidalMode, 'resonanceLinePoloidalMode');
    harmonic = validateIntegerScalar(harmonic, 'resonanceLineHarmonic');

    [phaseSpaceData, speciesLabel] = resolveSpeciesPhaseSpaceData(phaseSpaceResults, speciesName);
    dim = phaseCoordinateToDimension(fixedCoordinate);
    idx = parseOrbitFrequencySlice(sliceText, dim, phaseSpaceData);
    branchName = normalizeResonanceBranch(branchText);

    L0 = readPositiveScalar(normData, 'L0');
    VA0 = readPositiveScalar(normData, 'VA0');
    [residualHz, physicalToroidalMode] = calculateResonanceResidualField( ...
        phaseSpaceData, branchName, modeFrequencyHz, toroidalMode, poloidalMode, harmonic, L0, VA0);

    [Z, xVec, yVec, xlabelText, ylabelText, titleText] = sliceField( ...
        dim, idx, residualHz, phaseSpaceData.E1d, phaseSpaceData.Pphi1d, phaseSpaceData.Lambda1d);

    if ~any(isfinite(Z(:)))
        fprintf('[plot] %s %s has no resonance-residual data in this slice.\n', speciesLabel, branchName);
        return;
    end

    figureName = sprintf('%s %s resonance residual', speciesLabel, branchName);
    if strcmp(branchName, 'trapped')
        residualTitle = sprintf('%s %s, %s, $n=%d$, $l=%d$, $f=%.6g\\,\\mathrm{Hz}$', ...
            speciesLabel, branchName, titleText, physicalToroidalMode, harmonic, modeFrequencyHz);
    else
        residualTitle = sprintf('%s %s, %s, $n=%d$, $m=%d$, $l=%d$, $f=%.6g\\,\\mathrm{Hz}$', ...
            speciesLabel, branchName, titleText, physicalToroidalMode, poloidalMode, harmonic, modeFrequencyHz);
    end
    drawPhaseSpaceResonanceMap(Z, xVec, yVec, xlabelText, ylabelText, residualTitle, ...
        figureName, '$\Delta f_{\mathrm{res}}/\mathrm{Hz}$', colormapIndex);
end

function plotPhaseSpaceResonanceLineInteractive(phaseSpaceResults, normData, speciesName, fixedCoordinate, sliceText, branchText, ...
    modeFrequencyHz, modeFrequencyHzRange, toroidalMode, toroidalModeRange, poloidalMode, poloidalModeRange, ...
    harmonic, harmonicRange, colormapIndex, dynamicUpdate)

    colormapIndex = validateColormapIndex(colormapIndex);
    modeFrequencyHz = validateFiniteScalar(modeFrequencyHz, 'resonanceLineFrequencyHz');
    toroidalMode = readPositiveIntegerScalar(toroidalMode, 'resonanceLineToroidalMode');
    poloidalMode = readPositiveIntegerScalar(poloidalMode, 'resonanceLinePoloidalMode');
    harmonic = validateIntegerScalar(harmonic, 'resonanceLineHarmonic');

    modeFrequencyHzRange = normalizeSliderRange(modeFrequencyHzRange, modeFrequencyHz, 'resonanceLineFrequencyHzRange', false);
    toroidalModeRange = normalizeSliderRange(toroidalModeRange, toroidalMode, 'resonanceLineToroidalModeRange', true);
    poloidalModeRange = normalizeSliderRange(poloidalModeRange, poloidalMode, 'resonanceLinePoloidalModeRange', true);
    harmonicRange = normalizeSliderRange(harmonicRange, harmonic, 'resonanceLineHarmonicRange', true);
    toroidalModeRange(1) = max(1, toroidalModeRange(1));
    poloidalModeRange(1) = max(1, poloidalModeRange(1));
    if toroidalModeRange(1) == toroidalModeRange(2)
        toroidalModeRange(2) = toroidalModeRange(1) + 1;
    end
    if poloidalModeRange(1) == poloidalModeRange(2)
        poloidalModeRange(2) = poloidalModeRange(1) + 1;
    end

    [phaseSpaceData, speciesLabel] = resolveSpeciesPhaseSpaceData(phaseSpaceResults, speciesName);
    dim = phaseCoordinateToDimension(fixedCoordinate);
    idx0 = parseOrbitFrequencySlice(sliceText, dim, phaseSpaceData);
    branchName = normalizeResonanceBranch(branchText);

    L0 = readPositiveScalar(normData, 'L0');
    VA0 = readPositiveScalar(normData, 'VA0');
    nSlice = phaseDimensionSize(dim, phaseSpaceData);

    controls = [ ...
        integerSliderControl('sliceIndex', [char(fixedCoordinate) ' index'], idx0, 1, nSlice), ...
        numericSliderControl('modeFrequencyHz', 'f/Hz', modeFrequencyHz, modeFrequencyHzRange(1), modeFrequencyHzRange(2), 501), ...
        integerSliderControl('toroidalMode', 'n', toroidalMode, toroidalModeRange(1), toroidalModeRange(2)), ...
        integerSliderControl('poloidalMode', 'm', poloidalMode, poloidalModeRange(1), poloidalModeRange(2)), ...
        integerSliderControl('harmonic', 'l', harmonic, harmonicRange(1), harmonicRange(2))];

    figureName = sprintf('%s %s resonance residual', speciesLabel, branchName);
    plotInteractivePhaseSpaceMapWithSliders(figureName, controls, dynamicUpdate, @computePlotData, @renderPlotData);

    function plotData = computePlotData(values)
        [residualHz, physicalToroidalMode] = calculateResonanceResidualField( ...
            phaseSpaceData, branchName, values.modeFrequencyHz, values.toroidalMode, ...
            values.poloidalMode, values.harmonic, L0, VA0);

        [Z, xVec, yVec, xlabelText, ylabelText, sliceTitle] = sliceField( ...
            dim, values.sliceIndex, residualHz, phaseSpaceData.E1d, phaseSpaceData.Pphi1d, phaseSpaceData.Lambda1d);

        if strcmp(branchName, 'trapped')
            plotTitle = sprintf('%s %s, %s, $n=%d$, $l=%d$, $f=%.6g\\,\\mathrm{Hz}$', ...
                speciesLabel, branchName, sliceTitle, physicalToroidalMode, values.harmonic, values.modeFrequencyHz);
        else
            plotTitle = sprintf('%s %s, %s, $n=%d$, $m=%d$, $l=%d$, $f=%.6g\\,\\mathrm{Hz}$', ...
                speciesLabel, branchName, sliceTitle, physicalToroidalMode, values.poloidalMode, values.harmonic, values.modeFrequencyHz);
        end

        hasZeroContour = residualHasZeroContour(Z);
        if hasZeroContour
            contourText = 'yes';
        else
            contourText = 'no';
        end

        statusText = sprintf('%s %s | slice=%d/%d | f=%.6g Hz | n=%d | m=%d | l=%d | zero contour=%s', ...
            speciesLabel, branchName, values.sliceIndex, nSlice, values.modeFrequencyHz, ...
            physicalToroidalMode, values.poloidalMode, values.harmonic, contourText);

        plotData = phaseSpaceMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, plotTitle, ...
            '$\Delta f_{\mathrm{res}}/\mathrm{Hz}$', colormapIndex, statusText);
    end

    function renderPlotData(axHandle, plotData)
        renderPhaseSpaceResonanceMap(axHandle, plotData.Z, plotData.xVec, plotData.yVec, ...
            plotData.xlabelText, plotData.ylabelText, plotData.titleText, ...
            plotData.colorbarLabel, plotData.colormapIndex);
    end
end

function [residualHz, physicalToroidalMode] = calculateResonanceResidualField( ...
    phaseSpaceData, branchName, modeFrequencyHz, toroidalMode, poloidalMode, harmonic, L0, VA0)

    physicalToroidalMode = toroidalMode;
    classData = phaseSpaceData.(branchName);
    orbitTime = classData.dT .* (L0 / VA0);
    validOrbit = isfinite(orbitTime) & orbitTime > 0;

    % dT is normalized by L0/VA0; the residual below uses physical seconds.
    % Resonance is evaluated from the accumulated phase of n*zeta - m*theta - omega*t.
    phaseRate = nan(size(orbitTime));
    switch branchName
        case {'para', 'anti'}
            phaseAdvance = physicalToroidalMode .* classData.dphiTotal - poloidalMode .* classData.dtheta + 2 * pi * harmonic;
            phaseRate(validOrbit) = phaseAdvance(validOrbit) ./ orbitTime(validOrbit);
        case 'trapped'
            % Trapped particles use the precession phase; the poloidal mode m is averaged out.
            precessionFrequency = nan(size(orbitTime));
            bounceFrequency = nan(size(orbitTime));
            precessionPhase = classData.dphiTotal - classData.dphiVpara;
            precessionFrequency(validOrbit) = precessionPhase(validOrbit) ./ orbitTime(validOrbit);
            bounceFrequency(validOrbit) = 2 * pi ./ orbitTime(validOrbit);
            phaseRate(validOrbit) = physicalToroidalMode .* precessionFrequency(validOrbit) + ...
                harmonic .* bounceFrequency(validOrbit);
        otherwise
            error('Unsupported resonance branch "%s".', branchName);
    end

    modeAngularFrequency = 2 * pi * modeFrequencyHz;
    % Return signed residual/(2*pi) in Hz so residual = 0 can be contoured directly.
    residualHz = (modeAngularFrequency - phaseRate) ./ (2 * pi);
    residualHz(~isfinite(residualHz)) = NaN;
end

function drawPhaseSpaceResonanceMap(Z, xVec, yVec, xlabelText, ylabelText, titleText, figureName, colorbarLabel, colormapIndex)

    figHandle = figure('Name', figureName, 'Position', [100, 100, 900, 760]);
    axHandle = axes('Parent', figHandle);
    renderPhaseSpaceResonanceMap(axHandle, Z, xVec, yVec, xlabelText, ylabelText, titleText, colorbarLabel, colormapIndex);
end

function renderPhaseSpaceResonanceMap(axHandle, Z, xVec, yVec, xlabelText, ylabelText, titleText, colorbarLabel, colormapIndex)

    resetPhaseSpaceAxes(axHandle);
    [X, Y] = meshgrid(xVec, yVec);
    pcolor(axHandle, X, Y, Z);
    shading(axHandle, 'interp');
    hold(axHandle, 'on');

    [colorMin, colorMax] = symmetricFiniteColorLimits(Z);
    clim(axHandle, [colorMin, colorMax]);
    colormap(axHandle, mhdFieldColormap(colormapIndex, 256));
    cb = colorbar(axHandle);
    cb.FontName = 'Times New Roman';
    cb.FontSize = 14;
    ylabel(cb, colorbarLabel, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    if residualHasZeroContour(Z)
        contour(axHandle, X, Y, Z, [0, 0], 'EdgeColor', 'k', 'LineWidth', 1.5);
    end

    xlabel(axHandle, xlabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(axHandle, ylabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(axHandle, titleText, ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    axis(axHandle, 'tight');
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14, 'Color', 'white', 'Layer', 'top');
    hold(axHandle, 'off');
end

function drawPhaseSpaceColorMap(Z, xVec, yVec, xlabelText, ylabelText, titleText, figureName, colorbarLabel, colormapIndex, contourCount)

    figHandle = figure('Name', figureName, 'Position', [100, 100, 900, 760]);
    axHandle = axes('Parent', figHandle);
    renderPhaseSpaceColorMap(axHandle, Z, xVec, yVec, xlabelText, ylabelText, titleText, colorbarLabel, colormapIndex, contourCount);
end

function renderPhaseSpaceColorMap(axHandle, Z, xVec, yVec, xlabelText, ylabelText, titleText, colorbarLabel, colormapIndex, contourCount)

    resetPhaseSpaceAxes(axHandle);
    [X, Y] = meshgrid(xVec, yVec);
    pcolor(axHandle, X, Y, Z);
    shading(axHandle, 'interp');
    hold(axHandle, 'on');

    finiteZ = Z(isfinite(Z));
    if contourCount > 0 && ~isempty(finiteZ) && min(finiteZ) < max(finiteZ)
        contour(axHandle, X, Y, Z, contourCount, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end

    [colorMin, colorMax] = finiteColorLimits(Z);
    clim(axHandle, [colorMin, colorMax]);
    colormap(axHandle, mhdFieldColormap(colormapIndex, 256));
    cb = colorbar(axHandle);
    cb.FontName = 'Times New Roman';
    cb.FontSize = 14;
    ylabel(cb, colorbarLabel, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    xlabel(axHandle, xlabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(axHandle, ylabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(axHandle, titleText, ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    axis(axHandle, 'tight');
    set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 14, 'Color', 'white', 'Layer', 'top');
    hold(axHandle, 'off');
end

function resetPhaseSpaceAxes(axHandle)

    figHandle = ancestor(axHandle, 'figure');
    delete(findall(figHandle, 'Type', 'ColorBar'));
    cla(axHandle);
end

function plotData = phaseSpaceMapPlotData(Z, xVec, yVec, xlabelText, ylabelText, titleText, colorbarLabel, colormapIndex, statusText)

    plotData = struct( ...
        'Z', Z, ...
        'xVec', xVec, ...
        'yVec', yVec, ...
        'xlabelText', xlabelText, ...
        'ylabelText', ylabelText, ...
        'titleText', titleText, ...
        'colorbarLabel', colorbarLabel, ...
        'colormapIndex', colormapIndex, ...
        'status', statusText);
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

function plotInteractivePhaseSpaceMapWithSliders(figName, controls, dynamicUpdate, computePlotData, renderPlotData)

    figHandle = figure('Name', figName, 'Color', 'w', 'Position', [120, 80, 980, 760]);
    nControl = numel(controls);
    controlBottom = 0.030;
    controlSpacing = 0.041;
    statusBottom = controlBottom + controlSpacing * max(nControl, 1) + 0.010;
    axesBottom = statusBottom + 0.075;
    axHandle = axes('Parent', figHandle, 'Units', 'normalized', ...
        'Position', [0.10, axesBottom, 0.82, 0.94 - axesBottom]);
    statusText = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.10, statusBottom, 0.82, 0.035], 'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', 'FontName', 'Times New Roman', 'FontSize', 12);

    sliderLabels = gobjects(nControl, 1);
    sliders = gobjects(nControl, 1);
    for iControl = 1:nControl
        yPos = controlBottom + controlSpacing * (nControl - iControl);
        sliderLabels(iControl) = uicontrol(figHandle, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.10, yPos - 0.008, 0.19, 0.032], ...
            'String', sliderLabelText(controls(iControl), controls(iControl).value), ...
            'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', 12);
        sliders(iControl) = uicontrol(figHandle, 'Style', 'slider', 'Units', 'normalized', ...
            'Position', [0.30, yPos, 0.62, 0.024], ...
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
            warning('visualizePhaseSpaceOrbit:DynamicSliderUnavailable', ...
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

        plotData = computePlotData(values);
        renderPlotData(axHandle, plotData);
        set(statusText, 'String', plotData.status);
    end
end

function [phaseSpaceData, speciesLabel] = resolveSpeciesPhaseSpaceData(phaseSpaceResults, speciesName)

    speciesLabel = strtrim(char(speciesName));
    targetField = [speciesLabel 'PhaseSpaceData'];

    if isfield(phaseSpaceResults, targetField)
        phaseSpaceData = phaseSpaceResults.(targetField);
        return;
    end

    availableFields = fieldnames(phaseSpaceResults);
    for fieldIndex = 1:numel(availableFields)
        if strcmpi(availableFields{fieldIndex}, targetField)
            targetField = availableFields{fieldIndex};
            phaseSpaceData = phaseSpaceResults.(targetField);
            speciesLabel = strrep(targetField, 'PhaseSpaceData', '');
            return;
        end
    end

    error('No phase-space orbit data found for species "%s".', speciesLabel);
end

function dim = phaseCoordinateToDimension(fixedCoordinate)

    key = lower(strtrim(char(fixedCoordinate)));
    key = strrep(key, '_', '');
    key = strrep(key, '{', '');
    key = strrep(key, '}', '');
    key = strrep(key, '\', '');

    switch key
        case {'e', 'energy'}
            dim = 1;
        case {'pphi', 'pvarphi'}
            dim = 2;
        case {'lambda'}
            dim = 3;
        otherwise
            error('fixed coordinate must be one of "E", "Pphi", or "Lambda".');
    end
end

function idx = parseOrbitFrequencySlice(sliceText, dim, phaseSpaceData)

    if isnumeric(sliceText)
        idx = double(sliceText);
    else
        key = lower(strtrim(char(sliceText)));
        if ismember(key, {'middle', 'mid', 'center', 'centre'})
            idx = defaultSliceIndex(dim, phaseSpaceData.gridE, phaseSpaceData.gridPphi, phaseSpaceData.gridLambda);
        else
            idx = str2double(key);
        end
    end

    switch dim
        case 1
            idx = clampIndex(idx, phaseSpaceData.gridE, 'E');
        case 2
            idx = clampIndex(idx, phaseSpaceData.gridPphi, 'Pphi');
        case 3
            idx = clampIndex(idx, phaseSpaceData.gridLambda, 'Lambda');
    end
end

function nSlice = phaseDimensionSize(dim, phaseSpaceData)

    switch dim
        case 1
            nSlice = phaseSpaceData.gridE;
        case 2
            nSlice = phaseSpaceData.gridPphi;
        case 3
            nSlice = phaseSpaceData.gridLambda;
        otherwise
            error('slice dimension must be 1, 2, or 3.');
    end
end

function branchName = normalizeFrequencyBranch(branchText)

    branchName = lower(strtrim(char(branchText)));
    if ~ismember(branchName, {'para', 'anti'})
        error('frequency branch must be "para" or "anti".');
    end
end

function branchName = normalizeResonanceBranch(branchText)

    branchName = lower(strtrim(char(branchText)));
    if ~ismember(branchName, {'para', 'anti', 'trapped'})
        error('resonance branch must be "para", "anti", or "trapped".');
    end
end

function [unitScale, colorbarLabel] = frequencyUnitScale(unitText)

    key = lower(strtrim(char(unitText)));
    switch key
        case 'hz'
            unitScale = 1;
            colorbarLabel = '$f/\mathrm{Hz}$';
        case {'w', 'omega'}
            unitScale = 2 * pi;
            colorbarLabel = '$\omega/\mathrm{rad}$';
        otherwise
            error('frequency unit must be "Hz" or "w".');
    end
end

function value = readPositiveScalar(normData, fieldName)

    assert(isfield(normData, fieldName), 'normalization2D.mat is missing "%s".', fieldName);
    value = double(normData.(fieldName));
    assert(isscalar(value) && isfinite(value) && value > 0, '%s must be a positive scalar.', fieldName);
end

function value = validateFiniteScalar(rawValue, fieldName)

    value = double(rawValue);
    assert(isscalar(value) && isfinite(value), '%s must be a finite scalar.', fieldName);
end

function value = validateIntegerScalar(rawValue, fieldName)

    value = double(rawValue);
    assert(isscalar(value) && isfinite(value) && value == floor(value), ...
        '%s must be an integer scalar.', fieldName);
end

function range = normalizeSliderRange(rawRange, initialValue, fieldName, isInteger)

    range = reshape(double(rawRange), 1, []);
    assert(numel(range) == 2 && all(isfinite(range)), '%s must be [min max].', fieldName);
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

function value = clampSliderValue(value, control)

    if control.isInteger
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

    if control.isInteger
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

    errorE = mixedConservationError(diagnostic.E, baseE);
    errorPphi = mixedConservationError(diagnostic.Pphi, basePphi);
    errorLambda = mixedConservationError(diagnostic.Lambda, baseLambda);

    figure('Name', [speciesName ' ' classLabel ' 守恒量误差'], 'Position', [100, 100, 900, 760]);
    invariantLabels = {'E', 'Pphi', 'Lambda'};
    invariantErrors = {errorE, errorPphi, errorLambda};

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

function [Z, xVec, yVec, xlabelText, ylabelText, titleText] = sliceField(dim, idx, field, E1d, Pphi1d, Lambda1d)

    switch dim
        case 1
            idx = clampIndex(idx, numel(E1d), 'E');
            Z = reshape(field(idx, :, :), numel(Pphi1d), numel(Lambda1d));
            xVec = Lambda1d;
            yVec = -Pphi1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$P_{\varphi}$';
            titleText = sprintf('$E = %.6g$', E1d(idx));
        case 2
            idx = clampIndex(idx, numel(Pphi1d), 'Pphi');
            Z = reshape(field(:, idx, :), numel(E1d), numel(Lambda1d));
            xVec = Lambda1d;
            yVec = E1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$E$';
            titleText = sprintf('$P_{\\varphi} = %.6g$', -Pphi1d(idx));
        case 3
            idx = clampIndex(idx, numel(Lambda1d), 'Lambda');
            Z = reshape(field(:, :, idx), numel(E1d), numel(Pphi1d));
            xVec = -Pphi1d;
            yVec = E1d;
            xlabelText = '$P_{\varphi}$';
            ylabelText = '$E$';
            titleText = sprintf('$\\Lambda = %.6g$', Lambda1d(idx));
        otherwise
            error('slice dimension must be 1, 2, or 3.');
    end
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

function [colorMin, colorMax] = finiteColorLimits(data)

    finiteData = data(isfinite(data));
    if isempty(finiteData)
        colorMin = -1;
        colorMax = 1;
    else
        colorMin = min(finiteData);
        colorMax = max(finiteData);
        if colorMin == colorMax
            colorMin = colorMin - 1;
            colorMax = colorMax + 1;
        end
    end
end

function [colorMin, colorMax] = symmetricFiniteColorLimits(data)

    finiteData = data(isfinite(data));
    if isempty(finiteData)
        colorMin = -1;
        colorMax = 1;
        return;
    end

    colorLimit = max(abs(finiteData));
    if colorLimit <= 0
        colorLimit = 1;
    end

    colorMin = -colorLimit;
    colorMax = colorLimit;
end

function colormapIndex = validateColormapIndex(colormapIndex)

    assert(isscalar(colormapIndex) && isfinite(colormapIndex) && ...
        colormapIndex == floor(colormapIndex) && colormapIndex >= 1 && colormapIndex <= 2, ...
        'colormapIndex must be an integer from 1 to 2.');
end

function contourCount = validateContourCount(contourCount)

    assert(isscalar(contourCount) && isfinite(contourCount) && ...
        contourCount == floor(contourCount) && contourCount >= 0, ...
        'contourCount must be a nonnegative integer.');
end

function cmap = mhdFieldColormap(colormapIndex, nColor)

    switch colormapIndex
        case 1
            cmap = resampleColormap(customOrbitColormap(), nColor);
        case 2
            cmap = jet(nColor);
        otherwise
            error('colormapIndex must be an integer from 1 to 2.');
    end
end

function cmap = resampleColormap(baseCmap, nColor)

    if nColor == size(baseCmap, 1)
        cmap = baseCmap;
    else
        cmap = interp1(linspace(0, 1, size(baseCmap, 1)), baseCmap, linspace(0, 1, nColor), 'linear');
    end
end

function printSpeciesSummary(summary)

    totalRecords = summary.expectedRecords;
    initialized = summary.initializedRecords;

    fprintf('pad      : %d (%.2f%% of total)\n', summary.orbitCounts.pad, 100 * safeDivide(summary.orbitCounts.pad, totalRecords));
    fprintf('effective: %d (%.2f%% of total)\n', initialized, 100 * summary.initializedRatio);
    fprintf('\n');
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

function summary = initializeSummary(speciesName, orbitFile, numRecords, nPhase, validRecord, orbits)

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

function [gridE, gridPphi, gridLambda] = readPhaseGrid(normData)

    validateNormalizationMetadata(normData);
    gridE = readPositiveIntegerScalar(normData.gridE, 'gridE');
    gridPphi = readPositiveIntegerScalar(normData.gridPphi, 'gridPphi');
    gridLambda = readPositiveIntegerScalar(normData.gridLambda, 'gridLambda');
end

function validateNormalizationMetadata(normData)

    requiredFields = {'gridE', 'gridPphi', 'gridLambda'};
    for fieldIndex = 1:numel(requiredFields)
        assert(isfield(normData, requiredFields{fieldIndex}), ...
            'normalization2D.mat 缺少字段 "%s"。请先运行新版 generatePhaseSpaceMapping2D.m。', ...
            requiredFields{fieldIndex});
    end
end

function value = readPositiveIntegerScalar(rawValue, fieldName)

    value = double(rawValue);
    assert(isscalar(value) && isfinite(value) && value == floor(value) && value > 0, ...
        '%s 必须是正整数标量。', fieldName);
end

function phaseRange = readSpeciesRange(normData, rangeField)

    assert(isfield(normData, rangeField), ...
        'normalization2D.mat 缺少字段 "%s"。请先用 generatePhaseSpaceMapping2D.m 生成该物种映射。', rangeField);

    phaseRange = reshape(double(normData.(rangeField)), 1, []);
    assert(numel(phaseRange) == 6 && all(isfinite(phaseRange)), ...
        '%s 必须是包含 6 个有限数的数组：[minE maxE minPphi maxPphi minLambda maxLambda]。', rangeField);
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

function idx = defaultSliceIndex(dim, gridE, gridPphi, gridLambda)

    switch dim
        case 1
            idx = max(1, round(gridE / 2));
        case 2
            idx = max(1, round(gridPphi / 2));
        case 3
            idx = max(1, round(gridLambda / 2));
        otherwise
            error('slice dimension must be 1, 2, or 3.');
    end
end

function idx = clampIndex(idx, maxIndex, label)

    assert(isscalar(idx) && isfinite(idx) && idx == floor(idx), '%s 切片下标必须是整数。', label);
    if idx < 1 || idx > maxIndex
        error('%s 切片下标越界：idx=%d，有效范围为 [1, %d]。', label, idx, maxIndex);
    end
end

function cmap = customOrbitColormap()

    cmap = [ ...
        0.163302000000000, 0.119982000000000, 0.793530000000000; ...
        0.171177000000000, 0.136715000000000, 0.802131000000000; ...
        0.179052000000000, 0.153448000000000, 0.810732000000000; ...
        0.186927000000000, 0.170181000000000, 0.819333000000000; ...
        0.194801000000000, 0.186914000000000, 0.827934000000000; ...
        0.202676000000000, 0.203647000000000, 0.836535000000000; ...
        0.210551000000000, 0.220380000000000, 0.845136000000000; ...
        0.218426000000000, 0.237114000000000, 0.853737000000000; ...
        0.226301000000000, 0.253847000000000, 0.862338000000000; ...
        0.234176000000000, 0.270580000000000, 0.870939000000000; ...
        0.242051000000000, 0.287313000000000, 0.879540000000000; ...
        0.249926000000000, 0.304046000000000, 0.888142000000000; ...
        0.260241000000000, 0.322242000000000, 0.894622000000000; ...
        0.273484000000000, 0.342192000000000, 0.898558000000000; ...
        0.286727000000000, 0.362143000000000, 0.902494000000000; ...
        0.299970000000000, 0.382094000000000, 0.906430000000000; ...
        0.313213000000000, 0.402044000000000, 0.910366000000000; ...
        0.326456000000000, 0.421995000000000, 0.914302000000000; ...
        0.339699000000000, 0.441946000000000, 0.918238000000000; ...
        0.352943000000000, 0.461896000000000, 0.922174000000000; ...
        0.366186000000000, 0.481847000000000, 0.926109000000000; ...
        0.379429000000000, 0.501798000000000, 0.930045000000000; ...
        0.392672000000000, 0.521749000000000, 0.933981000000000; ...
        0.405915000000000, 0.541699000000000, 0.937917000000000; ...
        0.420158000000000, 0.558463000000000, 0.939249000000000; ...
        0.434501000000000, 0.574908000000000, 0.940321000000000; ...
        0.448844000000000, 0.591354000000000, 0.941393000000000; ...
        0.463187000000000, 0.607799000000000, 0.942465000000000; ...
        0.477530000000000, 0.624244000000000, 0.943537000000000; ...
        0.491873000000000, 0.640689000000000, 0.944609000000000; ...
        0.506216000000000, 0.657134000000000, 0.945681000000000; ...
        0.520559000000000, 0.673579000000000, 0.946752000000000; ...
        0.534902000000000, 0.690025000000000, 0.947824000000000; ...
        0.549245000000000, 0.706470000000000, 0.948896000000000; ...
        0.563588000000000, 0.722915000000000, 0.949968000000000; ...
        0.577366000000000, 0.737218000000000, 0.949957000000000; ...
        0.590155000000000, 0.747772000000000, 0.948050000000000; ...
        0.602944000000000, 0.758326000000000, 0.946143000000000; ...
        0.615734000000000, 0.768881000000000, 0.944236000000000; ...
        0.628523000000000, 0.779435000000000, 0.942329000000000; ...
        0.641312000000000, 0.789989000000000, 0.940423000000000; ...
        0.654102000000000, 0.800544000000000, 0.938516000000000; ...
        0.666891000000000, 0.811098000000000, 0.936609000000000; ...
        0.679681000000000, 0.821652000000000, 0.934702000000000; ...
        0.692470000000000, 0.832206000000000, 0.932795000000000; ...
        0.705259000000000, 0.842761000000000, 0.930889000000000; ...
        0.718049000000000, 0.853315000000000, 0.928982000000000; ...
        0.728215000000000, 0.858656000000000, 0.924361000000000; ...
        0.737798000000000, 0.862838000000000, 0.919138000000000; ...
        0.747381000000000, 0.867020000000000, 0.913914000000000; ...
        0.756965000000000, 0.871202000000000, 0.908690000000000; ...
        0.766548000000000, 0.875384000000000, 0.903467000000000; ...
        0.776131000000000, 0.879566000000000, 0.898243000000000; ...
        0.785714000000000, 0.883748000000000, 0.893019000000000; ...
        0.795298000000000, 0.887930000000000, 0.887796000000000; ...
        0.804881000000000, 0.892112000000000, 0.882572000000000; ...
        0.814464000000000, 0.896294000000000, 0.877349000000000; ...
        0.824047000000000, 0.900476000000000, 0.872125000000000; ...
        0.832515000000000, 0.902966000000000, 0.866069000000000; ...
        0.838010000000000, 0.900941000000000, 0.857794000000000; ...
        0.843504000000000, 0.898916000000000, 0.849520000000000; ...
        0.848999000000000, 0.896891000000000, 0.841245000000000; ...
        0.854493000000000, 0.894866000000000, 0.832970000000000; ...
        0.859987000000000, 0.892841000000000, 0.824695000000000; ...
        0.865482000000000, 0.890816000000000, 0.816421000000000; ...
        0.870976000000000, 0.888791000000000, 0.808146000000000; ...
        0.876470000000000, 0.886766000000000, 0.799871000000000; ...
        0.881965000000000, 0.884741000000000, 0.791596000000000; ...
        0.887459000000000, 0.882716000000000, 0.783322000000000; ...
        0.892954000000000, 0.880691000000000, 0.775047000000000; ...
        0.895305000000000, 0.874424000000000, 0.765238000000000; ...
        0.896479000000000, 0.866566000000000, 0.754854000000000; ...
        0.897652000000000, 0.858709000000000, 0.744470000000000; ...
        0.898825000000000, 0.850851000000000, 0.734086000000000; ...
        0.899999000000000, 0.842993000000000, 0.723702000000000; ...
        0.901172000000000, 0.835135000000000, 0.713318000000000; ...
        0.902346000000000, 0.827277000000000, 0.702935000000000; ...
        0.903519000000000, 0.819420000000000, 0.692551000000000; ...
        0.904692000000000, 0.811562000000000, 0.682167000000000; ...
        0.905866000000000, 0.803704000000000, 0.671783000000000; ...
        0.907039000000000, 0.795846000000000, 0.661399000000000; ...
        0.907472000000000, 0.787052000000000, 0.650848000000000; ...
        0.904570000000000, 0.774046000000000, 0.639547000000000; ...
        0.901669000000000, 0.761040000000000, 0.628246000000000; ...
        0.898768000000000, 0.748033000000000, 0.616944000000000; ...
        0.895867000000000, 0.735027000000000, 0.605643000000000; ...
        0.892966000000000, 0.722021000000000, 0.594342000000000; ...
        0.890065000000000, 0.709015000000000, 0.583040000000000; ...
        0.887164000000000, 0.696009000000000, 0.571739000000000; ...
        0.884263000000000, 0.683002000000000, 0.560438000000000; ...
        0.881362000000000, 0.669996000000000, 0.549136000000000; ...
        0.878461000000000, 0.656990000000000, 0.537835000000000; ...
        0.875560000000000, 0.643984000000000, 0.526534000000000; ...
        0.870352000000000, 0.628610000000000, 0.515234000000000; ...
        0.863825000000000, 0.611884000000000, 0.503935000000000; ...
        0.857298000000000, 0.595158000000000, 0.492636000000000; ...
        0.850771000000000, 0.578431000000000, 0.481336000000000; ...
        0.844244000000000, 0.561705000000000, 0.470037000000000; ...
        0.837718000000000, 0.544979000000000, 0.458738000000000; ...
        0.831191000000000, 0.528253000000000, 0.447439000000000; ...
        0.824664000000000, 0.511526000000000, 0.436140000000000; ...
        0.818137000000000, 0.494800000000000, 0.424841000000000; ...
        0.811610000000000, 0.478074000000000, 0.413542000000000; ...
        0.805083000000000, 0.461348000000000, 0.402243000000000; ...
        0.798257000000000, 0.444538000000000, 0.390997000000000; ...
        0.788430000000000, 0.426895000000000, 0.380282000000000; ...
        0.778603000000000, 0.409253000000000, 0.369567000000000; ...
        0.768776000000000, 0.391610000000000, 0.358852000000000; ...
        0.758949000000000, 0.373967000000000, 0.348137000000000; ...
        0.749123000000000, 0.356325000000000, 0.337422000000000; ...
        0.739296000000000, 0.338682000000000, 0.326707000000000; ...
        0.729469000000000, 0.321039000000000, 0.315992000000000; ...
        0.719642000000000, 0.303397000000000, 0.305277000000000; ...
        0.709815000000000, 0.285754000000000, 0.294562000000000; ...
        0.699989000000000, 0.268111000000000, 0.283847000000000; ...
        0.690162000000000, 0.250468000000000, 0.273131000000000; ...
        0.678532000000000, 0.235025000000000, 0.263462000000000; ...
        0.665400000000000, 0.221415000000000, 0.254665000000000; ...
        0.652268000000000, 0.207805000000000, 0.245867000000000; ...
        0.639136000000000, 0.194195000000000, 0.237070000000000; ...
        0.626004000000000, 0.180585000000000, 0.228272000000000; ...
        0.612873000000000, 0.166974000000000, 0.219475000000000; ...
        0.599741000000000, 0.153364000000000, 0.210677000000000; ...
        0.586609000000000, 0.139754000000000, 0.201880000000000; ...
        0.573477000000000, 0.126144000000000, 0.193082000000000; ...
        0.560345000000000, 0.112534000000000, 0.184285000000000; ...
        0.547213000000000, 0.0989234000000000, 0.175487000000000; ...
        0.534081000000000, 0.0853132000000000, 0.166690000000000];
end
