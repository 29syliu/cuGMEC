%% cuGMEC 相空间轨道频率可视化脚本

%{

脚本流程：

(1) 用户设置
设置输入目录, normalization2D.mat 路径, Orbit.bin 所在目录, 需要处理的物种，
以及是否打开守恒量诊断和旧的 dT 切片图。

(2) 读取 normalization2D.mat
读取相空间网格, 各物种 EPphiLambda 范围, 以及后续频率换算需要的归一化常数。

(3) 逐物种处理 Orbit.bin
读取 Ion/Alpha/BeamPhaseSpaceOrbit.bin, 剔除 pad 记录, 按 trapped, para,
anti 的轨道判据分类，并打印最终分类占比。

(4) 绘制相空间轨道频率
在本节选择绘图物种, 固定坐标, 切片, 频率单位和分支方向。

%}

%% (1) 用户设置

inputPath = 'C:\Users\Desktop\ITER\';
normalizationFile = fullfile(inputPath, 'normalization2D.mat');

orbitPath = inputPath;

speciesList = {'Ion', 'Alpha', 'Beam'};

plotInvariantDiagnostics = false;
plotDTSlices = false;

sliceDim = 1;       % 1: 固定 E, 2: 固定 Pphi, 3: 固定 Lambda
sliceIndex = [];    % 留空时自动取该方向中间切片；也可手动设为例如 50

%% (2) 读取 normalization2D.mat

normData = load(normalizationFile);
validateNormalizationMetadata(normData);

gridE = readPositiveIntegerScalar(normData.gridE, 'gridE');
gridPphi = readPositiveIntegerScalar(normData.gridPphi, 'gridPphi');
gridLambda = readPositiveIntegerScalar(normData.gridLambda, 'gridLambda');

if isempty(sliceIndex)
    sliceIndex = defaultSliceIndex(sliceDim, gridE, gridPphi, gridLambda);
end

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

    if plotInvariantDiagnostics
        plotConservationDiagnostics(speciesName, result.diagnostics, ...
            result.E1d, result.Pphi1d, result.Lambda1d);
    end

    if plotDTSlices
        plotSpeciesDTSlices(speciesName, result.phaseSpaceData, ...
            result.E1d, result.Pphi1d, result.Lambda1d, sliceDim, sliceIndex);
    end
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

plotPhaseSpaceOrbitFrequency(phaseSpaceResults, normData, ...
    orbitFrequencySpecies, orbitFrequencyFixedCoordinate, orbitFrequencySlice, ...
    orbitFrequencyUnit, orbitFrequencyBranch);

%% local functions

function result = analyzeSpeciesOrbit(speciesName, orbitFile, phaseRange, gridE, gridPphi, gridLambda)

    data = readOrbitBinary(orbitFile);

    Ids = int32(data(:, 1));
    orbits = data(:, 2);
    dtheta = data(:, 3);
    dphiTotal = data(:, 4);
    dphiVpara = data(:, 5);
    dTs = data(:, 6);
    Es = data(:, 7);
    Pphis = data(:, 8);
    Lambdas = data(:, 9);

    assertNoNaN(speciesName, Ids, orbits, dtheta, dphiTotal, dphiVpara, dTs, Es, Pphis, Lambdas);

    nPhase = gridE * gridPphi * gridLambda;
    expectedRecords = 2 * nPhase;
    numRecords = size(data, 1);
    assert(numRecords == expectedRecords, ...
        '%s 记录数为 %d，但 gridE*gridPphi*gridLambda*2 = %d。', ...
        orbitFile, numRecords, expectedRecords);

    E1d = linspace(phaseRange(1), phaseRange(2), gridE);
    Pphi1d = linspace(phaseRange(3), phaseRange(4), gridPphi);
    Lambda1d = linspace(phaseRange(5), phaseRange(6), gridLambda);

    phaseSpaceData = emptyPhaseSpaceData(gridE, gridPphi, gridLambda, E1d, Pphi1d, Lambda1d, phaseRange);

    recordBranch = ones(numRecords, 1);
    recordBranch(nPhase + 1:end) = -1;

    validRecord = ~isPadRecord(Ids);
    localId = abs(double(Ids));
    validLocalId = validRecord & localId >= 0 & localId < nPhase & localId == floor(localId);

    if any(validRecord & ~validLocalId)
        badIndex = find(validRecord & ~validLocalId, 1);
        error('%s 中存在非法相空间 id：record=%d, id=%d。', orbitFile, badIndex, Ids(badIndex));
    end

    warnIfSignedIdOrderLooksWrong(speciesName, Ids, recordBranch, validRecord);

    [plusIndex, minusIndex, duplicateCounts] = buildBranchIndex(localId, recordBranch, validLocalId, nPhase);

    diagnostics = struct();
    summary = initializeSummary(speciesName, orbitFile, numRecords, nPhase, validRecord, orbits);
    summary.duplicatePlusIds = duplicateCounts.plus;
    summary.duplicateMinusIds = duplicateCounts.minus;

    [phaseSpaceData.trapped, diagnostics.trapped, trappedSummary] = extractTrappedOrbits( ...
        phaseSpaceData.trapped, plusIndex, minusIndex, orbits, dtheta, dphiTotal, dphiVpara, dTs, ...
        Es, Pphis, Lambdas, gridE, gridPphi, gridLambda);
    summary.trapped = trappedSummary;

    [phaseSpaceData.para, diagnostics.para, paraSummary] = extractPassingOrbits( ...
        2.5, phaseSpaceData.para, validLocalId, localId, recordBranch, plusIndex, minusIndex, ...
        orbits, dtheta, dphiTotal, dphiVpara, dTs, Es, Pphis, Lambdas, ...
        gridE, gridPphi, gridLambda);
    summary.para = paraSummary;

    [phaseSpaceData.anti, diagnostics.anti, antiSummary] = extractPassingOrbits( ...
        3.5, phaseSpaceData.anti, validLocalId, localId, recordBranch, plusIndex, minusIndex, ...
        orbits, dtheta, dphiTotal, dphiVpara, dTs, Es, Pphis, Lambdas, ...
        gridE, gridPphi, gridLambda);
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

function [classData, diagnostic, summary] = extractTrappedOrbits(classData, plusIndex, minusIndex, ...
    orbits, dtheta, dphiTotal, dphiVpara, dTs, Es, Pphis, Lambdas, ...
    gridE, gridPphi, gridLambda)

    pairedLocalIndex = find(plusIndex > 0 & minusIndex > 0);
    plusRecord = plusIndex(pairedLocalIndex);
    minusRecord = minusIndex(pairedLocalIndex);

    trappedOrbit = isOrbit(orbits(plusRecord), 4.5) & isOrbit(orbits(minusRecord), 4.5);
    dTRelDiff = relativePairDifference(dTs(plusRecord), dTs(minusRecord));
    accepted = trappedOrbit & dTRelDiff < trappedOrbitRelativeTolerance();

    acceptedLocalId = pairedLocalIndex(accepted) - 1;
    acceptedPlusRecord = plusRecord(accepted);
    acceptedMinusRecord = minusRecord(accepted);

    classData = fillOrbitClass(classData, acceptedLocalId, ...
        averagePair(dtheta(acceptedPlusRecord), dtheta(acceptedMinusRecord)), ...
        averagePair(dphiTotal(acceptedPlusRecord), dphiTotal(acceptedMinusRecord)), ...
        averagePair(dphiVpara(acceptedPlusRecord), dphiVpara(acceptedMinusRecord)), ...
        averagePair(dTs(acceptedPlusRecord), dTs(acceptedMinusRecord)), ...
        gridE, gridPphi, gridLambda);

    diagnostic = struct( ...
        'localId', acceptedLocalId(:), ...
        'E', averagePair(Es(acceptedPlusRecord), Es(acceptedMinusRecord)), ...
        'Pphi', averagePair(Pphis(acceptedPlusRecord), Pphis(acceptedMinusRecord)), ...
        'Lambda', averagePair(Lambdas(acceptedPlusRecord), Lambdas(acceptedMinusRecord)), ...
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
    recordBranch, plusIndex, minusIndex, orbits, dtheta, dphiTotal, dphiVpara, dTs, Es, Pphis, Lambdas, ...
    gridE, gridPphi, gridLambda)

    targetRecord = find(validRecord & isOrbit(orbits, targetOrbit));
    targetLocalId = localId(targetRecord);
    targetBranch = recordBranch(targetRecord);

    counterpartRecord = zeros(size(targetRecord));
    plusTarget = targetBranch > 0;
    minusTarget = targetBranch < 0;

    counterpartRecord(plusTarget) = minusIndex(targetLocalId(plusTarget) + 1);
    counterpartRecord(minusTarget) = plusIndex(targetLocalId(minusTarget) + 1);

    hasCounterpart = counterpartRecord > 0;
    counterpartIsTrapped = false(size(targetRecord));
    counterpartIsTrapped(hasCounterpart) = isOrbit(orbits(counterpartRecord(hasCounterpart)), 4.5);

    keep = ~counterpartIsTrapped;
    keptRecord = targetRecord(keep);
    keptLocalId = targetLocalId(keep);
    keptHasCounterpart = hasCounterpart(keep);

    [uniqueLocalId, uniquePosition] = unique(keptLocalId, 'stable');
    selectedRecord = keptRecord(uniquePosition);
    selectedHasCounterpart = keptHasCounterpart(uniquePosition);

    classData = fillOrbitClass(classData, uniqueLocalId, ...
        dtheta(selectedRecord), dphiTotal(selectedRecord), dphiVpara(selectedRecord), dTs(selectedRecord), ...
        gridE, gridPphi, gridLambda);

    diagnostic = struct( ...
        'localId', uniqueLocalId(:), ...
        'E', Es(selectedRecord), ...
        'Pphi', Pphis(selectedRecord), ...
        'Lambda', Lambdas(selectedRecord), ...
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

function plotPhaseSpaceOrbitFrequency(phaseSpaceResults, normData, speciesName, fixedCoordinate, sliceText, unitText, branchText)

    [phaseSpaceData, speciesLabel] = resolveSpeciesPhaseSpaceData(phaseSpaceResults, speciesName);
    dim = phaseCoordinateToDimension(fixedCoordinate);
    idx = parseOrbitFrequencySlice(sliceText, dim, phaseSpaceData);
    branchName = normalizeFrequencyBranch(branchText);
    [unitScale, colorbarLabel] = frequencyUnitScale(unitText);

    L0 = readPositiveScalar(normData, 'L0');
    VA0 = readPositiveScalar(normData, 'VA0');

    dT = phaseSpaceData.(branchName).dT;
    trappedDT = phaseSpaceData.trapped.dT;
    fillFromTrapped = (dT == 0 | ~isfinite(dT)) & isfinite(trappedDT) & trappedDT > 0;
    dT(fillFromTrapped) = trappedDT(fillFromTrapped);
    dT(~isfinite(dT) | dT <= 0) = NaN;

    frequency = unitScale ./ (dT .* (L0 / VA0));
    frequency(~isfinite(frequency)) = NaN;

    [Z, xVec, yVec, xlabelText, ylabelText, titleText] = sliceField( ...
        dim, idx, frequency, phaseSpaceData.E1d, phaseSpaceData.Pphi1d, phaseSpaceData.Lambda1d);

    if all(isnan(Z(:)))
        fprintf('[plot] %s %s has no orbit-frequency data in this slice.\n', speciesLabel, branchName);
        return;
    end

    Z = fillIsolatedNaN(Z);

    figure('Name', sprintf('%s %s orbit frequency', speciesLabel, branchName), ...
        'Position', [100, 100, 900, 760]);
    [X, Y] = meshgrid(xVec, yVec);
    pcolor(X, Y, Z);
    shading interp;
    hold on;

    finiteZ = Z(isfinite(Z));
    if numel(unique(finiteZ)) > 1
        contour(X, Y, Z, 10, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end

    colormap(jet);
    cb = colorbar;
    cb.FontName = 'Times New Roman';
    cb.FontSize = 14;
    ylabel(cb, colorbarLabel, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    xlabel(xlabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(ylabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(titleText, ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    axis tight;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Color', 'white', 'Layer', 'top');
    hold off;
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

function branchName = normalizeFrequencyBranch(branchText)

    branchName = lower(strtrim(char(branchText)));
    if ~ismember(branchName, {'para', 'anti'})
        error('frequency branch must be "para" or "anti".');
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

function plotConservationDiagnostics(speciesName, diagnostics, E1d, Pphi1d, Lambda1d)

    classNames = {'trapped', 'para', 'anti'};
    classLabels = {'trapped', 'para', 'anti'};

    for classIndex = 1:numel(classNames)
        className = classNames{classIndex};
        if ~isfield(diagnostics, className) || isempty(diagnostics.(className).localId)
            fprintf('[plot] %s %s 无守恒量诊断点。\n', speciesName, classLabels{classIndex});
            continue;
        end

        plotInvariantErrors(speciesName, classLabels{classIndex}, diagnostics.(className), E1d, Pphi1d, Lambda1d);
    end
end

function plotInvariantErrors(speciesName, classLabel, diagnostic, E1d, Pphi1d, Lambda1d)

    [baseE, basePphi, baseLambda] = initialCoordinatesFromLocalId(diagnostic.localId, E1d, Pphi1d, Lambda1d);

    errorE = mixedConservationError(diagnostic.E, baseE);
    errorPphi = mixedConservationError(diagnostic.Pphi, basePphi);
    errorLambda = mixedConservationError(diagnostic.Lambda, baseLambda);

    figure('Name', [speciesName ' ' classLabel ' 守恒量误差'], 'Position', [100, 100, 900, 760]);

    subplot(3, 1, 1);
    plot(errorE, '.');
    grid on;
    ylabel('E', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(sprintf('%s %s relative error', speciesName, classLabel), ...
        'Interpreter', 'none', 'FontName', 'Times New Roman', 'FontSize', 14);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

    subplot(3, 1, 2);
    plot(errorPphi, '.');
    grid on;
    ylabel('Pphi', 'FontName', 'Times New Roman', 'FontSize', 14);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

    subplot(3, 1, 3);
    plot(errorLambda, '.');
    grid on;
    ylabel('Lambda', 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('particle index', 'FontName', 'Times New Roman', 'FontSize', 14);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

function errorValue = mixedConservationError(finalValue, baselineValue)

    finalValue = finalValue(:);
    baselineValue = baselineValue(:);

    errorValue = finalValue - baselineValue;
    relativeMask = abs(baselineValue) > 1e-12;
    errorValue(relativeMask) = errorValue(relativeMask) ./ baselineValue(relativeMask);
end

function plotSpeciesDTSlices(speciesName, phaseSpaceData, E1d, Pphi1d, Lambda1d, sliceDim, sliceIndex)

    plotDTSlice([speciesName ' trapped'], sliceDim, sliceIndex, phaseSpaceData.trapped.dT, E1d, Pphi1d, Lambda1d);
    plotDTSlice([speciesName ' para'], sliceDim, sliceIndex, phaseSpaceData.para.dT, E1d, Pphi1d, Lambda1d);
    plotDTSlice([speciesName ' anti'], sliceDim, sliceIndex, phaseSpaceData.anti.dT, E1d, Pphi1d, Lambda1d);
end

function plotDTSlice(titlePrefix, dim, idx, dT, E1d, Pphi1d, Lambda1d)

    [Z, xVec, yVec, xlabelText, ylabelText, titleText] = sliceField(dim, idx, dT, E1d, Pphi1d, Lambda1d);
    Z(Z == 0) = NaN;

    if all(isnan(Z(:)))
        fprintf('[plot] %s 在当前切片没有 dT 数据。\n', titlePrefix);
        return;
    end

    Z = fillIsolatedNaN(Z);

    figure('Name', [titlePrefix ' dT']);
    [X, Y] = meshgrid(xVec, yVec);
    pcolor(X, Y, Z);
    shading interp;
    hold on;
    contour(X, Y, Z, 10, 'EdgeColor', 'k', 'LineWidth', 0.5);
    colormap(jet);
    cb = colorbar;
    cb.FontName = 'Times New Roman';
    cb.FontSize = 12;
    ylabel(cb, 'dT');

    xlabel(xlabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(ylabelText, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title([titlePrefix ', ' titleText], 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'white');
end

function [Z, xVec, yVec, xlabelText, ylabelText, titleText] = sliceField(dim, idx, field, E1d, Pphi1d, Lambda1d)

    switch dim
        case 1
            idx = clampIndex(idx, numel(E1d), 'E');
            Z = squeeze(field(idx, :, :));
            xVec = Lambda1d;
            yVec = -Pphi1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$P_{\varphi}$';
            titleText = sprintf('$E = %.6g$', E1d(idx));
        case 2
            idx = clampIndex(idx, numel(Pphi1d), 'Pphi');
            Z = squeeze(field(:, idx, :));
            xVec = Lambda1d;
            yVec = E1d;
            xlabelText = '$\Lambda$';
            ylabelText = '$E$';
            titleText = sprintf('$P_{\\varphi} = %.6g$', -Pphi1d(idx));
        case 3
            idx = clampIndex(idx, numel(Lambda1d), 'Lambda');
            Z = squeeze(field(:, :, idx));
            xVec = -Pphi1d;
            yVec = E1d;
            xlabelText = '$P_{\varphi}$';
            ylabelText = '$E$';
            titleText = sprintf('$\\Lambda = %.6g$', Lambda1d(idx));
        otherwise
            error('sliceDim 必须是 1, 2 或 3。');
    end
end

function Z = fillIsolatedNaN(Z)

    [ny, nx] = size(Z);
    if ny < 3 || nx < 3
        return;
    end

    offsets = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 1; 1, -1; 1, 0; 1, 1];
    dist = sqrt(sum(offsets.^2, 2));
    weights = 1 ./ dist;
    weights = weights / sum(weights);

    for i = 2:ny - 1
        for j = 2:nx - 1
            if isnan(Z(i, j))
                neighborValues = zeros(1, 8);
                valid = true;
                for k = 1:8
                    ni = i + offsets(k, 1);
                    nj = j + offsets(k, 2);
                    val = Z(ni, nj);
                    if isnan(val)
                        valid = false;
                        break;
                    end
                    neighborValues(k) = val;
                end
                if valid
                    Z(i, j) = neighborValues * weights;
                end
            end
        end
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

function assertNoNaN(speciesName, Ids, orbits, dtheta, dphiTotal, dphiVpara, dTs, Es, Pphis, Lambdas)

    nanCount = sum(isnan(double(Ids))) + sum(isnan(orbits)) + sum(isnan(dtheta)) + ...
        sum(isnan(dphiTotal)) + sum(isnan(dphiVpara)) + sum(isnan(dTs)) + ...
        sum(isnan(Es)) + sum(isnan(Pphis)) + sum(isnan(Lambdas));

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

function idx = defaultSliceIndex(sliceDim, gridE, gridPphi, gridLambda)

    switch sliceDim
        case 1
            idx = max(1, round(gridE / 2));
        case 2
            idx = max(1, round(gridPphi / 2));
        case 3
            idx = max(1, round(gridLambda / 2));
        otherwise
            error('sliceDim 必须是 1, 2 或 3。');
    end
end

function idx = clampIndex(idx, maxIndex, label)

    assert(isscalar(idx) && isfinite(idx) && idx == floor(idx), '%s 切片下标必须是整数。', label);
    if idx < 1 || idx > maxIndex
        error('%s 切片下标越界：idx=%d，有效范围为 [1, %d]。', label, idx, maxIndex);
    end
end
