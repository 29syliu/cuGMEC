%% cuGMEC 相空间轨道映射生成脚本

%{

1. 脚本功能
读取 normalization2D.mat，生成 cuGMEC 相空间轨道映射文件。

2. 代码块 (1) 的参数含义
inputPath：normalization2D.mat 所在目录。当前脚本只需要该目录中已有 normalization2D.mat。
outputPath：相空间映射二进制文件的输出目录。默认与 inputPath 相同。
plotGeometry：是否画出平衡边界和当前 thetaIndex 对应的径向线。

3. 代码块 (2) 的参数含义
gridE，gridPphi，gridLambda 必须与 cuGMEC 参数文件中的相空间网格数量一致。
speciesList 用于设置需要生成映射文件的物种。
name 只能是 Ion、Alpha 或 Beam。
Mass，Char，Vmin，Vmax 必须与 cuGMEC 中对应物种参数一致。

4. 运行代码块 (3)，outputPath 下会生成：
IonPhaseSpaceMapping.bin，AlphaPhaseSpaceMapping.bin 或 BeamPhaseSpaceMapping.bin。

%}

%% (1)

inputPath = 'C:\Users\Desktop\ITER\';
outputPath = inputPath;

plotGeometry = true;
phaseSpaceFile = fullfile(inputPath, 'normalization2D.mat');

%% (2)

gridE = 64;
gridPphi = 64;
gridLambda = 64;

% 增加数组元素即可在一次运行中生成多个物种。
% name 只能是 Ion、Alpha 或 Beam。
speciesList = struct( ...
    'name', {'Alpha'}, ...
    'Mass', {4.0}, ...
    'Char', {2.0}, ...
    'Vmin', {0.0733}, ...
    'Vmax', {1.466});

%% (3)


phaseSpace = load(phaseSpaceFile);
validatePhaseSpaceData(phaseSpace);
phaseSpaceRanges = initializePhaseSpaceRangeMetadata(phaseSpace);

thetaIndex = size(phaseSpace.B, 2) / 2;

if plotGeometry
    figure;
    hold on;
    plot(phaseSpace.R(1, :), phaseSpace.Z(1, :), 'r');
    axis equal;
    plot(phaseSpace.R(end, :), phaseSpace.Z(end, :), 'r');
    axis equal;
    plot(phaseSpace.R(:, thetaIndex), phaseSpace.Z(:, thetaIndex), 'r');
    axis equal;
end

for speciesIndex = 1:numel(speciesList)
    summary = writePhaseSpaceMappingForSpecies(phaseSpace, speciesList(speciesIndex), ...
        thetaIndex, gridE, gridPphi, gridLambda, outputPath);
    phaseSpaceRanges.([summary.species 'EPphiLambda']) = summary.EPphiLambda;
end

savePhaseSpaceMappingMetadata(phaseSpaceFile, gridE, gridPphi, gridLambda, phaseSpaceRanges);

%%

function summary = writePhaseSpaceMappingForSpecies(eq, species, thetaIndex, sizeE, sizePphi, sizeLambda, outputPath)

    validateSpecies(species);
    padId = phaseSpacePadId();

    fprintf('\n============================================================\n');
    fprintf('物种：%s，相空间映射生成\n', species.name);
    fprintf('============================================================\n');

    [paraRho, paraVpara, paraMu, paraId, paraError, paraNonemptyRatio, ranges] = ...
        computePhaseSpaceMapping(eq, species, 1, thetaIndex, sizeE, sizePphi, sizeLambda, padId);

    [antiRho, antiVpara, antiMu, antiId, antiError, antiNonemptyRatio] = ...
        computePhaseSpaceMapping(eq, species, -1, thetaIndex, sizeE, sizePphi, sizeLambda, padId);

    % 用 id 符号区分反向分支，同时保留 pad 标记不变。
    antiMask = antiId ~= padId;
    antiId(antiMask) = -antiId(antiMask);

    allRho = [paraRho; antiRho];
    allVpara = [paraVpara; antiVpara];
    allMu = [paraMu; antiMu];
    allId = [paraId; antiId];
    allError = [paraError; antiError];

    validateMappingArrays(allRho, allVpara, allMu, allId, padId);

    fprintf('\n--- 相空间范围 ---\n');
    fprintf('E      : [%g, %g]\n', ranges.minE, ranges.maxE);
    fprintf('Pphi   : [%g, %g]\n', ranges.minPphi, ranges.maxPphi);
    fprintf('Lambda : [%g, %g]\n', ranges.minLambda, ranges.maxLambda);

    fprintf('\n--- 映射有效率 ---\n');
    fprintf('正向分支：%.2f%%\n', 100 * paraNonemptyRatio);
    fprintf('反向分支：%.2f%%\n', 100 * antiNonemptyRatio);
    fprintf('总体映射：%.2f%%\n', 100 * sum(allId ~= padId) / numel(allId));

    validErrors = allError(allId ~= padId);
    if ~isempty(validErrors)
        fprintf('\n--- Pphi 反解相对误差 ---\n');
        fprintf('最大值：%.3e\n', max(validErrors));
        fprintf('平均值：%.3e\n', mean(validErrors));
    end

    outputFile = fullfile(outputPath, [species.name 'PhaseSpaceMapping.bin']);
    writeMappingBinary(outputFile, allId, allRho, allVpara, allMu);

    fprintf('============================================================\n\n');

    summary = struct( ...
        'species', species.name, ...
        'outputFile', outputFile, ...
        'records', numel(allId), ...
        'validRecords', sum(allId ~= padId), ...
        'EPphiLambda', [ranges.minE, ranges.maxE, ranges.minPphi, ranges.maxPphi, ...
            ranges.minLambda, ranges.maxLambda]);
end

function phaseSpaceRanges = initializePhaseSpaceRangeMetadata(phaseSpace)

    speciesNames = {'Ion', 'Alpha', 'Beam'};
    phaseSpaceRanges = struct();

    for speciesIndex = 1:numel(speciesNames)
        fieldName = [speciesNames{speciesIndex} 'EPphiLambda'];
        if isfield(phaseSpace, fieldName) && numel(phaseSpace.(fieldName)) == 6
            phaseSpaceRanges.(fieldName) = reshape(double(phaseSpace.(fieldName)), 1, 6);
        else
            phaseSpaceRanges.(fieldName) = nan(1, 6);
        end
    end
end

function savePhaseSpaceMappingMetadata(phaseSpaceFile, gridE, gridPphi, gridLambda, phaseSpaceRanges)

    IonEPphiLambda = phaseSpaceRanges.IonEPphiLambda;
    AlphaEPphiLambda = phaseSpaceRanges.AlphaEPphiLambda;
    BeamEPphiLambda = phaseSpaceRanges.BeamEPphiLambda;

    save(phaseSpaceFile, 'gridE', 'gridPphi', 'gridLambda', ...
        'IonEPphiLambda', 'AlphaEPphiLambda', 'BeamEPphiLambda', '-append');

    fprintf('\n--- normalization2D.mat 相空间元数据 ---\n');
    fprintf('gridE/gridPphi/gridLambda: %d / %d / %d\n', gridE, gridPphi, gridLambda);
    fprintf('IonEPphiLambda   : [%s]\n', num2str(IonEPphiLambda, ' %.16g'));
    fprintf('AlphaEPphiLambda : [%s]\n', num2str(AlphaEPphiLambda, ' %.16g'));
    fprintf('BeamEPphiLambda  : [%s]\n', num2str(BeamEPphiLambda, ' %.16g'));
end

function padId = phaseSpacePadId()

    padId = int32(20251106);
end

function validateSpecies(species)

    validSpecies = {'Ion', 'Alpha', 'Beam'};

    assert(isfield(species, 'name') && any(strcmp(species.name, validSpecies)), ...
        'species.name 必须是 Ion、Alpha 或 Beam。');
    assert(isfield(species, 'Mass') && isscalar(species.Mass) && species.Mass > 0, ...
        'species.Mass 必须为正数。');
    assert(isfield(species, 'Char') && isscalar(species.Char) && species.Char ~= 0, ...
        'species.Char 不能为 0。');
    assert(isfield(species, 'Vmin') && isscalar(species.Vmin) && species.Vmin >= 0, ...
        'species.Vmin 必须为非负数。');
    assert(isfield(species, 'Vmax') && isscalar(species.Vmax) && species.Vmax > species.Vmin, ...
        'species.Vmax 必须大于 species.Vmin。');
end

function validatePhaseSpaceData(phaseSpace)

    requiredFields = {'MP', 'QE', 'B', 'J', 'psip', 'SFAcovyz', 'B0', 'L0', 'VA0', ...
        'RHO0', 'RHO1', 'PSITMAX', 'q', 'R', 'Z', 'rho', 'theta'};

    for fieldIndex = 1:numel(requiredFields)
        assert(isfield(phaseSpace, requiredFields{fieldIndex}), ...
            'normalization2D.mat 缺少字段 "%s"。', requiredFields{fieldIndex});
    end
end

function validateMappingArrays(rhoArr, vparaArr, muArr, idArr, padId)

    assert(numel(rhoArr) == numel(vparaArr) && numel(rhoArr) == numel(muArr) && numel(rhoArr) == numel(idArr), ...
        '映射数组长度不一致。');
    assert(all(isfinite(rhoArr)) && all(isfinite(vparaArr)) && all(isfinite(muArr)), ...
        '映射数组包含 NaN 或 Inf。');
    assert(all(rhoArr >= 0 & rhoArr < 1), '映射后的 rho 必须位于 [0, 1) 内。');
    assert(any(idArr ~= padId), '映射结果中没有有效相空间点。');
end

function writeMappingBinary(outputFile, idArr, rhoArr, vparaArr, muArr)

    fid = fopen(outputFile, 'wb');
    if fid == -1
        error('无法创建文件：%s', outputFile);
    end
    closeFile = onCleanup(@() fclose(fid));

    assert(fwrite(fid, idArr, 'int32') == numel(idArr), '写入 id 数组失败。');
    assert(fwrite(fid, rhoArr, 'double') == numel(rhoArr), '写入 rho 数组失败。');
    assert(fwrite(fid, vparaArr, 'double') == numel(vparaArr), '写入 vpara 数组失败。');
    assert(fwrite(fid, muArr, 'double') == numel(muArr), '写入 mu 数组失败。');

    clear closeFile;
    fileInfo = dir(outputFile);
    fprintf('\n--- 输出文件 ---\n');
    fprintf('文件：%s\n', outputFile);
    fprintf('记录数：%d\n', numel(idArr));
    fprintf('大小：%.3f MB\n', fileInfo.bytes / 1024^2);
end

function [rhoArr, vparaArr, muArr, idArr, PphiErrorArr, nonemptyRatio, ranges] = ...
    computePhaseSpaceMapping(eq, species, branchSign, thetaIndex, sizeE, sizePphi, sizeLambda, padId)

    B = eq.B;
    J = eq.J;
    psip = eq.psip;
    SFAcovyz = eq.SFAcovyz;
    q = eq.q;
    rho = eq.rho;

    drho = eq.RHO1 - eq.RHO0;
    psitmax = eq.PSITMAX / (eq.B0 * eq.L0 * eq.L0);
    cm = eq.VA0 / (eq.L0 * (eq.QE * eq.B0 / eq.MP));

    minE = 0.5 * species.Mass * species.Vmin^2;
    maxE = 0.5 * species.Mass * species.Vmax^2;

    Pphi0 = -cm * species.Mass * species.Vmax * 2 * psitmax * drho .* rho .* SFAcovyz ./ ...
        (q .* J .* B) - species.Char .* psip;
    Pphi1 = cm * species.Mass * species.Vmax * 2 * psitmax * drho .* rho .* SFAcovyz ./ ...
        (q .* J .* B) - species.Char .* psip;

    minPphi = min(Pphi0(:));
    maxPphi = max(Pphi1(:));

    minLambda = 0.0;
    maxLambda = max(max(1 ./ B));

    gridE = linspace(minE, maxE, sizeE);
    gridPphi = linspace(minPphi, maxPphi, sizePphi);
    gridLambda = linspace(minLambda, maxLambda, sizeLambda);

    ranges = struct( ...
        'minE', minE, ...
        'maxE', maxE, ...
        'minPphi', minPphi, ...
        'maxPphi', maxPphi, ...
        'minLambda', minLambda, ...
        'maxLambda', maxLambda);

    rhoRadial = rho(:, thetaIndex);
    BRadial = B(:, thetaIndex);
    SFAcovyzRadial = SFAcovyz(:, thetaIndex);
    qRadial = q(:, thetaIndex);
    JRadial = J(:, thetaIndex);
    psipRadial = psip(:, thetaIndex);

    countValid0 = 0;
    countValid1 = 0;
    minValidGt1 = inf;
    nonmonotonicCount = 0;
    singleExtremaCount = 0;
    multipleExtremaCount = 0;
    equalAdjacentCount = 0;

    validRhoPphi = cell(numel(gridE), numel(gridLambda));

    for i = 1:numel(gridE)
        E = gridE(i);
        for k = 1:numel(gridLambda)
            Lambda = gridLambda(k);
            sqrtArg = 1 - Lambda * BRadial;
            valid = sqrtArg >= 0;
            validCount = sum(valid);

            if validCount == 0
                countValid0 = countValid0 + 1;
                continue;
            elseif validCount == 1
                countValid1 = countValid1 + 1;
                continue;
            elseif validCount < minValidGt1
                minValidGt1 = validCount;
            end

            rhoValid = rhoRadial(valid);
            BValid = BRadial(valid);
            SFAcovyzValid = SFAcovyzRadial(valid);
            qValid = qRadial(valid);
            JValid = JRadial(valid);
            psipValid = psipRadial(valid);

            Vpara = branchSign * sqrt(2 * E * sqrtArg(valid) / species.Mass);
            Pphi = cm * species.Mass * Vpara .* 2 * psitmax * drho .* rhoValid .* SFAcovyzValid ./ ...
                (qValid .* JValid .* BValid) - species.Char .* psipValid;

            diffPphi = diff(Pphi);
            equalAdjacentCount = equalAdjacentCount + any(abs(diffPphi) < 1e-12);

            isStrictMonotonic = all(diffPphi > 0) || all(diffPphi < 0);
            isConstant = all(diffPphi == 0);

            needStore = false;
            storedRho = [];
            storedPphi = [];

            if isStrictMonotonic && ~isConstant
                needStore = true;
                storedRho = rhoValid;
                storedPphi = Pphi;
            elseif ~isStrictMonotonic && ~isConstant
                nonmonotonicCount = nonmonotonicCount + 1;

                signSeq = sign(diffPphi);
                filteredSign = signSeq(signSeq ~= 0);
                if isempty(filteredSign)
                    numExtrema = 0;
                else
                    numExtrema = sum(abs(diff(filteredSign)) > 0);
                end

                if numExtrema == 1
                    singleExtremaCount = singleExtremaCount + 1;
                    [~, idxMax] = max(Pphi);
                    [~, idxMin] = min(Pphi);
                    if idxMax > 1 && idxMax < numel(Pphi)
                        extremumIdx = idxMax;
                    elseif idxMin > 1 && idxMin < numel(Pphi)
                        extremumIdx = idxMin;
                    else
                        continue;
                    end

                    rhoLeft = rhoValid(extremumIdx) - rhoValid(1);
                    rhoRight = rhoValid(end) - rhoValid(extremumIdx);
                    if rhoLeft < rhoRight
                        storedRho = rhoValid(extremumIdx:end);
                        storedPphi = Pphi(extremumIdx:end);
                    else
                        storedRho = rhoValid(1:extremumIdx);
                        storedPphi = Pphi(1:extremumIdx);
                    end
                    needStore = true;
                elseif numExtrema >= 2
                    multipleExtremaCount = multipleExtremaCount + 1;
                end
            end

            if needStore
                validRhoPphi{i, k} = struct('rho', storedRho, 'Pphi', storedPphi);
            end
        end
    end

    if branchSign > 0
        branchName = '正向分支 vpara > 0';
    else
        branchName = '反向分支 vpara < 0';
    end

    fprintf('\n--- %s ---\n', branchName);
    fprintf('有效 rho 点数统计：无有效点=%d，仅 1 个有效点=%d，多于 1 个有效点时最少点数=%g\n', ...
        countValid0, countValid1, minValidGt1);
    fprintf('Pphi 单调性诊断：相邻近等值=%d，非单调=%d，单极值=%d，多极值=%d\n', ...
        equalAdjacentCount, nonmonotonicCount, singleExtremaCount, multipleExtremaCount);

    M = sizePphi * sizeLambda;
    N = sizeE * M;

    rhoCell = cell(sizeE, 1);
    vparaCell = cell(sizeE, 1);
    muCell = cell(sizeE, 1);
    idCell = cell(sizeE, 1);
    errorCell = cell(sizeE, 1);

    parfor i = 1:sizeE
        rhoLocal = 0.5 * ones(M, 1);
        vparaLocal = zeros(M, 1);
        muLocal = zeros(M, 1);
        idLocal = repmat(padId, M, 1);
        errorLocal = zeros(M, 1);

        Ei = gridE(i);
        for k = 1:sizeLambda
            if isempty(validRhoPphi{i, k})
                continue;
            end

            PphiVec = validRhoPphi{i, k}.Pphi;
            rhoVec = validRhoPphi{i, k}.rho;

            % pchip 需要有序网格；排序也能统一处理 Pphi 递减分支。
            [PphiInterp, sortIdx] = sort(PphiVec);
            rhoInterp = rhoVec(sortIdx);
            [PphiInterp, uniqueIdx] = unique(PphiInterp, 'stable');
            rhoInterp = rhoInterp(uniqueIdx);
            if numel(PphiInterp) < 2
                continue;
            end

            minPphiVec = PphiInterp(1);
            maxPphiVec = PphiInterp(end);
            Lambda_k = gridLambda(k);

            for j = 1:sizePphi
                pphiVal = gridPphi(j);
                if pphiVal < minPphiVec || pphiVal > maxPphiVec
                    continue;
                end

                rhoVal = interp1(PphiInterp, rhoInterp, pphiVal, 'pchip');
                if ~isfinite(rhoVal)
                    continue;
                end

                BInterp = interp1(rhoRadial, BRadial, rhoVal, 'pchip');
                SFAcovyzInterp = interp1(rhoRadial, SFAcovyzRadial, rhoVal, 'pchip');
                qInterp = interp1(rhoRadial, qRadial, rhoVal, 'pchip');
                JInterp = interp1(rhoRadial, JRadial, rhoVal, 'pchip');
                psipInterp = interp1(rhoRadial, psipRadial, rhoVal, 'pchip');

                sqrtArg = 1 - Lambda_k * BInterp;
                if sqrtArg < -1e-12
                    continue;
                end
                sqrtArg = max(sqrtArg, 0);

                vparaVal = branchSign * sqrt(2 * Ei * sqrtArg / species.Mass);
                PphiCalc = cm * species.Mass * vparaVal .* 2 * psitmax * drho .* rhoVal .* SFAcovyzInterp ./ ...
                    (qInterp .* JInterp .* BInterp) - species.Char .* psipInterp;

                % Pphi 远离 0 时用相对误差，接近 0 时用绝对误差。
                absErr = abs(PphiCalc - pphiVal);
                if abs(pphiVal) > 1e-12
                    pphiError = absErr / abs(pphiVal);
                    if pphiError > 1e-3
                        continue;
                    end
                else
                    pphiError = absErr;
                    if pphiError > 1e-12
                        continue;
                    end
                end

                rhoNorm = (rhoVal - eq.RHO0) / (eq.RHO1 - eq.RHO0);
                if abs(rhoNorm - 1) < 1e-8
                    rhoNorm = 1 - 1e-8;
                end
                if rhoNorm < 0 || rhoNorm >= 1
                    continue;
                end

                linearIdx = (j - 1) * sizeLambda + k;
                rhoLocal(linearIdx) = rhoNorm;
                vparaLocal(linearIdx) = vparaVal;
                muLocal(linearIdx) = Ei * Lambda_k;
                idLocal(linearIdx) = int32((i - 1) * M + linearIdx - 1);
                errorLocal(linearIdx) = pphiError;
            end
        end

        rhoCell{i} = rhoLocal;
        vparaCell{i} = vparaLocal;
        muCell{i} = muLocal;
        idCell{i} = idLocal;
        errorCell{i} = errorLocal;
    end

    rhoArr = vertcat(rhoCell{:});
    vparaArr = vertcat(vparaCell{:});
    muArr = vertcat(muCell{:});
    idArr = vertcat(idCell{:});
    PphiErrorArr = vertcat(errorCell{:});

    nonemptyRatio = sum(idArr ~= padId) / N;
    fprintf('该分支有效映射比例：%.2f%%\n', nonemptyRatio * 100);
end
