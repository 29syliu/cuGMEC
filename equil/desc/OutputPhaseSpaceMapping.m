
%%

MP = mi;
QE = e;
B = B / B0;
J = Jxyz / J0;
psip = psip / (B0*L0^2);
SFAcovyz = gcovSFA{2,3} / L0^2;
SFAcovyz = SFAcovyz(:,3:66);
theta = theta_pest;

save(fullfile(outputPath, 'phaseSpaceVariables.mat'), ...
    'MP', 'QE', 'B', 'J', 'psip', 'SFAcovyz', ...
    'B0', 'L0', 'VA0', 'RHO0', 'RHO1', 'PSITMAX', ...
    'q', 'R', 'Z', 'rho', 'theta');

%%

load('phaseSpaceVariables.mat')

%%

thetaIndex = 33;

figure;
hold on;
plot(R(1,:),Z(1,:),'r');axis equal;
plot(R(end,:),Z(end,:),'r');axis equal;
plot(R(:,thetaIndex),Z(:,thetaIndex),'r');axis equal;

%%

gridE = 64;
gridPphi = 128;
gridLambda = 128;

IonMass = 4;
IonChar = 2;
IonVmin = 0.0733;
IonVmax = 1.466;

%%

[paraRho, paraVpara, paraE, paraPphi, paraLambda, paraMu, paraId, paraError, paraNonemptyRatio, ...
    minE, maxE, minPphi, maxPphi, minLambda, maxLambda] = ...
    computePhaseSpaceMapping(IonMass, IonChar, IonVmin, IonVmax, 1, thetaIndex, gridE, gridPphi, gridLambda);

[antiRho, antiVanti, antiE, antiPphi, antiLambda, antiMu, antiId, antiError, antiNonemptyRatio] = ...
    computePhaseSpaceMapping(IonMass, IonChar, IonVmin, IonVmax, -1, thetaIndex, gridE, gridPphi, gridLambda);

%%

[allRho, allVpara, allE, allPphi, allLambda, allMu, allId, allError] = ...
    extractMappings(paraRho, paraVpara, paraE, paraPphi, paraLambda, paraMu, paraId, paraError, ...
    antiRho, antiVanti, antiE, antiPphi, antiLambda, antiMu, antiId, antiError);

%%

min([allRho]) > double(0)
max([allRho]) < double(1)

sum(isnan(allId))+sum(isnan(allRho))+sum(isnan(allVpara))+sum(isnan(allMu))+sum(isnan(allE))+sum(isnan(allPphi))+sum(isnan(allLambda))

fid = fopen('AlphaPhaseSpaceMapping.bin', 'wb');
if fid == -1
    error('无法创建文件');
end

fwrite(fid, allId, 'int32');
fwrite(fid, allRho, 'double');
fwrite(fid, allVpara, 'double');
fwrite(fid, allMu, 'double');

fclose(fid);

%%

function [rhoArr, vparaArr, EArr, PphiArr, LambdaArr, muArr, idArr, PphiErrorArr, nonemptyRatio, ...
    minE, maxE, minPphi, maxPphi, minLambda, maxLambda] = ...
    computePhaseSpaceMapping(Mass, Char, Vmin, Vmax, Sign, thetaIndex, sizeE, sizePphi, sizeLambda)
    
    B       = evalin('base', 'B');
    B0      = evalin('base', 'B0');
    J       = evalin('base', 'J');
    L0      = evalin('base', 'L0');
    MP      = evalin('base', 'MP');
    psip    = evalin('base', 'psip');
    PSITMAX = evalin('base', 'PSITMAX');
    q       = evalin('base', 'q');
    QE      = evalin('base', 'QE');
    rho     = evalin('base', 'rho');
    RHO0    = evalin('base', 'RHO0');
    RHO1    = evalin('base', 'RHO1');
    SFAcovyz= evalin('base', 'SFAcovyz');
    VA0     = evalin('base', 'VA0');

    drho = RHO1 - RHO0;
    psitmax = PSITMAX / (B0 * L0 * L0);
    cm = VA0 / (L0 * (QE * B0 / MP));
    
    minE = 0.5 * Mass * Vmin^2;
    maxE = 0.5 * Mass * Vmax^2;
    
    Pphi0 = -cm * Mass * Vmax * 2 * psitmax * drho .* rho .* SFAcovyz ./ (q .* J .* B) - Char .* psip;
    Pphi1 =  cm * Mass * Vmax * 2 * psitmax * drho .* rho .* SFAcovyz ./ (q .* J .* B) - Char .* psip;
    
    minPphi = min(min(Pphi0));
    maxPphi = max(max(Pphi1));
    
    minLambda = 0.0;
    maxLambda = max(max(1./B));
    
    gridE = linspace(minE, maxE, sizeE);
    gridPphi = linspace(minPphi, maxPphi, sizePphi);
    gridLambda = linspace(minLambda, maxLambda, sizeLambda);
    
    rho_radial = rho(:, thetaIndex);
    B_radial = B(:, thetaIndex);
    SFAcovyz_radial = SFAcovyz(:, thetaIndex);
    q_radial = q(:, thetaIndex);
    J_radial = J(:, thetaIndex);
    psip_radial = psip(:, thetaIndex);
    
    countValid0 = 0;
    countValid1 = 0;
    minValidGt1 = inf;
    nonmonotonicCount = 0;
    singleExtremaCount = 0;
    multipleExtremaCount = 0;
    percentArray = [];
    hasEqualAdjacent = false;
    equalAdjacentCount = 0;
    
    validRhoPphi = cell(length(gridE), length(gridLambda));
    
    for i = 1:length(gridE)
        E = gridE(i);
        for j = 1:length(gridLambda)
            Lambda = gridLambda(j);
            sqrtArg = 1 - Lambda * B_radial;
            valid = sqrtArg >= 0;
            validCount = sum(valid);
            
            if validCount == 0
                countValid0 = countValid0 + 1;
                continue;
            elseif validCount == 1
                countValid1 = countValid1 + 1;
                continue;
            elseif validCount > 1
                if validCount < minValidGt1
                    minValidGt1 = validCount;
                end
            end
            
            rhoValid = rho_radial(valid);
            BValid = B_radial(valid);
            SFAcovyzValid = SFAcovyz_radial(valid);
            qValid = q_radial(valid);
            JValid = J_radial(valid);
            psipValid = psip_radial(valid);
            
            Vpara = Sign * sqrt(2 * E * sqrtArg(valid) / Mass);
            Pphi = cm * Mass * Vpara .* 2 * psitmax * drho .* rhoValid .* SFAcovyzValid ./ ...
                   (qValid .* JValid .* BValid) - Char .* psipValid;
            
            diffPphi = diff(Pphi);
            if any(abs(diffPphi) < 1e-12)
                hasEqualAdjacent = true;
                equalAdjacentCount = equalAdjacentCount + 1;
            end
            
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
                    changes = sum(abs(diff(filteredSign)) > 0);
                    numExtrema = changes;
                end
                
                if numExtrema == 1
                    singleExtremaCount = singleExtremaCount + 1;
                    [~, idxMax] = max(Pphi);
                    [~, idxMin] = min(Pphi);
                    if idxMax > 1 && idxMax < length(Pphi)
                        extremumIdx = idxMax;
                    elseif idxMin > 1 && idxMin < length(Pphi)
                        extremumIdx = idxMin;
                    else
                        continue;
                    end
                    
                    totalRho = rhoValid(end) - rhoValid(1);
                    rhoLeft = rhoValid(extremumIdx) - rhoValid(1);
                    rhoRight = rhoValid(end) - rhoValid(extremumIdx);
                    percentVal = min(rhoLeft, rhoRight) / totalRho * 100;
                    percentArray(end+1) = percentVal;
                    
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
                validRhoPphi{i, j} = struct('rho', storedRho, 'Pphi', storedPphi);
            end
        end
    end
    
    fprintf('--- 有效rho统计 ---\n');
    fprintf('有效rho为0的(E, Lambda)组合数：%d\n', countValid0);
    fprintf('有效rho为1的(E, Lambda)组合数：%d\n', countValid1);
    if isinf(minValidGt1)
        fprintf('没有有效rho大于1的组合。\n');
    else
        fprintf('有效rho大于1时，最少的有效rho为：%d\n', minValidGt1);
    end
    
    fprintf('\n--- 相邻Pphi相等统计 ---\n');
    if hasEqualAdjacent
        fprintf('存在相邻两个有效rho点对应的Pphi相等的(E, Lambda)组合\n');
        fprintf('出现这种情况的组合总数为：%d\n', equalAdjacentCount);
    else
        fprintf('不存在相邻两个有效rho点对应的Pphi相等\n');
    end
    fprintf('\n');
    
    fprintf('发现%d个(E, Lambda)组合下，Pphi相对rho是非单调的\n', nonmonotonicCount);
    fprintf('其中只含一个极值点的有%d个，含多个极值点的有%d个\n', singleExtremaCount, multipleExtremaCount);
    if ~isempty(percentArray)
        maxPercent = max(percentArray);
        meanPercent = mean(percentArray);
        fprintf('单极值点中，两端较窄一侧的rho范围占总的有效rho范围的百分比：最大值为%.2f%%，平均值为%.2f%%\n', maxPercent, meanPercent);
    end
    
    function rhoVal = Pphi2Rho(i, j, pphiVal, validRhoPphi)
        PphiVec = validRhoPphi{i,j}.Pphi;
        rhoVec = validRhoPphi{i,j}.rho;
        if min(pphiVal) < min(PphiVec) || max(pphiVal) > max(PphiVec)
            rhoVal = [];
            return;
        end
        rhoVal = interp1(PphiVec, rhoVec, pphiVal, 'spline');
    end

    M = sizePphi * sizeLambda;
    N = sizeE * M;
    
    rhoCell = cell(sizeE, 1);
    vparaCell = cell(sizeE, 1);
    ECell = cell(sizeE, 1);
    PphiCell = cell(sizeE, 1);
    LambdaCell = cell(sizeE, 1);
    muCell = cell(sizeE, 1);
    idCell = cell(sizeE, 1);
    errorCell = cell(sizeE, 1);
    
    parfor i = 1:sizeE

        rhoLocal = zeros(M, 1);
        vparaLocal = zeros(M, 1);
        ELocal = zeros(M, 1);
        PphiLocal = zeros(M, 1);
        LambdaLocal = zeros(M, 1);
        muLocal = zeros(M, 1);
        idLocal = int32(zeros(M, 1));
        errorLocal = zeros(M, 1);
        
        defaultIdLocal = int32(20251106);
        rhoLocal(:) = 0.5;
        vparaLocal(:) = 0;
        ELocal(:) = 0;
        PphiLocal(:) = 0;
        LambdaLocal(:) = 0;
        muLocal(:) = 0;
        idLocal(:) = defaultIdLocal;
        errorLocal(:) = 0;
        
        Ei = gridE(i);
        for k = 1:sizeLambda
            if isempty(validRhoPphi{i, k})
                continue;
            end
            PphiVec = validRhoPphi{i, k}.Pphi;
            rhoVec = validRhoPphi{i, k}.rho;
            minPphiVec = min(PphiVec);
            maxPphiVec = max(PphiVec);
            Lambda_k = gridLambda(k);
            for j = 1:sizePphi
                pphiVal = gridPphi(j);
                if pphiVal < minPphiVec || pphiVal > maxPphiVec
                    continue;
                end
                rhoVal = interp1(PphiVec, rhoVec, pphiVal, 'spline');
                if isnan(rhoVal)
                    continue;
                end

                BInterp = interp1(rho_radial, B_radial, rhoVal, 'spline');
                SFAcovyzInterp = interp1(rho_radial, SFAcovyz_radial, rhoVal, 'spline');
                qInterp = interp1(rho_radial, q_radial, rhoVal, 'spline');
                JInterp = interp1(rho_radial, J_radial, rhoVal, 'spline');
                psipInterp = interp1(rho_radial, psip_radial, rhoVal, 'spline');
                
                sqrtArg = 1 - Lambda_k * BInterp;
                if sqrtArg < 0
                    continue;
                end

                vparaVal = Sign * sqrt(2 * Ei * sqrtArg / Mass);
                PphiCalc = cm * Mass * vparaVal .* 2 * psitmax * drho .* rhoVal .* SFAcovyzInterp ./ ...
                    (qInterp .* JInterp .* BInterp) - Char .* psipInterp;
                pphiError = abs(PphiCalc - pphiVal) / abs(pphiVal);
                if pphiError > 1e-3
                    continue;
                end

                rhoNorm = (rhoVal - RHO0) / (RHO1 - RHO0);
                if abs(rhoNorm-1) < 1e-8
                    rhoNorm = 1 - 1e-8;
                end
                
                linearIdx = (j-1) * sizeLambda + k;
                rhoLocal(linearIdx) = rhoNorm;
                vparaLocal(linearIdx) = vparaVal;
                ELocal(linearIdx) = Ei;
                PphiLocal(linearIdx) = pphiVal;
                LambdaLocal(linearIdx) = Lambda_k;
                muLocal(linearIdx) = Ei * Lambda_k;
                idLocal(linearIdx) = int32((i-1)*M + linearIdx - 1);
                errorLocal(linearIdx) = pphiError;
            end
        end

        rhoCell{i} = rhoLocal;
        vparaCell{i} = vparaLocal;
        ECell{i} = ELocal;
        PphiCell{i} = PphiLocal;
        LambdaCell{i} = LambdaLocal;
        muCell{i} = muLocal;
        idCell{i} = idLocal;
        errorCell{i} = errorLocal;
    end
    
    rhoArr = vertcat(rhoCell{:});
    vparaArr = vertcat(vparaCell{:});
    EArr = vertcat(ECell{:});
    PphiArr = vertcat(PphiCell{:});
    LambdaArr = vertcat(LambdaCell{:});
    muArr = vertcat(muCell{:});
    idArr = vertcat(idCell{:});
    PphiErrorArr = vertcat(errorCell{:});

    nonemptyRatio = sum(idArr ~= int32(20251106)) / N;
    fprintf('相空间映射非空比例: %.2f%%\n\n', nonemptyRatio * 100);

end

%%

function [rhoAll, vparaAll, EAll, PphiAll, LambdaAll, muAll, idAll, PphiErrorAll] = ...
    extractMappings(rhoArr1, vparaArr1, EArr1, PphiArr1, LambdaArr1, muArr1, idArr1, pphiErrorArr1, ...
                    rhoArr2, vparaArr2, EArr2, PphiArr2, LambdaArr2, muArr2, idArr2, pphiErrorArr2)

    N1 = length(rhoArr1);
    if nargin < 9
        N2 = N1;
        rhoArr2 = ones(N2,1) * 0.5;
        vparaArr2 = zeros(N2,1);
        EArr2 = zeros(N2,1);
        PphiArr2 = zeros(N2,1);
        LambdaArr2 = zeros(N2,1);
        muArr2 = zeros(N2,1);
        idArr2 = ones(N2,1) * 20251106;
        pphiErrorArr2 = zeros(N2,1);
    end
    
    rhoAll = [rhoArr1; rhoArr2];
    vparaAll = [vparaArr1; vparaArr2];
    EAll = [EArr1; EArr2];
    PphiAll = [PphiArr1; PphiArr2];
    LambdaAll = [LambdaArr1; LambdaArr2];
    muAll = [muArr1; muArr2];
    idAll = [idArr1; idArr2];
    PphiErrorAll = [pphiErrorArr1; pphiErrorArr2];
end