fid = fopen('AlphaPhaseSpaceOrbit.bin', 'rb');
if fid == -1
    error('无法打开文件 PhaseSpaceOrbit.bin');
end

fseek(fid, 0, 'eof');
fileSize = ftell(fid);
fseek(fid, 0, 'bof');

bytesPerParticle = 9 * 8;
if mod(fileSize, bytesPerParticle) ~= 0
    fclose(fid);
    error('文件大小不是 %d 的整数倍，数据可能不完整', bytesPerParticle);
end
numParticles = fileSize / bytesPerParticle;

data = fread(fid, numParticles * 9, 'double');
fclose(fid);

data = reshape(data, 9, numParticles)';

Ids       = data(:, 1);
orbits    = data(:, 2);
dtheta    = data(:, 3);
dphiTotal = data(:, 4);
dphiVpara = data(:, 5);
dTs       = data(:, 6);
Es        = data(:, 7);
Pphis     = data(:, 8);
Lambdas   = data(:, 9);
Ids = int32(Ids);

E1d = linspace(minE, maxE, gridE);
Pphi1d = linspace(minPphi, maxPphi, gridPphi);
Lambda1d = linspace(minLambda, maxLambda, gridLambda);

%%

sum(isnan(Ids))+sum(isnan(orbits))+sum(isnan(dtheta))+sum(isnan(dphiTotal))+sum(isnan(dphiVpara))+ ...
    sum(isnan(dTs))+sum(isnan(Es))+sum(isnan(Pphis))+sum(isnan(Lambdas))

%% 分类粒子

idPad = orbits<1;
idLoss = orbits>1&orbits<2;
idUnknown = orbits>5;

idPara = orbits>2&orbits<3;
idAnti = orbits>3&orbits<4;
idTrapped = orbits>4&orbits<5;

sum(idLoss)
sum(idUnknown)
sum(idPara)
sum(idAnti)
sum(idTrapped/2)

%% 捕获粒子提取

valid = Ids ~= 20251106;
tempIds = Ids(valid);
tempOrbits = orbits(valid);
tempDTs = dTs(valid);

[unique_ids, ~, ic] = unique(tempIds);
counts = accumarray(ic, 1);
twice_mask = counts == 2;
twice_ids = unique_ids(twice_mask);

idx_groups = accumarray(ic, (1:length(tempIds))', [], @(x) {x});
idx_groups = idx_groups(twice_mask);

n_twice = sum(twice_mask);
trapped_ids = int32(zeros(n_twice, 1));
trapped_orbits1 = zeros(n_twice, 1);
trapped_orbits2 = zeros(n_twice, 1);
trapped_dTs1 = zeros(n_twice, 1);
trapped_dTs2 = zeros(n_twice, 1);

tol = 1e-12;
cnt = 0;

for i = 1:n_twice

    idx = idx_groups{i};
    if abs(tempOrbits(idx(1)) - 4.5) < tol && abs(tempOrbits(idx(2)) - 4.5) < tol
        cnt = cnt + 1;
        trapped_ids(cnt) = twice_ids(i);
        trapped_orbits1(cnt) = tempOrbits(idx(1));
        trapped_orbits2(cnt) = tempOrbits(idx(2));
        trapped_dTs1(cnt) = tempDTs(idx(1));
        trapped_dTs2(cnt) = tempDTs(idx(2));
    end

end

trapped_ids = trapped_ids(1:cnt);
trapped_orbits1 = trapped_orbits1(1:cnt);
trapped_orbits2 = trapped_orbits2(1:cnt);
trapped_dTs1 = trapped_dTs1(1:cnt);
trapped_dTs2 = trapped_dTs2(1:cnt);

result_table = table(trapped_ids', trapped_orbits1', trapped_orbits2', trapped_dTs1', trapped_dTs2', ...
                         'VariableNames', {'ID', 'Orbit1', 'Orbit2', 'dT1', 'dT2'});

figure;
plot(abs(trapped_dTs1-trapped_dTs2)./trapped_dTs2)
sum(abs((trapped_dTs1-trapped_dTs2)./trapped_dTs2) < 0.05)/size(trapped_dTs2,1)


%% 捕获粒子诊断

validTrappedIndex = abs((trapped_dTs1 - trapped_dTs2) ./ trapped_dTs2) < 0.05;
validTrappedId = trapped_ids(validTrappedIndex);

validTrappedIndex2 = findIndex2(Ids, validTrappedId);
max(validTrappedIndex2(:,1)-validTrappedIndex2(:,2))
min(validTrappedIndex2(:,1)-validTrappedIndex2(:,2))

figure;
plot((Es(validTrappedIndex2(:,1))-allE(validTrappedIndex2(:,1)))./allE(validTrappedIndex2(:,1)));
figure;
plot((Es(validTrappedIndex2(:,2))-allE(validTrappedIndex2(:,1)))./allE(validTrappedIndex2(:,1)));

figure;
plot((Pphis(validTrappedIndex2(:,1))-allPphi(validTrappedIndex2(:,1)))./allPphi(validTrappedIndex2(:,1)));
figure;
plot((Pphis(validTrappedIndex2(:,2))-allPphi(validTrappedIndex2(:,1)))./allPphi(validTrappedIndex2(:,1)));

figure;
plot((Lambdas(validTrappedIndex2(:,1))-allLambda(validTrappedIndex2(:,1)))./allLambda(validTrappedIndex2(:,1)));
figure;
plot((Lambdas(validTrappedIndex2(:,2))-allLambda(validTrappedIndex2(:,1)))./allLambda(validTrappedIndex2(:,1)));

figure;
plot(dtheta(validTrappedIndex2(:,1)))
figure;
plot(dtheta(validTrappedIndex2(:,2)))

%% 捕获粒子相空间

validTrappedDT = (dTs(validTrappedIndex2(:,1))+dTs(validTrappedIndex2(:,2)))/2;
validTrappedBounceT = zeros(gridE,gridPphi,gridLambda);

for index = 1:size(validTrappedIndex2(:,1),1)
    [i,j,k] = idx2ijk(validTrappedIndex2(index,1),gridE,gridPphi,gridLambda);
    validTrappedBounceT(i,j,k) = validTrappedDT(index);
end

%% 捕获粒子绘图

plotSlice(1,1,50,validTrappedBounceT, E1d, Pphi1d, Lambda1d)

%% 同向通行粒子2提取

valid = Ids ~= 20251106;
tempIds = Ids(valid);
tempOrbits = orbits(valid);
tempDTs = dTs(valid);

[unique_ids, ~, ic] = unique(tempIds);
counts = accumarray(ic, 1);
twice_mask = counts == 2;
twice_ids = unique_ids(twice_mask);

idx_groups = accumarray(ic, (1:length(tempIds))', [], @(x) {x});
idx_groups = idx_groups(twice_mask);

n_twice = sum(twice_mask);
para_ids = int32(zeros(n_twice, 1));
para_orbits1 = zeros(n_twice, 1);
para_orbits2 = zeros(n_twice, 1);
para_dTs1 = zeros(n_twice, 1);
para_dTs2 = zeros(n_twice, 1);

tol = 1e-12;
cnt = 0;

pos_25 = [];

for i = 1:n_twice

    idx = idx_groups{i};
    o1 = tempOrbits(idx(1));
    o2 = tempOrbits(idx(2));
    is25_1 = abs(o1 - 2.5) < tol;
    is25_2 = abs(o2 - 2.5) < tol;
    isOther_1 = abs(o1 - 1.5) < tol || abs(o1 - 3.5) < tol || abs(o1 - 5.5) < tol;
    isOther_2 = abs(o2 - 1.5) < tol || abs(o2 - 3.5) < tol || abs(o2 - 5.5) < tol;
    
    if (is25_1 && isOther_2) || (is25_2 && isOther_1)
        cnt = cnt + 1;
        para_ids(cnt) = twice_ids(i);
        para_orbits1(cnt) = o1;
        para_orbits2(cnt) = o2;
        para_dTs1(cnt) = tempDTs(idx(1));
        para_dTs2(cnt) = tempDTs(idx(2));
        
        if is25_1
            pos_25(end+1) = 1;
        else
            pos_25(end+1) = 2;
        end
    end
end

para_ids = para_ids(1:cnt);
para_orbits1 = para_orbits1(1:cnt);
para_orbits2 = para_orbits2(1:cnt);
para_dTs1 = para_dTs1(1:cnt);
para_dTs2 = para_dTs2(1:cnt);

if cnt > 0
    result_table = table(para_ids', para_orbits1', para_orbits2', para_dTs1', para_dTs2', ...
                         'VariableNames', {'ID', 'Orbit1', 'Orbit2', 'dT1', 'dT2'});
    
    if all(pos_25 == 1)
        fprintf('2.5 全部出现在 Orbit1 列（共 %d 个）。\n', cnt);
    elseif all(pos_25 == 2)
        fprintf('2.5 全部出现在 Orbit2 列（共 %d 个）。\n', cnt);
    else
        num_col1 = sum(pos_25 == 1);
        num_col2 = sum(pos_25 == 2);
        fprintf('2.5 不是全部在同一列。在 Orbit1 列出现 %d 次 (%.1f%%)，在 Orbit2 列出现 %d 次 (%.1f%%)。\n', ...
                num_col1, 100*num_col1/cnt, num_col2, 100*num_col2/cnt);
    end
    
   
else
    disp('没有找到符合条件的配对（一个2.5，另一个为1.5/3.5/5.5）。');
end

%% 同向通行粒子2诊断

validParaId = para_ids;

validParaIndex2 = findIndex2(Ids, validParaId);
max(validParaIndex2(:,1)-validParaIndex2(:,2))
min(validParaIndex2(:,1)-validParaIndex2(:,2))

figure;
plot((Es(validParaIndex2(:,1))-allE(validParaIndex2(:,1)))./allE(validParaIndex2(:,1)));

figure;
plot((Pphis(validParaIndex2(:,1))-allPphi(validParaIndex2(:,1)))./allPphi(validParaIndex2(:,1)));

figure;
plot((Lambdas(validParaIndex2(:,1))-allLambda(validParaIndex2(:,1)))./allLambda(validParaIndex2(:,1)));

figure;
plot(dtheta(validParaIndex2(:,1)))

figure;
plot(dphiTotal(validParaIndex2(:,1)))

%% 同向通行粒子1提取

[unique_ids, ~, ic] = unique(tempIds);
counts = accumarray(ic, 1);
single_ids = unique_ids(counts == 1);

tol = 1e-12;
mask = ismember(tempIds, single_ids) & (abs(tempOrbits - 2.5) < tol);
orbit25_single_ids = tempIds(mask);

fprintf('只出现一次且 orbit 为 2.5 的 ID 数量：%d\n', length(orbit25_single_ids));

validParaIndex1 = findIndex1(Ids, orbit25_single_ids);

%% 同向通行粒子1诊断

figure;
plot((Es(validParaIndex1)-allE(validParaIndex1))./allE(validParaIndex1));

figure;
plot((Pphis(validParaIndex1)-allPphi(validParaIndex1))./allPphi(validParaIndex1));

figure;
plot((Lambdas(validParaIndex1)-allLambda(validParaIndex1))./allLambda(validParaIndex1));

figure;
plot(dtheta(validParaIndex1))

figure;
plot(dphiTotal(validParaIndex1))

%% 同向通行粒子相空间

validPara2DT = dTs(validParaIndex2(:,1));
validParaTransitT = zeros(gridE,gridPphi,gridLambda);

for index = 1:size(validParaIndex2(:,1),1)
    [i,j,k] = idx2ijk(validParaIndex2(index,1),gridE,gridPphi,gridLambda);
    validParaTransitT(i,j,k) = validPara2DT(index);
end

validPara1DT = dTs(validParaIndex1);

for index = 1:size(validParaIndex1,1)
    [i,j,k] = idx2ijk(validParaIndex1(index),gridE,gridPphi,gridLambda);
    validParaTransitT(i,j,k) = validPara1DT(index);
end

%% 同向通行粒子绘图

plotSlice(1,1,50,validParaTransitT, E1d, Pphi1d, Lambda1d)

%% 反向通行粒子2提取

valid = Ids ~= 20251106;
tempIds = Ids(valid);
tempOrbits = orbits(valid);
tempDTs = dTs(valid);
tempEs = Es(valid);
tempPphis = Pphis(valid);
tempLambdas = Lambdas(valid);

[unique_ids, ~, ic] = unique(tempIds);
counts = accumarray(ic, 1);
twice_mask = counts == 2;
twice_ids = unique_ids(twice_mask);

idx_groups = accumarray(ic, (1:length(tempIds))', [], @(x) {x});
idx_groups = idx_groups(twice_mask);

n_twice = sum(twice_mask);
anti_ids = zeros(n_twice, 1);
anti_orbits1 = zeros(n_twice, 1);
anti_orbits2 = zeros(n_twice, 1);
anti_dTs1 = zeros(n_twice, 1);
anti_dTs2 = zeros(n_twice, 1);

tol = 1e-12;
cnt = 0;
pos_35 = [];   % 记录3.5出现的位置

for i = 1:n_twice
    idx = idx_groups{i};
    o1 = tempOrbits(idx(1));
    o2 = tempOrbits(idx(2));
    is35_1 = abs(o1 - 3.5) < tol;
    is35_2 = abs(o2 - 3.5) < tol;
    isOther_1 = abs(o1 - 1.5) < tol || abs(o1 - 2.5) < tol || abs(o1 - 5.5) < tol;
    isOther_2 = abs(o2 - 1.5) < tol || abs(o2 - 2.5) < tol || abs(o2 - 5.5) < tol;
    
    if (is35_1 && isOther_2) || (is35_2 && isOther_1)
        cnt = cnt + 1;
        anti_ids(cnt) = twice_ids(i);
        anti_orbits1(cnt) = o1;
        anti_orbits2(cnt) = o2;
        anti_dTs1(cnt) = tempDTs(idx(1));
        anti_dTs2(cnt) = tempDTs(idx(2));
        
        if is35_1
            pos_35(end+1) = 1;
        else
            pos_35(end+1) = 2;
        end
    end
end

anti_ids = anti_ids(1:cnt);
anti_orbits1 = anti_orbits1(1:cnt);
anti_orbits2 = anti_orbits2(1:cnt);
anti_dTs1 = anti_dTs1(1:cnt);
anti_dTs2 = anti_dTs2(1:cnt);

if cnt > 0
    result_table = table(anti_ids', anti_orbits1', anti_orbits2', anti_dTs1', anti_dTs2', ...
                         'VariableNames', {'ID', 'Orbit1', 'Orbit2', 'dT1', 'dT2'});
    
    if all(pos_35 == 1)
        fprintf('3.5 全部出现在 Orbit1 列（共 %d 个）。\n', cnt);
    elseif all(pos_35 == 2)
        fprintf('3.5 全部出现在 Orbit2 列（共 %d 个）。\n', cnt);
    else
        num_col1 = sum(pos_35 == 1);
        num_col2 = sum(pos_35 == 2);
        fprintf('3.5 不是全部在同一列。在 Orbit1 列出现 %d 次 (%.1f%%)，在 Orbit2 列出现 %d 次 (%.1f%%)。\n', ...
                num_col1, 100*num_col1/cnt, num_col2, 100*num_col2/cnt);
    end
else
    disp('没有找到符合条件的配对（一个3.5，另一个为1.5/2.5/5.5）。');
end

%% 反向通行粒子2诊断

validAntiId = anti_ids;

validAntiIndex2 = findIndex2(Ids, validAntiId);
max(validAntiIndex2(:,1)-validAntiIndex2(:,2))
min(validAntiIndex2(:,1)-validAntiIndex2(:,2))

figure;
plot((Es(validAntiIndex2(:,2))-allE(validAntiIndex2(:,1)))./allE(validAntiIndex2(:,1)));

figure;
plot((Pphis(validAntiIndex2(:,2))-allPphi(validAntiIndex2(:,1)))./allPphi(validAntiIndex2(:,1)));

figure;
plot((Lambdas(validAntiIndex2(:,2))-allLambda(validAntiIndex2(:,1)))./allLambda(validAntiIndex2(:,1)));

figure;
plot(dtheta(validAntiIndex2(:,2)))

figure;
plot(dphiTotal(validAntiIndex2(:,2)))

%% 反向通行粒子1提取

[unique_ids, ~, ic] = unique(tempIds);
counts = accumarray(ic, 1);
single_ids = unique_ids(counts == 1);

tol = 1e-12;
mask = ismember(tempIds, single_ids) & (abs(tempOrbits - 3.5) < tol);
orbit35_single_ids = tempIds(mask);

fprintf('只出现一次且 orbit 为 3.5 的 ID 数量：%d\n', length(orbit35_single_ids));

validAntiIndex1 = findIndex1(Ids, orbit35_single_ids);

%% 反向通行粒子1诊断

figure;
plot((Es(validAntiIndex1)-allE(validAntiIndex1))./allE(validAntiIndex1));

figure;
plot((Pphis(validAntiIndex1)-allPphi(validAntiIndex1))./allPphi(validAntiIndex1));

figure;
plot((Lambdas(validAntiIndex1)-allLambda(validAntiIndex1))./allLambda(validAntiIndex1));

figure;
plot(dtheta(validAntiIndex1))

figure;
plot(dphiTotal(validAntiIndex1))

%% 反向通行粒子相空间

validAnti2DT = dTs(validAntiIndex2(:,2));
validAntiTransitT = zeros(gridE,gridPphi,gridLambda);

for index = 1:size(validAntiIndex2(:,1),1)
    [i,j,k] = idx2ijk(validAntiIndex2(index,1),gridE,gridPphi,gridLambda);
    validAntiTransitT(i,j,k) = validAnti2DT(index);
end

validAnti1DT = dTs(validAntiIndex1);

for index = 1:size(validAntiIndex1,1)
    [i,j,k] = idx2ijk(validAntiIndex1(index)-gridE*gridPphi*gridLambda,gridE,gridPphi,gridLambda);
    validAntiTransitT(i,j,k) = validAnti1DT(index);
end

%% 反向通行粒子绘图

plotSlice(-1,1,50,validAntiTransitT, E1d, Pphi1d, Lambda1d)

%% 相空间有效粒子统计

(nnz(validAntiTransitT)+nnz(validParaTransitT)+2*nnz(validTrappedBounceT))/(gridE*gridPphi*gridLambda*2)

sum(idLoss) / (gridE*gridPphi*gridLambda*2)

nnz((validAntiTransitT ~= 0) & (validTrappedBounceT ~= 0))

nnz((validParaTransitT ~= 0) & (validTrappedBounceT ~= 0))

load("C:\Users\ALFVEN\Desktop\cuGMEC 1.1\MATLAB\colorZQ.mat")

[nnz(validParaTransitT), sum(idPara)]
[nnz(validAntiTransitT), sum(idAnti)]
[nnz(validTrappedBounceT), sum(idTrapped/2)]

%% vpara > 0

validT = (validParaTransitT + validTrappedBounceT)*L0/VA0;

plotSlice(1,1,40,validT, E1d, Pphi1d, Lambda1d)

colormap(colorZQ);

%% vpara < 0

validT = (validAntiTransitT + validTrappedBounceT)*L0/VA0;

plotSlice(-1,1,40,validT, E1d, Pphi1d, Lambda1d)

colormap(colorZQ);

%% valid

valid = paraId ~= 20251106;
valid = paraId(valid) + 1;

validPara = zeros(gridE,gridPphi,gridLambda);

for index = 1:size(valid,1)
    [i,j,k] = idx2ijk(int32(valid(index)),gridE,gridPphi,gridLambda);
    validPara(i,j,k) = 1;
end

%%
int32(valid(end))
ans = idx2ijk(1048332,gridE,gridPphi,gridLambda)
ans = idx2ijk(int32(valid(end)),gridE,gridPphi,gridLambda)


%%

valid = allId((allId ~= 20251106) & (allVpara < 0));

hhh2 = findIndex1(allId, valid);

%%

validAnti = zeros(gridE,gridPphi,gridLambda);

for index = 1:size(hhh2,1)
    [i,j,k] = idx2ijk(hhh2(index),gridE,gridPphi,gridLambda);
    validAnti(i,j,k) = 1;
end

%%
plotSlice(1,1,40,validPara, E1d, Pphi1d, Lambda1d)

%%

function result = findIndex2(a, b)
    [unique_vals, ~, ic] = unique(a);
    indices = accumarray(ic, (1:length(a))', [], @(x) {int32(x)});
    [~, loc] = ismember(b, unique_vals);
    n = length(b);
    result = zeros(n, 2, 'int32');
    for i = 1:n
        idx_list = indices{loc(i)};
        result(i, :) = idx_list(1:2);
    end
end

function result = findIndex1(a, b)
    [unique_vals, ~, ic] = unique(a);
    firstIdx = accumarray(ic, (1:length(a))', [], @min);
    [~, loc] = ismember(b, unique_vals);
    result = firstIdx(loc);
end

%%

function [i, j, k] = idx2ijk(Id, sizeE, sizePphi, sizeLambda)

    A = sizePphi * sizeLambda;
    B = sizeLambda;

    i = floor(Id / A) + 1;
    remainder = mod(Id, A);
    j = floor(remainder / B) + 1;
    k = mod(remainder, B) + 1;

end

%%

function plotSlice(sign, dim, idx, bounceT, E1d, Pphi1d, Lambda1d)

    figure;

    switch dim
        case 1  % 固定 E
            Z = squeeze(1./bounceT(idx, :, :));
            x_vec = Lambda1d;
            y_vec = -Pphi1d;          % 取反
            xlabel_str = '$\Lambda$';
            ylabel_str = '$P_{\varphi}$';
            title_str = sprintf('$E = %.3f$, $\\sigma = %+d$', E1d(idx), sign);
        case 2  % 固定 Pphi
            Z = squeeze(1./bounceT(:, idx, :));
            x_vec = Lambda1d;
            y_vec = E1d;
            xlabel_str = '$\Lambda$';
            ylabel_str = '$E$';
            title_str = sprintf('$P_{\\varphi} = %.3f$, $\\sigma = %+d$', -Pphi1d(idx), sign);
        case 3  % 固定 Lambda
            Z = squeeze(1./bounceT(:, :, idx));
            x_vec = -Pphi1d;          % 取反
            y_vec = E1d;
            xlabel_str = '$P_{\varphi}$';
            ylabel_str = '$E$';
            title_str = sprintf('$\\Lambda = %.3f$, $\\sigma = %+d$', Lambda1d(idx), sign);
        otherwise
            error('dim 必须是 1、2 或 3');
    end

    % 填充无穷大或NaN的点
    Z(isinf(Z)) = NaN;
    [ny, nx] = size(Z);
    offsets = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 1; 1, -1; 1, 0; 1, 1];
    dist = sqrt(sum(offsets.^2, 2));
    weights = 1 ./ dist;
    weights = weights / sum(weights);
    for i = 2:ny-1
        for j = 2:nx-1
            if isnan(Z(i,j))
                neighbor_vals = zeros(1,8);
                valid = true;
                for k = 1:8
                    ni = i + offsets(k,1);
                    nj = j + offsets(k,2);
                    val = Z(ni, nj);
                    if isnan(val)
                        valid = false;
                        break;
                    end
                    neighbor_vals(k) = val;
                end
                if valid
                    Z(i,j) = neighbor_vals * weights;
                end
            end
        end
    end

    [X, Y] = meshgrid(x_vec, y_vec);
    pcolor(X, Y, Z);
    shading interp;
    hold on;
    contour(X, Y, Z, 10, 'EdgeColor', 'k', 'LineWidth', 0.5);
    colormap(jet);
    colorbar;
    xlabel(xlabel_str, 'Interpreter', 'latex');
    ylabel(ylabel_str, 'Interpreter', 'latex');
    title(title_str, 'Interpreter', 'latex');

    cb = colorbar;
    cb.FontName = 'Times New Roman';
    cb.FontSize = 12;
    
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    
    xlabel(xlabel_str, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(ylabel_str, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(title_str, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

    set(gca, 'Color', 'white');
end