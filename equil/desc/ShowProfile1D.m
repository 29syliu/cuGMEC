%%

[ne,dne_dr] = sliding_poly_derivative(rhoSample, neSample, rho, 'ne', neWindow ,neOrder);
[Te,dTe_dr] = sliding_poly_derivative(rhoSample, TeSample, rho, 'Te', TeWindow ,TeOrder);


if IonType~=3
    [ni,dni_dr] = sliding_poly_derivative(rhoSample, niSample, rho, 'ni', niWindow ,niOrder);
    [Ti,dTi_dr] = sliding_poly_derivative(rhoSample, TiSample, rho, 'Ti', TiWindow ,TiOrder);
    if IonType==2
        [Pi,dPi_dr] = sliding_poly_derivative(rhoSample, PiSample, rho, 'Pi', PiWindow ,PiOrder);
    end
end


if AlphaType~=3
    [na,dna_dr] = sliding_poly_derivative(rhoSample, naSample, rho, 'na', naWindow ,naOrder);
    [Ta,dTa_dr] = sliding_poly_derivative(rhoSample, TaSample, rho, 'Ta', TaWindow ,TaOrder);
    if AlphaType==2
        [Pa,dPa_dr] = sliding_poly_derivative(rhoSample, PaSample, rho, 'Pa', PaWindow ,PaOrder);
    end
end


if BeamType~=3
    [nb,dnb_dr] = sliding_poly_derivative(rhoSample, nbSample, rho, 'nb', nbWindow ,nbOrder);
    [Tb,dTb_dr] = sliding_poly_derivative(rhoSample, TbSample, rho, 'Tb', TbWindow ,TbOrder);
    if BeamType==2
        [Pb,dPb_dr] = sliding_poly_derivative(rhoSample, PbSample, rho, 'Pb', PbWindow ,PbOrder);
    end
end

ni_c = niSample(1,1)*1e19;

if IonType~=3
    if min(ni)<0
        disp('Warning: min(ni) < 0');
    end
    if min(Ti)<0
        disp('Warning: min(Ti) < 0');
    end
    if IonType==2
        if min(Pi)<0
            disp('Warning: min(Pi) < 0');
        end
    end
end

if AlphaType~=3
    if min(na)<0
        disp('Warnang: min(na) < 0');
    end
    if min(Ta)<0
        disp('Warnang: min(Ta) < 0');
    end
    if AlphaType==2
        if min(Pa)<0
            disp('Warnang: min(Pa) < 0');
        end
    end
end

if BeamType~=3
    if min(nb)<0
        disp('Warnbng: min(nb) < 0');
    end
    if min(Tb)<0
        disp('Warnbng: min(Tb) < 0');
    end
    if BeamType==2
        if min(Pb)<0
            disp('Warnbng: min(Pb) < 0');
        end
    end
end

%%


function [y2, dy2dx2] = sliding_poly_derivative(x, y, x2, label, window_len, poly_deg, plot_fig1, plot_fig2)
% SLIDING_POLY_DERIVATIVE  基于滑动窗口多项式拟合计算任意点上的函数值及导数
%   [y2, dy2dx2] = SLIDING_POLY_DERIVATIVE(x, y, x2, label, window_len, poly_deg)
%   使用长度为 window_len 的滑动窗口（步长=1），每个窗口拟合 poly_deg 次多项式，
%   计算目标点 x2 上的插值函数值 y2 和导数值 dydx2。同时绘制两幅图：
%       图1：原始数据（蓝色）和拟合值（红色），图例分别为 'original label' 和 'fitted label'，标题为 label。
%       图2：原始数据点上的导数（蓝色）与 x2 上导数值连线（红色），图例分别为 'original range' 和 'simulated range'，标题为 fitted dlabel/drho。
%
%   插值策略：对于每个待求点，选择所有覆盖该点的窗口中，窗口中心最靠近该点的多项式进行插值。
%
%   输入：
%       x, y     - 原始数据向量（1×N 或 N×1）
%       x2       - 目标插值点，可以是向量或二维数组（Nx×Ny）
%       label    - 字符串，用于图例和标题
%       window_len - 窗口大小（整数，2 <= window_len <= length(x)）
%       poly_deg   - 多项式次数（整数，0 <= poly_deg < window_len）
%       plot_fig1  - 是否绘制图1，逻辑值，默认 true
%       plot_fig2  - 是否绘制图2，逻辑值，默认 true
%   输出：
%       y2       - 在 x2 上的插值函数值（与 x2 同形状）
%       dy2dx2   - 在 x2 上的导数值（与 x2 同形状）

    % 设置可选参数默认值
    if nargin < 7 || isempty(plot_fig1)
        plot_fig1 = true;
    end
    if nargin < 8 || isempty(plot_fig2)
        plot_fig2 = true;
    end

    % 统一 x, y 为行向量
    x = x(:).';
    y = y(:).';
    N = length(x);

    % 检查输入参数有效性
    if N < 2
        error('数据点数量至少为2');
    end
    if window_len < 2
        error('窗口长度至少为2');
    end
    if window_len > N
        error('窗口长度不能大于数据点数量');
    end
    if poly_deg < 0
        error('多项式次数不能小于0');
    end
    if poly_deg >= window_len
        error('多项式次数必须小于窗口长度');
    end

    % 存储所有窗口的多项式及其区间和中心
    polys = {};      % 多项式系数
    ranges = {};     % 窗口的 x 区间 [min, max]
    centers = [];    % 窗口中心

    % 滑动窗口：窗口大小 = window_len，步长 = 1
    for i = 1:N - window_len + 1
        idx = i:i+window_len-1;
        x_win = x(idx);
        y_win = y(idx);
        p = polyfit(x_win, y_win, poly_deg);
        polys{end+1} = p;
        ranges{end+1} = [min(x_win), max(x_win)];
        centers(end+1) = (min(x_win) + max(x_win)) / 2;
    end

    % ---------- 计算原始 x 点上的导数和拟合值（用于绘图） ----------
    dydx_orig = zeros(1, N);
    y_fit = zeros(1, N);
    for i = 1:N
        xi = x(i);
        valid_idx = [];
        for j = 1:length(polys)
            if xi >= ranges{j}(1) && xi <= ranges{j}(2)
                valid_idx = [valid_idx, j];
            end
        end
        if isempty(valid_idx)
            error('点 %g 没有被任何窗口覆盖，请检查数据分布', xi);
        end
        if length(valid_idx) == 1
            best_idx = valid_idx(1);
        else
            dist = abs(centers(valid_idx) - xi);
            [~, min_pos] = min(dist);
            best_idx = valid_idx(min_pos);
        end
        p = polys{best_idx};
        y_fit(i) = polyval(p, xi);
        if poly_deg == 0
            dydx_orig(i) = 0;
        else
            dp = polyder(p);
            dydx_orig(i) = polyval(dp, xi);
        end
    end

    % ---------- 处理 x2：支持二维数组，但仅沿第一维变化 ----------
    x2_orig_shape = size(x2);
    % 确保 x2 为行向量或二维矩阵（如果是一维，则 reshape 为列向量方便处理）
    if ndims(x2) == 2
        % 如果是二维，取第一列作为代表计算，然后扩展
        x2_col = x2(:, 1);        % 第一列，列向量
        % 对第一列进行插值计算
        [y2_col, dydx2_col] = interpolate_on_x2(x2_col);
        % 扩展到与 x2 相同的列数
        y2 = repmat(y2_col, 1, size(x2, 2));
        dy2dx2 = repmat(dydx2_col, 1, size(x2, 2));
    else
        % 一维或更高维，直接计算（但通常只有一维）
        x2_vec = x2(:).';
        [y2_vec, dydx2_vec] = interpolate_on_x2(x2_vec);
        % 恢复形状
        y2 = reshape(y2_vec, x2_orig_shape);
        dy2dx2 = reshape(dydx2_vec, x2_orig_shape);
    end

    % 嵌套函数：对给定 x2 向量进行插值（返回列向量）
    function [y2_vec, dydx2_vec] = interpolate_on_x2(x2_in)
        x2_in = x2_in(:);   % 确保列向量
        n2 = length(x2_in);
        y2_vec = zeros(n2, 1);
        dydx2_vec = zeros(n2, 1);
        for k = 1:n2
            xk = x2_in(k);
            valid_idx = [];
            for j = 1:length(polys)
                if xk >= ranges{j}(1) && xk <= ranges{j}(2)
                    valid_idx = [valid_idx, j];
                end
            end
            if isempty(valid_idx)
                error('x2 中的点 %g 不在任何窗口范围内', xk);
            end
            if length(valid_idx) == 1
                best_idx = valid_idx(1);
            else
                dist = abs(centers(valid_idx) - xk);
                [~, min_pos] = min(dist);
                best_idx = valid_idx(min_pos);
            end
            p = polys{best_idx};
            y2_vec(k) = polyval(p, xk);
            if poly_deg == 0
                dydx2_vec(k) = 0;
            else
                dp = polyder(p);
                dydx2_vec(k) = polyval(dp, xk);
            end
        end
    end

    % ---------- 绘图1：原始数据（蓝）与拟合值（红） ----------
    if plot_fig1
        figure;
        plot(x, y, 'b-', 'LineWidth', 1.5, 'DisplayName', ['original ', label]);
        hold on;
        plot(x, y_fit, 'r-', 'LineWidth', 1, 'DisplayName', ['fitted ', label]);
        hold off;
        legend('show');
        title(label);
        grid on;
    end

    % ---------- 绘图2：原始导数（蓝）与 x2 导数（红） ----------
    if plot_fig2
        figure;
        hold on;
        % 原始数据点上的导数（蓝色曲线）
        plot(x, dydx_orig, 'b-', 'LineWidth', 1.5, 'DisplayName', 'original range');
        % x2上的导数值连线（红色曲线）
        x2_flat = x2(:);
        dydx2_flat = dy2dx2(:);
        if ~isempty(x2_flat)
            [x2_sorted, sort_idx] = sort(x2_flat);
            dydx2_sorted = dydx2_flat(sort_idx);
            plot(x2_sorted, dydx2_sorted, 'r-', 'LineWidth', 2.5, 'DisplayName', 'simulated range');
        end
        hold off;
        legend('show');
        title(['fitted d', label, '/drho']);
        grid on;
    end
end