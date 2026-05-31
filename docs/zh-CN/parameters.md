# cuGMEC 参数表

本文档对应 `src/cuGMEC_param.h`。修改这些参数后必须重新编译。

示范列只给出一种写法，不表示推荐值。

记号说明：

| 记号 | 代表 |
|---|---|
| `<Species>` | `Ion`、`Alpha`、`Beam` |
| `<MHDField>` | `Phi`、`A`、`dNe`、`dTe` |
| `<PICField>` | `dPi`、`dPa`、`dPb` |

## 数据类型和文件路径

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `mhdReal` | `float` / `double` | MHD 计算精度。 | MHD 精度应高于或等于 PIC 精度。 | `using mhdReal = double;` |
| `picReal` | `float` / `double` | PIC 计算精度。 | PIC 精度应低于或等于 MHD 精度。 | `using picReal = float;` |
| `inputDir` | 字符串路径 | 输入文件目录。 | 从头算时，仅需提供 `MHDCollocated` 以及可选的 `MHDStaggered` 文件。<br>继续算时，还需要提供 `MHDContinue` 和 `PICContinue` 文件，它们来自上一次任务 `outputDir/final`。 | `const std::string inputDir = "/path/to/input";` |
| `outputDir` | 字符串路径 | 输出文件目录。 |  | `const std::string outputDir = "/path/to/output";` |
| `MHDCollocated` | `.bin` 文件名 | collocated MHD 平衡文件。 |  | `const std::string MHDCollocated = "MHDCollocated.bin";` |
| `MHDStaggered` | `.bin` 文件名 | staggered MHD 平衡文件。 | 仅 `ifStaggered=trueType` 时需要。 | `const std::string MHDStaggered = "MHDStaggered.bin";` |
| `<Species>PhaseSpaceMapping` | `.bin` 文件名 | 相空间轨道映射文件。 | 仅 `ifOutputPhaseSpaceOrbit=trueType` 时需要。 | `const std::string IonPhaseSpaceMapping = "IonPhaseSpaceMapping.bin";` |

## 归一化参数

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `B0` | 正实数 | 磁场归一化基准。 |  | `const double B0 = 4.921751144619735;` |
| `L0` | 正实数 | 长度归一化基准。 |  | `const double L0 = 6.595629295925759;` |
| `VA0` | 正实数 | 速度归一化基准。 |  | `const double VA0 = 8.864164667700194e+06;` |
| `RHO0` | 正实数 | 径向区间左端。 |  | `const double RHO0 = 0.08;` |
| `RHO1` | 正实数 | 径向区间右端。 |  | `const double RHO1 = 0.90;` |
| `PSITMAX` | 正实数 | 最外层环向磁通（Wb/rad）。 |  | `const double PSITMAX = 18.868213504765762;` |

## 网格和并行规模

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `hostNums` | 正整数 | MPI 进程数。 |  | `const int hostNums = 4;` |
| `devNums` | 正整数 | 每个 MPI 进程使用的 GPU 数。 |  | `const int devNums = 4;` |
| `gridNx` | 正整数 | 径向网格数。 |  | `const int gridNx = 256;` |
| `gridNy` | 正整数 | 沿场线方向网格数。 | 需能被 `hostNums * devNums` 整除。 | `const int gridNy = 32;` |
| `gridNz` | 正整数 | 环向网格数。 |  | `const int gridNz = 96;` |
| `ppcNums` | 正整数 | 每种离子在每一个网格内的平均粒子数。 |  | `const int ppcNums = 256;` |

## 环向范围和初始扰动

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `tubes` | 正整数 | 模拟 `1/tubes` 个环向区域。 |  | `const int tubes = 6;` |
| `leftN` | 非负整数 | `1/tubes` 个环向区域里保留的环向模数下限。 | 实际环向模数范围为 `tubes*[leftN,rightN]`，间隔为 `tubes`。 | `const int leftN = 1;` |
| `rightN` | 非负整数 | `1/tubes` 个环向区域里保留的环向模数上限。 |  | `const int rightN = 6;` |
| `refinedTimes` | 正整数 | `selectNM_*` 极向滤波使用的沿场线方向插值细分倍数。 | 细分后网格数为 `gridNy * refinedTimes`。 | `const int refinedTimes = 32;` |
| `perturbLeftN` | 非负整数 | `1/tubes` 个环向区域里初始扰动环向模数下限。 | 实际环向模数范围为 `tubes*[perturbLeftN,perturbRightN]`，间隔为 `tubes`。 | `const int perturbLeftN = 1;` |
| `perturbRightN` | 非负整数 | `1/tubes` 个环向区域里初始扰动环向模数上限。 |  | `const int perturbRightN = 6;` |
| `perturbRadialIndex` | `1` 到 `gridNx` | 初始高斯扰动的径向峰值网格点。 |  | `const int perturbRadialIndex = 85;` |
| `perturbWidth` | 正实数 | 初始高斯扰动的径向宽度。 |  | `const mhdReal perturbWidth = 0.04;` |
| `perturbAmplitude` | 实数 | 初始高斯扰动的归一化 `Phi` 幅值。 |  | `const mhdReal perturbAmplitude = 2.5e-9;` |

## 滤波

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `ifFilterN_<MHDField>` | `trueType` / `falseType` | 是否对 MHD 场进行环向滤波。 |  | `using ifFilterN_Phi = trueType;` |
| `ifFilterN_dP` | `trueType` / `falseType` | 是否对 PIC 扰动压强进行环向滤波。 | 作用于 `dPi/dPa/dPb`。 | `using ifFilterN_dP = trueType;` |
| `removeN_<MHDField>` | `{}` 或环向模数列表 | 从 MHD 场中移除指定环向模数。 |  | `constexpr std::array<int, 1> removeN_Phi = {7};` |
| `removeN_dP` | `{}` 或环向模数列表 | 从 PIC 扰动压强中移除指定环向模数。 | 作用于 `dPi/dPa/dPb`。 | `constexpr std::array<int, 0> removeN_dP = {};` |
| `selectNM_<MHDField>` | `{{{N, leftM, rightM}, ...}}` | 对 MHD 场的指定环向模数进行极向滤波。 |  | `constexpr std::array<std::tuple<int, int, int>, 1> selectNM_Phi = {{{0, 0, 0}}};` |
| `selectNM_dP` | `{{{N, leftM, rightM}, ...}}` | 对 PIC 扰动压强的指定环向模数进行极向滤波。 | 作用于 `dPi/dPa/dPb`。 | `constexpr std::array<std::tuple<int, int, int>, 1> selectNM_dP = {{{0, 0, 1}}};` |

## MHD 物理开关

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `ifStaggered` | `trueType` / `falseType` | 是否包含交错网格。 | 需要 `MHDStaggered` 文件。 | `using ifStaggered = falseType;` |
| `ifNonlinearMHD` | `trueType` / `falseType` | 是否包含 MHD 非线性项。 |  | `using ifNonlinearMHD = trueType;` |
| `ifEparallel` | `trueType` / `falseType` | 是否包含平行电场项。 |  | `using ifEparallel = trueType;` |
| `ifFLRMHD` | `trueType` / `falseType` | 是否包含泊松方程中由背景热离子带来的 FLR 效应。 |  | `using ifFLRMHD = falseType;` |
| `ifQNeutrality` | `trueType` / `falseType` | 是否由准中性条件计算电子密度扰动。 | 开启时，`dNe` 由涡量项和启用物种的扰动密度共同决定。 | `using ifQNeutrality = trueType;` |
| `ifMaxwellStress` | `trueType` / `falseType` | 是否包含 Maxwell stress。 | 仅 `ifNonlinearMHD=trueType` 时有效。 | `using ifMaxwellStress = trueType;` |
| `ifReynoldsStress` | `trueType` / `falseType` | 是否包含 Reynolds stress。 | 仅 `ifNonlinearMHD=trueType` 时有效。 | `using ifReynoldsStress = trueType;` |
| `MaxwellStressCoef` | 实数 | Maxwell stress 系数。 | 仅 `ifMaxwellStress=trueType` 时有效。 | `const mhdReal MaxwellStressCoef = 1.0;` |
| `ReynoldsStressCoef` | 实数 | Reynolds stress 系数。 | 仅 `ifReynoldsStress=trueType` 时有效。 | `const mhdReal ReynoldsStressCoef = 1.0;` |

## 耗散

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `ifNablaPerp2<MHDField>` | `trueType` / `falseType` | 是否对 MHD 场施加垂直方向二阶耗散。 |  | `using ifNablaPerp2Phi = trueType;` |
| `perp2<MHDField>` | 正实数 | MHD 场的垂直二阶耗散系数。 |  | `const mhdReal perp2Phi = 1.0e-7;` |
| `ifNablaPara4<MHDField>` | `trueType` / `falseType` | 是否对 MHD 场施加平行方向四阶耗散。 |  | `using ifNablaPara4Phi = trueType;` |
| `para4<MHDField>` | 正实数 | MHD 场的平行四阶耗散系数。 |  | `const mhdReal para4Phi = 1.0e-6;` |
| `ifNablaPerp2dP` | `trueType` / `falseType` | 是否对 PIC 扰动压强施加垂直方向二阶耗散。 | 作用于 `dPi/dPa/dPb`。 | `using ifNablaPerp2dP = trueType;` |
| `perp2dP` | 正实数 | PIC 扰动压强的垂直二阶耗散系数。 | 作用于 `dPi/dPa/dPb`。 | `const mhdReal perp2dP = 1.0e-6;` |
| `ifNablaPara4dP` | `trueType` / `falseType` | 是否对 PIC 扰动压强施加平行方向四阶耗散。 | 作用于 `dPi/dPa/dPb`。 | `using ifNablaPara4dP = falseType;` |
| `para4dP` | 正实数 | PIC 扰动压强的平行四阶耗散系数。 | 作用于 `dPi/dPa/dPb`。 | `const mhdReal para4dP = 0.0;` |

## PIC 物理开关

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `ifFLRPIC` | `trueType` / `falseType` | 是否开启 PIC 部分的 FLR 效应。 |  | `using ifFLRPIC = trueType;` |
| `ifNonlinearPIC` | `trueType` / `falseType` | 是否开启 PIC 非线性项。 |  | `using ifNonlinearPIC = trueType;` |
| `gyroNums` | 正整数 | gyro-average 点数。 |  | `const int gyroNums = 4;` |

## PIC 物种参数

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `if<Species>` | `trueType` / `falseType` | 是否开启对应离子。 |  | `using ifIon = trueType;` |
| `<Species>Type` | `Maxwell` / `Slowing0` / `Slowing1` / `Slowing2` / `Slowing3` | 速度分布类型。 |  | `const disType IonType = Maxwell;` |
| `<Species>Space` | `spaceReal` / `spaceUniform` | 空间采样方式（marker 位置空间真实分布或者均匀分布）。 |  | `const spaceType IonSpace = spaceReal;` |
| `<Species>Velocity` | `velocityReal` / `velocityUniform` | 速度采样方式（marker 速度空间真实分布或者均匀分布）。 |  | `const velocityType IonVelocity = velocityUniform;` |
| `<Species>Mass` | 正实数 | 离子质量，按质子质量归一化。 |  | `const picReal IonMass = 2.5;` |
| `<Species>Char` | 正实数 | 离子电荷，按电子电荷归一化。 |  | `const picReal IonChar = 1.0;` |
| `<Species>Beta` | 正实数 | 离子磁轴处比压（`P/(B0^2/(2*mu0))`）。 |  | `const picReal IonBeta = 0.037793721898356;` |
| `<Species>Vmin` | 正实数 | 离子速度下限，按 `VA0` 归一化。 |  | `const picReal IonVmin = 0.0135;` |
| `<Species>Vmax` | 正实数 | 离子速度上限，按 `VA0` 归一化。 |  | `const picReal IonVmax = 0.54;` |
| `<Species>Vb` | 正实数 | slowing-down 分布的截断速度。 |  | `const picReal BeamVb = 0.1;` |
| `<Species>DeltaV` | 正实数 | slowing-down 分布的截断宽度。 |  | `const picReal BeamDeltaV = 0.1;` |
| `<Species>Lambda0` | 正实数 | 各向异性分布的 `Lambda` 中心。 |  | `const picReal BeamLambda0 = 0.4;` |
| `<Species>DeltaLambda2` | 正实数 | 各向异性分布的 `Lambda` 展宽平方。 |  | `const picReal BeamDeltaLambda2 = 1.0 / (4.5 * 4.5);` |

| 分布类型 | 公式 |
|---|---|
| `Maxwell` | `f ~ n * T^(-3/2) * exp[-m*v^2/(2*T)]` |
| `Slowing0` | `f ~ n / (v^3 + v_c^3)` |
| `Slowing1` | `f ~ n * [1 + erf((V_b - v)/DeltaV)] / (v^3 + v_c^3)` |
| `Slowing2` | `f ~ n * exp[-(Lambda - Lambda_0)^2 / DeltaLambda^2] / (v^3 + v_c^3)` |
| `Slowing3` | `f ~ n * [1 + erf((V_b - v)/DeltaV)] * exp[-(Lambda - Lambda_0)^2 / DeltaLambda^2] / (v^3 + v_c^3)` |

## 诊断开关

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `ifDiagAmplitude` | `trueType` / `falseType` | 是否诊断模式振幅。 |  | `using ifDiagAmplitude = trueType;` |
| `ifDiagFrequency` | `trueType` / `falseType` | 是否诊断模式频率。 |  | `using ifDiagFrequency = trueType;` |
| `ifDiagEparallel` | `trueType` / `falseType` | 是否诊断平行电场。 |  | `using ifDiagEparallel = trueType;` |
| `ifDiagDensity` | `trueType` / `falseType` | 是否诊断离子扰动密度。 |  | `using ifDiagDensity = trueType;` |
| `ifDiagDiffusivity` | `trueType` / `falseType` | 是否诊断离子扩散系数。 |  | `using ifDiagDiffusivity = trueType;` |
| `ifDiagZFDrive` | `trueType` / `falseType` | 是否诊断 zonal-flow 驱动源。 |  | `using ifDiagZFDrive = falseType;` |
| `ifCheckNAN` | `trueType` / `falseType` | 是否检查 NaN。 | 如果诊断到 NaN，程序会立刻停止并进行输出。 | `using ifCheckNAN = trueType;` |

## MHD 场输出开关

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `ifOutputPhi` | `trueType` / `falseType` | 是否输出多个时间点的 `Phi`。 |  | `using ifOutputPhi = trueType;` |
| `ifOutputA` | `trueType` / `falseType` | 是否输出多个时间点的 `A`。 |  | `using ifOutputA = falseType;` |
| `ifOutputdNe` | `trueType` / `falseType` | 是否输出多个时间点的 `dNe`。 |  | `using ifOutputdNe = falseType;` |
| `ifOutputdTe` | `trueType` / `falseType` | 是否输出多个时间点的 `dTe`。 |  | `using ifOutputdTe = falseType;` |
| `ifOutputdPi` | `trueType` / `falseType` | 是否输出多个时间点的 `dPi`。 |  | `using ifOutputdPi = falseType;` |
| `ifOutputdPa` | `trueType` / `falseType` | 是否输出多个时间点的 `dPa`。 |  | `using ifOutputdPa = falseType;` |
| `ifOutputdPb` | `trueType` / `falseType` | 是否输出多个时间点的 `dPb`。 |  | `using ifOutputdPb = falseType;` |

## 相空间输出参数

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `gridE` | 正整数 | `E` 方向 phase-space 网格数。 |  | `const int gridE = 96;` |
| `gridPphi` | 正整数 | `Pphi` 方向 phase-space 网格数。 |  | `const int gridPphi = 128;` |
| `gridLambda` | 正整数 | `Lambda` 方向 phase-space 网格数。 |  | `const int gridLambda = 48;` |
| `ppcPhase` | 正整数 | 每种离子在每一个 phase-space 网格内的平均粒子数。 |  | `const int ppcPhase = 2048;` |
| `ifOutputPhaseSpaceJacobian` | `trueType` / `falseType` | 是否输出 phase-space Jacobian。 | 初始阶段输出。 | `using ifOutputPhaseSpaceJacobian = trueType;` |
| `ifOutputPhaseSpaceOrbit` | `trueType` / `falseType` | 是否输出 phase-space 轨道频率。 | 需要 `<Species>PhaseSpaceMapping` 文件。 | `using ifOutputPhaseSpaceOrbit = falseType;` |
| `ifOutputPhaseSpaceF0` | `trueType` / `falseType` | 是否输出 phase-space 平衡分布函数。 | 初始阶段输出。 | `using ifOutputPhaseSpaceF0 = falseType;` |
| `ifOutputPhaseSpaceDeltaF` | `trueType` / `falseType` | 是否输出 phase-space 扰动分布函数。 |  | `using ifOutputPhaseSpaceDeltaF = falseType;` |
| `ifOutputPhaseSpacePower` | `trueType` / `falseType` | 是否输出 phase-space 波粒作用功率。 |  | `using ifOutputPhaseSpacePower = falseType;` |
| `gridVpara` | 正整数 | `vpara` 方向 pitch-space 网格数。 |  | `const int gridVpara = 128;` |
| `gridVperp` | 正整数 | `vperp` 方向 pitch-space 网格数。 |  | `const int gridVperp = 64;` |
| `ppcPitch` | 正整数 | 每种离子在每一个 pitch-space 网格内的平均粒子数。 |  | `const int ppcPitch = 2048;` |
| `ifOutputPitchSpaceJacobian` | `trueType` / `falseType` | 是否输出 pitch-space Jacobian。 | 初始阶段输出。 | `using ifOutputPitchSpaceJacobian = trueType;` |
| `ifOutputPitchSpaceF0` | `trueType` / `falseType` | 是否输出 pitch-space 平衡分布函数。 | 初始阶段输出。 | `using ifOutputPitchSpaceF0 = falseType;` |
| `ifOutputPitchSpaceDeltaF` | `trueType` / `falseType` | 是否输出 pitch-space 扰动分布函数。 |  | `using ifOutputPitchSpaceDeltaF = falseType;` |
| `ifOutputPitchSpacePower` | `trueType` / `falseType` | 是否输出 pitch-space 波粒作用功率。 |  | `using ifOutputPitchSpacePower = falseType;` |

## 时间推进

| 参数名 | 允许的值 | 含义 | 注意事项 | 示范 |
|---|---|---|---|---|
| `ifContinue` | `trueType` / `falseType` | 是否从继续算文件启动。 | `inputDir` 下需要有 `MHDContinue` 和 `PICContinue` 继续算文件。它们来自上一次任务 `outputDir/final`。 | `using ifContinue = falseType;` |
| `continueSteps` | 非负整数 | 继续算起始步数。 | 需与继续算文件名后缀一致。 | `const int continueSteps = 0;` |
| `dt` | 正实数 | MHD 时间步长。 |  | `const double dt = 0.02;` |
| `totalSteps` | 正整数 | 本次运行总 MHD 步数。 |  | `const int totalSteps = 20000;` |
| `ratioDt` | 正整数 | PIC 步长相对 MHD 步长的比例。 | PIC 每次推进使用 `dt * ratioDt`。 | `const int ratioDt = 1;` |
| `sortSteps` | 正整数 | 粒子排序间隔。 |  | `const int sortSteps = 25;` |
| `diagSteps` | 正整数 | 诊断采样间隔。 |  | `const int diagSteps = 1;` |
| `outputSteps` | 正整数 | 场和相空间输出间隔。 |  | `const int outputSteps = 2500;` |

