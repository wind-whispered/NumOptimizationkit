# NumOptimizationkit

**Mathematica 数值计算工具包**

一个统一的 Mathematica Paclet，将数值优化、方程求根、数值积分、常微分方程、偏微分方程、线性方程组、矩阵特征分析、多项式插值与逼近等算法整合为 **20 个公开入口函数**、**67 种算法**，通过 `Method` 选项切换。

---

## 目录

1. [安装方法](#安装方法)
2. [包结构](#包结构)
3. [函数与算法参考](#函数与算法参考)
   - [无约束极值求解](#无约束极值求解)
   - [方程求根](#方程求根)
   - [非线性方程组](#非线性方程组)
   - [数值积分](#数值积分)
   - [常微分方程初值问题](#常微分方程初值问题)
   - [边值问题](#边值问题)
   - [偏微分方程数值解](#偏微分方程数值解)
   - [线性方程组](#线性方程组)
   - [矩阵特征分析](#矩阵特征分析)
   - [多项式插值](#多项式插值)
   - [多项式逼近](#多项式逼近)
4. [快速示例](#快速示例)
5. [设计规范](#设计规范)

---

## 安装方法

```mathematica
(* 方式一：直接加载 *)
Get["D:\\path\\to\\NumOptimizationkit\\NumOptimizationkit.wl"]

(* 方式二：安装为 Paclet 后按名称加载 *)
PacletInstall["D:\\path\\to\\NumOptimizationkit"]
Needs["NumOptimizationkit`"]
```

---

## 包结构

```
NumOptimizationkit/
├── NumOptimizationkit.wl                    ← 公开声明 + 模块加载
├── PacletInfo.wl                            ← Paclet 元数据（版本 1.0.0）
├── Kernel/init.m                            ← Paclet 入口
├── Modules/
│   ├── Internal.wl                          ← 私有工具函数（所有模块共享）
│   ├── UnconstrainedMinimization.wl         ← FindMinimum1D, FindMinimumND
│   ├── NumericalQuadrature.wl               ← NumericalQuadrature
│   ├── EquationRootFinding.wl               ← FindEquationRoot
│   ├── NonlinearEquationSystems.wl          ← SolveNonlinearEquationSystem
│   ├── LinearEquationSystems.wl             ← SolveLinearEquationSystem + 矩阵分解
│   ├── MatrixEigenanalysis.wl               ← FindDominantEigenvalue, FindMatrixEigenvalues
│   ├── PolynomialInterpolation.wl           ← PolynomialInterpolate, HermiteInterpolate,
│   │                                           SplineInterpolate
│   ├── PolynomialApproximation.wl           ← ChebyshevApproximate, LeastSquaresApproximate,
│   │                                           PadeApproximate, RemezApproximate
│   ├── OrdinaryDifferentialEquations.wl     ← SolveInitialValueProblem
│   ├── BoundaryValueProblems.wl             ← SolveBoundaryValueProblem
│   └── PartialDifferentialEquations.wl      ← SolveParabolicPDE, SolveHyperbolicPDE
├── Tests/
│   ├── BasicTests.wl                        ← 自动化测试套件
│   └── 测试脚本.nb                           ← 可交互运行的示例 Notebook
└── Docs/
    ├── README.md                            ← 英文文档
    └── README_CN.md                         ← 本文件（中文文档）
```

### 私有工具函数（`Internal.wl`）

所有模块共享，位于 `NumOptimizationkit`Private`` 上下文：

| 符号 | 功能 |
|---|---|
| `numGrad[f, x]` | 向量点 x 处 f 的中心差分梯度 |
| `numDeriv[f, x, k]` | 标量函数 k 阶中心差分导数 |
| `numJacobian[F, x]` | 向量值函数 F 的数值 Jacobi 矩阵 |
| `numHessian[f, x]` | 标量函数 f 的数值 Hessian 矩阵 |
| `armijoLineSearch[f,x,g,d]` | Armijo 回溯线搜索 |
| `quadLineSearch[f,x,g,d]` | 二次插值线搜索 |
| `qrHouseholder[A]` | Householder QR 分解（线性代数、特征分析、逼近模块共用）|
| `laThomasAlgorithm[l,d,u,b]` | 三对角方程组追赶法（线性代数、边值问题模块共用）|
| `gaussLegendreNodesWeights[n]` | Newton 迭代法计算 Gauss-Legendre 节点和权重 |

---

## 函数与算法参考

### 无约束极值求解

#### `FindMinimum1D[f, {a, b}]`

在区间 [a, b] 上寻找单峰一元纯函数 f 的局部极小值。  
f 的调用形式为 `f[x]`，返回实数。

**返回值** `<| "Point" -> xp, "Value" -> fp, "Iterations" -> k |>`

**选项** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 200`

| Method | 算法 | 默认 |
|---|---|:---:|
| `"GoldenSection"` | 黄金搜索法 | ★ |
| `"Fibonacci"` | Fibonacci 搜索法 | |
| `"QuadraticInterpolation"` | 二次插值（抛物线）法 | |
| `"Newton"` | Newton 法（解 f'=0）| |

---

#### `FindMinimumND[f, x0]`

从初始点 x0（列表）出发，寻找多元纯函数 f 的局部极小值。  
f 的调用形式为 `f[{x1, x2, ...}]`，返回实数。

**返回值** `<| "Point" -> xp, "Value" -> fp, "Iterations" -> k, "Converged" -> True/False |>`

**选项** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`, `Gradient -> Automatic`

| Method | 算法 | 类型 | 默认 |
|---|---|---|:---:|
| `"BFGS"` | Broyden-Fletcher-Goldfarb-Shanno 变尺度法 | 拟 Newton | ★ |
| `"DFP"` | Davidon-Fletcher-Powell 变尺度法 | 拟 Newton | |
| `"GradientDescent"` | 最速下降法（梯度法）| 一阶 | |
| `"ConjugateGradientPR"` | 共轭梯度法（Polak-Ribière 公式）| 一阶 | |
| `"ConjugateGradientFR"` | 共轭梯度法（Fletcher-Reeves 公式）| 一阶 | |
| `"Newton"` | Newton 法（数值 Hessian）| 二阶 | |
| `"NelderMead"` | Nelder-Mead 单纯形法（无导数）| 直接搜索 | |
| `"SimulatedAnnealing"` | 模拟退火算法 | 全局 | |
| `"GeneticAlgorithm"` | 实数编码遗传算法 | 全局 | |

---

### 方程求根

#### `FindEquationRoot[f, {a, b}]`

寻找单变量纯函数 f 的一个根。

- **括号法**（Bisection、RegulaFalsi、Brent、Muller）：[a, b] 须含根，即 f(a)·f(b) < 0。
- **开放法**（Newton、Halley、FixedPoint、Steffensen、Aitken）：a 为初始迭代点，b 不使用。
- **割线法**（Secant）：a = x₀，b = x₁（两个不同的初始点）。

**返回值** `<| "Root" -> xp, "Residual" -> |f(xp)|, "Iterations" -> k, "Converged" -> True/False |>`

**选项** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 200`

| Method | 算法 | 收敛阶 | 类型 | 默认 |
|---|---|---|---|:---:|
| `"Bisection"` | 二分法 | 线性 | 括号法 | ★ |
| `"RegulaFalsi"` | 试位法 | 线性 | 括号法 | |
| `"Brent"` | Brent 法 | 超线性 | 括号法 | |
| `"Newton"` | Newton-Raphson 迭代法 | 二阶 | 开放法 | |
| `"Secant"` | 割线法 | ≈ 1.618 阶 | 两点法 | |
| `"Halley"` | Halley 法 | 三阶 | 开放法 | |
| `"Muller"` | 抛物线法（Müller）| ≈ 1.839 阶 | 括号法 | |
| `"FixedPoint"` | 不动点迭代法（f 为迭代函数 g(x)）| 线性 | 开放法 | |
| `"Steffensen"` | Steffensen 加速法 | 二阶 | 开放法 | |
| `"Aitken"` | Aitken Δ² 加速法 | 超线性 | 开放法 | |

---

### 非线性方程组

#### `SolveNonlinearEquationSystem[F, x0]`

从初始向量 x0 出发，求解 F(x) = 0（F : ℝⁿ → ℝⁿ）。  
F 调用形式为 `F[{x1,...,xn}]`，返回 n 个残差组成的列表。  
`Jacobian -> function` 可提供解析 Jacobi 矩阵（否则自动数值计算）。

**返回值** `<| "Solution" -> xp, "Residual" -> ‖F(xp)‖, "Iterations" -> k, "Converged" -> True/False |>`

**选项** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`, `Jacobian -> Automatic`

| Method | 算法 | 默认 |
|---|---|:---:|
| `"Newton"` | 多元 Newton 法（二阶收敛）| ★ |
| `"Broyden"` | Broyden 秩 1 拟 Newton 法 | |
| `"FixedPoint"` | 向量不动点迭代法 | |
| `"Continuation"` | 数值延拓法（同伦方法）| |

---

### 数值积分

#### `NumericalQuadrature[f, {x, a, b}]`

数值计算 ∫ₐᵇ f(x) dx。f 调用形式为 `f[x]`，返回实数。

**返回值**：实数（积分近似值）

**选项** — `Method`, `Points -> 5`（Gauss 型节点数）, `Intervals -> 100`（复合公式子区间数）, `Tolerance -> 1*^-6`

| Method | 算法 | 说明 | 默认 |
|---|---|---|:---:|
| `"GaussLegendre"` | Gauss-Legendre 求积公式 | `Points` 控制节点数 | ★ |
| `"Trapezoidal"` | 复合梯形公式 | `Intervals` 控制子区间数 | |
| `"Simpson"` | 复合 Simpson 公式 | | |
| `"AdaptiveSimpson"` | 自适应 Simpson 公式 | `Tolerance` 控制误差 | |
| `"Romberg"` | Romberg 公式（Richardson 外推）| | |
| `"GaussChebyshev"` | Gauss-Chebyshev I 型 | 积分 f(x)/√(1−x²) | |
| `"GaussHermite"` | Gauss-Hermite | 积分 f(x)·e^{-x²}，(−∞,∞) | |
| `"GaussLaguerre"` | Gauss-Laguerre | 积分 f(x)·e^{-x}，[0,∞) | |
| `"GaussLobatto"` | Gauss-Lobatto | 含两个端点 | |
| `"GaussRadau"` | Gauss-Radau | 含左端点 | |

---

### 常微分方程初值问题

#### `SolveInitialValueProblem[f, {t, t0, tEnd}, y0]`

求解一阶常微分方程初值问题 y'(t) = f(t, y)，y(t₀) = y₀。  
方程组形式：`f[t, {y1,...}]` 返回列表，`y0` 为初始向量。

**返回值** `<| "Grid" -> {t₀,...,tₙ}, "Solution" -> {y₀,...,yₙ} |>`

**选项** — `Method`, `StepSize -> Automatic`（= (tEnd−t₀)/100）

| Method | 算法 | 精度阶 | 默认 |
|---|---|---|:---:|
| `"RungeKutta4"` | 经典四阶 Runge-Kutta | 4 | ★ |
| `"Euler"` | 前向 Euler 法 | 1 | |
| `"ImplicitEuler"` | 隐式（后向）Euler 法 | 1 | |
| `"Heun"` | Heun 法（显式梯形）| 2 | |
| `"RungeKutta3"` | 三阶 Kutta 法 | 3 | |
| `"RungeKuttaFehlberg"` | RKF45 自适应步长控制 | 4/5 自适应 | |
| `"AdamsBashforth"` | Adams-Bashforth 四步显式法 | 4 | |
| `"AdamsPC4"` | Adams 四阶预测-校正法 | 4 | |
| `"HammingPC"` | Hamming 预测-校正法 | 4 | |

---

### 边值问题

#### `SolveBoundaryValueProblem[f, {x, a, b}, {ya, yb}]`
#### `SolveBoundaryValueProblem[f, {x, a, b}, {ya, yb}, n]`

求解二阶常微分方程边值问题 y''(x) = f(x, y, y')，y(a) = yₐ，y(b) = y_b。  
n 为内部节点数（默认 100）。

**返回值** `<| "Grid" -> {x₀,...,xₙ}, "Solution" -> {y₀,...,yₙ} |>`

**选项** — `Method`

| Method | 算法 | 默认 |
|---|---|:---:|
| `"Shooting"` | 打靶法（叠加两个 RK4 初值问题）| ★ |
| `"FiniteDifference"` | 中心差分法（三对角方程组，追赶法求解）| |

---

### 偏微分方程数值解

#### `SolveParabolicPDE[alpha2, {x,xL,xR,nx}, {t,t0,tEnd,nt}, ic, {bc0,bcL}]`

求解抛物型方程（热方程）u_t = alpha2·u_xx，x ∈ [xL,xR]，t ∈ [t0,tEnd]。  
共 nx+2 个空间节点，nt+1 个时间层。
- `ic[x]` — 初始条件 u(x, t₀)
- `bc0[t]`, `bcL[t]` — 左/右边界条件

**返回值** `<| "SpatialGrid" -> xs, "TimeGrid" -> ts, "Solution" -> 矩阵 |>`（矩阵`[[j,i]] = u(xᵢ, tⱼ)`）

| Method | 算法 | 稳定性条件 | 默认 |
|---|---|---|:---:|
| `"CrankNicolson"` | Crank-Nicolson 法 | 无条件稳定；O(Δt²,Δx²) | ★ |
| `"ForwardDifference"` | 显式向前差分法 | r = alpha2·Δt/Δx² ≤ 0.5 | |
| `"BackwardDifference"` | 隐式向后差分法 | 无条件稳定；O(Δt,Δx²) | |

#### `SolveHyperbolicPDE[c2, {x,xL,xR,nx}, {t,t0,tEnd,nt}, ic, vic, {bc0,bcL}]`

求解双曲型方程（波动方程）u_tt = c2·u_xx。  
`ic[x]` — 初始位移；`vic[x]` — 初始速度。  
CFL 稳定条件：√c2·Δt/Δx ≤ 1。

| Method | 算法 | 默认 |
|---|---|:---:|
| `"ExplicitDifference"` | 显式二阶有限差分法 | ★ |

---

### 线性方程组

#### `SolveLinearEquationSystem[A, b]`

求解线性方程组 A·x = b。返回解向量 x。

**选项** — `Method`, `Omega -> 1.5`（SOR 松弛因子）, `Tolerance -> 1*^-8`, `MaxIterations -> 1000`

| Method | 算法 | 类型 | 适用条件 | 默认 |
|---|---|---|---|:---:|
| `"GaussianEliminationPivot"` | Gauss 列主元消去法 | 直接法 | 通用 | ★ |
| `"GaussianElimination"` | Gauss 顺序消去法 | 直接法 | 通用 | |
| `"GaussJordan"` | Gauss-Jordan 全消去法 | 直接法 | 通用 | |
| `"LU"` | 列主元 LU 分解法 | 直接法 | 通用 | |
| `"Cholesky"` | Cholesky（LLᵀ）分解法 | 直接法 | 对称正定 | |
| `"LDLT"` | LDLᵀ 分解法 | 直接法 | 对称 | |
| `"Tridiagonal"` | 追赶法（Thomas 算法）| 直接法 | 三对角矩阵 | |
| `"Jacobi"` | Jacobi 迭代法 | 迭代法 | 对角占优 | |
| `"GaussSeidel"` | Gauss-Seidel 迭代法 | 迭代法 | 对角占优 | |
| `"SOR"` | 逐次超松弛法（ω = Omega）| 迭代法 | 对角占优 | |

#### 矩阵分解独立函数

| 函数 | 返回值 | 分解形式 |
|---|---|---|
| `LUDecompose[A]` | `{L, U, P}` | P·A = L·U（列主元置换）|
| `QRDecompose[A]` | `{Q, R}` | A = Q·R；`Method -> "Householder"` ★ 或 `"GramSchmidt"` |
| `CholeskyDecompose[A]` | `L` | A = L·Lᵀ（对称正定）|
| `LDLTDecompose[A]` | `{L, d}` | A = L·diag(d)·Lᵀ（对称）|

---

### 矩阵特征分析

#### `FindDominantEigenvalue[A]`

求矩阵 A 模最大的特征值（主特征值）及其特征向量。

**返回值** `<| "Eigenvalue" -> λ, "Eigenvector" -> v, "Iterations" -> k, "Converged" -> True/False |>`

**选项** — `Method`, `Shift -> 0`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`

| Method | 算法 | 默认 |
|---|---|:---:|
| `"PowerMethod"` | 乘幂法 | ★ |
| `"InversePowerMethod"` | 移位反幂法（求最靠近 Shift 的特征值）| |

#### `FindMatrixEigenvalues[A]`

求方阵 A 的全部特征值及特征向量。

**返回值** `<| "Eigenvalues" -> {λ₁,...}, "Eigenvectors" -> {v₁,...}, "Iterations" -> k |>`

**选项** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`

| Method | 算法 | 适用范围 | 默认 |
|---|---|---|:---:|
| `"QRIteration"` | 带原点位移的 QR 迭代法 | 一般实矩阵 | ★ |
| `"JacobiMethod"` | Jacobi 旋转法 | 仅实对称矩阵 | |

---

### 多项式插值

#### `PolynomialInterpolate[xs, ys, xq]`

通过数据节点 `{xs[[i]], ys[[i]]}` 构造插值多项式，在查询点 xq 处求值。

**选项** — `Method`

| Method | 算法 | 默认 |
|---|---|:---:|
| `"Newton"` | Newton 差商法（Horner 求值，效率高）| ★ |
| `"Lagrange"` | Lagrange 插值法 | |

#### `HermiteInterpolate[xs, ys, dys, xq]`

利用节点 xs 处的函数值 ys 和导数值 dys 构造 Hermite 插值多项式，在 xq 处求值。

#### `SplineInterpolate[xs, ys, xq]`

构造样条插值函数，在 xq 处求值。

**选项** — `Degree -> 3`, `BoundaryCondition -> "Natural"`

| Degree | 算法 | 默认 |
|---|---|:---:|
| `1` | 分段线性样条 | |
| `2` | 分段二次样条 | |
| `3` | 自然三次样条 | ★ |

---

### 多项式逼近

| 函数 | 算法 | 返回值 |
|---|---|---|
| `ChebyshevApproximate[f, {x,a,b}, n]` | [a,b] 上的 n 次 Chebyshev 展开 | 纯函数（直接调用 `fApprox[xq]`）|
| `LeastSquaresApproximate[xs, ys, basis]` | 法方程最小二乘拟合 | 系数向量 `{c₁,...,cₘ}` |
| `LeastSquaresApproximate[xs, ys, basis, "Orthogonal"]` | QR 分解最小二乘（更稳定）| 系数向量 |
| `PadeApproximate[f, x0, {p, q}]` | [p/q] Padé 有理逼近 | 纯函数 `R[xq]` |
| `RemezApproximate[f, {a,b}, n]` | Remez 交换算法（minimax 最佳一致逼近）| `<|"Coefficients", "Error", "Nodes"|>` |

---

## 快速示例

```mathematica
(* 加载程序包 *)
Get["D:\\path\\to\\NumOptimizationkit\\NumOptimizationkit.wl"]

(* 1. 黄金搜索法：x^2 - Cos[x] 在 [-2,2] 上的极小值 *)
FindMinimum1D[#^2 - Cos[#] &, {-2., 2.}]
(* <| "Point" -> 0.4502, "Value" -> -0.7969, "Iterations" -> 85 |> *)

(* 2. BFGS 求 Rosenbrock 函数最小值（全局最小值在 (1,1)）*)
rb = Function[x, (1 - x[[1]])^2 + 100 (x[[2]] - x[[1]]^2)^2];
FindMinimumND[rb, {-1., 0.5}, Method -> "BFGS"]

(* 3. Brent 法求 x^3 - x - 2 = 0 的根 *)
FindEquationRoot[Function[x, x^3 - x - 2], {1., 2.}, Method -> "Brent"]
(* <| "Root" -> 1.5214, "Residual" -> 0., "Iterations" -> 9, "Converged" -> True |> *)

(* 4. Newton 法求解非线性方程组 {x^2+y^2=5, x-y=1} *)
SolveNonlinearEquationSystem[
  Function[v, {v[[1]]^2 + v[[2]]^2 - 5, v[[1]] - v[[2]] - 1}], {1., 2.}]
(* <| "Solution" -> {2., 1.}, ... |> *)

(* 5. Gauss-Legendre 求积：∫₀^π sin(x)dx = 2 *)
NumericalQuadrature[Sin, {x, 0., Pi}, Method -> "GaussLegendre", Points -> 5]
(* 2.0（机器精度）*)

(* 6. 经典 RK4：求解 y'=-y，y(0)=1，精确解 Exp[-t] *)
sol = SolveInitialValueProblem[Function[{t, y}, -y], {t, 0., 5.}, 1., StepSize -> 0.05];
ListLinePlot[Transpose[{sol["Grid"], sol["Solution"]}]]

(* 7. Crank-Nicolson 热方程 *)
u = SolveParabolicPDE[1., {x, 0., 1., 49}, {t, 0., 0.1, 100},
      Function[x, Sin[Pi x]], {Function[t, 0.], Function[t, 0.]}];
MatrixPlot[u["Solution"], ColorFunction -> "TemperatureMap"]

(* 8. LU 分解 *)
{L, U, P} = LUDecompose[{{4.,2.,1.},{2.,5.,3.},{1.,3.,6.}}];
Norm[P . A - L . U]   (* 应为 0 *)

(* 9. Jacobi 法求实对称矩阵的全部特征值 *)
FindMatrixEigenvalues[{{4.,1.,2.},{1.,3.,0.},{2.,0.,5.}}, Method -> "JacobiMethod"]

(* 10. Chebyshev 逼近 Exp[x] 在 [0,2] 上 *)
fApprox = ChebyshevApproximate[Exp, {x, 0., 2.}, 10];
Plot[Exp[x] - fApprox[x], {x, 0., 2.}, PlotLabel -> "逼近误差"]
```

---

## 设计规范

### 函数输入约定

所有数值函数接受**纯函数**作为主要参数：

```mathematica
(* 单变量 *)
f = #^2 - Cos[#] &
f = Function[x, x^2 - Cos[x]]

(* 多变量（参数为向量列表）*)
F = Function[x, {x[[1]]^2 + x[[2]]^2 - 1, x[[1]] - x[[2]]}]

(* ODE 右端函数 f(t, y) *)
f = Function[{t, y}, -y]

(* 方程组 f(t, {y1,y2,...}) *)
f = Function[{t, y}, {y[[2]], -y[[1]]}]
```

### 返回值结构

| 函数类型 | 返回类型 | 主要键名 |
|---|---|---|
| 极值求解 | `Association` | `"Point"`, `"Value"`, `"Iterations"`, `"Converged"` |
| 方程求根 | `Association` | `"Root"`, `"Residual"`, `"Iterations"`, `"Converged"` |
| 非线性方程组 | `Association` | `"Solution"`, `"Residual"`, `"Iterations"`, `"Converged"` |
| 初值/边值问题 | `Association` | `"Grid"`, `"Solution"` |
| PDE | `Association` | `"SpatialGrid"`, `"TimeGrid"`, `"Solution"` |
| 线性方程组 | `List` | 直接返回解向量 |
| 矩阵分解 | `List` | 因子，如 `{L, U, P}` |
| 逼近函数（Chebyshev、Padé）| `Function` | 直接调用 `fApprox[xq]` |

### Method 选项

所有迭代函数均接受 `Method -> "方法名"` 选项。  
未知方法名触发 `::badmethod` 警告并自动回退到默认方法。  
★ 标记各函数的默认算法。

### 数值微分

未提供解析梯度/Jacobian 时，所有方法自动使用**中心差分**：

| 量 | 步长 |
|---|---|
| 梯度 / Jacobian | h = `$MachineEpsilon`^(1/3) |
| Hessian | h = `$MachineEpsilon`^(1/4) |

### 私有函数命名规范

所有内部实现函数采用**camelCase** 加模块前缀，**不使用下划线**连接词（Mathematica 中下划线是模式语法，会与受保护符号如 `Fibonacci`、`FixedPoint` 冲突）：

```
min1DGoldenSection   minNDQuasiNewton    quadGaussLegendre
rootBisection        rootFixedPointIter  nlsysNewton
laGaussElim          qrHouseholder       eigPowerMethod
ivpRungeKutta4       bvpShooting         pdeCrankNicolson
```
