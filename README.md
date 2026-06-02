# NumOptimizationkit

**Numerical Computation Toolkit for Mathematica**

A unified Mathematica paclet that unifies numerical algorithms across optimization, root-finding, quadrature, ODE/PDE solvers, linear algebra, eigenanalysis, interpolation, and function approximation into **20 public entry-point functions** with **67 algorithm choices** selectable via the `Method` option.

---

## Table of Contents

1. [Installation](#installation)
2. [Package Structure](#package-structure)
3. [Function & Algorithm Reference](#function--algorithm-reference)
   - [Unconstrained Minimization](#unconstrained-minimization)
   - [Equation Root Finding](#equation-root-finding)
   - [Nonlinear Equation Systems](#nonlinear-equation-systems)
   - [Numerical Quadrature](#numerical-quadrature)
   - [ODE Initial Value Problems](#ode-initial-value-problems)
   - [Boundary Value Problems](#boundary-value-problems)
   - [Partial Differential Equations](#partial-differential-equations)
   - [Linear Equation Systems](#linear-equation-systems)
   - [Matrix Eigenanalysis](#matrix-eigenanalysis)
   - [Polynomial Interpolation](#polynomial-interpolation)
   - [Polynomial Approximation](#polynomial-approximation)
4. [Quick Examples](#quick-examples)
5. [Design Conventions](#design-conventions)

---

## Installation

```mathematica
(* Load directly by path *)
Get["D:\\path\\to\\NumOptimizationkit\\NumOptimizationkit.wl"]

(* Or install as a paclet, then load by name *)
PacletInstall["D:\\path\\to\\NumOptimizationkit"]
Needs["NumOptimizationkit`"]
```

---

## Package Structure

```
NumOptimizationkit/
├── NumOptimizationkit.wl                    ← public declarations + module loader
├── PacletInfo.wl                            ← paclet metadata (v1.0.0)
├── Kernel/init.m                            ← paclet entry point
├── Modules/
│   ├── Internal.wl                          ← shared private utilities
│   ├── UnconstrainedMinimization.wl         ← FindMinimum1D, FindMinimumND
│   ├── NumericalQuadrature.wl               ← NumericalQuadrature
│   ├── EquationRootFinding.wl               ← FindEquationRoot
│   ├── NonlinearEquationSystems.wl          ← SolveNonlinearEquationSystem
│   ├── LinearEquationSystems.wl             ← SolveLinearEquationSystem + decompositions
│   ├── MatrixEigenanalysis.wl               ← FindDominantEigenvalue, FindMatrixEigenvalues
│   ├── PolynomialInterpolation.wl           ← PolynomialInterpolate, HermiteInterpolate,
│   │                                           SplineInterpolate
│   ├── PolynomialApproximation.wl           ← ChebyshevApproximate, LeastSquaresApproximate,
│   │                                           PadeApproximate, RemezApproximate
│   ├── OrdinaryDifferentialEquations.wl     ← SolveInitialValueProblem
│   ├── BoundaryValueProblems.wl             ← SolveBoundaryValueProblem
│   └── PartialDifferentialEquations.wl      ← SolveParabolicPDE, SolveHyperbolicPDE
├── Tests/
│   ├── BasicTests.wl                        ← automated test suite
│   └── 测试脚本.nb                           ← interactive demo notebook (English content)
└── Docs/
    ├── README.md                            ← this file
    └── README_CN.md                         ← Chinese documentation
```

### Shared Private Utilities (`Internal.wl`)

All modules share the following private helpers (in `NumOptimizationkit`Private`` context):

| Symbol | Description |
|---|---|
| `numGrad[f, x]` | Central-difference gradient of `f` at vector `x` |
| `numDeriv[f, x, k]` | k-th order central-difference derivative (scalar) |
| `numJacobian[F, x]` | Numerical Jacobian matrix of vector-valued `F` |
| `numHessian[f, x]` | Numerical Hessian matrix of scalar `f` |
| `armijoLineSearch[f,x,g,d]` | Armijo backtracking line search |
| `quadLineSearch[f,x,g,d]` | Quadratic-interpolation line search |
| `qrHouseholder[A]` | Householder QR factorisation (used by linear algebra, eigenanalysis, and approximation modules) |
| `laThomasAlgorithm[l,d,u,b]` | Thomas algorithm for tridiagonal systems (used by linear algebra and BVP modules) |
| `gaussLegendreNodesWeights[n]` | Newton-iteration computation of Gauss-Legendre nodes and weights |

---

## Function & Algorithm Reference

### Unconstrained Minimization

#### `FindMinimum1D[f, {a, b}]`

Finds a local minimum of a univariate pure function `f` on `[a, b]`.  
`f` is called as `f[x]` and must return a real number.

**Returns** `<| "Point" -> xp, "Value" -> fp, "Iterations" -> k |>`

**Options** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 200`

| Method | Algorithm | Default |
|---|---|:---:|
| `"GoldenSection"` | Golden-section search | ★ |
| `"Fibonacci"` | Fibonacci search | |
| `"QuadraticInterpolation"` | Quadratic (parabolic) interpolation search | |
| `"Newton"` | Newton's method applied to f' = 0 | |

#### `FindMinimumND[f, x0]`

Finds a local minimum of a multivariate pure function starting from `x0` (a list).  
`f` is called as `f[{x1, x2, ...}]` and must return a real number.

**Returns** `<| "Point" -> xp, "Value" -> fp, "Iterations" -> k, "Converged" -> True/False |>`

**Options** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`, `Gradient -> Automatic`

| Method | Algorithm | Type | Default |
|---|---|---|:---:|
| `"BFGS"` | Broyden-Fletcher-Goldfarb-Shanno | Quasi-Newton | ★ |
| `"DFP"` | Davidon-Fletcher-Powell | Quasi-Newton | |
| `"GradientDescent"` | Steepest descent with quadratic line search | First-order | |
| `"ConjugateGradientPR"` | Conjugate gradient, Polak-Ribière formula | First-order | |
| `"ConjugateGradientFR"` | Conjugate gradient, Fletcher-Reeves formula | First-order | |
| `"Newton"` | Newton's method with numerical Hessian | Second-order | |
| `"NelderMead"` | Nelder-Mead simplex (derivative-free) | Direct search | |
| `"SimulatedAnnealing"` | Simulated annealing | Global | |
| `"GeneticAlgorithm"` | Real-coded genetic algorithm | Global | |

---

### Equation Root Finding

#### `FindEquationRoot[f, {a, b}]`

Finds a root of a univariate pure function `f`.

- **Bracket methods** (Bisection, RegulaFalsi, Brent, Muller): `[a, b]` must bracket a root, i.e. `f(a)*f(b) < 0`.
- **Open methods** (Newton, Halley, FixedPoint, Steffensen, Aitken): `a` is the starting iterate; `b` is ignored.
- **Secant**: `a` = x₀, `b` = x₁ (two distinct starting points).

**Returns** `<| "Root" -> xp, "Residual" -> |f(xp)|, "Iterations" -> k, "Converged" -> True/False |>`

**Options** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 200`

| Method | Algorithm | Order | Type | Default |
|---|---|---|---|:---:|
| `"Bisection"` | Bisection method | Linear | Bracket | ★ |
| `"RegulaFalsi"` | Regula falsi (false position) | Linear | Bracket | |
| `"Brent"` | Brent's method | Superlinear | Bracket | |
| `"Newton"` | Newton-Raphson iteration | Quadratic | Open | |
| `"Secant"` | Secant method | ≈ 1.618 | Two-point | |
| `"Halley"` | Halley's method | Cubic | Open | |
| `"Muller"` | Muller's parabolic method | ≈ 1.839 | Bracket | |
| `"FixedPoint"` | Fixed-point iteration (`f` is `g(x)`) | Linear | Open | |
| `"Steffensen"` | Steffensen acceleration | Quadratic | Open | |
| `"Aitken"` | Aitken Δ² acceleration | Superlinear | Open | |

---

### Nonlinear Equation Systems

#### `SolveNonlinearEquationSystem[F, x0]`

Solves the system F(x) = 0 for F : ℝⁿ → ℝⁿ, starting from initial vector `x0`.  
`F` is called as `F[{x1,...,xn}]` and must return a list of n residuals.

**Returns** `<| "Solution" -> xp, "Residual" -> ‖F(xp)‖, "Iterations" -> k, "Converged" -> True/False |>`

**Options** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`, `Jacobian -> Automatic`

| Method | Algorithm | Default |
|---|---|:---:|
| `"Newton"` | Multivariate Newton's method (quadratic convergence) | ★ |
| `"Broyden"` | Broyden rank-1 quasi-Newton update | |
| `"FixedPoint"` | Vector fixed-point iteration | |
| `"Continuation"` | Numerical continuation / homotopy method | |

---

### Numerical Quadrature

#### `NumericalQuadrature[f, {x, a, b}]`

Numerically evaluates `Integral[f[x], {x, a, b}]`.  
`f` is called as `f[x]` and must return a real number.

**Returns** a real number.

**Options** — `Method`, `Points -> 5`, `Intervals -> 100`, `Tolerance -> 1*^-6`

| Method | Algorithm | Notes | Default |
|---|---|---|:---:|
| `"GaussLegendre"` | Gauss-Legendre quadrature | `Points` controls nodes | ★ |
| `"Trapezoidal"` | Composite trapezoidal rule | `Intervals` controls subintervals | |
| `"Simpson"` | Composite Simpson's rule | | |
| `"AdaptiveSimpson"` | Recursive adaptive Simpson | Error-controlled via `Tolerance` | |
| `"Romberg"` | Romberg's method (Richardson extrapolation) | | |
| `"GaussChebyshev"` | Gauss-Chebyshev type I | Integrates f(x)/√(1−x²) | |
| `"GaussHermite"` | Gauss-Hermite | Integrates f(x)·Exp[−x²] on (−∞,∞) | |
| `"GaussLaguerre"` | Gauss-Laguerre | Integrates f(x)·Exp[−x] on [0,∞) | |
| `"GaussLobatto"` | Gauss-Lobatto | Includes both endpoints | |
| `"GaussRadau"` | Gauss-Radau | Includes left endpoint | |

---

### ODE Initial Value Problems

#### `SolveInitialValueProblem[f, {t, t0, tEnd}, y0]`

Solves y'(t) = f(t, y),  y(t₀) = y₀.  
For systems: `f[t, {y1,...}]` returns a list; `y0` is a list.

**Returns** `<| "Grid" -> {t₀,...,tₙ}, "Solution" -> {y₀,...,yₙ} |>`

**Options** — `Method`, `StepSize -> Automatic` (= (tEnd−t0)/100)

| Method | Algorithm | Order | Default |
|---|---|---|:---:|
| `"RungeKutta4"` | Classical 4th-order Runge-Kutta | 4 | ★ |
| `"Euler"` | Forward Euler | 1 | |
| `"ImplicitEuler"` | Implicit (backward) Euler | 1 | |
| `"Heun"` | Heun's method (explicit trapezoidal) | 2 | |
| `"RungeKutta3"` | Classical 3rd-order Kutta method | 3 | |
| `"RungeKuttaFehlberg"` | RKF45 adaptive step-size control | 4/5 adaptive | |
| `"AdamsBashforth"` | 4-step Adams-Bashforth explicit | 4 | |
| `"AdamsPC4"` | Adams 4th-order predictor-corrector | 4 | |
| `"HammingPC"` | Hamming predictor-corrector | 4 | |

---

### Boundary Value Problems

#### `SolveBoundaryValueProblem[f, {x, a, b}, {ya, yb}]`
#### `SolveBoundaryValueProblem[f, {x, a, b}, {ya, yb}, n]`

Solves y''(x) = f(x, y, y'),  y(a) = yₐ,  y(b) = y_b  
using `n` interior grid points (default `n = 100`).

**Returns** `<| "Grid" -> {x₀,...,xₙ}, "Solution" -> {y₀,...,yₙ} |>`

**Options** — `Method`

| Method | Algorithm | Default |
|---|---|:---:|
| `"Shooting"` | Linear shooting method (superposition of two IVPs via RK4) | ★ |
| `"FiniteDifference"` | Central finite differences, Thomas algorithm | |

---

### Partial Differential Equations

#### `SolveParabolicPDE[alpha2, {x,xL,xR,nx}, {t,t0,tEnd,nt}, ic, {bc0,bcL}]`

Solves the heat equation  u_t = alpha2 · u_xx  on [xL,xR] × [t0,tEnd].  
`nx+2` spatial nodes; `nt+1` time levels.

- `ic[x]` — initial condition u(x, t₀)
- `bc0[t]`, `bcL[t]` — left/right boundary conditions

**Returns** `<| "SpatialGrid" -> xs, "TimeGrid" -> ts, "Solution" -> matrix |>`  
where `matrix[[j, i]] = u(xᵢ, tⱼ)`.

| Method | Algorithm | Stability | Default |
|---|---|---|:---:|
| `"CrankNicolson"` | Crank-Nicolson | Unconditional; O(Δt²,Δx²) | ★ |
| `"ForwardDifference"` | Explicit forward difference | r = alpha2·Δt/Δx² ≤ 0.5 | |
| `"BackwardDifference"` | Implicit backward difference | Unconditional; O(Δt,Δx²) | |

#### `SolveHyperbolicPDE[c2, {x,xL,xR,nx}, {t,t0,tEnd,nt}, ic, vic, {bc0,bcL}]`

Solves the wave equation  u_tt = c2 · u_xx.  
`ic[x]` — initial displacement; `vic[x]` — initial velocity.  
CFL stability condition: √c2 · Δt/Δx ≤ 1.

| Method | Algorithm | Default |
|---|---|:---:|
| `"ExplicitDifference"` | Explicit second-order finite difference | ★ |

---

### Linear Equation Systems

#### `SolveLinearEquationSystem[A, b]`

Solves A·x = b.  Returns solution vector x.

**Options** — `Method`, `Omega -> 1.5` (SOR relaxation factor), `Tolerance -> 1*^-8`, `MaxIterations -> 1000`

| Method | Algorithm | Type | Requirement | Default |
|---|---|---|---|:---:|
| `"GaussianEliminationPivot"` | Gaussian elimination with partial pivoting | Direct | General | ★ |
| `"GaussianElimination"` | Gaussian elimination without pivoting | Direct | General | |
| `"GaussJordan"` | Gauss-Jordan full elimination | Direct | General | |
| `"LU"` | LU factorisation with partial pivoting | Direct | General | |
| `"Cholesky"` | Cholesky (LL^T) factorisation | Direct | Sym. pos. def. | |
| `"LDLT"` | LDL^T factorisation | Direct | Symmetric | |
| `"Tridiagonal"` | Thomas algorithm | Direct | Tridiagonal | |
| `"Jacobi"` | Jacobi iteration | Iterative | Diag. dominant | |
| `"GaussSeidel"` | Gauss-Seidel iteration | Iterative | Diag. dominant | |
| `"SOR"` | Successive over-relaxation (ω = Omega) | Iterative | Diag. dominant | |

#### Matrix Decompositions

Standalone factorisation functions (also used internally):

| Function | Returns | Factorisation |
|---|---|---|
| `LUDecompose[A]` | `{L, U, P}` | P·A = L·U, partial pivoting |
| `QRDecompose[A]` | `{Q, R}` | A = Q·R; `Method -> "Householder"` ★ or `"GramSchmidt"` |
| `CholeskyDecompose[A]` | `L` | A = L·Lᵀ (symmetric positive-definite) |
| `LDLTDecompose[A]` | `{L, d}` | A = L·diag(d)·Lᵀ (symmetric) |

---

### Matrix Eigenanalysis

#### `FindDominantEigenvalue[A]`

Finds the eigenvalue with the largest absolute value (dominant eigenvalue) and its eigenvector.

**Returns** `<| "Eigenvalue" -> λ, "Eigenvector" -> v, "Iterations" -> k, "Converged" -> True/False |>`

**Options** — `Method`, `Shift -> 0`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`

| Method | Algorithm | Default |
|---|---|:---:|
| `"PowerMethod"` | Power iteration | ★ |
| `"InversePowerMethod"` | Shifted inverse power iteration (finds eigenvalue nearest `Shift`) | |

#### `FindMatrixEigenvalues[A]`

Finds all eigenvalues (and eigenvectors) of a square matrix.

**Returns** `<| "Eigenvalues" -> {λ₁,...}, "Eigenvectors" -> {v₁,...}, "Iterations" -> k |>`

**Options** — `Method`, `Tolerance -> 1*^-8`, `MaxIterations -> 500`

| Method | Algorithm | Restriction | Default |
|---|---|---|:---:|
| `"QRIteration"` | QR iteration with origin shift | General real matrices | ★ |
| `"JacobiMethod"` | Jacobi rotation method | Symmetric matrices only | |

---

### Polynomial Interpolation

#### `PolynomialInterpolate[xs, ys, xq]`

Evaluates the interpolating polynomial through data nodes `{xs[[i]], ys[[i]]}` at query point `xq`.

**Options** — `Method`

| Method | Algorithm | Default |
|---|---|:---:|
| `"Newton"` | Newton divided-difference form, Horner evaluation | ★ |
| `"Lagrange"` | Lagrange interpolation | |

#### `HermiteInterpolate[xs, ys, dys, xq]`

Evaluates the Hermite interpolating polynomial matching both function values `ys` and derivative values `dys` at nodes `xs`, at query point `xq`.

#### `SplineInterpolate[xs, ys, xq]`

Evaluates a piecewise spline interpolant at query point `xq`.

**Options** — `Degree -> 3`, `BoundaryCondition -> "Natural"`

| Degree | Algorithm | Default |
|---|---|:---:|
| `1` | Piecewise linear spline | |
| `2` | Piecewise quadratic spline | |
| `3` | Natural cubic spline | ★ |

---

### Polynomial Approximation

| Function | Algorithm | Returns |
|---|---|---|
| `ChebyshevApproximate[f, {x,a,b}, n]` | Degree-n Chebyshev expansion on [a,b] | Pure function (call as `fApprox[xq]`) |
| `LeastSquaresApproximate[xs, ys, basis]` | Normal equations least-squares fit | Coefficient vector `{c₁,...,cₘ}` |
| `LeastSquaresApproximate[xs, ys, basis, "Orthogonal"]` | QR-based least-squares (more stable) | Coefficient vector |
| `PadeApproximate[f, x0, {p, q}]` | [p/q] Padé rational approximant at x0 | Pure function `R[xq]` |
| `RemezApproximate[f, {a,b}, n]` | Remez exchange algorithm (minimax) | `<|"Coefficients", "Error", "Nodes"|>` |

---

## Quick Examples

```mathematica
(* Load the package *)
Get["D:\\path\\to\\NumOptimizationkit\\NumOptimizationkit.wl"]

(* 1. Golden-section search: minimum of x^2 - Cos[x] on [-2, 2] *)
FindMinimum1D[#^2 - Cos[#] &, {-2., 2.}]
(* <| "Point" -> 0.4502, "Value" -> -0.7969, "Iterations" -> 85 |> *)

(* 2. BFGS on Rosenbrock function, global minimum at (1, 1) *)
rb = Function[x, (1 - x[[1]])^2 + 100 (x[[2]] - x[[1]]^2)^2];
FindMinimumND[rb, {-1., 0.5}, Method -> "BFGS"]

(* 3. Brent root-finding for x^3 - x - 2 = 0 *)
FindEquationRoot[Function[x, x^3 - x - 2], {1., 2.}, Method -> "Brent"]
(* <| "Root" -> 1.5214, "Residual" -> 0., "Iterations" -> 9, "Converged" -> True |> *)

(* 4. Newton's method for the nonlinear system x^2+y^2=5, x-y=1 *)
SolveNonlinearEquationSystem[
  Function[v, {v[[1]]^2 + v[[2]]^2 - 5, v[[1]] - v[[2]] - 1}], {1., 2.}]
(* <| "Solution" -> {2., 1.}, "Residual" -> 0., ... |> *)

(* 5. Gauss-Legendre: integral of Sin from 0 to Pi = 2 *)
NumericalQuadrature[Sin, {x, 0., Pi}, Method -> "GaussLegendre", Points -> 5]
(* 2.0 (to machine precision) *)

(* 6. RK4: solve y' = -y, y(0) = 1; exact solution Exp[-t] *)
sol = SolveInitialValueProblem[Function[{t, y}, -y], {t, 0., 5.}, 1., StepSize -> 0.05];
ListLinePlot[Transpose[{sol["Grid"], sol["Solution"]}]]

(* 7. Crank-Nicolson heat equation *)
u = SolveParabolicPDE[1., {x, 0., 1., 49}, {t, 0., 0.1, 100},
      Function[x, Sin[Pi x]], {Function[t, 0.], Function[t, 0.]}];
MatrixPlot[u["Solution"], ColorFunction -> "TemperatureMap"]

(* 8. LU decomposition *)
{L, U, P} = LUDecompose[{{4.,2.,1.},{2.,5.,3.},{1.,3.,6.}}];
Norm[P . A - L . U]   (* should be 0 *)

(* 9. All eigenvalues via Jacobi method *)
FindMatrixEigenvalues[{{4.,1.,2.},{1.,3.,0.},{2.,0.,5.}}, Method -> "JacobiMethod"]

(* 10. Chebyshev approximation of Exp on [0,2] *)
fApprox = ChebyshevApproximate[Exp, {x, 0., 2.}, 10];
Plot[Exp[x] - fApprox[x], {x, 0., 2.}, PlotLabel -> "Approximation error"]
```

---

## Design Conventions

### Input

All numerical functions accept **pure functions** as their primary argument:

```mathematica
f = #^2 - Cos[#] &                               (* univariate *)
F = Function[x, {x[[1]]^2 + x[[2]]^2 - 1, ...}] (* vector-valued, R^n -> R^n *)
f = Function[{t, y}, -y]                          (* ODE: f(t, y) *)
f = Function[{t, y}, {y[[2]], -y[[1]]}]           (* ODE system *)
```

### Return Values

| Function class | Return type | Keys |
|---|---|---|
| Minimisation | `Association` | `"Point"`, `"Value"`, `"Iterations"`, `"Converged"` |
| Root finding | `Association` | `"Root"`, `"Residual"`, `"Iterations"`, `"Converged"` |
| Nonlinear system | `Association` | `"Solution"`, `"Residual"`, `"Iterations"`, `"Converged"` |
| IVP / BVP | `Association` | `"Grid"`, `"Solution"` |
| PDE | `Association` | `"SpatialGrid"`, `"TimeGrid"`, `"Solution"` |
| Linear system | `List` | solution vector directly |
| Decomposition | `List` | factors, e.g. `{L, U, P}` |
| Approximation (Chebyshev, Padé) | `Function` | call as `fApprox[xq]` |

### Method Option

All iterative functions accept `Method -> "MethodName"` as a string.  
Unknown method names trigger a `::badmethod` warning and fall back to the default.  
★ marks the default method for each function.

### Numerical Differentiation

When `Gradient -> Automatic` or `Jacobian -> Automatic`, all methods use **central finite differences**:

| Quantity | Step size |
|---|---|
| Gradient / Jacobian | h = `$MachineEpsilon`^(1/3) |
| Hessian | h = `$MachineEpsilon`^(1/4) |

### Private Naming Convention

All internal implementation functions use **camelCase** with a descriptive module prefix — no underscores in the middle of names (underscores are pattern syntax in Mathematica and would conflict with protected symbols such as `Fibonacci` and `FixedPoint`):

```
min1DGoldenSection   minNDQuasiNewton    quadGaussLegendre
rootBisection        rootFixedPointIter  nlsysNewton
laGaussElim          qrHouseholder       eigPowerMethod
ivpRungeKutta4       bvpShooting         pdeCrankNicolson
```
