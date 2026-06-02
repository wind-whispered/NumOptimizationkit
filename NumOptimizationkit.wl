(* ::Package:: *)

(* NumOptimizationkit -- Numerical Computation Toolkit for Mathematica
   Provides unified entry-point functions for unconstrained optimization,
   equation root-finding, nonlinear systems, numerical quadrature,
   ODE/PDE solvers (including stiff), 2D PDE solvers, linear algebra
   (dense and sparse), eigenanalysis, polynomial interpolation,
   scattered-data interpolation and function approximation.
   The Method option selects the specific algorithm.

   Batch-3 additions
     Stiff ODE    : TrapezoidalRule, BDF2, BDF4
     Constrained  : FindConstrainedMinimum
     2D PDE       : SolveEllipticPDE, SolveParabolicPDE2D, SolveHyperbolicPDE2D
     Sparse linear: ConjugateGradient, SparseDirectLU
     Scattered    : ScatteredInterpolate
*)

BeginPackage["NumOptimizationkit`"]

(* ============================================================
   Unconstrained Minimization
   ============================================================ *)

FindMinimum1D::usage =
"FindMinimum1D[f, {a, b}] locates a local minimum of the univariate pure \
function f on [a, b].  f is called as f[x] and must return a real number.
Returns <|\"Point\"->xp, \"Value\"->fp, \"Iterations\"->k|>.

Options
  Method        -> \"GoldenSection\" (default) | \"Fibonacci\" |
                   \"QuadraticInterpolation\" | \"Newton\"
  Tolerance     -> 1*^-8
  MaxIterations -> 200

Algorithms
  \"GoldenSection\"          golden-section search
  \"Fibonacci\"              Fibonacci search
  \"QuadraticInterpolation\" quadratic (parabolic) interpolation search
  \"Newton\"                 Newton's method applied to f' = 0"

FindMinimumND::usage =
"FindMinimumND[f, x0] locates a local minimum of the multivariate pure \
function f starting from initial point x0 (a list).
f is called as f[{x1, x2, ...}] and must return a real number.
Returns <|\"Point\"->xp, \"Value\"->fp, \"Iterations\"->k, \"Converged\"->True/False|>.

Options
  Method             -> \"BFGS\" (default) | \"DFP\" | \"GradientDescent\" |
                        \"ConjugateGradientPR\" | \"ConjugateGradientFR\" |
                        \"Newton\" | \"NelderMead\" | \"SimulatedAnnealing\" |
                        \"GeneticAlgorithm\"
  Tolerance          -> 1*^-8
  MaxIterations      -> 500
  Gradient           -> Automatic | gradient-function
  ConvergenceHistory -> False  (when True, result includes \"ValueHistory\" and \"PointHistory\")
  Compiled           -> Automatic  (Batch-4: Automatic tries to Compile f for speed;
                                   typical 10-50x faster for large numerical problems)

Algorithms
  \"GradientDescent\"     steepest-descent with quadratic line search
  \"ConjugateGradientPR\" conjugate gradient, Polak-Ribiere formula
  \"ConjugateGradientFR\" conjugate gradient, Fletcher-Reeves formula
  \"Newton\"              Newton's method using numerical Hessian
  \"BFGS\"                Broyden-Fletcher-Goldfarb-Shanno quasi-Newton
  \"DFP\"                 Davidon-Fletcher-Powell quasi-Newton
  \"NelderMead\"          Nelder-Mead simplex (derivative-free)
  \"SimulatedAnnealing\"  simulated annealing (global)
  \"GeneticAlgorithm\"    real-coded genetic algorithm (global)"

(* ============================================================
   Equation Root Finding
   ============================================================ *)

FindEquationRoot::usage =
"FindEquationRoot[f, {a, b}] finds a root of the univariate pure function f.

  Bracket methods (Bisection, RegulaFalsi, Brent, Muller):
    [a, b] must bracket a root, i.e. f(a)*f(b) < 0.
  Open methods (Newton, Halley, FixedPoint, Steffensen, Aitken):
    a is the starting iterate; b is ignored.
  Secant method: a = x0, b = x1 (two distinct starting points).

Returns <|\"Root\"->xp, \"Residual\"->|f(xp)|, \"Iterations\"->k, \"Converged\"->True/False|>.

Options
  Method             -> \"Bisection\" (default) | \"RegulaFalsi\" | \"Brent\" |
                        \"Newton\" | \"Secant\" | \"Halley\" | \"Muller\" |
                        \"FixedPoint\" | \"Steffensen\" | \"Aitken\"
  Tolerance          -> 1*^-8
  MaxIterations      -> 200
  ConvergenceHistory -> False  (when True, result includes \"ResidualHistory\")
  WorkingPrecision   -> MachinePrecision  (Batch-4: e.g. 30 for 30-digit roots)

Algorithms
  \"Bisection\"   bisection method (linear convergence, guaranteed)
  \"RegulaFalsi\" false-position / regula falsi (linear)
  \"Brent\"       Brent's method (superlinear, guaranteed)
  \"Newton\"      Newton-Raphson (quadratic convergence)
  \"Secant\"      secant method (order ~1.618)
  \"Halley\"      Halley's method (cubic convergence)
  \"Muller\"      Muller's parabolic method (order ~1.839)
  \"FixedPoint\"  fixed-point iteration (f is the iteration function g(x))
  \"Steffensen\"  Steffensen acceleration (quadratic)
  \"Aitken\"      Aitken delta-squared acceleration (superlinear)"

(* ============================================================
   Nonlinear Equation Systems
   ============================================================ *)

SolveNonlinearEquationSystem::usage =
"SolveNonlinearEquationSystem[F, x0] solves the nonlinear system F(x) = 0
starting from initial vector x0.
F is a pure function called as F[{x1,...,xn}] returning a list of n residuals.
Returns <|\"Solution\"->xp, \"Residual\"->\[LeftDoubleBracketingBar]F(xp)\[RightDoubleBracketingBar], \"Iterations\"->k, \"Converged\"->True/False|>.

Options
  Method             -> \"Newton\" (default) | \"Broyden\" | \"FixedPoint\" | \"Continuation\"
  Tolerance          -> 1*^-8
  MaxIterations      -> 500
  Jacobian           -> Automatic | analytical-Jacobian-function
  ConvergenceHistory -> False  (when True, result includes \"ResidualHistory\")

Algorithms
  \"Newton\"       multivariate Newton's method (quadratic)
  \"Broyden\"      Broyden rank-1 quasi-Newton update
  \"FixedPoint\"   vector fixed-point iteration
  \"Continuation\" numerical continuation / homotopy method"

(* ============================================================
   Numerical Quadrature
   ============================================================ *)

NumericalQuadrature::usage =
"NumericalQuadrature[f, {x, a, b}] numerically evaluates Integral[f[x], {x, a, b}].
f is a pure function called as f[x] and must return a real number.
Returns a real number.

Options
  Method            -> \"GaussLegendre\" (default) | \"Trapezoidal\" | \"Simpson\" |
                       \"AdaptiveSimpson\" | \"Romberg\" | \"GaussChebyshev\" |
                       \"GaussHermite\" | \"GaussLaguerre\" | \"GaussLobatto\" | \"GaussRadau\"
  Points            -> 5    (nodes for Gauss-type rules)
  Intervals         -> 100  (subintervals for composite rules)
  Tolerance         -> 1*^-6
  MaxRecursionDepth -> 60   (maximum recursion depth for AdaptiveSimpson)
  Weighted          -> True (GaussChebyshev/Hermite/Laguerre: True = weighted integral,
                             False = standard integral of f[x])

Algorithms
  \"Trapezoidal\"    composite trapezoidal rule
  \"Simpson\"        composite Simpson's rule
  \"AdaptiveSimpson\" recursive adaptive Simpson (error-controlled)
  \"Romberg\"        Romberg's method via Richardson extrapolation
  \"GaussLegendre\"  Gauss-Legendre quadrature on [a, b]
  \"GaussChebyshev\" Gauss-Chebyshev (type I); integrates f(x)/Sqrt[1-x^2]
  \"GaussHermite\"   Gauss-Hermite; integrates f(x)*Exp[-x^2] on (-Inf,Inf)
  \"GaussLaguerre\"  Gauss-Laguerre; integrates f(x)*Exp[-x] on [0, Inf)
  \"GaussLobatto\"   Gauss-Lobatto (includes endpoints)
  \"GaussRadau\"     Gauss-Radau (includes left endpoint)"

(* ============================================================
   ODE Initial Value Problems
   ============================================================ *)

SolveInitialValueProblem::usage =
"SolveInitialValueProblem[f, {t, t0, tEnd}, y0] solves the first-order ODE
initial value problem y'(t) = f(t, y), y(t0) = y0.
For systems: f[t, {y1,...}] returns a list; y0 is a list.
Returns <|\"Grid\"->{t0,...,tn}, \"Solution\"->{y0,...,yn}|>.

Options
  Method   -> \"RungeKutta4\" (default) | \"Euler\" | \"ImplicitEuler\" |
              \"Heun\" | \"RungeKutta3\" | \"RungeKuttaFehlberg\" |
              \"AdamsBashforth\" | \"AdamsPC4\" | \"HammingPC\" |
              \"TrapezoidalRule\" | \"BDF2\" | \"BDF4\"
  StepSize -> Automatic (= (tEnd - t0)/100)
  Compiled -> Automatic  (Batch-4: Automatic tries to Compile f[t,y]; skip with False)
Note (Batch-4): if y0 or f[t0,y0] is complex-valued, the system is automatically
split into a 2n real system, solved, and reconstructed as complex output.

Algorithms
  \"Euler\"              forward Euler (order 1)
  \"ImplicitEuler\"      implicit / backward Euler (order 1)
  \"Heun\"               Heun's method (order 2)
  \"RungeKutta3\"        classical 3rd-order Kutta method
  \"RungeKutta4\"        classical 4th-order Runge-Kutta
  \"RungeKuttaFehlberg\" RKF45 adaptive step-size control
  \"AdamsBashforth\"     4-step Adams-Bashforth explicit (order 4)
  \"AdamsPC4\"           Adams 4th-order predictor-corrector
  \"HammingPC\"          Hamming predictor-corrector
  \"TrapezoidalRule\"    implicit trapezoidal / Crank-Nicolson for ODE (A-stable, order 2)
  \"BDF2\"               2nd-order backward differentiation formula (A-stable, order 2)
  \"BDF4\"               4th-order backward differentiation formula (order 4)"

(* ============================================================
   Boundary Value Problems
   ============================================================ *)

SolveBoundaryValueProblem::usage =
"SolveBoundaryValueProblem[f, {x, a, b}, {ya, yb}] solves the second-order
ODE boundary value problem y''(x) = f(x, y, y'),  y(a) = ya,  y(b) = yb.
SolveBoundaryValueProblem[f, {x,a,b}, {ya,yb}, n, opts] specifies n interior nodes.
Returns <|\"Grid\"->{x0,...,xn}, \"Solution\"->{y0,...,yn}|>.

Options
  Method -> \"Shooting\" (default) | \"FiniteDifference\"

Algorithms
  \"Shooting\"        linear shooting method (superposition of two IVPs)
  \"FiniteDifference\" central finite differences (tridiagonal system)"

(* ============================================================
   Constrained Optimization  (Batch-3)
   ============================================================ *)

FindConstrainedMinimum::usage =
"FindConstrainedMinimum[f, x0] finds a local minimum of f subject to
box and/or equality constraints specified via options.
f is called as f[{x1,...,xn}] and must return a real number.
Returns <|\"Point\"->xp, \"Value\"->fp, \"Iterations\"->k, \"Converged\"->True/False|>.

Options
  LowerBounds              -> Automatic | list  (Automatic = -Inf)
  UpperBounds              -> Automatic | list  (Automatic = +Inf)
  EqualityConstraints      -> None | {h1,h2,...}  (h_i[x] = 0)
  Method                   -> Automatic | \"ProjectedGradient\" | \"AugmentedLagrangian\"
  Tolerance                -> 1*^-8
  MaxIterations            -> 1000
  AugmentedLagrangianPenalty -> 10.

Auto-selection: box only -> ProjectedGradient; equality present -> AugmentedLagrangian."

FindConstrainedMinimum::badmethod =
  "Unknown Method \"`1`\". Falling back to \"ProjectedGradient\"."

(* ============================================================
   Partial Differential Equations
   ============================================================ *)

SolveParabolicPDE::usage =
"SolveParabolicPDE[alpha2, {x, xL, xR, nx}, {t, t0, tEnd, nt}, ic, {bc0, bcL}]
solves the parabolic PDE (heat equation)  u_t = alpha2 * u_xx
on [xL, xR] x [t0, tEnd] with nx+2 spatial nodes and nt+1 time levels.
  ic[x]   -- initial condition u(x, t0)
  bc0[t]  -- left  boundary u(xL, t)
  bcL[t]  -- right boundary u(xR, t)
Returns <|\"SpatialGrid\"->xs, \"TimeGrid\"->ts, \"Solution\"->matrix|>
where matrix[[j, i]] = u(x_i, t_j).

Options
  Method -> \"CrankNicolson\" (default) | \"ForwardDifference\" | \"BackwardDifference\"

Algorithms
  \"CrankNicolson\"      Crank-Nicolson (2nd-order in t and x; unconditionally stable)
  \"ForwardDifference\"  explicit forward difference (stable when r = alpha2*dt/dx^2 <= 0.5)
  \"BackwardDifference\" implicit backward difference (unconditionally stable)"

SolveHyperbolicPDE::usage =
"SolveHyperbolicPDE[c2, {x, xL, xR, nx}, {t, t0, tEnd, nt}, ic, vic, {bc0, bcL}]
solves the hyperbolic PDE (wave equation)  u_tt = c2 * u_xx.
  ic[x]   -- initial displacement u(x, t0)
  vic[x]  -- initial velocity u_t(x, t0)
  bc0[t]  -- left  boundary u(xL, t)
  bcL[t]  -- right boundary u(xR, t)
Returns <|\"SpatialGrid\"->xs, \"TimeGrid\"->ts, \"Solution\"->matrix|>.
CFL stability condition:  Sqrt[c2] * dt / dx <= 1.

Options
  Method -> \"ExplicitDifference\" (default)

Algorithms
  \"ExplicitDifference\" explicit second-order finite difference"

(* ============================================================
   2D PDE Solvers  (Batch-3)
   ============================================================ *)

SolveEllipticPDE::usage =
"SolveEllipticPDE[f, {x,xL,xR,nx}, {y,yL,yR,ny}, {bcLeft,bcRight,bcBottom,bcTop}]
solves the Poisson equation  nabla^2 u = f(x,y)  on a uniform rectangular grid.
Dirichlet BCs on all four edges.
Returns <|\"SpatialGridX\"->xs, \"SpatialGridY\"->ys, \"Solution\"->matrix|>
where matrix[[j,i]] = u(x_i, y_j).  Uses a sparse 5-point stencil."

SolveParabolicPDE2D::usage =
"SolveParabolicPDE2D[alpha2, {x,xL,xR,nx}, {y,yL,yR,ny}, {t,t0,tEnd,nt},
  ic, {bcLeft,bcRight,bcBottom,bcTop}]
solves the 2D heat equation  u_t = alpha2*(u_xx+u_yy)  using the
Peaceman-Rachford ADI method (alternating direction implicit; order 2 in t and space).
ic[x,y] -- initial condition.
Returns <|\"SpatialGridX\"->xs, \"SpatialGridY\"->ys, \"TimeGrid\"->ts, \"Solution\"->{u0,...}|>."

SolveHyperbolicPDE2D::usage =
"SolveHyperbolicPDE2D[c2, {x,xL,xR,nx}, {y,yL,yR,ny}, {t,t0,tEnd,nt},
  ic, vic, {bcLeft,bcRight,bcBottom,bcTop}]
solves the 2D wave equation  u_tt = c2*(u_xx+u_yy)  using an explicit
second-order finite difference scheme.
ic[x,y] -- initial displacement;  vic[x,y] -- initial velocity.
CFL stability: Sqrt[c2]*dt*Sqrt[1/hx^2+1/hy^2] <= 1.
Returns <|\"SpatialGridX\"->xs, \"SpatialGridY\"->ys, \"TimeGrid\"->ts, \"Solution\"->{u0,...}|>."

SolveHyperbolicPDE2D::cfl =
  "CFL number `1` > 1. ExplicitDifference may be unstable."

(* ============================================================
   Scattered Data Interpolation  (Batch-3)
   ============================================================ *)

ScatteredInterpolate::usage =
"ScatteredInterpolate[points, values, q] interpolates scattered 2D data at
query point q = {xq, yq}.  points is a list of {xi,yi} locations.
Returns the interpolated scalar value.

Options
  Method     -> \"ThinPlateSpline\" (default) | \"NearestNeighbor\" | \"InverseDistance\"
  IDWPower   -> 2   (inverse-distance exponent p)
  Regularise -> 0.  (Tikhonov regularisation for ThinPlateSpline)

Algorithms
  \"NearestNeighbor\"  O(n) per query -- return value at closest data point
  \"InverseDistance\"  O(n) per query -- IDW with weight 1/||x-xi||^p
  \"ThinPlateSpline\"  O(n^3) setup, O(n) per query -- global RBF phi(r) = r^2 log(r)"

ScatteredInterpolate::badmethod =
  "Unknown Method \"`1`\". Falling back to \"ThinPlateSpline\"."
ScatteredInterpolate::baddims =
  "points and values must have the same length. Got `1` and `2`."
ScatteredInterpolate::toofewpoints =
  "ThinPlateSpline requires at least 3 data points. Got `1`."

(* ============================================================
   Linear Equation Systems
   ============================================================ *)

SolveLinearEquationSystem::usage =
"SolveLinearEquationSystem[A, b] solves the linear system A.x = b.
Returns the solution vector x.

Options
  Method        -> \"GaussianEliminationPivot\" (default) | \"GaussianElimination\" |
                   \"GaussJordan\" | \"LU\" | \"Cholesky\" | \"LDLT\" |
                   \"Tridiagonal\" | \"Jacobi\" | \"GaussSeidel\" | \"SOR\"
  Omega         -> 1.5   (relaxation factor omega for SOR)
  Tolerance     -> 1*^-8
  MaxIterations -> 1000

Direct methods
  \"GaussianEliminationPivot\" Gaussian elimination with partial pivoting
  \"GaussianElimination\"      Gaussian elimination without pivoting
  \"GaussJordan\"              Gauss-Jordan full elimination
  \"LU\"                       LU factorisation with partial pivoting
  \"Cholesky\"                 Cholesky (LL^T) for symmetric positive-definite A
  \"LDLT\"                     LDL^T factorisation for symmetric A
  \"Tridiagonal\"              Thomas algorithm (tridiagonal A only)
Iterative methods
  \"Jacobi\"             Jacobi iteration
  \"GaussSeidel\"        Gauss-Seidel iteration
  \"SOR\"                successive over-relaxation with factor omega
  \"ConjugateGradient\"  conjugate gradient for symmetric positive-definite A;
                        accepts both dense Matrix and SparseArray (Batch-3)
Sparse direct
  \"SparseDirectLU\"     converts A to SparseArray then uses Mathematica's
                        sparse direct solver (SuperLU/UMFPACK); best for large
                        sparse systems such as PDE-discretised matrices (Batch-3)"

LUDecompose::usage =
"LUDecompose[A] computes the LU factorisation of A with partial pivoting.
Returns {L, U, P} satisfying P.A = L.U,
where L is unit lower-triangular, U is upper-triangular, P is a permutation matrix."

QRDecompose::usage =
"QRDecompose[A] computes the QR factorisation of A.
Returns {Q, R} satisfying A = Q.R, Q orthogonal, R upper-triangular.
Option Method -> \"Householder\" (default) | \"GramSchmidt\""

CholeskyDecompose::usage =
"CholeskyDecompose[A] computes the Cholesky factorisation of symmetric positive-definite A.
Returns lower-triangular L satisfying A = L.Transpose[L]."

LDLTDecompose::usage =
"LDLTDecompose[A] computes the LDL^T factorisation of symmetric A.
Returns {L, d} satisfying A = L.DiagonalMatrix[d].Transpose[L],
where L is unit lower-triangular and d is the diagonal vector."

(* ============================================================
   Matrix Eigenanalysis
   ============================================================ *)

FindDominantEigenvalue::usage =
"FindDominantEigenvalue[A] finds the eigenvalue of A with the largest absolute value
(the dominant eigenvalue) and its eigenvector.
Returns <|\"Eigenvalue\"->lambda, \"Eigenvector\"->v, \"Iterations\"->k, \"Converged\"->True/False|>.

Options
  Method        -> \"PowerMethod\" (default) | \"InversePowerMethod\"
  Shift         -> 0  (shift mu; InversePowerMethod finds eigenvalue nearest Shift)
  Tolerance     -> 1*^-8
  MaxIterations -> 500

Algorithms
  \"PowerMethod\"        power iteration (finds dominant eigenvalue)
  \"InversePowerMethod\" shifted inverse power iteration"

FindMatrixEigenvalues::usage =
"FindMatrixEigenvalues[A] finds all eigenvalues (and eigenvectors) of square matrix A.
Returns <|\"Eigenvalues\"->{lambda1,...}, \"Eigenvectors\"->{v1,...}, \"Iterations\"->k|>.

Options
  Method        -> \"QRIteration\" (default) | \"JacobiMethod\"
  Tolerance     -> 1*^-8
  MaxIterations -> 500

Algorithms
  \"QRIteration\"  QR iteration with origin shift (general real matrices)
  \"JacobiMethod\" Jacobi rotation method (symmetric matrices only)"

(* ============================================================
   Polynomial Interpolation
   ============================================================ *)

PolynomialInterpolate::usage =
"PolynomialInterpolate[xs, ys, xq] evaluates the interpolating polynomial through
data nodes {xs[[i]], ys[[i]]} at query point xq.
Returns the interpolated value.

Options
  Method -> \"Newton\" (default) | \"Lagrange\"

Algorithms
  \"Newton\"   Newton divided-difference form, Horner evaluation
  \"Lagrange\" Lagrange interpolation"

HermiteInterpolate::usage =
"HermiteInterpolate[xs, ys, dys, xq] evaluates the Hermite interpolating polynomial
through nodes xs with function values ys and derivative values dys at query point xq.
Returns the interpolated value."

SplineInterpolate::usage =
"SplineInterpolate[xs, ys, xq] evaluates a spline interpolant at query point xq.
Returns the interpolated value.

Options
  Degree            -> 3       (1 = linear, 2 = quadratic, 3 = cubic)
  BoundaryCondition -> \"Natural\" | \"Clamped\" | \"Periodic\"

Algorithms
  Degree->1  piecewise linear spline
  Degree->2  piecewise quadratic spline
  Degree->3  natural cubic spline (default)"

(* ============================================================
   Polynomial Approximation
   ============================================================ *)

ChebyshevApproximate::usage =
"ChebyshevApproximate[f, {x, a, b}, n] approximates f on [a, b] with a
degree-n Chebyshev polynomial expansion.
Returns a pure function representing the approximant."

LeastSquaresApproximate::usage =
"LeastSquaresApproximate[xs, ys, basis] fits discrete data {xs[[i]], ys[[i]]}
to a linear combination Sum[c[[k]] * basis[[k]][x]].
Returns the coefficient vector {c1, ..., cm}.
LeastSquaresApproximate[xs, ys, basis, \"Orthogonal\"] uses QR factorisation."

PadeApproximate::usage =
"PadeApproximate[f, x0, {p, q}] computes the [p/q] Pade rational approximant of f
at expansion point x0.
Returns a pure function R[x] = numerator(x) / denominator(x)."

RemezApproximate::usage =
"RemezApproximate[f, {a, b}, n] finds the best uniform (minimax) polynomial
approximation of degree n to f on [a, b] via the Remez exchange algorithm.
Returns <|\"Coefficients\"->{c0,...,cn}, \"Error\"->Emax, \"Nodes\"->{x0,...}|>."

(* ============================================================
   Error messages
   ============================================================ *)

(* ── Method routing errors ───────────────────────────────────────────────── *)
FindMinimum1D::badmethod              = "Unknown Method \"`1`\". Falling back to \"GoldenSection\"."
FindMinimumND::badmethod              = "Unknown Method \"`1`\". Falling back to \"BFGS\"."
FindEquationRoot::badmethod           = "Unknown Method \"`1`\". Falling back to \"Bisection\"."
FindEquationRoot::bracket             = "Bracket methods require f(a)*f(b) < 0. Check the input interval."
SolveNonlinearEquationSystem::badmethod = "Unknown Method \"`1`\". Falling back to \"Newton\"."
NumericalQuadrature::badmethod        = "Unknown Method \"`1`\". Falling back to \"GaussLegendre\"."
NumericalQuadrature::maxdepth         = "AdaptiveSimpson reached maximum recursion depth (`1`). Result may be inaccurate. Increase MaxRecursionDepth or switch to \"Romberg\"."
SolveInitialValueProblem::badmethod   = "Unknown Method \"`1`\". Falling back to \"RungeKutta4\"."
SolveInitialValueProblem::shorttspan  = "Time span too short for multi-step method \"`1`\" (requires >= 3 steps). Falling back to RungeKutta4."
SolveBoundaryValueProblem::badmethod  = "Unknown Method \"`1`\". Falling back to \"Shooting\"."
SolveBoundaryValueProblem::nonlinear  = "FiniteDifference detected possible nonlinearity in f. Results may be inaccurate. Consider Method -> \"Shooting\"."
SolveParabolicPDE::badmethod          = "Unknown Method \"`1`\". Falling back to \"CrankNicolson\"."
SolveParabolicPDE::stiff              = "Stability parameter r = `1` > 0.5. ForwardDifference may be unstable."
SolveHyperbolicPDE::cfl               = "CFL number `1` > 1. ExplicitDifference may be unstable."
SolveLinearEquationSystem::badmethod  = "Unknown Method \"`1`\". Falling back to \"GaussianEliminationPivot\"."
FindDominantEigenvalue::badmethod     = "Unknown Method \"`1`\". Falling back to \"PowerMethod\"."
FindMatrixEigenvalues::badmethod      = "Unknown Method \"`1`\". Falling back to \"QRIteration\"."
PolynomialInterpolate::badmethod      = "Unknown Method \"`1`\". Falling back to \"Newton\"."
CholeskyDecompose::notspdefinite      = "Matrix is not symmetric positive-definite."

(* ── Input validation errors (Batch 2) ──────────────────────────────────── *)
FindMinimum1D::badinterval     = "Require a < b. Got a = `1`, b = `2`."
FindMinimum1D::badfunc         = "f does not return a real number at the left endpoint. Check f."
FindMinimumND::badx0           = "x0 must be a numeric vector."
FindMinimumND::badfunc         = "f[x0] does not return a real number. Check that f accepts a list."
FindEquationRoot::badfunc      = "f does not evaluate to a number at the starting point."
SolveNonlinearEquationSystem::baddims  = "F[x0] has length `1` but x0 has length `2`. They must match."
SolveNonlinearEquationSystem::badfunc  = "F[x0] does not return a numeric list. Check F."
SolveNonlinearEquationSystem::singular = "Jacobian is singular at iteration `1`. Returning current best estimate."
NumericalQuadrature::equalendpoints   = "a == b: returning 0 (integral over an empty interval)."
NumericalQuadrature::weighted         = "GaussChebyshev/GaussHermite/GaussLaguerre integrate weighted forms (see usage). Set Weighted -> False to integrate f[x] directly."
SolveInitialValueProblem::badtspan    = "Require t0 < tEnd. Got t0 = `1`, tEnd = `2`."
SolveInitialValueProblem::badstep     = "StepSize must be positive. Got `1`."
SolveLinearEquationSystem::notsquare  = "Matrix A must be square. Got `1` x `2`."
SolveLinearEquationSystem::baddims    = "A has `1` rows but b has length `2`. They must match."
SolveLinearEquationSystem::singular   = "Matrix is singular or nearly singular. No solution returned."
CholeskyDecompose::nonsymmetric       = "Cholesky requires a symmetric matrix. Input is not symmetric."
LDLTDecompose::nonsymmetric           = "LDL^T requires a symmetric matrix. Input is not symmetric."
FindMatrixEigenvalues::nonsymmetric   = "JacobiMethod requires a symmetric matrix. Use QRIteration for general matrices."
PolynomialInterpolate::baddims        = "xs and ys must have the same length. Got `1` and `2`."
PolynomialInterpolate::dupnodes       = "xs contains duplicate values. Interpolation nodes must be distinct."
SplineInterpolate::toofewpoints       = "Need at least `1` data points for Degree `2` spline. Got `3`."
SolveBoundaryValueProblem::badinterval = "Require a < b. Got a = `1`, b = `2`."

Begin["`Private`"]

With[{base = DirectoryName[$InputFileName]},
  (* Core utilities -- must load first *)
  Get[FileNameJoin[{base, "Modules", "Internal.wl"}]];
  (* Unconstrained optimisation and root finding *)
  Get[FileNameJoin[{base, "Modules", "UnconstrainedMinimization.wl"}]];
  Get[FileNameJoin[{base, "Modules", "NumericalQuadrature.wl"}]];
  Get[FileNameJoin[{base, "Modules", "EquationRootFinding.wl"}]];
  Get[FileNameJoin[{base, "Modules", "NonlinearEquationSystems.wl"}]];
  (* Linear algebra (defines qrHouseholder, laThomasAlgorithm) *)
  Get[FileNameJoin[{base, "Modules", "LinearEquationSystems.wl"}]];
  Get[FileNameJoin[{base, "Modules", "MatrixEigenanalysis.wl"}]];
  (* Interpolation and approximation *)
  Get[FileNameJoin[{base, "Modules", "PolynomialInterpolation.wl"}]];
  Get[FileNameJoin[{base, "Modules", "PolynomialApproximation.wl"}]];
  (* ODE (defines ivpRungeKutta4, ivpRK4StartUp, ivpMultiStepCheck) *)
  Get[FileNameJoin[{base, "Modules", "OrdinaryDifferentialEquations.wl"}]];
  Get[FileNameJoin[{base, "Modules", "StiffODESolvers.wl"}]];         (* Batch-3 *)
  Get[FileNameJoin[{base, "Modules", "BoundaryValueProblems.wl"}]];
  (* PDE: 1D then 2D *)
  Get[FileNameJoin[{base, "Modules", "PartialDifferentialEquations.wl"}]];
  Get[FileNameJoin[{base, "Modules", "PartialDifferentialEquations2D.wl"}]]; (* Batch-3 *)
  (* Batch-3: constrained optimisation and scattered interpolation *)
  Get[FileNameJoin[{base, "Modules", "ConstrainedOptimization.wl"}]];
  Get[FileNameJoin[{base, "Modules", "ScatteredInterpolation.wl"}]]
]

End[]

EndPackage[]
