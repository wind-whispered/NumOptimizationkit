(* BoundaryValueProblems.wl -- ODE boundary value problem solvers
   SolveBoundaryValueProblem solves y''(x) = f(x, y, y'), y(a)=ya, y(b)=yb.
   Methods: Shooting, FiniteDifference

   Batch-1 improvement
     * bvpFiniteDifference now handles the GENERAL LINEAR BVP
       y''(x) = p(x)·y'(x) + q(x)·y(x) + r(x)
       by extracting p, q, r via central finite differences of f,
       then building the correct non-symmetric tridiagonal system.
       The previous implementation incorrectly assumed p=q=0.
     * Nonlinearity detection: warns when f appears nonlinear in y or y'
       and recommends the Shooting method instead.

   Private helper naming convention: camelCase prefixed by "bvp"
   Depends on: ivpRungeKutta4 (OrdinaryDifferentialEquations.wl)
               laThomasAlgorithm (LinearEquationSystems.wl)
*)

Options[SolveBoundaryValueProblem] = {
  Method -> "Shooting"
}

SolveBoundaryValueProblem::nonlinear =
  "FiniteDifference detected that f may be nonlinear in y or y' \
(second-order finite-difference check failed at some grid points). \
Results may be inaccurate for nonlinear BVPs. Consider Method -> \"Shooting\"."

SolveBoundaryValueProblem[f_, {x_, a_?NumericQ, b_?NumericQ},
    {ya_?NumericQ, yb_?NumericQ}, n_Integer:100, opts:OptionsPattern[]] :=
  Module[{method},
    If[N@a >= N@b,
      Message[SolveBoundaryValueProblem::badinterval, a, b];
      Return[$Failed, Module]];
    method = OptionValue[Method];
    Switch[method,
      "Shooting",         bvpShooting[f, N@a, N@b, N@ya, N@yb, n],
      "FiniteDifference", bvpFiniteDifference[f, N@a, N@b, N@ya, N@yb, n],
      _,
        Message[SolveBoundaryValueProblem::badmethod, method];
        bvpShooting[f, N@a, N@b, N@ya, N@yb, n]
    ]
  ]

(* ── Shooting Method (linear superposition of two IVPs) ─────────────────── *)
(* Converts y'' = f(x, y, y') to {y1'=y2, y2'=f}; solves two IVPs and
   combines via superposition so that y(b) = yb. Works for linear and
   nonlinear BVPs (nonlinear: superposition is approximate; use Newton
   shooting for full nonlinear accuracy). *)
bvpShooting[f_, a_, b_, ya_, yb_, n_] :=
  Module[{h = (b-a)/n, F, sol1, sol2, u1, u2, c},
    F = Function[{t, vec}, {vec[[2]], f[t, vec[[1]], vec[[2]]]}];
    sol1 = ivpRungeKutta4[F, a, b, {ya, 0.}, h];
    u1   = sol1["Solution"];
    sol2 = ivpRungeKutta4[F, a, b, {0., 1.}, h];
    u2   = sol2["Solution"];
    c = (yb - u1[[-1, 1]]) / (u2[[-1, 1]] + $MachineEpsilon);
    <|"Grid"     -> sol1["Grid"],
      "Solution" -> (u1[[All, 1]] + c u2[[All, 1]])|>
  ]

(* ── Finite Difference Method ───────────────────────────────────────────── *)
(*
   For the GENERAL LINEAR BVP  y'' = p(x)·y' + q(x)·y + r(x)
   the central-difference discretisation (uniform grid, step h) is:

     (y_{i+1} - 2y_i + y_{i-1})/h^2 = p_i·(y_{i+1}-y_{i-1})/(2h) + q_i·y_i + r_i

   Rearranging into tridiagonal form:
     lower_i = 1/h^2 - p_i/(2h)      (coefficient of y_{i-1}, i=2..n)
     diag_i  = -2/h^2 + q_i          (coefficient of y_i,     i=1..n)
     upper_i = 1/h^2 + p_i/(2h)      (coefficient of y_{i+1}, i=1..n-1)
     rhs_i   = r_i

   p_i, q_i, r_i are extracted via numerical differentiation of f:
     r_i =  f(x_i, 0, 0)
     q_i = (f(x_i, eps, 0) - f(x_i, -eps, 0)) / (2*eps)   [df/dy]
     p_i = (f(x_i, 0, eps) - f(x_i, 0, -eps)) / (2*eps)   [df/dy']

   Nonlinearity check: if the second mixed difference of f w.r.t. y or y'
   is large relative to the first, the problem is nonlinear and a warning
   is issued recommending the Shooting method.
*)
bvpFiniteDifference[f_, a_, b_, ya_, yb_, n_] :=
  Module[{h = (b-a)/(n+1), xs, eps = 1*^-6, p, q, r,
          lower, diag, upper, rhs, nonlinearFlag = False, yi},

    xs = a + Range[1, n] h;

    (* Extract linear coefficients at each interior node *)
    r = Table[f[xs[[i]], 0., 0.], {i, n}] // N;
    q = Table[(f[xs[[i]], eps, 0.] - f[xs[[i]], -eps, 0.])/(2 eps), {i, n}] // N;
    p = Table[(f[xs[[i]], 0., eps] - f[xs[[i]], 0., -eps])/(2 eps), {i, n}] // N;

    (* Nonlinearity check: second-order differences should be ~0 for linear f *)
    Do[
      With[{d2y  = (f[xs[[i]], eps, 0.] - 2 f[xs[[i]], 0., 0.] + f[xs[[i]], -eps, 0.]) / eps^2,
            d2dy = (f[xs[[i]], 0., eps] - 2 f[xs[[i]], 0., 0.] + f[xs[[i]], 0., -eps]) / eps^2},
        If[Abs[d2y] > 1*^-4 || Abs[d2dy] > 1*^-4,
          nonlinearFlag = True; Break[]
        ]
      ],
      {i, n}
    ];
    If[nonlinearFlag, Message[SolveBoundaryValueProblem::nonlinear]];

    (* Build tridiagonal system *)
    lower = Table[1./h^2 - p[[i+1]]/(2h), {i, n-1}];   (* rows 2..n *)
    diag  = Table[-2./h^2 + q[[i]],       {i, n}];
    upper = Table[1./h^2 + p[[i]]/(2h),   {i, n-1}];   (* rows 1..n-1 *)
    rhs   = N@r;

    (* Apply boundary conditions *)
    rhs[[1]]  -= (1./h^2 - p[[1]]/(2h))  * ya;
    rhs[[-1]] -= (1./h^2 + p[[-1]]/(2h)) * yb;

    yi = laThomasAlgorithm[lower, diag, upper, rhs];
    <|"Grid"     -> Join[{a}, xs, {b}],
      "Solution" -> Join[{ya}, yi, {yb}]|>
  ]
