(* Internal.wl -- Private utilities shared across all modules
   All symbols defined here live in NumOptimizationkit`Private`

   Batch-4 additions
   4.1  WorkingPrecision support
        nkWP           -- session-level precision setting (Block-scoped by dispatchers)
        nkN[x]         -- precision-aware N[x]
        nkEps[]        -- precision-aware machine epsilon
        nkZero[]       -- precision-aware zero
        Numerical differentiation helpers updated to use nkEps[]
   4.2  Compile acceleration helpers
        tryCompileIVP  -- attempt to Compile an ODE RHS; falls back silently
        tryCompileND   -- attempt to Compile a multivariate objective
   4.4  Complex ODE detection
        ivpIsComplex   -- True when y0 or f[t,y] is complex-valued
        ivpComplexSplit -- real/imaginary split wrapper for complex IVPs
*)

(* ══════════════════════════════════════════════════════════════════════════
   4.1  Working-precision infrastructure
   ══════════════════════════════════════════════════════════════════════════ *)

(* Session-level precision; overridden per-call via Block[{nkWP = wp}, ...]. *)
nkWP = MachinePrecision

(* Precision-aware numeric conversion *)
nkN[x_]      := If[nkWP === MachinePrecision, N[x], N[x, nkWP]]
nkN[x_, p_]  := N[x, p]

(* Precision-aware convergence epsilon.
   For machine precision: $MachineEpsilon (~2e-16).
   For exact precision p: 10^(2 - p)  (two guard digits).  *)
nkEps[] := If[nkWP === MachinePrecision, $MachineEpsilon, 10^(2 - nkWP)]

(* Precision-aware zero scalar *)
nkZero[] := If[nkWP === MachinePrecision, 0., N[0, nkWP]]

(* ── Numerical differentiation  (precision-aware step sizes) ─────────────── *)

(* Optimal central-difference step for gradient/Jacobian: eps^(1/3) *)
nkGradStep[]  := nkEps[]^(1/3)
(* Optimal step for Hessian: eps^(1/4) *)
nkHessStep[]  := nkEps[]^(1/4)

numGrad[f_, x_?VectorQ] :=
  Module[{h = nkGradStep[], n = Length[x]},
    Table[
      (f[ReplacePart[x, i -> x[[i]] + h]] - f[ReplacePart[x, i -> x[[i]] - h]]) / (2 h),
      {i, n}
    ]
  ]

numDeriv[f_, x_?NumericQ, order_:1] :=
  Module[{h = nkGradStep[]},
    If[order == 1,
      (f[x + h] - f[x - h]) / (2 h),
      (f[x + h] - 2 f[x] + f[x - h]) / h^2
    ]
  ]

numJacobian[F_, x_?VectorQ] :=
  Module[{h = nkGradStep[], n = Length[x]},
    Transpose@Table[
      (F[ReplacePart[x, i -> x[[i]] + h]] - F[ReplacePart[x, i -> x[[i]] - h]]) / (2 h),
      {i, n}
    ]
  ]

numHessian[f_, x_?VectorQ] :=
  Module[{h = nkHessStep[], n = Length[x], fx = f[x], ei, ej},
    Table[
      ei = UnitVector[n, i]; ej = UnitVector[n, j];
      If[i == j,
        (f[x + h ei] - 2 fx + f[x - h ei]) / h^2,
        (f[x + h (ei + ej)] - f[x + h (ei - ej)] -
         f[x + h (ej - ei)] + f[x - h (ei + ej)]) / (4 h^2)
      ],
      {i, n}, {j, n}
    ]
  ]

(* ── Line search ─────────────────────────────────────────────────────────── *)

armijoLineSearch[f_, x0_, grad0_, dir_, alpha0_:1., rho_:0.01, beta_:0.5] :=
  Module[{alpha = nkN[alpha0], fx0 = f[x0], slope = grad0 . dir, k = 0},
    While[f[x0 + alpha dir] > fx0 + rho alpha slope && k++ < 50,
      alpha *= beta
    ];
    alpha
  ]

quadLineSearch[f_, x0_, grad0_, dir_, alpha0_:1.] :=
  Module[{l = 0.15, u = 0.85, rho = 0.01,
          alpha = nkN[alpha0], fx0 = f[x0], gd = grad0 . dir,
          x1 = x0, alpha1, fx, k = 0},
    Catch[
      While[k++ <= 50,
        fx = f[x1 + alpha dir];
        If[fx < fx0 + alpha rho gd, Throw[{x1 + alpha dir, fx, alpha, True}]];
        alpha1 = -gd alpha^2 0.5 / (fx - fx0 - alpha gd);
        alpha1 = Max[alpha1, l alpha];
        alpha  = Min[alpha1, u alpha]
      ];
      {x0, fx0, 0., False}
    ]
  ]

(* ── Safe LinearSolve wrapper ────────────────────────────────────────────── *)

safeLinearSolve[A_, b_] :=
  Quiet[
    Check[LinearSolve[A, b], $Failed,
      {LinearSolve::sing, LinearSolve::luc, LinearSolve::nosol}],
    {LinearSolve::sing, LinearSolve::luc, LinearSolve::nosol}
  ]

(* ── Input validation helpers ─────────────────────────────────────────────── *)

numericVectorQ[x_] := VectorQ[x, NumericQ]

realValuedAtQ[f_, x_] :=
  Quiet[With[{v = f[x]}, NumericQ[v] && Im[v] == 0], $MessageGroups]

toRealVec[x_List]     := N[x]
toRealVec[x_?NumericQ] := {N[x]}

(* ══════════════════════════════════════════════════════════════════════════
   4.2  Compile acceleration helpers
   ══════════════════════════════════════════════════════════════════════════ *)

(* Try to Compile an ODE RHS f[t, y] or f[t, {y1,...}].
   Returns CompiledFunction on success, original f on failure (silent fallback).
   isScalar: True when y0 is a scalar. *)
tryCompileIVP[f_, t0_, y0_] :=
  Module[{isScalar = !VectorQ[y0], cf},
    cf = Quiet@Check[
      If[isScalar,
        Compile[{{t, _Real}, {y, _Real}}, f[t, y],
          RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}],
        Compile[{{t, _Real}, {y, _Real, 1}}, f[t, y],
          RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}]
      ],
      $Failed
    ];
    If[Head[cf] === CompiledFunction, cf, f]
  ]

(* Try to Compile a multivariate objective f[{x1,...}] -> Real.
   Used by FindMinimumND when Compiled -> True. *)
tryCompileND[f_, x0_] :=
  Module[{n = Length[x0], cf},
    cf = Quiet@Check[
      Compile[{{x, _Real, 1}}, f[x],
        RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}],
      $Failed
    ];
    If[Head[cf] === CompiledFunction, cf, f]
  ]

(* ══════════════════════════════════════════════════════════════════════════
   4.4  Complex ODE helpers
   ══════════════════════════════════════════════════════════════════════════ *)

(* Test whether an IVP involves complex arithmetic. *)
ivpIsComplex[f_, t0_, y0_] :=
  Module[{fval},
    If[!FreeQ[y0, Complex], Return[True]];
    fval = Quiet[f[N@t0, N@y0]];
    !FreeQ[fval, Complex]
  ]

(* Split a complex IVP into a real 2n system, solve, then reconstruct.
   For scalar y0 = x+Iy, splits to {Re, Im}.
   For vector y0 = {z1,...,zn}, splits to {Re[z1],...,Re[zn], Im[z1],...,Im[zn]}. *)
ivpComplexSplit[f_, t0_, tEnd_, y0_, h_, method_] :=
  Module[{isVec = VectorQ[y0], n, y0r, fReal, sol, tsOut, zsOut},
    n = If[isVec, Length[y0], 1];
    (* Real initial vector: [Re part ... | Im part ...] *)
    y0r = If[isVec,
      Join[Re[y0], Im[y0]],
      {Re[y0], Im[y0]}
    ] // N;
    (* Real-valued surrogate RHS *)
    fReal = Function[{tt, yy},
      Module[{zc, res},
        zc = If[isVec,
          yy[[1;;n]] + I yy[[n+1;;2n]],
          yy[[1]] + I yy[[2]]
        ];
        res = f[tt, zc];
        If[isVec,
          Join[Re[res], Im[res]],
          {Re[res], Im[res]}
        ] // N
      ]
    ];
    (* Solve with the existing real dispatcher (avoid infinite recursion) *)
    sol = SolveInitialValueProblem[fReal, {t, t0, tEnd}, y0r,
      Method -> method, StepSize -> h];
    If[sol === $Failed, Return[$Failed, Module]];
    tsOut = sol["Grid"];
    zsOut = If[isVec,
      sol["Solution"][[All, 1;;n]] + I sol["Solution"][[All, n+1;;2n]],
      sol["Solution"][[All, 1]]    + I sol["Solution"][[All, 2]]
    ];
    <|"Grid" -> tsOut, "Solution" -> zsOut|>
  ]

(* ══════════════════════════════════════════════════════════════════════════
   4.3  Namespace documentation & protection
   Cross-module symbol dependency map (for reference and protection):
     LinearEquationSystems -> qrHouseholder, laThomasAlgorithm
     OrdinaryDifferentialEquations -> ivpRungeKutta4, ivpRK4StartUp,
                                      ivpMultiStepCheck
   These are protected after definition to prevent accidental redefinition.
   ══════════════════════════════════════════════════════════════════════════ *)
