(* NonlinearEquationSystems.wl -- Nonlinear equation system solvers F(x)=0
   Methods: Newton, Broyden, FixedPointIter, Continuation

   Batch-2 improvements
     * Validates F[x0] returns a list of the same length as x0
     * safeLinearSolve replaces raw LinearSolve; singular Jacobian issues
       a domain message and returns the current best estimate rather than crashing
     * ConvergenceHistory option: when True, result includes "ResidualHistory"
       (list of ||F(x_k)|| per iteration)

   Private helper naming convention: camelCase prefixed by "nlsys"
   Note: "FixedPoint" is a protected Mathematica symbol; helper is nlsysFixedPointIter.
*)

Options[SolveNonlinearEquationSystem] = {
  Method             -> "Newton",
  Tolerance          -> 1*^-8,
  MaxIterations      -> 500,
  Jacobian           -> Automatic,
  ConvergenceHistory -> False
}

SolveNonlinearEquationSystem[F_, x0_?VectorQ, opts:OptionsPattern[]] :=
  Module[{method, tol, maxIter, jacFn, track, x0r, Fx0, raw, result},
    (* Validation *)
    x0r = Quiet[N@x0];
    If[!numericVectorQ[x0r],
      Message[SolveNonlinearEquationSystem::badx0]; Return[$Failed, Module]];
    Fx0 = Quiet[F[x0r]];
    If[!VectorQ[Fx0, NumericQ],
      Message[SolveNonlinearEquationSystem::badfunc]; Return[$Failed, Module]];
    If[Length[Fx0] =!= Length[x0r],
      Message[SolveNonlinearEquationSystem::baddims, Length[Fx0], Length[x0r]];
      Return[$Failed, Module]];
    method  = OptionValue[Method];
    tol     = OptionValue[Tolerance];
    maxIter = OptionValue[MaxIterations];
    track   = TrueQ[OptionValue[ConvergenceHistory]];
    jacFn   = With[{j = OptionValue[Jacobian]},
                If[j === Automatic, numJacobian[F, #]&, j]];
    If[track,
      raw = Reap[
        result = Switch[method,
          "Newton",       nlsysNewton[F, jacFn, x0r, tol, maxIter, True],
          "Broyden",      nlsysBroyden[F, x0r, tol, maxIter, True],
          "FixedPoint",   nlsysFixedPointIter[F, x0r, tol, maxIter, True],
          "Continuation", nlsysContinuation[F, jacFn, x0r, tol, maxIter],
          _,
            Message[SolveNonlinearEquationSystem::badmethod, method];
            nlsysNewton[F, jacFn, x0r, tol, maxIter, True]
        ]
      , "nkhist"];
      If[result === $Failed, Return[$Failed, Module]];
      Append[result,
        "ResidualHistory" -> If[Last[raw] === {}, {}, First[Last[raw]]]],
      Switch[method,
        "Newton",       nlsysNewton[F, jacFn, x0r, tol, maxIter],
        "Broyden",      nlsysBroyden[F, x0r, tol, maxIter],
        "FixedPoint",   nlsysFixedPointIter[F, x0r, tol, maxIter],
        "Continuation", nlsysContinuation[F, jacFn, x0r, tol, maxIter],
        _,
          Message[SolveNonlinearEquationSystem::badmethod, method];
          nlsysNewton[F, jacFn, x0r, tol, maxIter]
      ]
    ]
  ]

(* ── Multivariate Newton's Method ───────────────────────────────────────── *)
nlsysNewton[F_, jacFn_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = x0, J, Fx, dx, err, k = 0, converged = False},
    Fx = F[x];
    If[track, Sow[Norm[Fx], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        err = Norm[Fx];
        If[err < tol, converged = True; Throw[Null]];
        J  = jacFn[x];
        dx = safeLinearSolve[J, -Fx];
        If[dx === $Failed,
          Message[SolveNonlinearEquationSystem::singular, k];
          Throw[Null]
        ];
        x  += dx; Fx = F[x];
        If[track, Sow[Norm[Fx], "nkhist"]];
        If[Norm[dx]/(Norm[x] + $MachineEpsilon) < tol,
          converged = True; Throw[Null]]
      ]
    ];
    <|"Solution" -> x, "Residual" -> Norm[F[x]],
      "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Broyden Rank-1 Quasi-Newton ────────────────────────────────────────── *)
nlsysBroyden[F_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = x0, J = numJacobian[F, x0], Fx = F[x0],
          dx, xNew, FxNew, s, y, k = 0, converged = False},
    If[track, Sow[Norm[Fx], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        If[Norm[Fx] < tol, converged = True; Throw[Null]];
        dx    = safeLinearSolve[J, -Fx];
        If[dx === $Failed,
          Message[SolveNonlinearEquationSystem::singular, k]; Throw[Null]];
        xNew  = x + dx; FxNew = F[xNew];
        s = xNew - x; y = FxNew - Fx;
        J = J + Outer[Times, y - J . s, s] / (s . s + $MachineEpsilon);
        x = xNew; Fx = FxNew;
        If[track, Sow[Norm[Fx], "nkhist"]]
      ]
    ];
    <|"Solution" -> x, "Residual" -> Norm[F[x]],
      "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Vector Fixed-Point Iteration ──────────────────────────────────────── *)
nlsysFixedPointIter[G_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = x0, xp, k = 0, converged = False, diverged = False},
    If[track, Sow[Norm[G[x] - x], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        xp = G[x];
        If[track, Sow[Norm[xp - x], "nkhist"]];
        If[Norm[xp - x] < tol, converged = True; x = xp; Throw[Null]];
        x = xp;
        (* Bail out once the iterate runs away, before the next G[x]
           evaluation would overflow machine arithmetic and spam messages *)
        If[!(VectorQ[x, NumericQ] && Max[Abs[x]] < 1.*^150),
          diverged = True; Throw[Null]]
      ]
    ];
    <|"Solution" -> x,
      "Residual" -> If[diverged, Infinity, Norm[G[x] - x]],
      "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Numerical Continuation / Homotopy Method ──────────────────────────── *)
nlsysContinuation[F_, jacFn_, x0_, tol_, maxIter_, nSteps_:20] :=
  Module[{x = x0, lambda = 0., dlambda = 1./nSteps, n = Length[x0], G, jacG, result},
    While[lambda < 1. - $MachineEpsilon,
      lambda = Min[lambda + dlambda, 1.];
      G      = Function[xx, lambda F[xx] + (1 - lambda)(xx - x0)];
      (* Newton on the homotopy G needs the Jacobian of G itself,
         not of the original F: J_G = lambda J_F + (1-lambda) I *)
      jacG   = Function[xx, lambda jacFn[xx] + (1 - lambda) IdentityMatrix[n]];
      result = nlsysNewton[G, jacG, x, tol, maxIter];
      If[!result["Converged"],
        Return[<|"Solution" -> x, "Residual" -> Norm[F[x]],
                 "Iterations" -> result["Iterations"], "Converged" -> False|>, Module]];
      x = result["Solution"]
    ];
    nlsysNewton[F, jacFn, x, tol, maxIter]
  ]
