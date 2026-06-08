(* StiffODESolvers.wl -- Implicit multi-step ODE solvers for stiff problems
   New methods for SolveInitialValueProblem:
     "TrapezoidalRule" -- implicit trapezoidal (Crank-Nicolson for ODE), order 2
     "BDF2"            -- 2nd-order backward differentiation formula
     "BDF4"            -- 4th-order backward differentiation formula

   All three reduce to the same implicit residual form
     G(y) = y - rhs - coeff * f(t_{n+1}, y) = 0
   solved per step by Newton iteration with numerical Jacobian.

   Depends on: ivpRK4StartUp, ivpRungeKutta4 (OrdinaryDifferentialEquations.wl)
               safeLinearSolve, numJacobian, numDeriv (Internal.wl)

   Private helper naming convention: camelCase prefixed by "ivp"
*)

(* ── Shared implicit-step Newton solver ─────────────────────────────────── *)
(* Solves  G(y) = y - rhs - coeff*f(tNew, y) = 0  by Newton iteration.
   Works for both scalar and vector y. *)
ivpImplicitSolve[f_, tNew_, rhs_, coeff_, yPred_, maxIter_:15, tol_:1*^-10] :=
  Module[{y = yPred, Gval, JG, dY, isVec = VectorQ[yPred]},
    Catch[
      Do[
        Gval = y - rhs - coeff f[tNew, y];
        JG   = If[isVec,
                 IdentityMatrix[Length@y]*1. - coeff numJacobian[f[tNew, #]&, y],
                 1. - coeff numDeriv[f[tNew, #]&, y, 1]
               ];
        dY   = If[isVec,
                 safeLinearSolve[JG, -Gval],
                 If[Abs[JG] < $MachineEpsilon, Throw[$Failed], -Gval/JG]
               ];
        If[dY === $Failed, Throw[$Failed]];
        y += dY;
        If[Norm[Flatten[{dY}]] < tol, Throw[y]],
        {maxIter}
      ];
      y  (* return best estimate even if not fully converged *)
    ]
  ]

(* ── Implicit Trapezoidal Rule (order 2, A-stable) ──────────────────────── *)
(*   y_{n+1} = y_n + h/2 (f(t_n,y_n) + f(t_{n+1},y_{n+1}))
     G(y)    = y - [y_n + h/2 f(t_n,y_n)] - h/2 f(t_{n+1},y) = 0 *)
ivpTrapezoidalRule[f_, t0_, tEnd_, y0_, h_] :=
  Module[{nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t = t0, y = y0, rhs, yNew, i = 1},
    ts = ConstantArray[0.,  nMax];
    ys = ConstantArray[y0, nMax];
    ts[[1]] = t0; ys[[1]] = y0;
    While[t + h/2 <= tEnd,
      rhs  = y + h/2 f[t, y];   (* explicit half-step *)
      yNew = ivpImplicitSolve[f, t+h, rhs, h/2, rhs];
      If[yNew === $Failed, yNew = rhs];  (* fallback on failure *)
      t += h; i++;
      y = yNew; ts[[i]] = t; ys[[i]] = y
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── BDF-2 (order 2, A-stable) ──────────────────────────────────────────── *)
(*   (3y_{n+1} - 4y_n + y_{n-1}) / (2h) = f(t_{n+1}, y_{n+1})
     y_{n+1} = (4/3)y_n - (1/3)y_{n-1} + (2h/3) f(t_{n+1}, y_{n+1}) *)
ivpBDF2[f_, t0_, tEnd_, y0_, h_] :=
  Module[{check, startup, nStart, nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t, i, rhs, yNew},
    check = ivpMultiStepCheck["BDF2", t0, tEnd, y0, h,
              ivpRungeKutta4[f, t0, tEnd, y0, h]];
    If[check =!= Null, Return[check, Module]];
    startup = ivpRK4StartUp[f, t0, 2, y0, h];
    nStart  = Length[startup["Grid"]];
    ts = ConstantArray[0.,  nMax]; ys = ConstantArray[y0, nMax];
    ts[[1;;nStart]] = startup["Grid"];
    ys[[1;;nStart]] = startup["Solution"];
    i = nStart; t = ts[[i]];
    While[t + h/2 <= tEnd,
      rhs  = (4/3) ys[[i]] - (1/3) ys[[i-1]];   (* extrapolated predictor *)
      yNew = ivpImplicitSolve[f, t+h, rhs, 2h/3, rhs];
      If[yNew === $Failed, yNew = rhs];
      t += h; i++;
      ts[[i]] = t; ys[[i]] = yNew
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── BDF-4 (order 4, A-stable up to a sector) ───────────────────────────── *)
(*   (25y_{n+1} - 48y_n + 36y_{n-1} - 16y_{n-2} + 3y_{n-3}) / (12h) = f(t_{n+1}, y_{n+1})
     y_{n+1} = (48y_n - 36y_{n-1} + 16y_{n-2} - 3y_{n-3})/25 + (12h/25) f(t_{n+1}, y_{n+1}) *)
ivpBDF4[f_, t0_, tEnd_, y0_, h_] :=
  Module[{check, startup, nStart, nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t, i, rhs, yNew},
    check = ivpMultiStepCheck["BDF4", t0, tEnd, y0, h,
              ivpRungeKutta4[f, t0, tEnd, y0, h]];
    If[check =!= Null, Return[check, Module]];
    startup = ivpRK4StartUp[f, t0, 4, y0, h];
    nStart  = Length[startup["Grid"]];
    ts = ConstantArray[0.,  nMax]; ys = ConstantArray[y0, nMax];
    ts[[1;;nStart]] = startup["Grid"];
    ys[[1;;nStart]] = startup["Solution"];
    i = nStart; t = ts[[i]];
    While[t + h/2 <= tEnd,
      rhs  = (48 ys[[i]] - 36 ys[[i-1]] + 16 ys[[i-2]] - 3 ys[[i-3]]) / 25;
      yNew = ivpImplicitSolve[f, t+h, rhs, 12h/25, rhs];
      If[yNew === $Failed, yNew = rhs];
      t += h; i++;
      ts[[i]] = t; ys[[i]] = yNew
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]
