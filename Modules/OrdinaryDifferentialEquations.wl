(* OrdinaryDifferentialEquations.wl -- ODE initial value problem solvers
   SolveInitialValueProblem solves y'(t) = f(t, y), y(t0) = y0.
   Methods: Euler, ImplicitEuler, Heun, RungeKutta3, RungeKutta4,
            RungeKuttaFehlberg, AdamsBashforth, AdamsPC4, HammingPC
            TrapezoidalRule, BDF2, BDF4  (Batch-3, in StiffODESolvers.wl)

   Batch-1: pre-allocated arrays, Reap/Sow for RKF45, safe multi-step indexing
   Batch-4 additions
     4.2 Compiled -> Automatic option: when set, attempts to Compile f[t,y]
         before the integration loop; typically 10-50x faster for numeric f.
     4.4 Complex ODE support: if y0 or f[t0,y0] involves complex numbers,
         the system is transparently split into a 2n real system and solved
         with the requested method; the complex solution is reconstructed.

   Private helper naming convention: camelCase prefixed by "ivp"
   Cross-module exports: ivpRungeKutta4, ivpRK4StartUp, ivpMultiStepCheck
*)

Options[SolveInitialValueProblem] = {
  Method   -> "RungeKutta4",
  StepSize -> Automatic,
  Compiled -> Automatic   (* Batch-4: Automatic tries to Compile f; False skips *)
}

SolveInitialValueProblem::shorttspan =
  "Time span too short for multi-step method \"`1`\" (requires >= 3 steps). \
Falling back to RungeKutta4."

SolveInitialValueProblem[f_, {t_, t0_?NumericQ, tEnd_?NumericQ}, y0_,
    opts:OptionsPattern[]] :=
  Module[{method, h, compiled, fUse},
    (* Validation *)
    If[N@t0 >= N@tEnd,
      Message[SolveInitialValueProblem::badtspan, t0, tEnd];
      Return[$Failed, Module]];
    method   = OptionValue[Method];
    compiled = OptionValue[Compiled];
    h        = OptionValue[StepSize];
    If[h === Automatic, h = (tEnd - t0)/100];
    h = N@h;
    If[h <= 0,
      Message[SolveInitialValueProblem::badstep, h];
      Return[$Failed, Module]];
    (* 4.4  Complex ODE: transparently split into real/imag system *)
    If[ivpIsComplex[f, t0, y0],
      Return[ivpComplexSplit[f, t0, tEnd, y0, h, method], Module]];
    (* 4.2  Compile: attempt to compile f for speed *)
    fUse = If[TrueQ[compiled] || compiled === Automatic,
      tryCompileIVP[f, N@t0, N@y0],
      f
    ];
    Switch[method,
      "Euler",              ivpEuler[fUse, N@t0, tEnd, N@y0, h],
      "ImplicitEuler",      ivpImplicitEuler[fUse, N@t0, tEnd, N@y0, h],
      "Heun",               ivpHeun[fUse, N@t0, tEnd, N@y0, h],
      "RungeKutta3",        ivpRungeKutta3[fUse, N@t0, tEnd, N@y0, h],
      "RungeKutta4",        ivpRungeKutta4[fUse, N@t0, tEnd, N@y0, h],
      "RungeKuttaFehlberg", ivpRungeKuttaFehlberg[fUse, N@t0, tEnd, N@y0, h],
      "AdamsBashforth",     ivpAdamsBashforth[fUse, N@t0, tEnd, N@y0, h],
      "AdamsPC4",           ivpAdamsPC4[fUse, N@t0, tEnd, N@y0, h],
      "HammingPC",          ivpHammingPC[fUse, N@t0, tEnd, N@y0, h],
      (* Batch-3: stiff solvers use original f (compilation less useful for Newton) *)
      "TrapezoidalRule",    ivpTrapezoidalRule[f, N@t0, tEnd, N@y0, h],
      "BDF2",               ivpBDF2[f, N@t0, tEnd, N@y0, h],
      "BDF4",               ivpBDF4[f, N@t0, tEnd, N@y0, h],
      _,
        Message[SolveInitialValueProblem::badmethod, method];
        ivpRungeKutta4[fUse, N@t0, tEnd, N@y0, h]
    ]
  ]

ivpResult[ts_, ys_] := <|"Grid" -> ts, "Solution" -> ys|>

(* ══════════════════════════════════════════════════════════════════════════
   Fixed-step methods — pre-allocated array pattern
   nMax = Ceiling[(tEnd-t0)/h] + 2 gives a safe upper bound.
   ConstantArray[y0, nMax] works for both scalar and vector y0:
     scalar y0  -> flat list  {y0, y0, ...}
     vector y0  -> 2-D list   {{y0...}, {y0...}, ...}
   Both support direct indexed assignment ys[[i]] = newY.
   ══════════════════════════════════════════════════════════════════════════ *)

(* ── Forward Euler (order 1) ────────────────────────────────────────────── *)
ivpEuler[f_, t0_, tEnd_, y0_, h_] :=
  Module[{nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t = t0, y = y0, i = 1},
    ts = ConstantArray[0.,  nMax];
    ys = ConstantArray[y0, nMax];
    ts[[1]] = t0; ys[[1]] = y0;
    While[t + h/2 <= tEnd,
      y = y + h f[t, y]; t += h; i++;
      ts[[i]] = t; ys[[i]] = y
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── Implicit Euler (backward Euler, fixed-point corrector) ─────────────── *)
ivpImplicitEuler[f_, t0_, tEnd_, y0_, h_] :=
  Module[{nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t = t0, y = y0, yNew, i = 1},
    ts = ConstantArray[0.,  nMax];
    ys = ConstantArray[y0, nMax];
    ts[[1]] = t0; ys[[1]] = y0;
    While[t + h/2 <= tEnd,
      yNew = y + h f[t, y];
      Do[yNew = y + h f[t + h, yNew], {20}];
      t += h; y = yNew; i++;
      ts[[i]] = t; ys[[i]] = y
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── Heun's Method (explicit trapezoidal, order 2) ──────────────────────── *)
ivpHeun[f_, t0_, tEnd_, y0_, h_] :=
  Module[{nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t = t0, y = y0, k1, k2, i = 1},
    ts = ConstantArray[0.,  nMax];
    ys = ConstantArray[y0, nMax];
    ts[[1]] = t0; ys[[1]] = y0;
    While[t + h/2 <= tEnd,
      k1 = f[t, y]; k2 = f[t+h, y+h k1];
      y = y + h/2 (k1+k2); t += h; i++;
      ts[[i]] = t; ys[[i]] = y
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── 3rd-Order Kutta Method ─────────────────────────────────────────────── *)
ivpRungeKutta3[f_, t0_, tEnd_, y0_, h_] :=
  Module[{nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t = t0, y = y0, k1, k2, k3, i = 1},
    ts = ConstantArray[0.,  nMax];
    ys = ConstantArray[y0, nMax];
    ts[[1]] = t0; ys[[1]] = y0;
    While[t + h/2 <= tEnd,
      k1 = f[t, y];
      k2 = f[t + h/2, y + h/2 k1];
      k3 = f[t + h,   y - h k1 + 2h k2];
      y = y + h/6 (k1 + 4k2 + k3); t += h; i++;
      ts[[i]] = t; ys[[i]] = y
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* Unprotect in case this module is being reloaded (cf. Protect[] at end of file) *)
Unprotect[ivpRungeKutta4, ivpRK4StartUp, ivpMultiStepCheck];

(* ── Classical 4th-Order Runge-Kutta (also used by BoundaryValueProblems) ── *)
ivpRungeKutta4[f_, t0_, tEnd_, y0_, h_] :=
  Module[{nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, t = t0, y = y0, k1, k2, k3, k4, i = 1},
    ts = ConstantArray[0.,  nMax];
    ys = ConstantArray[y0, nMax];
    ts[[1]] = t0; ys[[1]] = y0;
    While[t + h/2 <= tEnd,
      k1 = f[t,       y];
      k2 = f[t + h/2, y + h/2 k1];
      k3 = f[t + h/2, y + h/2 k2];
      k4 = f[t + h,   y + h k3];
      y = y + h/6 (k1 + 2k2 + 2k3 + k4); t += h; i++;
      ts[[i]] = t; ys[[i]] = y
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── RKF45 Adaptive Step-Size Control ──────────────────────────────────── *)
(* Step count is not known in advance, so Reap/Sow is used instead of
   pre-allocation. Reap is significantly faster than repeated AppendTo. *)
ivpRungeKuttaFehlberg[f_, t0_, tEnd_, y0_, h0_, tol_:1*^-6] :=
  Module[{t = t0, y = y0, h = h0,
          k1, k2, k3, k4, k5, k6, y4, y5, err, hNew, raw},
    raw = Last@Reap[
      Sow[{t, y}];
      While[t < tEnd,
        h = Min[h, tEnd - t];
        k1 = h f[t,            y];
        k2 = h f[t + h/4,      y + k1/4];
        k3 = h f[t + 3h/8,     y + 3k1/32  + 9k2/32];
        k4 = h f[t + 12h/13,   y + 1932k1/2197 - 7200k2/2197 + 7296k3/2197];
        k5 = h f[t + h,        y + 439k1/216 - 8k2 + 3680k3/513 - 845k4/4104];
        k6 = h f[t + h/2,      y - 8k1/27 + 2k2 - 3544k3/2565 + 1859k4/4104 - 11k5/40];
        y4 = y + 25k1/216 + 1408k3/2565 + 2197k4/4104 - k5/5;
        y5 = y + 16k1/135 + 6656k3/12825 + 28561k4/56430 - 9k5/50 + 2k6/55;
        err = Norm[Flatten[{y5 - y4}]];
        If[err < tol || h < 1*^-12,
          t += h; y = y4;
          Sow[{t, y}]
        ];
        hNew = 0.84 (tol / (err + $MachineEpsilon))^(1/4) h;
        h = Clip[hNew, {h/5, 5h}]
      ]
    ];
    (* raw is either {} (no steps accepted) or a list of {t,y} pairs *)
    With[{pairs = If[raw === {}, {{t0, y0}}, First[raw]]},
      ivpResult[pairs[[All, 1]], pairs[[All, 2]]]
    ]
  ]

(* ── RK4 startup for multi-step methods ────────────────────────────────── *)
ivpRK4StartUp[f_, t0_, nSteps_, y0_, h_] :=
  ivpRungeKutta4[f, t0, t0 + nSteps h, y0, h]

(* ══════════════════════════════════════════════════════════════════════════
   Multi-step methods — pre-allocated ts / ys / fvals
   Startup (ivpRK4StartUp) produces 4 initial points (indices 1..nStart).
   Subsequent steps fill in forward using positive indexing i, i-1, i-2, i-3.
   This eliminates negative-tail indexing (fvals[[-1]] etc.) which was
   fragile whenever startup size changed.
   A short-span guard falls back to RK4 when the interval holds fewer than
   3 steps, making it impossible to build the required startup history.
   ══════════════════════════════════════════════════════════════════════════ *)

(* Helper: multi-step short-interval check *)
ivpMultiStepCheck[method_, t0_, tEnd_, y0_, h_, fallback_] :=
  If[(tEnd - t0) < 3 h,
    Message[SolveInitialValueProblem::shorttspan, method];
    fallback,
    Null   (* caller proceeds *)
  ]

(* ── Adams-Bashforth 4-step Explicit (order 4) ──────────────────────────── *)
ivpAdamsBashforth[f_, t0_, tEnd_, y0_, h_] :=
  Module[{check, startup, nStart, nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, fvals, t, i, y},
    check = ivpMultiStepCheck["AdamsBashforth", t0, tEnd, y0, h,
              ivpRungeKutta4[f, t0, tEnd, y0, h]];
    If[check =!= Null, Return[check, Module]];
    startup = ivpRK4StartUp[f, t0, 3, y0, h];
    nStart  = Length[startup["Grid"]];
    ts    = ConstantArray[0.,         nMax];
    ys    = ConstantArray[y0,         nMax];
    fvals = ConstantArray[f[t0, y0],  nMax];
    ts[[1;;nStart]]    = startup["Grid"];
    ys[[1;;nStart]]    = startup["Solution"];
    fvals[[1;;nStart]] = MapThread[f, {startup["Grid"], startup["Solution"]}];
    i = nStart; t = ts[[i]];
    While[t + h/2 <= tEnd,
      y = ys[[i]] + h/24 (55 fvals[[i]] - 59 fvals[[i-1]] + 37 fvals[[i-2]] - 9 fvals[[i-3]]);
      t += h; i++;
      ts[[i]] = t; ys[[i]] = y; fvals[[i]] = f[t, y]
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── Adams 4th-Order Predictor-Corrector (order 4) ─────────────────────── *)
ivpAdamsPC4[f_, t0_, tEnd_, y0_, h_] :=
  Module[{check, startup, nStart, nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, fvals, t, i, y, yPred, fPred},
    check = ivpMultiStepCheck["AdamsPC4", t0, tEnd, y0, h,
              ivpRungeKutta4[f, t0, tEnd, y0, h]];
    If[check =!= Null, Return[check, Module]];
    startup = ivpRK4StartUp[f, t0, 3, y0, h];
    nStart  = Length[startup["Grid"]];
    ts    = ConstantArray[0.,         nMax];
    ys    = ConstantArray[y0,         nMax];
    fvals = ConstantArray[f[t0, y0],  nMax];
    ts[[1;;nStart]]    = startup["Grid"];
    ys[[1;;nStart]]    = startup["Solution"];
    fvals[[1;;nStart]] = MapThread[f, {startup["Grid"], startup["Solution"]}];
    i = nStart; t = ts[[i]];
    While[t + h/2 <= tEnd,
      (* Adams-Bashforth predictor *)
      yPred = ys[[i]] + h/24 (55 fvals[[i]] - 59 fvals[[i-1]] + 37 fvals[[i-2]] - 9 fvals[[i-3]]);
      fPred = f[t+h, yPred];
      (* Adams-Moulton corrector *)
      y = ys[[i]] + h/24 (9 fPred + 19 fvals[[i]] - 5 fvals[[i-1]] + fvals[[i-2]]);
      t += h; i++;
      ts[[i]] = t; ys[[i]] = y; fvals[[i]] = f[t, y]
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── Hamming Predictor-Corrector ────────────────────────────────────────── *)
ivpHammingPC[f_, t0_, tEnd_, y0_, h_] :=
  Module[{check, startup, nStart, nMax = Ceiling[(tEnd-t0)/h] + 2,
          ts, ys, fvals, t, i, y, yPred, fPred, yCorr},
    check = ivpMultiStepCheck["HammingPC", t0, tEnd, y0, h,
              ivpRungeKutta4[f, t0, tEnd, y0, h]];
    If[check =!= Null, Return[check, Module]];
    startup = ivpRK4StartUp[f, t0, 3, y0, h];
    nStart  = Length[startup["Grid"]];
    ts    = ConstantArray[0.,         nMax];
    ys    = ConstantArray[y0,         nMax];
    fvals = ConstantArray[f[t0, y0],  nMax];
    ts[[1;;nStart]]    = startup["Grid"];
    ys[[1;;nStart]]    = startup["Solution"];
    fvals[[1;;nStart]] = MapThread[f, {startup["Grid"], startup["Solution"]}];
    i = nStart; t = ts[[i]];
    While[t + h/2 <= tEnd,
      (* Milne open predictor (needs ys 3 steps back: index i-3) *)
      yPred = ys[[i-3]] + 4h/3 (2 fvals[[i]] - fvals[[i-1]] + 2 fvals[[i-2]]);
      fPred = f[t+h, yPred];
      (* Hamming corrector (needs ys 1 and 2 steps back: i and i-2) *)
      yCorr = (9 ys[[i]] - ys[[i-2]])/8 + 3h/8 (fPred + 2 fvals[[i]] - fvals[[i-1]]);
      t += h; i++;
      ts[[i]] = t; ys[[i]] = yCorr; fvals[[i]] = f[t, yCorr]
    ];
    ivpResult[ts[[;;i]], ys[[;;i]]]
  ]

(* ── 4.3  Protect cross-module exports ──────────────────────────────────── *)
Protect[ivpRungeKutta4, ivpRK4StartUp, ivpMultiStepCheck]
