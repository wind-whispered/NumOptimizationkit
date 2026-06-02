(* EquationRootFinding.wl -- Univariate equation root-finding
   Bracket methods : Bisection, RegulaFalsi, Brent, Muller
   Open methods    : Newton, Halley, FixedPointIter, Steffensen, Aitken
   Two-point       : Secant

   Batch-2: input validation, ConvergenceHistory option
   Batch-4: WorkingPrecision option -- inner methods use nkEps[] for all
            convergence and division-guard thresholds; useful for high-precision
            root finding when f returns exact symbolic values.

   Private helper naming convention: camelCase prefixed by "root"
   Note: "FixedPoint" is a protected Mathematica symbol; helper is rootFixedPointIter.
*)

Options[FindEquationRoot] = {
  Method             -> "Bisection",
  Tolerance          -> 1*^-8,
  MaxIterations      -> 200,
  ConvergenceHistory -> False,
  WorkingPrecision   -> MachinePrecision  (* Batch-4 *)
}

FindEquationRoot[f_, {a_?NumericQ, b_?NumericQ}, opts:OptionsPattern[]] :=
  Module[{method, tol, maxIter, track, wp, fa, raw, result},
    wp  = OptionValue[WorkingPrecision];
    (* Validation: f must evaluate at a *)
    fa = Quiet[f[nkN[a]]];
    If[!NumericQ[fa],
      Message[FindEquationRoot::badfunc]; Return[$Failed, Module]];
    method  = OptionValue[Method];
    tol     = OptionValue[Tolerance];
    maxIter = OptionValue[MaxIterations];
    track   = TrueQ[OptionValue[ConvergenceHistory]];
    (* 4.1  Set session-level working precision for inner methods *)
    Block[{nkWP = wp},
    If[track,
      raw = Reap[
        result = Switch[method,
          "Bisection",   rootBisection[f, N@a, N@b, tol, maxIter, True],
          "RegulaFalsi", rootRegulaFalsi[f, N@a, N@b, tol, maxIter, True],
          "Brent",       rootBrent[f, N@a, N@b, tol, maxIter, True],
          "Newton",      rootNewton[f, N@a, tol, maxIter, True],
          "Secant",      rootSecant[f, N@a, N@b, tol, maxIter, True],
          "Halley",      rootHalley[f, N@a, tol, maxIter, True],
          "Muller",      rootMuller[f, N@a, N@b, N@(a+b)/2, tol, maxIter, True],
          "FixedPoint",  rootFixedPointIter[f, N@a, tol, maxIter, True],
          "Steffensen",  rootSteffensen[f, N@a, tol, maxIter, True],
          "Aitken",      rootAitken[f, N@a, tol, maxIter, True],
          _,
            Message[FindEquationRoot::badmethod, method];
            rootBisection[f, N@a, N@b, tol, maxIter, True]
        ]
      , "nkhist"];
      If[result === $Failed, Return[$Failed, Module]];
      Append[result,
        "ResidualHistory" -> If[Last[raw] === {}, {}, First[Last[raw]]]],
      (* No history *)
      result = Switch[method,
        "Bisection",   rootBisection[f, N@a, N@b, tol, maxIter],
        "RegulaFalsi", rootRegulaFalsi[f, N@a, N@b, tol, maxIter],
        "Brent",       rootBrent[f, N@a, N@b, tol, maxIter],
        "Newton",      rootNewton[f, N@a, tol, maxIter],
        "Secant",      rootSecant[f, N@a, N@b, tol, maxIter],
        "Halley",      rootHalley[f, N@a, tol, maxIter],
        "Muller",      rootMuller[f, N@a, N@b, N@(a+b)/2, tol, maxIter],
        "FixedPoint",  rootFixedPointIter[f, N@a, tol, maxIter],
        "Steffensen",  rootSteffensen[f, N@a, tol, maxIter],
        "Aitken",      rootAitken[f, N@a, tol, maxIter],
        _,
          Message[FindEquationRoot::badmethod, method];
          rootBisection[f, N@a, N@b, tol, maxIter]
      ];
      result
    ]   (* end If[track] *)
    ]   (* end Block[{nkWP}] *)
  ]     (* end Module *)

(* 鈹€鈹€ Bisection Method 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootBisection[f_, a0_, b0_, tol_, maxIter_, track_:False] :=
  Module[{a = a0, b = b0, c = a0, fa = f[a0], fb = f[b0], fc, err, k = 0,
          converged = False},
    If[fa fb > 0, Message[FindEquationRoot::bracket]; Return[$Failed, Module]];
    If[track, Sow[Abs[fa], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        c = (a + b)/2; fc = f[c]; err = (b - a)/2;
        If[track, Sow[Abs[fc], "nkhist"]];
        If[Abs[fc] < nkEps[] || err < tol,
          converged = True; Throw[Null]];
        If[fc fa > 0, a = c; fa = fc, b = c]
      ]
    ];
    <|"Root" -> c, "Residual" -> Abs[fc], "Iterations" -> k, "Converged" -> converged|>
  ]

(* 鈹€鈹€ Regula Falsi 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootRegulaFalsi[f_, a0_, b0_, tol_, maxIter_, track_:False] :=
  Module[{a = a0, b = b0, c = a0, fa = f[a0], fb = f[b0], fc, k = 0,
          converged = False},
    If[fa fb > 0, Message[FindEquationRoot::bracket]; Return[$Failed, Module]];
    If[track, Sow[Abs[fa], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        c  = (a fb - b fa) / (fb - fa);
        fc = f[c];
        If[track, Sow[Abs[fc], "nkhist"]];
        If[Abs[fc] < tol, converged = True; Throw[Null]];
        If[fc fa > 0, a = c; fa = fc, b = c; fb = fc]
      ]
    ];
    <|"Root" -> c, "Residual" -> Abs[fc], "Iterations" -> k, "Converged" -> converged|>
  ]

(* 鈹€鈹€ Brent's Method 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootBrent[f_, a0_, b0_, tol_, maxIter_, track_:False] :=
  Module[{a = a0, b = b0, c = a0, d, s, fa, fb, fc, mflag = True, k = 0,
          converged = False},
    fa = f[a]; fb = f[b];
    If[fa fb > 0, Message[FindEquationRoot::bracket]; Return[$Failed, Module]];
    If[Abs[fa] < Abs[fb], {a, b} = {b, a}; {fa, fb} = {fb, fa}];
    c = a; fc = fa;
    If[track, Sow[Abs[fb], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        If[Abs[fb] < tol, converged = True; Throw[Null]];
        If[fa != fc && fb != fc,
          s = a fb fc/((fa-fb)(fa-fc)) + b fa fc/((fb-fa)(fb-fc)) + c fa fb/((fc-fa)(fc-fb)),
          s = b - fb (b-a)/(fb-fa)
        ];
        If[! (((3a+b)/4 < s < b || b < s < (3a+b)/4) &&
              !(mflag  && Abs[s-b] < Abs[b-c]/2) &&
              !(!mflag && Abs[s-b] < Abs[c-d]/2) &&
              !(mflag  && Abs[b-c] > tol) &&
              !(!mflag && Abs[c-d] > tol)),
          s = (a+b)/2; mflag = True, mflag = False
        ];
        d = c; c = b; fc = fb;
        If[fa f[s] < 0, b = s; fb = f[s], a = s; fa = f[s]];
        If[Abs[fa] < Abs[fb], {a, b} = {b, a}; {fa, fb} = {fb, fa}];
        If[track, Sow[Abs[fb], "nkhist"]]
      ]
    ];
    <|"Root" -> b, "Residual" -> Abs[fb], "Iterations" -> k, "Converged" -> converged|>
  ]

(* 鈹€鈹€ Newton-Raphson 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootNewton[f_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = N@x0, fx, dfx, err, k = 0, converged = False},
    If[track, Sow[Abs[f[x]], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        fx  = f[x]; dfx = numDeriv[f, x, 1];
        If[Abs[dfx] < nkEps[], Throw[Null]];
        x  -= fx/dfx;
        err = Abs[fx/(Abs[x] + nkEps[])];
        If[track, Sow[Abs[f[x]], "nkhist"]];
        If[err < tol || Abs[f[x]] < nkEps[],
          converged = True; Throw[Null]]
      ]
    ];
    <|"Root" -> x, "Residual" -> Abs[f[x]], "Iterations" -> k, "Converged" -> converged|>
  ]

(* 鈹€鈹€ Secant Method 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootSecant[f_, x0_, x1_, tol_, maxIter_, track_:False] :=
  Module[{pts = {N@x0, N@x1}, fp = {f[N@x0], f[N@x1]}, x2, fx2, err,
          k = 0, converged = False},
    If[track, Sow[Abs[fp[[1]]], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        x2  = pts[[2]] - fp[[2]] (pts[[2]] - pts[[1]]) / (fp[[2]] - fp[[1]] + nkEps[]);
        fx2 = f[x2];
        err = Abs[x2 - pts[[2]]] / (Abs[x2] + nkEps[]);
        pts = {pts[[2]], x2}; fp = {fp[[2]], fx2};
        If[track, Sow[Abs[fx2], "nkhist"]];
        If[err < tol || Abs[fx2] < nkEps[], converged = True; Throw[Null]]
      ]
    ];
    <|"Root" -> pts[[2]], "Residual" -> Abs[fp[[2]]], "Iterations" -> k,
      "Converged" -> converged|>
  ]

(* 鈹€鈹€ Halley's Method 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootHalley[f_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = N@x0, fx, dfx, d2fx, dx, k = 0, converged = False},
    If[track, Sow[Abs[f[x]], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        fx   = f[x]; dfx = numDeriv[f, x, 1]; d2fx = numDeriv[f, x, 2];
        dx   = 2 fx dfx / (2 dfx^2 - fx d2fx + nkEps[]);
        x   -= dx;
        If[track, Sow[Abs[f[x]], "nkhist"]];
        If[Abs[dx] < tol, converged = True; Throw[Null]]
      ]
    ];
    <|"Root" -> x, "Residual" -> Abs[f[x]], "Iterations" -> k, "Converged" -> converged|>
  ]

(* 鈹€鈹€ Muller's Parabolic Method 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootMuller[f_, x0_, x2_, x1_, tol_, maxIter_, track_:False] :=
  Module[{pts = {N@x0, N@x1, N@x2}, fp = f /@ {N@x0, N@x1, N@x2},
          h0, h1, del0, del1, a, b, c, disc, xNew, k = 0, converged = False},
    If[track, Sow[Abs[fp[[3]]], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        h0   = pts[[2]] - pts[[1]]; h1 = pts[[3]] - pts[[2]];
        del0 = (fp[[2]] - fp[[1]])/h0; del1 = (fp[[3]] - fp[[2]])/h1;
        a    = (del1 - del0)/(h1 + h0); b = a h1 + del1; c = fp[[3]];
        disc = Sqrt[Max[0, b^2 - 4 a c]];
        xNew = pts[[3]] - 2c / (If[Abs[b+disc] > Abs[b-disc], b+disc, b-disc] + nkEps[]);
        pts  = {pts[[2]], pts[[3]], xNew};
        fp   = {fp[[2]],  fp[[3]],  f[xNew]};
        If[track, Sow[Abs[fp[[3]]], "nkhist"]];
        If[Abs[fp[[3]]] < tol, converged = True; Throw[Null]]
      ]
    ];
    <|"Root" -> pts[[3]], "Residual" -> Abs[fp[[3]]], "Iterations" -> k,
      "Converged" -> converged|>
  ]

(* 鈹€鈹€ Fixed-Point Iteration 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootFixedPointIter[g_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = N@x0, xp, k = 0, converged = False},
    If[track, Sow[Abs[g[x] - x], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        xp = g[x];
        If[track, Sow[Abs[xp - x], "nkhist"]];
        If[Abs[xp - x] < tol, converged = True; x = xp; Throw[Null]];
        x = xp
      ]
    ];
    <|"Root" -> x, "Residual" -> Abs[g[x] - x], "Iterations" -> k, "Converged" -> converged|>
  ]

(* 鈹€鈹€ Steffensen's Acceleration 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootSteffensen[f_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = N@x0, x1, x2, k = 0, converged = False},
    If[track, Sow[Abs[f[x]], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        x1 = x + f[x]; x2 = x1 + f[x1];
        x  = x - (x1 - x)^2 / (x2 - 2 x1 + x + nkEps[]);
        If[track, Sow[Abs[f[x]], "nkhist"]];
        If[Abs[f[x]] < tol, converged = True; Throw[Null]]
      ]
    ];
    <|"Root" -> x, "Residual" -> Abs[f[x]], "Iterations" -> k, "Converged" -> converged|>
  ]

(* 鈹€鈹€ Aitken Delta-Squared Acceleration 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
rootAitken[f_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = N@x0, x1, x2, xHat, k = 0, converged = False},
    If[track, Sow[Abs[f[x] - x], "nkhist"]];
    Catch[
      While[k++ < maxIter,
        x1   = f[x]; x2 = f[x1];
        xHat = x - (x1 - x)^2 / (x2 - 2 x1 + x + nkEps[]);
        If[track, Sow[Abs[f[xHat] - xHat], "nkhist"]];
        If[Abs[xHat - x] < tol, converged = True; x = xHat; Throw[Null]];
        x = xHat
      ]
    ];
    <|"Root" -> x, "Residual" -> Abs[f[x] - x], "Iterations" -> k, "Converged" -> converged|>
  ]

