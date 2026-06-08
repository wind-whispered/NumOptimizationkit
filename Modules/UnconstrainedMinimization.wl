(* UnconstrainedMinimization.wl -- Unconstrained extremum solvers
   Univariate : GoldenSection, Fibonacci, QuadraticInterpolation, Newton
   Multivariate: GradientDescent, ConjugateGradient{PR,FR}, Newton,
                 BFGS, DFP, NelderMead, SimulatedAnnealing, GeneticAlgorithm

   Batch-2 improvements
     * FindMinimum1D: validates a < b and that f evaluates to a real number
     * FindMinimumND: validates x0 is numeric and f[x0] is real
     * ConvergenceHistory option: when True, result Association includes
       "PointHistory" (list of iterates) and "ValueHistory" (list of f values).
       Implemented via tagged Sow/Reap; inner methods accept optional track_ flag.

   Private helper naming convention: camelCase prefixed by scope
     min1D*  -- univariate minimisation helpers
     minND*  -- multivariate minimisation helpers
*)

Options[FindMinimum1D] = {
  Method             -> "GoldenSection",
  Tolerance          -> 1*^-8,
  MaxIterations      -> 200
}

Options[FindMinimumND] = {
  Method             -> "BFGS",
  Tolerance          -> 1*^-8,
  MaxIterations      -> 500,
  Gradient           -> Automatic,
  ConvergenceHistory -> False,
  Compiled           -> Automatic   (* Batch-4: Automatic tries to Compile f *)
}

(* ── FindMinimum1D dispatcher ──────────────────────────────────────────── *)

FindMinimum1D[f_, {a_?NumericQ, b_?NumericQ}, opts:OptionsPattern[]] :=
  Module[{method, tol, maxIter, fa},
    (* Validation *)
    If[N@a >= N@b,
      Message[FindMinimum1D::badinterval, a, b]; Return[$Failed, Module]];
    fa = Quiet[f[N@a]];
    If[!NumericQ[fa] || !RealValuedQ[N@fa],
      Message[FindMinimum1D::badfunc]; Return[$Failed, Module]];
    method  = OptionValue[Method];
    tol     = OptionValue[Tolerance];
    maxIter = OptionValue[MaxIterations];
    Switch[method,
      "GoldenSection",          min1DGoldenSection[f, N@a, N@b, tol, maxIter],
      "Fibonacci",              min1DFibonacci[f, N@a, N@b, tol, maxIter],
      "QuadraticInterpolation", min1DQuadInterp[f, N@a, N@b, tol, maxIter],
      "Newton",                 min1DNewton[f, N@(a+b)/2, tol, maxIter],
      _,
        Message[FindMinimum1D::badmethod, method];
        min1DGoldenSection[f, N@a, N@b, tol, maxIter]
    ]
  ]

(* ── Golden-Section Search ──────────────────────────────────────────────── *)
min1DGoldenSection[f_, a0_, b0_, tol_, maxIter_] :=
  Module[{a = a0, b = b0, r = (Sqrt[5] - 1)/2, c, d, fc, fd, k = 0},
    c = a + (1 - r)(b - a); d = a + r(b - a);
    fc = f[c]; fd = f[d];
    While[(b - a) > tol && k++ < maxIter,
      If[fc < fd,
        b = d; d = c; fd = fc;
        c = a + (1 - r)(b - a); fc = f[c],
        a = c; c = d; fc = fd;
        d = a + r(b - a); fd = f[d]
      ]
    ];
    If[fc <= fd,
      <|"Point" -> c, "Value" -> fc, "Iterations" -> k|>,
      <|"Point" -> d, "Value" -> fd, "Iterations" -> k|>
    ]
  ]

(* ── Fibonacci Search ───────────────────────────────────────────────────── *)
min1DFibonacci[f_, a0_, b0_, tol_, maxIter_] :=
  Module[{a = a0, b = b0, fibs, n, k = 0, c, d, fc, fd},
    fibs = NestList[{#[[2]], #[[1]] + #[[2]]}&, {1, 1}, maxIter + 1][[All, 1]];
    n = SelectFirst[Range[Length[fibs]], fibs[[#]] > (b - a)/tol &, Length[fibs]];
    c = a + fibs[[n - 2]]/fibs[[n]] (b - a);
    d = a + fibs[[n - 1]]/fibs[[n]] (b - a);
    fc = f[c]; fd = f[d];
    While[k < n - 2,
      k++;
      If[fc < fd,
        b = d; d = c; fd = fc;
        c = a + fibs[[n - k - 2]]/fibs[[n - k]] (b - a); fc = f[c],
        a = c; c = d; fc = fd;
        d = a + fibs[[n - k - 1]]/fibs[[n - k]] (b - a); fd = f[d]
      ]
    ];
    If[fc <= fd,
      <|"Point" -> c, "Value" -> fc, "Iterations" -> k|>,
      <|"Point" -> d, "Value" -> fd, "Iterations" -> k|>
    ]
  ]

(* ── Quadratic Interpolation Search ────────────────────────────────────── *)
min1DQuadInterp[f_, a0_, b0_, tol_, maxIter_] :=
  Module[{a = a0, b = b0, m = (a0+b0)/2, fa, fb, fm, xp, denom, k = 0},
    fa = f[a]; fb = f[b]; fm = f[m];
    Catch[
      While[k++ < maxIter,
        denom = 2((m-a)(fb-fm) - (m-b)(fm-fa));
        If[Abs[denom] < $MachineEpsilon,
          Throw[<|"Point" -> m, "Value" -> fm, "Iterations" -> k|>]];
        xp = m - ((m-a)^2 (fb-fm) - (m-b)^2 (fm-fa)) / denom;
        xp = Clip[xp, {a, b}];
        If[Abs[xp - m] < tol,
          Throw[<|"Point" -> xp, "Value" -> f[xp], "Iterations" -> k|>]];
        Which[
          fa == Max[fa, fm, fb], a = xp; fa = f[xp],
          fb == Max[fa, fm, fb], b = xp; fb = f[xp],
          True,                   m = xp; fm = f[xp]
        ];
        m = (a + b)/2; fm = f[m]
      ];
      <|"Point" -> m, "Value" -> fm, "Iterations" -> k|>
    ]
  ]

(* ── Newton's Method for Minimisation (solves f' = 0) ──────────────────── *)
min1DNewton[f_, x0_, tol_, maxIter_] :=
  Module[{x = N@x0, fp, fpp, dx, k = 0},
    Catch[
      While[k++ < maxIter,
        fp  = numDeriv[f, x, 1];
        fpp = numDeriv[f, x, 2];
        If[Abs[fpp] < $MachineEpsilon,
          Throw[<|"Point" -> x, "Value" -> f[x], "Iterations" -> k|>]];
        dx = -fp/fpp;
        x += dx;
        If[Abs[dx] < tol,
          Throw[<|"Point" -> x, "Value" -> f[x], "Iterations" -> k|>]]
      ];
      <|"Point" -> x, "Value" -> f[x], "Iterations" -> k|>
    ]
  ]

(* ══════════════════════════════════════════════════════════════════════════
   FindMinimumND dispatcher
   ConvergenceHistory: inner methods accept track_:False; when True they
   Sow {x, f[x]} pairs tagged "nkhist" once per major iteration.
   The dispatcher wraps with Reap and appends history to the result.
   ══════════════════════════════════════════════════════════════════════════ *)

FindMinimumND[f_, x0_?VectorQ, opts:OptionsPattern[]] :=
  Module[{method, tol, maxIter, gradFn, track, compiled, x0r, fx0, raw, result, fUse},
    (* Validation *)
    x0r = Quiet[N@x0];
    If[!numericVectorQ[x0r],
      Message[FindMinimumND::badx0]; Return[$Failed, Module]];
    fx0 = Quiet[f[x0r]];
    If[!NumericQ[fx0] || !RealValuedQ[N@fx0],
      Message[FindMinimumND::badfunc]; Return[$Failed, Module]];
    method   = OptionValue[Method];
    tol      = OptionValue[Tolerance];
    maxIter  = OptionValue[MaxIterations];
    track    = TrueQ[OptionValue[ConvergenceHistory]];
    compiled = OptionValue[Compiled];
    (* 4.2  Attempt compilation of f for numerical methods that call f many times *)
    fUse = If[(TrueQ[compiled] || compiled === Automatic) &&
              MemberQ[{"BFGS","DFP","GradientDescent","ConjugateGradientPR","ConjugateGradientFR","NelderMead"}, method],
      tryCompileND[f, x0r],
      f
    ];
    gradFn  = With[{g = OptionValue[Gradient]},
                If[g === Automatic, numGrad[fUse, #]&, g]];
    If[track,
      raw = Reap[
        result = Switch[method,
          "GradientDescent",     minNDGradDesc[f, gradFn, x0r, tol, maxIter, True],
          "ConjugateGradientPR", minNDConjGrad[f, gradFn, x0r, tol, maxIter, "PR", True],
          "ConjugateGradientFR", minNDConjGrad[fUse, gradFn, x0r, tol, maxIter, "FR", True],
          "Newton",              minNDNewton[fUse, gradFn, x0r, tol, maxIter, True],
          "BFGS",                minNDQuasiNewton[fUse, gradFn, x0r, tol, maxIter, "BFGS", True],
          "DFP",                 minNDQuasiNewton[fUse, gradFn, x0r, tol, maxIter, "DFP", True],
          "NelderMead",          minNDNelderMead[fUse, x0r, tol, maxIter, True],
          "SimulatedAnnealing",  minNDSimAnnealing[f, x0r, tol, maxIter, True],  (* no compile benefit *)
          "GeneticAlgorithm",    minNDGenetic[f, x0r, tol, maxIter, True],
          _,
            Message[FindMinimumND::badmethod, method];
            minNDQuasiNewton[fUse, gradFn, x0r, tol, maxIter, "BFGS", True]
        ]
      , "nkhist"];
      With[{hist = If[Last[raw] === {}, {}, First[Last[raw]]]},
        (* Append only takes a single element; merge both new keys with Join *)
        Join[result,
          <|"PointHistory" -> hist[[All, 1]],
            "ValueHistory" -> hist[[All, 2]]|>]
      ],
      (* No history *)
      Switch[method,
        "GradientDescent",     minNDGradDesc[fUse, gradFn, x0r, tol, maxIter],
        "ConjugateGradientPR", minNDConjGrad[fUse, gradFn, x0r, tol, maxIter, "PR"],
        "ConjugateGradientFR", minNDConjGrad[fUse, gradFn, x0r, tol, maxIter, "FR"],
        "Newton",              minNDNewton[fUse, gradFn, x0r, tol, maxIter],
        "BFGS",                minNDQuasiNewton[fUse, gradFn, x0r, tol, maxIter, "BFGS"],
        "DFP",                 minNDQuasiNewton[fUse, gradFn, x0r, tol, maxIter, "DFP"],
        "NelderMead",          minNDNelderMead[fUse, x0r, tol, maxIter],
        "SimulatedAnnealing",  minNDSimAnnealing[f, x0r, tol, maxIter],
        "GeneticAlgorithm",    minNDGenetic[f, x0r, tol, maxIter],
        _,
          Message[FindMinimumND::badmethod, method];
          minNDQuasiNewton[fUse, gradFn, x0r, tol, maxIter, "BFGS"]
      ]
    ]
  ]

(* ── Gradient Descent ───────────────────────────────────────────────────── *)
minNDGradDesc[f_, gradFn_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = x0, xp, fx, g, res, k = 0, converged = False},
    fx = f[x];
    If[track, Sow[{x, fx}, "nkhist"]];
    Catch[
      While[k++ < maxIter,
        g = gradFn[x];
        If[Norm[g] < tol, converged = True; Throw[Null]];
        g = g / Norm[g];
        res = quadLineSearch[f, x, g, -g];
        If[!res[[4]], converged = False; Throw[Null]];
        xp = res[[1]];
        If[Norm[xp - x] < tol && Abs[res[[2]] - fx] < tol,
          converged = True; Throw[Null]];
        x = xp; fx = res[[2]];
        If[track, Sow[{x, fx}, "nkhist"]]
      ]
    ];
    <|"Point" -> x, "Value" -> f[x], "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Conjugate Gradient (PR or FR) ─────────────────────────────────────── *)
minNDConjGrad[f_, gradFn_, x0_, tol_, maxIter_, variant_, track_:False] :=
  Module[{x = x0, fx = f[x0], g, gp, s, xp, res, k = 0, beta, converged = False},
    g = gradFn[x]; s = -g;
    If[track, Sow[{x, fx}, "nkhist"]];
    Catch[
      While[k++ < maxIter,
        res = quadLineSearch[f, x, g, s];
        If[!res[[4]], converged = False; Throw[Null]];
        xp = res[[1]];
        If[Norm[xp - x] < tol && Abs[res[[2]] - fx] < tol,
          converged = True; Throw[Null]];
        gp = gradFn[xp];
        beta = If[variant === "PR",
          Max[0., (gp - g) . gp / (g . g + $MachineEpsilon)],
          gp . gp / (g . g + $MachineEpsilon)
        ];
        s = -gp + beta s;
        x = xp; fx = res[[2]]; g = gp;
        If[track, Sow[{x, fx}, "nkhist"]]
      ]
    ];
    <|"Point" -> x, "Value" -> f[x], "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Newton's Method (multivariate) ────────────────────────────────────── *)
minNDNewton[f_, gradFn_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = x0, g, H, d, k = 0, converged = False},
    If[track, Sow[{x, f[x]}, "nkhist"]];
    Catch[
      While[k++ < maxIter,
        g = gradFn[x];
        If[Norm[g] < tol, converged = True; Throw[Null]];
        H = numHessian[f, x];
        d = safeLinearSolve[H, -g];
        If[d === $Failed, converged = False; Throw[Null]];
        x += d;
        If[track, Sow[{x, f[x]}, "nkhist"]];
        If[Norm[d] < tol, converged = True; Throw[Null]]
      ]
    ];
    <|"Point" -> x, "Value" -> f[x], "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── BFGS / DFP Quasi-Newton ────────────────────────────────────────────── *)
minNDQuasiNewton[f_, gradFn_, x0_, tol_, maxIter_, variant_, track_:False] :=
  Module[{x = x0, n = Length[x0], H = IdentityMatrix[Length[x0]]*1.,
          E = IdentityMatrix[Length[x0]], g, d, res, xp, gp, s, y, sy,
          k = 0, converged = False},
    g = gradFn[x];
    If[track, Sow[{x, f[x]}, "nkhist"]];
    Catch[
      While[k++ < maxIter,
        If[Norm[g] < tol, converged = True; Throw[Null]];
        d = -H . g;
        res = quadLineSearch[f, x, g, d];
        If[!res[[4]], H = E; res = quadLineSearch[f, x, g, d]];
        xp = res[[1]]; gp = gradFn[xp];
        s = xp - x; y = gp - g; sy = s . y;
        If[Abs[sy] > $MachineEpsilon,
          H = Switch[variant,
            "BFGS",
              (E - Outer[Times, s, y]/sy) . H . (E - Outer[Times, y, s]/sy) +
              Outer[Times, s, s]/sy,
            "DFP",
              H + Outer[Times, s, s]/sy -
              (H . Outer[Times, y, y] . H) / (y . H . y)
          ]
        ];
        x = xp; g = gp;
        If[track, Sow[{x, res[[2]]}, "nkhist"]]
      ]
    ];
    <|"Point" -> x, "Value" -> f[x], "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Nelder-Mead Simplex ────────────────────────────────────────────────── *)
minNDNelderMead[f_, x0_, tol_, maxIter_, track_:False] :=
  Module[{n = Length[x0], simplex, fvals, order, best, worst, centroid,
          xr, xe, xc, fr, fe, fc, k = 0, converged = False},
    simplex = Prepend[Table[x0 + IdentityMatrix[n][[i]], {i, n}], x0];
    fvals   = f /@ simplex;
    If[track, Sow[{simplex[[1]], fvals[[1]]}, "nkhist"]];
    Catch[
      While[k++ < maxIter,
        order   = Ordering[fvals];
        simplex = simplex[[order]]; fvals = fvals[[order]];
        best    = simplex[[1]]; worst = simplex[[-1]];
        If[Norm[worst - best] < tol && Abs[fvals[[-1]] - fvals[[1]]] < tol,
          converged = True; Throw[Null]];
        centroid = Mean[Most[simplex]];
        xr = 2 centroid - worst; fr = f[xr];
        If[fr < fvals[[1]],
          xe = 3 centroid - 2 worst; fe = f[xe];
          {simplex[[-1]], fvals[[-1]]} = If[fe < fr, {xe, fe}, {xr, fr}],
          If[fr < fvals[[-2]],
            simplex[[-1]] = xr; fvals[[-1]] = fr,
            xc = (centroid + worst)/2; fc = f[xc];
            If[fc < fvals[[-1]],
              simplex[[-1]] = xc; fvals[[-1]] = fc,
              simplex = Join[{best}, (best + (# - best)/2 &) /@ Rest[simplex]];
              fvals = f /@ simplex
            ]
          ]
        ];
        If[track, Sow[{simplex[[Ordering[fvals, 1][[1]]]], Min[fvals]}, "nkhist"]]
      ]
    ];
    order = Ordering[fvals];
    <|"Point" -> simplex[[order[[1]]]], "Value" -> fvals[[order[[1]]]],
      "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Simulated Annealing ────────────────────────────────────────────────── *)
minNDSimAnnealing[f_, x0_, tol_, maxIter_, track_:False] :=
  Module[{x = x0, xBest = x0, fx = f[x0], fBest = f[x0],
          T0 = 100., Tmin = 1*^-8, cool = 0.95, stepSize = 0.5,
          xNew, fNew, dE, k = 0, converged = False},
    SeedRandom[42];
    If[track, Sow[{x, fx}, "nkhist"]];
    While[k++ < maxIter && (T0 > Tmin || Norm[x - xBest] > tol),
      xNew = x + stepSize RandomVariate[NormalDistribution[], Length[x]];
      fNew = f[xNew]; dE = fNew - fx;
      (* As T0 cools toward Tmin, Exp[-dE/T0] correctly underflows to 0.
         (rejecting the worse move) -- silence the resulting General::munfl *)
      If[dE < 0 || RandomReal[] < Quiet[Exp[-dE/T0], General::munfl],
        x = xNew; fx = fNew;
        If[fx < fBest, xBest = x; fBest = fx]
      ];
      T0 *= cool;
      If[track, Sow[{xBest, fBest}, "nkhist"]]
    ];
    converged = Norm[numGrad[f, xBest]] < tol;
    <|"Point" -> xBest, "Value" -> fBest,
      "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Genetic Algorithm ──────────────────────────────────────────────────── *)
minNDGenetic[f_, x0_, tol_, maxIter_, track_:False] :=
  Module[{n = Length[x0], popSize = 30, elite = 4, bound = 5.,
          pop, fit, order, parents, child, mutRate = 0.1, mutSigma = 0.5, k = 0},
    SeedRandom[123];
    pop = Join[{x0}, Table[x0 + bound RandomVariate[NormalDistribution[], n], {popSize - 1}]];
    While[k++ < maxIter,
      fit   = f /@ pop; order = Ordering[fit];
      pop   = pop[[order]]; fit = fit[[order]];
      If[track, Sow[{pop[[1]], fit[[1]]}, "nkhist"]];
      If[Norm[numGrad[f, pop[[1]]]] < tol, Break[]];
      pop = Join[Take[pop, elite],
        Table[
          parents = RandomSample[Range[popSize], 2];
          child   = 0.5 (pop[[parents[[1]]]] + pop[[parents[[2]]]]);
          child + If[RandomReal[] < mutRate, mutSigma RandomVariate[NormalDistribution[], n], 0 x0],
          {popSize - elite}
        ]
      ]
    ];
    fit = f /@ pop; order = Ordering[fit];
    <|"Point" -> pop[[order[[1]]]], "Value" -> fit[[order[[1]]]],
      "Iterations" -> k,
      "Converged"  -> Norm[numGrad[f, pop[[order[[1]]]]]] < tol|>
  ]
