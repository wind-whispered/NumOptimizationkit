(* ConstrainedOptimization.wl -- Constrained minimisation
   FindConstrainedMinimum[f, x0] minimises f subject to constraints
   specified through options.

   Phase A: Box constraints   (LowerBounds, UpperBounds options)
     Method "ProjectedGradient" -- gradient projection onto the feasible box
   Phase B: Equality constraints  (EqualityConstraints option)
     Method "AugmentedLagrangian" -- multiplier method; inner loop uses BFGS
   Combined: both box and equality constraints work together.

   Return value: same as FindMinimumND — Association with keys
     "Point", "Value", "Iterations", "Converged"

   Private helper naming convention: camelCase prefixed by "cnstr"
   Depends on: numGrad, quadLineSearch, minNDQuasiNewton (via direct call)
*)

Options[FindConstrainedMinimum] = {
  Method                   -> Automatic,  (* Automatic selects based on constraints *)
  LowerBounds              -> Automatic,  (* Automatic = -Infinity for all components *)
  UpperBounds              -> Automatic,  (* Automatic = +Infinity for all components *)
  EqualityConstraints      -> None,       (* list of h_i functions: h_i[x]=0 *)
  Tolerance                -> 1*^-8,
  MaxIterations            -> 1000,
  AugmentedLagrangianPenalty -> 10.      (* initial penalty mu *)
}

FindConstrainedMinimum::usage =
"FindConstrainedMinimum[f, x0] finds a local minimum of f subject to constraints
specified through options. f is called as f[{x1,...,xn}].

Options
  LowerBounds         -> Automatic | list  (box lower bounds; Automatic = -Inf)
  UpperBounds         -> Automatic | list  (box upper bounds; Automatic = +Inf)
  EqualityConstraints -> None | {h1, h2,...}  (equality constraints h_i[x]=0)
  Method              -> Automatic | \"ProjectedGradient\" | \"AugmentedLagrangian\"
  Tolerance           -> 1*^-8
  MaxIterations       -> 1000
  AugmentedLagrangianPenalty -> 10.

Auto-selection: if only box constraints -> ProjectedGradient;
                if equality constraints present -> AugmentedLagrangian."

FindConstrainedMinimum[f_, x0_?VectorQ, opts:OptionsPattern[]] :=
  Module[{method, lb, ub, eqs, tol, maxIter, mu, n, x0r,
          hasBox, hasEq},
    x0r    = N@x0;
    n      = Length[x0r];
    lb     = OptionValue[LowerBounds];
    ub     = OptionValue[UpperBounds];
    eqs    = OptionValue[EqualityConstraints];
    tol    = OptionValue[Tolerance];
    maxIter = OptionValue[MaxIterations];
    mu     = N[OptionValue[AugmentedLagrangianPenalty]];
    (* Normalise bounds: Automatic -> +-Infinity *)
    lb     = If[lb === Automatic, ConstantArray[-Infinity, n], N@lb];
    ub     = If[ub === Automatic, ConstantArray[+Infinity, n], N@ub];
    hasBox = Or @@ Table[lb[[i]] > -Infinity || ub[[i]] < Infinity, {i, n}];
    hasEq  = eqs =!= None && eqs =!= {};
    method = OptionValue[Method];
    If[method === Automatic,
      method = If[hasEq, "AugmentedLagrangian", "ProjectedGradient"]
    ];
    Switch[method,
      "ProjectedGradient",
        cnstrProjectedGradient[f, x0r, lb, ub, tol, maxIter],
      "AugmentedLagrangian",
        cnstrAugmentedLagrangian[f, x0r, lb, ub, eqs, mu, tol, maxIter],
      _,
        Message[FindConstrainedMinimum::badmethod, method];
        cnstrProjectedGradient[f, x0r, lb, ub, tol, maxIter]
    ]
  ]

FindConstrainedMinimum::badmethod =
  "Unknown Method \"`1`\". Falling back to \"ProjectedGradient\"."

(* ── Box projection ──────────────────────────────────────────────────────── *)
(* Project x onto the box [lb, ub] element-wise. *)
cnstrProjectBox[x_, lb_, ub_] :=
  MapThread[Clip[#1, {#2, #3}]&, {x, lb, ub}]

(* ── Projected Gradient Method (Phase A: box constraints only) ───────────── *)
(* Gradient projection onto feasible box with Armijo backtracking. *)
cnstrProjectedGradient[f_, x0_, lb_, ub_, tol_, maxIter_] :=
  Module[{x = cnstrProjectBox[x0, lb, ub], g, xNew, alpha, fx, fNew,
          k = 0, converged = False},
    fx = f[x];
    Catch[
      While[k++ < maxIter,
        g      = numGrad[f, x];
        (* Armijo backtracking along projected gradient direction *)
        alpha = 1.;
        Do[
          xNew = cnstrProjectBox[x - alpha g, lb, ub];
          fNew = f[xNew];
          If[fNew < fx - 1*^-4 alpha (g . (x - xNew)),
            Break[]   (* Armijo condition satisfied *)
          ];
          alpha *= 0.5,
          {30}
        ];
        If[Norm[xNew - x] < tol && Abs[fNew - fx] < tol,
          converged = True; Throw[Null]];
        x = xNew; fx = fNew
      ]
    ];
    <|"Point" -> x, "Value" -> f[x], "Iterations" -> k, "Converged" -> converged|>
  ]

(* ── Augmented Lagrangian Method (Phase B: equality + optional box) ──────── *)
(*   Minimise L_mu(x, lambda) = f(x) + lambda.h(x) + mu/2 ||h(x)||^2
     over the feasible box (projected BFGS for inner loop),
     then update lambda <- lambda + mu * h(x),
     optionally increase mu until convergence. *)
cnstrAugmentedLagrangian[f_, x0_, lb_, ub_, eqs_, mu0_, tol_, maxIter_] :=
  Module[{x = cnstrProjectBox[x0, lb, ub], n = Length[x0], m,
          lambda, mu = mu0, hVals, augF, result,
          k = 0, converged = False},
    m      = Length[eqs];
    lambda = ConstantArray[0., m];
    Catch[
      While[k++ < maxIter,
        hVals = Through[eqs[x]];   (* evaluate all constraint functions *)
        If[Norm[hVals] < tol, converged = True; Throw[Null]];
        (* Build augmented Lagrangian as a pure function *)
        With[{lam = lambda, muC = mu},
          augF = Function[xx,
            f[xx] + lam . Through[eqs[xx]] + muC/2 Norm[Through[eqs[xx]]]^2
          ]
        ];
        (* Inner minimisation with projected gradient (box constraints) *)
        result = cnstrProjectedGradient[augF, x, lb, ub, tol/10, 200];
        x = result["Point"];
        (* Multiplier update *)
        hVals  = Through[eqs[x]];
        lambda = lambda + mu hVals;
        (* Penalty increase for slow feasibility improvement *)
        If[Norm[hVals] > 0.25 Norm[hVals], mu *= 2.]
      ]
    ];
    <|"Point" -> x, "Value" -> f[x], "Iterations" -> k, "Converged" -> converged|>
  ]
