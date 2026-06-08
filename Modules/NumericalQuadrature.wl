(* NumericalQuadrature.wl -- Numerical Quadrature
   Composite Newton-Cotes: Trapezoidal, Simpson
   Adaptive:               AdaptiveSimpson (depth-controlled), Romberg
   Gauss rules:            GaussLegendre, GaussChebyshev, GaussHermite,
                           GaussLaguerre, GaussLobatto, GaussRadau

   Batch-1 improvements
     * Gauss node/weight tables memoized per n (computed once per session)
     * AdaptiveSimpson: explicit depth counter replaces unbounded recursion;
       MaxRecursionDepth option (default 60) prevents stack overflow
*)

Options[NumericalQuadrature] = {
  Method            -> "GaussLegendre",
  Points            -> 5,
  Intervals         -> 100,
  Tolerance         -> 1*^-6,
  MaxRecursionDepth -> 60,
  Weighted          -> True
}

NumericalQuadrature::maxdepth =
  "AdaptiveSimpson reached maximum recursion depth (`1`). \
Result may be inaccurate. Increase MaxRecursionDepth or switch to \"Romberg\"."

NumericalQuadrature::equalendpoints =
  "a == b: returning 0 (integral over an empty interval)."

NumericalQuadrature::weighted =
  "GaussChebyshev/GaussHermite/GaussLaguerre integrate weighted forms by default \
(see ::usage). Use Weighted -> False to integrate f[x] directly."

(* ── Primary form: NumericalQuadrature[f, {a, b}, opts] ─────────────────── *)
(* Batch-2: the {x,a,b} dummy-variable form is kept for backward compatibility;
   both forms dispatch to the same internal nqDispatch. *)

NumericalQuadrature[f_, {a_?NumericQ, b_?NumericQ}, opts:OptionsPattern[]] :=
  nqDispatch[f, N@a, N@b, opts]

(* Legacy form: NumericalQuadrature[f, {x, a, b}, opts] — x is ignored *)
NumericalQuadrature[f_, {_, a_?NumericQ, b_?NumericQ}, opts:OptionsPattern[]] :=
  nqDispatch[f, N@a, N@b, opts]

nqDispatch[f_, a_, b_, opts:OptionsPattern[NumericalQuadrature]] :=
  Module[{method, n, m, tol, maxDepth, weighted},
    If[a == b, Message[NumericalQuadrature::equalendpoints]; Return[0., Module]];
    method   = OptionValue[NumericalQuadrature, {opts}, Method];
    n        = OptionValue[NumericalQuadrature, {opts}, Points];
    m        = OptionValue[NumericalQuadrature, {opts}, Intervals];
    tol      = OptionValue[NumericalQuadrature, {opts}, Tolerance];
    maxDepth = OptionValue[NumericalQuadrature, {opts}, MaxRecursionDepth];
    weighted = OptionValue[NumericalQuadrature, {opts}, Weighted];
    Switch[method,
      "Trapezoidal",     quadTrapezoidal[f, a, b, m],
      "Simpson",         quadSimpson[f, a, b, m],
      "AdaptiveSimpson", quadAdaptiveSimpson[f, a, b, tol, maxDepth],
      "Romberg",         quadRomberg[f, a, b, tol],
      "GaussLegendre",   quadGaussLegendre[f, a, b, n],
      "GaussChebyshev",  quadGaussChebyshev[f, a, b, n, weighted],
      "GaussHermite",    quadGaussHermite[f, n, weighted],
      "GaussLaguerre",   quadGaussLaguerre[f, n, weighted],
      "GaussLobatto",    quadGaussLobatto[f, a, b, n],
      "GaussRadau",      quadGaussRadau[f, a, b, n],
      _,
        Message[NumericalQuadrature::badmethod, method];
        quadGaussLegendre[f, a, b, n]
    ]
  ]

(* ── Composite Trapezoidal Rule ─────────────────────────────────────────── *)
quadTrapezoidal[f_, a_, b_, m_] :=
  Module[{h = (b - a)/m, nodes = a + Range[0, m] (b - a)/m},
    h/2 (f[a] + 2 Total[f /@ nodes[[2 ;; -2]]] + f[b])
  ]

(* ── Composite Simpson's Rule ───────────────────────────────────────────── *)
quadSimpson[f_, a_, b_, m_] :=
  Module[{n = 2 Ceiling[m/2], h = (b - a)/(2 Ceiling[m/2]),
          nodes = a + Range[0, 2 Ceiling[m/2]] (b - a)/(2 Ceiling[m/2])},
    (* Use Range[...] (not a literal Span) so the "odd interior nodes" term
       gracefully evaluates to an empty selection when n == 2 (m == 1),
       rather than the invalid part-spec n-1 < 3 *)
    h/3 (f[a] + f[b] +
      4 Total[f /@ nodes[[Range[2, n, 2]]]] +
      2 Total[f /@ nodes[[Range[3, n-1, 2]]]])
  ]

(* ── Adaptive Simpson's Rule (depth-controlled) ─────────────────────────── *)
(* Entry point: seeds the depth counter and passes maxDepth through. *)
quadAdaptiveSimpson[f_, a_, b_, tol_, maxDepth_] :=
  quadAdaptHelper[f, a, b, tol, quadSimpson[f, a, b, 1], 0, maxDepth]

(* Recursive helper carrying current depth d and ceiling maxD. *)
quadAdaptHelper[f_, a_, b_, tol_, Sab_, d_, maxD_] :=
  Module[{c = (a + b)/2, Sac, Scb, S12},
    If[d >= maxD,
      Message[NumericalQuadrature::maxdepth, maxD];
      Return[Sab, Module]   (* return best estimate without further splitting *)
    ];
    Sac = quadSimpson[f, a, c, 1];
    Scb = quadSimpson[f, c, b, 1];
    S12 = Sac + Scb;
    If[Abs[S12 - Sab] / 15 < tol,
      S12 + (S12 - Sab) / 15,           (* Richardson correction *)
      quadAdaptHelper[f, a, c, tol/2, Sac, d+1, maxD] +
      quadAdaptHelper[f, c, b, tol/2, Scb, d+1, maxD]
    ]
  ]

(* ── Romberg Integration (Richardson extrapolation) ────────────────────── *)
quadRomberg[f_, a_, b_, tol_] :=
  Module[{T, h = b - a, m = 1, n = 0, maxRow = 25},
    T = ConstantArray[0., {maxRow + 1, 5}];
    T[[1, 1]] = h/2 (f[a] + f[b]);
    Catch[
      While[n++ < maxRow,
        h /= 2; m *= 2;
        T[[n+1, 1]] = T[[n, 1]]/2 + h Total[f /@ (a + (2 Range[m/2] - 1) h)];
        Do[
          T[[n+1, j+1]] = (4^j T[[n+1, j]] - T[[n, j]]) / (4^j - 1),
          {j, Min[n, 4]}
        ];
        If[n >= 4 && Abs[T[[n+1, Min[n+1, 5]]] - T[[n, Min[n, 5]]]] < tol,
          Throw[T[[n+1, Min[n+1, 5]]]]
        ]
      ];
      T[[n, Min[n, 5]]]
    ]
  ]

(* ══════════════════════════════════════════════════════════════════════════
   Gauss Quadrature — memoized node/weight tables
   Each gaussXxxNodesWeights[n] stores its result on first call;
   subsequent calls with the same n return the cached value instantly.
   ══════════════════════════════════════════════════════════════════════════ *)

(* ── Gauss-Legendre nodes and weights on [-1, 1] ────────────────────────── *)
gaussLegendreNodesWeights[n_Integer] := gaussLegendreNodesWeights[n] =
  Module[{k, x, xNew, Pm2, Pm1, P, Pp, itTol = 1*^-14,
          nodes = ConstantArray[0., n], weights = ConstantArray[0., n]},
    Do[
      x = Cos[Pi (k - 0.25) / (n + 0.5)];
      Do[
        Pm2 = 1.; Pm1 = x;
        Do[P = ((2j-1) x Pm1 - (j-1) Pm2)/j; Pm2 = Pm1; Pm1 = P, {j, 2, n}];
        Pp   = n (Pm2 - x Pm1) / (1 - x^2);
        xNew = x - Pm1/Pp;
        If[Abs[xNew - x] < itTol, x = xNew; Break[]];
        x = xNew,
        {50}
      ];
      nodes[[k]]   = x;
      weights[[k]] = 2 / ((1 - x^2) Pp^2),
      {k, n}
    ];
    {nodes, weights}
  ]

quadGaussLegendre[f_, a_, b_, n_] :=
  Module[{pts, wts},
    {pts, wts} = gaussLegendreNodesWeights[n];
    (b - a)/2 Total[wts * (f /@ ((b - a)/2 pts + (a + b)/2))]
  ]

(* ── Gauss-Chebyshev Type-I ──────────────────────────────────────────────── *)
(* Weighted -> True (default): integrates f(x)/Sqrt[1-x^2].
   Weighted -> False          : integrates f(x) directly by multiplying
                               each sample by Sqrt[1-t_k^2] (the weight). *)
quadGaussChebyshev[f_, a_, b_, n_, weighted_:True] :=
  Module[{t = Cos[(2 Range[n] - 1) Pi / (2n)], xPts, vals},
    xPts = (b - a)/2 t + (a + b)/2;
    vals = f /@ xPts;
    If[TrueQ[weighted],
      Pi/n (b - a)/2 Total[vals],
      Pi/n (b - a)/2 Total[Sqrt[1 - t^2] * vals]   (* cancel out weight *)
    ]
  ]

(* ── Gauss-Hermite nodes and weights (integrates f(x)*Exp[-x^2]) ────────── *)
(* Physicist's Hermite recurrence: H_0=1, H_1=2x,
   H_{k+1} = 2x H_k - 2k H_{k-1}
   Derivative: H_n'(x) = 2n H_{n-1}(x)
   Weight:     w_k = 2^{n-1} n! sqrt(pi) / (n H_{n-1}(x_k))^2          *)
gaussHermiteNodesWeights[n_Integer] := gaussHermiteNodesWeights[n] =
  Module[{itTol = 1*^-14, nodes = ConstantArray[0., n], weights = ConstantArray[0., n],
          x, xNew, Hm1, H, jj},
    Do[
      x = N[Cos[Pi (k - 0.5) / n] Sqrt[2 n + 1]];  (* Chebyshev-based initial guess *)
      Do[
        (* Evaluate H_n(x) and H_{n-1}(x) via recurrence *)
        Hm1 = 1.; H = 2. x;
        Do[{Hm1, H} = {H, 2 x H - 2(jj-1) Hm1}, {jj, 2, n}];
        (* H_n'(x) = 2n H_{n-1}(x);  Newton step *)
        xNew = x - H / (2 n Hm1);
        If[Abs[xNew - x] < itTol, x = xNew; Break[]];
        x = xNew,
        {50}
      ];
      nodes[[k]]   = x;
      weights[[k]] = 2^(n-1) Factorial[n] Sqrt[Pi] / (n Hm1)^2,
      {k, n}
    ];
    {Sort[nodes], weights[[Ordering[nodes]]]}
  ]

(* Weighted -> True (default): integrates f(x)*Exp[-x^2] on (-Inf,Inf).
   Weighted -> False          : integrates f(x) by dividing out Exp[-x^2].
   Note: Weighted -> False can be numerically unstable for large |x|. *)
quadGaussHermite[f_, n_, weighted_:True] :=
  Module[{pts, wts},
    {pts, wts} = gaussHermiteNodesWeights[n];
    If[TrueQ[weighted],
      Total[wts * (f /@ pts)],
      Total[wts * Exp[pts^2] * (f /@ pts)]   (* divide out Exp[-x^2] weight *)
    ]
  ]

(* ── Gauss-Laguerre nodes and weights (integrates f(x)*Exp[-x] on [0,Inf)) *)
(* Laguerre recurrence: L_0=1, L_1=1-x,
   (k+1) L_{k+1} = (2k+1-x) L_k - k L_{k-1}
   Derivative:    x L_n'(x) = n(L_n(x) - L_{n-1}(x))
   Weight:        w_k = x_k / (n L_{n-1}(x_k))^2
   (Derived by substituting L_n(x_k)=0 into the (n+1)L_{n+1} recurrence.) *)
gaussLaguerreNodesWeights[n_Integer] := gaussLaguerreNodesWeights[n] =
  Module[{itTol = 1*^-14, nodes = ConstantArray[0., n], weights = ConstantArray[0., n],
          x, xNew, Lm1, L, jj, Lp},
    Do[
      (* Gatteschi initial guess for k-th Laguerre root *)
      x = N[(4k - 1) Pi^2 / (8 n)];
      Do[
        (* Evaluate L_n(x) and L_{n-1}(x) via recurrence *)
        Lm1 = 1.; L = 1 - x;
        Do[{Lm1, L} = {L, ((2jj-1-x) L - (jj-1) Lm1) / jj}, {jj, 2, n}];
        (* x L_n'(x) = n(L_n - L_{n-1}) => L_n'(x) = n(L-Lm1)/x;  Newton step *)
        Lp   = n (L - Lm1) / (x + 1*^-100);  (* guard against x≈0 *)
        If[Abs[L] < itTol, Break[]];
        xNew = Max[x - L / Lp, 1*^-10];  (* keep x strictly positive *)
        If[Abs[xNew - x] < itTol, x = xNew; Break[]];
        x = xNew,
        {50}
      ];
      nodes[[k]]   = x;
      weights[[k]] = x / (n Lm1)^2,   (* w_k = x_k/(n L_{n-1}(x_k))^2 *)
      {k, n}
    ];
    {Sort[nodes], weights[[Ordering[nodes]]]}
  ]

(* Weighted -> True (default): integrates f(x)*Exp[-x] on [0,Inf).
   Weighted -> False          : integrates f(x) by multiplying by Exp[x].
   Note: Weighted -> False can overflow for large x nodes. *)
quadGaussLaguerre[f_, n_, weighted_:True] :=
  Module[{pts, wts},
    {pts, wts} = gaussLaguerreNodesWeights[n];
    If[TrueQ[weighted],
      Total[wts * (f /@ pts)],
      Total[wts * Exp[pts] * (f /@ pts)]   (* divide out Exp[-x] weight *)
    ]
  ]

(* ── Gauss-Lobatto nodes and weights (includes both endpoints) ──────────── *)
(* Interior Lobatto nodes are roots of P'_{n-1}(x) in (-1,1).
   Legendre recurrence: P_0=1, P_1=x,  k P_k = (2k-1)x P_{k-1} - (k-1)P_{k-2}
   P'_{n-1}(x) = (n-1)(P_{n-2}(x) - x P_{n-1}(x))/(1-x^2)
   P''_{n-1}(x) = (2x P'_{n-1}(x) - (n-1)n P_{n-1}(x))/(1-x^2)   [from Legendre ODE]
   Newton step:  x <- x - P'_{n-1}(x)/P''_{n-1}(x)
   Weight:       w_k = 2/(n(n-1) P_{n-1}(x_k)^2)  for ALL nodes incl. endpoints.
   At endpoints: P_{n-1}(±1) = (±1)^{n-1} => w = 2/(n(n-1)).              *)
gaussLobattoNodesWeights[n_Integer] := gaussLobattoNodesWeights[n] =
  Module[{itTol = 1*^-14, interior, allPts, wts,
          x, xNew, Pm2, Pm1, P, Pp, Ppp, jj},
    interior = If[n > 2,
      Sort[Table[
        x = N[Cos[Pi k / (n-1)]];  (* initial guess for k-th interior node *)
        Do[
          (* Evaluate P_{n-1}(x) and P_{n-2}(x) via recurrence *)
          Pm2 = 1.; Pm1 = x;
          Do[P = ((2jj-1) x Pm1 - (jj-1) Pm2)/jj; Pm2 = Pm1; Pm1 = P, {jj, 2, n-1}];
          (* Pm1 = P_{n-1}(x),  Pm2 = P_{n-2}(x) *)
          Pp  = (n-1)(Pm2 - x Pm1) / (1 - x^2 + itTol);  (* P'_{n-1} *)
          Ppp = (2 x Pp - (n-1) n Pm1) / (1 - x^2 + itTol);  (* P''_{n-1} *)
          If[Abs[Ppp] < itTol, Break[]];
          xNew = x - Pp/Ppp;
          If[Abs[xNew - x] < itTol, x = xNew; Break[]];
          x = xNew,
          {50}
        ];
        x,
        {k, n-2}
      ]] // N,
      {}
    ];
    allPts = Join[{-1.}, interior, {1.}];
    (* Compute weights: evaluate P_{n-1} at each node via recurrence *)
    wts = Table[
      Module[{xi = allPts[[ii]]},
        Pm2 = 1.; Pm1 = xi;
        Do[P = ((2jj-1) xi Pm1 - (jj-1) Pm2)/jj; Pm2 = Pm1; Pm1 = P, {jj, 2, n-1}];
        2. / ((n (n-1)) Pm1^2)  (* w = 2/(n(n-1) P_{n-1}(x)^2) *)
      ],
      {ii, Length[allPts]}
    ];
    {allPts, wts}
  ]

quadGaussLobatto[f_, a_, b_, n_] :=
  Module[{pts, wts},
    {pts, wts} = gaussLobattoNodesWeights[n];
    (b-a)/2 Total[wts * (f /@ ((b-a)/2 pts + (a+b)/2))]
  ]

(* ── Gauss-Radau nodes and weights (includes left endpoint) ─────────────── *)
(* Interior Radau nodes are roots of g(x) = P_{n-1}(x) + P_n(x) in (-1,1).
   Three-register Legendre recurrence simultaneously provides P_{n-2}, P_{n-1}, P_n:
     Pm3 = 1 (P_0),  Pm2 = x (P_1), then update to degree n.
   g(x)  = P_{n-1} + P_n  (= Pm2 + Pm1 after loop)
   g'(x) = [(n-1)(P_{n-2}-x P_{n-1}) + n(P_{n-1}-x P_n)] / (1-x^2)
          = [(n-1)(Pm3-x Pm2) + n(Pm2-x Pm1)] / (1-x^2)
   Newton: x <- x - g/g'
   Weight at x=-1:      w = 2/n^2
   Weight at x_k>-1:    w_k = (1-x_k)/(n P_{n-1}(x_k))^2 = (1-x_k)/(n Pm2_k)^2 *)
gaussRadauNodesWeights[n_Integer] := gaussRadauNodesWeights[n] =
  Module[{itTol = 1*^-14, interior, allPts, wts,
          x, xNew, Pm3, Pm2, Pm1, P, g, gp, jj},
    interior = Sort[Table[
      x = N[Cos[Pi (2k-1) / (2n-1)]];
      Do[
        (* Three-level Legendre recurrence.
           Initialise: Pm2=P_0=1, Pm1=P_1=x, Pm3 is a carry slot.
           After loop to j=n: Pm3=P_{n-2}, Pm2=P_{n-1}, Pm1=P_n. *)
        Pm3 = 1.; Pm2 = 1.; Pm1 = x;
        Do[
          P = ((2jj-1) x Pm1 - (jj-1) Pm2)/jj;
          Pm3 = Pm2; Pm2 = Pm1; Pm1 = P,
          {jj, 2, n}
        ];
        g  = Pm2 + Pm1;                              (* P_{n-1}+P_n *)
        gp = ((n-1)(Pm3 - x Pm2) + n(Pm2 - x Pm1)) / (1 - x^2 + itTol);
        If[Abs[gp] < itTol, Break[]];
        xNew = x - g/gp;
        If[Abs[xNew - x] < itTol, x = xNew; Break[]];
        x = xNew,
        {50}
      ];
      x,
      {k, n-1}
    ]] // N;
    allPts = Join[{-1.}, interior];
    (* Weights: re-evaluate P_{n-1} at each interior node with a clean recurrence *)
    wts = Join[{2./n^2},
      Table[
        Module[{xi = interior[[k]], Qm2, Qm1, Q},
          Qm2 = 1.; Qm1 = xi;
          Do[Q = ((2jj-1) xi Qm1 - (jj-1) Qm2)/jj; Qm2 = Qm1; Qm1 = Q, {jj, 2, n-1}];
          (* Qm1 = P_{n-1}(xi) *)
          (1 - xi) / (n Qm1)^2
        ],
        {k, n-1}
      ]
    ];
    {allPts, wts}
  ]

quadGaussRadau[f_, a_, b_, n_] :=
  Module[{pts, wts},
    {pts, wts} = gaussRadauNodesWeights[n];
    (b-a)/2 Total[wts * (f /@ ((b-a)/2 pts + (a+b)/2))]
  ]
