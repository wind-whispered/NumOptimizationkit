(* PolynomialApproximation.wl -- Function approximation
   ChebyshevApproximate    : near-best uniform approximation via Chebyshev expansion
   LeastSquaresApproximate : discrete least-squares polynomial fitting
   PadeApproximate         : Pade rational approximant
   RemezApproximate        : Remez exchange algorithm (minimax)

   All linear systems are solved through our own implementations:
     Normal equations  -> laGaussElim  (partial-pivoting Gaussian elimination)
     QR back-sub       -> laGaussElim  (upper-triangular system, no pivot needed)
     Pade system       -> laGaussElim
     Remez system      -> laGaussElim
   Depends on: qrHouseholder (LinearEquationSystems.wl)
               laGaussElim   (LinearEquationSystems.wl)
               min1DGoldenSection (UnconstrainedMinimization.wl)
*)

(* ── Chebyshev Polynomial Expansion ────────────────────────────────────── *)
(* Uses ChebyshevT for evaluation only (built-in polynomial definition),
   not as a solver.  The DCT-like coefficient formula is implemented directly. *)
ChebyshevApproximate[f_, {x_, a_?NumericQ, b_?NumericQ}, n_Integer] :=
  Module[{nodes, vals, c},
    nodes = Cos[(2 Range[n+1] - 1) Pi / (2(n+1))];
    vals  = f /@ ((b-a)/2 nodes + (a+b)/2);
    c = Table[
      2/(n+1) Sum[vals[[k+1]] ChebyshevT[j, nodes[[k+1]]], {k, 0, n}],
      {j, 0, n}];
    c[[1]] /= 2;
    Function[Evaluate[x],
      Evaluate[Total[MapIndexed[#1 ChebyshevT[#2[[1]]-1, (2(x)-a-b)/(b-a)]&, c]]]
    ]
  ]

(* ── Least-Squares Fitting ──────────────────────────────────────────────── *)
(* Normal-equations variant: solve (A^T A) c = A^T y using our Gaussian elimination. *)
LeastSquaresApproximate[xs_?VectorQ, ys_?VectorQ, basis_List] :=
  Module[{A = Table[basis[[j]][xs[[i]]], {i, Length[xs]}, {j, Length[basis]}],
          ATA, ATy},
    ATA = Transpose[A] . A;
    ATy = Transpose[A] . N@ys;
    laGaussElim[N@ATA, N@ATy, True]   (* our partial-pivoting Gaussian elimination *)
  ]

(* QR-based variant: solve R c = Q^T y by back-substitution via our Gaussian elimination.
   R is upper-triangular, so laGaussElim with no pivoting is exact for back-sub. *)
LeastSquaresApproximate[xs_?VectorQ, ys_?VectorQ, basis_List, "Orthogonal"] :=
  Module[{A = Table[basis[[j]][xs[[i]]], {i, Length[xs]}, {j, Length[basis]}],
          Q, R, n = Length[basis]},
    {Q, R} = qrHouseholder[N@A];
    (* qrHouseholder returns the FULL QR (Q: m x m, R: m x n); take the thin
       factors (Q1: m x n, R1: n x n upper-triangular) so laGaussElim sees a
       square system, per the QR least-squares normal form R1 c = Q1^T y *)
    laGaussElim[R[[1;;n, 1;;n]], Transpose[Q[[All, 1;;n]]] . N@ys, False]
  ]

(* ── Pade Rational Approximant [p/q] ────────────────────────────────────── *)
PadeApproximate[f_, x0_?NumericQ, {p_Integer, q_Integer}] :=
  Module[{h = 1*^-5, taylor, A, rhs, bCoeffs, aCoeffs, xp},
    taylor = Table[
      Sum[(-1)^j Binomial[k, j] f[x0 + (k - 2j) h] / (2^k h^k), {j, 0, k}] / k! // N,
      {k, 0, p+q}
    ];
    If[q == 0,
      Function[Evaluate[xp], Evaluate[Sum[taylor[[k+1]] (xp - x0)^k, {k, 0, p}]]],
      A   = Table[If[k-j >= 0, taylor[[k-j+1]], 0.], {k, p+1, p+q}, {j, 1, q}];
      rhs = -Table[taylor[[k+1]], {k, p+1, p+q}];
      (* Solve Pade denominator system with our Gaussian elimination *)
      bCoeffs = Prepend[laGaussElim[N@A, N@rhs, True], 1.];
      aCoeffs = Table[
        Sum[If[j+1 <= Length[bCoeffs],
          bCoeffs[[j+1]] (If[k-j >= 0, taylor[[k-j+1]], 0.]), 0.],
          {j, 0, k}],
        {k, 0, p}];
      Function[Evaluate[xp],
        Evaluate[Sum[aCoeffs[[k+1]] (xp-x0)^k, {k, 0, p}] /
                 Sum[bCoeffs[[k+1]] (xp-x0)^k, {k, 0, q}]]]
    ]
  ]

(* ── Remez Exchange Algorithm (minimax polynomial) ──────────────────────── *)
RemezApproximate[f_, {a_?NumericQ, b_?NumericQ}, n_Integer,
                 maxIter_Integer:30, tol_:1*^-10] :=
  Module[{xs, A, rhs, sol, c, E, Eprev = Infinity, xExtr, k = 0, sgn},
    xs = Table[(a+b)/2 + (b-a)/2 Cos[(n+1-j) Pi/(n+1)], {j, 0, n+1}] // N;
    Catch[
      While[k++ < maxIter,
        A   = Table[Append[Table[If[j == 0, 1., xs[[i]]^j], {j, 0, n}], (-1.)^i], {i, n+2}];
        rhs = f /@ xs // N;
        (* Solve the (n+2)x(n+2) Remez system with our Gaussian elimination *)
        sol = laGaussElim[N@A, N@rhs, True];
        E   = Last[sol]; c = Most[sol];
        (* Convergence = the leveled (equioscillation) error has stabilized
           between sweeps -- NOT that it is near zero. For non-polynomial f
           the minimax error itself converges to a small but generally
           nonzero constant, so testing Abs[E] < tol would never trigger. *)
        If[Abs[E - Eprev] < tol, Throw[<|"Coefficients" -> c, "Error" -> Abs[E], "Nodes" -> xs|>]];
        Eprev = E;
        sgn = Sign[f[xs[[1]]] - c . Table[If[j == 0, 1., xs[[1]]^j], {j, 0, n}]];
        (* The next reference set must again hold n+2 alternation points.
           Keep both endpoints (which generically belong to the
           equioscillation set on a closed interval, matching the initial
           Chebyshev-node guess above) and locate the n interior extrema
           by searching the n gaps between them -- this keeps xs at a
           constant length n+2 across sweeps instead of shrinking to n+1. *)
        xExtr = Join[{a},
          Table[
            Module[{r = min1DGoldenSection[
              Function[xx, sgn (-1.)^i (f[xx] - c . Table[If[j == 0, 1., xx^j], {j, 0, n}])],
              xs[[i]], xs[[i+1]], tol, 200]},
              r["Point"]],
            {i, 2, n+1}],
          {b}];
        xs = xExtr
      ]
    ];
    <|"Coefficients" -> c, "Error" -> Abs[E], "Nodes" -> xs|>
  ]
