(* Convergence_RootFinding.wl -- Per-iteration convergence order for root finders *)

NKBeginSuite["Root Finding Convergence Order", "B", "RootFinding_Orders"]

With[{f=Function[x,x^3-x-2], bracket={1.,2.}},
  With[{rh = FindEquationRoot[f,bracket,Method->"Bisection",MaxIterations->30,ConvergenceHistory->True]["ResidualHistory"]},
    NKTest["RootFinding/Bisection: residuals decrease monotonically",
      And @@ (Drop[rh,-1] >= Rest[rh])];
    NKTest["RootFinding/Bisection: 20 iterations improve by factor ~2^20",
      rh[[1]]/rh[[-1]] > 1*^5]];
  With[{rh = FindEquationRoot[f,{1.5,0.},Method->"Newton",MaxIterations->20,ConvergenceHistory->True]["ResidualHistory"]},
    NKTest["RootFinding/Newton: converges in <= 10 iterations",
      Length[rh] <= 10 || rh[[10]] < 1*^-10]];
  With[{rh = FindEquationRoot[f,bracket,Method->"Brent",MaxIterations->50,ConvergenceHistory->True]["ResidualHistory"]},
    NKTest["RootFinding/Brent: final residual < 1e-8", Last[rh] < 1*^-8]]
]

With[{errors = Table[
        Abs[FindEquationRoot[Sin,{3,4},Method->"Newton",
                             WorkingPrecision->p, Tolerance->10^(-p+2)]["Root"]-Pi],
        {p,{15,20,30}}]},
  NKTest["RootFinding/WorkingPrecision: higher precision -> smaller error",
    errors[[1]] > errors[[2]] > errors[[3]]];
  NKTest["RootFinding/WorkingPrecision 30: root of Sin accurate to 25 digits",
    errors[[3]] < 10^-25]
]

NKEndSuite[]
