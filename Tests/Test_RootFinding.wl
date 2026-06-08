(* Test_RootFinding.wl -- Correctness tests for FindEquationRoot *)

NKBeginSuite["Equation Root Finding Correctness", "A", "RootFinding"]

With[{f = Function[x, x^3-x-2],
      exact = 1.5213797068},
  Scan[Function[m,
    NKTestNear["FindEquationRoot/" <> m <> ": x^3-x-2 in [1,2]",
      FindEquationRoot[f, {1., 2.}, Method -> m]["Root"],
      exact, "AbsoluteTolerance" -> 1*^-6]],
  {"Bisection", "RegulaFalsi", "Brent", "Secant", "Muller"}];
  Scan[Function[m,
    NKTestNear["FindEquationRoot/" <> m <> ": x^3-x-2 from x0=1.5",
      FindEquationRoot[f, {1.5, 0.}, Method -> m]["Root"],
      exact, "AbsoluteTolerance" -> 1*^-6]],
  {"Newton", "Halley", "Steffensen"}];
  NKTest["FindEquationRoot/Bisection: Converged -> True",
    TrueQ[FindEquationRoot[f, {1., 2.}]["Converged"]]];
  NKTest["FindEquationRoot/Newton: Residual < 1e-10",
    FindEquationRoot[f, {1.5, 0.}, Method -> "Newton"]["Residual"] < 1*^-10]
]

NKTestNear["FindEquationRoot/FixedPoint: (x+2)^(1/3) fixed point",
  FindEquationRoot[Function[x, (x+2.)^(1/3)], {1.5, 0.}, Method -> "FixedPoint"]["Root"],
  1.5213797068, "AbsoluteTolerance" -> 1*^-4]

NKTestNear["FindEquationRoot/Aitken: (x+2)^(1/3) fixed point",
  FindEquationRoot[Function[x, (x+2.)^(1/3)], {1.5, 0.}, Method -> "Aitken"]["Root"],
  1.5213797068, "AbsoluteTolerance" -> 1*^-4]

With[{r = FindEquationRoot[Sin, {3, 4}, Method -> "Newton",
                            WorkingPrecision -> 30, Tolerance -> 10^-28]},
  NKTest["FindEquationRoot/WorkingPrecision 30: root of Sin near Pi",
    Precision[r["Root"]] >= 25 && Abs[r["Root"] - Pi] < 10^-25]
]

With[{r = FindEquationRoot[Function[x, x^3-x-2], {1., 2.},
           Method -> "Newton", ConvergenceHistory -> True]},
  NKTest["FindEquationRoot/ConvergenceHistory: key present",
    KeyExistsQ[r, "ResidualHistory"]];
  NKTest["FindEquationRoot/ConvergenceHistory: residuals decrease",
    With[{rh = r["ResidualHistory"]}, Length[rh] > 1 && Last[rh] < First[rh]]]
]

NKTestFails["FindEquationRoot: no sign change returns $Failed",
  FindEquationRoot[Function[x, x^2+1], {-2., 2.}, Method -> "Bisection"]]

NKEndSuite[]
