(* Test_Optimization.wl -- Correctness tests for FindMinimum1D and FindMinimumND *)

NKBeginSuite["Unconstrained Minimization Correctness", "A", "Optimization"]

With[{f = #^2 - Cos[#] &, exact = 0.45018361129763},
  Scan[Function[m,
    NKTestNear["FindMinimum1D/" <> m <> ": x^2-Cos[x] on [-2,2]",
      FindMinimum1D[f, {-2., 2.}, Method -> m]["Point"],
      exact, "RelativeTolerance" -> 1*^-4]],
  {"GoldenSection", "Fibonacci", "QuadraticInterpolation", "Newton"}]
]

Scan[Function[m,
  NKTestNear["FindMinimum1D/" <> m <> ": x^2 on [-1,1] -> 0",
    FindMinimum1D[#^2 &, {-1., 1.}, Method -> m]["Point"],
    0., "AbsoluteTolerance" -> 1*^-6]],
{"GoldenSection", "Fibonacci", "Newton"}]

NKTest["FindMinimum1D: result has Point/Value/Iterations keys",
  With[{r = FindMinimum1D[#^2 - Cos[#] &, {-2., 2.}]},
    KeyExistsQ[r,"Point"] && KeyExistsQ[r,"Value"] && KeyExistsQ[r,"Iterations"]]]

With[{rb = Function[x, (1-x[[1]])^2 + 100(x[[2]]-x[[1]]^2)^2]},
  Scan[Function[m,
    NKTestNear["FindMinimumND/" <> m <> ": Rosenbrock -> (1,1)",
      FindMinimumND[rb, {-1., 0.5}, Method -> m, MaxIterations -> 5000]["Point"],
      {1., 1.}, "AbsoluteTolerance" -> 0.01]],
  {"BFGS", "DFP", "ConjugateGradientPR"}]
]

With[{qf = Function[x, (x[[1]]-2)^2 + (x[[2]]+1)^2]},
  NKTestNear["FindMinimumND/GradientDescent: quadratic min at (2,-1)",
    FindMinimumND[qf, {0., 0.}, Method -> "GradientDescent"]["Point"],
    {2., -1.}, "AbsoluteTolerance" -> 1*^-4];
  With[{r = FindMinimumND[qf, {0., 0.}, Method -> "Newton"]},
    NKTestNear["FindMinimumND/Newton: quadratic exact min",
      r["Point"], {2., -1.}, "AbsoluteTolerance" -> 1*^-6];
    NKTest["FindMinimumND/Newton: Converged -> True", TrueQ[r["Converged"]]]
  ]
]

NKTestNear["FindMinimumND/NelderMead: Rosenbrock 2D",
  FindMinimumND[Function[x,(1-x[[1]])^2+100(x[[2]]-x[[1]]^2)^2],
    {-1., 0.5}, Method -> "NelderMead", MaxIterations -> 10000]["Point"],
  {1., 1.}, "AbsoluteTolerance" -> 0.05]

With[{r = FindMinimumND[Function[x,(x[[1]]-2)^2+(x[[2]]+1)^2],
          {0., 0.}, Method -> "BFGS", ConvergenceHistory -> True]},
  NKTest["FindMinimumND/ConvergenceHistory: ValueHistory key present",
    KeyExistsQ[r, "ValueHistory"]];
  NKTest["FindMinimumND/ConvergenceHistory: values non-increasing",
    With[{vh = r["ValueHistory"]}, Length[vh] > 1 && Last[vh] <= First[vh]]]
]

With[{r = FindConstrainedMinimum[Function[x, x[[1]]^2+x[[2]]^2], {1.,1.},
           LowerBounds -> {1., 1.}, UpperBounds -> {5., 5.}]},
  NKTestNear["FindConstrainedMinimum: box constraint forces (1,1)",
    r["Point"], {1., 1.}, "AbsoluteTolerance" -> 0.01]
]

NKEndSuite[]
