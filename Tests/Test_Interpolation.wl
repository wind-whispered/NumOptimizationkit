(* Test_Interpolation.wl -- Correctness tests for interpolation functions *)

NKBeginSuite["Polynomial Interpolation Correctness", "A", "Interpolation"]

With[{xs = {0.,1.,2.,3.,4.}, ys = Sin /@ {0.,1.,2.,3.,4.}},
  Scan[Function[xi,
    NKTestNear["PolynomialInterpolate/Newton: exact at node x=" <> ToString[xi],
      PolynomialInterpolate[xs, ys, N@xi], Sin[N@xi], "AbsoluteTolerance" -> 1*^-12]],
  {0,1,2,3,4}];
  Scan[Function[xq,
    NKTestNear["PolynomialInterpolate: Newton vs Lagrange at x=" <> ToString[xq],
      PolynomialInterpolate[xs,ys,xq,"Method"->"Newton"] -
      PolynomialInterpolate[xs,ys,xq,"Method"->"Lagrange"],
      0., "AbsoluteTolerance" -> 1*^-12]],
  {0.5, 1.5, 2.5}]
]

With[{xs = {0.,Pi/2,Pi}, ys = Sin/@{0.,Pi/2,Pi}, dys = Cos/@{0.,Pi/2,Pi}},
  Scan[Function[xi,
    NKTestNear["HermiteInterpolate: exact at node",
      HermiteInterpolate[xs,ys,dys,N@xi], Sin[N@xi], "AbsoluteTolerance" -> 1*^-10]],
  {0,Pi/2,Pi}];
  NKTestNear["HermiteInterpolate: accurate at Pi/4",
    HermiteInterpolate[xs,ys,dys,Pi/4], Sin[Pi/4], "AbsoluteTolerance" -> 0.01]
]

With[{xs = Range[0.,5.]*1., ys = Range[0.,5.]*1.},
  NKTestNear["SplineInterpolate/Cubic: exact on linear data at 2.7",
    SplineInterpolate[xs,ys,2.7], 2.7, "AbsoluteTolerance" -> 1*^-10]
]

With[{xs = {0.,1.,2.,3.}, ys = {0.,1.,0.,1.}},
  NKTestNear["SplineInterpolate/Degree 1: at midpoint",
    SplineInterpolate[xs,ys,0.5,Degree->1], 0.5, "AbsoluteTolerance" -> 1*^-10]
]

NKTestFails["PolynomialInterpolate: xs/ys length mismatch",
  PolynomialInterpolate[{1.,2.,3.},{1.,2.},1.5]]
NKTestFails["PolynomialInterpolate: duplicate nodes",
  PolynomialInterpolate[{1.,1.,3.},{1.,2.,3.},1.5]]

NKEndSuite[]
