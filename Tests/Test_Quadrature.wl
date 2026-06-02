(* Test_Quadrature.wl -- Correctness tests for NumericalQuadrature *)

NKBeginSuite["Numerical Quadrature Correctness", "A", "Quadrature"]

Scan[Function[m,
  NKTestNear["NumericalQuadrature/" <> m <> ": Sin on [0,Pi] = 2",
    NumericalQuadrature[Sin, {x, 0., Pi}, Method -> m],
    2., "RelativeTolerance" -> 1*^-5]],
{"Trapezoidal", "Simpson", "AdaptiveSimpson", "Romberg",
 "GaussLegendre", "GaussLobatto", "GaussRadau"}]

NKTestNear["NumericalQuadrature/GaussLegendre: Exp on [0,1] = e-1",
  NumericalQuadrature[Exp, {x, 0., 1.}, Method -> "GaussLegendre", Points -> 8],
  N[E - 1], "RelativeTolerance" -> 1*^-12]

With[{n = 5},
  NKTestNear["NumericalQuadrature/GaussLegendre: x^9 exact with 5 nodes",
    NumericalQuadrature[#^9 &, {x, 0., 1.}, Method -> "GaussLegendre", Points -> n],
    1./10, "RelativeTolerance" -> 1*^-12]
]

NKTestNear["NumericalQuadrature/GaussHermite (Weighted): x^2 * Exp[-x^2]",
  NumericalQuadrature[#^2 &, {x, -1., 1.}, Method -> "GaussHermite",
    Points -> 5, Weighted -> True],
  N[Sqrt[Pi]/2], "RelativeTolerance" -> 1*^-8]

NKTestNear["NumericalQuadrature/GaussLaguerre (Weighted): x * Exp[-x]",
  NumericalQuadrature[Function[x, x], {x, 0., 1.}, Method -> "GaussLaguerre",
    Points -> 5, Weighted -> True],
  1., "RelativeTolerance" -> 1*^-8]

NKTestNear["NumericalQuadrature: a=b returns 0",
  Quiet[NumericalQuadrature[Sin, {x, 1., 1.}]], 0., "AbsoluteTolerance" -> 1*^-15]

NKEndSuite[]
