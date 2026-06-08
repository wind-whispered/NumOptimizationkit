(* Test_Approximation.wl -- Correctness tests for polynomial approximation *)

NKBeginSuite["Polynomial Approximation Correctness", "A", "Approximation"]

Table[
  With[{fApprox = ChebyshevApproximate[Exp, {x,0.,2.}, n],
        errMax = Max[Abs[Exp[#] - ChebyshevApproximate[Exp,{x,0.,2.},n][#] & /@ Range[0.,2.,0.1]]]},
    NKTest["ChebyshevApproximate: degree " <> ToString[n] <> " error < 10^{-" <> ToString[n-2] <> "}",
      errMax < 10^-(n-2)]],
  {n, {4, 8, 12}}]

With[{xs = Range[-2.,2.,0.5], ys = 2+3#-#^2 & /@ Range[-2.,2.,0.5],
      basis = {Function[x,1.],Function[x,x],Function[x,x^2]}},
  With[{c = LeastSquaresApproximate[xs,ys,basis]},
    NKTestNear["LeastSquaresApproximate: c[1]=2", c[[1]], 2., "AbsoluteTolerance" -> 1*^-8];
    NKTestNear["LeastSquaresApproximate: c[2]=3", c[[2]], 3., "AbsoluteTolerance" -> 1*^-8];
    NKTestNear["LeastSquaresApproximate: c[3]=-1", c[[3]],-1., "AbsoluteTolerance" -> 1*^-8]]
]

With[{R = PadeApproximate[Exp, 0., {2,2}]},
  NKTestNear["PadeApproximate [2/2]: Exp(0.5)",
    R[0.5], Exp[0.5], "RelativeTolerance" -> 0.001];
  NKTestNear["PadeApproximate [2/2]: exact at expansion point",
    R[0.], 1., "AbsoluteTolerance" -> 1*^-10]
]

With[{res = RemezApproximate[Sin, {0.,Pi}, 3]},
  NKTest["RemezApproximate: returns Coefficients/Error/Nodes",
    KeyExistsQ[res,"Coefficients"] && KeyExistsQ[res,"Error"] && KeyExistsQ[res,"Nodes"]];
  NKTest["RemezApproximate: Error > 0", res["Error"] > 0];
  NKTest["RemezApproximate: 5 equioscillation nodes for degree 3",
    Length[res["Nodes"]] == 5]
]

SeedRandom[42];
With[{pts = RandomReal[{0.,1.},{15,2}], vals = Sin[Pi #[[1]]] Cos[Pi #[[2]]] & /@ RandomReal[{0.,1.},{15,2}]},
  NKTestNear["ScatteredInterpolate/TPS: exact at data point",
    ScatteredInterpolate[pts,vals,pts[[1]]], vals[[1]], "AbsoluteTolerance" -> 0.01];
  NKTestNear["ScatteredInterpolate/InverseDistance: exact at data point",
    ScatteredInterpolate[pts,vals,pts[[1]],Method->"InverseDistance"],
    vals[[1]], "AbsoluteTolerance" -> 1*^-6]
]

NKEndSuite[]
