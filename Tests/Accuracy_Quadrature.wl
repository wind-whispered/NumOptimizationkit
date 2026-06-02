(* Accuracy_Quadrature.wl -- Quadrature accuracy benchmarks vs reference solutions *)

NKBeginSuite["Quadrature Accuracy Benchmarks", "D", "Quad_Reference"]

NKTestNear["NumericalQuadrature/GaussLegendre: Gaussian integral ~ sqrt(pi)",
  NumericalQuadrature[Exp[-#^2]&,{x,-5.,5.},Method->"GaussLegendre",Points->10],
  N[Sqrt[Pi]],"AbsoluteTolerance"->1*^-8]

With[{poly=Function[x,x^5-x^3/2+x/6]},
  NKTestNear["NumericalQuadrature/GaussLegendre: degree-5 poly exact with 3 nodes",
    NumericalQuadrature[poly,{x,0.,1.},Method->"GaussLegendre",Points->3],
    N[1/6-1/8+1/12],"AbsoluteTolerance"->1*^-14]
]

NKTestNear["NumericalQuadrature/GaussHermite: x^2*Exp[-x^2] = sqrt(pi)/2",
  NumericalQuadrature[#^2&,{x,-1.,1.},Method->"GaussHermite",Points->5,Weighted->True],
  N[Sqrt[Pi]/2],"RelativeTolerance"->1*^-8]

With[{ref=Quiet@NIntegrate[Sin[x]/(1+x^2),{x,0,3}]},
  If[NumericQ[ref],
    NKTestNear["NumericalQuadrature/GaussLegendre vs NIntegrate: sin(x)/(1+x^2)",
      NumericalQuadrature[Function[x,Sin[x]/(1+x^2)],{x,0.,3.},
        Method->"GaussLegendre",Points->10],
      ref,"RelativeTolerance"->1*^-4],
    NKTest["Accuracy_Quadrature: NIntegrate unavailable",True]]
]

NKEndSuite[]
