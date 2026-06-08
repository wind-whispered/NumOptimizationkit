(* Performance_Quadrature.wl -- Performance regression for quadrature (tests memoization) *)

NKBeginSuite["Quadrature Performance", "E", "Quad_Timing"]

NKTestPerformanceRatio["Performance/GaussLegendre: 1000 calls Sin on [0,Pi] n=5",
  Table[NumericalQuadrature[Sin,{x,0.,Pi},Points->5],{1000}],
  "Gauss_Quad_1000calls","MaxRatio"->3.0]

With[{t1=First[AbsoluteTiming[gaussLegendreNodesWeights[10]]],
      t2=First[AbsoluteTiming[gaussLegendreNodesWeights[10]]]},
  NKTest["Performance/GaussLegendre memoization: second call < 0.1x first",
    t2 < Max[t1*0.1, 1*^-4]]
]

With[{xs=Range[0.,10.,0.01],ys=Sin/@Range[0.,10.,0.01]},
  NKTestPerformanceRatio["Performance/SplineInterpolate: 1000 data points",
    Table[SplineInterpolate[xs,ys,x],{x,0.5,9.5,0.1}],
    "Spline_1000pts","MaxRatio"->3.0]
]

NKEndSuite[]
