(* Convergence_Quadrature.wl -- Convergence order verification for quadrature *)

NKBeginSuite["Quadrature Convergence Order Verification", "B", "Quadrature_Orders"]

With[{f=Sin, a=0., b=Pi, exact=2., mVals={10,20,40,80}},
  NKTestConvergenceOrder["Quadrature/Trapezoidal: order 2",
    1.0/mVals,
    Table[Abs[NumericalQuadrature[f,{x,a,b},Method->"Trapezoidal",Intervals->m]-exact],{m,mVals}],
    2, "OrderTolerance"->0.2];
  NKTestConvergenceOrder["Quadrature/Simpson: order 4",
    1.0/mVals,
    Table[Abs[NumericalQuadrature[f,{x,a,b},Method->"Simpson",Intervals->m]-exact],{m,mVals}],
    4, "OrderTolerance"->0.3]
]

NKTest["Quadrature/Romberg: more accurate than Simpson at 50 intervals",
  Abs[NumericalQuadrature[Sin,{x,0.,Pi},Method->"Romberg"]-2.] <
  Abs[NumericalQuadrature[Sin,{x,0.,Pi},Method->"Simpson",Intervals->50]-2.]]

With[{err5 = Abs[NumericalQuadrature[Sin,{x,0.,Pi},Points->5]-2.],
      err15= Abs[NumericalQuadrature[Sin,{x,0.,Pi},Points->15]-2.]},
  NKTest["Quadrature/GaussLegendre: n=5 more accurate than n=3", err5 < Abs[NumericalQuadrature[Sin,{x,0.,Pi},Points->3]-2.]];
  NKTest["Quadrature/GaussLegendre: n=15 reaches machine precision", err15 < 1*^-12]
]

NKEndSuite[]
