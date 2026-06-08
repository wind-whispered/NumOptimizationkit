(* Test_PDE.wl -- Correctness tests for SolveParabolicPDE, SolveHyperbolicPDE, SolveEllipticPDE *)

NKBeginSuite["PDE Numerical Solution Correctness", "A", "PDE"]

With[{ic = Function[x, Sin[Pi x]], bc = {Function[t,0.], Function[t,0.]},
      T = 0.05, exactC = N[Sin[Pi 0.5] Exp[-Pi^2 0.05]]},
  Scan[Function[m,
    With[{u = SolveParabolicPDE[1., {x,0.,1.,49}, {t,0.,T,50}, ic, bc, Method -> m]},
      NKTestNear["SolveParabolicPDE/" <> m <> ": u(0.5,0.05) vs exact",
        u["Solution"][[-1, 26]], exactC, "AbsoluteTolerance" -> 0.01]]],
  {"CrankNicolson", "ForwardDifference", "BackwardDifference"}]
]

With[{u = SolveParabolicPDE[1., {x,0.,1.,9}, {t,0.,0.01,10},
           Function[x,Sin[Pi x]], {Function[t,0.],Function[t,0.]}]},
  NKTest["SolveParabolicPDE: Solution matrix dims match grids",
    Length[u["Solution"]] == Length[u["TimeGrid"]] &&
    Length[First[u["Solution"]]] == Length[u["SpatialGrid"]]]
]

With[{uw = SolveHyperbolicPDE[1., {x,0.,1.,49}, {t,0.,1.,100},
           Function[x,Sin[Pi x]], Function[x,0.],
           {Function[t,0.],Function[t,0.]}],
      exactH = N[Sin[Pi 0.5] Cos[Pi 0.5]]},
  With[{tIdx = First[Nearest[uw["TimeGrid"] -> Range[Length[uw["TimeGrid"]]], 0.5]]},
    NKTestNear["SolveHyperbolicPDE: u(0.5,0.5) vs exact",
      uw["Solution"][[tIdx, 26]], exactH, "AbsoluteTolerance" -> 0.05]
  ]
]

With[{uE = SolveEllipticPDE[Function[{x,y},-2.],
           {x,0.,1.,9},{y,0.,1.,9},
           {Function[y,0.],Function[y,0.],Function[x,0.],Function[x,0.]}]},
  NKTestNear["SolveEllipticPDE: nabla^2 u=-2, center ~ 0.5",
    uE["Solution"][[6,6]], 0.5, "AbsoluteTolerance" -> 0.02]
]

NKEndSuite[]
