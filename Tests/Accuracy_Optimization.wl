(* Accuracy_Optimization.wl -- Optimization accuracy benchmarks vs reference solutions *)

NKBeginSuite["Optimization Accuracy Benchmarks", "D", "Opt_Reference"]

With[{rb=Function[x,(1-x[[1]])^2+100(x[[2]]-x[[1]]^2)^2]},
  Scan[Function[m,
    NKTestNear["Optimization/" <> m <> " vs exact: Rosenbrock -> (1,1)",
      Norm[FindMinimumND[rb,{-1.,0.5},Method->m,MaxIterations->5000]["Point"]-{1.,1.}],
      0.,"AbsoluteTolerance"->0.01]],
  {"BFGS","DFP","ConjugateGradientPR"}]
]

With[{sphere=Function[x,Total[x^2]],n=5},
  With[{r=FindMinimumND[sphere,ConstantArray[1.,n],Method->"BFGS"]},
    NKTestNear["Optimization/BFGS: sphere 5D -> origin",
      Norm[r["Point"]],0.,"AbsoluteTolerance"->1*^-6]
  ]
]

With[{booth=Function[x,(x[[1]]+2x[[2]]-7)^2+(2x[[1]]+x[[2]]-5)^2]},
  NKTestNear["Optimization/BFGS: Booth function -> (1,3)",
    Norm[FindMinimumND[booth,{0.,0.},Method->"BFGS"]["Point"]-{1.,3.}],
    0.,"AbsoluteTolerance"->1*^-5]
]

NKEndSuite[]
