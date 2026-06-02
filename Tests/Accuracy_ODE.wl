(* Accuracy_ODE.wl -- ODE accuracy benchmarks vs reference solutions *)

NKBeginSuite["ODE Accuracy Benchmarks", "D", "ODE_Reference"]

With[{mu=0.1, f=Function[{t,y},{y[[2]],mu(1-y[[1]]^2)y[[2]]-y[[1]]}], y0={2.,0.}},
  With[{refSol=Quiet@NDSolve[{y1'[t]==y2[t],y2'[t]==mu(1-y1[t]^2)y2[t]-y1[t],
                               y1[0]==2,y2[0]==0},{y1,y2},{t,0,5},
                               WorkingPrecision->30,MaxStepSize->0.001]},
    If[Head[refSol]===List,
      With[{ref={y1[5],y2[5]}/.First[refSol]//N,
            ours=Last[SolveInitialValueProblem[f,{t,0.,5.},y0,
                        Method->"RungeKutta4",StepSize->0.01]["Solution"]]},
        NKTestNear["ODE/RK4 vs NDSolve: van der Pol mu=0.1 at t=5",
          Norm[ours-ref],0.,"AbsoluteTolerance"->0.01]],
      NKTest["ODE/van der Pol: NDSolve reference available",False]]]
]

With[{rob=Function[{t,y},{-0.04y[[1]]+1e4 y[[2]]y[[3]],
                            0.04y[[1]]-1e4 y[[2]]y[[3]]-3e7 y[[2]]^2,
                            3e7 y[[2]]^2}], y0={1.,0.,0.}},
  With[{sol=SolveInitialValueProblem[rob,{t,0.,0.01},y0,Method->"BDF2",StepSize->1*^-4]},
    NKTestNear["ODE/BDF2 Robertson: conservation y1+y2+y3=1 at t=0.01",
      Plus@@Last[sol["Solution"]],1.,"AbsoluteTolerance"->1*^-6]]
]

With[{sol=SolveInitialValueProblem[Function[{t,y},-I y],{t,0.,2Pi},1.+0.I,
           Method->"RungeKutta4",StepSize->0.02]},
  NKTestNear["ODE/Complex RK4: |y(2pi)| = 1 (periodicity)",
    Abs[Last[sol["Solution"]]],1.,"AbsoluteTolerance"->0.001]
]

NKEndSuite[]
