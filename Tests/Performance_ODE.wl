(* Performance_ODE.wl -- Performance regression tests for ODE solvers *)

NKBeginSuite["ODE Performance", "E", "ODE_Timing"]

NKTestPerformanceRatio["Performance/ODE RK4: 10000-step scalar y'=-y",
  SolveInitialValueProblem[Function[{t,y},-y],{t,0.,100.},1.,
    Method->"RungeKutta4",StepSize->0.01],
  "ODE_RK4_10000steps","MaxRatio"->2.0]

NKTestPerformanceRatio["Performance/ODE RK4: 1000-step 2D oscillator",
  SolveInitialValueProblem[Function[{t,y},{y[[2]],-y[[1]]}],
    {t,0.,10.},{0.,1.},Method->"RungeKutta4",StepSize->0.01],
  "ODE_RK4_System_1000steps","MaxRatio"->2.0]

NKTest["Performance/ODE RKF45: adaptive integration < 5s",
  First[AbsoluteTiming[
    SolveInitialValueProblem[Function[{t,y},{y[[2]],-y[[1]]}],
      {t,0.,100.},{0.,1.},Method->"RungeKuttaFehlberg",StepSize->0.1]]] < 5.]

With[{sol100  = SolveInitialValueProblem[Function[{t,y},-y],{t,0.,1.},1.,StepSize->0.01],
      sol1000 = SolveInitialValueProblem[Function[{t,y},-y],{t,0.,1.},1.,StepSize->0.001]},
  NKTestNear["Performance/ODE: 100-step vs 1000-step results consistent",
    Abs[Last[sol100["Solution"]]-Last[sol1000["Solution"]]],
    0.,"AbsoluteTolerance"->1*^-3]
]

NKEndSuite[]
