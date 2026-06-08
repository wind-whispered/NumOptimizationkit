(* Convergence_ODE.wl -- Convergence order verification for ODE methods *)

NKBeginSuite["ODE Convergence Order Verification", "B", "ODE_Orders"]

With[{f = Function[{t,y},-y], exact = N[Exp[-1]], hVals = {0.1,0.05,0.025,0.0125}},
  methodErrors[m_] := Table[
    Abs[Last[SolveInitialValueProblem[f,{t,0.,1.},1.,Method->m,StepSize->h]["Solution"]]-exact],
    {h, hVals}];
  NKTestConvergenceOrder["ODE/Euler order 1",    hVals, methodErrors["Euler"],         1, "OrderTolerance"->0.2];
  NKTestConvergenceOrder["ODE/Heun order 2",     hVals, methodErrors["Heun"],          2, "OrderTolerance"->0.2];
  NKTestConvergenceOrder["ODE/RK3 order 3",      hVals, methodErrors["RungeKutta3"],   3, "OrderTolerance"->0.3];
  NKTestConvergenceOrder["ODE/RK4 order 4",      hVals, methodErrors["RungeKutta4"],   4, "OrderTolerance"->0.3];
  NKTestConvergenceOrder["ODE/Trapezoidal order 2", hVals, methodErrors["TrapezoidalRule"], 2, "OrderTolerance"->0.3];
  With[{hBDF={0.1,0.05,0.025}},
    NKTestConvergenceOrder["ODE/BDF2 order 2", hBDF,
      Table[Abs[Last[SolveInitialValueProblem[f,{t,0.,1.},1.,Method->"BDF2",StepSize->h]["Solution"]]-exact],
            {h,hBDF}],
      2, "OrderTolerance"->0.4]]
]

NKTest["ODE/RKF45: error < 1e-4 at default settings",
  Abs[Last[SolveInitialValueProblem[Function[{t,y},-y],{t,0.,1.},1.,
             Method->"RungeKuttaFehlberg",StepSize->0.1]["Solution"]] - N[Exp[-1]]] < 1*^-4]

NKEndSuite[]
