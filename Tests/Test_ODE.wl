(* Test_ODE.wl -- Correctness tests for SolveInitialValueProblem *)

NKBeginSuite["ODE Initial Value Problem Correctness", "A", "ODE"]

With[{f = Function[{t, y}, -y], exact = N[Exp[-1]]},
  Scan[Function[m,
    NKTestNear["SolveInitialValueProblem/" <> m <> ": y'=-y, y(1)=Exp[-1]",
      Last[SolveInitialValueProblem[f, {t, 0., 1.}, 1.,
            Method -> m, StepSize -> 0.01]["Solution"]],
      exact, "AbsoluteTolerance" -> 0.002]],
  {"Euler", "Heun", "RungeKutta3", "RungeKutta4", "RungeKuttaFehlberg"}];
  Scan[Function[m,
    NKTestNear["SolveInitialValueProblem/" <> m <> ": y'=-y, y(1)=Exp[-1]",
      Last[SolveInitialValueProblem[f, {t, 0., 1.}, 1.,
            Method -> m, StepSize -> 0.01]["Solution"]],
      exact, "AbsoluteTolerance" -> 0.002]],
  {"AdamsBashforth", "AdamsPC4", "HammingPC"}]
]

With[{stiff = Function[{t, y}, -100. y], exact = N[Exp[-10.]]},
  Scan[Function[m,
    NKTestNear["SolveInitialValueProblem/" <> m <> " (stiff): y'=-100y, y(0.1)",
      Last[SolveInitialValueProblem[stiff, {t, 0., 0.1}, 1.,
            Method -> m, StepSize -> 0.01]["Solution"]],
      exact, "AbsoluteTolerance" -> 0.01]],
  {"TrapezoidalRule", "BDF2", "BDF4"}]
]

With[{sol = SolveInitialValueProblem[Function[{t,y},{y[[2]],-y[[1]]}],
             {t, 0., 2 Pi}, {0., 1.}, Method -> "RungeKutta4", StepSize -> 0.05]},
  NKTestNear["SolveInitialValueProblem/RK4 system: y1(pi/2) = 1",
    sol["Solution"][[Round[Pi/2/0.05]+1, 1]], 1., "AbsoluteTolerance" -> 1*^-4];
  NKTestNear["SolveInitialValueProblem/RK4 system: energy at t=2pi",
    Plus @@ (Last[sol["Solution"]]^2), 1., "AbsoluteTolerance" -> 1*^-4]
]

With[{sol = SolveInitialValueProblem[Function[{t,y}, -I y], {t, 0., Pi}, 1.+0.I,
             Method -> "RungeKutta4", StepSize -> 0.05]},
  NKTest["SolveInitialValueProblem: complex ODE returns complex solution",
    !FreeQ[sol["Solution"], Complex]];
  NKTestNear["SolveInitialValueProblem: complex ODE |y(pi)| ~ 1",
    Abs[Last[sol["Solution"]]], 1., "AbsoluteTolerance" -> 0.01];
  NKTestNear["SolveInitialValueProblem: complex ODE Re[y(pi)] ~ -1",
    Re[Last[sol["Solution"]]], -1., "AbsoluteTolerance" -> 0.01]
]

With[{sol = SolveInitialValueProblem[Function[{t,y},-y], {t,0.,1.}, 1.]},
  NKTest["SolveInitialValueProblem: Grid and Solution same length",
    Length[sol["Grid"]] == Length[sol["Solution"]]];
  NKTest["SolveInitialValueProblem: Grid starts at t0",
    First[sol["Grid"]] == 0.]
]

NKTestFails["SolveInitialValueProblem: t0 >= tEnd returns $Failed",
  SolveInitialValueProblem[Function[{t,y},-y], {t, 1., 0.}, 1.]]

NKEndSuite[]
