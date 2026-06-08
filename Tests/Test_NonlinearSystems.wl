(* Test_NonlinearSystems.wl -- Correctness tests for SolveNonlinearEquationSystem *)

NKBeginSuite["Nonlinear Equation Systems Correctness", "A", "NonlinearSystems"]

With[{F = Function[v, {v[[1]]^2 + v[[2]]^2 - 5, v[[1]] - v[[2]] - 1}]},
  NKTestNear["SolveNonlinearEquationSystem/Newton: circle-line at (2,1)",
    SolveNonlinearEquationSystem[F, {1., 2.}, Method -> "Newton"]["Solution"],
    {2., 1.}, "AbsoluteTolerance" -> 1*^-8];
  NKTestNear["SolveNonlinearEquationSystem/Broyden: circle-line at (2,1)",
    SolveNonlinearEquationSystem[F, {1., 2.}, Method -> "Broyden"]["Solution"],
    {2., 1.}, "AbsoluteTolerance" -> 1*^-6];
  NKTestNear["SolveNonlinearEquationSystem/Continuation: circle-line",
    SolveNonlinearEquationSystem[F, {1., 2.}, Method -> "Continuation"]["Solution"],
    {2., 1.}, "AbsoluteTolerance" -> 1*^-6];
  NKTest["SolveNonlinearEquationSystem: Residual < 1e-10",
    SolveNonlinearEquationSystem[F, {1., 2.}]["Residual"] < 1*^-10];
  NKTest["SolveNonlinearEquationSystem: Converged -> True",
    TrueQ[SolveNonlinearEquationSystem[F, {1., 2.}]["Converged"]]]
]

With[{G = Function[v, {v[[1]]^2+v[[2]]^2+v[[3]]^2-3,
                       v[[1]]+v[[2]]-2v[[3]],
                       v[[1]] v[[2]]-v[[3]]^2}]},
  NKTestNear["SolveNonlinearEquationSystem/Newton: 3D system -> (1,1,1)",
    SolveNonlinearEquationSystem[G, {0.9, 0.9, 0.9}]["Solution"],
    {1., 1., 1.}, "AbsoluteTolerance" -> 1*^-6]
]

With[{r = SolveNonlinearEquationSystem[
            Function[v, {v[[1]]^2+v[[2]]^2-5, v[[1]]-v[[2]]-1}],
            {1., 2.}, Method -> "Newton", ConvergenceHistory -> True]},
  NKTest["SolveNonlinearEquationSystem/ConvergenceHistory: key present",
    KeyExistsQ[r, "ResidualHistory"]];
  NKTest["SolveNonlinearEquationSystem/ConvergenceHistory: residuals decrease",
    With[{rh = r["ResidualHistory"]},
      Length[rh] > 1 && Last[rh] < First[rh]]]
]

NKEndSuite[]
