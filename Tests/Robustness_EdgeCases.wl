(* Robustness_EdgeCases.wl -- Edge case and boundary input tests *)

NKBeginSuite["Edge Cases", "C", "EdgeCases"]

NKTestFails["FindEquationRoot: positive f on bracket",
  FindEquationRoot[Function[x,x^2+1],{-2.,2.},Method->"Bisection"]]
With[{r=FindEquationRoot[Function[x,x^3-x-2],{1.,2.},MaxIterations->1]},
  NKTest["FindEquationRoot: MaxIterations=1 returns Association",AssociationQ[r]];
  NKTest["FindEquationRoot: MaxIterations=1 gives Converged->False",!TrueQ[r["Converged"]]]
]
NKTestNear["NumericalQuadrature: a=b returns 0",
  Quiet[NumericalQuadrature[Sin,{x,1.,1.}]],0.,"AbsoluteTolerance"->1*^-15]
NKTestFails["SolveLinearEquationSystem: singular matrix",
  SolveLinearEquationSystem[{{1.,2.},{2.,4.}},{1.,3.}]]
NKTestFails["SolveLinearEquationSystem: non-square",
  SolveLinearEquationSystem[{{1.,2.,3.},{4.,5.,6.}},{1.,2.}]]
NKTestFails["SolveLinearEquationSystem: dim mismatch",
  SolveLinearEquationSystem[{{1.,2.},{3.,4.}},{1.,2.,3.}]]
NKTestFails["CholeskyDecompose: non-symmetric",
  CholeskyDecompose[{{1.,2.},{3.,4.}}]]
NKTestFails["CholeskyDecompose: indefinite",
  CholeskyDecompose[{{1.,0.},{0.,-1.}}]]
NKTestFails["FindMatrixEigenvalues/Jacobi: non-symmetric",
  FindMatrixEigenvalues[{{1.,2.},{3.,4.}},Method->"JacobiMethod"]]
NKTestFails["SolveInitialValueProblem: t0>=tEnd",
  SolveInitialValueProblem[Function[{t,y},-y],{t,1.,0.},1.]]
NKTestFails["PolynomialInterpolate: length mismatch",
  PolynomialInterpolate[{1.,2.,3.},{1.,2.},1.5]]
NKTestFails["PolynomialInterpolate: duplicate xs",
  PolynomialInterpolate[{1.,1.,3.},{1.,2.,3.},1.5]]
NKTestFails["SolveBoundaryValueProblem: a>=b",
  SolveBoundaryValueProblem[Function[{x,y,dy},-y],{x,Pi,0.},{0.,0.}]]
With[{r=SolveInitialValueProblem[Function[{t,y},-y],{t,0.,0.001},1.,
           Method->"AdamsBashforth",StepSize->0.001]},
  NKTest["SolveInitialValueProblem/AdamsBashforth: short interval fallback",
    AssociationQ[r] && KeyExistsQ[r,"Grid"]]
]

NKEndSuite[]
