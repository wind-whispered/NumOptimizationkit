(* Test_LinearAlgebra.wl -- Correctness tests for SolveLinearEquationSystem and decompositions *)

NKBeginSuite["Linear Equation Systems Correctness", "A", "LinearAlgebra"]

With[{A = {{4.,2.,1.},{2.,5.,3.},{1.,3.,6.}}, b = {1.,2.,3.}},
  Scan[Function[m,
    NKTestNear["SolveLinearEquationSystem/" <> m <> ": residual",
      Norm[A . SolveLinearEquationSystem[A, b, Method -> m] - b],
      0., "AbsoluteTolerance" -> 1*^-10]],
  {"GaussianEliminationPivot","GaussianElimination","GaussJordan",
   "LU","Cholesky","LDLT"}];
  Scan[Function[m,
    NKTestNear["SolveLinearEquationSystem/" <> m <> ": residual (iterative)",
      Norm[A . SolveLinearEquationSystem[A, b, Method -> m, MaxIterations->2000] - b],
      0., "AbsoluteTolerance" -> 1*^-6]],
  {"Jacobi","GaussSeidel","SOR","ConjugateGradient"}];
  With[{T = {{4.,-1.,0.},{-1.,4.,-1.},{0.,-1.,4.}}, bt = {1.,1.,1.}},
    NKTestNear["SolveLinearEquationSystem/Tridiagonal: residual",
      Norm[T . SolveLinearEquationSystem[T, bt, Method->"Tridiagonal"] - bt],
      0., "AbsoluteTolerance" -> 1*^-12]]
]

With[{A = {{4.,2.,1.},{2.,5.,3.},{1.,3.,6.}}},
  With[{f = LUDecompose[A]},
    NKTestNear["LUDecompose: ||P.A - L.U||",
      Norm[f[[3]].A - f[[1]].f[[2]]], 0., "AbsoluteTolerance" -> 1*^-12]];
  With[{QR = QRDecompose[A]},
    NKTestNear["QRDecompose: ||A - Q.R||", Norm[A - QR[[1]].QR[[2]]],
      0., "AbsoluteTolerance" -> 1*^-12];
    NKTestNear["QRDecompose: ||Q^T.Q - I||",
      Norm[Transpose[QR[[1]]].QR[[1]] - IdentityMatrix[3]],
      0., "AbsoluteTolerance" -> 1*^-12]];
  With[{L = CholeskyDecompose[A]},
    NKTestNear["CholeskyDecompose: ||A - L.L^T||",
      Norm[A - L.Transpose[L]], 0., "AbsoluteTolerance" -> 1*^-12]];
  With[{ld = LDLTDecompose[A]},
    NKTestNear["LDLTDecompose: ||A - L.D.L^T||",
      Norm[A - ld[[1]].DiagonalMatrix[ld[[2]]].Transpose[ld[[1]]]],
      0., "AbsoluteTolerance" -> 1*^-12]]
]

NKTestFails["SolveLinearEquationSystem: non-square A returns $Failed",
  SolveLinearEquationSystem[{{1.,2.},{3.,4.},{5.,6.}},{1.,2.,3.}]]
NKTestFails["SolveLinearEquationSystem: dim mismatch returns $Failed",
  SolveLinearEquationSystem[{{1.,2.},{3.,4.}},{1.,2.,3.}]]

NKEndSuite[]
