(* Test_Eigenanalysis.wl -- Correctness tests for eigenvalue functions *)

NKBeginSuite["Matrix Eigenanalysis Correctness", "A", "Eigenanalysis"]

With[{M = {{3.,1.,0.},{1.,3.,1.},{0.,1.,3.}}},
  NKTestNear["FindDominantEigenvalue/PowerMethod: dominant eigenvalue = 3+sqrt(2)",
    FindDominantEigenvalue[M]["Eigenvalue"], N[3+Sqrt[2]], "AbsoluteTolerance" -> 1*^-5];
  NKTestNear["FindDominantEigenvalue: eigenvector unit norm",
    Norm[FindDominantEigenvalue[M]["Eigenvector"]], 1., "AbsoluteTolerance" -> 1*^-8];
  With[{all = FindMatrixEigenvalues[M, Method -> "JacobiMethod"]},
    NKTestNear["FindMatrixEigenvalues/Jacobi: eigenvalue sum = trace = 9",
      Total[all["Eigenvalues"]], 9., "AbsoluteTolerance" -> 1*^-8]];
  NKTestNear["FindMatrixEigenvalues/QRIteration: eigenvalue sum = trace = 9",
    Total[FindMatrixEigenvalues[M, Method -> "QRIteration"]["Eigenvalues"]],
    9., "AbsoluteTolerance" -> 1*^-5]
]

With[{M = {{2.,1.},{1.,3.}}},
  With[{res = FindDominantEigenvalue[M]},
    NKTestNear["FindDominantEigenvalue: Rayleigh quotient check",
      res["Eigenvector"].M.res["Eigenvector"] / (res["Eigenvector"].res["Eigenvector"]),
      res["Eigenvalue"], "AbsoluteTolerance" -> 1*^-8]]
]

NKTestFails["FindMatrixEigenvalues/JacobiMethod: non-symmetric returns $Failed",
  FindMatrixEigenvalues[{{1.,2.},{3.,4.}}, Method -> "JacobiMethod"]]

NKEndSuite[]
