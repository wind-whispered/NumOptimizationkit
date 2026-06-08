(* Performance_LinearAlgebra.wl -- Performance regression for linear algebra *)

NKBeginSuite["Linear Algebra Performance", "E", "LA_Timing"]

With[{A=N@RandomReal[{0.,1.},{200,200}]},
  NKTestPerformanceRatio["Performance/LUDecompose: 200x200",
    LUDecompose[A],"LU_200x200","MaxRatio"->3.0]
]

With[{n=1000,
      A=SparseArray[{Band[{1,1}]->4.,Band[{1,2}]->-1.,Band[{2,1}]->-1.},{n,n}],
      b=ConstantArray[1.,n]},
  NKTestPerformanceRatio["Performance/ConjugateGradient: 1000x1000 tridiagonal",
    SolveLinearEquationSystem[Normal@A,b,Method->"ConjugateGradient",Tolerance->1*^-8],
    "CG_1000x1000","MaxRatio"->3.0]
]

NKTest["Performance/GaussSeidel: 100x100 diagonal dominant < 3s",
  First[AbsoluteTiming[
    SolveLinearEquationSystem[
      DiagonalMatrix[ConstantArray[4.,100]]+
      DiagonalMatrix[ConstantArray[-1.,99],1]+
      DiagonalMatrix[ConstantArray[-1.,99],-1],
      ConstantArray[1.,100],
      Method->"GaussSeidel",Tolerance->1*^-8,MaxIterations->2000]]] < 3.]

NKEndSuite[]
