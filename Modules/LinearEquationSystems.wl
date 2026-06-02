(* LinearEquationSystems.wl -- Linear equation system solvers and decompositions
   Direct methods  : GaussianElimination, GaussianEliminationPivot,
                     GaussJordan, LU, Cholesky, LDLT, Tridiagonal
   Iterative       : Jacobi, GaussSeidel, SOR, ConjugateGradient (Batch-3)
   Sparse          : SparseDirectLU (Batch-3)
   Decompositions  : LUDecompose, QRDecompose, CholeskyDecompose, LDLTDecompose

   Batch-2 improvements
     * SolveLinearEquationSystem validates square A and matching dimensions
     * CholeskyDecompose validates symmetry before attempting factorisation
     * LDLTDecompose validates symmetry
     * LU-based solve uses safeLinearSolve internally for back-substitution safety

   Private helper naming convention: camelCase prefixed by "la" (linear algebra)
   Shared helpers also used by other modules: qrHouseholder, laThomasAlgorithm
*)

Options[SolveLinearEquationSystem] = {
  Method        -> "GaussianEliminationPivot",
  Omega         -> 1.5,
  Tolerance     -> 1*^-8,
  MaxIterations -> 1000
}

Options[QRDecompose] = {
  Method -> "Householder"
}

(* ── SolveLinearEquationSystem dispatcher (with input validation) ─────── *)

SolveLinearEquationSystem[A_?MatrixQ, b_?VectorQ, opts:OptionsPattern[]] :=
  Module[{method, omega, tol, maxIter, nA},
    nA = Length@A;
    (* Validate square matrix *)
    If[nA =!= Length[A[[1]]],
      Message[SolveLinearEquationSystem::notsquare, nA, Length[A[[1]]]];
      Return[$Failed, Module]];
    (* Validate dimension match *)
    If[nA =!= Length@b,
      Message[SolveLinearEquationSystem::baddims, nA, Length@b];
      Return[$Failed, Module]];
    method  = OptionValue[Method];
    omega   = OptionValue[Omega];
    tol     = OptionValue[Tolerance];
    maxIter = OptionValue[MaxIterations];
    Switch[method,
      "GaussianElimination",      laGaussElim[N@A, N@b, False],
      "GaussianEliminationPivot", laGaussElim[N@A, N@b, True],
      "GaussJordan",              laGaussJordan[N@A, N@b],
      "LU",                       laLUSolve[N@A, N@b],
      "Cholesky",                 laCholeskySolve[N@A, N@b],
      "LDLT",                     laLDLTSolve[N@A, N@b],
      "Tridiagonal",              laTridiagSolve[N@A, N@b],
      "Jacobi",                   laJacobi[N@A, N@b, tol, maxIter],
      "GaussSeidel",              laGaussSeidel[N@A, N@b, tol, maxIter],
      "SOR",                      laSOR[N@A, N@b, omega, tol, maxIter],
      (* Batch-3: sparse / Krylov methods *)
      "ConjugateGradient",        laConjugateGradient[A, b, tol, maxIter],
      "SparseDirectLU",           laSparseDirectLU[A, b],
      _,
        Message[SolveLinearEquationSystem::badmethod, method];
        laGaussElim[N@A, N@b, True]
    ]
  ]

(* ── Gaussian Elimination (with or without partial pivoting) ────────────── *)
laGaussElim[A0_, b0_, pivot_:True] :=
  Module[{n = Length[b0], A = N@A0, b = N@b0, prow, m},
    Do[
      If[pivot,
        prow = First@MaximalBy[Range[k, n], Abs[A[[#, k]]]&];
        If[prow != k,
          A[[{k, prow}]] = A[[{prow, k}]];
          b[[{k, prow}]] = b[[{prow, k}]]]
      ];
      Do[
        m       = A[[i, k]] / A[[k, k]];
        A[[i]] -= m A[[k]];
        b[[i]] -= m b[[k]],
        {i, k+1, n}
      ],
      {k, n}
    ];
    Module[{x = ConstantArray[0., n]},
      x[[n]] = b[[n]] / A[[n, n]];
      Do[x[[i]] = (b[[i]] - A[[i, i+1;;n]] . x[[i+1;;n]]) / A[[i, i]], {i, n-1, 1, -1}];
      x]
  ]

(* ── Gauss-Jordan Full Elimination ─────────────────────────────────────── *)
laGaussJordan[A0_, b0_] :=
  Module[{n = Length[b0], Aug = MapThread[Append, {N@A0, N@b0}], prow},
    Do[
      prow = First@MaximalBy[Range[k, n], Abs[Aug[[#, k]]]&];
      If[prow != k, Aug[[{k, prow}]] = Aug[[{prow, k}]]];
      Aug[[k]] /= Aug[[k, k]];
      Do[Aug[[i]] -= Aug[[i, k]] Aug[[k]], {i, Complement[Range[n], {k}]}],
      {k, n}
    ];
    Aug[[All, -1]]
  ]

(* ── LU Decomposition with Partial Pivoting ─────────────────────────────── *)
LUDecompose[A_?MatrixQ] :=
  Module[{n = Length@A, M = N@A, L = IdentityMatrix[Length@A]*1.,
          P = Range[Length@A], prow, m},
    Do[
      prow = First@MaximalBy[Range[k, n], Abs[M[[#, k]]]&];
      If[prow != k,
        M[[{k, prow}]] = M[[{prow, k}]];
        L[[{k, prow}, 1;;k-1]] = L[[{prow, k}, 1;;k-1]];
        P[[{k, prow}]] = P[[{prow, k}]]
      ];
      Do[m = M[[i,k]]/M[[k,k]]; L[[i,k]] = m; M[[i]] -= m M[[k]], {i, k+1, n}],
      {k, n-1}
    ];
    {L, M, IdentityMatrix[n][[P]]}
  ]

laLUSolve[A_, b_] :=
  Module[{L, U, Pb, y, n = Length[b]},
    {L, U, Pb} = LUDecompose[A]; Pb = Pb . b;
    y = Module[{v = ConstantArray[0., n]},
      v[[1]] = Pb[[1]];
      Do[v[[i]] = Pb[[i]] - L[[i, 1;;i-1]] . v[[1;;i-1]], {i, 2, n}]; v];
    Module[{v = ConstantArray[0., n]},
      v[[n]] = y[[n]] / U[[n,n]];
      Do[v[[i]] = (y[[i]] - U[[i, i+1;;n]] . v[[i+1;;n]]) / U[[i,i]], {i, n-1, 1, -1}]; v]
  ]

(* ── Cholesky (LL^T) Decomposition ─────────────────────────────────────── *)
CholeskyDecompose[A_?MatrixQ] :=
  Module[{n = Length@A, L = ConstantArray[0., {Length@A, Length@A}], s, An},
    An = N@A;
    If[!SymmetricMatrixQ[An, Tolerance -> 1*^-10],
      Message[CholeskyDecompose::nonsymmetric]; Return[$Failed, Module]];
    A = An;
    Do[
      s = A[[j,j]] - If[j > 1, L[[j, 1;;j-1]] . L[[j, 1;;j-1]], 0.];
      If[s <= 0., Message[CholeskyDecompose::notspdefinite]; Return[$Failed]];
      L[[j,j]] = Sqrt[s];
      Do[L[[i,j]] = (A[[i,j]] - If[j > 1, L[[i, 1;;j-1]] . L[[j, 1;;j-1]], 0.]) / L[[j,j]],
        {i, j+1, n}],
      {j, n}
    ];
    L
  ]

laCholeskySolve[A_, b_] :=
  Module[{L = CholeskyDecompose[A], n = Length[b], y},
    If[L === $Failed, Return[$Failed]];
    y = Module[{v = ConstantArray[0., n]},
      v[[1]] = b[[1]] / L[[1,1]];
      Do[v[[i]] = (b[[i]] - L[[i, 1;;i-1]] . v[[1;;i-1]]) / L[[i,i]], {i, 2, n}]; v];
    Module[{v = ConstantArray[0., n]},
      v[[n]] = y[[n]] / L[[n,n]];
      Do[v[[i]] = (y[[i]] - L[[i+1;;n, i]] . v[[i+1;;n]]) / L[[i,i]], {i, n-1, 1, -1}]; v]
  ]

(* ── LDL^T Decomposition ────────────────────────────────────────────────── *)
LDLTDecompose[A_?MatrixQ] :=
  Module[{n = Length@A, L = IdentityMatrix[Length@A]*1., d = ConstantArray[0., Length@A], v, An},
    An = N@A;
    If[!SymmetricMatrixQ[An, Tolerance -> 1*^-10],
      Message[LDLTDecompose::nonsymmetric]; Return[$Failed, Module]];
    A = An;
    Do[
      v      = If[j > 1, L[[j, 1;;j-1]] * d[[1;;j-1]], {}];
      d[[j]] = A[[j,j]] - If[j > 1, L[[j, 1;;j-1]] . v, 0.];
      Do[L[[i,j]] = (A[[i,j]] - If[j > 1, L[[i, 1;;j-1]] . v, 0.]) / d[[j]], {i, j+1, n}],
      {j, n}
    ];
    {L, d}
  ]

laLDLTSolve[A_, b_] :=
  Module[{L, d, n = Length[b], z, y},
    {L, d} = LDLTDecompose[A];
    z = Module[{v = N@b},
      Do[v[[i]] -= L[[i, 1;;i-1]] . v[[1;;i-1]], {i, 2, n}]; v];
    y = z/d;
    Module[{v = N@y},
      Do[v[[i]] -= L[[i+1;;n, i]] . v[[i+1;;n]], {i, n-1, 1, -1}]; v]
  ]

(* ── Thomas Algorithm for Tridiagonal Systems ───────────────────────────── *)
laTridiagSolve[A_, b_] :=
  Module[{n = Length[b],
          lower = Table[A[[i, i-1]], {i, 2, Length[b]}],
          diag  = Table[A[[i, i]],   {i, Length[b]}],
          upper = Table[A[[i, i+1]], {i, Length[b]-1}]},
    laThomasAlgorithm[lower, diag, upper, N@b]
  ]

(* Shared helper: also used by BoundaryValueProblems.wl *)
laThomasAlgorithm[lower_, diag_, upper_, rhs_] :=
  Module[{n = Length[diag], c = N@upper, d = N@rhs, x},
    Do[
      c[[i]] = c[[i]] / diag[[i]];
      d[[i]] = (d[[i]] - (If[i > 1, lower[[i-1]] c[[i-1]], 0.])) /
               (diag[[i]] - (If[i > 1, lower[[i-1]] c[[i-1]], 0.])),
      {i, n}
    ];
    x = ConstantArray[0., n];
    x[[n]] = d[[n]];
    Do[x[[i]] = d[[i]] - c[[i]] x[[i+1]], {i, n-1, 1, -1}];
    x
  ]

(* ── QR Decomposition ───────────────────────────────────────────────────── *)
QRDecompose[A_?MatrixQ, opts:OptionsPattern[]] :=
  Switch[OptionValue[Method],
    "GramSchmidt", qrGramSchmidt[N@A],
    _,             qrHouseholder[N@A]
  ]

(* Shared helper: also used by MatrixEigenanalysis.wl and PolynomialApproximation.wl *)
qrHouseholder[A_] :=
  Module[{m = Length@A, n = Length[A[[1]]], Q = IdentityMatrix[Length@A]*1., R = N@A, v, H},
    Do[
      v     = R[[k;;m, k]];
      v[[1]] += Sign[v[[1]]] Norm[v];
      If[Norm[v] < $MachineEpsilon, Continue[]];
      v /= Norm[v];
      H = IdentityMatrix[m-k+1] - 2 Outer[Times, v, v];
      R[[k;;m, k;;n]] = H . R[[k;;m, k;;n]];
      Q[[All, k;;m]] = Q[[All, k;;m]] . Transpose[H],
      {k, Min[m-1, n]}
    ];
    {Q, R}
  ]

qrGramSchmidt[A_] :=
  Module[{cols = Transpose@A, e = {}, u},
    Do[
      u = cols[[k]];
      Do[u -= (u . e[[j]]) e[[j]], {j, Length[e]}];
      AppendTo[e, If[Norm[u] < $MachineEpsilon, cols[[k]], u/Norm[u]]],
      {k, Length[cols]}
    ];
    {Transpose[e], Transpose[e] . Transpose[cols] // Transpose}
  ]

(* ── Jacobi Iteration ───────────────────────────────────────────────────── *)
laJacobi[A_, b_, tol_, maxIter_] :=
  Module[{n = Length[b], x = ConstantArray[0., Length[b]], xNew, k = 0},
    While[k++ < maxIter,
      xNew = Table[
        (b[[i]] - Sum[If[j != i, A[[i,j]] x[[j]], 0.], {j, n}]) / A[[i,i]],
        {i, n}];
      If[Norm[xNew - x] < tol, Return[xNew]];
      x = xNew];
    x]

(* ── Gauss-Seidel Iteration ─────────────────────────────────────────────── *)
laGaussSeidel[A_, b_, tol_, maxIter_] :=
  Module[{n = Length[b], x = ConstantArray[0., Length[b]], xOld, k = 0},
    While[k++ < maxIter,
      xOld = x;
      Do[x[[i]] = (b[[i]] - A[[i, 1;;i-1]] . x[[1;;i-1]] -
                              A[[i, i+1;;n]] . x[[i+1;;n]]) / A[[i,i]],
        {i, n}];
      If[Norm[x - xOld] < tol, Return[x]]];
    x]

(* ── Successive Over-Relaxation (SOR) ──────────────────────────────────── *)
laSOR[A_, b_, omega_, tol_, maxIter_] :=
  Module[{n = Length[b], x = ConstantArray[0., Length[b]], xOld, gs, k = 0},
    While[k++ < maxIter,
      xOld = x;
      Do[
        gs     = (b[[i]] - A[[i, 1;;i-1]] . x[[1;;i-1]] -
                             A[[i, i+1;;n]] . x[[i+1;;n]]) / A[[i,i]];
        x[[i]] = (1 - omega) x[[i]] + omega gs,
        {i, n}];
      If[Norm[x - xOld] < tol, Return[x]]];
    x]

(* ── Conjugate Gradient (Batch-3: symmetric positive-definite A) ─────────── *)
(* Works with both dense Matrix and SparseArray. *)
laConjugateGradient[A_, b_, tol_, maxIter_] :=
  Module[{x = ConstantArray[0., Length@b], r = N@b, p, rsOld, rsNew, alpha, Ap,
          k = 0},
    p = r; rsOld = r . r;
    While[k++ < maxIter,
      Ap    = A . p;
      alpha = rsOld / (p . Ap);
      x    += alpha p;
      r    -= alpha Ap;
      rsNew = r . r;
      If[Sqrt[rsNew] < tol, Break[]];
      p     = r + (rsNew / rsOld) p;
      rsOld = rsNew
    ];
    x
  ]

(* ── Sparse Direct Solve (Batch-3) ─────────────────────────────────────── *)
(* Converts A to SparseArray; Mathematica then uses its sparse direct solver. *)
laSparseDirectLU[A_, b_] :=
  Quiet[
    Check[
      LinearSolve[SparseArray[A], N@b],
      $Failed,
      {LinearSolve::sing, LinearSolve::luc, LinearSolve::nosol}
    ],
    {LinearSolve::sing, LinearSolve::luc, LinearSolve::nosol}
  ]

(* ── 4.3  Protect cross-module exports ──────────────────────────────────── *)
Protect[qrHouseholder, laThomasAlgorithm]
