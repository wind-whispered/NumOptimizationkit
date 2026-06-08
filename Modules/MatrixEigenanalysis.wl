(* MatrixEigenanalysis.wl -- Matrix eigenvalue computation
   FindDominantEigenvalue  : PowerMethod, InversePowerMethod
   FindMatrixEigenvalues   : JacobiMethod, QRIteration

   Private helper naming convention: camelCase prefixed by "eig"
   Depends on: qrHouseholder (defined in LinearEquationSystems.wl)
*)

Options[FindDominantEigenvalue] = {
  Method        -> "PowerMethod",
  Shift         -> 0,
  Tolerance     -> 1*^-8,
  MaxIterations -> 500
}

Options[FindMatrixEigenvalues] = {
  Method        -> "QRIteration",
  Tolerance     -> 1*^-8,
  MaxIterations -> 500
}

(* ── FindDominantEigenvalue dispatcher ──────────────────────────────────── *)

FindDominantEigenvalue[A_?MatrixQ, opts:OptionsPattern[]] :=
  Module[{method, tol, maxIter, shift},
    method  = OptionValue[Method];
    tol     = OptionValue[Tolerance];
    maxIter = OptionValue[MaxIterations];
    shift   = OptionValue[Shift];
    Switch[method,
      "PowerMethod",        eigPowerMethod[N@A, tol, maxIter],
      "InversePowerMethod", eigInversePowerMethod[N@A, N@shift, tol, maxIter],
      _,
        Message[FindDominantEigenvalue::badmethod, method];
        eigPowerMethod[N@A, tol, maxIter]
    ]
  ]

(* ── Power Method ───────────────────────────────────────────────────────── *)
eigPowerMethod[A_, tol_, maxIter_] :=
  Module[{n = Length@A, v = N@UnitVector[Length@A, 1], vNew,
          lambda = 0., lambdaOld, k = 0, converged = False},
    Catch[
      While[k++ < maxIter,
        vNew      = A . v;
        lambdaOld = lambda;
        lambda    = vNew[[First@Ordering[Abs@vNew, -1]]];
        vNew     /= lambda;
        If[Abs[lambda - lambdaOld] < tol, converged = True; v = vNew; Throw[Null]];
        v = vNew
      ]
    ];
    lambda = (A . v) . v / (v . v);   (* Rayleigh quotient for sign *)
    <|"Eigenvalue"  -> lambda,
      "Eigenvector" -> v / Norm[v],
      "Iterations"  -> k,
      "Converged"   -> converged|>
  ]

(* ── Shifted Inverse Power Method (finds eigenvalue nearest Shift) ───────── *)
eigInversePowerMethod[A_, shift_, tol_, maxIter_] :=
  Module[{B = A - shift IdentityMatrix[Length@A], result},
    result = eigPowerMethod[Inverse[B], tol, maxIter];
    <|"Eigenvalue"  -> shift + 1./result["Eigenvalue"],
      "Eigenvector" -> result["Eigenvector"],
      "Iterations"  -> result["Iterations"],
      "Converged"   -> result["Converged"]|>
  ]

(* ── FindMatrixEigenvalues dispatcher (with symmetry validation) ─────────── *)

FindMatrixEigenvalues[A_?MatrixQ, opts:OptionsPattern[]] :=
  Module[{method, tol, maxIter, An},
    An = N@A;
    method  = OptionValue[Method];
    tol     = OptionValue[Tolerance];
    maxIter = OptionValue[MaxIterations];
    (* JacobiMethod requires a symmetric matrix *)
    If[method === "JacobiMethod" && !SymmetricMatrixQ[An, Tolerance -> 1*^-10],
      Message[FindMatrixEigenvalues::nonsymmetric];
      Return[$Failed, Module]];
    Switch[method,
      "JacobiMethod", eigJacobi[An, tol, maxIter],
      "QRIteration",  eigQRIteration[An, tol, maxIter],
      _,
        Message[FindMatrixEigenvalues::badmethod, method];
        eigQRIteration[An, tol, maxIter]
    ]
  ]

(* ── Jacobi Rotation Method (symmetric matrices only) ───────────────────── *)
eigJacobi[A0_, tol_, maxIter_] :=
  Module[{n = Length@A0, A = N@A0, V = IdentityMatrix[Length@A0]*1.,
          p, q, theta, c, s, G, k = 0, offMax},
    Catch[
      While[k++ < maxIter,
        offMax = Max[Abs@Flatten@Table[If[i != j, A[[i,j]], 0.], {i,n},{j,n}]];
        If[offMax < tol, Throw[Null]];
        {p, q} = First@MaximalBy[
          Flatten[Table[{i,j}, {i,n},{j,i+1,n}], 1],
          Abs[A[[#[[1]], #[[2]]]]]&];
        theta = If[Abs[A[[p,p]] - A[[q,q]]] < $MachineEpsilon,
          Pi/4,
          ArcTan[2 A[[p,q]] / (A[[p,p]] - A[[q,q]])] / 2];
        c = Cos[theta]; s = Sin[theta];
        G = IdentityMatrix[n]*1.;
        G[[p,p]] = c;  G[[q,q]] = c;
        G[[p,q]] = s;  G[[q,p]] = -s;
        A = Transpose[G] . A . G;
        V = V . G
      ]
    ];
    <|"Eigenvalues"  -> Diagonal[A],
      "Eigenvectors" -> Transpose[V],
      "Iterations"   -> k|>
  ]

(* ── QR Iteration with Origin Shift ─────────────────────────────────────── *)
eigQRIteration[A0_, tol_, maxIter_] :=
  Module[{n = Length@A0, A = eigHessenberg[N@A0],
          Q0 = IdentityMatrix[Length@A0]*1., Q, R, shift, k = 0},
    Catch[
      While[k++ < maxIter,
        shift = A[[-1,-1]];
        {Q, R} = qrHouseholder[A - shift IdentityMatrix[n]];
        A  = R . Q + shift IdentityMatrix[n];
        Q0 = Q0 . Q;
        If[Max[Abs@Flatten@Table[If[i>j, A[[i,j]], 0.], {i,n},{j,n}]] < tol,
          Throw[Null]]
      ]
    ];
    <|"Eigenvalues"  -> Diagonal[A],
      "Eigenvectors" -> Transpose[Q0],
      "Iterations"   -> k|>
  ]

(* ── Hessenberg Reduction (Householder similarity) ──────────────────────── *)
eigHessenberg[A_] :=
  Module[{n = Length@A, M = N@A, v, H},
    Do[
      v     = M[[k+1;;n, k]];
      v[[1]] += Sign[v[[1]]] Norm[v];
      If[Norm[v] < $MachineEpsilon, Continue[]];
      v /= Norm[v];
      H = IdentityMatrix[n-k] - 2 Outer[Times, v, v];
      M[[k+1;;n, k;;n]] = H . M[[k+1;;n, k;;n]];
      M[[1;;n, k+1;;n]] = M[[1;;n, k+1;;n]] . Transpose[H],
      {k, n-2}
    ];
    M
  ]
