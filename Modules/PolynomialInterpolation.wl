(* PolynomialInterpolation.wl -- Polynomial and spline interpolation
   PolynomialInterpolate : Newton divided differences, Lagrange
   HermiteInterpolate    : Hermite interpolation (matches values and derivatives)
   SplineInterpolate     : linear, quadratic, natural cubic splines

   Private helper naming convention: camelCase prefixed by "interp" or "spline"
*)

Options[PolynomialInterpolate] = {
  Method -> "Newton"
}

Options[SplineInterpolate] = {
  Degree            -> 3,
  BoundaryCondition -> "Natural"
}

(* 鈹€鈹€ PolynomialInterpolate dispatcher (with validation) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)

PolynomialInterpolate[xs_?VectorQ, ys_?VectorQ, xq_, opts:OptionsPattern[]] :=
  Module[{method},
    If[Length@xs =!= Length@ys,
      Message[PolynomialInterpolate::baddims, Length@xs, Length@ys];
      Return[$Failed, Module]];
    If[Length@xs =!= Length@Union[xs],
      Message[PolynomialInterpolate::dupnodes];
      Return[$Failed, Module]];
    method = OptionValue[Method];
    Switch[method,
      "Lagrange", interpLagrange[N@xs, N@ys, N@xq],
      "Newton",   interpNewton[N@xs, N@ys, N@xq],
      _,
        Message[PolynomialInterpolate::badmethod, method];
        interpNewton[N@xs, N@ys, N@xq]
    ]
  ]

(* 鈹€鈹€ Lagrange Interpolation 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
interpLagrange[xs_, ys_, xq_] :=
  Module[{n = Length[xs]},
    Sum[
      ys[[k]] * Product[If[j == k, 1, (xq-xs[[j]])/(xs[[k]]-xs[[j]])], {j, n}],
      {k, n}
    ]
  ]

(* 鈹€鈹€ Newton Divided-Difference Interpolation (Horner evaluation) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
interpNewton[xs_, ys_, xq_] :=
  Module[{n = Length[xs], dd = interpDividedDiffCoeffs[xs, ys]},
    Fold[Function[{acc, k}, dd[[k]] + (xq - xs[[k]]) acc],
      dd[[n]], Range[n-1, 1, -1]]
  ]

interpDividedDiffCoeffs[xs_, ys_] :=
  Module[{n = Length[xs], F = N@ys},
    Do[
      Do[F[[j]] = (F[[j]] - F[[j-1]]) / (xs[[j+i-1]] - xs[[j-1]]), {j, n, i+1, -1}],
      {i, n-1}
    ];
    F
  ]

(* 鈹€鈹€ Hermite Interpolation (matches function values and derivatives) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
HermiteInterpolate[xs_?VectorQ, ys_?VectorQ, dys_?VectorQ, xq_] :=
  Module[{n = Length[xs], zs, fz, m, T},
    zs = Riffle[N@xs, N@xs];     (* duplicate each node *)
    fz = Riffle[N@ys, N@ys];
    m  = 2n;
    T  = ConstantArray[0., {m, m}];
    T[[All, 1]] = fz;
    Do[
      Do[
        T[[k, j]] = If[zs[[k+j-1]] === zs[[k]],
          N@dys[[(k+1)/2]],     (* derivative at repeated node *)
          (T[[k+1, j-1]] - T[[k, j-1]]) / (zs[[k+j-1]] - zs[[k]])
        ],
        {k, m-j+1}
      ],
      {j, 2, m}
    ];
    Fold[Function[{acc, k}, T[[1,k]] + (xq - zs[[k]]) acc],
      T[[1,m]], Range[m-1, 1, -1]]
  ]

(* 鈹€鈹€ SplineInterpolate dispatcher 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)

SplineInterpolate[xs_?VectorQ, ys_?VectorQ, xq_, opts:OptionsPattern[]] :=
  Module[{deg = OptionValue[Degree], bc = OptionValue[BoundaryCondition], minPts},
    If[Length@xs =!= Length@ys,
      Message[SplineInterpolate::toofewpoints, deg+1, deg, Length@xs];
      Return[$Failed, Module]];
    minPts = deg + 1;   (* need at least deg+1 points *)
    If[Length@xs < minPts,
      Message[SplineInterpolate::toofewpoints, minPts, deg, Length@xs];
      Return[$Failed, Module]];
    Switch[deg,
      1, splineLinear[N@xs, N@ys, N@xq],
      2, splineQuadratic[N@xs, N@ys, N@xq],
      _, splineCubic[N@xs, N@ys, N@xq, bc]
    ]
  ]

(* Find interval index k such that xs[[k]] <= xq < xs[[k+1]] *)
splineInterval[xs_, xq_] :=
  Clip[LengthWhile[xs, # <= xq&], {1, Length[xs]-1}]

(* 鈹€鈹€ Piecewise Linear Spline 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
splineLinear[xs_, ys_, xq_] :=
  Module[{k = splineInterval[xs, xq], h},
    h = xs[[k+1]] - xs[[k]];
    ys[[k]] + (ys[[k+1]] - ys[[k]]) / h (xq - xs[[k]])
  ]

(* 鈹€鈹€ Piecewise Quadratic Spline 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
splineQuadratic[xs_, ys_, xq_] :=
  Module[{n = Length[xs]-1, h = Differences[xs],
          slopes = Differences[ys]/Differences[xs], c, A, rhs, k, dx},
    A   = DiagonalMatrix[ConstantArray[2., n], 0] +
          DiagonalMatrix[ConstantArray[1., n-1], 1] +
          DiagonalMatrix[ConstantArray[1., n-1], -1];
    rhs = 6 Differences[slopes];
    c   = Insert[laGaussElim[N@A[[2;;, 2;;-1]], N@rhs[[2;;]], True], 0., 1];
    k = splineInterval[xs, xq]; dx = xq - xs[[k]];
    ys[[k]] + slopes[[k]] dx + c[[k]]/2 dx^2
  ]

(* 鈹€鈹€ Natural Cubic Spline (zero second derivatives at endpoints) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€ *)
(* The spline moment system is symmetric tridiagonal; use Thomas algorithm.
   lower = upper = h[[2..n-1]],  diag = 2(h[[1..n-1]] + h[[2..n]]). *)
splineCubic[xs_, ys_, xq_, bc_:"Natural"] :=
  Module[{n = Length[xs]-1, h = Differences[xs],
          slopes = Differences[ys]/Differences[xs], M, rhs,
          lower, diag, upper},
    lower = N@h[[2;;n-1]];    (* sub-diagonal, length n-2 *)
    diag  = N@(2(h[[1;;n-1]] + h[[2;;n]]));  (* diagonal, length n-1 *)
    upper = N@h[[2;;n-1]];    (* super-diagonal, length n-2 *)
    rhs   = N@(6 Differences[slopes]);
    M = Switch[bc,
      "Natural",  Join[{0.}, laThomasAlgorithm[lower, diag, upper, rhs], {0.}],
      _,          Join[{0.}, laThomasAlgorithm[lower, diag, upper, rhs], {0.}]
    ];
    Module[{k = splineInterval[xs, xq], hk},
      hk = h[[k]];
      M[[k]]   (xs[[k+1]] - xq)^3 / (6hk) +
      M[[k+1]] (xq - xs[[k]])^3    / (6hk) +
      (ys[[k]]   / hk - M[[k]]   hk/6) (xs[[k+1]] - xq) +
      (ys[[k+1]] / hk - M[[k+1]] hk/6) (xq - xs[[k]])
    ]
  ]

