(* ScatteredInterpolation.wl -- Interpolation from scattered 2D data
   ScatteredInterpolate[points, values, q] evaluates an interpolant at
   query point q = {xq, yq} given data pairs {points[[i]], values[[i]]}.

   Methods
     "NearestNeighbor"  -- return the value of the closest data point (O(n))
     "InverseDistance"  -- IDW weighted average with 1/r^p weights (O(n))
     "ThinPlateSpline"  -- global RBF interpolation phi(r) = r^2 log(r) (O(n^3) setup)

   Options
     Method     -> "ThinPlateSpline" (default)
     IDWPower   -> 2   (exponent p for InverseDistance)
     Regularise -> 0.  (Tikhonov regularisation for ThinPlateSpline; helps noisy data)

   Private helper naming convention: camelCase prefixed by "si"
*)

Options[ScatteredInterpolate] = {
  Method     -> "ThinPlateSpline",
  IDWPower   -> 2,
  Regularise -> 0.
}

ScatteredInterpolate::usage =
"ScatteredInterpolate[points, values, q] interpolates scattered 2D data at
query point q = {xq, yq}.
  points  -- list of {xi, yi} data locations
  values  -- list of function values f(xi,yi)
  q       -- query point {xq, yq}
Returns the interpolated scalar value.

Options
  Method     -> \"ThinPlateSpline\" (default) | \"NearestNeighbor\" | \"InverseDistance\"
  IDWPower   -> 2   (inverse-distance exponent)
  Regularise -> 0.  (Tikhonov regularisation for ThinPlateSpline)"

ScatteredInterpolate::baddims =
  "points and values must have the same length. Got `1` and `2`."

ScatteredInterpolate::toofewpoints =
  "ThinPlateSpline requires at least 3 data points. Got `1`."

ScatteredInterpolate[points_?MatrixQ, values_?VectorQ, q_?VectorQ,
    opts:OptionsPattern[]] :=
  Module[{method, p, reg, pts, vals},
    If[Length@points =!= Length@values,
      Message[ScatteredInterpolate::baddims, Length@points, Length@values];
      Return[$Failed, Module]];
    method = OptionValue[Method];
    p      = OptionValue[IDWPower];
    reg    = OptionValue[Regularise];
    pts    = N@points; vals = N@values;
    Switch[method,
      "NearestNeighbor", siNearestNeighbor[pts, vals, N@q],
      "InverseDistance", siInverseDistance[pts, vals, N@q, N@p],
      "ThinPlateSpline", siThinPlateSpline[pts, vals, N@q, N@reg],
      _,
        Message[ScatteredInterpolate::badmethod, method];
        siThinPlateSpline[pts, vals, N@q, N@reg]
    ]
  ]

ScatteredInterpolate::badmethod =
  "Unknown Method \"`1`\". Falling back to \"ThinPlateSpline\"."

(* Vectorised form: evaluate at multiple query points *)
ScatteredInterpolate[points_?MatrixQ, values_?VectorQ, qs_?MatrixQ,
    opts:OptionsPattern[]] :=
  ScatteredInterpolate[points, values, #, opts]& /@ qs

(* ── Nearest-Neighbour Interpolation ────────────────────────────────────── *)
siNearestNeighbor[pts_, vals_, q_] :=
  Module[{dists = Norm[# - q]& /@ pts},
    vals[[First@Ordering[dists]]]
  ]

(* ── Inverse Distance Weighting ─────────────────────────────────────────── *)
siInverseDistance[pts_, vals_, q_, p_] :=
  Module[{dists = Norm[# - q]& /@ pts, weights, exactIdx},
    exactIdx = FirstPosition[dists, 0., {0}][[1]];
    If[exactIdx =!= 0, Return[vals[[exactIdx]], Module]];
    weights = 1 / dists^p;
    (weights . vals) / Total[weights]
  ]

(* ── Thin Plate Spline Interpolation ────────────────────────────────────── *)
(*   Interpolant: s(x) = a0 + a1*x + a2*y + Sum_i w_i phi(r_i)
     Radial basis: phi(r) = r^2 log(r)  (thin plate spline kernel)
     Fit solves the (n+3) x (n+3) symmetric system:
       [Phi + lambda*I   P ] [w]   [vals]
       [P^T              0 ] [a] = [0   ]
     where lambda = Regularise parameter for Tikhonov regularisation. *)

siThinPlateSpline[pts_, vals_, q_, reg_:0.] :=
  Module[{n, Phi, P, A, rhs, coeff, w, a, phiQ, pQ},
    n = Length[pts];
    If[n < 3,
      Message[ScatteredInterpolate::toofewpoints, n]; Return[$Failed, Module]];
    (* Thin plate spline kernel phi(r) = r^2 log(r), phi(0) = 0 *)
    tpsKernel[r_] := If[r < $MachineEpsilon, 0., r^2 Log[r]];
    (* Gram matrix Phi_ij = phi(||xi - xj||) *)
    Phi = Table[tpsKernel[Norm[pts[[i]] - pts[[j]]]], {i, n}, {j, n}];
    If[reg > 0., Phi += reg IdentityMatrix[n]];  (* regularisation *)
    (* Polynomial part P: columns [1, x, y] for each data point *)
    P = Prepend[#, 1.]& /@ pts;  (* n x 3 matrix *)
    (* Full system *)
    A = ArrayFlatten[{{Phi, P}, {Transpose@P, ConstantArray[0., {3, 3}]}}];
    rhs = Join[vals, {0., 0., 0.}];
    (* Solve *)
    (* Solve TPS system using our Gaussian elimination (general square matrix). *)
    (* TPS matrix is symmetric but indefinite, so we use partial-pivoting GE. *)
    coeff = Quiet[laGaussElim[N@A, N@rhs, True]];
    If[!VectorQ[coeff, NumericQ], Return[vals[[1]], Module]];
    w = coeff[[1;;n]]; a = coeff[[n+1;;n+3]];  (* weights and polynomial coeffs *)
    (* Evaluate at query point q = {xq, yq} *)
    phiQ = Table[tpsKernel[Norm[q - pts[[i]]]], {i, n}];
    pQ   = Prepend[q, 1.];
    w . phiQ + a . pQ
  ]
