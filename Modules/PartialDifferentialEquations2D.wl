(* PartialDifferentialEquations2D.wl -- 2D PDE finite difference solvers
   SolveEllipticPDE    -- Poisson equation  nabla^2 u = f(x,y)
   SolveParabolicPDE2D -- Heat equation     u_t = alpha^2 (u_xx + u_yy)
   SolveHyperbolicPDE2D -- Wave equation    u_tt = c^2 (u_xx + u_yy)

   All methods use uniform Cartesian grids with Dirichlet boundary conditions.
   The 5-point Laplacian stencil is assembled as a SparseArray for the elliptic
   solve. The 2D parabolic uses Peaceman-Rachford ADI (alternating direction
   implicit), which decomposes into independent tridiagonal sweeps. The 2D
   hyperbolic uses an explicit second-order scheme.

   Private helper naming convention: camelCase prefixed by "pde2"
*)

(* ══════════════════════════════════════════════════════════════════════════
   SolveEllipticPDE  --  nabla^2 u = f(x,y)
   ══════════════════════════════════════════════════════════════════════════ *)

Options[SolveEllipticPDE] = {
  Method -> "DirectSparse"
}

SolveEllipticPDE::usage =
"SolveEllipticPDE[f, {x,xL,xR,nx}, {y,yL,yR,ny}, {bcLeft,bcRight,bcBottom,bcTop}]
solves the Poisson equation  nabla^2 u = f(x,y)  on [xL,xR]*[yL,yR].
  f[x,y]          -- forcing function (RHS)
  bcLeft[y]       -- u(xL, y)    bcRight[y] -- u(xR, y)
  bcBottom[x]     -- u(x, yL)   bcTop[x]   -- u(x, yR)
Returns <|\"SpatialGridX\"->xs, \"SpatialGridY\"->ys, \"Solution\"->matrix|>
where matrix[[j,i]] = u(x_i, y_j) on the full (nx+2)*(ny+2) grid."

SolveEllipticPDE[f_, {x_, xL_:0, xR_?NumericQ, nx_Integer},
    {y_, yL_:0, yR_?NumericQ, ny_Integer},
    {bcLeft_, bcRight_, bcBottom_, bcTop_}, opts:OptionsPattern[]] :=
  Module[{hx, hy, xs, ys, N2, idx, coeffRules, A, b, k, xi, yj,
          uInner, uFull, row},
    hx = (xR - xL) / (nx + 1);
    hy = (yR - yL) / (ny + 1);
    xs = xL + Range[0, nx+1] hx;   (* full x-grid including boundaries *)
    ys = yL + Range[0, ny+1] hy;   (* full y-grid including boundaries *)
    N2 = nx * ny;
    idx[i_, j_] := (j-1)*nx + i;   (* linear index for interior node (i,j) *)
    (* ── Build sparse coefficient matrix ── *)
    coeffRules = Flatten[Table[
      Module[{k = idx[i, j]},
        {
          {k, k} -> -2./hx^2 - 2./hy^2,
          If[i > 1,  {k, idx[i-1, j]} -> 1./hx^2, Nothing],
          If[i < nx, {k, idx[i+1, j]} -> 1./hx^2, Nothing],
          If[j > 1,  {k, idx[i, j-1]} -> 1./hy^2, Nothing],
          If[j < ny, {k, idx[i, j+1]} -> 1./hy^2, Nothing]
        }
      ],
      {j, ny}, {i, nx}
    ], 2];
    A = SparseArray[coeffRules, {N2, N2}];
    (* ── Build right-hand side ── *)
    b = Flatten[Table[
      Module[{xi = xs[[i+1]], yj = ys[[j+1]], rhs},
        rhs = N[f[xi, yj]];
        If[i == 1,   rhs -= N[bcLeft[yj]]   / hx^2];
        If[i == nx,  rhs -= N[bcRight[yj]]  / hx^2];
        If[j == 1,   rhs -= N[bcBottom[xi]] / hy^2];
        If[j == ny,  rhs -= N[bcTop[xi]]    / hy^2];
        rhs
      ],
      {j, ny}, {i, nx}
    ]];
    (* ── Solve sparse system via our own Conjugate Gradient method ── *)
    (* The 5-point Laplacian (with Dirichlet BCs) produces a symmetric
       negative-definite matrix A.  Negate to make it positive-definite for CG.
       For large problems, use laSparseDirectLU (our sparse interface) instead. *)
    uInner = If[N2 <= 400,
      laConjugateGradient[-Normal@A, -N@b, 1*^-10, 2 N2],  (* dense CG on -A *)
      laSparseDirectLU[A, N@b]  (* our sparse interface for large systems *)
    ];
    (* ── Reshape and add boundaries ── *)
    uFull = Table[
      Which[
        i == 0,    N[bcLeft[ys[[j+1]]]],
        i == nx+1, N[bcRight[ys[[j+1]]]],
        j == 0,    N[bcBottom[xs[[i+1]]]],
        j == ny+1, N[bcTop[xs[[i+1]]]],
        True,      uInner[[idx[i, j]]]
      ],
      {j, 0, ny+1}, {i, 0, nx+1}
    ];
    <|"SpatialGridX" -> xs, "SpatialGridY" -> ys, "Solution" -> uFull|>
  ]

(* ══════════════════════════════════════════════════════════════════════════
   SolveParabolicPDE2D  --  u_t = alpha^2 (u_xx + u_yy)
   Peaceman-Rachford ADI (alternating direction implicit, order 2 in t and x,y)
   Each half-step reduces to nx (or ny) independent tridiagonal systems.
   ══════════════════════════════════════════════════════════════════════════ *)

Options[SolveParabolicPDE2D] = {
  Method -> "ADI"
}

SolveParabolicPDE2D::usage =
"SolveParabolicPDE2D[alpha2, {x,xL,xR,nx}, {y,yL,yR,ny}, {t,t0,tEnd,nt},
  ic, {bcLeft,bcRight,bcBottom,bcTop}]
solves the 2D heat equation  u_t = alpha2*(u_xx+u_yy)  using the
Peaceman-Rachford ADI method.
  ic[x,y]    -- initial condition u(x,y,t0)
  bcXxx[...]  -- boundary conditions (same convention as SolveEllipticPDE)
Returns <|\"SpatialGridX\"->xs, \"SpatialGridY\"->ys, \"TimeGrid\"->ts, \"Solution\"->{u0,u1,...}|>
where each u_k is a (ny+2)*(nx+2) matrix."

SolveParabolicPDE2D[alpha2_?NumericQ,
    {x_, xL_:0, xR_?NumericQ, nx_Integer},
    {y_, yL_:0, yR_?NumericQ, ny_Integer},
    {t_, t0_:0, tEnd_?NumericQ, nt_Integer},
    ic_, {bcLeft_, bcRight_, bcBottom_, bcTop_}, opts:OptionsPattern[]] :=
  Module[{hx, hy, dt, xs, ys, ts, rx, ry, u, U},
    hx = (xR - xL)/(nx+1); hy = (yR - yL)/(ny+1);
    dt = (tEnd - t0)/nt;
    xs = xL + Range[0, nx+1] hx;
    ys = yL + Range[0, ny+1] hy;
    ts = t0 + Range[0, nt] dt;
    rx = alpha2 dt / (2 hx^2);   (* half-step diffusion numbers *)
    ry = alpha2 dt / (2 hy^2);
    (* Initialise full grid (including boundaries) *)
    u = Table[N[ic[xs[[i+1]], ys[[j+1]]]], {j, 0, ny+1}, {i, 0, nx+1}];
    (* Set boundary rows/columns in initial condition *)
    u = pde2SetBC2D[u, xs, ys, ts[[1]], bcLeft, bcRight, bcBottom, bcTop, nx, ny];
    U = {u};
    Do[
      u = pde2ADIStep[u, xs, ys, ts[[step]], dt, rx, ry, nx, ny,
                      bcLeft, bcRight, bcBottom, bcTop];
      AppendTo[U, u],
      {step, 2, nt+1}
    ];
    <|"SpatialGridX" -> xs, "SpatialGridY" -> ys,
      "TimeGrid"     -> ts, "Solution"     -> U|>
  ]

(* Set Dirichlet boundary values on the 4 edges *)
pde2SetBC2D[u_, xs_, ys_, t_, bcL_, bcR_, bcB_, bcT_, nx_, ny_] :=
  Module[{v = u},
    Do[v[[j+1, 1]]    = N[bcL[ys[[j+1]]]]; (* left:   i=0 *)
       v[[j+1, nx+2]] = N[bcR[ys[[j+1]]]]; (* right:  i=nx+1 *)
       v[[1,    i+1]] = N[bcB[xs[[i+1]]]]; (* bottom: j=0 *)
       v[[ny+2, i+1]] = N[bcT[xs[[i+1]]]], (* top:    j=ny+1 *)
       {j, 0, ny+1}, {i, 0, nx+1}
    ];
    v
  ]

(* One Peaceman-Rachford ADI step:
   x-sweep: implicit in x (tridiagonal per row j=1..ny),  explicit in y
   y-sweep: explicit in x,                                 implicit in y (tridiagonal per col i=1..nx)
   Index convention: u[[j+1, i+1]] = u(x_i, y_j), i/j are 0-based for the grid. *)
pde2ADIStep[u_, xs_, ys_, t_, dt_, rx_, ry_, nx_, ny_,
            bcL_, bcR_, bcB_, bcT_] :=
  Module[{uStar = u, uNew = u, lower, diag, upper, rhs, val},

    (* ── x-sweep: for each row j (interior j=1..ny) ── *)
    Do[
      lower = ConstantArray[-rx, nx-1];
      diag  = ConstantArray[1 + 2rx, nx];
      upper = ConstantArray[-rx, nx-1];
      (* RHS: explicit y-direction *)
      rhs = Table[
        ry u[[j+1, i]] + (1 - 2ry) u[[j+1, i+1]] + ry u[[j+1, i+2]],
        {i, nx}   (* interior x: i=1..nx, index i+1 in u *)
      ] // N;
      (* Apply x-boundary contributions *)
      rhs[[1]]  += rx N[bcL[ys[[j+1]]]];
      rhs[[-1]] += rx N[bcR[ys[[j+1]]]];
      uStar[[j+1, 2;;nx+1]] = laThomasAlgorithm[lower, diag, upper, rhs],
      {j, ny}
    ];
    (* Fill x-boundaries in uStar *)
    Do[
      uStar[[j+1, 1]]    = N[bcL[ys[[j+1]]]];
      uStar[[j+1, nx+2]] = N[bcR[ys[[j+1]]]],
      {j, 0, ny+1}
    ];

    (* ── y-sweep: for each column i (interior i=1..nx) ── *)
    Do[
      lower = ConstantArray[-ry, ny-1];
      diag  = ConstantArray[1 + 2ry, ny];
      upper = ConstantArray[-ry, ny-1];
      (* RHS: explicit x-direction from uStar *)
      rhs = Table[
        rx uStar[[j+1, i]] + (1 - 2rx) uStar[[j+1, i+1]] + rx uStar[[j+1, i+2]],
        {j, ny}
      ] // N;
      rhs[[1]]  += ry N[bcB[xs[[i+1]]]];
      rhs[[-1]] += ry N[bcT[xs[[i+1]]]];
      uNew[[2;;ny+1, i+1]] = laThomasAlgorithm[lower, diag, upper, rhs],
      {i, nx}
    ];
    (* Fill y-boundaries in uNew *)
    Do[
      uNew[[1,    i+1]] = N[bcB[xs[[i+1]]]];
      uNew[[ny+2, i+1]] = N[bcT[xs[[i+1]]]],
      {i, 0, nx+1}
    ];
    uNew
  ]

(* ══════════════════════════════════════════════════════════════════════════
   SolveHyperbolicPDE2D  --  u_tt = c^2 (u_xx + u_yy)
   Explicit second-order finite difference.
   CFL condition: c^2 * dt^2 * (1/hx^2 + 1/hy^2) <= 1
   ══════════════════════════════════════════════════════════════════════════ *)

Options[SolveHyperbolicPDE2D] = {
  Method -> "ExplicitDifference"
}

SolveHyperbolicPDE2D::cfl =
  "CFL number `1` > 1. ExplicitDifference may be unstable."

SolveHyperbolicPDE2D::usage =
"SolveHyperbolicPDE2D[c2, {x,xL,xR,nx}, {y,yL,yR,ny}, {t,t0,tEnd,nt},
  ic, vic, {bcLeft,bcRight,bcBottom,bcTop}]
solves the 2D wave equation  u_tt = c2*(u_xx+u_yy).
  ic[x,y]   -- initial displacement  u(x,y,t0)
  vic[x,y]  -- initial velocity      u_t(x,y,t0)
CFL stability condition: c2*dt^2*(1/hx^2 + 1/hy^2) <= 1.
Returns <|\"SpatialGridX\"->xs, \"SpatialGridY\"->ys, \"TimeGrid\"->ts, \"Solution\"->{u0,u1,...}|>"

SolveHyperbolicPDE2D[c2_?NumericQ,
    {x_, xL_:0, xR_?NumericQ, nx_Integer},
    {y_, yL_:0, yR_?NumericQ, ny_Integer},
    {t_, t0_:0, tEnd_?NumericQ, nt_Integer},
    ic_, vic_, {bcLeft_, bcRight_, bcBottom_, bcTop_}, opts:OptionsPattern[]] :=
  Module[{hx, hy, dt, xs, ys, ts, r2x, r2y, cfl, u0, u1, u2, U},
    hx = (xR - xL)/(nx+1); hy = (yR - yL)/(ny+1);
    dt = (tEnd - t0)/nt;
    xs = xL + Range[0, nx+1] hx;
    ys = yL + Range[0, ny+1] hy;
    ts = t0 + Range[0, nt] dt;
    r2x = c2 dt^2 / hx^2;
    r2y = c2 dt^2 / hy^2;
    cfl = Sqrt[c2] dt Sqrt[1/hx^2 + 1/hy^2];
    If[cfl > 1., Message[SolveHyperbolicPDE2D::cfl, N[cfl]]];
    (* Initial displacement *)
    u0 = Table[N[ic[xs[[i+1]], ys[[j+1]]]], {j, 0, ny+1}, {i, 0, nx+1}];
    u0 = pde2SetBC2D[u0, xs, ys, ts[[1]], bcLeft, bcRight, bcBottom, bcTop, nx, ny];
    (* First step uses initial velocity (Taylor expansion) *)
    u1 = Table[
      Which[
        i == 0,    N[bcLeft[ys[[j+1]]]],
        i == nx+1, N[bcRight[ys[[j+1]]]],
        j == 0,    N[bcBottom[xs[[i+1]]]],
        j == ny+1, N[bcTop[xs[[i+1]]]],
        True,
          u0[[j+1,i+1]] + dt N[vic[xs[[i+1]],ys[[j+1]]]] +
          r2x/2 (u0[[j+1,i+2]] - 2u0[[j+1,i+1]] + u0[[j+1,i]]) +
          r2y/2 (u0[[j+2,i+1]] - 2u0[[j+1,i+1]] + u0[[j,i+1]])
      ],
      {j, 0, ny+1}, {i, 0, nx+1}
    ];
    U = {u0, u1};
    Do[
      u2 = Table[
        Which[
          i == 0,    N[bcLeft[ys[[j+1]]]],
          i == nx+1, N[bcRight[ys[[j+1]]]],
          j == 0,    N[bcBottom[xs[[i+1]]]],
          j == ny+1, N[bcTop[xs[[i+1]]]],
          True,
            2 u1[[j+1,i+1]] - u0[[j+1,i+1]] +
            r2x (u1[[j+1,i+2]] - 2u1[[j+1,i+1]] + u1[[j+1,i]]) +
            r2y (u1[[j+2,i+1]] - 2u1[[j+1,i+1]] + u1[[j,i+1]])
        ],
        {j, 0, ny+1}, {i, 0, nx+1}
      ];
      u0 = u1; u1 = u2; AppendTo[U, u2],
      {3, nt+1}
    ];
    <|"SpatialGridX" -> xs, "SpatialGridY" -> ys,
      "TimeGrid"     -> ts, "Solution"     -> U|>
  ]
