(* PartialDifferentialEquations.wl -- PDE finite difference solvers
   SolveParabolicPDE  : heat equation  u_t = alpha^2 * u_xx
     Methods: ForwardDifference, BackwardDifference, CrankNicolson
   SolveHyperbolicPDE : wave equation  u_tt = c^2 * u_xx
     Method : ExplicitDifference

   Private helper naming convention: camelCase prefixed by "pde"
*)

Options[SolveParabolicPDE] = {
  Method -> "CrankNicolson"
}

Options[SolveHyperbolicPDE] = {
  Method -> "ExplicitDifference"
}

(* ── SolveParabolicPDE dispatcher ───────────────────────────────────────── *)

SolveParabolicPDE[alpha2_?NumericQ,
    {x_, xL_:0, xR_?NumericQ, nx_Integer},
    {t_, t0_:0, tEnd_?NumericQ, nt_Integer},
    ic_, {bc0_, bcL_}, opts:OptionsPattern[]] :=
  Module[{method, dx, dt, xs, ts},
    method = OptionValue[Method];
    dx = (xR - xL)/(nx + 1);
    dt = (tEnd - t0)/nt;
    xs = xL + Range[0, nx+1] dx;
    ts = t0 + Range[0, nt] dt;
    Switch[method,
      "ForwardDifference",  pdeForwardDifference[alpha2, xs, ts, ic, bc0, bcL, dx, dt],
      "BackwardDifference", pdeBackwardDifference[alpha2, xs, ts, ic, bc0, bcL, dx, dt],
      "CrankNicolson",      pdeCrankNicolson[alpha2, xs, ts, ic, bc0, bcL, dx, dt],
      _,
        Message[SolveParabolicPDE::badmethod, method];
        pdeCrankNicolson[alpha2, xs, ts, ic, bc0, bcL, dx, dt]
    ]
  ]

pdeResult[xs_, ts_, U_] :=
  <|"SpatialGrid" -> xs, "TimeGrid" -> ts, "Solution" -> U|>

(* ── Explicit Forward Difference (stable when r = alpha2*dt/dx^2 <= 0.5) ── *)
pdeForwardDifference[alpha2_, xs_, ts_, ic_, bc0_, bcL_, dx_, dt_] :=
  Module[{r = alpha2 dt/dx^2, nx = Length[xs]-2, nt = Length[ts]-1, U, u, uNew},
    If[r > 0.5, Message[SolveParabolicPDE::stiff, r]];
    u = ic /@ xs; U = {u};
    Do[
      uNew = u;
      uNew[[1]]  = bc0[ts[[j]]];
      uNew[[-1]] = bcL[ts[[j]]];
      Do[uNew[[i]] = u[[i]] + r (u[[i+1]] - 2u[[i]] + u[[i-1]]), {i, 2, nx+1}];
      u = uNew; AppendTo[U, u],
      {j, 2, nt+1}
    ];
    pdeResult[xs, ts, U]]

(* ── Implicit Backward Difference (unconditionally stable) ──────────────── *)
pdeBackwardDifference[alpha2_, xs_, ts_, ic_, bc0_, bcL_, dx_, dt_] :=
  Module[{r = alpha2 dt/dx^2, nx = Length[xs]-2, nt = Length[ts]-1,
          U, u, uNew, rhs, lower, diag, upper},
    (* tridiagonal coefficients — constant for all steps; computed here
       (not in the localization header) so nx is already concretely bound *)
    lower = ConstantArray[-r, nx-1];
    diag  = ConstantArray[1+2r, nx];
    upper = ConstantArray[-r, nx-1];
    u = ic /@ xs; U = {u};
    Do[
      rhs = u[[2;;nx+1]];
      rhs[[1]]  += r bc0[ts[[j]]];
      rhs[[-1]] += r bcL[ts[[j]]];
      (* Use our Thomas algorithm (laThomasAlgorithm) instead of LinearSolve *)
      uNew = Join[{bc0[ts[[j]]]}, laThomasAlgorithm[lower, diag, upper, rhs], {bcL[ts[[j]]]}];
      u = uNew; AppendTo[U, u],
      {j, 2, nt+1}
    ];
    pdeResult[xs, ts, U]]

(* ── Crank-Nicolson (2nd-order in both time and space) ──────────────────── *)
pdeCrankNicolson[alpha2_, xs_, ts_, ic_, bc0_, bcL_, dx_, dt_] :=
  Module[{r = alpha2 dt/(2 dx^2), nx = Length[xs]-2, nt = Length[ts]-1,
          U, u, uNew, rhs, lowerA, diagA, upperA, B},
    lowerA = ConstantArray[-r, nx-1];
    diagA  = ConstantArray[1+2r, nx];
    upperA = ConstantArray[-r, nx-1];
    u = ic /@ xs; U = {u};
    B = DiagonalMatrix[ConstantArray[1-2r, nx], 0] +
        DiagonalMatrix[ConstantArray[r, nx-1], 1] +
        DiagonalMatrix[ConstantArray[r, nx-1], -1];
    Do[
      rhs = B . u[[2;;nx+1]];
      rhs[[1]]  += r (bc0[ts[[j]]] + bc0[ts[[j+1]]]);
      rhs[[-1]] += r (bcL[ts[[j]]] + bcL[ts[[j+1]]]);
      (* Use our Thomas algorithm instead of LinearSolve *)
      uNew = Join[{bc0[ts[[j+1]]]}, laThomasAlgorithm[lowerA, diagA, upperA, rhs], {bcL[ts[[j+1]]]}];
      u = uNew; AppendTo[U, u],
      {j, 1, nt}
    ];
    pdeResult[xs, ts, U]]

(* ── SolveHyperbolicPDE dispatcher ─────────────────────────────────────── *)

SolveHyperbolicPDE[c2_?NumericQ,
    {x_, xL_:0, xR_?NumericQ, nx_Integer},
    {t_, t0_:0, tEnd_?NumericQ, nt_Integer},
    ic_, vic_, {bc0_, bcL_}, opts:OptionsPattern[]] :=
  Module[{method, dx, dt, xs, ts},
    method = OptionValue[Method];
    dx = (xR - xL)/(nx + 1);
    dt = (tEnd - t0)/nt;
    xs = xL + Range[0, nx+1] dx;
    ts = t0 + Range[0, nt] dt;
    pdeWaveExplicit[c2, xs, ts, ic, vic, bc0, bcL, dx, dt]
  ]

(* ── Explicit Finite Difference (stable when CFL = Sqrt[c2]*dt/dx <= 1) ── *)
pdeWaveExplicit[c2_, xs_, ts_, ic_, vic_, bc0_, bcL_, dx_, dt_] :=
  Module[{r2 = c2 dt^2/dx^2, nx = Length[xs]-2, nt = Length[ts]-1, U, u0, u1, u2},
    If[r2 > 1., Message[SolveHyperbolicPDE::cfl, Sqrt[r2]]];
    u0 = ic /@ xs;
    (* First step incorporates initial velocity via Taylor expansion *)
    u1 = Table[
      Which[
        i == 1,      bc0[ts[[2]]],
        i == nx+2,   bcL[ts[[2]]],
        True, u0[[i]] + dt vic[xs[[i]]] + r2/2 (u0[[i+1]] - 2u0[[i]] + u0[[i-1]])
      ], {i, nx+2}];
    U = {u0, u1};
    Do[
      u2 = Table[
        Which[
          i == 1,    bc0[ts[[j]]],
          i == nx+2, bcL[ts[[j]]],
          True, 2u1[[i]] - u0[[i]] + r2 (u1[[i+1]] - 2u1[[i]] + u1[[i-1]])
        ], {i, nx+2}];
      u0 = u1; u1 = u2; AppendTo[U, u2],
      {j, 3, nt+1}
    ];
    pdeResult[xs, ts, U]]
