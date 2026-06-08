(* Robustness_ErrorMessages.wl -- Verify correct error messages are emitted *)

NKBeginSuite["Error Message Verification", "C", "ErrorMessages"]

NKTestMessage["FindEquationRoot::bracket fires for no sign change",
  FindEquationRoot[Function[x,x^2+1],{-2.,2.}],
  "FindEquationRoot::bracket"]
NKTestMessage["SolveLinearEquationSystem::badmethod fires",
  SolveLinearEquationSystem[{{1.,2.},{3.,4.}},{1.,2.},Method->"Unknown"],
  "SolveLinearEquationSystem::badmethod"]
NKTestMessage["FindMinimum1D::badmethod fires",
  FindMinimum1D[#^2&,{-1.,1.},Method->"Unknown"],
  "FindMinimum1D::badmethod"]
NKTest["FindMinimum1D: badmethod falls back and returns valid result",
  AssociationQ[Quiet[FindMinimum1D[#^2&,{-1.,1.},Method->"Unknown"]]]]
NKTestMessage["SolveParabolicPDE::stiff fires when r>0.5",
  SolveParabolicPDE[1.,{x,0.,1.,5},{t,0.,0.1,1},
    Function[x,Sin[Pi x]],{Function[t,0.],Function[t,0.]},Method->"ForwardDifference"],
  "SolveParabolicPDE::stiff"]
NKTestMessage["NumericalQuadrature::badmethod fires",
  NumericalQuadrature[Sin,{x,0.,Pi},Method->"BadMethod"],
  "NumericalQuadrature::badmethod"]
NKTest["FindMinimum1D::badmethod can be silenced with Off[]",
  Module[{fired=False},
    Off[FindMinimum1D::badmethod];
    Check[FindMinimum1D[#^2&,{-1.,1.},Method->"X"],fired=True,FindMinimum1D::badmethod];
    On[FindMinimum1D::badmethod];
    !fired]]

NKEndSuite[]
