(* Test_BVP.wl -- Correctness tests for SolveBoundaryValueProblem *)

NKBeginSuite["Boundary Value Problem Correctness", "A", "BVP"]

With[{bvpF = Function[{x, y, dy}, -y]},
  With[{sol = SolveBoundaryValueProblem[bvpF, {x, 0., Pi}, {0., 0.}, 100,
               Method -> "Shooting"]},
    NKTestNear["SolveBoundaryValueProblem/Shooting: y''=-y, y(pi/2) ~ 1",
      sol["Solution"][[51]], 1., "AbsoluteTolerance" -> 0.01];
    NKTest["SolveBoundaryValueProblem/Shooting: boundary y(0) = 0",
      Abs[First[sol["Solution"]]] < 1*^-10]
  ];
  With[{sol = SolveBoundaryValueProblem[bvpF, {x, 0., Pi}, {0., 0.}, 100,
               Method -> "FiniteDifference"]},
    NKTestNear["SolveBoundaryValueProblem/FiniteDiff: y''=-y, y(pi/2) ~ 1",
      sol["Solution"][[51]], 1., "AbsoluteTolerance" -> 0.01]
  ]
]

With[{sol = SolveBoundaryValueProblem[Function[{x,y,dy}, 0.],
             {x, 0., 1.}, {1., 3.}, 50, Method -> "Shooting"]},
  NKTestNear["SolveBoundaryValueProblem/Shooting: y''=0, y(0.5) = 2",
    sol["Solution"][[26]], 2., "AbsoluteTolerance" -> 0.01]
]

With[{sol = SolveBoundaryValueProblem[Function[{x,y,dy},-y], {x,0.,Pi},{0.,0.}]},
  NKTest["SolveBoundaryValueProblem: has Grid and Solution keys",
    KeyExistsQ[sol,"Grid"] && KeyExistsQ[sol,"Solution"]]
]

NKTestFails["SolveBoundaryValueProblem: a >= b returns $Failed",
  SolveBoundaryValueProblem[Function[{x,y,dy},-y], {x, Pi, 0.}, {0.,0.}]]

NKEndSuite[]
