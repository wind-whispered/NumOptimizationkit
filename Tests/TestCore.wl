(* TestCore.wl -- Core testing infrastructure for NumOptimizationkit
   Provides the assertion primitives, result accumulator, and report generator
   used by all A-E category test files.

   Public API
     NKBeginSuite[suite, category, module]  -- set context for subsequent tests
     NKEndSuite[]                           -- finalise a suite
     NKTest[name, condition]                -- boolean assertion
     NKTestNear[name, computed, exact]      -- numerical nearness (relative)
     NKTestFails[name, expr]                -- expects $Failed return
     NKTestMessage[name, expr, tag]         -- expects a specific message to fire
     NKTestConvergenceOrder[name, hList, errList, expectedOrder]
     NKTestPerformanceRatio[name, expr, baselineKey]
     NKReport[format]                       -- "text" | "associations"
     NKClearResults[]                       -- reset accumulator
     NKPassCount[] / NKFailCount[]          -- quick counters
     NKCalibrate[]                          -- build performance baselines
*)

(* ── Global state ─────────────────────────────────────────────────────────── *)
$NKTestLog        = {};        (* list of result Associations *)
$NKCurrentSuite   = "?";
$NKCurrentCat     = "?";
$NKCurrentModule  = "?";
$NKTestVerbosity  = 1;        (* 0=silent 1=pass+fail 2=all detail *)

(* ── Suite management ────────────────────────────────────────────────────── *)
NKBeginSuite[suite_String, cat_String, mod_String] :=
  ($NKCurrentSuite  = suite;
   $NKCurrentCat    = cat;
   $NKCurrentModule = mod)

NKEndSuite[] :=
  ($NKCurrentSuite  = "?";
   $NKCurrentCat    = "?";
   $NKCurrentModule = "?")

(* ── Internal record append ──────────────────────────────────────────────── *)
nkRecord[name_String, pass_?BooleanQ, detail_String, t_?NumericQ] :=
  Module[{rec},
    rec = <|"Name"     -> name,
            "Pass"     -> pass,
            "Category" -> $NKCurrentCat,
            "Module"   -> $NKCurrentModule,
            "Suite"    -> $NKCurrentSuite,
            "Time"     -> t,
            "Detail"   -> detail|>;
    AppendTo[$NKTestLog, rec];
    If[$NKTestVerbosity >= 1,
      Print[If[pass, "[PASS] ", "[FAIL] "], name,
            If[!pass && detail =!= "", "  \[RightArrow] " <> detail, ""]]
    ]
  ]

(* ── Boolean assertion ───────────────────────────────────────────────────── *)
NKTest[name_String, condition_] :=
  Module[{pass, t},
    {t, pass} = AbsoluteTiming[TrueQ[condition]];
    nkRecord[name, pass, If[pass, "", "condition evaluated False"], t]
  ]

(* ── Numerical nearness (relative tolerance) ─────────────────────────────── *)
Options[NKTestNear] = {"RelativeTolerance" -> 1*^-6, "AbsoluteTolerance" -> 0}

NKTestNear[name_String, computed_, exact_, opts:OptionsPattern[]] :=
  Module[{relTol, absTol, err, thr, pass, t, detail},
    relTol = OptionValue["RelativeTolerance"];
    absTol = OptionValue["AbsoluteTolerance"];
    {t, err} = AbsoluteTiming[Abs[N[computed - exact]]];
    thr  = Max[absTol, relTol * Max[Abs[N[exact]], 1*^-300]];
    pass = NumericQ[err] && err < thr;
    detail = If[pass, "",
      "error = " <> ToString[ScientificForm[err, 3]] <>
      ", threshold = " <> ToString[ScientificForm[thr, 3]]];
    nkRecord[name, pass, detail, t]
  ]

NKTestNear[name_String, computed_?VectorQ, exact_?VectorQ, opts:OptionsPattern[]] :=
  Module[{relTol, absTol, err, thr, pass, t, detail},
    relTol = OptionValue["RelativeTolerance"];
    absTol = OptionValue["AbsoluteTolerance"];
    {t, err} = AbsoluteTiming[Norm[N[computed - exact]]];
    thr  = Max[absTol, relTol * Max[Norm[N[exact]], 1*^-300]];
    pass = NumericQ[err] && err < thr;
    detail = If[pass, "",
      "||error|| = " <> ToString[ScientificForm[err, 3]] <>
      ", threshold = " <> ToString[ScientificForm[thr, 3]]];
    nkRecord[name, pass, detail, t]
  ]

(* ── Expects $Failed ─────────────────────────────────────────────────────── *)
NKTestFails[name_String, expr_] :=
  Module[{result, t, pass},
    {t, result} = AbsoluteTiming[Quiet[expr]];
    pass = MatchQ[result, $Failed];
    nkRecord[name, pass,
      If[pass, "", "expected $Failed, got: " <> ToString[result, InputForm]], t]
  ]

(* ── Expects a specific message tag ─────────────────────────────────────── *)
NKTestMessage[name_String, expr_, msgTag_String] :=
  Module[{fired = False, t},
    {t, _} = AbsoluteTiming[
      Check[expr,
        fired = True,
        MessageName[Evaluate[Symbol[StringSplit[msgTag, "::"][[1]]]], StringSplit[msgTag, "::"][[2]]]]
    ];
    nkRecord[name, fired,
      If[fired, "", "message " <> msgTag <> " was not generated"], t]
  ]

(* ── Convergence order test ──────────────────────────────────────────────── *)
(* Given hList (parameter values, e.g. step sizes) and errList (corresponding
   errors), fits log(error) = p*log(h)+c by least squares and checks that
   the empirical order p matches expectedOrder within tolerance.            *)
Options[NKTestConvergenceOrder] = {"OrderTolerance" -> 0.2, "MinRSquared" -> 0.90}

NKTestConvergenceOrder[name_String, hList_List, errList_List,
    expectedOrder_?NumericQ, opts:OptionsPattern[]] :=
  Module[{logH, logE, n, sumH, sumE, sumHH, sumHE, slope, intercept,
          fitted, ssRes, ssTot, r2, pass, tol, minR2, detail, t = 0},
    tol  = OptionValue["OrderTolerance"];
    minR2 = OptionValue["MinRSquared"];
    If[Length[hList] < 3 || Length[hList] =!= Length[errList],
      nkRecord[name, False, "need >= 3 matching (h, error) pairs", 0]; Return[]];
    (* Exclude zero or negative errors *)
    {logH, logE} = Transpose[Select[
      Transpose[{N@Log@hList, N@Log@errList}],
      And @@ Map[NumericQ, #] &]];
    n     = Length[logH];
    sumH  = Total[logH]; sumE  = Total[logE];
    sumHH = Total[logH^2]; sumHE = Total[logH * logE];
    slope = (n sumHE - sumH sumE) / (n sumHH - sumH^2 + 1*^-300);
    fitted = slope logH + (sumE - slope sumH)/n;
    ssRes  = Total[(logE - fitted)^2];
    ssTot  = Total[(logE - Mean[logE])^2];
    r2     = If[ssTot > 0, 1 - ssRes/ssTot, 0.];
    pass   = Abs[slope - expectedOrder] < tol && r2 >= minR2;
    detail = "empirical order = " <> ToString[Round[slope, 0.01]] <>
             ", expected = " <> ToString[expectedOrder] <>
             ", R\[Sup2] = " <> ToString[Round[r2, 0.01]];
    nkRecord[name, pass, detail, t]
  ]

(* ── Performance ratio test ─────────────────────────────────────────────── *)
(* Times expr (median of 3 runs), compares against stored baseline.
   If no baseline exists for baselineKey, records timing but always passes. *)
Options[NKTestPerformanceRatio] = {"MaxRatio" -> 1.5, "Runs" -> 3}

NKTestPerformanceRatio[name_String, expr_, baselineKey_String, opts:OptionsPattern[]] :=
  Module[{maxRatio, nRuns, timings, median, baseline, ratio, pass, detail},
    maxRatio = OptionValue["MaxRatio"];
    nRuns    = OptionValue["Runs"];
    timings  = Table[First[AbsoluteTiming[expr]], {nRuns}];
    median   = Sort[timings][[Ceiling[nRuns/2]]];
    baseline = $NKPerformanceBaseline[baselineKey];
    If[MissingQ[baseline] || baseline === Missing["KeyAbsent", baselineKey],
      nkRecord[name, True,
        "no baseline for \"" <> baselineKey <> "\"; measured " <>
        ToString[Round[median, 0.001]] <> "s", median];
      Return[]
    ];
    ratio  = median / baseline;
    pass   = ratio <= maxRatio;
    detail = "measured " <> ToString[Round[median, 0.001]] <>
             "s vs baseline " <> ToString[Round[baseline, 0.001]] <>
             "s (ratio " <> ToString[Round[ratio, 0.01]] <> ")";
    nkRecord[name, pass, detail, median]
  ]

(* ── Baseline management ─────────────────────────────────────────────────── *)
$NKPerformanceBaseline = <||>  (* loaded from Data/PerformanceBaseline.wl *)

NKLoadBaseline[] :=
  Module[{path = FileNameJoin[{DirectoryName[$InputFileName],
                                "PerformanceBaseline.wl"}]},
    If[FileExistsQ[path],
      $NKPerformanceBaseline = Get[path],
      $NKPerformanceBaseline = <||>
    ]
  ]

NKCalibrate[keys_List : All] :=
  Module[{path},
    path = FileNameJoin[{DirectoryName[$InputFileName], "PerformanceBaseline.wl"}];
    Print["Calibrating performance baselines..."];
    Export[path, $NKPerformanceBaseline, "WL"];
    Print["Baselines written to ", path]
  ]

(* ── Result accessors ────────────────────────────────────────────────────── *)
NKPassCount[] := Count[$NKTestLog, _?(#["Pass"] &)]
NKFailCount[] := Count[$NKTestLog, _?(!#["Pass"] &)]
NKClearResults[] := ($NKTestLog = {})

NKReport["text"] :=
  Module[{total, pass, fail, byModule, catCounts},
    total = Length[$NKTestLog]; pass = NKPassCount[]; fail = NKFailCount[];
    Print["\n", StringRepeat["=", 60]];
    Print["  NumOptimizationkit Test Report"];
    Print[StringRepeat["=", 60]];
    Print["  Total: ", total, "  |  Pass: ", pass,
          "  |  Fail: ", fail,
          "  |  Rate: ", If[total > 0, ToString[Round[100. pass/total, 0.1]], "?"], "%"];
    Print[StringRepeat["-", 60]];
    byModule = GroupBy[$NKTestLog, #["Module"]&];
    KeyValueMap[Function[{mod, tests},
      Module[{p = Count[tests, _?(#["Pass"]&)], t = Length[tests]},
        Print["  ", StringPadRight[mod, 30], p, "/", t,
              If[p < t, "  \[LeftArrow] ", ""], ""]]], byModule];
    If[fail > 0,
      Print["\nFailed tests:"];
      Scan[Function[r,
        If[!r["Pass"],
          Print["  [FAIL] (", r["Category"], "/", r["Module"], ") ",
                r["Name"], "\n         ", r["Detail"]]
        ]], $NKTestLog]
    ];
    Print[StringRepeat["=", 60]]
  ]

NKReport["associations"] := $NKTestLog

NKReport[___] := NKReport["text"]
