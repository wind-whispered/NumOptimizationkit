(* TestRunner.wl -- Test scheduler and dispatcher for NumOptimizationkit

   All test files live flat in the Tests/ directory alongside this file.

   Naming convention (no single-letter prefixes to avoid conflicts):
     Test_<Module>.wl           Correctness tests  (Category A)
     Convergence_<Module>.wl    Convergence order  (Category B)
     Robustness_<Type>.wl       Robustness         (Category C)
     Accuracy_<Module>.wl       Accuracy benchmarks(Category D)
     Performance_<Module>.wl    Performance timing  (Category E)

   Scheduling design (dependency-aware)
   ─────────────────────────────────────
     Category A  runs first (foundation for all other categories)
     Category C  runs alongside A  (independent of A's results)
     Category E  runs alongside A  (independent of A's results)
     Category B  runs after A, only if A pass rate >= 50%
     Category D  runs after A, only if A pass rate >= 50%

   Public API
     RunTests[opts]              master entry point
     RunTestCategory[cat, opts]  run one category (A/B/C/D/E)
     RunTestModule[mod, opts]    run tests for one module name
     RunQuickTests[]             fast smoke test (Category A only, silent)
     NKTestSummary[]             per-category pass/fail table
*)

(* ── Test registry (flat file paths) ────────────────────────────────────── *)
(* IMPORTANT: $NKTestsBaseDir is set EAGERLY at load time (when $InputFileName
   is valid). $NKTestRegistry uses :=/:= memoisation but reads $NKTestsBaseDir,
   not $InputFileName, so the path is correct even when called later.          *)
$NKTestsBaseDir = DirectoryName[$InputFileName];   (* evaluated immediately on Get *)

$NKTestRegistry := $NKTestRegistry =
  Module[{base = $NKTestsBaseDir},
    {
      (* ── Category A: Correctness tests ── *)
      <|"File" -> FileNameJoin[{base, "Test_Optimization.wl"}],
        "Category" -> "A", "Module" -> "Optimization",      "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_RootFinding.wl"}],
        "Category" -> "A", "Module" -> "RootFinding",       "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_NonlinearSystems.wl"}],
        "Category" -> "A", "Module" -> "NonlinearSystems",  "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_Quadrature.wl"}],
        "Category" -> "A", "Module" -> "Quadrature",        "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_ODE.wl"}],
        "Category" -> "A", "Module" -> "ODE",               "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_BVP.wl"}],
        "Category" -> "A", "Module" -> "BVP",               "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_PDE.wl"}],
        "Category" -> "A", "Module" -> "PDE",               "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_LinearAlgebra.wl"}],
        "Category" -> "A", "Module" -> "LinearAlgebra",     "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_Eigenanalysis.wl"}],
        "Category" -> "A", "Module" -> "Eigenanalysis",     "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_Interpolation.wl"}],
        "Category" -> "A", "Module" -> "Interpolation",     "Tags" -> {"correctness"}|>,
      <|"File" -> FileNameJoin[{base, "Test_Approximation.wl"}],
        "Category" -> "A", "Module" -> "Approximation",     "Tags" -> {"correctness"}|>,
      (* ── Category B: Convergence order ── *)
      <|"File" -> FileNameJoin[{base, "Convergence_ODE.wl"}],
        "Category" -> "B", "Module" -> "ODE_Orders",        "Tags" -> {"convergence"}|>,
      <|"File" -> FileNameJoin[{base, "Convergence_Quadrature.wl"}],
        "Category" -> "B", "Module" -> "Quadrature_Orders", "Tags" -> {"convergence"}|>,
      <|"File" -> FileNameJoin[{base, "Convergence_RootFinding.wl"}],
        "Category" -> "B", "Module" -> "RootFinding_Orders","Tags" -> {"convergence"}|>,
      (* ── Category C: Robustness ── *)
      <|"File" -> FileNameJoin[{base, "Robustness_EdgeCases.wl"}],
        "Category" -> "C", "Module" -> "EdgeCases",         "Tags" -> {"robustness"}|>,
      <|"File" -> FileNameJoin[{base, "Robustness_ErrorMessages.wl"}],
        "Category" -> "C", "Module" -> "ErrorMessages",     "Tags" -> {"robustness"}|>,
      (* ── Category D: Accuracy benchmarks ── *)
      <|"File" -> FileNameJoin[{base, "Accuracy_ODE.wl"}],
        "Category" -> "D", "Module" -> "ODE_Reference",     "Tags" -> {"accuracy"}|>,
      <|"File" -> FileNameJoin[{base, "Accuracy_Quadrature.wl"}],
        "Category" -> "D", "Module" -> "Quad_Reference",    "Tags" -> {"accuracy"}|>,
      <|"File" -> FileNameJoin[{base, "Accuracy_Optimization.wl"}],
        "Category" -> "D", "Module" -> "Opt_Reference",     "Tags" -> {"accuracy"}|>,
      (* ── Category E: Performance regression ── *)
      <|"File" -> FileNameJoin[{base, "Performance_ODE.wl"}],
        "Category" -> "E", "Module" -> "ODE_Timing",        "Tags" -> {"performance"}|>,
      <|"File" -> FileNameJoin[{base, "Performance_LinearAlgebra.wl"}],
        "Category" -> "E", "Module" -> "LA_Timing",         "Tags" -> {"performance"}|>,
      <|"File" -> FileNameJoin[{base, "Performance_Quadrature.wl"}],
        "Category" -> "E", "Module" -> "Quad_Timing",       "Tags" -> {"performance"}|>
    }
  ]

(* ── Internal helpers ────────────────────────────────────────────────────── *)

nkRunFile[entry_Association] :=
  Module[{path = entry["File"], t, p0, f0, p1, f1},
    If[!FileExistsQ[path],
      Print["[SKIP] file not found: ", FileNameTake[path]]; Return[{0, 0, 0}]];
    p0 = NKPassCount[]; f0 = NKFailCount[];
    {t, _} = AbsoluteTiming[Get[path]];
    p1 = NKPassCount[]; f1 = NKFailCount[];
    {p1 - p0, f1 - f0, t}
  ]

nkRunFileList[entries_List, parallel_] :=
  Module[{results},
    results = If[TrueQ[parallel] && $KernelCount > 1,
      ParallelMap[nkRunFile, entries, DistributedContexts -> Automatic],
      Map[nkRunFile, entries]
    ];
    If[results === {}, {0, 0, 0}, Total /@ Transpose[results]]
  ]

nkCategoryPassRate[cat_String] :=
  Module[{tests = Select[$NKTestLog, #["Category"] === cat &]},
    If[Length[tests] == 0, 1.,
      N[Count[tests, _?(#["Pass"] &)] / Length[tests]]]
  ]

(* ── Main options ────────────────────────────────────────────────────────── *)
Options[RunTests] = {
  "Categories"           -> All,
  "Modules"              -> All,
  "Parallel"             -> False,
  "Verbosity"            -> 1,
  "SkipBDonLowA"         -> True,
  "LowPassRateThreshold" -> 0.50,
  "ReportFormat"         -> "text"
}

(* ── Master entry point ──────────────────────────────────────────────────── *)
RunTests[opts:OptionsPattern[]] :=
  Module[{cats, mods, par, verb, skipBD, thresh, fmt,
          aE, bE, cE, dE, eE, aRate, skip, startTime, filter},

    cats   = OptionValue["Categories"];
    mods   = OptionValue["Modules"];
    par    = TrueQ[OptionValue["Parallel"]];
    verb   = OptionValue["Verbosity"];
    skipBD = TrueQ[OptionValue["SkipBDonLowA"]];
    thresh = OptionValue["LowPassRateThreshold"];
    fmt    = OptionValue["ReportFormat"];

    $NKTestVerbosity = verb;
    NKClearResults[];
    NKLoadBaseline[];

    filter = Function[e,
      (cats === All || MemberQ[cats, e["Category"]]) &&
      (mods === All || MemberQ[mods, e["Module"]])
    ];
    aE = Select[$NKTestRegistry, filter[#] && #["Category"] === "A" &];
    bE = Select[$NKTestRegistry, filter[#] && #["Category"] === "B" &];
    cE = Select[$NKTestRegistry, filter[#] && #["Category"] === "C" &];
    dE = Select[$NKTestRegistry, filter[#] && #["Category"] === "D" &];
    eE = Select[$NKTestRegistry, filter[#] && #["Category"] === "E" &];

    startTime = AbsoluteTime[];

    (* ── Phase 1: A — correctness foundation ── *)
    If[verb >= 1, Print["\n>> [A] Correctness Tests"]];
    nkRunFileList[aE, par];
    aRate = nkCategoryPassRate["A"];
    If[verb >= 1,
      Print["   Pass rate: ", Round[100 aRate, 0.1], "% (",
            NKPassCount[], " / ", NKPassCount[] + NKFailCount[], ")"]];

    (* ── Phase 2: C + E — independent, run regardless of A ── *)
    If[Length[cE] > 0,
      If[verb >= 1, Print["\n>> [C] Robustness Tests"]];
      nkRunFileList[cE, par]];
    If[Length[eE] > 0,
      If[verb >= 1, Print["\n>> [E] Performance Tests"]];
      nkRunFileList[eE, par]];

    (* ── Phase 3: B + D — conditional on A pass rate ── *)
    skip = skipBD && aRate < thresh;
    If[skip,
      If[verb >= 1,
        Print["\n[WARNING] A pass rate (", Round[100 aRate, 0.1],
              "%) < threshold (", Round[100 thresh, 0.1],
              "%): skipping B (Convergence) and D (Accuracy)."]],
      If[Length[bE] > 0,
        If[verb >= 1, Print["\n>> [B] Convergence Order Tests"]];
        nkRunFileList[bE, par]];
      If[Length[dE] > 0,
        If[verb >= 1, Print["\n>> [D] Accuracy Benchmarks"]];
        nkRunFileList[dE, par]]
    ];

    If[verb >= 1,
      Print["\n   Total elapsed: ", Round[AbsoluteTime[] - startTime, 0.1], " s"]];
    NKReport[fmt]
  ]

(* ── Convenience runners ─────────────────────────────────────────────────── *)

RunTestCategory[cat_String, opts:OptionsPattern[RunTests]] :=
  RunTests["Categories" -> {cat}, opts]

RunTestModule[mod_String, opts:OptionsPattern[RunTests]] :=
  RunTests["Modules" -> {mod}, opts]

RunQuickTests[] :=
  Module[{p, f},
    RunTests["Categories" -> {"A"}, "Verbosity" -> 0,
             "ReportFormat" -> "associations"];
    p = NKPassCount[]; f = NKFailCount[];
    Print["Quick check: ", p, " passed / ", f, " failed  (",
          Round[100. p / Max[p + f, 1], 0.1], "%)"];
    f == 0
  ]

(* ── Per-category summary table ─────────────────────────────────────────── *)
NKTestSummary[] :=
  Module[{byCat},
    byCat = GroupBy[$NKTestLog, #["Category"] &];
    Print["Cat  Pass  Fail  Rate"];
    Print[StringRepeat["-", 28]];
    KeyValueMap[Function[{cat, tests},
      Module[{p = Count[tests, _?(#["Pass"] &)], t = Length[tests]},
        Print[StringPadRight[cat, 5],
              StringPadLeft[ToString[p], 5],
              StringPadLeft[ToString[t - p], 6], "  ",
              Round[100. p / Max[t, 1], 0.1], "%"]
      ]], KeySort[byCat]];
    Print[StringRepeat["-", 28]];
    Print[StringPadRight["ALL", 5],
          StringPadLeft[ToString[NKPassCount[]], 5],
          StringPadLeft[ToString[NKFailCount[]], 6], "  ",
          Round[100. NKPassCount[] / Max[Length[$NKTestLog], 1], 0.1], "%"]
  ]

(* ── Auto-initialise on load ────────────────────────────────────────────── *)
(* All paths below use $NKTestsBaseDir (captured eagerly above) so they
   remain correct whether the function is called at load time or later.     *)

(* Package is one directory above Tests/ *)
$NKPackageDir = ParentDirectory[$NKTestsBaseDir];

nkEnsurePackageLoaded[] :=
  If[!MemberQ[$Packages, "NumOptimizationkit`"],
    Quiet@Get[FileNameJoin[{$NKPackageDir, "NumOptimizationkit.wl"}]]
  ]

(* Load TestCore.wl from the captured Tests/ base directory *)
Get[FileNameJoin[{$NKTestsBaseDir, "TestCore.wl"}]];

(* Set the Data/ path eagerly so TestCore's NKLoadBaseline / NKCalibrate work *)
$NKTestDataDir = FileNameJoin[{$NKTestsBaseDir, "Data"}];

nkEnsurePackageLoaded[];
NKLoadBaseline[];
Print["NumOptimizationkit Test Runner ready.  Use RunTests[] to begin."];
