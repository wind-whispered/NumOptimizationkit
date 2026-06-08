(* ::Package:: *)

(* Kernel/init.m -- Paclet bootstrap: loads NumOptimizationkit once,
   then reports basic runtime diagnostics (version / OS / front-end). *)

(* ---- Load guard: skip if something already pulled in the package ---- *)
If[! MemberQ[$Packages, "NumOptimizationkit`"],
  Get[FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]], "NumOptimizationkit.wl"}]]
];

(* ---- Version compatibility notice (soft warning, not a hard stop) --- *)
(* Association literals, SelectFirst, Nothing etc. used throughout the
   package only require WL >= 10.1, but 12.0 is the floor this package
   has actually been built and exercised against -- warn below that.    *)
If[$VersionNumber < 12.0,
  Print[Style[
    "NumOptimizationkit: running on Wolfram Language " <> ToString[$VersionNumber] <>
    " -- this package is developed and tested on 12.0 or later; some \
functions may not behave as expected on older kernels.",
    Orange]]
];

(* ---- Self-description, handy for bug reports ------------------------ *)
NumOptimizationkit`PackageInformation::usage =
  "PackageInformation[] returns an Association describing the package \
version and the runtime it is loaded into (Wolfram Language version, \
operating system, system ID, kernel count) -- include its output when \
reporting issues.";
NumOptimizationkit`PackageInformation[] := <|
  "Name"            -> "NumOptimizationkit",
  "Version"         -> "1.0.0",
  "WolframVersion"  -> $VersionNumber,
  "OperatingSystem" -> $OperatingSystem,   (* "Windows" | "MacOSX" | "Unix" *)
  "SystemID"        -> $SystemID,          (* e.g. "Windows-x86-64" *)
  "KernelCount"     -> $KernelCount
|>;

(* ---- Front-end-aware "ready" banner ---------------------------------- *)
With[{msg = "NumOptimizationkit 1.0.0 ready (" <> $OperatingSystem <> ", WL " <>
            ToString[$VersionNumber] <> ") -- see PackageInformation[]."},
  If[$FrontEnd =!= $Failed,
    PrintTemporary[Style[msg, Darker[Green]]],
    Print[msg]
  ]
];
