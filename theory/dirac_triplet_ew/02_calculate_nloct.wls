(* QuirkSim experimental NLOCT electroweak one-loop calculation, stage 2.
   MUST use a fresh kernel, not the kernel in which FeynRules was loaded.
   Environments:
        FEYNRULES_HOME=/path/to/feynrules
        FEYNARTS_HOME=/path/to/feynarts
   Wolfram invocation:
        wolframscript -file theory/dirac_triplet_ew/02_calculate_nloct.wls
   This stage requires FeynArts + NLOCT and their compatibility.
   NLOCT default renormalization schemes must be reviewed for complete EW
   renormalization (e.g. tadpoles, Gmu charge renormalization, wave-function
   mixing, photon, IR, complex-mass scheme).
*)
frHome=Environment["FEYNRULES_HOME"];
faHome=Environment["FEYNARTS_HOME"];
If[!StringQ[frHome]||StringLength[frHome]==0||
   !StringQ[faHome]||StringLength[faHome]==0,
   Print["ERROR: FEYNRULES_HOME and FEYNARTS_HOME required"];Exit[2]];
frHome=ExpandFileName[frHome];
faHome=ExpandFileName[faHome];
scriptDir=DirectoryName[ExpandFileName[$InputFileName]];
out=FileNameJoin[{scriptDir,"generated"}];
If[!DirectoryQ[out],Print["ERROR: Run stage 01 first"];Exit[2]];
AppendTo[$Path,faHome];
SetDirectory[faHome];
Get["FeynArts.m"];
AppendTo[$Path,frHome];
SetDirectory[frHome];
Get["NLOCT.m"];

SetDirectory[out];
Print["Starting NLOCT *full* QCD+EW UV/R2 construction. No QCDOnly shortcut ..."];
WriteCT[
    "DiracTripletEW_QCDEW_FA",
    "Lorentz",
    Output -> "DiracTripletEW_QCDEW_CT",
    QCDOnly -> False,
    Exclude4ScalarsCT -> False,
    Assumptions -> {MQ > 0, MZ > 0, MW > 0}
];
Print["NLOCT output attempted in: ",out];
Print["WARNING: Counterterms are unvalidated and MUST NOT yet be used for NLO-EW predictions."];
