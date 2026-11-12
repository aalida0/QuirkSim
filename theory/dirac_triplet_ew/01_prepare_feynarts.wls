(* QuirkSim experimental NLO QCD+EW renormalization pipeline: step 1 of 3.
   Run as a fresh Mathematica kernel:
      FEYNRULES_HOME=/path/to/feynrules \
      wolframscript -file theory/dirac_triplet_ew/01_prepare_feynarts.wls

   FEYNRULES_HOME must contain FeynRules.m and Models/SM/SM.fr.
   The output is NOT yet an NLO EW UFO; it is a renormalized FeynArts
   input for a subsequent NLOCT calculation.
*)

frHome=Environment["FEYNRULES_HOME"];
If[!StringQ[frHome] || StringLength[frHome]==0,
   Print["ERROR: set FEYNRULES_HOME to a FeynRules installation"];
   Exit[2]
];
frHome=ExpandFileName[frHome];
If[!FileExistsQ[FileNameJoin[{frHome,"FeynRules.m"}]],
   Print["ERROR: no FeynRules.m at ",frHome];Exit[2]
];
scriptDir=DirectoryName[$InputFileName];
modelFile=FileNameJoin[{scriptDir,"TripletDirac.fr"}];
If[!FileExistsQ[modelFile],Print["Missing ",modelFile];Exit[2]];
out=FileNameJoin[{scriptDir,"generated"}];
If[!DirectoryQ[out],CreateDirectory[out]];
AppendTo[$Path,frHome];
SetDirectory[frHome];
Get["FeynRules.m"];
LoadModel["SM.fr"];
Get[modelFile];
FeynmanGauge=True;
Print["Checking new triplet Lagrangian's Hermiticity ..."];
CheckHermiticity[LTripletDirac];
Print["Generating EW-field and mass counterterm Lagrangian (QCDOnly=False) ..."];
LTripletRen = OnShellRenormalization[
    LTripletFull, QCDOnly -> False, FlavorMixing -> False
];
SetDirectory[out];
WriteFeynArtsOutput[
    LTripletRen,
    Output -> "DiracTripletEW_QCDEW_FA",
    GenericFile -> False,
    FlavorExpand -> True
];
Print["FeynArts export attempted at: ",out];
Print["NOTE: model is NOT yet validated at NLO EW. Inspect generated files, UV finite-parts and renormalization scheme."];
