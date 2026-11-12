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
scriptDir=DirectoryName[ExpandFileName[$InputFileName]];
modelFile=FileNameJoin[{scriptDir,"TripletDirac.fr"}];
If[!FileExistsQ[modelFile],Print["Missing ",modelFile];Exit[2]];
qsGeneratedDir=FileNameJoin[{scriptDir,"generated"}];
If[!DirectoryQ[qsGeneratedDir],CreateDirectory[qsGeneratedDir]];
AppendTo[$Path,frHome];
SetDirectory[frHome];
Get["FeynRules.m"];
LoadModel[FileNameJoin[{frHome,"Models","SM","SM.fr"}], modelFile];
Print["FR_TRIPLET_MERGED_MODEL_LOADED"];
FeynmanGauge=True;
Print["Checking new triplet Lagrangian's Hermiticity ..."];
CheckHermiticity[LTripletDirac];
Print["Generating EW-field and mass counterterm Lagrangian (QCDOnly=False) ..."];
LTripletRen = OnShellRenormalization[
    LTripletFull, QCDOnly -> False, FlavorMixing -> False
];
SetDirectory[qsGeneratedDir];
WriteFeynArtsOutput[
    LTripletRen,
    Output -> "DiracTripletEW_QCDEW_FA",
    GenericFile -> False,
    FlavorExpand -> True
];
qsModelMod=FileNameJoin[{qsGeneratedDir,"DiracTripletEW_QCDEW_FA","DiracTripletEW_QCDEW_FA.mod"}];
If[!FileExistsQ[qsModelMod],
   Print["ERROR: FeynArts .mod missing after WriteFeynArtsOutput: ",qsModelMod];
   Print["Stage 01 did not finish. Do not run stage 02."];
   Exit[3]
];
Print["FR_STAGE1_FEYNARTS_MODEL_OK: ", qsModelMod];
Print["NOTE: model is NOT yet validated at NLO EW. Inspect generated files, UV finite-parts and renormalization scheme."];
