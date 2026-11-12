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
qsSourceLagrangian=LTripletDirac;
If[Head[qsSourceLagrangian] === Symbol || qsSourceLagrangian === 0,
   Print["ERROR: LTripletDirac is undefined or zero. Do NOT renormalize."];
   Exit[5]
];
qsTreeVertices=FeynmanRules[qsSourceLagrangian,
                          ScreenOutput->False,MaxParticles->3];
Print["FR_STAGE1_TRIPLET_TREE_VERTEX_COUNT: ",Length[qsTreeVertices]];
If[!ListQ[qsTreeVertices]||Length[qsTreeVertices]<8,
   Print["ERROR: fewer than 8 tree vertices; refusing false-positive hermiticity test."];
   Exit[5]
];
CheckHermiticity[LTripletDirac];
Print["Generating EW-field and mass counterterm Lagrangian (QCDOnly=False) ..."];
LTripletRen = Check[
    OnShellRenormalization[
        LTripletFull, QCDOnly -> False, FlavorMixing -> False
    ],
    $Failed,
    {Part::partd, Part::partw, Set::partd, Set::partw,
     First::nofirst, Function::slotn, StringJoin::string}
];
If[LTripletRen === $Failed,
   Print["FR_STAGE1_RENORM_INVALID: an unresolved Mathematica error occurred ",
         "inside OnShellRenormalization. FeynArts export was NOT attempted."];
   Print["Do NOT run NLOCT until the EW renormalization scheme is repaired."];
   Exit[4]
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
Print["FR_STAGE1_FEYNARTS_FILES_WRITTEN_UNVALIDATED: ", qsModelMod];
Print["STAGE1_REQUIRES_LOG_AUDIT: files written; an external fail-closed log check is mandatory before NLOCT."];
Print["NOTE: model is NOT yet validated at NLO EW. Inspect generated files, UV finite-parts and renormalization scheme."];
