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
qsGeneratedRoot=FileNameJoin[{scriptDir,"generated"}];
If[!DirectoryQ[qsGeneratedRoot],Print["ERROR: Run stage 01 first"];Exit[2]];
AppendTo[$Path,faHome];
SetDirectory[faHome];
Get["FeynArts.m"];
AppendTo[$Path,frHome];
SetDirectory[frHome];
Get["NLOCT.m"];

SetDirectory[qsGeneratedRoot];
If[!FileExistsQ[FileNameJoin[{qsGeneratedRoot,"DiracTripletEW_QCDEW_FA",
                                  "DiracTripletEW_QCDEW_FA.mod"}]],
   Print["ERROR: stage 01 must produce a FeynArts .mod model file in generated/DiracTripletEW_QCDEW_FA"];Exit[2]];
If[!FileExistsQ[FileNameJoin[{qsGeneratedRoot,"DiracTripletEW_QCDEW_FA","DiracTripletEW_QCDEW_FA.gen"}]],
   Print["ERROR: stage 01 must produce a model-specific .gen file for EW CTs; Lorentz.gen is insufficient."];
   Exit[2]
];
qsAuditFile=FileNameJoin[{qsGeneratedRoot,"stage1_log_audit.json"}];
qsModelStem=FileNameJoin[{qsGeneratedRoot,"DiracTripletEW_QCDEW_FA",
                          "DiracTripletEW_QCDEW_FA"}];
If[!FileExistsQ[qsAuditFile],
   Print["ERROR: stage 01 has no successful log-audit manifest. ",
         "Use scripts/run_triplet_ew_stage1.sh and fix errors first."];
   Exit[3]
];
qsAudit=Import[qsAuditFile,"RawJSON"];
If[!TrueQ[Lookup[qsAudit,"accepted_for_nloct_preflight",False]],
   Print["ERROR: stage 01 FeynArts model failed fail-closed log audit; ",
         "NLOCT would consume malformed vertices."];
   Exit[3]
];
qsHashes=Lookup[qsAudit,"model_sha256",<||>];
If[!AssociationQ[qsHashes] ||
   !And@@Table[
      FileHash[qsModelStem<>ext,"SHA256","HexString"]===Lookup[qsHashes,ext,""],
      {ext,{".mod",".pars",".gen"}}],
   Print["ERROR: FeynArts .mod/.pars changed since the stage-1 log audit."];
   Exit[3]
];
Print["FR_STAGE1_LOG_AUDIT_PASSED_ONLY: no known rejected-vertex patterns."];
Print["CAUTION: physics-level NLO EW counterterms are not validated."];
Print["Starting NLOCT *full* QCD+EW UV/R2 construction. No QCDOnly shortcut ..."];
qsNLOOutputFile=FileNameJoin[{qsGeneratedRoot,"DiracTripletEW_QCDEW_CT.nlo"}];
If[FileExistsQ[qsNLOOutputFile],
   Print["ERROR: existing NLOCT output found: ",qsNLOOutputFile,
         ". Archive it explicitly before starting a new integration."];
   Exit[3]
];
qsNLOCTResult=Check[
    WriteCT[
        "DiracTripletEW_QCDEW_FA/DiracTripletEW_QCDEW_FA",
        "DiracTripletEW_QCDEW_FA/DiracTripletEW_QCDEW_FA",
        Output -> "DiracTripletEW_QCDEW_CT",
        QCDOnly -> False,
        Exclude4ScalarsCT -> False,
        Assumptions -> {MQ > 0, MZ > 0, MW > 0}
    ],
    $Failed,
    {Part::partd, Part::partw, Set::partd, Set::partw,
     First::nofirst, Last::nolast, Function::slotn,
     StringJoin::string, Power::infy, Infinity::indet}
];
If[qsNLOCTResult === $Failed,
   Print["NLOCT_STAGE2_REJECTED_INVALID_MATHEMATICA_EXPRESSION: ",
         "one-loop counterterm generation emitted structural errors."];
   Exit[5]
];
If[!FileExistsQ[qsNLOOutputFile] ||
   FileByteCount[qsNLOOutputFile] < 200,
   Print["ERROR: NLOCT did not produce a sufficiently populated .nlo file at ",
         qsNLOOutputFile];
   Exit[6]
];
Print["NLOCT_STAGE2_FILE_WRITTEN_UNVALIDATED: ",qsNLOOutputFile,
      " (",FileByteCount[qsNLOOutputFile]," bytes)"];
Print["NLOCT_STAGE2_REQUIRES_LOG_AND_CONTENT_AUDIT"];
Print["NLOCT output attempted in: ",qsGeneratedRoot];
Print["WARNING: Counterterms are unvalidated and MUST NOT yet be used for NLO-EW predictions."];
