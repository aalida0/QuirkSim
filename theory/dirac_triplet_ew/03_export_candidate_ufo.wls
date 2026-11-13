(* QuirkSim experimental NLO QCD+EW UFO export, stage 3.
   Run in a third, fresh Mathematica kernel after stage 02:
      FEYNRULES_HOME=/path/to/feynrules \
      wolframscript -file theory/dirac_triplet_ew/03_export_candidate_ufo.wls
   Candidate UFO MUST pass independent one-loop QED UV/IR, gauge, SM and
   NLO-QCD regression checks before being eligible for physics runs.
*)
frHome=Environment["FEYNRULES_HOME"];
If[!StringQ[frHome]||StringLength[frHome]==0,
   Print["ERROR: FEYNRULES_HOME required"];Exit[2]];
frHome=ExpandFileName[frHome];
scriptDir=DirectoryName[ExpandFileName[$InputFileName]];
qsStage3GeneratedRoot=FileNameJoin[{scriptDir,"generated"}];
nlo=FileNameJoin[{qsStage3GeneratedRoot,"DiracTripletEW_QCDEW_CT.nlo"}];
If[!FileExistsQ[nlo],
   Print["ERROR: Missing ",nlo, ". Run stage 02 and inspect errors."];
   Exit[2]
];
qsStage2AuditFile=FileNameJoin[{qsStage3GeneratedRoot,"stage2_log_audit.json"}];
If[!FileExistsQ[qsStage2AuditFile],
   Print["ERROR: stage-2 NLOCT audit manifest missing; use ",
         "scripts/run_triplet_ew_stage2.sh before UFO export."];
   Exit[3]
];
qsStage2Audit=Import[qsStage2AuditFile,"RawJSON"];
If[!TrueQ[Lookup[qsStage2Audit,"preflight_pass",False]],
   Print["ERROR: stage-2 audit failed; candidate UFO export refused."];
   Exit[3]
];
If[FileHash[nlo,"SHA256","HexString"]=!=Lookup[qsStage2Audit,"nlo_sha256",""],
   Print["ERROR: NLOCT .nlo file has changed since stage-2 audit."];
   Exit[3]
];
qsProjectRoot=DirectoryName[DirectoryName[scriptDir]];
qsContentReportFile=FileNameJoin[{qsProjectRoot,"results",
                                  "triplet_ew_nloct_content_report.json"}];
If[!FileExistsQ[qsContentReportFile],
   Print["ERROR: large NLOCT content report is missing; run ",
         "python3 scripts/inspect_triplet_nloct_large.py first."];
   Exit[4]
];
qsContentReport=Import[qsContentReportFile,"RawJSON"];
If[Lookup[qsContentReport,"status",""]=!=
     "CONTENT_STRUCTURAL_SCREEN_PASS_NOT_PHYSICS_CERTIFIED" ||
   Lookup[qsContentReport,"sha256",""]=!=FileHash[nlo,"SHA256","HexString"],
   Print["ERROR: NLOCT .nlo content inspection failed, or saved SHA256 ",
         "does not match the file being exported."];
   Exit[4]
];
qsHeaderAuditFile=FileNameJoin[{qsProjectRoot,"results",
                                 "triplet_ew_nloct_header_audit.json"}];
If[!FileExistsQ[qsHeaderAuditFile],
    Print["ERROR: missing UV/R2 external-leg header audit. Run ",
          "python3 scripts/inspect_triplet_nloct_vertex_headers.py ",
          "before candidate UFO export."];
    Exit[4]
];
qsHeaderAudit=Import[qsHeaderAuditFile,"RawJSON"];
If[Lookup[qsHeaderAudit,"status",""]=!=
     "HEADERS_PARSED_REQUIRE_PHYSICS_VALIDATION" ||
   Lookup[qsHeaderAudit,"file_sha256",""]=!=FileHash[nlo,"SHA256","HexString"],
    Print["ERROR: UV/R2 vertex headers are malformed or have an ",
          "input hash mismatch; do not export the UFO."];
    Exit[4]
];
Print["NLOCT_UV_R2_EXTERNAL_HEADERS_PARSED_ONLY"];
Print["NLOCT_CONTENT_STRUCTURAL_SCREEN_PASS_NOT_PHYSICS_CERTIFIED"];
Print["FR_STAGE2_PREFLIGHT_PASS_NOT_PHYSICS_CERTIFIED"];
AppendTo[$Path,frHome];
SetDirectory[frHome];
Get["FeynRules.m"];
LoadModel[FileNameJoin[{frHome,"Models","SM","SM.fr"}],
          FileNameJoin[{scriptDir,"TripletDirac.fr"}]];
FeynmanGauge=True;
(* Fail before reading the 136 MiB CT definitions if FeynRules shadows
   a short path-variable name. Directory[] at this moment is frHome. *)
If[!DirectoryQ[qsStage3GeneratedRoot],
   Print["ERROR: generated output root does not exist: ",qsStage3GeneratedRoot];
   Exit[5]
];
qsStage3BeforeLoadDir=Check[SetDirectory[qsStage3GeneratedRoot],$Failed,
                             {SetDirectory::badfile,SetDirectory::cdir}];
If[qsStage3BeforeLoadDir===$Failed ||
   Directory[]=!=qsStage3GeneratedRoot,
   Print["NLOCT_STAGE3_REJECTED_BAD_OUTPUT_DIRECTORY_PRELOAD: ",Directory[]];
   Exit[5]
];
Print["NLOCT_STAGE3_OUTPUT_DIRECTORY_PREFLIGHT_OK: ",Directory[]];
If[Environment["QUIRKSIM_STAGE3_PREFLIGHT_ONLY"]==="1",
   Print["NLOCT_STAGE3_PREFLIGHT_ONLY_PASS: no heavy counterterms loaded; exiting."];
   Exit[0]
];
Print["NLOCT_STAGE3_LOADING_FULL_COUNTERTERMS: large input file ",
      FileByteCount[nlo], " bytes. Inspect available memory."];
qsReadNLOCT=Check[Get[nlo],$Failed,
    {Part::partd,Part::partw,Set::partd,Function::slotn,StringJoin::string,
     SetDirectory::badfile,SetDirectory::cdir,Power::infy,Infinity::indet}];
If[qsReadNLOCT===$Failed,
   Print["ERROR: NLOCT counterterm definitions failed during Get[nlo]."];
   Exit[5]
];
If[!ValueQ[UV$vertlist]||!ValueQ[R2$vertlist]||
   !ListQ[UV$vertlist]||!ListQ[R2$vertlist],
   Print["ERROR: NLOCT output lacks parseable UV$vertlist or R2$vertlist."];
   Exit[5]
];
If[Length[UV$vertlist]==0||Length[R2$vertlist]==0,
   Print["ERROR: NLOCT output has no UV or R2 vertices."];
   Exit[5]
];
If[!ValueQ[FR$CTparam] || !ListQ[FR$CTparam],
   Print["ERROR: NLOCT source lacks its FR$CTparam list; ",
         "UFO export would drop required auxiliary counterterm definitions."];
   Exit[5]
];
Print["NLOCT_STAGE3_CT_PARAMETER_COUNT: ",Length[FR$CTparam]];

Print["NLOCT_STAGE3_COUNTERTERMS_LOADED: UV=",Length[UV$vertlist],
      ", R2=",Length[R2$vertlist]];
qsCandidateName="DiracTripletEW_QCDEW_CANDIDATE_UFO";
qsCandidatePath=FileNameJoin[{qsStage3GeneratedRoot,qsCandidateName}];
If[DirectoryQ[qsCandidatePath],
   Print["ERROR: candidate UFO directory already exists; do not overwrite ",
         qsCandidatePath];
   Exit[6]
];
qsStage3SetDirectory=Check[
    SetDirectory[qsStage3GeneratedRoot],
    $Failed,
    {SetDirectory::badfile, SetDirectory::cdir}
];
If[qsStage3SetDirectory===$Failed ||
   Directory[]=!=qsStage3GeneratedRoot,
   Print["NLOCT_STAGE3_REJECTED_BAD_OUTPUT_DIRECTORY: expected ",
         qsStage3GeneratedRoot, " but got ", Directory[]];
   Exit[6]
];
Print["NLOCT_STAGE3_OUTPUT_DIRECTORY_OK: ",Directory[]];
qsExportResult=Check[
    WriteUFO[
      LTripletFull,
      UVCounterterms -> UV$vertlist,
      R2Vertices -> R2$vertlist,
      CTParameters -> FR$CTparam,
      AddDecays -> False,
      Output -> qsCandidateName
    ],
    $Failed,
    {Part::partd,Part::partw,Set::partd,Function::slotn,StringJoin::string,
     Power::infy,Infinity::indet}
];
If[qsExportResult===$Failed,
   Print["NLOCT_STAGE3_REJECTED: Mathematica reported a structural error ",
         "in UFO generation; the partial directory, if any, is unvalidated."];
   Exit[7]
];
(* FeynRules 2.3.x only needs a separate CT_parameters.py if NLOCT
   supplied auxiliary CT parameters. A valid FR$CTparam == {} cannot
   require an otherwise unused file. *)
qsMandatory=Join[
    {"particles.py","vertices.py","parameters.py","couplings.py",
     "CT_vertices.py","CT_couplings.py"},
    If[Length[FR$CTparam]>0,{"CT_parameters.py"},{}]
];
qsMissing=Select[qsMandatory,
    !FileExistsQ[FileNameJoin[{qsCandidatePath,#}]] &];
If[Length[qsMissing]>0,
   Print["NLOCT_STAGE3_REJECTED: missing UFO files: ",qsMissing];
   Exit[8]
];
Print["NLOCT_STAGE3_CANDIDATE_FILES_PRESENT_UNVALIDATED: ",
      qsCandidatePath];
Print["UFO candidate exported. Its name is deliberate: NOT CERTIFIED."];
Print["STOP: Do not quote NLO EW production or use as an approved UFO before
 UV/IR pole cancellation, photon-PDF, SM reduction, and SU2 Ward identity tests."];
