(* A fast FeynRules tree-level model sanity probe.
   Run before EW counterterm renormalization:
     wolframscript -file theory/dirac_triplet_ew/00_tree_vertex_probe.wls

   This stage does NOT calculate any UV/R2. It prints vertices
   for manual checks against the existing NLO-QCD UFO and charge
   conservation. A positive run does not certify loop renormalization.
*)

qsFrHome=Environment["FEYNRULES_HOME"];
If[!StringQ[qsFrHome] || StringLength[qsFrHome]==0,
    Print["ERROR: FEYNRULES_HOME is unset"];Exit[2]];
qsFrHome=ExpandFileName[qsFrHome];
qsHere=DirectoryName[ExpandFileName[$InputFileName]];
qsSmFr=FileNameJoin[{qsFrHome,"Models","SM","SM.fr"}];
qsTripletFr=FileNameJoin[{qsHere,"TripletDirac.fr"}];
If[!FileExistsQ[qsSmFr] || !FileExistsQ[qsTripletFr],
   Print["ERROR: Missing SM or triplet FeynRules source."];Exit[2]];
AppendTo[$Path,qsFrHome];
SetDirectory[qsFrHome];
Get["FeynRules.m"];
LoadModel[qsSmFr,qsTripletFr];
FeynmanGauge=True;
qsTripletClassCheck=Table[
    !FreeQ[M$ClassesDescription,F[qsIndex]],{qsIndex,{61,62,63}}];
Print["FR_TRIPLET_CLASS_DECLARATIONS_PRESENT: ",
      InputForm[qsTripletClassCheck]];
If[!And@@qsTripletClassCheck,
   Print["ERROR: Missing one or more Dirac triplet classes after LoadModel."];
   Exit[3]];
qsTripletLagrangian=LTripletDirac;
qsTripletLagrangianText=ToString[InputForm[qsTripletLagrangian]];
Print["FR_TRIPLET_LAGRANGIAN_HEAD: ", Head[qsTripletLagrangian]];
Print["FR_TRIPLET_LAGRANGIAN_PREVIEW: ",
      StringTake[qsTripletLagrangianText,UpTo[1000]]];
If[Head[qsTripletLagrangian] === Symbol ||
   qsTripletLagrangian === 0 ||
   StringContainsQ[qsTripletLagrangianText,"bar[Qp1]"] ||
   StringContainsQ[qsTripletLagrangianText,"bar[Qzero]"] ||
   StringContainsQ[qsTripletLagrangianText,"bar[Qm2]"],
   Print["ERROR: triplet Lagrangian missing, zero, or contains unrecognized bar[psi] antifields."];
   Exit[3]
];
Print["FR_TRIPLET_LAGRANGIAN_DECLARATION_OK"];
Print["FR_TREE_PROBE_START"];
qsV=FeynmanRules[qsTripletLagrangian,ScreenOutput->False,MaxParticles->3];
Print["FR_TREE_VERTEX_COUNT: ",Length[qsV]];
If[!ListQ[qsV] || Length[qsV]<8,
   Print["ERROR: fewer than 8 triplet tree vertices; cannot proceed to EW NLO."];
   Exit[4]];
Do[
 Print["FR_TREE_VERTEX_",qsK,": ",InputForm[qsV[[qsK]]]],
 {qsK,Length[qsV]}];
Print["FR_TREE_PROBE_COMPLETED_UNVALIDATED: inspect FFV vertices against QCD UFO."];
