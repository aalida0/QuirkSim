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
out=FileNameJoin[{scriptDir,"generated"}];
nlo=FileNameJoin[{out,"DiracTripletEW_QCDEW_CT.nlo"}];
If[!FileExistsQ[nlo],
   Print["ERROR: Missing ",nlo, ". Run stage 02 and inspect errors."];
   Exit[2]
];
AppendTo[$Path,frHome];
SetDirectory[frHome];
Get["FeynRules.m"];
LoadModel[FileNameJoin[{frHome,"Models","SM","SM.fr"}],
          FileNameJoin[{scriptDir,"TripletDirac.fr"}]];
FeynmanGauge=True;
Get[nlo];
If[!ValueQ[UV$vertlist]||!ValueQ[R2$vertlist],
   Print["ERROR: NLOCT output lacks UV$vertlist or R2$vertlist."];Exit[2]
];
If[Length[UV$vertlist]==0||Length[R2$vertlist]==0,
   Print["ERROR: NLOCT output has no UV or R2 vertices."];Exit[2]
];
SetDirectory[out];
WriteUFO[
  LTripletFull,
  UVCounterterms -> UV$vertlist,
  R2Vertices -> R2$vertlist,
  Output -> "DiracTripletEW_QCDEW_CANDIDATE_UFO"
];
Print["UFO candidate exported. Its name is deliberate: NOT CERTIFIED."];
Print["STOP: Do not quote NLO EW production or use as an approved UFO before
 UV/IR pole cancellation, photon-PDF, SM reduction, and SU2 Ward identity tests."];
