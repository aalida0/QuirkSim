(* Set FEYNRULES_HOME and execute with Mathematica or wolframscript. *)
ClearAll[frHome, scriptDir];
frHome = Environment["FEYNRULES_HOME"];
If[!StringQ[frHome] || StringLength[frHome]==0,
 Print["Set FEYNRULES_HOME"]; Abort[]];
scriptDir = DirectoryName[ExpandFileName[$InputFileName]];
SetDirectory[frHome];
Get["FeynRules`"];
LoadModel[FileNameJoin[{frHome,"Models","SM","SM.fr"}],
  FileNameJoin[{scriptDir,"DiracTripletEW.fr"}]];
FeynmanGauge = False;
Print[CheckHermiticity[LTriplet]];
Print[FeynmanRules[LTriplet, FlavorExpand -> True]];
SetDirectory[scriptDir];
WriteUFO[L, Output -> "DiracTripletEW_FR_UFO", AddDecays -> False];
Print["Runtime validation of generated vertices remains required."];
