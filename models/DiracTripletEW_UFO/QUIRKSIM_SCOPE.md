# QuirkSim DiracTripletEW_UFO (LO electroweak spectator-infracolor benchmark)

Independently added SU(2)L Dirac triplet Y=0 to the canonical MadGraph SM UFO.
Three Dirac fields q1+ (PDG 9000111), q0 (9000112), q2- (9000113) and three separate antiparticles.
Tree-level vectorlike A/Z/W couplings; full T+=sqrt(2) ladder => W coefficient e/sin(thetaW).
Eight FFV new vertices. Neutral q0 has no diagonal A or Z coupling.
Tree-level gauge-invariant mass choice is MQC=MQ0=225 GeV, the default.
The loop-improved charged mass may be set to 225.15717 GeV at production-level accuracy.
No dynamic SU(N_IC) gauge gluons or SU(N_IC) index. For pure LO Drell–Yan, multiply cross sections by N_IC if spectator-color assumptions apply.
No weak pion decays, string dynamics, hidden radiation, or annihilation. New particles width=0 so propagation/decay is done by QuirkSim.
Do not confuse this with a full quirk UFO. Matrix-element-level MadGraph numerical checks are still pending.

Provenance: Base SM UFO from mg5amcnlo/mg5amcnlo '3.x' models/sm, extended by QuirkSim.
Companion model builder, tests and source are in scripts/ and FeynRules/ of this branch.

Example (inside MadGraph):
import model DiracTripletEW_UFO
generate u u~ > q1+ q1- QED=2 QCD=0
output triplet_uu -f

Charge-conserving W+ example:
generate u d~ > q1+ q0~ QED=2 QCD=0
