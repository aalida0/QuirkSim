# Dirac triplet LO UFO v0.2

Complete standalone SM+Dirac-triplet UFO at models/DiracTripletEW_UFO.
Three independent Dirac fields and three distinct antiparticles: q1+, q0, q2-.
Tree-level SU(2)L triplet Y=0, with vectorlike W+/- coupling e/sw.
Photon and Z couple only to the charged components (T3=+1,-1).
Tree-invariant default: MQ0=MQC=225 GeV.
Loop-split mass benchmark is a pole-mass approximation.
No dynamical infracolor, strings, radiation or weak pion decays in the UFO.
Apply N_IC multiplicity to LO electroweak cross sections, not amplitudes.
MadGraph and FeynRules runtime checks are still pending.

Example in MG5:
import model /full/path/to/QuirkSim/models/DiracTripletEW_UFO
generate u u~ > q1+ q1- QED=2 QCD=0
output dirac_triplet_pair -f

Associated W+ process:
import model /full/path/to/QuirkSim/models/DiracTripletEW_UFO
generate u d~ > q1+ q0~ QED=2 QCD=0
add process u d~ > q0 q2+ QED=2 QCD=0
output dirac_triplet_associated -f
