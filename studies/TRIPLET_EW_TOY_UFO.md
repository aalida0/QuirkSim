# Dirac triplet EW toy UFO — first implementation

This is a **tree-level electroweak-only production model**, NOT a full quirks UFO.

## Fields
A vectorlike SU(2)_L triplet with Y=0:

- `q1+` and its antiparticle `q1-`: PDG ±9000111, charge ±1
- `q0` and `q0~`: PDG ±9000112, charge 0
- `q2-` and its antiparticle `q2+`: PDG ±9000113, charge ∓1

The label q1- denotes an *antiparticle*, whereas q2- is a *different component* of the original triplet.

Inputs: m(Q0)=225 GeV, m(charged)=225.15717 GeV. Loop splitting is input, not calculated.

Gauge couplings from the triplet covariant derivative: photon: ±e; Z: ±e cot(theta_W); W transition: g=e/sin(theta_W). Both chiralities participate.

## Building against MadGraph SM UFO

```bash
git switch feature/dirac-triplet-ew-toy-ufo
python scripts/build_triplet_toy_ufo.py \
  --sm-ufo ~/codes/Herwig/herwigpp/opt/MG5_aMC_v3_5_1/models/sm \
  --out ~/codes/Herwig/herwigpp/opt/MG5_aMC_v3_5_1/models/TripletEW_UFO

~/codes/Herwig/herwigpp/opt/MG5_aMC_v3_5_1/bin/mg5_aMC
```

Inside MadGraph:

```text
import model TripletEW_UFO
display particles
generate u u~ > q1+ q1- QED=2 QCD=0
output triplet_ew_uu -f
launch triplet_ew_uu
```

To test a charged-current process:

```text
import model TripletEW_UFO
generate u d~ > q1+ q0~ QED=2 QCD=0
output triplet_ew_ud -f
```

For proton-proton production at 13.6 TeV, set ebeam1 = ebeam2 = 6800 GeV in the run card. The model contains `q1+`, `q1-`, `q2+`, and `q2-` as different species with distinct PDG codes.

## Scope / limitations

1. **No infracolor representations**. The `color=1` assignment is an explicitly unphysical EW production-only proxy, and absolute cross sections cannot yet be interpreted as physical quirks rates. An IC fundamental introduces additional multiplicity and hidden interactions.
2. **No quirk strings, weak pion decays, capture, radiation, annihilation or ATLAS detector**. Those require QuirkSim or later simulation.
3. **No automatic one-loop splitting**. The masses are benchmark values. Gauge invariance of a split multiplet is approximate at this order, so compare with a degenerate benchmark when testing the model.
4. **Validation status**. Offline structural tests pass. The installed MG5 runtime was unavailable in the authoring environment, so MG5 import, numerical production rates, and process-level tests remain to be run.

Do not claim any physics cross section, ATLAS sensitivity, or UFO validation without the MadGraph execution and independent reference comparisons.
