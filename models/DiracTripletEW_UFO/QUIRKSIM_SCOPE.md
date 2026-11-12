# QuirkSim Dirac triplet: LO and NLO-QCD production model scope

This directory (`DiracTripletEW_UFO`) is the **original LO electroweak
prototype**. It is validated for six hard Drell--Yan channels and their
charge assignments, including MadGraph LO amplitudes and 1,000-event
production. It has **no one-loop SM-QCD UV/R2 counterterms**.

For **NLO QCD**, use the reproducible model generator in:

```bash
python3 scripts/build_dirac_triplet_nlo_ufo.py \
  --mg5-dir /path/to/MG5_aMC --force
python3 scripts/validate_dirac_triplet_nlo_ufo.py \
  --mg5-dir /path/to/MG5_aMC
```

It produces `models/DiracTripletEW_NLO_UFO/`, based on MadGraph's NLO
`loop_sm` and retaining its validated SM-QCD UV/R2 sector. QuirkSim
adds the complete eight tree-level vectorlike triplet electroweak vertices
and correct charge-conjugate particles; the new triplet does not couple
to SU(3)c, so no additional BSM QCD counterterms are required.

### Lagrangian / mass convention

Three Dirac fields q1+ (PDG 9000111), q0 (PDG 9000112), q2-
(PDG 9000113), with three separate antiparticles. They form **one**
vectorlike SU(2)L triplet with Y=0; there are not six triplets.

Tree-level A/Z/W interactions have T3=(+1,0,-1), W coefficient g2
and identical vector couplings to left- and right-handed components.
The neutral field has no diagonal photon/Z tree coupling.

Gauge invariant tree mass is MQC=MQ0=M (default 225 GeV).
The 157 MeV charged-neutral mass splitting used in QuirkSim is a
one-loop *electroweak pole-mass correction*, not an independent
gauge-invariant tree-level mass term. For NLO-QCD production use
degenerate electroweak masses and apply splitting in the subsequent
weak-decay propagation. A split-mass hard UFO is only an approximation.

### What is not included

The triplet is represented as an SU(3)c singlet with an **unrepresented
infracolor index**. The NLO-QCD partonic cross section can be multiplied
by N_IC for degenerate spectator hidden colors. There are no hidden
gluons, infracolor radiation, hadronization, string worldsheets,
quirkonium/annihilation or weak pion decay vertices inside either UFO.
These require the independent QuirkSim simulation stage.

Neither UFO contains NLO electroweak corrections or a gauge-invariant
dynamical infracolor completion. The NLO-QCD UFO is complete **within
its explicitly restricted perturbative production scope**, not a full
physical model of all quirk processes.

See `studies/TRIPLET_NLO_QCD_UFO.md` for physics scope, source
provenance, NLO validation and numerical cross-section requirements.
