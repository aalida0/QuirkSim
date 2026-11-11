# First public-backend weak-recoil result — 2026-11-11

Verdict: deterministic directional recoil is resolved on the micrometre scale;
selected local handoff checks pass. Weak event-time precision and benchmark input
precision require correction before statistical ensembles. No annihilation or
survival inference, no general convergence certificate.

Input: quirkpheno_weak.tar.gz, SHA-256
98630d08a5b126753c0e71294516cf40afe65792a05ec331e0d5da82dbfd2587.
Upstream d91cb27ed3786bc63900faad892513841b5b6ed4; adapter matches d846f752
before the subsequent UI precision change. Source patch hash matches manifest.
All eight compiled/running cases exit zero with 130 ps coverage, one event and
both endpoints, expected charged/neutral final species, and no G4Exception marker.

## Local handoff evidence

Each weak run has exactly one pion decay and one consumed endpoint-ownership
transfer. Both history fingerprints match across handoff. The neutral daughter
has mass 225 GeV and zero charge; partner remains charged. Parent/daughter and
pion identities agree, including the negative endpoint. Observed four-vectors
satisfy the local sum with maximum component residual 1.461e-11 MeV, and all
logged mass-shell checks pass. The pion rest-frame direction matches the case.
These checks do not establish total endpoint-plus-string invariants or complete
causal evolution. A fingerprint match is a targeted state-preservation check.

## Measured angular recoil effect

Both pion directions are physical conditional two-body decays, at a fixed
proper threshold 50 ps. They share exactly the pre-decay parent position,
four-vector and event time within a given level. Only the rest-frame pion
+y/-y direction changes; charge and species changes are shared.
This isolates sensitivity to recoil direction, not recoil vs zero recoil.

At lab time 130 ps, magnitudes of y(+y pion)-y(-y pion):

| Resolution | Decaying endpoint (micrometres) | Partner (micrometres) |
|---|---|---|
| 7 | 12.05966746 | 0.69928700 |
| 8 | 12.06201656 | 0.70006584 |

Level-7/8 changes of these angular differences are 2.3491 nm and 0.77884 nm:
about 0.0195% and 0.111% of the respective signals. The direction effect is
large relative to this selected numerical comparison. A single fixed direction
pair does not give an angular distribution, efficiency, or survival probability.
The plotted partner direction difference is zero at all 2001 checked samples
outside a conservative Euclidean future light cone. The observed later response
is consistent with retardation; this is not a proof of propagation along the
actual worldsheet or exhaustive absence of acausal error.

## Precision limits exposed by this pilot

Individual weak trajectories differ by as much as 1.59242 nm between levels,
about 6.25 confinement reference lengths. The stable control remains near
0.02198 nm (0.0862 reference lengths). Therefore the earlier stable criterion
must not be transferred to the weak trajectories. Micrometre directional
sensitivity is supported; quantum-scale close encounters remain unvalidated.

Endpoint-0 weak events occur at lab 54.46270293 ps (level 7) and 54.44528739 ps
(level 8). Proper-clock overshoots above the fixed 50 ps threshold are
0.01837920 and 0.00220040 ps. Endpoint-1 overshoots are 0.01119739 and
0.00877633 ps. First-post-step event selection and upstream approximate clock
integration need better localization before lifetime ensembles. Never claim
these clocks are exact proper-time roots.

The observed charged mass is 225.157 GeV, while the manifest requests
225.1571702763212 GeV. The actual gap is 157.0 rather than 157.170276 MeV.
The backend pion mass is 139.5701 MeV, compared with QuirkSim's reference
139.57039 MeV. The corresponding rest recoil is 71.87199 vs 72.24235 MeV,
approximately 0.51% difference. Local conservation uses actual backend masses
and is valid; the intended benchmark parameters were not realized exactly.

Geant4 v11.4.1 G4UIcmdWithADoubleAndUnit::DoIt normalizes units through
G4UIcommand::ConvertToString; the latter uses stream default precision unless
G4UImanager::UseDoublePrecisionStr(true) is enabled. The adapter now enables
that documented switch before macro commands. The same change also preserves
full requested initial momentum precision. Nine exact source preimages and
patch application are checked; this additional UI change has not run on AF yet.
Primary API sources:
- https://github.com/Geant4/geant4/blob/v11.4.1/source/intercoms/src/G4UIcmdWithADoubleAndUnit.cc
- https://github.com/Geant4/geant4/blob/v11.4.1/source/intercoms/src/G4UIcommand.cc
- https://github.com/Geant4/geant4/blob/v11.4.1/source/intercoms/include/G4UImanager.hh

## Reproduction and next step

    python scripts/audit_quirkpheno_weak_pilot.py quirkpheno_weak.tar.gz --output audit.json --plot recoil.svg

Plot uses separate pre/post-decay interpolants; it never interpolates momentum
through the recoil jump. Solid curves are level 8 and dashed level 7; lower
panel shows the change in the angular difference, in nm. No uncertainty band
or statistical confidence interval is inferred from two resolutions.

Next: localized decay time and consistent high-precision inputs, then sampled
proper lifetimes/channels and two decays. Keep physical opposite-direction
controls; a stable comparison also changes mass and charge. Doublet stays paused.
No additional AF run was submitted during this audit. Radiation/annihilation
kernels, encounter matching and physical survival remain unresolved.
