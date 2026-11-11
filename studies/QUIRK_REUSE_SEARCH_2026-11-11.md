# Literature and code that can help QuirkSim

Search date: 11 November 2026. Scope: electroweak, QCD-color-neutral heavy
fermionic quirks, charged-to-neutral endpoint decays, classical string
evolution, radiation and quantum annihilation. Current cards are doublet
150 GeV / 2154 eV and triplet 225 GeV / 774 eV. This is a targeted reuse
search, not an exhaustive novelty certification. No external code was compiled
or shown to reproduce our benchmark in this search. No authors were contacted.

## Prioritized reuse matrix

| Source and inspected locator | What it supplies | Applicability / evidence | Priority |
|---|---|---|---|
| [Athena Quirks source](https://github.com/rolandjansky/athena/tree/master/Simulation/G4Extensions/Quirks/src), `InfracolorForce.cxx`, `HyperbolaStep.cxx`, `HyperbolaStepper.cxx` | Deque of string four-vectors, incoming/outgoing reflection, segment-wise hyperbolic propagation, rapidity limiter, midpoint external field update, string angular momentum and energy moment | Source inspected. Most concrete independent numerical route found. Geant4/CLHEP integration; not a standalone library or a matched radiation/annihilation generator. Existing user Athena workflow gives an integration route, but its version must be pinned separately. | 1, numerical |
| G. S. Sharov, [hep-ph/0004003](https://arxiv.org/pdf/hep-ph/0004003), Sec. II, Eq. (3) and initial-boundary construction | Relativistic string with massive ends; continuation of characteristic data through massive boundary conditions; perturbed rotations | Full text inspected. Useful analytic benchmark/alternate parametrization. Does not include our magnetic field, charge-changing weak event, or infracolor radiation. No accompanying executable verified. | 1, formulation |
| Knapen, Lou, Papucci, Setford, [1708.02243](https://arxiv.org/pdf/1708.02243), Secs. II–III | Phenomenological trajectory simulation and planar-hit strategy for vectorlike fermions | Full text inspected. Good trajectory benchmark lead, not a certified internal-string solver following arbitrary weak recoil. Macroscopic tracking accuracy is not core-scale encounter accuracy. | 1, benchmark/code lead |
| Sha et al., [2410.00269v2](https://arxiv.org/html/2410.00269v2), Sec. III Eq. (2), acknowledgements; journal DOI [10.1103/dvzc-lgj1](https://doi.org/10.1103/dvzc-lgj1) | Explicit acknowledgement of trajectory code shared by Knapen/Papucci | Verified code-existence lead, NOT a verified public download or reuse license. Its central-force trajectory model supplies a useful control but does not establish arbitrary internal-string dynamics. Consider an author request later; none sent. | 1, acquisition lead |
| Li, Li, Pei, Zhang, [2002.07503](https://arxiv.org/abs/2002.07503), [EPJC 80, 651](https://doi.org/10.1140/epjc/s10052-020-8209-y), Secs. 2–6 | Approximate analytic trajectory, oscillation scales and external-force corrections | Publisher full text inspected. Useful zero/small-field controls and warnings about string-direction assumptions. No public implementation verified. No direct weak-impulse or quantum matching solution. | 2, benchmark |
| Kang, Luty, [0805.4642v3](https://arxiv.org/pdf/0805.4642), Sec. 5, especially Eq. (5.29) | WKB encounter annihilation prescription and discussion of angular momentum/radiation | Full text inspected. Useful conditional framework. Neither the populated quantum density nor arbitrary post-decay mixed-flavor rates are supplied by this equation. Small-loss per-crossing probabilities need a stated survival prescription. | 1, quantum matching |
| Evans, Luty, [1811.08903v1](https://arxiv.org/pdf/1811.08903), Sec. 4 Eq. (4.4), App. B.1–B.2 | Explicit phenomenological glueball emission probability, angular-momentum random walk, electromagnetic cooling regimes | Full text inspected. Assumes emission probability 0.1 for its study. Provides a declared comparison model, not a measured/calibrated three-dimensional recoil kernel. Its material-bound regime is not our vacuum evolution. | 1, uncertainty model |
| Li, Liao, Ni, Pei, [2311.15486](https://arxiv.org/abs/2311.15486), [PRD 109, 095005](https://scoap3-prod-backend.s3.cern.ch/media/files/85154/10.1103/PhysRevD.109.095005.pdf), Sec. III A Eq. (7), Table II, Eq. (26) | Assumed glueball emission per oscillation, linear/Coulomb cooling, partial s-wave hard annihilation coefficients, radial Schrodinger shooting prescription | Full published text inspected. Separate fermion/scalar coefficients; Table II is explicitly same-flavor. Does not close W/Z channels for all our multiplet/flavor states or give their populations. Period conventions and scheme/color factors require matching. | 1, conditional rates |
| Fok, Kribs, [1106.3101v1](https://arxiv.org/pdf/1106.3101), Sec. V, Apps. A–C | Spin/projector formalism, charged and neutral quirkonium decays, vectorlike versus chiral comparison, transitions | Full text inspected. Important source for missing electroweak hard channels and transition checks. Sec. V vectorlike example has different electroweak charges from our charged/neutral doublet; triplet is not a plug-in replacement. Re-derive vertices/group factors and check normalizations before reuse. Does not determine our incoming quantum state. | 1, hard channels |
| Harnik, Kribs, Martin, [1106.2569v1](https://arxiv.org/pdf/1106.2569), Sec. IV and App. B | Weak decay competing with spin-down, altered flavor content and subsequent annihilation | Full text inspected. Direct conceptual predecessor: weak decay before/after cooling is already studied. Model uses QCD-colored scalars, so numerical widths/radiation are not our fermion model. | 1, scope/novelty |
| Forsyth, Low, Tenney, Verhaaren, [2504.02940v2](https://arxiv.org/html/2504.02940v2), Secs. 2.1, 4.1.1, App. A.1; [2026 erratum](https://doi.org/10.1103/lwnc-cmpn) | Electroweak scalar-quirk dynamics, beta decay, photon deexcitation and explicit discussion of uncertain glueball radiation | Full original HTML inspected; erratum publisher record/search text checked. Original omits deexcitation glueball radiation under a stated assumption. Scalar formulas do not transfer to fermions; conservative signal reasoning is specific to their observable. Erratum changes must be reconciled before reproducing plots or glueball lifetimes. | 1, modeling scope |

## Code audit and concrete implementation idea

The inspected public Athena mirror uses:

- `InfracolorForce::PushStringVector`, `PopTo`, `CombineStringVector` to
  exchange and optionally merge finite string segments between endpoints;
- `HyperbolaStep::PrepareNextStep` and `Step` to reflect incoming string
  four-vectors, update endpoint four-position/momentum, then add the field;
- a rapidity-based step limiter, rather than relying solely on an outgoing
  angular threshold as our current solver does.

This suggests an **independent piecewise-segment boundary propagator** as a
promising control for the doublet stiffness. That is a proposal from source
inspection, not evidence that it will fix our present error. The inspected
`HyperbolaStepper` sets error estimates to zero with a TODO; `DistChord` also
returns zero with a TODO. `HyperbolaStep` explicitly notes possible instability
in its mass-shell renormalization. Therefore it should not be treated as an
infallible reference for our sub-confinement-length closure gate.

Pinned file blob hashes in the mirror:

| File | Blob SHA |
|---|---|
| `InfracolorForce.cxx` | `35c91cefe63f8b38a4a694656eb4c9bd865d22fc` |
| `HyperbolaStep.cxx` | `d6f21f4c089cad414a16140341ea486d9a4583d5` |
| `HyperbolaStepper.cxx` | `589ae9adf263c14e5edddb7f4734ade6586303a2` |

Latest path commit returned for `HyperbolaStep.cxx`:
`9790641f007cf0dded380e3b4e7f6a5ebd7fda29` (2016-04-08).
Mirror root LICENSE states Apache 2.0 except where other licenses apply;
these source files retain CERN/ATLAS copyright headers. No source was copied
into QuirkSim. The current CERN GitLab directory could not be inspected with
the available web tool; this is an older public mirror, not verification of
current Athena behavior.

[TurinLatticeFieldTheoryGroup/NambuGotoCNF](https://github.com/TurinLatticeFieldTheoryGroup/NambuGotoCNF)
is public and MIT licensed, with a README linking [2307.01107](https://arxiv.org/abs/2307.01107).
Its purpose is lattice Nambu–Goto sampling with normalizing flows. It is not
a real-time open string/massive-endpoint evolution code, so it is not a fix
for our solver. README inspected, notebook not executed. Other generic
Nambu–Goto repository hits were not promoted to usable references without
a verified paper, boundary model and validation history.

The user's known quirk-chain workflow remains relevant, but a direct GitHub
README lookup at `Quirk-Analysis/quirk-chain` returned 404 in this session;
public access, license and exact revision were not verified here.

## Scientific decision

We do not need to invent a nonperturbative glueball calculation to make any
progress. Published analyses use explicitly assumed prescriptions, including
nonzero emission and omitted emission. Those can define **conditional model
comparisons**, provided energy/momentum accounting, period conventions and
the observable's sensitivity to recoil are checked. They are not calibrated
confidence intervals or universal upper/lower bounds on suppression.

For our central annihilation-suppression claim, the decisive question is
whether the weak-recoil effect survives defensible variations of these
assumptions and quantum-state matching. A detector-scale trajectory reference
cannot validate that by itself. A conditional classical-recoil study is a
narrower possible scope; its publishability and novelty remain unestablished.
Weak decay versus spin-down, by itself, is not new (2011 and 2025 predecessors).

Next sequence: (1) inspect stop diagnostics without loosening gates; (2)
compare against a source-matched segment propagator or analytic benchmark;
(3) map the Fok–Kribs hard channels to our exact representations; (4) run a
declared radiation/matching uncertainty comparison. No collider sensitivity
or unconditional physical survival claim is released.

## Search coverage ledger

Searches used web engines and GitHub repository metadata on 2026-11-11.
Queries included `quirks simulation code github Nambu Goto radiation annihilation`,
`quirks weak decay recoil electroweak multiplet strings annihilation`,
`site:github.com quirks infracolor`, `Learning to reconstruct quirky tracks github code`,
`Visible collider signals of natural quirks github`, `Stopping Quirks code simulation github`,
`Sharov string massive ends numerical initial boundary problem orthonormal gauge`,
and `Chiral quirkonium decays code supplementary Mathematica Fok Kribs`.
GitHub metadata queries included `quirk physics` (4 mostly unrelated hits),
`quirks particle` (0), `infracolor` (0), and `Nambu Goto` (first page inspected).
Repository discovery results are not a global absence proof.

Backward citation paths followed selectively: 2410.00269 -> 1708.02243;
1106.3101 -> earlier quirkonium work; the cooling discussions in 1811.08903,
2311.15486 and 2504.02940 -> 0805.4642. The 2025 paper's linked erratum was
checked. No full forward-citation crawl, INSPIRE graph enumeration, systematic
supplement download or private author-code acquisition was completed. These
limits prevent an exhaustive novelty or public-code-absence claim.
