# Reuse existing calculations: frozen decisions and executable controls

The immediate development policy is to reuse the public transport, hard
amplitudes already calculated, and published radiation prescriptions. No
additional independent Sommerfeld engine or precision annihilation project
is needed for a first declared-model phenomenology study.

## Source and code reuse matrix

| Ingredient | Source | Decision and scope |
| --- | --- | --- |
| Geant4 quirk transport | tvami/QuirkPheno, pinned d91cb27ed3786bc63900faad892513841b5b6ed4 | Keep the existing Apache-2.0 backend and experimental weak-handoff adapter. Its README identifies transport/recasting functionality; it does not establish our radiation or quantum-state implementation. |
| Mean glueball energy loss | Li et al., arXiv:2311.15486v2, Eq. (7), Sec. III.A.4 and Fig. 9 | Reuse mean loss epsilon*Lambda/T. Use epsilon=0.1 and the paper's broad 0.01--1 variations as model scenarios, plus epsilon=0. These are not statistical errors. |
| Radial period | Forsyth et al., arXiv:2504.02940v2, Eq. (6) | Reuse T=4p/sigma for a full signed radial period. A factor-two mean-count scenario exposes encounter/period convention dependence. No silent conversion of source probabilities. |
| No-glueball baseline | Same source, Sec. 2.1 | A published modeling choice. Our neutral reference has no leading point-charge photon cooling; persistent charged photon cooling must not be assigned to neutral daughters. |
| Competition of beta decay, spin-down and annihilation | Harnik, Kribs and Martin, arXiv:1106.2569 | Reuse the workflow structure. Their colored scalar rates are not transferred to this color-neutral Dirac triplet. |
| Analytic Sommerfeld/partial-wave machinery | El Hedri, Kaminska and de Vries, arXiv:1612.02825; EPJC 77 (2017) 622; author code below | Use the published framework for future partial-wave dressing; retain spin/color projections. Do not build another general implementation. |
| Coupled electroweak capture | Asadi et al., arXiv:1610.07617, with erratum | Method reference. Majorana identical-particle normalization and spin restrictions require translation for our Dirac model. No verified public code integration is claimed here. |

## Author code actually executed

The paper's reference 44 explicitly links
[MDT-Maikel/Sommerfeld](https://github.com/MDT-Maikel/Sommerfeld).
Pinned commit: f0c97f32d10bb0f23f791765e652426117211955.
Upstream README declares CC-BY-4.0. We do not vendor the large notebook or
generated Python source. The upstream source is pinned by git blob
982db5d710283acced8e4d5b580784fd6807d93b in the validation adapter.

`check_author_sommerfeld.py` executes the author's Python function definitions
without rewriting them. Only the Python-2 CLI is excluded; its original
global `alphasommerfeld` is explicitly supplied because the generated
functions use that global despite a similarly named formal argument.
It checks a **color/spin-averaged SU(3) fundamental, l=0 control**:

* Bare leading sigma*v_Moller = 7*pi*alpha_s^2/(27*m^2).
* Coulomb enhancement = (2/7) S(4 alpha_s/3) + (5/7) S(-alpha_s/6).
* The author convention uses p/m=beta/sqrt(1-beta^2) in S, rather than
  substituting beta directly into the nonrelativistic momentum formula.

All checked normalizations agree to better than 1e-12 relative. At
alpha_s=alpha_Sommerfeld=0.04 the enhancement is 1.566384 at beta=0.05
and 1.041343 at beta=0.4. Results and source SHA-256 are in
`data/author_sommerfeld_control_2026-11-11.json`.

The SU(2) IC singlet is not the SU(3) color-averaged state. The generated
Python CLI supports representations 3,6,8; the notebook and stored results
include general SU(N). General-N notebook execution and a singlet-projected
transfer remain unvalidated. The l=0 control is not a reproduction of the
entire finite-velocity spin-dependent calculation or a confinement/capture
calculation. This code is a useful reusable component with a scoped passing
external control, not a turnkey simulator of the whole project.

Reproduce after obtaining the pinned `sommerfeld.py` from the author repo:

```bash
python scripts/check_author_sommerfeld.py /path/to/author/sommerfeld.py \
  --output author_sommerfeld_control.json
```

## Literature-based baseline now connected to existing code

`literature_deexcitation_baseline.py` reuses `OrbitAveragedModel` and
`linear_cooling_time_ps`; it does not implement another solver. The cards
cover charged-stable and neutral radial references, epsilon=0,0.01,0.1,1,
and a factor-two mean-count convention variation. The point-production
reference has beta_each=0.4; it is not the old pre-stretched conformal pilot.
The linear reference stops at the declared energy 8 Lambda, without
identifying that stop with capture or annihilation.

The shortest reference cooling time across these 16 cards is about
7.55e9 ps (7.55 ms), versus a charged endpoint weak lifetime of about
258.5 ps. This is a **counterfactual radial all-encounters-eligible reference**,
not the post-decay evolution of the G4 event. The mean emitted energy is a
loss parameter, not an asserted on-shell glueball mass or spectrum.

The script also consumes the existing first-return speed bounds. For radii
1,2,8 times hbar*c/Lambda_tension, the four recoil cases are outside the
gate throughout their uploaded common windows. The two stable controls
remain geometrically unresolved. One epsilon trial per eligible encounter
is a declared additional choice; the output gives eligibility/probability
bounds for that choice, not a physical annihilation probability. No fake
glueball recoil or energy loss is applied after the fact to unchanged tracks.
The saved bounds are in the lab frame. This gate is an explicitly labeled
diagnostic; translating the source's CM encounter criterion is still needed
before using it as a physical event-level emission prescription.

Important convention checks:

* Our Lambda_tension=sqrt(sigma). Forsyth et al. use sigma approximately
  3.6 Lambda_prime^2; these confinement scales must not be equated.
* The published photon coefficient and our rationalized-EM implementation
  have a normalization discrepancy already recorded in
  `RADIATION_ANNIHILATION_SOURCE_AUDIT_2026-11-10.md`. This reuse step does
  not silently switch the existing EM convention.
* Li et al.'s estimate of T and E0 is not equated to the exact relativistic
  full-signed-period formula without a declared energy/period convention.
* A zero-glueball baseline is not guaranteed conservative for every recoil
  observable; its direction of conservativeness must be established.

Reproduce:

```bash
python scripts/literature_deexcitation_baseline.py \
  studies/data/triplet_first_return_window_2026-11-11.json \
  --output literature_deexcitation_baseline.json
PYTHONPATH=. python -m unittest discover -s research/quirkpheno_recoil -p test_literature_baseline.py
```

## Remaining project-specific implementation

The public calculations cover much of the generic machinery. The novel
work is the charge-changing recoil interface, evolution of angular momentum
and spin between weak decays and repeated encounters, and their impact on
a stated de-excitation/capture model. Those modeling choices must be varied
alongside epsilon. The current executable reference study does not yet close
that full feedback loop or produce physical annihilation probabilities.

Release: use these source-backed model cards and reusable components for
the scoped phenomenology development. A publication claim of robust
recoil-induced annihilation suppression still requires an ensemble and
explicit tests of angular-momentum, capture and spin assumptions.

Sources inspected:

- https://arxiv.org/html/2311.15486v2
- https://arxiv.org/html/2504.02940v2
- https://arxiv.org/pdf/1106.2569
- https://link.springer.com/article/10.1140/epjc/s10052-017-5168-z
- https://github.com/MDT-Maikel/Sommerfeld
- https://github.com/tvami/QuirkPheno

This is a targeted reuse audit, not a new exhaustive novelty review.
