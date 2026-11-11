# Literature-motivated angular emission budgets

## Result

Implemented a necessary low-angular-momentum reachability check for a class
of discrete emission prescriptions. This avoids simulating up to 1e8 events
per trajectory when angular capture is already excluded by a triangle
inequality. It is not a newly derived nonperturbative emission rate.

With emission energy >=Lambda and angular step <=1 hbar, none of the 14
isotropic fits or 22 collinear 2 T fits can reach L<=1. All six zero-field
collinear fits already have L=0. Even the broader energy >=0.1 Lambda and
step <=10 scenario leaves those same 36 high-L cases inaccessible. This is
stronger than the previous gate-closed statement within this declared class:
no gate or time limit is imposed, and compact binding energy is generously
credited to the budget before capture.

## Primary-source reuse matrix

| Source / locator | Verified content | Transfer decision |
|---|---|---|
| Kang and Luty, arXiv:0805.4642v3, Sec.5.5 Eq.(5.56)--(5.58), Fig.8, PDF pp.31--33 | Quantum toy emissions motivate angular randomization, not deterministic downward drift | Motivation only; toy QCD process and low-wave initial conditions are not our neutral high-L triplet evolution |
| Same source, Sec.5.7, PDF pp.35--36; Sec.7.2, PDF p.41 | IC glueball emission at close separation; order-Lambda energy, order-one total angular momentum; efficient/inefficient scenarios | Reuse parametric scales; hard floors/ceilings below are additional declared assumptions, not source theorems |
| Li et al., arXiv:2311.15486v2, Sec.III.A.4 | Geometry-dependent glueball eligibility, magnetic/material effects; no calibrated high-L drift | Retain the contact motivation; no new large directed torque inferred |
| Curtin et al., PRD113 (2026)015010, DOI10.1103/g7c7-6qh2, Sec.II, PDF p.2 | Glueball deexcitation rate remains difficult to determine; their study uses microscopic, SM-colored quirks | Their hadronic shower setup does not supply a neutral 774-eV-triplet radiation law |

Links: https://arxiv.org/pdf/0805.4642 ;
https://arxiv.org/html/2311.15486v2 ;
https://journals.aps.org/prd/pdf/10.1103/g7c7-6qh2 .

Scoped search on 2026-11-12: queries `quirks glueball radiation angular
momentum random walk annihilation Kang Luty`, `quirkonium deexcitation angular
momentum radiation random walk 0805.4642`, and `quirks glueball angular
momentum random walk`; checked these primary texts and their immediate
relevant discussion. No exhaustive citation-graph or repository survey was
performed in this step. No externally matched emission kernel was found in
these inspected sources.

## Necessary-condition calculation

Let each emission consume at least omega_min and change the angular-action
vector by at most j_max. Credit energy E_available=E_fit+B_compact. Then

    N_max=floor(E_available/omega_min)
    L_final >= max(0,L_initial-N_max*j_max).

The triangle inequality holds for every sequence, including maximally aligned
downward changes. If this lower bound exceeds L_capture, low L cannot be
reached at any emission within the budget. The result does not require an
unbiased random walk. In the unit-step/unit-energy card N_max is about
5.3e7--1.46e8, while isotropic angular gaps are 3.61e11--5.65e12 hbar.
Collinear 2 T gaps are 2.10e10--3.41e11 hbar.

An additional diagnostic assumes zero conditional mean of vector increments.
For their sum M_k, E|M_N|^2<=N*j_max^2. The nonnegative-submartingale maximal
inequality bounds the probability of ever reaching the low-L ball by
min(1,N*j_max^2/(L_initial-L_capture)^2). Deterministically inaccessible cases
have probability zero in this prescription. This is our mathematical bound,
not a probability quoted from the quirk literature. It is not the final-wave
1/sqrt(N) estimate for a source initially near low partial waves.

## Validation and scope

378 cards cover energy floors .1,1,10 Lambda and step ceilings 1,10,100 hbar
for each of the 42 level-8 neutral fits. At the extreme .1-Lambda/100-hbar
card, 11 collinear 2 T cases become kinematically reachable, whereas all
isotropic cases remain inaccessible. Reachable does not mean captured.

Three new tests check deterministic thresholds, zero-energy/invalid limits,
and compare the maximal bound with exhaustive first-passage enumeration of
all 1024 ten-step unbiased walks. Together with seven existing evolution
tests, all ten pass. This is internal mathematical validation; there is no
external physical reproduction of our hard emission bounds.

Release decision: conditional low-L inaccessibility can be reported for this
specified class. No unique physical annihilation probability is released.
Soft massless radiation, high-angular-momentum glueball orbital emission,
coherent torque and exchange of angular momentum with string modes are not
excluded. Endpoint fitted central L is not independently verified as the
total string action or a quantum occupation distribution. These limitations
cannot be resolved by asserting a larger directed drift.

Reproduce:

```bash
python scripts/check_triplet_angular_budget.py \
  studies/data/triplet_ensemble_capture_models_2026-11-11.json.gz \
  --output angular_budget.json
```

The saved JSON records input and original AF archive hashes and every card.
