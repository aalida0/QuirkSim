# Charge-changing recoil connected to a declared repeated-encounter model

## Implementation and scope

The model now connects the uploaded post-weak trajectories to repeated
central-orbit encounters, mean radiative cooling, an orbital-angular-momentum
loss prescription, a capture criterion and competing annihilation channels.
It produces **conditional model probabilities**, with energy, angular-momentum
and probability ledgers. It does not claim to be a quantum or Nambu--Goto
continuation of the uploaded trajectories.

This implements the phenomenology route proposed in
`LITERATURE_REUSE_DECISIONS_2026-11-11.md`, reusing the existing radial cooling
integral and the previously calculated compact neutral spin-zero width.
Published motivation is Li et al. arXiv:2311.15486v2 Eq. (7) for mean loss,
and Forsyth et al. arXiv:2504.02940v2 Eq. (6) for the explicitly chosen
full signed period. Directed angular drift, capture and spin reset below are
additional declared closure choices, not results supplied by those papers.

## Handoff from charge-changing transport

Only the four weak-decay cases are fitted; all their window endpoints are
neutral. Their saved worldlines and momenta are Lorentz transformed to the
fixed production CM (beta_z=0.25) before finding simultaneous closest approach.
Interpolate both momenta at the same fitted time and define

r=(x0-x1), p_rel=(p0-p1)/2,
L=|r cross p_rel|/hbar,
E=2[sqrt(m0^2+p_rel^2)-m0]+sigma*|r|.

These define an **endpoint central-Hamiltonian reduction**, not the full
string energy or angular momentum. Residual pair endpoint momentum is
retained in the fit report. The model assumes no subsequent external fields,
material kicks or string-mode transfer and begins after the weak decays;
pre-handoff radiation is absent from this fit.

| Case | Level | Fitted L/hbar | Central pericentre (micrometres) |
| --- | ---: | ---: | ---: |
| Sampled recoil | 7 | 2.08426e12 | 4.188400 |
| Sampled recoil | 8 | 2.08428e12 | 4.188433 |
| Opposite recoil | 7 | 2.02251e12 | 4.063294 |
| Opposite recoil | 8 | 2.02254e12 | 4.063334 |

The frame and reduction differ from the previous lab-frame estimates. These
numbers must not be identified with total string angular momentum.

## Explicit evolution prescriptions

1. Turning geometry uses the central Hamiltonian
   H=2 sqrt(m0^2+p_r^2+L^2/r^2)+sigma*r. The inner turning root is solved
   from this Hamiltonian at each updated (E,L).
2. A hard emission gate requires r_peri <= 1/Lambda_tension. Outside it,
   mean cooling, angular loss and the supplied excited-annihilation hazard
   are disabled. This is a model cutoff, not proof of zero quantum tunnelling.
3. Within the gate, mean loss is epsilon*Lambda per full signed radial
   period, epsilon=0.01,0.1,1. The reference period is 4p/sigma. It has two
   radial encounters. The clock uses a radial approximation even for nonzero
   L; it becomes least controlled near circular or quantum-scale orbits.
4. Two angular prescriptions are tested: conserved L, and directed loss
   dL/dN_effective=-1 until zero. The latter is not an unbiased angular random
   walk. It parametrizes efficient damping; no real glueball momenta or spin
   spectrum are sampled. The secular update changes both energy and the next
   turning geometry, unlike a hazard evaluated on an unchanged trajectory.
5. Capture requires E <= 8 Lambda and L <= 1. It transfers the pair into the
   declared compact state sector and assigns E+binding_energy to a capture
   radiation bath, plus the remaining orbital L. This scalar energy/action
   ledger is not a local four-momentum-conserving microscopic capture process.
6. At capture, a supplied fraction f=0,1/4,1 enters the previously calculated
   neutral-connected spin-zero eigenstate. The rest remains in an unresolved
   captured sector. f=0 is a spin-preserving stress limit, 1/4 is a completely
   randomized-spin ansatz and 1 is an all-singlet stress limit. None is inferred
   from the spinless histories. The singlet width and binding energy are read
   from the pinned model result rather than recalculated.
7. A separate optional excited-state hazard has optical depth h=0 or 0.01
   per eligible radial encounter. This is a supplied overlap/hazard scenario,
   not a rate derived from the finite-velocity hard matrix. The default h=0
   tests capture-mediated annihilation only. Survival to capture is included
   before the compact decay contribution is added, avoiding double counting.

The secular clock permits fractional encounter counts and does not resolve
the phase of individual passages. The energy/action ledger follows the
conditional surviving orbit; it is not an ensemble-averaged radiated-energy
ledger weighted by survival. Histories are censored if an update leaves the
allowed central-orbit domain, rather than assigning a complete-horizon
probability in that case.

## Results of the connected study

The driver evaluates 96 orbit scenarios: four fitted cases, their zero-L
action references, three epsilon values, two angular-loss prescriptions and
two time horizons. Each orbit is reweighted over six capture-spin/excited-
hazard choices. Reweighting reuses the same secular path and does not rerun
the expensive compact eigenstate or hard-scattering calculations.

The zero-L reference has the same fitted excess energy and neutral mass but
artificially sets L=0. It is **not a physical no-recoil weak-decay control**;
it isolates dependence on the initial action in the declared model.

* All four fitted recoil orbits are outside the emission gate. Their actions
  remain fixed while encounters repeat. Conditional model capture and
  annihilation probabilities are zero at both 1000 ps and 1e13 ps after
  handoff. This zero is a consequence of the hard cutoff and isolated
  central-reduction assumptions, not a model-independent lifetime claim.
* For the zero-L action reference matched to sampled seed 1, level 8,
  epsilon=0.1 gives capture at 1.51092e11 ps, about 0.151 seconds. With h=0
  and f=1/4, the declared-model annihilation probability by ten seconds is
  1/4. With f=0, it is zero; the captured sector is retained explicitly.
* With h=0.01, the same zero-L action reference has an excited-annihilation
  probability of about 4.5% during the first 1000 ps of this secular clock.
  The fitted recoil orbit has no eligible encounter under the gate.
* A synthetic eligible-orbit probe, E=0.001 GeV and L=100, reaches the capture
  energy but fails the L threshold if L is conserved. With directed angular
  loss, it is captured. An additional E=40 GeV, L=1e7 test verifies that
  conserved L closes the emission gate during cooling, while directed loss
  permits capture. This exercises feedback that the fitted large-miss cases
  cannot constrain because their gates start closed.

The zero-L reference spans annihilation probabilities 0--1 over the full
closure set at ten seconds. Thus the integrated model exposes substantial
capture/spin/hazard dependence; it does not establish a robust physical
annihilation suppression factor for a realistic ensemble.

## Validation and release

Seven tests check the Hamiltonian turning root, the existing exact radial
cooling-time reference, repeated gated-off passages, angular-loss-dependent
capture, exact competing exponential hazards, finite observation windows,
spin-sector accounting, ledger closure and step refinement. The comparison
driver validates the compact result's mass, N_IC, coupling and tension card,
and records hashes of both its source archives and compact-state input.

Release: the requested project-specific connection is implemented as a
declared phenomenological model. It is suitable for sensitivity studies of
these assumptions. Physical with/without-recoil probabilities still need
a representative ensemble, an appropriate physical control and tests of the
central reduction, hard gate, orbital/spin loss and capture closures.
No new AF run was performed or submitted.

The full 96-path result is saved losslessly in
`data/repeated_encounter_comparison_2026-11-11.json.gz`. Decompress it to read
the model cards, action ledgers and closure scenarios. The uncompressed
JSON SHA-256 is
`2a85b44f3348bdf6fa1bceb0a487c84e27f15ac4d7c24d4e51142544a22b0b47`;
both input-file hashes are also retained inside the result itself.

```bash
python scripts/compare_repeated_encounters.py triplet_return_window.tar.gz \
  --output repeated_encounter_comparison.json
PYTHONPATH=. python -m unittest discover -s research/quirkpheno_recoil -p test_repeated_encounters.py
```
