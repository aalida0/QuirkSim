# Compact-state triplet annihilation: scoped LO result

For the 225 GeV neutral-mass triplet benchmark, the existing LO split gives
m_charged = 225.157170276 GeV. A **specified pure same-flavor charged 1S0
IC-singlet state** has the following leading nonrelativistic two-boson widths:

Gamma(gamma gamma) = 4 pi N_IC alpha_EM^2 Q^4 |psi(0)|^2 / m^2,

Gamma(g' g') = pi (N_IC^2-1) alpha_IC,hard^2 |psi(0)|^2 / (N_IC m^2).

These are the spin-zero projections of the spin-averaged coefficients in
[Li et al., Table II and Eq. (25)](https://arxiv.org/html/2311.15486v2).
The factor four follows because only one of the four initial spin states
contributes to these leading S-wave two-boson channels. The photon expression
reduces to the LO parapositronium width alpha^5 m_e/2 after removing the
spectator IC color factor. The SU(3) gluon expression gives the standard
eta_Q coefficient 8 pi alpha_s^2/(3 m^2) multiplying |psi(0)|^2. These checks
do not validate formation of such a state in the transported triplet ensemble.

## Numerical benchmark and approximations

We declare N_IC=2, alpha_IC,potential=alpha_IC,hard=0.04 and alpha_EM=alpha(0).
The couplings are fixed LO inputs, with no inferred relation between the
macroscopic tension and Lambda_MS. For this single-channel illustration,
V(r)=-(C_F alpha_IC,potential + alpha_EM)/r and the reduced mass is m/2.
The analytic Coulomb state gives |psi(0)|^2 = 23.5641 GeV^3, a binding energy
of 0.0783036 GeV and a Bohr momentum of 4.19888 GeV.

| Computed channel | LO partial width |
| --- | ---: |
| Two photons | 0.622086 keV |
| Two infragluons | 3.50462 keV |

hbar divided by their sum is 1.595e-19 s. This quantity is not reported as
the physical total lifetime or used to assign branching fractions.

The expectation of the linear potential at Lambda_tension=774 eV is only
2.73e-12 of this Coulomb binding energy. This check establishes the smallness
of that particular omitted potential term. It does not control the omitted
running coupling, finite-range EW potential, coupled-flavor mixing,
relativistic or higher-order corrections. No numerical uncertainty band is
claimed. Coherent flavor components can interfere in annihilation; the
same-flavor coefficients must not be summed as incoherent flavor probabilities.

## Release decision

The hard partial-width function and its normalization checks are implemented.
The single-channel numerical illustration is approximate. Its formation
probability and the physical annihilation probability of the original
triplet transport benchmark remain UNKNOWN. No compact state is inserted
into a classical trajectory or assigned a population by this script.

Reproduce with:

```bash
python scripts/triplet_compact_annihilation.py --output compact_annihilation.json
PYTHONPATH=. python -m unittest discover -s research/quirkpheno_recoil -p test_compact_annihilation.py
```

The first-return AF transport study can now extend the existing paired
pilot. Run from the repository on the node; the cases are sequential and use
no Condor submission:

```bash
python scripts/af_quirkpheno_weak_pilot.py \
  --output-dir /data/aalida/QuirkSim/triplet_first_return_260ps \
  --sampled-pion-seeds 1 --duration-ps 260 \
  --case-timeout-seconds 1800 --audit-close-return
```

Use a fresh output directory; `uproot` must be installed in this Python
environment. This prepares six cases (two stable levels and two paired
directions at two levels for seed 1). The compact manifest contains the
post-run close-return classifications. Return that manifest and failed logs
first; retain the ROOT files on AF for later inspection. A censored interval
does not by itself establish a contact or annihilation. No AF execution of
this extended pilot has been performed in this calculation.
