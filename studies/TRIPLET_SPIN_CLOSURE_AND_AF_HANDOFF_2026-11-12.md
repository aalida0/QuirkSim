# Spin-sector calculation, population accounting and remaining AF execution

## New physical subcalculation

The missing specified compact spin-one contribution has been calculated at
LO in the IC-singlet two-vector sector. It uses derivatives of the complete
t/u amplitude at zero incoming momentum, with explicit orbital x spin
projectors for 3P0, 3P1 and 3P2. It does not turn a finite-velocity cross
section into an S-wave width. A coupled three-channel L=1 radial problem is
solved in the same fixed-coupling potential as the previous spin-zero result.

| Specified lowest compact state | Sum of modeled two-vector widths (GeV) | Proper lifetime from modeled channels (ps) |
|---|---:|---:|
| 1S0, previous calculation | 1.97500e-5 | 3.33271e-8 |
| 3P0 | 3.65728e-10 | 0.00179973 |
| 3P1 | 4.14623e-12 | 0.158749 |
| 3P2 | 1.04946e-10 | 0.00627192 |

The modeled channels are WW, IC gg, gamma gamma, gamma Z and ZZ. Massive
electroweak vectors allow a nonzero 3P1 contribution; the two-identical-
massless-vector channels vanish. These are widths of specified compact
states, not their formation probabilities or an inclusive transport rate.
Additional higher-multiplicity channels and radiative corrections are not
claimed calculated. Mesh stability is not a physical uncertainty estimate.

## Why a neutral 3S1 width is not the default missing channel

Bai and Hill, [1005.0008v2](https://arxiv.org/pdf/1005.0008), Eqs.(1)--(3),
give a generalized charge-conjugation symmetry for vectorlike Y=0 weak
multiplets. Their hidden gauge representation is complex, with N>=3; their
light-pion spectrum is not transferred to this two-color heavy-quark model.
Instead we check the two-color extension explicitly. For

    A = epsilon_IC tensor S_EW,

the six gauge generators satisfy A*T.conjugate()*A.dagger()=-T and
A*A.conjugate()=-I. The symmetric charged-neutral mass splitting also
respects this transformation. Thus quirk charge conjugation combined with
these internal transformations can leave both sets of gauge fields fixed.
The transformation acts as a Z2 on mesons, while its square on a fundamental
fermion is minus one. This is a derivation for the minimal gauge Lagrangian,
not a reuse of an unrelated lattice spectrum. The pseudoreal two-color
structure is also discussed in [1809.09117](https://arxiv.org/pdf/1809.09117).

For a neutral particle-antiparticle partial wave, G=(-1)^(L+S). The SM vector
production current is G-even; SM emissions and pure two-color glue states
are G-even. Therefore the lowest neutral G-even compact sectors include
1S0 and 3P_J, whereas a neutral 3S1 is G-odd. A single-pair G-even source
with G-even emissions cannot populate that 3S1 state. This statement assumes
Y=0, vectorlike couplings and no G-breaking higher-dimensional interactions.
It is not a claim that arbitrary externally prepared spin-one states decay.

The AF gun records individual classical charge/history branches and no
coherent production amplitude. Applying this selection at quantum matching
is an explicit physical source model, not information recovered from ROOT.
The original statement that a heavy vector current has no initial 1S0
population remains valid; the P-wave calculation supplies a spin-preserving
annihilation channel without imposing an artificial 1/4 singlet reset.

## Validation

In the massless abelian limit the derivative calculation gives
Gamma(3P0)=9*alpha^2*|Rprime(0)|^2/m^4 and
Gamma(3P2)=(12/5)*alpha^2*|Rprime(0)|^2/m^4, with Gamma(3P1)=0.
Multiplication by N_c=3 reproduces the Born quarkonium normalization and
the 4/15 ratio given by Sang et al.,
[1511.06288v2](https://arxiv.org/pdf/1511.06288), Eqs.(3)--(8).
Only the Born hard reference is reused; QCD corrections are not transferred.

Hard matrices change by less than 2e-8 relative when both angular quadrature
and the momentum-derivative step are refined. The last 8000-to-16000-point
radial width changes are below 0.027%. Equality of ordered t/u P projections
checks the charged/neutral WW amplitude factors (.5,1,.5); the vector
s-channel has no linear-momentum P-wave projection. Kernel positivity,
J projectors, symmetry algebra, existing S-wave/finite-velocity references
and population/packaging checks are tested. All this validates the stated
LO calculation, not a nonperturbative formation law or an AF quantum state.

## Capture-to-annihilation model is now executable without a fixed spin fraction

Within the declared compact support {lowest 1S0, lowest 3P0,3P1,3P2}, any
population mixture has its post-capture annihilation probability between
C*(1-exp(-Gamma_min*t/hbar)) and C*(1-exp(-Gamma_max*t/hbar)). No spin
occupation needs to be invented. The slowest modeled lifetime is 0.159 ps,
so at macroscopic post-capture ages these endpoints both approach C.

The new output evaluates matching fractions eta=0,.1,1, where C=eta times
the earlier eligible-capture probability. Eta is a compact-subspace
occupation scenario, not a per-encounter transition rate. The eta=1 ideal
absorbing prescription therefore closes the conditional model; it does not
prove all physical captures immediately reach the lowest radial states.
Excited radial populations, transitions and off-gate radiation remain
explicit limitations. The earlier broad spin interval remains appropriate
outside this specified compact support.

The existing recoil fits remain capture-blocked in the selected family.
The zero-field collinear controls reach the efficient boundary at about
.151 seconds and then annihilate rapidly in the eta=1 compact-subspace
model. None of these statements establishes a detector-time or inclusive
annihilation contrast. The one-nanosecond result applies after the fitted
both-neutral starting point, not to the whole event since production.

## Population and branching accounting

All eight isotropic seeds per point remain in the accounting. The original
three/five both-neutral fits are not treated as the whole event population.
For the remaining zero/one-decay categories, the code retains outcome
intervals [0,1]. Exact binomial intervals quantify the very small fixed-point
pilot; points share seeds and are not pooled as independent observations.

The pion branching fraction is .9760480104 per endpoint. Hence a fraction
1-B_pi^2=.04733028135 of all two-decay channel assignments contains at least
one leptonic decay. Without claiming a leptonic transport calculation, an
interval [lo,hi] in the forced-pion ensemble becomes
[B_pi^2*lo, B_pi^2*hi+1-B_pi^2]. This is a missing-channel envelope, not a
statistical confidence interval. With the presently unresolved histories,
the all-category pilot upper endpoints remain .643 (benchmark/tilted/0T)
or .405 (faster). They must not be reported as measured physical rates.

## AF continuation: prepared, not executed here

The wrapper can now run independent one-CPU cases concurrently. A 32-seed,
four-point suite has 404 transport cases, including stable/refinement
controls. Four-nanosecond macros and hashes were checked in preparation
mode. No Condor job was submitted. Native execution of this new parallel
campaign has not been tested here because AF/CVMFS execution is unavailable.

On an AF node with four CPUs available:

```bash
cd ~/QuirkSim
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
out=/data/aalida/QuirkSim/triplet_late32_$(date -u +%Y%m%dT%H%M%SZ)
python scripts/af_triplet_recoil_ensemble.py --seeds 32 --workers 4 \
  --duration-ps 4000 --case-timeout-seconds 600 --output-dir "$out"
```

The wrapper builds once, then runs four independent processes and exports
compact post-neutral windows, logs and provenance; ROOT files stay on AF.
The standalone post-neutral packager can recover completed runs without
rerunning transport. Failed transport/audit cases stop the procedure. Histories
with no resolved post-neutral minimum remain explicitly censored. Synthetic
ROOT round-trip tests check that the compact archive preserves those counts.

## Remaining-task release decision

| Task | Actual status |
|---|---|
| Lowest G-even compact spin-one two-vector widths | Implemented and checked at declared LO |
| Conditional compact capture-to-annihilation solver | Implemented, including occupation envelopes |
| Zero/one-decay and missing branching accounting | Implemented as explicit bounds; no silent dropping |
| Larger and later-return transport | Prepared; requires AF execution |
| Production-weighted inclusive probability | Not established by these unweighted fixed-point pilots |
| Unique nonperturbative capture/transition law | Not supplied by the inspected literature or classical histories |
| Submission-ready physical-annihilation paper | Not finished; conditional phenomenology is the supported scope |

The software and compact-state calculations are available for review. It
would be incorrect to label every remaining task complete before the AF
sample is run or to replace the remaining physical assumptions with test
counts. The dominant missing physical input is now formation/deexcitation,
not a missing lowest-P-wave spin-one annihilation width.
