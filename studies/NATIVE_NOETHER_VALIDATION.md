# Direct endpoint + string conservation validation

Status: implementation and standalone invariant checks complete; native Geant4
transport conservation remains UNKNOWN until the AF diagnostic is returned.
The earlier 40/40 inclusive audit establishes vertices/history handoffs only.
It did not save string queues, so full field charges cannot be recovered from
that archive or from endpoint momenta alone. No survival/capture probability
is inferred here.

## Charge definition and analytic checks

Use the pinned public QuirkPheno backend
`tvami/QuirkPheno@d91cb27ed3786bc63900faad892513841b5b6ed4`, specifically
`QuirkHyperbolaStep.cc` and `QuirkInfracolorForce.cc`. The incoming queue at
endpoint i consists of ordered four-displacements s_ik. All spacetime
coordinates below use X=(ct,x,y,z) in mm; P=(E,p) in MeV; T in MeV/mm.

Let S_i=sum_k s_ik. The two future-directed queue paths must join:

    C = X_0 + S_0 = X_1 + S_1.
    P_string = T (S_0 + S_1).
    M_string = T sum_i sum_k (X_i + sum_{l<k} s_il) wedge s_ik.

Add endpoint P_i and X_i wedge P_i. This is a causal-contour charge evaluated
on the actual two asynchronous endpoint events, not a sum of endpoints at
unmatched times with an assumed straight instantaneous string. The queue
joining residual is measured separately. No finite closure defect is
silently projected away. For an exactly static straight string length L,
this gives E_string=T L and M^{01}=-T L^2/2, agreeing with a direct uniform
stress-energy integral. Generic queue masses/merges are retained; merge
artifacts are tested by refinement, not declared absent.

The backend consumes incoming a, appends reflected b to the partner, and
updates X by a+b and P by T(a-b). Thus the change of T(S_0+S_1) cancels
endpoint momentum exactly. Its reflection law implies

    (a+b) wedge P = T a wedge b,

which cancels the angular-tensor update, including partial queue consumption.
The standalone test exercises 20,000 alternating reflections on nontrivial
curved queues and checks all tensor components and contour closure. It
also checks the independently integrated static-string reference, pointlike
emission flux, and the failure caused by omitting a queue. These are analytic
and internal numerical checks, not external certification of Geant4.

At a weak vertex add each emitted product's P and X_vertex wedge P to a
removed-flux ledger. The surviving string history is unchanged. Emitted
particles are accounted at their birth vertices and removed from the audited
subsystem; their later trajectories and interactions are not included.

At B=0, test full four-momentum and all six classical orbital Lorentz
generators. At uniform Bz, only energy and torque-corrected Jz are gated:

    Delta Jz_external = - (q/e) 0.299792458 Bz Delta(r_perp^2)/2 [MeV mm].

This follows directly from integrating r cross q(v cross B); no approximate
orbit fitting is needed. The end-point charge is used only until it
neutralizes. Full momentum and transverse angular components are not
conserved in the imposed external field. Spin and quantum occupations are
outside this classical check. Boosts along Bz induce no electric field.

## Initialization and numerical scope

The backend uses a finite initial null seed and a borrowed-string debt.
Snapshots explicitly skip first-step/debt states, rather than treating those
as a physical field. Production four-momentum and angular tensor remain the
fixed references: later snapshots are not rebased to hide drift. The finite
seed can leave a contour/charge defect; it is measured and must pass the
budget. The two numerical settings reduce both transport tolerances and
seed length (1e-10 mm and 5e-11 mm). This is a regulator refinement, not a
claim of exact quantum-scale initialization. No O(hbar) validation follows.

Declared maximum sampled residual budgets:

| Quantity | Budget |
|---|---:|
| Total energy | 1e-3 MeV |
| Each total momentum component, B=0 | 1e-3 MeV |
| Each orbital tensor component at B=0; corrected Jz at Bz | 1e-3 MeV mm |
| Each queue-contour joining component | 1e-8 mm |
| Discontinuity across a weak vertex, native charge components | 1e-6 |

Snapshots are saved every 2000 accepted endpoint stepping callbacks, before
and after each weak vertex, and at event end. A pass means those sampled
charges meet the declared budgets; it is not a rigorous bound between samples.
A numerical fail is retained and investigated, never repaired by inferring
missing string energy from the expected conserved total. Missing/incomplete
records, nonfinite values, reversed steps, changed production references,
non-antisymmetric tensors and premature stops are rejected.

## AF run

The direct worker campaign has 16 independent cases: stable parent and the
three forced weak channels, isotropic directions, B=0/2 T, levels 7/8;
m_neutral=225 GeV, Lambda=4642 eV, the same boosted initial kinematics as the
40-case suite. It stops at 150 ps, after the forced 50/100 ps proper thresholds
and initial post-neutral returns. It does not measure physical branching
fractions or lifetimes. No radiation/capture model is enabled. Diagnostics
require one event per process in vacuum and a declared uniform Bz. The
runner supplies the matching card and environment; other geometries are
outside this contract.

```bash
cd ~/QuirkSim
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
out=/data/aalida/QuirkSim/noether_$(date -u +%Y%m%dT%H%M%SZ)
nohup python3 -u scripts/af_validate_native_noether.py \
  --output-dir "$out" --workers 4 >"${out}.driver.log" 2>&1 &
echo "$out"
```

Run on a compute node with four CPUs. No Condor submission. The build is
shared; each native process uses one worker. ROOT event files remain on AF,
but ROOT step recording is disabled to avoid another large transport archive.
The driver produces `<out>_summary.tar.gz` with the audit, input manifest,
macros/cards/plans and compact raw Noether snapshots. Failed cases are packed
as UNKNOWN/FAIL; no missing files are promoted to a pass.

Packing recovery (does not rerun transport):

```bash
python3 scripts/af_validate_native_noether.py --output-dir "$out" --audit-only
```

After this conservation gate passes, the next physics task is a sampled
lifetime/branching ensemble with matched direction controls. Radiation,
capture and annihilation still require declared and varied physical models;
a conservation pass alone does not justify a unique annihilation prediction.
