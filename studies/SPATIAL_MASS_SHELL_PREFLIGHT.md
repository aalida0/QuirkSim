# Spatial-impulse mass-shell precision proposal

Status: local counterexample/regression checks pass. Native AF verification of
this change is pending; existing boost-generator failures remain failures.
No conservation tolerance or production reference is changed.

## Mechanism and scope

In the pinned public QuirkPheno analytic string step, incoming a and reflected b
update endpoint P by T(a-b) and X by a+b. Without electromagnetic action, the
exact reflection preserves P^2=m^2. The legacy code subsequently rescales *all*
four momentum components by m/sqrt(P^2). Repeating that floating-point scale
factor can add a biased spatial impulse at each step. The outgoing string
vectors do not include this additional numerical impulse.

The proposed opt-in change retains the computed spatial string impulse and
sets E=hypot(m,|p|) when no electromagnetic kick acts. It preserves the legacy
projection if a nonzero external field acts on a charged endpoint. A neutral
endpoint has no EM kick even in a nonzero imposed field. Material/other processes
are not modified. This changes the numerical mass-shell projection; it does not
add a physical radiation/capture force or correct measured charges to their
expected targets. The kinetic trajectory and Noether diagnostic both see the
same updated particle state.

QP_SPATIAL_MASS_SHELL=1 enables the proposed method. Other native runners keep
their legacy default. The physical-ensemble runner records and enables it for
this new preflight. Exact string-reflection momentum balance still has to pass
the native AF audit after the Geant4 energy/direction round-trip.

## Local evidence and limits

The standalone native_shell_check.cc reproduces the original projection and a
Geant4-style energy/direction momentum round-trip over two million alternating
null-segment impulses. This is a controlled local numerical counterexample, not
a full native string run or an external physics benchmark. The legacy maximum
accumulated four-momentum residual is about 7.4e-5 MeV; the proposed projection
gives about 4.8e-8 MeV for that fixture. The test requires at least 100-fold
improvement and a new residual below 1e-6 MeV. Configuration parsing and a
bounded two-case manifest are checked separately. Field-kick integration is
unchanged; passing the no-kick test does not certify electromagnetic transport.

A separate analytic diagnostic explains the boost-error amplification: a small
momentum nonconservation source contributes X wedge d(delta P) to orbital
Lorentz-generator drift. In a predominantly longitudinally boosted history,
M03 approximately accumulates integral ct d(delta Pz - beta delta E).
This distinguishes a momentum-source contribution from an independent clock
error. It is a diagnostic approximation, not a subtraction applied to the audit.
All original energy, momentum, rotational and boost-generator budgets remain.
The summary now reports the individual conserved tensor-component maxima.

## Two-case AF rerun only

On a compute node with two available CPUs:

```bash
cd ~/QuirkSim
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
out=/data/aalida/QuirkSim/shell_preflight_$(date -u +%Y%m%dT%H%M%SZ)
nohup python3 -u scripts/af_physical_weak_ensemble.py \
  --output-dir "$out" --workers 2 --preflight-only >"${out}.driver.log" 2>&1 &
echo "$out"
```

This runs only the same seed-1, zero-field recoil event at levels 7/8 to 1000 ps.
It does not release the 56 remaining population/control runs. Return
<out>_summary.tar.gz, including failures. ROOT steps remain disabled. Repacking
uses --audit-only with the same output directory; it does not rerun transport.
Only after both native conservation cases pass should the larger pilot resume.
The retained first-return geometry is still classical and does not establish
annihilation or quantum capture.
