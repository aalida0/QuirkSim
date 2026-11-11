# Precision confirmation and weak-recoil handoff — 2026-11-11

Status: precision suite prepared; no new AF result or Geant4 decay implementation
is claimed. Doublet paused. This supplements the 130 ps refinement result.

## Precision confirmation

Use --suite precision: eight sequential stable controls, Bz=0/2 T,
levels 5,6,7,8 at 130 ps, same initial state and upstream pinned revision.

| Level | maxBoost | maxMerge (mm) |
|---|---|---|
| 5 | 0.0003125 | 3.125e-10 |
| 6 | 0.00015625 | 1.5625e-10 |
| 7 | 0.000078125 | 7.8125e-11 |
| 8 | 0.0000390625 | 3.90625e-11 |

Per-case 120 s plus 5 s kill grace; one CPU, direct node execution,
no Condor, no set -e/-u. The build retains its separate 600 s limit.
The prior approx factor-two decrease predicts a level-7/8 difference near
0.022 nm (0.086 reference lengths). This is a forecast, not a certified
error estimate. Review new data before advancing the physics claim.
The existing generalized ROOT audit already handles these levels.

Local verification: generated all eight macro/card pairs; checked exact controls,
130 ps unit conversion, level parameters and disabled crossing restart.
Python compilation passes. C++ simulation cannot run in this local environment.

## Recoil extension: source-grounded implementation constraints

Inspected pinned upstream:
- QuirkStringStore::stringFor selects the endpoint by PDG sign.
- QuirkTransportation::EndTracking marks every nonsuspended track's string
  endpoint killed. It cannot treat a weak-decay parent as an ordinary dead end.
- QuirkInfracolorForce::StartTracking treats an initialized first-step track
  as stale state and can clear both endpoint histories.
- QuirkWatcher also has orphan/death handling.
- Particle processes are explicitly replaced in QuirkPhysics.

A neutral daughter therefore needs the same endpoint identity, an explicit
one-time handoff token, and a neutral Quirk particle definition with transport.
A generic Geant4 decay table alone is insufficient.

Proposed narrow implementation sequence (not implemented):
1. Register charged/neutral Quirk definitions with consistent particle/antiparticle
   signs and physical split masses. Map both to the same endpoint slot.
2. Add a local weak-decay process with sampled proper-time threshold. Limit the
   step to decay, record actual time and four-momenta, and transfer endpoint
   ownership to exactly one daughter track.
3. During handoff, suppress parent endpoint-death handling in transportation
   and watcher. On daughter start, consume the token and preserve incoming
   string vectors, borrowed/first-string state, coupling and partner state.
   Ordinary starts/resumptions and genuine track kills retain normal behavior.
4. Update cached endpoint velocity/kinetic state to the daughter. The next
   outgoing string segment uses its new velocity; do not edit the partner
   state instantaneously or clear/restart the string.
5. Use the existing physical decay conventions in quirksim/decays.py as the
   reference (sample_decay, pion_decay_rest, leptonic_decay_rest, boost).
   Test local on-shell daughter/products and four-momentum conservation.
   Lifetime/channel inputs must match the selected triplet model.
6. Log complete string/endpoint state across handoff, proper clock, daughter
   identity and emitted products. Test one and two decays, unchanged partner
   before causal arrival, no extra string creation, and nondecay controls.
   Test both charge signs and both endpoint slots.

Do not substitute zero dynamic charge on a charged particle definition as
a validated neutral species: attached material processes and particle identity
also need to be consistent. A diagnostic no-recoil control must be labeled
as such rather than presented as a physical two-body pion decay.

The adapter can be developed while precision data is pending, but cannot be
called compiled or validated until executed in Geant4. Absolute annihilation
and radiation conclusions remain outside these transport gates.

## Direct AF compute-node command

```bash
cd ~/QuirkSim
git pull --ff-only origin research/variational-worldsheet
batch=/data/aalida/QuirkSim/quirkpheno_precision_$(date +%Y%m%dT%H%M%S)
python3 scripts/af_quirkpheno_transport.py --suite precision --output-dir "$batch"
tar -czf ~/quirkpheno_precision.tar.gz -C /data/aalida/QuirkSim "$(basename "$batch")"
```
