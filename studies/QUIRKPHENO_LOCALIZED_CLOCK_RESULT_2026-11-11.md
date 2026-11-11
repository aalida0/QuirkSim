# Localized-clock AF result — 2026-11-11

Verdict: PASS_SELECTED_FIXED_THRESHOLD_TRIPLET_PILOT. All ten cases compile,
complete 130 ps, preserve the requested parent mass, and pass local handoff,
kinematic and event-clock checks. This is not a stochastic lifetime ensemble,
full worldsheet validation or annihilation/survival result.

Input quirkpheno_clock.tar.gz SHA-256:
3f08df392f1044f4101feef1826c30ff60889ad15a60c194e7711d8e0da37ccf.
Upstream d91cb27ed3786bc63900faad892513841b5b6ed4; adapter 66eef0b1.
Actual parent mass 225.1571702763212 GeV, neutral 225 GeV, pion 139.5701 MeV
from Geant4. Earlier UI rounding is corrected.

## Clock and controls

Six localized weak cases have |proper time - 50 ps| <=5.571e-12 ps, comfortably
within the computed-clock tolerance. Endpoint-0 lab times are 54.44266345716183
and 54.44266345716375 ps at levels 7/8; endpoint-1 times differ by <5e-12 ps.
The same-input legacy weak control overshoots by 0.00436336 ps and decays at
54.44757553619 ps. The new clock removes that selected event-step bias.
The localized and legacy stable trajectories coincide at every compared sample
(maximum difference zero for both endpoints), supporting that replay sampling
has not changed this stable geometry. No continuum clock error bound is inferred.
All local four-momentum component residuals are below 6e-11 MeV. Neutral charge,
mass, endpoint identity and history-handoff checks pass, including negative PDG.

## Trajectory refinement and directional effect

| Mode, worst endpoint | Maximum level-7/8 position difference (nm) |
|---|---:|
| Stable | 0.0219840 |
| Pion +y, endpoint 0 decays | 0.0218572 |
| Pion -y, endpoint 0 decays | 0.0218572 |
| Pion +y, endpoint 1 decays | 0.0218591 |

These are below 0.1 confinement reference lengths (reference 0.254944419 nm)
over the 13001 sampled common times, using piecewise Hermite interpolation.
No interpolation error certificate or general convergence order is established.

At 130 ps, reversing the endpoint-0 pion rest direction gives y differences:
- decaying endpoint: 12.1246217401 micrometres;
- partner: 0.70379608156 micrometres.

The level-7/8 angular differences agree much more closely than the individual
trajectories, a cancellation that must not be advertised as an absolute error
bound. Pre-decay parent states/times match within each angular pair. No sampled
partner difference appears outside the conservative Euclidean future cone.
This supports directional recoil sensitivity with selected retarded behavior,
not a physical survival probability or annihilation suppression.

Next: small complete fixed-seed sample of independent proper lifetimes and
isotropic pion directions at both endpoints, before inclusive channels. Keep
it conditional on the pion channel, fixed initial kinematics and a finite
130 ps window. Doublet remains paused; no simulation was launched by this audit.
