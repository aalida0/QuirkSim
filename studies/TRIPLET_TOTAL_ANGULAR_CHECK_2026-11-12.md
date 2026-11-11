# Endpoint versus pair-plus-string angular momentum

## Result

The available conservation checks do not reveal a large string contribution
cancelling the fitted recoil angular momentum. This is a conservation-target
comparison, not a direct numerical integration of the string stress tensor.

At zero field, reconstructing the angular tensor carried away by the two
pions determines the full surviving pair-plus-string target. Boosting that
target to the surviving system rest frame gives:

| Isotropic seed, level 8 | Fitted relative L / hbar | Full rest-frame target / hbar |
|---:|---:|---:|
| 1 | 2.05387696e12 | 2.05389895e12 |
| 4 | 4.44376924e11 | 4.44395291e11 |
| 7 | 2.70968616e12 | 2.70972993e12 |

The largest relative discrepancy of these magnitudes is 4.133e-5 (0.00414%).
These large angular scales therefore also occur in the full conservation
target; they are not merely an extrinsic center-of-mass term in the endpoint
fit. The equality is not asserted at the quantum one-hbar scale.

For all 33 level-8 neutral-qualified cases at 2 T, the axial total target
agrees with the endpoint axial angular momentum within 0.151%. Its absolute
value is at least 1.82e10 hbar. Only this component is reconstructed at 2 T;
the other components are not conserved by the external solenoid. These
fixed-production-frame axial targets are not silently identified with
rest-frame quantum partial waves.

## Distinct conservation method

Production is at the origin with zero classical orbital angular momentum.
At B=0 the surviving system has

    M_pair+string = -sum_vertices(x_pion wedge p_pion).

Both x and p are boosted to the fixed production CM, including emission
times in x^0=ct. The tensor is then boosted to the surviving system rest
frame using initial pair four-momentum minus pion four-momenta. The tiny
initial numerical string seed is omitted from this energy convention.

The interpolated endpoints independently give

    J_end = r_relative cross p_relative + R_midpoint cross P_end.

The second, extrinsic term is explicitly retained. Subtracting J_end from
the full target leaves a required string contribution OR numerical error;
it does not measure either separately. For zero-field isotropic fine cases,
the norm of this vector residual is at most 0.0402% of the fixed-frame target.

At uniform Bz, the integrated external axial torque before the endpoints
neutralize can be evaluated without trajectory integration:

    Delta Jz_ext = -sum_endpoints(q*Bz/2)*r_perp(decay)^2
    Jz_pair+string = Delta Jz_ext - sum_vertices(r cross p_pion)_z.

The conversion for q=+e is 0.000299792458 GeV/mm/T, followed by division by
hbar*c in GeV mm. Initial positions vanish. Removing each pion at its weak
vertex avoids including its subsequent magnetic motion in the subsystem.
The frame boost is along B, so it induces no electric field. The pilot field
is uniform over the reachable region: the card spans 100 cm radially and
longitudinally, beyond the light-travel distance of these sub-ns histories.

This axial flux identity follows from Lorentz torque; the canonical-angular
momentum relation is also given in Groening, Xiao and Chung,
[arXiv:1709.09538v1, Eq.(3)](https://arxiv.org/pdf/1709.09538).
Only that single-particle identity is reused; their paraxial beam extension
is not needed. Extension to the pair includes its internal string stresses.

## Code and numerical evidence

The pinned QuirkPheno backend already has `GetAngMomentum`, `GetMomentOfE`
and `GetSumStrings` in `quirkcore/src/QuirkInfracolorForce.cc`. They operate
on causal queues, rather than exporting a complete common-time string slice.
The current compact archive contains endpoints and weak vertices, not those
queues or a directly measured field angular tensor. Summing arbitrary
asynchronous queue values would not establish the required closed-system
charge. No direct field measurement is claimed here.

The new script evaluates all 54 neutral-qualified windows including
refinements. Six tests pass: three new tests verify removal of extrinsic
orbital action for a point particle, Lorentz transformation of a rest angular
tensor, and axial torque against an exact helical Lorentz orbit; three
existing angular-budget tests also pass. These are internal mathematical
checks, not an external certification of Geant4 string charges.

## Decision

The available classical conservation targets support the large fitted
angular-momentum scale, especially in the fully reconstructible zero-field
cases. The bounded-emission reachability argument remains conditional on its
emission bounds and on this classical action interpretation. A unique
annihilation probability is still not established. Direct common-slice
string charges and quantum state occupations remain unmeasured; collinear
zero-field numerical residuals of order a few to hundreds of hbar are not
promoted into physical orbital excitation.

Source archive SHA256:
`ca6ff21c5447efa156a65732430be467d27ba993ab8218558e329d6be6773cf1`.

```bash
python scripts/audit_triplet_total_angular.py triplet_ensemble_v1_audit.tar.gz \
  --output total_angular.json
```
