# Radiation and annihilation: source-matched scope

Date: 2026-11-10. Model in the simulator: color-neutral vectorlike Dirac
doublet/triplet, charged-to-neutral weak transitions; string tension is
sigma=Lambda_tension^2. Neither NIC nor its hard coupling is fixed by those
electroweak representations. The conditional study varies NIC=2,3 and
explicitly supplies alpha_IC=0.04. This is not a fit of Lambda_MS to the
macroscopic string tension.

| Ingredient | Inspected primary source | Transfer and implementation |
|---|---|---|
| Glueball energy loss | Li, Liao, Ni, Pei, PRD109 (2024) 095005, III.A, Eq.(7) | Published epsilon times an assumed emission energy per oscillation is a phenomenological input. Implemented as a declared average loss, with period convention explicit. It does not supply a local momentum/angle/recoil kernel. |
| Hard same-flavor Dirac annihilation | Same paper, Table II and Eq.(25), PDF page 7 | IC-singlet, spin-averaged QED and IC threshold coefficients implemented. Electron coefficient extended to other photon-coupled fermions with threshold factors. No complete doublet/triplet EW total is claimed. |
| Excited-state annihilation | Kang and Luty, 0805.4642v3, Eqs.(5.21), (5.28), (5.29), PDF pp.24–26 | WKB short-distance loss requires the quantum s-wave overlap and the hard coefficient. The orbital model accepts the overlap explicitly; geometric proximity does not determine it. |
| Independent treatment of uncertain glueball radiation | Forsyth et al., 2504.02940v2, section 2.1, PDF pp.5–7 | Scalar-model study discusses uncertain glueball deexcitation and chooses to omit it. It supplies no calibrated Dirac-doublet/triplet emission kernel. Its scalar hard widths are not reused. |
| Smooth electromagnetic radiation | PDG, Accelerator physics of colliders (2002), radiation section; rationalized SI-to-natural-unit conversion | Implemented single-charge Lienard power and checked against its angular integral and the constant-longitudinal-force limit. Smooth radiation excludes instantaneous weak-vertex bremsstrahlung and pair interference. |

Links:

- https://doi.org/10.1103/PhysRevD.109.095005
- https://arxiv.org/pdf/0805.4642
- https://arxiv.org/pdf/2504.02940
- https://pdg.lbl.gov/2002/accelrpp.pdf

## Normalization checks and conventions

The photon coefficient in Li et al. Table II reproduces the textbook
parapositronium width after projecting the spin average onto its singlet.
The IC/photon ratio is independently checked by explicit SU(N) generator
traces, including the singlet initial color projection. The tests do not
merely compare two copies of the same scalar formula.

Kang-Luty's displayed two-boson widths and later compact per-crossing
expressions do not reproduce all the Table-II normalizations by naive
substitution. In particular, the displayed IC width in Eq.(5.35), combined
with Eq.(5.29), differs by a factor two from Eq.(5.39); Table II of the
2024 paper has a further factor-two difference from the displayed 2009
width. This comparison was checked on rendered PDF equations. We do not
mix those compact formulas with the selected hard coefficients. The
WKB flux relation and selected hard kernels are recorded separately.
This is a convention/normalization discrepancy register, not a claim
that every source uses identical state/velocity definitions.

In rationalized EM with alpha=e^2/(4 pi), a single charge has
P=(2 alpha/3) hbar gamma^6 [beta_dot^2-|beta cross beta_dot|^2]
in GeV/ps. Some inspected quirk estimates display an 8 pi alpha/3
coefficient. We retain our explicitly derived convention and independent
angular-integral check rather than importing a cooling-time normalization
without matching its electromagnetic definitions.

One full signed radial period is 4p/f; successive zero-length radial
crossings are separated by 2p/f, with f=sigma/hbar. The conditional
glueball input is specified per **full signed period**. A probability
per encounter is not silently substituted for it.

The average loss epsilon*Lambda does not itself describe an on-shell
glueball. Its mass spectrum and relation to sqrt(sigma) require a scale
match; an actual radiation event must respect its mass threshold and
four-momentum balance. The current worldsheet never emits an arbitrary
on-shell glueball from an isolated on-shell endpoint of the same mass:
that process is kinematically forbidden without energy/momentum transfer
from the string or another degree of freedom.

## Reuse decision

The selected threshold kernels and conditional master equations are
reviewable calculations. The radiation transition/recoil law, density
matrix at the classical-to-quantum matching surface, complete EW
annihilation matrix, and quantum capture/return map remain missing
physical inputs. A finite-level master-equation solver is implemented
but its physical rate matrix is not supplied by a classical NG trajectory.

The literature search here was scoped to these formula sources and their
model applicability. It is not a new exhaustive novelty or citation-graph
audit, and it does not certify a detector sensitivity.
