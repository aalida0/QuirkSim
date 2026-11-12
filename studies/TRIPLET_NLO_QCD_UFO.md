# Dirac triplet electroweak production UFO at NLO QCD

## Model definition and what "complete" means here

The physical Lagrangian used for the **production** model is

```
L = L_SM + bar(Q)(i gamma_mu D_EW^mu - M) Q
Q = (Q1+, Q0, Q2-)^T, Dirac (3,0) under SU(2)_L x U(1)_Y,
singlet under SU(3)_c, with a spectator hidden color index.
```

The vectorlike EW vertices are all eight `bar(Q) gamma_mu Q A/Z/W`
conjugate/charged-current structures:

- `A Q1bar Q1`: +e and `A Q2bar Q2`: -e;
- `Z Q1bar Q1`: +e cot(theta_W), `Z Q2bar Q2`: -e cot(theta_W);
- four charged-current W+ and W- transitions with full-strength
  `g_2=ee/sw`, both fermion chiralities; no diagonal neutral A/Z vertex.

There are three **distinct Dirac quirk fields and their antiparticles**,
so six particles in the UFO, not six unrelated triplet representations.

### Gauge-invariant mass interpretation

At tree level the SM gauge symmetry forces the common mass `MQC=MQ0=M`.
One-loop EW splitting (~157 MeV at m=225 GeV) is a physical
**pole-mass correction**, not a new renormalizable gauge-breaking
tree mass term. For gauge-consistent NLO QCD hard scattering use
`MQC=MQ0=225 GeV`. In QuirkSim weak decay propagation use the
loop-improved charged mass. Using unequal masses directly in the
tree-level UFO is only a small phenomenological approximation and
does not preserve the exact SU(2)_L Ward identities.

### Which NLO calculation is covered?

The triplet is an SU(3)_c singlet. At NLO QCD in Drell--Yan production,
gluon virtual loops and real emission couple only to initial-state SM
quarks/gluons. The renormalized SM QCD UV and rational R2 counterterms
are therefore sufficient, without introducing a spurious
`QCD`-charged triplet. The generated production model uses MG5
`models/loop_sm` and inherits those CT files **without modification**.

This provides the relevant *NLO QCD* order `O(alpha_s alpha_EW^2)`
for the cross section. It does not include NLO electroweak corrections
`O(alpha_EW^3)`, infracolor gauge loops, hidden gluon radiation,
confined strings or weak pion decays. The latter are deliberately
applied by the independent QuirkSim propagation stage.

Because the hard amplitude is independent of a spectator hidden-color
index, `sigma_NIC = N_IC * sigma_color_singlet` **at LO and NLO QCD**,
provided the `N_IC` hidden copies are degenerate, have identical EW
couplings and are not subject to hidden-radiation effects at production.

This is **not a full dynamical infracolor or all-orders NLO model**.

## Reproducible UFO build

Use a local MadGraph clone which has `models/loop_sm`
(e.g. 3.8.1, and check your local 3.5.1 installation):

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
MG5="$HOME/codes/Herwig/herwigpp/opt/MG5_aMC_v3_5_1"
test -f "$MG5/models/loop_sm/CT_vertices.py"
python3 scripts/build_dirac_triplet_nlo_ufo.py --mg5-dir "$MG5" --force
python3 scripts/validate_dirac_triplet_nlo_ufo.py --mg5-dir "$MG5"
```

Output: `models/DiracTripletEW_NLO_UFO/`.

The builder preserves the original `loop_sm` CT sector and adds
only two external masses, six quirk/antiquirk particles, five electroweak
couplings, and eight vectorlike vertices. It also fixes all inherited
`restrict_*.dat` cards to include the triplet masses and selects
the massless 5-flavour SM restriction by default for NLO PDFs.
The generated `QUIRKSIM_NLO_PROVENANCE.json` records the MG5 source SHA.

The NLO UFO is assembled from the base code at build time instead of
checking a copied SM model into QuirkSim's repository.

## NLO process commands

```bash
"$MG5/bin/mg5_aMC" examples/pp_dirac_triplet_13p6TeV_NLO.mg5
```

The card generates six inclusive final states with `[QCD]`,
`QCD=0`, `QED=2`, including both charged-pair and charged--neutral
associated hard production. No tree-level `q0 q0~` Drell--Yan exists.

For a **physical NLO cross section** after process generation, run
an aMC@NLO fixed-order integration with a documented NLO PDF
(e.g. an LHAPDF NNPDF set), `sqrt(s)=13.6 TeV`, central scales,
and seven-point `muR/muF` variations. Match LO to the same PDF,
mass, scale, EW input, and flavour convention before reporting K factors.
Do not reuse the old 2.779 pb LO result as an NLO normalization.
Do not call a mere NLO matrix-element `output` a completed NLO rate.

For true NLO+PS generation, use MC@NLO matching and propagate LHE
weights, including negative weights, correctly; the current
LO Pythia8 JSONL adapter is not yet NLO-weight aware.

## Required physics release gates

1. Python import, six species, eight vector vertices and charge/quirk-number conservation.
2. CT UV/R2 inheritance and a valid restriction/parameter card.
3. MG5 true NLO `[QCD]` generation for `u u~ > q1+ q1-` and
   `u d~ > q1+ q0~` (real+virtual+FKS, not LO only).
4. UV/IR pole cancellation and numerical stable finite virtuals.
5. Inclusive NLO cross sections with correct PDFs, scale variation,
   5F flavour scheme, mass convention, all six channels, and hidden-index
   multiplicity checks.
6. Independent gauge and charge conjugation checks, plus at least one
   external electroweak DY reference calculation for the NLO K factor.

Items 4--6 are **not certified by structural tests alone**.

## Provenance

- MadGraph5_aMC source `models/loop_sm` (standard QCD NLO UV/R2),
  upstream: https://github.com/mg5amcnlo/mg5amcnlo
- MadGraph NLO usage: Alwall et al. (2014), arXiv:1405.0301.
- Vectorlike lepton NLO QCD reference templates:
  https://github.com/prudhvibhattiprolu/VLL-UFOs
- FeynRules UFO NLO examples: https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki/NLOModels

No hidden SU(N_IC) gauge interactions or worldsheet dynamics are
misrepresented as SM QCD corrections.
