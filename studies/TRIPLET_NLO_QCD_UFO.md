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

## Verified NLO QCD numerical pilot (MG5_aMC 3.8.1)

GitHub Actions [run 43255774024](https://github.com/aalida0/QuirkSim/actions/runs/43255774024)
successfully completed a fixed-order NLO integration on November 12, 2026
for **only** `u u~ > q1+ q1-` at 13.6 TeV using the internal pilot PDF settings
and `req_acc_FO=0.05`.

| Diagnostic | Result |
|---|---|
| Born/virtual/real/FKS generation | PASS |
| Pole cancellation | 20/20 phase-space points at tolerance 1e-5 |
| Soft/collinear tests | PASS |
| MadLoop evaluations | 1,030 stable / 1,030 |
| Single-channel NLO hadronic contribution | 0.1757 +/- 0.00093 pb (Monte Carlo error) |
| MadLoop reduction | CutTools, double precision |

**This is not the inclusive pp triplet production cross section.** It omits
other initial flavours and five remaining charge channels; scale/PDF
uncertainties and the spectator `N_IC` multiplicity are not represented in
the quoted pilot number.

### MG5_aMC 3.5.1 compatibility

The older `loop_sm/write_param_card.py` contains the Python-2-only
`need_writing.sort(self.order_param)`. QuirkSim's builder now patches
this in its *copied* model using
`sort(key=__import__('functools').cmp_to_key(self.order_param))`.
The original installed MadGraph files are not modified.

The user confirmed local structural validation with MG5 3.5.1:
6 states, 8 EW vertices, 103 CT vertices and a generated parameter card.
Local MG5_aMC 3.5.1 has now **passed** NLO process generation, compilation,
20/20 UV/IR pole-cancellation checks, and one-channel fixed-order integration.
The local and 3.8.1 pilot rates agree: **0.1757 ± 0.00093 pb** (internal PDF,
only \`u u~ > q1+ q1-\`, no N_IC factor). This does not determine the
six-channel inclusive rate or full NLO PDF/scale uncertainties.

All six charged-pair and charged--neutral NLO processes successfully
generated Born, real-emission, virtual and FKS contributions in CI
[run 43270607758](https://github.com/aalida0/QuirkSim/actions/runs/43270607758).
(The earlier failed CI was due solely to its assertion expecting one aggregate
six-subprocess message rather than six one-subprocess messages.)

## Inclusive six-final-state fixed-order NLO calculation

The command-card builder `scripts/prepare_triplet_nlo_inclusive.py` creates
**two separate MG5 command cards** so first-time MadLoop library installation
cannot consume commands intended for a later launch. It includes five-flavour
incoming partons and all six physical Drell--Yan final states.

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
export MG5="$HOME/codes/Herwig/herwigpp/opt/MG5_aMC_v3_5_1"

# PILOT (diagnostic internal PDFs; choose a new --output-dir if it exists)
python3 scripts/prepare_triplet_nlo_inclusive.py --pdf internal \
  --output-dir results/pp_dirac_triplet_13p6TeV_NLO \
  --card-dir results/triplet_nlo_cards --accuracy 0.02
"$MG5/bin/mg5_aMC" results/triplet_nlo_cards/triplet_six_nlo_generate.mg5 \
  2>&1 | tee results/triplet_six_nlo_generation.log
"$MG5/bin/mg5_aMC" results/triplet_nlo_cards/triplet_six_nlo_integrate.mg5 \
  2>&1 | tee results/triplet_six_nlo_integration.log
```

Do not interpret this internal-PDF result as the final physical total rate.
For final NLO QCD cross sections, install/configure LHAPDF in MG5 first,
choose an existing **NLO** 5F PDF set and central-member ID, and generate a
**different** output directory:

```bash
python3 scripts/prepare_triplet_nlo_inclusive.py \
  --pdf lhapdf --lhaid YOUR_NLO_LHAPDF_ID \
  --output-dir results/pp_dirac_triplet_13p6TeV_NLO_LHAPDF \
  --card-dir results/triplet_nlo_cards_lhapdf \
  --accuracy 0.01
```

The builder does not install PDFs. Verify `lhapdf-config`, actual PDF
membership, its alpha_s(MZ), and central 5F scheme first. Propagate PDF
uncertainty and an appropriate 7-point scale envelope separately. Default MG5
scale summaries may use a **nine-point** envelope rather than the standard
seven-point one.

The inclusive MG5 result still represents **one hidden-color copy**.
If hidden color is a spectator and there are N_IC degenerate copies, use
`sigma_physical = N_IC * sigma_NLO`. Keep this normalization **outside**
the UFO to avoid accidentally double-counting its hidden degrees of freedom.

## Inclusive NLO internal-PDF integration (local MG5 3.5.1)

Local six-channel `pp` NLO run at `sqrt(s)=13.6 TeV` and equal
tree-level quirk masses of 225 GeV completed successfully:

| Quantity | Local result |
|---|---:|
| Inclusive NLO partonic hard-process sum (one spectator hidden-color copy) | **3.746 +/- 0.0055 pb** |
| MG5 reported 9-point scale envelope | **+1.5% / -1.1%** |
| Grid-setup estimate | 3.752 +/- 0.0098 pb |
| UV/IR poles | 16 successful 20/20 phase-space checks listed in local log |
| N_IC=2 spectator-rescaled diagnostic rate | 7.492 +/- 0.011 pb |

These numbers are from MG5's **internal PDF setup** and do NOT give
an externally verified NLO PDF/alpha_s-consistent prediction.
The numerical error is an integration uncertainty, *not* a PDF or
perturbative theory uncertainty. The nine-point envelope is not necessarily
the standard seven-point muR/muF envelope, and its small magnitude must
not be interpreted as complete theoretical accuracy.

Do **not** compare this rate to the prior 2.779 pb LO run as an official
K factor, because the old LO sample used different PDF/scale/mass
settings. A valid K factor requires matched electroweak inputs, scales,
masses, PDFs and flavour scheme.

### Explicit NLO LHAPDF production configuration

For an example modern 5-flavour NLO set use
`NNPDF40_nlo_as_01180` (LHAPDF **331700**, alpha_s(MZ)=0.118).
A user reported multiple existing `lhapdf-config` installations.
Choose a single working executable and verify `--version`, `--libdir`,
`--datadir`, and that the PDF files are installed and visible.

```bash
export LHACFG="$HOME/codes/tools/MG5_aMC/HEPTools/lhapdf6_py3/bin/lhapdf-config"
"$LHACFG" --version
"$LHACFG" --libdir
"$LHACFG" --datadir
# If the chosen data directory does not contain the set, install it
# with that LHAPDF installation's lhapdf manager before continuing.
PDFDIR="$("$LHACFG" --datadir)"
test -f "$PDFDIR/NNPDF40_nlo_as_01180/NNPDF40_nlo_as_01180.info"
```

To configure **the exact MG5 instance used for QuirkSim**, run a
dedicated configuration script before generating the new NLO directory:

```bash
cat > /tmp/quirksim_set_lhapdf.mg5 <<EOF
set lhapdf $LHACFG
quit
EOF
"$MG5/bin/mg5_aMC" /tmp/quirksim_set_lhapdf.mg5
python3 scripts/prepare_triplet_nlo_inclusive.py \
  --pdf lhapdf --lhaid 331700 --accuracy 0.01 \
  --output-dir results/pp_dirac_triplet_13p6TeV_NLO_NNPDF40 \
  --card-dir results/triplet_nlo_cards_NNPDF40
"$MG5/bin/mg5_aMC" results/triplet_nlo_cards_NNPDF40/triplet_six_nlo_generate.mg5 \
  2>&1 | tee results/triplet_six_nlo_NNPDF40_generate.log
"$MG5/bin/mg5_aMC" results/triplet_nlo_cards_NNPDF40/triplet_six_nlo_integrate.mg5 \
  2>&1 | tee results/triplet_six_nlo_NNPDF40_integrate.log
```

Verify that the actual run card records `pdlabel=lhapdf`, `lhaid=331700`,
the correct LHAPDF prefix and the exact PDF set name. Reweight or
rerun to evaluate PDF replica uncertainty and a standard seven-point
scale envelope before release. N_IC should be multiplied **once** outside
the UFO.

## NNPDF4.0 NLO central-rate integration (user local MG5 3.5.1)

On 2026-11-12 the six-final-state inclusive fixed-order NLO run completed
with `pdlabel=lhapdf`, requested set ID `331700`
(`NNPDF40_nlo_as_01180`), LHAPDF 6.5.5, `sqrt(s)=13.6 TeV`,
degenerate triplet masses `M=225 GeV`, one spectator hidden-color copy.

| Quantity | Result |
|---|---:|
| Fixed-order NLO central rate, 1 hidden color | **3.991 ± 0.0057 pb** |
| Grid-setup estimate | 3.997 ± 0.010 pb |
| For `N_IC=2` under spectator-index assumptions | **7.982 ± 0.0114 pb** |
| Pole cancellation | 16 reported checks, each 20/20 points at 1e-5 |
| PDF set chosen in MG5 run card | `pdlabel=lhapdf`, `lhaid=331700` |
| PDF provenance from installed .info | `SetIndex=331700`, `NumMembers=101`, `OrderQCD=1`, `AlphaS_MZ=0.1180024` |
| Compared with internal-PDF NLO central | 3.746 ± 0.0055 pb (NOT a K factor) |

The MG5 log prints `Failed to access python version of LHAPDF`
but then `INFO: Using LHAPDF v6.5.5 interface for PDFs`. This indicates
the compiled LHAPDF interface is available; missing Python bindings
do not invalidate the central fixed-order calculation. Before final
release, cross-check the stored `Cards/run_card.dat` and banner
to confirm `lhaid 331700` was used, and exclude any PDF-loading
error messages.

**This is a central NLO QCD result, not a fully uncertainty-qualified
prediction.** Integration error is not a theory uncertainty. The printed
9-point scale envelope is not necessarily the standard 7-point envelope.
NLO PDF replica uncertainty and a scale-envelope prescription need to be
reported. A true NLO K factor requires an LO run with matched NLO PDFs,
alpha_s, electroweak parameters, scales, final-state masses and flavour
scheme (in addition to the conventional LO-PDF LO baseline as appropriate).

For NLO fixed-order scale/PDF uncertainties, MadGraph's on-the-fly
reweighting options `reweight_scale`, `rw_rscale`, `rw_fscale`
and `reweight_PDF` are documented in its NLO run card. Use the
compiled LHAPDF interface; the separate Python-based systematics
module can require LHAPDF Python bindings.

## NNPDF4.0 replica-uncertainty recovery (November 12, 2026)

The second six-channel NLO run with `reweight_PDF=True` returned

```text
Total cross section: 3.989e+00 +- 5.7e-03 pb
Dynamical_scale_choice -1 (envelope of 9 values): +1.5% -1.1%
PDF variation: 331700 (101 members; combination method unknown): 3.989e+00 pb
```

The central value agrees with the previous 3.991 pb integration.
The absence of a PDF uncertainty does **not** mean 0% PDF error.

**Root cause verified from MG5_aMC 3.5.1 source:**
`madgraph/interface/amcatnlo_run_interface.py`,
`pdf_scale_from_reweighting()`. If
`misc.import_python_lhapdf()` fails, the code marks the PDF uncertainty
type as `unknown` instead of computing `getPDFSet(...).uncertainty(...)`.
The underlying 101 member-weighted fixed-order cross sections are still
written to `SubProcesses/P*/all_G*/scale_pdf_dependence.dat`. This
does not affect fixed-order integration through the compiled LHAPDF
interface.

**Recovery without rebuilding LHAPDF or rerunning MG5:**

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
python3 scripts/extract_triplet_nlo_pdf_replicas.py \
  --process results/pp_dirac_triplet_13p6TeV_NLO_NNPDF40 \
  --run-log results/triplet_nlo_uncertainties.log \
  --nic 2
```

The Python script reconstructs each PDF member cross section by
**exactly** replicating MG5's summation over job directories weighted by
`wgt_frac` from the trusted local `SubProcesses/job_status.pkl`, then
uses the sample standard deviation across the 100 NNPDF replicas for a
68%-CL PDF-only uncertainty. It asserts that the reconstructed central
matches the original MG5 summary within 1%. The 100-replica
population-RMS is printed as a convention cross-check.

This method is restricted to `NNPDF40_nlo_as_01180` and assumes its
standard 101-member central-plus-replicas organization. It does **not**
include an alpha_s uncertainty. A full scale uncertainty should use the
appropriate seven-point muR/muF envelope rather than the default
nine-point envelope printed by MG5.

**Preservation caution:** MG5 overwrites `SubProcesses/job_status.pkl`
and scale/PDF variations when you rerun the same process directory.
Extract the replica statistics now, before performing another integration.

## Confirmed 100-replica PDF RMS and 7-point scale extractor

User's November 12, 2026 local output from 17 saved NLO integration jobs:

| Observable (one hidden-color copy) | Value |
|---|---:|
| NNPDF40 central | **3.9893368 pb** |
| 100-replica mean | 3.9892335 pb |
| 68%-CL PDF uncertainty (sample RMS) | **0.0573033 pb = 1.436%** |
| Population-RMS cross-check | 0.0570161 pb |
| Replica mean minus central | -0.0001033 pb |
| N_IC=2 central | **7.9786736 pb** |
| N_IC=2 PDF-only uncertainty | **0.1146067 pb** |

This replica RMS is a PDF-only uncertainty, and **does not** include
alpha_s variation or the perturbative scale envelope.

To extract a **standard seven-point** scale envelope from the **same**
17 sets of saved per-job reweighting tables without a new NLO run:

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
python3 scripts/extract_triplet_nlo_7point.py \
  --process results/pp_dirac_triplet_13p6TeV_NLO_NNPDF40 \
  --run-log results/triplet_nlo_uncertainties.log \
  --nic 2 \
  --json-out results/triplet_nlo_normalisation.json \
  | tee results/triplet_nlo_7point.txt
```

This script replicates MG5's `wgt_frac`-weighted aggregation of the
`scale_pdf_dependence.dat` files and verifies that the reconstructed
central cross section is consistent with the MG5 log and the independent
PDF-weight central. MadGraph 3.5.1 swaps the scale factor `1.0` to
the **first element** of each `rw_(r/f)scale` list. The reweighting
array is ordered with muR fastest, muF slowest. The seven-point envelope
excludes `(muR,muF)=(0.5,2)` and `(2,0.5)` while the nine-point
envelope is printed as a consistency check against MG5's
`+1.5%/-1.1%` summary.

The script reports the central and all nine points for auditing, the
seven-point absolute/percentage envelope, and the PDF RMS with correct
`N_IC` normalization. The JSON record captures the numbers **before
a new integration overwrites the active MG5 job-weight files**.

### Matched LO baseline after preserving the scale/PDF result

MG5_aMC's NLO output directory supports `order=LO` and
`fixed_order=on`. Only **after** writing the
`triplet_nlo_normalisation.json` file above, run a separate LO
integration using the *same* NLO UFO and PDF set, equal quirk
masses, 13.6 TeV beams and central scale choice:

```bash
cd ~/QuirkSim
PROC="$PWD/results/pp_dirac_triplet_13p6TeV_NLO_NNPDF40"
cat > /tmp/triplet_matched_lo.mg5 <<EOF
launch $PROC
order=LO
fixed_order=on
shower=off
done
set pdlabel lhapdf
set lhaid 331700
set iseed 7482
set req_acc_FO 0.01
set nevents 0
set ebeam1 6800
set ebeam2 6800
set fixed_ren_scale False
set fixed_fac_scale False
set dynamical_scale_choice -1
set reweight_scale False
set reweight_PDF False
done
quit
EOF
"$MG5/bin/mg5_aMC" /tmp/triplet_matched_lo.mg5 \
  2>&1 | tee results/triplet_matched_lo.log
```

If the LO run succeeds, the matched perturbative
`K_NLO/LO = 3.9893368 pb / sigma_LO_same_NLO_PDF` is meaningful
for the chosen **same-NLO-PDF** convention (PDF and scale inputs held
fixed). A separate usual PDF-order-matched LO baseline can be given
with an LO PDF set if phenomenologically useful.

## Matched LO result and provisional NLO/LO K-factor (November 12, 2026)

The user completed `order=LO`, `fixed_order=on` on the *same six-channel NLO
process output* at `sqrt(s)=13.6 TeV`, using the previously prepared
NNPDF40 NLO input and `M_Q=225 GeV` (the exact saved run card should be
retained with the result):

| Quantity | One spectator hidden color | N_IC=2 spectator |
|---|---:|---:|
| Matched LO total | **3.437 ± 0.007 pb** | **6.874 ± 0.014 pb** |
| NLO central (PDF replicas extraction) | **3.9893368 ± 0.0057 pb** | **7.9786736 ± 0.0114 pb** |
| NLO/LO fixed-NLO-PDF K-factor | **1.16070 ± 0.00289 (MC only)** | **1.16070 ± 0.00289 (MC only)** |
| NLO enhancement vs matched LO | **+16.07%** | **+16.07%** |
| NLO PDF uncertainty | **±0.0573033 pb (1.436%)** | **±0.1146067 pb (1.436%)** |

The ratio uncertainty is obtained by propagating the independently quoted
LO and NLO integration errors; this is only a diagnostic estimate and
ignores any statistical correlations. This K-factor uses **the same NLO PDF**
for numerator and denominator: it is not the conventional ratio between a
full NLO calculation with NLO PDFs and a LO prediction using LO PDFs.
The hidden-color multiplicity cancels in the ratio.

The seven-point NLO scale uncertainty was subsequently **recovered
and verified** from the NLO run (see numerical result below). Outstanding
items for a broader publication-level theory assessment include independent
electroweak-production cross-checks, alpha_s/PDF-set variations, potential
NLO EW and infracolor effects, and additional quirk mass benchmarks.
The active `SubProcesses/job_status.pkl` may now reflect the later LO run:
use the *already-saved NLO JSON*, never the overwritten LO job weights.

## Completed 7-point NLO scale-envelope validation — 225 GeV benchmark

On November 12, 2026 the user independently ran
`scripts/extract_triplet_nlo_7point.py` on the saved **17 NLO integration jobs**
*before* the matched LO integration. The output was preserved as
`results/triplet_nlo_normalisation.json` (1.3 kB, local machine) and
`results/triplet_nlo_7point.txt`.

The reconstructed central `3.9893368 pb` exactly matches the PDF-member
central and agrees with MadGraph's rounded `3.989 pb` summary.

| muR/mu0 | muF/mu0 | Cross section (pb) | Standard 7-point? |
|---:|---:|---:|---|
| 0.5 | 0.5 | 4.0477817 | Yes — **maximum** |
| 0.5 | 1.0 | 4.0380961 | Yes |
| 0.5 | 2.0 | 4.0438984 | No — antipodal |
| 1.0 | 0.5 | 4.0021631 | Yes |
| 1.0 | 1.0 | 3.9893368 | Yes — **central** |
| 1.0 | 2.0 | 3.9898075 | Yes |
| 2.0 | 0.5 | 3.9634451 | No — antipodal |
| 2.0 | 1.0 | 3.9480101 | Yes |
| 2.0 | 2.0 | 3.9440466 | Yes — **minimum** |

The two antipodal pairs are *not* extrema; hence the standard 7-point
and complete 9-point envelopes happen to coincide:

- **7-point NLO scale:** +0.058445/-0.045290 pb
  (**+1.465%/-1.135%**) for one hidden-color copy.
- **PDF-only NNPDF40 68% RMS:** +/-0.0573033 pb (**1.436%**).
- **MC integration:** +/-0.0057 pb.
- For **N_IC=2** spectators the central and errors are doubled:
  `7.9786736 pb`, `+0.116890/-0.090580 pb` (scale),
  `+/-0.1146067 pb` (PDF), `+/-0.0114 pb` (MC).

One-hidden-color final benchmark (independent errors **quoted separately**):

```text
sigma_NLO(pp -> all six Dirac-triplet DY channels; sqrt(s)=13.6 TeV,
          mQ=225 GeV, NNPDF40_nlo_as_01180, N_IC spectator=1)
   = 3.9893 (+0.0584/-0.0453) [7pt scale]
            +/- 0.0573 [68%-CL PDF replica RMS]
            +/- 0.0057 [MC integration] pb

sigma_LO (same NLO PDF and central scale) = 3.437 +/- 0.007 pb
K_NLO/LO(same NLO PDF) = 1.1607 +/- 0.0029 (MC integration only)
```

The seven-point envelope does **not** include PDF or alpha_s uncertainty.
The fixed-NLO-PDF K factor is not the conventional NLO(PDF-NLO)/LO(PDF-LO)
ratio. The result is a **validated fixed-order NLO QCD production
benchmark**, not yet a completed NLO+PS QuirkSim event-production chain or
a full infracolor calculation.

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

Items 1--5 are satisfied for the **M=225 GeV fixed-order spectator production
benchmark** to the stated level: MG5 3.5.1 locally generated, integrated,
checked NLO poles, reconstructed NNPDF replica RMS and the standard
7-point scale envelope. Independent external comparisons and wider
mass/scale/PDF scheme checks in item 6 are still open, as are all
NLO+PS and quirk-string/decay interfaces.

## Provenance

- MadGraph5_aMC source `models/loop_sm` (standard QCD NLO UV/R2),
  upstream: https://github.com/mg5amcnlo/mg5amcnlo
- MadGraph NLO usage: Alwall et al. (2014), arXiv:1405.0301.
- Vectorlike lepton NLO QCD reference templates:
  https://github.com/prudhvibhattiprolu/VLL-UFOs
- FeynRules UFO NLO examples: https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki/NLOModels

No hidden SU(N_IC) gauge interactions or worldsheet dynamics are
misrepresented as SM QCD corrections.
