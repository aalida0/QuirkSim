# Dirac SU(2)L Y=0 triplet: NLO QCD + EW development

**STATUS (2026-11-12): SOURCE + NLOCT BUILD WORKFLOW PREPARED; NO FULL NLO-EW UFO CERTIFIED.**

This is **not** a claim of successful full NLO QCD+EW implementation.
Do **not** launch physical `[QED]` calculations using the validated
`models/DiracTripletEW_NLO_UFO`, which has only inherited SM **QCD** UV/R2.

## Physics definition

The model contains one vectorlike Dirac fermion triplet
`Q=(Q1+,Q0,Q2-)`, with `Y=0`, mass `M_Q`, and conserved quirk number:

```text
L = L_SM + bar Q (i gamma.mu D_mu - M_Q) Q
Q ~ (3,0) of SU(2)L x U(1)Y, singlet of SM color
T3 = diag(+1,0,-1), g(Q1-Q0-W) = g(Q0-Q2-W) = g2.
```

It contains THREE distinct Dirac species and their antiparticles; the
neutral Q0 is **Dirac**, unlike the **Majorana neutral** component of
minimal type-III seesaw models. No mixing with SM leptons, no Yukawas
or tree-level splitting is assumed. At tree level Q1/Q0/Q2 share one
gauge-invariant `M_Q` parameter. In the physical basis Goldstone
couplings to the vectorlike fermion vanish for this gauge-invariant
mass, but the SM ghost and Goldstone sector is retained in Feynman gauge.

A perturbative fixed-order `NLO QCD+EW` prediction means

```text
sigma = sigma_LO[alpha_EW^2]
      + sigma_NLOQCD[alpha_s alpha_EW^2]
      + sigma_NLOEW[alpha_EW^3]
```

for pure DY topology, **not** a mixed two-loop
`alpha_s*alpha_EW^3` calculation.

NLO-EW Drell-Yan involves virtual heavy-fermion and SM gauge
self-energies, UV parameter/wave-function counterterms, one-loop
rational R2 terms, photon emission and photon-induced production.
Its renormalization scheme (on-shell/complex-mass, alpha(0)/Gmu),
weak-mass pole renormalization, and QED initial-state factorization
must be fixed consistently; an SM-only NLO EW CT model is insufficient.

### Hidden-color issue unique to NLO EW

At **NLO QCD**, under the spectator hidden-color assumption, the
external hidden index supplies a simple factor `N_IC` to the inclusive
triplet production rate.

At **NLO EW**, closed virtual quirk loops in gauge-boson self-energies
also carry the hidden-color trace `N_IC`. The NLO EW correction can
therefore depend on `N_IC` *internally* as well as via external
multiplicity. A result calculated with one virtual triplet and then
multiplied by `N_IC` is **not** generally the exact NLO EW prediction
for `N_IC>1`. To include this term consistently, model all hidden-index
copies in virtual loops or derive and validate explicit trace factors
within the EW renormalization scheme. Do not double-count the
external hidden-index multiplicity.

## Existing public resources audited

1. **MadGraph5_aMC `loop_qcd_qed_sm_Gmu`**: includes UV/R2
   counterterms for **SM** NLO QCD+EW, with mixed-coupling
   process syntax `[QCD QED]`. This is the right *SM* subtraction,
   factorization and renormalization reference; **not** a complete
   BSM-triplet EW model. Frederix et al., *JHEP* 07 (2018) 185,
   [doi:10.1007/JHEP07(2018)185](https://doi.org/10.1007/JHEP07(2018)185).
   Official source: https://github.com/mg5amcnlo/UFOmodeldb .
2. **FeynRules Type-III seesaw and vectorlike lepton models**:
   Ajjath, Fuks, Shao, Simon,
   *Phys. Rev. D* **107** (2023) 075011,
   [arXiv:2301.03640](https://arxiv.org/abs/2301.03640).
   Closely related SU(2)L triplet production and NLO+PS **QCD**
   benchmarks. Not a complete electroweak-renormalized Dirac
   quirk triplet. Public source and model files:
   https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki/NLOModels .
   The original Type III seesaw model is Majorana neutral and
   includes Yukawa mixing; not directly identical to our Dirac model:
   https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki/TypeIIISeeSaw .
3. **FeynRules, FeynArts and NLOCT**:
   C. Degrande, *Comput. Phys. Commun.* **197** (2015), 239,
   [doi:10.1016/j.cpc.2015.08.015](https://doi.org/10.1016/j.cpc.2015.08.015).
   Example renormalization pipeline:
   https://indico.ihep.ac.cn/event/5301/contributions/74634/attachments/38345/44341/Tutorial_NLO.pdf .
   More detailed EW-loop application: *MadDM v3.2*,
   [doi:10.1140/epjc/s10052-023-11377-2](https://doi.org/10.1140/epjc/s10052-023-11377-2).
   `QCDOnly -> False` enables the EW construction, but software
   support and renormalization conditions **must be validated**;
   it does not guarantee a successful completed counterterm
   export in an arbitrary gauge theory.
4. **EW pole-mass difference**:
   Cirelli, Fornengo, Strumia, *Nucl. Phys. B* **753** (2006),
   178, [hep-ph/0512090](https://arxiv.org/abs/hep-ph/0512090).
   Ibe, Matsumoto, Sato, *Phys. Lett. B* **721** (2013),
   252, [arXiv:1212.5989](https://arxiv.org/abs/1212.5989)
   for two-loop corrections to the wino-like triplet splitting.
   The vectorlike Dirac triplet has the same minimal one-loop
   charged-neutral gauge-boson self-energy difference as one
   triplet copy with corresponding gauge couplings (but has
   additional particle content and independent global charge).
5. **UFO 2.0 format:** Degrande et al.,
   [doi:10.1140/epjc/s10052-023-11780-9](https://doi.org/10.1140/epjc/s10052-023-11780-9).
   Genuine NLO UFOs need `CT_vertices.py`,
   `CT_couplings.py` and `CT_parameters.py` in addition
   to the tree-level files, and their content must match the
   virtual loop particle spectrum and renormalization scheme.

**Photon PDFs:** A full photon-initiated NLO EW setup needs
a PDF set with a photon component. The original QCD-only
`NNPDF40_nlo_as_01180`, ID **331700**, is *not* that dataset.
Use `NNPDF40_nlo_as_01180_qed`, LHAPDF ID **335900** as
a candidate input and validate its member definitions, QED
factorization scheme and SM alpha configuration.
Official catalogue:
https://www.lhapdf.org/pdfsets.html .
The SM / quirk amplitudes must include photon-induced Born
channels consistently if they are of the chosen coupling order.

## Workstation setup: FeynArts is the missing dependency

On the current workstation:
- FeynRules is installed at `$HOME/.Wolfram/Applications/FeynRules`.
- NLOCT.m is installed inside the same FeynRules directory.
- `find $HOME/.Wolfram $HOME/.Mathematica $HOME/codes -name FeynArts.m`
  returned no results as of 2026-11-12.

Use FeynArts 3.11 for an initial compatibility test with the
FeynRules 2.3.x/NLOCT generation; its use in a FeynRules/NLOCT
1.02 setup is documented in the public Zee-Babu implementation.
The official download URL is
`https://feynarts.de/FeynArts-3.11.tar.gz`; the FeynArts
author also distributes 3.12, which can be tested separately.

```bash
cd ~/QuirkSim
export FEYNRULES_HOME="$HOME/.Wolfram/Applications/FeynRules"
export WOLFRAM_APP_ROOT="$HOME/.Wolfram/Applications"
mkdir -p "$WOLFRAM_APP_ROOT"
curl -fL --retry 3 -o /tmp/FeynArts-3.11.tar.gz \
  https://feynarts.de/FeynArts-3.11.tar.gz
tar -tzf /tmp/FeynArts-3.11.tar.gz >/dev/null
tar -xzf /tmp/FeynArts-3.11.tar.gz -C "$WOLFRAM_APP_ROOT"
export FEYNARTS_HOME="$WOLFRAM_APP_ROOT/FeynArts-3.11"
test -f "$FEYNARTS_HOME/FeynArts.m"
test -f "$FEYNRULES_HOME/NLOCT.m"
wolframscript -code 'AppendTo[$Path,Environment["FEYNARTS_HOME"]]; Get["FeynArts.m"]; Print["FEYNARTS_IMPORT_OK"]'
```

If the tarball has a different top-level directory name, set
`FEYNARTS_HOME` to the directory that directly contains
`FeynArts.m`; do not assume FeynArts is installed elsewhere.
Run stage 01 first; stage 02 requires its exported model, and stage
03 requires the `.nlo` file emitted by stage 02.

**Risk**: Stage 01 and later Wolfram scripts are source-level
prototypes and have not passed Mathematica's vertex or NLOCT
calculations. Record the first-stage log and correct model-specific
FeynRules syntax/renormalization choices before trying stage 02.
Do not interpret installed dependencies as a completed EW UFO.

## Implemented in QuirkSim

### A. One-loop EW pole-mass difference — functional, tested

```bash
cd ~/QuirkSim
python3 scripts/triplet_ew_pole_spectrum.py --mass 225 \
  --output results/triplet_ew_225_spectrum.json
python3 -m unittest discover -s tests \
  -p 'test_triplet_ew_pole_spectrum.py' -v
```

Implements the full finite-mass one-loop integral

```text
Delta M = alpha2 M/(4 pi) *
   [ f(MW/M) - cos^2(thetaW) f(MZ/M) - sin^2(thetaW) f(0) ]
f(r) = 2 int_0^1 dx (1+x) log[x^2 + (1-x) r^2]
```

One **input neutral pole mass** defines the charged spectrum in a
specific renormalization convention; `MQ0=MQC` remains in the
gauge-invariant **tree-level NLO production UFO**.

The default spectrum benchmark uses `alphaEW(MZ)^-1=127.9`,
`GF=1.1663787e-5 GeV^-2`, `MZ=91.1876 GeV` with MW
**derived** as in the SM.fr tree relation. At `M=225 GeV`,
the result is approximately **161.844 MeV** in this input
convention. With `MW=80.379 GeV` supplied explicitly, the
answer changes at higher-order due to the different input
scheme; record that choice. The heavy-mass asymptotic
one-loop limit in the default scheme is about 166.4 MeV.
This is an EW *pole-mass difference*, **not** a complete
renormalized NLO EW cross-section calculation.

### B. FeynRules source Lagrangian — source prepared, NOT run

`theory/dirac_triplet_ew/TripletDirac.fr` adds three Dirac
fields and their eight vectorlike electroweak gauge currents
to a stock FeynRules `SM.fr` Lagrangian, with **one**
mass parameter `MQ`. The source is defined in the broken
mass basis, preserving gauge coupling relations and the
SM Feynman-gauge Goldstones. Requires a Mathematica-backed
vertex extraction and comparison to our working UFO before
using it for precision calculations.

### C. Three *separate-kernel* NLOCT stages — not executed

NLOCT's official tutorial explicitly requires separate
Mathematica kernels when switching between FeynRules and
FeynArts. Three driver scripts are supplied accordingly:

```bash
export FEYNRULES_HOME=/actual/path/to/FeynRules
export FEYNARTS_HOME=/actual/path/to/FeynArts

wolframscript -file theory/dirac_triplet_ew/01_prepare_feynarts.wls
wolframscript -file theory/dirac_triplet_ew/02_calculate_nloct.wls
wolframscript -file theory/dirac_triplet_ew/03_export_candidate_ufo.wls
```

These are **unvalidated driver scripts**, not proof that
NLOCT can successfully renormalize the entire model in
the chosen scheme or that the exported counterterms are
complete. They may require version-specific adaptation
of FeynRules/FeynArts/NLOCT and renormalization conditions.
The candidate UFO is named
`generated/DiracTripletEW_QCDEW_CANDIDATE_UFO`,
not the final `DiracTripletEW_NLO_QCDEW_UFO`.

### D. Fail-closed structural screening

```bash
python3 scripts/audit_triplet_nlo_ew_candidate.py --ufo \
  theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CANDIDATE_UFO
```

The check requires a UFO with full NLO files and
triplet-specific UV and R2 vertices. The existing NLO
QCD-only model fails this screen by design. Even if a
candidate passes the static screen, it is **not** certified:
only actual NLO virtual amplitude/UV-IR numerical tests,
independent amplitude checks and EW Ward identities can
support a full NLO-EW claim.

## Required hard physics acceptance gates

1. **Lagrangian and gauge invariance:** Independent
   FeynRules vertex comparison with `DiracTripletEW_NLO_UFO`
   at LO (two charged Dirac species, neutral Dirac Q0, 8
   gauge vertices, full-strength W couplings, no spurious
   tree Goldstone/Yukawa, conservation of quirk number).
2. **Renormalization prescription:** Explicitly document
   OS masses, EW charge scheme, tadpole scheme, W/Z mixing,
   electric charge, photon IR, field/mass CTs, the quirk
   on-shell mass convention, and 5-flavour QED MSbar
   factorization and photon PDF. Verify both SM and new
   triplet counterterms, including virtual quirk loops.
3. **UV and rational terms:** Produce complete nonzero
   `CT_vertices.py`, `CT_couplings.py`, `CT_parameters.py`
   for both SM and triplet EW with scheme-consistent
   rational pieces. Check all Lorentz and electric charge
   structures. Test that no SM counterterms were lost.
4. **MG5 fixed-order validation:** Compare the new UFO
   against SM `loop_qcd_qed_sm_Gmu` for SM-only processes.
   Check UV/IR pole cancellation for
   `u u~ > q1+ q1- [QED]`,
   `d d~ > q2- q2+ [QED]`,
   and charged-neutral `u d~ > q1+ q0~ [QED]`.
   Verify the full FKS QED soft/collinear checks, virtual
   stability, finite cross sections, and gauge independence.
5. **Mixed-coupling and PDFs:** Verify `[QCD QED]` includes
   requested Born powers and additive O(alpha_s alphaEW^2)
   and O(alphaEW^3) terms. Include photon-induced Born
   channels, QED PDFs (LHAPDF 335900 example), and
   the correct factorization and photon recombination
   conventions for IR-safe fiducial observables.
   Compare against the known NLO-QCD-only results and
   SM limits; independently validate NLO EW against
   an analytical or separate-loop-program result.
6. **Hidden multiplicity:** Derive all `N_IC` factors
   for external triplet lines and closed virtual
   triplet-loop EW counterterms; do not naïvely multiply
   the full NLO-EW coefficient obtained for one copy
   by two. Quirk hidden-sector dynamics remains absent.
7. **Pole masses:** Independently validate the finite
   one-loop splitting and test its value under the
   same electroweak input scheme. The mass splitting
   is applied to physical decays downstream and not
   inserted as an arbitrary SU(2)-breaking bare mass
   in the gauge-symmetric theory.

**No full NLO QCD+EW cross section has been generated
or verified in this project as of this document.**

## Existing production benchmark must not be overwritten

Our fully validated result is still *NLO QCD only* for
all six triplet channels at 13.6 TeV and M=225 GeV:

```text
sigma(N_IC=1) = 3.9893 (+0.0584/-0.0453) [7-point QCD scale]
                       +/- 0.0573 [NNPDF40 QCD PDF] pb
K(same NLO QCD PDF and central scales) = 1.1607
```

These numbers are from **NNPDF40_nlo_as_01180
(LHAPDF 331700)**, not the QED PDF.
They must not be relabeled as NLO EW or as a QED-PDF
prediction. See `studies/TRIPLET_NLO_QCD_UFO.md`.

## Status at a glance

| Artifact | Status |
|---|---|
| Validated six-channel NLO QCD UFO | **PASS** (existing model; untouched) |
| Charged-pair-only LO/NLO QCD cards | Prepared, numerical run pending |
| One-loop EW Dirac-triplet pole-mass difference | **Implemented; tested** |
| FeynRules Dirac source Lagrangian | Written, Mathematica vertex audit pending |
| NLOCT full EW UV/R2 driver scripts | Written; **not executed** |
| Full BSM EW UV/R2 | **NOT CALCULATED / NOT VERIFIED** |
| NLO QCD+EW UFO | **NOT VALIDATED / NOT RELEASED** |
| NLO EW xsec at 13.6 TeV | **NOT CALCULATED** |

### Resolved first-stage LoadModel abort (2026-11-12)

FeynArts 3.11 imports successfully on the user workstation; its
`SetDelayed::wrsym: Symbol Discard is Protected` import warning did not
prevent `FEYNARTS_IMPORT_OK`. The first stage originally failed **before**
any physics generation because `TripletDirac.fr` tried to read and append
to `M$ClassesDescription` during FeynRules' intermediate per-file loading,
when the SM file is not yet merged. The extension has been corrected to
declare only the three triplet classes and let `LoadModel[SM.fr,
TripletDirac.fr]` perform the merge.

The first stage now prints `FR_TRIPLET_MERGED_MODEL_LOADED` after loading,
uses `qsGeneratedDir` instead of the shadowing global name `out`,
and exits nonzero if the expected FeynArts `.mod` output is absent.
It should **not** be mistaken for a completed NLOCT model.

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
export FEYNRULES_HOME="$HOME/.Wolfram/Applications/FeynRules"
# retain FEYNARTS_HOME from actual install
set -o pipefail
wolframscript -file theory/dirac_triplet_ew/01_prepare_feynarts.wls \
  2>&1 | tee results/triplet_ew_stage1.log
grep -E 'FR_TRIPLET_MERGED_MODEL_LOADED|FR_STAGE1_FEYNARTS_MODEL_OK|ERROR|Abort|::' \
  results/triplet_ew_stage1.log | tail -n 60
```

This resolves the identified loading misuse in source; the Mathematica
renormalization and output have **not** been executed or checked here.
