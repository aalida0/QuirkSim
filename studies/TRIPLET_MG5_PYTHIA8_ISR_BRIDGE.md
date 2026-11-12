# QuirkSim: MG5 -> Pythia8 QCD ISR -> charged/neutral string propagation

**Status (2026-11-12): six-channel synthetic end-to-end smoke PASS,
GitHub Actions run `43245033841`. Real 1000-event sample is user-local and
has not been run through Pythia8 here.**

## Why

The LO 2->2 MG5 sample has zero pair pT. QCD ISR gives the pair realistic
transverse recoil, especially important for detector acceptance. The primary
deliverable is *initial post-shower quirk four-momenta*, not a precision jet
cross section. The PYTHIA chain is unmerged LO+PS. Hard-ISR tails need
explicit 0+1-jet matching/merging and uncertainty studies.

## Local setup

From the `feature/dirac-triplet-ew-lo` branch:

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
# Check whether Pythia8 is already available:
command -v pythia8-config || true
find ~/codes ~/Codes -path '*/include/Pythia8/Pythia.h' 2>/dev/null | head
```

If installed, expose its prefix by `export PYTHIA8_DIR=/path/to/pythia8`
(where include/Pythia8/Pythia.h and lib/libpythia8.so live).
If you use conda, a reproducible environment is:

```bash
conda create -n quirksim-pythia -c conda-forge python=3.11 pythia8 numpy scipy pip
conda activate quirksim-pythia
python -m pip install -e .
export PYTHIA8_DIR="$CONDA_PREFIX"
```

## Shower your generated 1000-event LHE sample

```bash
python scripts/run_triplet_pythia8.py \
  results/pp_dirac_triplet_13p6TeV_LO/Events/run_01/unweighted_events.lhe.gz \
  --out results/triplet_pythia_isr.jsonl --nevents 0 --seed 12345
python scripts/summarize_triplet_isr.py results/triplet_pythia_isr.jsonl
```

`--nevents 0` reads all events. The wrapper decompresses LHE transparently
and compiles the Pythia8 runner from source. Pythia8 recognizes custom PDGs
9000111, 9000112 and 9000113 and their distinct antiparticles.
All quirks are stable in Pythia8. QCD ISR/FSR enabled; QED shower,
MPI and hadronization disabled for this controlled study.

## Feed momenta into QuirkSim

```bash
python scripts/propagate_triplet_isr.py results/triplet_pythia_isr.jsonl \
  --nevents 3 --lambda-ev 774 --b-tesla 2 \
  --force-model ng_endpoint --mode recoil --duration-ps 500 \
  --output results/triplet_pythia_recoil_pilot.jsonl
```

Do not start large propagation jobs on shared login nodes.
Each event retains both PDGs and four-momenta, with the true electric charge
per endpoint (including initially neutral Q0 and Q0bar). The string persists
after a charged endpoint changes to neutral in a weak decay.

**Interpretation:** initially one charged endpoint is a production topology;
it is NOT a one-decayed event. The report separately records
`initial_charged_endpoints`, `n_weak_decays`, and `final_charged_endpoints`.

**Physics limitations:** no dynamical SU(N_IC), infraglue emission,
hadronization, detector response, material effects, full Nambu-Goto
worldsheet, quirkonium, or annihilation. A straight-string endpoint-force
approximation is used. The original LO MG5 cross section should NOT
be read as a precision LHC prediction. Hidden color multiplicity N_IC
is separate from the Pythia8 shower.

**Validation:** GitHub Actions passes six synthetic channels through
Pythia8 8.312 and QuirkSim, including 4 one-charged and 2 two-charged
topologies. The synthetic test is NOT a cross-section or realistic-recoil
benchmark. The actual 1000-event LHE file must be checked independently.

Refs: https://pythia.org/latest-manual/LHEF.html and
https://pythia.org/latest-manual/ParticleDataScheme.html
