# Triplet weak-decay radius and ATLAS pixel-barrel cylinder diagnostic

**Status:** endpoint propagation, decay-position serialization, first-outward cylinder
crossings, and summary tests pass on GitHub Actions. User's original 100-event
propagation records lack vertex coordinates; regenerate them with the updated
`quirksim/production_dynamics.py` before using the new summary.

## Physics meaning

- Initial charge and final charge are not enough to infer track visibility.
- Weak decay is sampled in proper time and recorded at its true simulated
  lab-frame time and 3D position, including transverse decay radius.
- For each produced quirk endpoint and each cylinder radius, integrate the
  same ODE with a nonterminal outward root condition
  `sqrt(x(t)^2+y(t)^2)-R = 0`, direction `+1`.
- Record the **first outward crossing** whether charged or neutral.
- Classify *initially charged* endpoints as: first crossing still charged,
  first crossing only after neutralization, or no crossing within the window.
- Initially neutral endpoints have their own reached/unreached counts.
- The curves are **geometric proxies, not detector hits or efficiencies**.
  No endcaps, z acceptance, finite sensor lengths, material, charge deposition,
  hit time windows, track fitting, reconstruction or triggering are modeled.
- Charge-changing pions and leptons are not propagated in detector material.

## Barrel cylinder radii

| Layer | Approx. radius (cm) |
|---|---:|
| IBL | 3.325 |
| Pixel B | 5.05 |
| Pixel 1 | 8.85 |
| Pixel 2 | 12.25 |

Reference: ATLAS Collaboration, *Alignment of the ATLAS Inner Detector in
Run 2*, Eur. Phys. J. C 80 (2020) 1194, Table 1.
https://link.springer.com/article/10.1140/epjc/s10052-020-08700-6

## Recompute real MG5+Pythia8 events with new observables

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
python3 -m pytest -q tests/test_production_dynamics.py tests/test_triplet_barrel_summary.py
python3 scripts/propagate_triplet_isr.py results/triplet_pythia_isr.jsonl \
  --nevents 100 --lambda-ev 774 --b-tesla 2 \
  --force-model ng_endpoint --mode recoil --duration-ps 500 \
  --output results/triplet_isr_recoil100_geometry.jsonl
python3 scripts/summarize_triplet_propagation.py \
  results/triplet_isr_recoil100_geometry.jsonl
python3 scripts/summarize_triplet_barrel.py \
  results/triplet_isr_recoil100_geometry.jsonl \
  --json-out results/triplet_barrel100_geometry.json
```

Records contain `decays[].position_cm` and `decays[].radius_cm`,
and `first_outward_barrel_crossings[]` with endpoint, time, z, azimuth,
charge at crossing and radius, along with `barrel_layer_radii_cm`.

## Interpretation and limitations

The physical path may encounter a cylinder multiple times. Only the *first*
outward crossing is retained, so this is not a repeated-hit estimate.
The 500 ps integration window censors both late decays and late crossings.
At the maximum speed, a particle travels only 15 cm in 500 ps; this is
especially relevant when interpreting non-crossers of the outer pixel barrel.
Post-decay string propagation remains in the approximate straight-string
endpoint-force model, not a full worldsheet. No B-field map or actual finite
barrel support is included. Validate worldsheet and detector-level effects
before using these fractions in an ATLAS sensitivity claim.
