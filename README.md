# QuirkSim: charged-to-neutral macroscopic quirks

Independent public research on electroweak quirk multiplets, weak pion emission and recoil during macroscopic string oscillations. **Not an ATLAS analysis.** No unpublished internal ATLAS code is used.

## Current implementation status

Development branch: `development/standalone-simulator`.

- Dirac SU(2)L doublet Y=-1/2 and Dirac SU(2)L triplet Y=0, each in fundamental hidden SU(2)IC.
- One-loop charged-neutral mass splittings and leading weak-decay widths (reference formulas, not precision electroweak predictions).
- Exponentially sampled charged endpoint decay clocks in **proper time**, exact two-body pion recoil in the parent rest frame and Lorentz boost to the lab.
- Relativistic two-endpoint ODE with a smooth **central straight-string proxy** `V(r) = sigma (sqrt(r² + a²) - a)`.
- Three modes: `stable`, `no_recoil` (unphysical diagnostic) and `recoil`. Leptonic decays are explicitly marked **unresolved** rather than simulated incorrectly.

**Scientific limitations:** This is **not** a full moving Nambu–Goto worldsheet, and the string proxy omits radiation, annihilation, detector material, magnetic fields and LHC production. Any minimum separation measured from creation includes a trivial near-contact point, so it cannot be called an annihilation probability. The first goal is to validate physical dynamics against published calculations with nonradial initial conditions and matched decay-clock replay.

## Setup

```bash
git clone https://github.com/aalida0/QuirkSim.git
cd QuirkSim
git checkout development/standalone-simulator
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
pytest -q
python -m quirksim.cli --model triplet --mass 225 --lambda-ev 774 --mode recoil --duration-ps 500 --seed 8 --output result.json
```

## Physics roadmap

1. Verify one-loop splittings and LO pion/leptonic decay widths against original publications, scheme by scheme.
2. Validate exact undecayed radial solutions, energy conservation and convergence with solver and regulator variation.
3. Simulate both endpoint decays including leptonic channels and exact matched-clock replay.
4. Validate an independent moving Nambu–Goto string treatment for nonradial trajectories.
5. Add an uncertainty-bracketed contact-dependent infracolor radiation model; keep geometric contacts separate from annihilation.
6. Add a separately generated doublet/triplet UFO and realistic pp initial conditions.

## References

- Thomas & Wells (1998), arXiv:hep-ph/9804359.
- Kang & Luty (2009), arXiv:0805.4642.
- Martin (2011), arXiv:1012.2072.
- *The Quirk Trajectory* (2020), arXiv:2002.07503.

No results from this pilot should be interpreted as collider exclusions or validated physical annihilation rates.
