# QuirkPheno transport reuse pilot (2026-11-11)

Status: PREPARED, CXX BUILD BLOCKED BY LOCAL ENVIRONMENT; NO TRANSPORT OR RECOIL VALIDATION.

Upstream pinned source: https://github.com/tvami/QuirkPheno/commit/d91cb27ed3786bc63900faad892513841b5b6ed4
Apache-2.0; quirkcore derived from ATLAS transport, with separate Geant4 notice.
No upstream code copied into QuirkSim.

## Observed evidence

Downloaded 38 build/core/application files locally. CMake configuration with
HepMC3/Pythia8 disabled stops at CMakeLists.txt:11: geant4-config unavailable.
CVMFS and apptainer also unavailable here. No C++ build or simulation completed.
No AF node access is available in this session.

The upstream gun creates both endpoints at one vertex.
QuirkInfracolorForce::StartTracking builds the initial string from endpoint
four-velocities and firstStringLength. Our validated conformal pilot starts
with a 0.1 cm stretched string. Matching mass and Lambda, or setting that length
parameter to 1 mm, does not establish the same spatial worldsheet.

The upstream particle definitions are stable. ConstructProcess removes generic
processes. Adding only a decay table will not implement charge-changing weak
decay with retained string and partner.

## Deliverable and limits

scripts/af_quirkpheno_transport.py pins the clean upstream checkout, builds
with one CPU in upstream LCG_110, and runs four sequential one-event stable
vacuum controls: Bz=0/2 T and two maxBoost/merge settings. Uses mass 225 GeV,
Lambda 774 eV, beta each=.4 in pair CM, beta CM z=.25.
This is a transport control, not the loop-split charged triplet benchmark.

Nominal stop 20 ps, per-case timeout 120 s plus 5 s kill grace, build timeout
600 s. The upstream maxTime kills after an exceeding step: actual end times
must be checked. Every transport step saved. Experimental crossing restart
disabled. No Condor, no shell set flags, no decay or annihilation claims.
The runner requires an EL9 node with upstream LCG_110 CVMFS. It fails clearly
if setup is missing; it does not select another untested Geant4 version.
--prepare-only needs neither Geant4 nor network.

Python compilation, generated macro/card inspection against upstream command
definitions, and missing-CVMFS failure path checked locally. Runtime results,
refinement, and speed remain unknown. Manifest numerical status stays pending
ROOT audit even after successful execution.

## Next checks

1. Inspect ROOT end reasons, time coverage and aborts. Upstream maxSep/nCross
   use asynchronous last endpoint positions; do not interpret as equal-time minima.
2. Compare each backend under refinement, then implement a shared initial-state
   adapter before any cross-backend overlay. Interpolation errors need bounds.
3. Extend weak decay while preserving the partner/string history; test causal
   recoil propagation and local four-momentum balance without string restart.
4. Only then evaluate matched recoil controls. Glueball radiation and quantum
   annihilation remain independent physical inputs. Doublet stays paused.

## AF compute-node command

Run directly, not through Condor:

```bash
cd ~/QuirkSim
git pull --ff-only origin research/variational-worldsheet
python3 scripts/af_quirkpheno_transport.py --output-dir /data/aalida/QuirkSim/quirkpheno_transport_$(date +%Y%m%dT%H%M%S)
```

Default separate checkout: ~/QuirkPheno-public-pilot.
Return the complete output directory (manifest, build/case logs and ROOT files).
