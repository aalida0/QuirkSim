"""Bounded, unmodified QuirkPheno transport smoke test; not a recoil result."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess

PIN = "d91cb27ed3786bc63900faad892513841b5b6ed4"
URL = "https://github.com/tvami/QuirkPheno.git"
LCG = Path("/cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc13-opt/setup.sh")

def run_logged(cmd, cwd, log, seconds):
    with log.open("w") as stream:
        try:
            p = subprocess.run(cmd, cwd=cwd, stdout=stream,
                               stderr=subprocess.STDOUT, timeout=seconds)
            return {"returncode": p.returncode, "timeout": p.returncode in (124, 137)}
        except subprocess.TimeoutExpired:
            return {"returncode": None, "timeout": True}

def prepare(out, suite="smoke"):
    out.mkdir(parents=True, exist_ok=False)
    # Stable equal-mass transport control only. Multiplet decay is absent upstream.
    mass = 225.
    energy = mass / math.sqrt(1 - .4**2)
    px = .4 * energy
    pz = .25 * energy / math.sqrt(1 - .25**2)
    cases = []
    duration_ps = 20 if suite == "smoke" else 130
    levels = ((0, .01, 1e-8), (1, .005, 5e-9)) if suite == "smoke" else tuple(
        (level, .01 / 2**level, 1e-8 / 2**level) for level in range(1, 5))
    for field in (0, 2):
        for level, boost, merge in levels:
            name = f"stable_B{field}_level{level}"
            card = out / (name + ".card")
            card.write_text(f"world 100 100 G4_Galactic\nfield {field} 0 100 0 100\n")
            macro = out / (name + ".mac")
            macro.write_text(f"""/quirk/mass {mass} GeV
/quirk/lambda 774 eV
/quirk/charge 1
/quirk/firstStringLength 1e-6 mm
/quirk/maxBoost {boost}
/quirk/maxMerge {merge} mm
/quirk/crossingLength 0 mm
/quirk/maxRadius 0 cm
/quirk/maxZ 0 cm
/qsim/card {card}
/qsim/gen/mode gun
/qsim/gun/p1 {px:.17g} 0 {pz:.17g} GeV
/qsim/gun/p2 {-px:.17g} 0 {pz:.17g} GeV
/qsim/maxTime {duration_ps / 1000:.17g} ns
/qsim/out/stepPrescale 1
/qsim/out/towers false
/qsim/out/file {out / (name + '.root')}
/run/initialize
""")
            cases.append({"name": name, "macro": str(macro),
                          "macro_sha256": hashlib.sha256(macro.read_bytes()).hexdigest(),
                          "field_tesla": field, "level": level,
                          "max_boost": boost, "max_merge_mm": merge})
    return {"upstream": URL, "source_commit": PIN,
            "scope": "STABLE_POINT_PRODUCTION_TRANSPORT_ONLY_NO_RECOIL",
            "suite": suite,
            "mass_gev": mass, "lambda_ev": 774, "nominal_stop_ps": duration_ps,
            "mass_convention": "stable transport control; not loop-split charged parent",
            "beta_each_in_pair_cm": .4, "beta_pair_cm_z": .25,
            "first_string_vector_length_mm": 1e-6,
            "crossing_restart_enabled": False,
            "timeout_seconds_per_case": 120,
            "initial_state_matches_stretched_conformal_pilot": False,
            "numerical_validation": "NOT_RUN", "physical_survival_probability": None,
            "cases": cases}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--source-dir", default=str(Path.home() / "QuirkPheno-public-pilot"))
    p.add_argument("--prepare-only", action="store_true")
    p.add_argument("--suite", choices=("smoke", "refinement"), default="smoke",
                   help="refinement: 130 ps, B=0/2 T, four tighter levels, no decay")
    a = p.parse_args()
    out = Path(a.output_dir).expanduser().resolve()
    state = prepare(out, a.suite)
    manifest = out / "manifest.json"
    def save():
        manifest.write_text(json.dumps(state, indent=2) + "\n")
    save()
    print(f"Manifest: {manifest}", flush=True)
    if a.prepare_only:
        return
    # Run directly on the node: this script never submits Condor jobs.
    if not LCG.is_file():
        state["bootstrap_error"] = f"Required LCG setup unavailable: {LCG}"
        save()
        raise SystemExit(state["bootstrap_error"])
    source = Path(a.source_dir).expanduser().resolve()
    if not source.exists():
        result = run_logged(["git", "clone", URL, str(source)], out, out / "clone.log", 120)
        if result["returncode"] != 0:
            state["bootstrap_error"] = result
            save()
            raise SystemExit("Clone failed; see clone.log")
    # An existing checkout must already be pinned and clean; never overwrite user changes.
    head = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    dirty = subprocess.check_output(["git", "-C", str(source), "status", "--porcelain"], text=True).strip()
    if dirty:
        raise SystemExit("Source checkout has changes; use a fresh --source-dir")
    if head != PIN:
        result = run_logged(["git", "checkout", "--detach", PIN], source, out / "checkout.log", 60)
        if result["returncode"] != 0:
            raise SystemExit("Pinned checkout failed; see checkout.log")
    state["observed_source_commit"] = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    build = out / "build"
    def lcg(command):
        # subprocess timeout bounds the whole process tree via GNU timeout.
        return ["timeout", "--kill-after=10s", "600s", "bash", "-c",
                "source " + shlex.quote(str(LCG)) + " && " + command]
    cmd = "cmake -S " + shlex.quote(str(source)) + " -B " + shlex.quote(str(build))
    cmd += " -DQSIM_WITH_HEPMC3=OFF -DQSIM_WITH_PYTHIA8=OFF"
    cmd += " && cmake --build " + shlex.quote(str(build)) + " -j1"
    state["build"] = run_logged(lcg(cmd), source, out / "build.log", 620)
    save()
    if state["build"]["returncode"] != 0:
        raise SystemExit("Build failed; see build.log")
    exe = build / "quirksim" / "quirksim"
    for case in state["cases"]:
        print("Running " + case["name"], flush=True)
        command = "timeout --kill-after=5s 120s " + shlex.quote(str(exe))
        command += " --mt -t 1 -n 1 -s 1 -m " + shlex.quote(case["macro"])
        case["execution"] = run_logged(lcg(command), source,
                                       out / (case["name"] + ".log"), 145)
        root = out / (case["name"] + ".root")
        case["root_exists"] = root.is_file()
        # Exit zero and a ROOT file are smoke evidence, never convergence certification.
        save()
    state["numerical_validation"] = "PENDING_ROOT_TRAJECTORY_AUDIT"
    save()
    failed = [c["name"] for c in state["cases"] if c["execution"]["returncode"] != 0 or not c["root_exists"]]
    if failed:
        raise SystemExit("Incomplete cases: " + ", ".join(failed) + "; retain all logs for audit")
    print("Finished. Return the output directory for trajectory and stop-status audit.")

if __name__ == "__main__":
    main()

