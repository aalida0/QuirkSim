"""Pinned-source HTCondor pilots with explicit input/output transfers.

No physics ensemble is implied: these are selected numerical validation cases.
Preparation uses only the standard library and does not run on the login CPU.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


BRANCH = "research/variational-worldsheet"


def cases(suite):
    selected = []
    for model in ("triplet", "doublet"):
        for seed in (1, 25):
            for level in (0, 1):
                selected.append(dict(
                    name=f"{model}_recoil_seed{seed}_level{level}",
                    model=model, mode="recoil", seed=seed,
                    max_step_ps=.025 / 2**level,
                    ray_angle=.0001953125 / 2**level,
                    causal_fraction=.000390625 / 2**level))
    if suite == "remaining":
        return selected[:1]  # interrupted triplet seed 1, strict configuration
    return selected


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def safe_path(path):
    # Condor queue/executable syntax is deliberately kept simple and unambiguous.
    if any(c.isspace() or c in '\"\\$;' for c in str(path)):
        raise ValueError(f"Use a path without whitespace or Condor metacharacters: {path}")
    return path


def prepare(a):
    root = Path(__file__).resolve().parents[1]
    if git(root, "rev-parse", "--show-toplevel") != str(root):
        raise ValueError("Run from a real QuirkSim git checkout")
    if git(root, "branch", "--show-current") != BRANCH:
        raise ValueError(f"Switch to {BRANCH} before preparing the pilot")
    if git(root, "status", "--porcelain", "--untracked-files=no"):
        raise ValueError("Tracked changes must be committed or restored before source pinning")
    if not git(root, "ls-files", "scripts/af_conformal_validation.py"):
        raise ValueError("The AF helper must be present in the pinned git commit")
    python = safe_path(root / ".venv" / "bin" / "python")
    if not python.is_file():
        raise ValueError(f"Missing {python}; install QuirkSim in .venv first")
    if a.submit and not shutil.which("condor_submit"):
        raise ValueError("condor_submit is unavailable in this shell")
    if a.submit and a.suite == "remaining":
        raise ValueError("Single pilots run directly on the compute node; Condor is for the parallel convergence suite")
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    sha = git(root, "rev-parse", "HEAD")
    name = f"conformal_{a.suite}_{stamp}_{sha[:8]}"
    work = safe_path(root / "condor" / name)
    output = safe_path(Path(a.output_root).expanduser().resolve() / name)
    work.mkdir(parents=True, exist_ok=False)
    output.mkdir(parents=True, exist_ok=False)
    source = work / "source"
    source.mkdir()
    archive = work / "source.tar"
    subprocess.run(["git", "-C", str(root), "archive", "--format=tar",
                    f"--output={archive}", sha], check=True)
    subprocess.run(["tar", "-xf", str(archive), "-C", str(source)], check=True)
    freeze = subprocess.check_output([str(python), "-m", "pip", "freeze"], text=True)
    (work / "requirements.freeze.txt").write_text(freeze)
    manifest = dict(schema=1, source_commit=sha, branch=BRANCH, suite=a.suite,
                    created_utc=stamp, source_dir=str(source), python=str(python),
                    output_dir=str(output), work_dir=str(work),
                    duration_ps=1000., cases=cases(a.suite),
                    scope="SELECTED_NUMERICAL_PILOTS_PHYSICAL_SURVIVAL_UNKNOWN")
    manifest['transport'] = 'CONDOR_TRANSFER_WITH_SHARED_AF_PYTHON_RUNTIME'
    manifest_path = work / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.copy2(manifest_path, output / "manifest.json")
    shutil.copy2(work / "requirements.freeze.txt", output / "requirements.freeze.txt")
    (work / "cases.txt").write_text("".join(c["name"] + "\n" for c in manifest["cases"]))
    for case in manifest['cases']:
        (output / case['name']).mkdir()
    worker = work / "worker.sh"
    worker.write_text(f'''#!/bin/bash
set -eu
mkdir -p result source
if [ ! -x {python} ]; then
    echo "AF Python runtime unavailable: {python}. Source transfer is enabled, but this job still requires the AF virtualenv on the worker." | tee result/bootstrap_error.txt >&2
    exit 69
fi
tar -xf source.tar -C source
exec {python} -u source/scripts/af_conformal_validation.py run manifest.json "$1" --sandbox
''')
    worker.chmod(0o755)
    submit = work / "pilot.sub"
    submit.write_text(f"""universe = vanilla
executable = {worker}
arguments = $(case)
initialdir = {output}/$(case)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {archive},{manifest_path}
transfer_output_files = result
stream_output = True
stream_error = True
getenv = False
request_cpus = 1
request_memory = 4GB
request_disk = 2GB
+JobBatchName = "{name}"
output = {work}/$(case).$(ClusterId).$(ProcId).out
error = {work}/$(case).$(ClusterId).$(ProcId).err
log = {work}/scheduler.log
periodic_remove = (JobStatus == 2) && ((time() - JobCurrentStartDate) > 28800)
queue case from {work}/cases.txt
""")
    print(f"Pinned source: {sha}\nManifest: {manifest_path}\nResults: {output}", flush=True)
    print(f"Each case uses 1 CPU; eight-hour running limit. Cases: {len(manifest['cases'])}")
    if a.submit:
        dry = subprocess.run(["condor_submit", "-dry-run", str(work / "job.ads"), str(submit)],
                             text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (work / "dry-run.txt").write_text(dry.stdout)
        if dry.returncode:
            print(dry.stdout, end="")
            raise SystemExit(dry.returncode)
        # condor_submit validates its own syntax. Do not parse rendered ClassAds:
        # whitespace/serialization vary across versions, and site-added policy
        # must be assessed with condor_q -better-analyze after submission.
        result = subprocess.run(["condor_submit", str(submit)], text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (work / "submission.txt").write_text(result.stdout)
        print(result.stdout, end="")
        if result.returncode:
            raise SystemExit(result.returncode)
    else:
        print(f"Prepared only. Submit with: condor_submit {submit}")


def run(a):
    manifest = json.loads(Path(a.manifest).read_text())
    case = next(c for c in manifest["cases"] if c["name"] == a.case)
    destination = (Path.cwd() / "result" if a.sandbox
                   else Path(manifest["output_dir"]) / case["name"])
    destination.mkdir(exist_ok=True)
    result_path = destination / "evolution.json"
    if (destination / "completion.json").exists() or result_path.exists():
        raise ValueError("Result already exists; prepare a new submission instead of overwriting")
    source = Path.cwd() / "source" if a.sandbox else Path(manifest["source_dir"])
    env = dict(os.environ, PYTHONPATH=str(source), OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONUNBUFFERED="1")
    command = [manifest["python"], "-u", str(source / "scripts/validate_conformal_evolution.py"),
               "--model", case["model"], "--mode", case["mode"],
               "--seed", str(case["seed"]), "--duration-ps", str(manifest["duration_ps"]),
               "--max-step-ps", str(case["max_step_ps"]),
               "--ray-angle", str(case["ray_angle"]),
               "--causal-fraction", str(case["causal_fraction"]),
               "--max-refinements", "0", "--output", str(result_path)]
    provenance = dict(source_commit=manifest["source_commit"], case=case,
                      command=command, python_version=sys.version,
                      started_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    with (destination / "progress.log").open("w") as log:
        process = subprocess.Popen(command, cwd=source, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        code = process.wait()
    completion = dict(returncode=code, finished_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    (destination / "completion.json").write_text(json.dumps(completion, indent=2) + "\n")
    raise SystemExit(code)


def summarize(a):
    output = Path(a.output_dir)
    manifest = json.loads((output / "manifest.json").read_text())
    rows = []
    for case in manifest["cases"]:
        directory = output / case["name"]
        if (directory / 'result').is_dir():
            directory = directory / 'result'
        completion = directory / "completion.json"
        path = directory / "evolution.json"
        row = dict(case=case["name"], status="INCOMPLETE", source_commit=manifest["source_commit"])
        if completion.exists():
            code = json.loads(completion.read_text())["returncode"]
            row.update(returncode=code, status="PROCESS_FAILED" if code else "MISSING_RESULT")
            if path.exists():
                r = json.loads(path.read_text())
                row.update(termination=r["termination"], validation=r["numerical_validation"],
                           status="PASS" if code == 0 and r["numerical_validation"]["status"] == "PASS"
                           else "FAIL")
        rows.append(row)
    summary = dict(scope=manifest["scope"], cases=rows,
                   all_cases_pass=all(r["status"] == "PASS" for r in rows),
                   physical_survival_probability=None,
                   note="Passing both resolutions does not establish observable convergence; compare trajectories and matching events separately.")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if summary["all_cases_pass"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--suite", choices=["remaining", "convergence"], default="remaining")
    p.add_argument("--output-root", default="/data/aalida/QuirkSim/conformal_validation")
    p.add_argument("--submit", action="store_true")
    p = sub.add_parser("run")
    p.add_argument("manifest")
    p.add_argument("case")
    p.add_argument("--sandbox", action="store_true")
    p = sub.add_parser("summarize")
    p.add_argument("output_dir")
    a = parser.parse_args()
    if a.action == "prepare":
        prepare(a)
    elif a.action == "run":
        run(a)
    else:
        raise SystemExit(summarize(a))


if __name__ == "__main__":
    main()
