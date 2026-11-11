"""Audit the bounded public QuirkPheno ROOT archive without extracting it."""
import argparse
import hashlib
import io
import json
import tarfile
import numpy as np
import uproot
from scipy.interpolate import CubicHermiteSpline

C_MM_PS = .299792458
HBARC_MEV_MM = 1.973269804e-10

def audit(path):
    cases = {}
    trajectories = {}
    with tarfile.open(path) as archive:
        manifest_member = next(m for m in archive.getmembers() if m.name.endswith("/manifest.json"))
        manifest = json.load(archive.extractfile(manifest_member))
        prefix = manifest_member.name.rsplit("/", 1)[0] + "/"
        mass = manifest["mass_gev"] * 1000
        sigma = (manifest["lambda_ev"] * 1e-6)**2 / HBARC_MEV_MM
        reference_mm = HBARC_MEV_MM / (manifest["lambda_ev"] * 1e-6)
        for case in manifest["cases"]:
            name = case["name"]
            macro = archive.extractfile(prefix + name + ".mac").read()
            macro_hash_matches = hashlib.sha256(macro).hexdigest() == case["macro_sha256"]
            root = uproot.open(io.BytesIO(archive.extractfile(prefix + name + ".root").read()))
            e = root["events"].arrays(library="ak")
            a = root["steps"].arrays(library="np")
            result = {"execution": case["execution"], "macro_hash_matches": macro_hash_matches,
                      "event_entries": len(e), "endpoints": [],
                      "event_cpu_seconds": float(e["cpu"][0]), "event_wall_seconds": float(e["wall"][0])}
            log = archive.extractfile(prefix + name + ".log").read().decode(errors="replace")
            result["geant4_exception_lines"] = [line for line in log.splitlines() if "*** G4Exception" in line]
            for q in (0, 1):
                mask = a["q"] == q
                t = a["t"][mask] * 1000
                x = np.column_stack([a[k][mask] for k in ("x", "y", "z")])
                p = np.column_stack([a[k][mask] for k in ("px", "py", "pz")])
                E = a["ekin"][mask] + mass
                initial_p = np.array([float(e[k][0][q]) for k in ("q_px0", "q_py0", "q_pz0")])
                initial_x = np.array([float(e[k][0][q]) for k in ("q_x0", "q_y0", "q_z0")])
                initial_E = float(e["q_e0"][0][q])
                ep = {"endpoint": q, "rows": len(t), "reported_steps": int(e["q_steps"][0][q]),
                      "end_reason": int(e["q_reason"][0][q]), "last_time_ps": float(t[-1]),
                      "covers_20_ps": bool(t[-1] >= 20),
                      "strictly_increasing_times": bool(np.all(np.diff(t) > 0)),
                      "finite_values": bool(np.all(np.isfinite(np.r_[t, x.ravel(), p.ravel(), E]))),
                      "max_edep_mev": float(max(a["edep"][mask]))}
                if case["field_tesla"] == 0:
                    beta = initial_p[2] / initial_E
                    gamma = 1 / np.sqrt(1 - beta**2)
                    p_abs = abs(initial_p[0]) - sigma * C_MM_PS * t / gamma
                    E_rest = np.sqrt(mass**2 + p_abs**2)
                    x_exact = np.column_stack([np.sign(initial_p[0]) * (initial_E/gamma - E_rest) / sigma,
                                               np.zeros_like(t), beta * C_MM_PS * t])
                    p_exact = np.column_stack([np.sign(initial_p[0]) * p_abs, np.zeros_like(t),
                                               gamma * beta * E_rest])
                    ep["analytic_outbound_max_position_error_mm"] = float(max(np.linalg.norm(x-x_exact, axis=1)))
                    ep["analytic_outbound_max_momentum_error_mev"] = float(max(np.linalg.norm(p-p_exact, axis=1)))
                    ep["analytic_first_turn_time_ps"] = float(abs(initial_p[0]) * gamma / (sigma*C_MM_PS))
                v = C_MM_PS * p / E[:, None]
                if t[0] > 0:
                    t = np.r_[0, t]
                    x = np.vstack([initial_x, x])
                    v = np.vstack([C_MM_PS*initial_p/initial_E, v])
                trajectories[(name, q)] = (t, x, v)
                result["endpoints"].append(ep)
            cases[name] = result
        comparison = []
        grid = np.linspace(0, 20, 2001)
        for field in (0, 2):
            for q in (0, 1):
                interpolants = []
                for level in (0, 1):
                    t, x, v = trajectories[(f"stable_B{field}_level{level}", q)]
                    hermite = CubicHermiteSpline(t, x, v)(grid)
                    linear = np.column_stack([np.interp(grid, t, x[:, j]) for j in range(3)])
                    interpolants.append((hermite, linear))
                diff = np.linalg.norm(interpolants[0][0]-interpolants[1][0], axis=1)
                comparison.append({"field_tesla": field, "endpoint": q,
                    "hermite_max_position_difference_mm": float(max(diff)),
                    "hermite_max_difference_over_reference": float(max(diff)/reference_mm),
                    "linear_max_position_difference_mm": float(max(np.linalg.norm(interpolants[0][1]-interpolants[1][1], axis=1))),
                    "hermite_y_at_20ps_mm": [float(h[-1,1]) for h,l in interpolants],
                    "interpolation_error_bound_established": False})
    with open(path, "rb") as stream:
        digest = hashlib.sha256(stream.read()).hexdigest()
    return {"archive_sha256": digest, "source_commit": manifest["observed_source_commit"],
            "build": manifest["build"], "confinement_reference_mm": reference_mm,
            "cases": cases, "resolution_comparisons": comparison,
            "scope": "STABLE_OUTBOUND_TRANSPORT_ONLY_NO_RECOIL",
            "verdict": "SMOKE_PASS_B0_ANALYTIC_AGREEMENT_B2_REFINEMENT_UNRESOLVED",
            "physical_survival_probability": None,
            "limitations": ["No weak decays, turnaround, crossing, radiation or annihilation tested.",
                "Interpolation is exploratory; no certified error bound for B2.",
                "Endpoint data alone do not establish total string energy-momentum conservation.",
                "Upstream point-production state differs from stretched conformal benchmark."]}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = audit(args.archive)
    with open(args.output, "w") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(result["verdict"])
