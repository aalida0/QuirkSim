#!/usr/bin/env python3
"""Prepare reproducible six-channel triplet fixed-order NLO-QCD MG5 cards.

This generates only MG5 command files; it does not run MG5. The internal-PDF
mode is a diagnostic. LHAPDF mode requires an NLO PDF set already installed.
No hidden-color multiplicity is included inside MadGraph.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path


PROCESSES = (
    ("q1+ q1-", "charged_plus_pair"),
    ("q2- q2+", "charged_minus_pair"),
    ("q1+ q0~", "wplus_q1"),
    ("q0 q2+", "wplus_q2"),
    ("q1- q0", "wminus_q1"),
    ("q0~ q2-", "wminus_q2"),
)


def make_cards(
    *,
    model: Path,
    output_dir: Path,
    card_dir: Path,
    accuracy: float = 0.02,
    pdf: str = "internal",
    lhaid: int | None = None,
    seed: int = 7481,
    ebeam_gev: float = 6800,
    with_launch: bool = True,
) -> tuple[Path, Path]:
    if pdf not in ("internal", "lhapdf"):
        raise ValueError("pdf must be 'internal' or 'lhapdf'")
    if (pdf == "lhapdf") != (lhaid is not None):
        raise ValueError("Specify --lhaid only with --pdf lhapdf")
    if lhaid is not None and lhaid <= 0:
        raise ValueError("lhaid must be a positive LHAPDF set ID")
    if accuracy <= 0 or not math.isfinite(accuracy):
        raise ValueError("accuracy must be a finite positive fraction")
    if not (ebeam_gev > 0 and math.isfinite(ebeam_gev)):
        raise ValueError("ebeam_gev must be positive and finite")
    if not 0 < seed < 30081:
        raise ValueError("seed must be in [1,30080] for MadGraph")
    model = model.resolve()
    output_dir = output_dir.resolve()
    card_dir = card_dir.resolve()
    if not (model / "CT_vertices.py").is_file():
        raise FileNotFoundError(f"Missing NLO UFO CT vertices: {model}")
    if output_dir.exists():
        raise FileExistsError(f"Output already exists: {output_dir}; choose a new --output-dir")
    card_dir.mkdir(parents=True, exist_ok=True)
    generate = card_dir / "triplet_six_nlo_generate.mg5"
    launch = card_dir / "triplet_six_nlo_integrate.mg5"
    generate_lines = [
        f"import model {model}",
        "define p = g u u~ d d~ s s~ c c~ b b~",
    ]
    for i, (final_state, _) in enumerate(PROCESSES):
        op = "generate" if i == 0 else "add process"
        generate_lines.append(f"{op} p p > {final_state} QCD=0 QED=2 [QCD] @{i}")
    generate_lines += [f"output {output_dir} -f", "quit"]
    generate.write_text("\n".join(generate_lines) + "\n", encoding="utf-8")

    # Separate stages so a first-time loop-library installation prompt
    # cannot consume 'launch' or parameter-setting commands.
    launch_lines = [
        f"launch {output_dir}",
        "order=NLO",
        "fixed_order=on",
        "shower=off",
        "done",
        f"set iseed {seed}",
        f"set req_acc_FO {accuracy}",
        "set nevents 0",
        f"set ebeam1 {ebeam_gev:g}",
        f"set ebeam2 {ebeam_gev:g}",
        "set fixed_ren_scale False",
        "set fixed_fac_scale False",
        "set dynamical_scale_choice -1",
    ]
    if pdf == "lhapdf":
        launch_lines.extend(["set pdlabel lhapdf", f"set lhaid {lhaid}"])
    launch_lines += ["done", "quit"]
    launch.write_text("\n".join(launch_lines) + "\n", encoding="utf-8")
    return generate, launch


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model", type=Path, default=Path("models/DiracTripletEW_NLO_UFO")
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/pp_dirac_triplet_13p6TeV_NLO")
    )
    parser.add_argument(
        "--card-dir", type=Path, default=Path("results/triplet_nlo_cards")
    )
    parser.add_argument("--accuracy", type=float, default=0.02)
    parser.add_argument("--pdf", choices=["internal", "lhapdf"], default="internal")
    parser.add_argument("--lhaid", type=int)
    parser.add_argument("--seed", type=int, default=7481)
    parser.add_argument("--ebeam-gev", type=float, default=6800)
    args = parser.parse_args()
    generate, launch = make_cards(
        model=args.model, output_dir=args.output_dir,
        card_dir=args.card_dir, accuracy=args.accuracy,
        pdf=args.pdf, lhaid=args.lhaid, seed=args.seed,
        ebeam_gev=args.ebeam_gev,
    )
    print("NLO generation card:", generate)
    print("NLO integration card:", launch)
    if args.pdf == "internal":
        print("CAUTION: internal-PDF run is a numerical pilot, not a final physical prediction.")
    else:
        print("LHAPDF mode: confirm installed NLO set, alpha_s and 5F scheme before interpreting rates.")
    print("Run generation first, then integration with your MG5_aMC executable.")


if __name__ == "__main__":
    main()
