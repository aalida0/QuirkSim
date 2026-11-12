#!/usr/bin/env python3
"""Write MG5 3.5.1-compatible matched LO/NLO charged-pair production cards.

Only electrically charged pair channels, no charged-neutral associated production.
The default --species both adds q1+ q1- and q2- q2+ as independent final
states. --species q1 or q2 requests just one particular Dirac charged field.

The LO and NLO runs use the same NLO LHAPDF set, masses, central scale,
model, 5-flavour scheme and beam energies: a fixed-PDF-order K factor.
Run LO first, then NLO, to leave the final saved SubProcesses job weights
in the NLO state for scale and PDF uncertainty extraction.

No infracolor degrees are introduced into the UFO. Apply N_IC externally.
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path


CHANNELS = {
    "q1": "q1+ q1-",
    "q2": "q2- q2+",
}


def make_cards(
    model: Path,
    output_dir: Path,
    card_dir: Path,
    *,
    species: str = "both",
    lhaid: int = 331700,
    beam_energy: float = 6800.0,
    accuracy: float = 0.01,
) -> tuple[Path, Path, Path]:
    if species not in {"both", "q1", "q2"}:
        raise ValueError("species must be one of: both, q1, q2")
    if lhaid <= 0:
        raise ValueError("LHAPDF ID must be positive")
    if not 0.0 < accuracy < 1.0:
        raise ValueError("accuracy must be in (0,1)")
    if not math.isfinite(beam_energy) or beam_energy <= 0:
        raise ValueError("beam energy must be finite and positive")

    model, output_dir, card_dir = (
        model.resolve(), output_dir.resolve(), card_dir.resolve()
    )
    if not (model / "CT_vertices.py").is_file():
        raise FileNotFoundError("Build NLO UFO first (missing CT_vertices.py): " + str(model))
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if card_dir.exists() and any(card_dir.iterdir()):
        raise FileExistsError(f"Card directory is not empty: {card_dir}")
    card_dir.mkdir(parents=True, exist_ok=True)

    generate = card_dir / "charged_pair_generate.mg5"
    lo = card_dir / "charged_pair_lo.mg5"
    nlo = card_dir / "charged_pair_nlo.mg5"

    selected = list(CHANNELS.values()) if species == "both" else [CHANNELS[species]]
    gen = [
        f"import model {model}",
        "define p = g u u~ d d~ s s~ c c~ b b~",
    ]
    for idx, pair in enumerate(selected):
        command = "generate" if idx == 0 else "add process"
        gen.append(f"{command} p p > {pair} QCD=0 QED=2 [QCD] @{idx}")
    gen.extend([f"output {output_dir} -f", "quit"])
    generate.write_text("\n".join(gen) + "\n", encoding="utf-8")

    def launch(order: str, seed: int, with_reweights: bool) -> str:
        lines = [
            f"launch {output_dir}",
            f"order={order}",
            "fixed_order=on",
            "shower=off",
            "done",
            "set pdlabel lhapdf",
            f"set lhaid {lhaid}",
            f"set iseed {seed}",
            f"set req_acc_FO {accuracy:g}",
            "set nevents 0",
            f"set ebeam1 {beam_energy:g}",
            f"set ebeam2 {beam_energy:g}",
            "set fixed_ren_scale False",
            "set fixed_fac_scale False",
            "set dynamical_scale_choice -1",
            f"set reweight_scale {'True' if with_reweights else 'False'}",
            f"set reweight_PDF {'True' if with_reweights else 'False'}",
        ]
        if with_reweights:
            lines.extend([
                "set rw_rscale 1.0,0.5,2.0",
                "set rw_fscale 1.0,0.5,2.0",
            ])
        lines.extend(["done", "quit"])
        return "\n".join(lines) + "\n"

    lo.write_text(launch("LO",7482,False), encoding="utf-8")
    nlo.write_text(launch("NLO",7481,True), encoding="utf-8")
    return generate, lo, nlo


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", choices=["both","q1","q2"], default="both")
    ap.add_argument("--model", type=Path, default=Path("models/DiracTripletEW_NLO_UFO"))
    ap.add_argument("--output-dir", type=Path,
        default=Path("results/pp_dirac_triplet_charged_pair_matched_225"))
    ap.add_argument("--card-dir", type=Path,
        default=Path("results/triplet_charged_pair_cards_225"))
    ap.add_argument("--lhaid", type=int, default=331700)
    ap.add_argument("--beam-energy", type=float, default=6800.0)
    ap.add_argument("--accuracy", type=float, default=0.01)
    args=ap.parse_args()
    paths = make_cards(
        model=args.model, output_dir=args.output_dir,
        card_dir=args.card_dir, species=args.species,
        lhaid=args.lhaid,beam_energy=args.beam_energy,
        accuracy=args.accuracy,
    )
    print(f"Charged species: {args.species} (one hidden-color copy)")
    print("Generate:",paths[0])
    print("Run LO  :",paths[1])
    print("Run NLO :",paths[2])
    print("Run these cards in that order with the SAME MadGraph instance.")
    print("Use the final NLO log to extract 7-point/PDF uncertainties.")
    print("Masses are read from model/param_card.dat: check MQ0=MQC=225 GeV.")


if __name__=="__main__":
    main()
