#!/usr/bin/env python3
"""EXPERIMENTAL (unvalidated) MG5 NLO QCD+EW smoke cards.

A hard gate rejects a QCD-only UFO. A *passing* static counterterm gate
is not physics validation. These cards must not be used for cross-section
results until explicit NLO EW UV/IR, gauge, and QED factorization checks.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from audit_triplet_nlo_ew_candidate import inspect_ufo


def prepare(ufo: Path, output: Path, cards: Path, *,
            pdfid: int = 335900, acknowledge_experimental: bool=False):
    if not acknowledge_experimental:
        raise ValueError("Acknowledge this is experimental and *not* a valid EW result")
    result=inspect_ufo(ufo)
    if not result["passes_basic_triplet_ct_gate"]:
        raise ValueError("NLO EW UFO lacks triplet-specific UV/R2 counterterms")
    output=output.resolve()
    cards=cards.resolve()
    if output.exists() or (cards.exists() and any(cards.iterdir())):
        raise FileExistsError("Use a new process and empty cards directory")
    if pdfid != 335900:
        raise ValueError("Smoke test explicitly calibrated only for photon-PDF LHAPDF 335900")
    cards.mkdir(parents=True,exist_ok=True)
    g=cards/"qcdew_charged_pair_generate.mg5"
    l=cards/"qcdew_charged_pair_integrate.mg5"
    g.write_text("\n".join([
        "set complex_mass_scheme True",
        f"import model {ufo.resolve()}",
        # include the photon in initial partons and keep all photon-induced
        # Born terms of the chosen QCD/QED squared coupling order
        "define p = g u u~ d d~ s s~ c c~ b b~ a",
        "generate p p > q1+ q1- QCD=0 QED=2 [QCD QED] @0",
        "add process p p > q2- q2+ QCD=0 QED=2 [QCD QED] @1",
        f"output {output} -f",
        "quit",
    ])+"\n")
    l.write_text("\n".join([
        f"launch {output}",
        "order=NLO","fixed_order=on","shower=off","done",
        "set pdlabel lhapdf",f"set lhaid {pdfid}",
        "set ebeam1 6800","set ebeam2 6800",
        "set req_acc_FO 0.05","set nevents 0",
        "set fixed_ren_scale False","set fixed_fac_scale False",
        "set dynamical_scale_choice -1","done","quit",
    ])+"\n")
    return g,l


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ufo",type=Path,
        default=Path("theory/dirac_triplet_ew/generated/"
                     "DiracTripletEW_QCDEW_CANDIDATE_UFO"))
    p.add_argument("--output",type=Path,
        default=Path("results/pp_dirac_triplet_13p6TeV_QCDEW_EXPERIMENTAL"))
    p.add_argument("--cards",type=Path,
        default=Path("results/triplet_qcdew_EXPERIMENTAL_cards"))
    p.add_argument("--pdfid",type=int,default=335900)
    p.add_argument("--acknowledge-experimental",action="store_true")
    a=p.parse_args()
    g,l=prepare(a.ufo,a.output,a.cards,pdfid=a.pdfid,
        acknowledge_experimental=a.acknowledge_experimental)
    print("EXPERIMENTAL NLO QCD+EW cards; not a validated model:")
    print("Generate:",g)
    print("Integrate:",l)
    print("Required: QED PDF 335900 installed, UV+IR pole cancellation,"
          " soft tests, gauge-check cross sections, SM limit, N_IC audit.")


if __name__=="__main__":
    main()
