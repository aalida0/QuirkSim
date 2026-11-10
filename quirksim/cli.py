"""CLI for independent quirk-string research prototype."""
import argparse
import json
from pathlib import Path
from .dynamics import Parameters, propagate


def main():
    parser = argparse.ArgumentParser(description="Standalone weakly decaying quirk pair")
    parser.add_argument("--model", choices=["triplet","doublet"], default="triplet")
    parser.add_argument("--mass", type=float, default=225)
    parser.add_argument("--lambda-ev", type=float, default=774)
    parser.add_argument("--beta", type=float, default=0.4)
    parser.add_argument("--impact-cm", type=float, default=0)
    parser.add_argument("--duration-ps", type=float, default=500)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--mode", choices=["stable","no_recoil","recoil"], default="recoil")
    parser.add_argument("--output", default="result.json")
    a = parser.parse_args()
    p = Parameters(representation=a.model, mass_neutral_gev=a.mass, lambda_ev=a.lambda_ev,
                   beta=a.beta, impact_cm=a.impact_cm, duration_ps=a.duration_ps,
                   seed=a.seed, mode=a.mode)
    result = propagate(p)
    Path(a.output).write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="trace_time_ps_separation_cm"},indent=2))


if __name__ == "__main__":
    main()
