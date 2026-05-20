#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.fit import fit_example_cme_parameters
from ultrasens.utils import load_data_struct, save_data_struct, write_json


def parse_vector(value: str) -> list[float]:
    return [float(x) for x in value.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser(description="Python rewrite of ExampleCMEFit/CallFitting.m")
    parser.add_argument("--data-struct", required=True, help=".mat/.npz/.json containing DataStruct fields.")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--n-cpg", type=int, default=27)
    parser.add_argument("--iterations", type=int, default=60)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--parameters",
        default="0.1,0.1,4,0.14,0,3.5,0,0.73,0.06,0,25,0,0.0416666667,0.0208333333",
    )
    parser.add_argument(
        "--dist-length",
        default="0,0,0,0,33.63,33.63,33.63,33.63,43.38,43.38,43.38,43.38",
    )
    args = parser.parse_args()

    data = load_data_struct(args.data_struct)
    parameters = np.asarray(parse_vector(args.parameters), dtype=float)
    dist_length = np.asarray(parse_vector(args.dist_length), dtype=float)
    best, ssd, fit_params = fit_example_cme_parameters(
        parameters,
        dist_length,
        data,
        n_cpg=args.n_cpg,
        iterations=args.iterations,
        seed=args.seed,
    )
    save_data_struct(args.output_prefix, Parameters=fit_params, DistLength=dist_length, BestFitVector=best, SSD=np.asarray([ssd]))
    write_json(
        Path(args.output_prefix).with_suffix(".summary.json"),
        {"SSD": ssd, "BestFitVector": best.tolist(), "Parameters": fit_params.tolist(), "DistLength": dist_length.tolist()},
    )
    print(f"wrote example CME fit result to {args.output_prefix}.npz/.json")


if __name__ == "__main__":
    main()
