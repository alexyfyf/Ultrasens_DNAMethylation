#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cme import cme_fit_comparison, cme_fit_ssd, plot_cme_fit_comparison_panel
from ultrasens.utils import load_data_struct, write_json


def parse_vector(value: str) -> np.ndarray:
    return np.asarray([float(x) for x in value.split(",")], dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot CME fitted collaborative model versus no-collaboration standard model."
    )
    parser.add_argument("--data-struct", required=True, help=".npz/.json DataStruct from read_plot_data_wt.py.")
    parser.add_argument(
        "--fit-result",
        required=True,
        help=".npz/.json fit result with Parameters and DistLength from fit_cme_methylation.py.",
    )
    parser.add_argument("--output", required=True, help="Output figure path.")
    parser.add_argument("--summary-json", default=None)
    parser.add_argument("--label", default="HUES8")
    parser.add_argument("--n-cpg", type=int, default=27)
    parser.add_argument("--ds", default="2,3,4,5,6,7,8,9,11,13,17,26,110")
    args = parser.parse_args()

    data = load_data_struct(args.data_struct)
    fit = load_data_struct(args.fit_result)
    parameters = np.asarray(fit["Parameters"], dtype=float)
    dist_length = np.asarray(fit["DistLength"], dtype=float)
    ds = [int(x) for x in args.ds.split(",")]

    plot_cme_fit_comparison_panel(
        parameters,
        dist_length,
        data,
        args.output,
        n_cpg=args.n_cpg,
        label=args.label,
        ds=ds,
    )
    if args.summary_json:
        comparison = cme_fit_comparison(parameters, dist_length, data, n_cpg=args.n_cpg, ds=ds)
        write_json(
            args.summary_json,
            {
                "collaborative_ssd": cme_fit_ssd(comparison["collaborative"], data),
                "standard_ssd": cme_fit_ssd(comparison["standard"], data),
                "Parameters": parameters.tolist(),
                "DistLength": dist_length.tolist(),
            },
        )
    print(f"wrote CME fit comparison plot to {args.output}")


if __name__ == "__main__":
    main()
