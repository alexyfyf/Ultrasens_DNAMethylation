#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cpg_analysis import plot_figure5a_density_panels, summarize_cpg_density_file
from ultrasens.utils import write_json


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot CpG-density methylation summary panels from WGBS/density intersection inputs."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Intersected WGBS/density files with columns: chr start end WGBS Density.",
    )
    parser.add_argument("--labels", nargs="+", default=None, help="Panel labels, one per input.")
    parser.add_argument("--output", required=True, help="Output figure path.")
    parser.add_argument("--summary-json", default=None, help="Optional JSON summary of fitted values.")
    parser.add_argument(
        "--max-neighbors-per-100bp",
        type=float,
        default=50.0,
        help="Conversion from normalized density to # CpG neighbors in 100 bp.",
    )
    parser.add_argument("--hypo-cutoff", type=float, default=0.2)
    args = parser.parse_args()

    input_paths = [Path(p) for p in args.inputs]
    labels = args.labels or [p.stem.replace("_CpGsOnly_Chr1", "") for p in input_paths]
    if len(labels) != len(input_paths):
        raise SystemExit("--labels must have the same number of values as --inputs")

    summaries = [summarize_cpg_density_file(path, hypo_cutoff=args.hypo_cutoff) for path in input_paths]
    plot_figure5a_density_panels(
        summaries,
        labels,
        args.output,
        max_neighbors_per_100bp=args.max_neighbors_per_100bp,
    )

    if args.summary_json:
        serial = {}
        for label, summary in zip(labels, summaries):
            serial[label] = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in summary.items()}
        write_json(args.summary_json, serial)
    print(f"wrote CpG-density methylation summary panel to {args.output}")


if __name__ == "__main__":
    main()
