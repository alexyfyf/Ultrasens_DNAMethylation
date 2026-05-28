#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cpg_analysis import plot_cpg_matlab_style, plot_cpg_summary, summarize_cpg_density_file
from ultrasens.utils import write_json


def main() -> None:
    parser = argparse.ArgumentParser(description="Python rewrite of ReadPlotData_WT_example.m")
    parser.add_argument("--input", required=True, help="Intersected WGBS/density file: chr,start,end,WGBS,Density")
    parser.add_argument("--output-prefix", required=True, help="Writes .npz and .json DataStruct-compatible summaries.")
    parser.add_argument("--summary-json", default=None)
    parser.add_argument("--plot", default=None)
    parser.add_argument("--matlab-style-plot", default=None, help="Optional ReadPlotData_WT_example.m-style two-panel plot.")
    parser.add_argument("--label", default=None, help="Optional plot title label. Defaults to the input filename stem.")
    parser.add_argument("--hypo-cutoff", type=float, default=0.2)
    args = parser.parse_args()
    summary = summarize_cpg_density_file(args.input, hypo_cutoff=args.hypo_cutoff, output_prefix=args.output_prefix)
    serial = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in summary.items()}
    if args.summary_json:
        write_json(args.summary_json, serial)
    if args.plot:
        plot_cpg_summary(summary, args.plot)
    if args.matlab_style_plot:
        plot_cpg_matlab_style(args.input, summary, args.matlab_style_plot, label=args.label, hypo_cutoff=args.hypo_cutoff)
    print(f"wrote CpG-density DataStruct to {args.output_prefix}.npz and {args.output_prefix}.json")


if __name__ == "__main__":
    main()
