#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cgi_analysis import cgi_bivariate_summary, plot_cgi_bivariate, plot_cgi_bivariate_matlab_style
from ultrasens.utils import write_json


def main() -> None:
    parser = argparse.ArgumentParser(description="Python rewrite of BivariateHistogram_HumanWT_example.m")
    parser.add_argument("--island-files", nargs="+", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-plot", default=None)
    parser.add_argument("--matlab-style-plot", default=None)
    parser.add_argument("--hypo-cutoff", type=float, default=0.2)
    args = parser.parse_args()

    summary = cgi_bivariate_summary(args.island_files, hypo_cutoff=args.hypo_cutoff)
    write_json(args.output_json, summary)
    if args.output_plot:
        plot_cgi_bivariate(summary, args.output_plot)
    if args.matlab_style_plot:
        plot_cgi_bivariate_matlab_style(args.island_files, args.matlab_style_plot, hypo_cutoff=args.hypo_cutoff)
    print(f"wrote CGI bivariate summary to {args.output_json}")


if __name__ == "__main__":
    main()
