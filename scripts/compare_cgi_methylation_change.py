#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cgi_analysis import compare_cgi_methylation, plot_cgi_methylation_change


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare CGI-level methylation between two aggregate files.")
    parser.add_argument("--file-a", required=True)
    parser.add_argument("--file-b", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--plot", default=None)
    parser.add_argument("--title", default="")
    args = parser.parse_args()
    summary = compare_cgi_methylation(args.file_a, args.file_b)
    summary.to_csv(args.output, index=False)
    if args.plot:
        plot_cgi_methylation_change(summary, args.plot, title=args.title)
    print(f"wrote CGI methylation-change summary to {args.output}")


if __name__ == "__main__":
    main()
