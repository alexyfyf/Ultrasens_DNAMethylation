#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cpg_density import compute_density_table, load_cpg_positions, write_density_bed


def main() -> None:
    parser = argparse.ArgumentParser(description="Python rewrite of CpGDensity_Calc.m")
    parser.add_argument("--cpg-csv", required=True, help="CSV from seqkit locate with a start column.")
    parser.add_argument("--window", type=int, default=1000)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--output", required=True, help="BED-like output path.")
    parser.add_argument("--start-column", default="start")
    parser.add_argument("--trim-edges-like-matlab", action="store_true")
    args = parser.parse_args()

    positions = load_cpg_positions(args.cpg_csv, args.start_column)
    table = compute_density_table(positions, args.window, args.chrom, trim_edges_like_matlab=args.trim_edges_like_matlab)
    write_density_bed(table, args.output)
    print(f"wrote {len(table)} CpG density rows to {args.output}")


if __name__ == "__main__":
    main()

