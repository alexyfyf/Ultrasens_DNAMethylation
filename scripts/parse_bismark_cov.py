#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.wgbs import parse_bismark_coverage


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert Bismark .cov/.cov.gz to processed WGBS BED.")
    parser.add_argument("--cov", required=True, help="Bismark coverage file: .cov, .cov.gz, or equivalent TSV.")
    parser.add_argument("--output", required=True, help="Output 4-column BED: chr,start,end,WGBS.")
    parser.add_argument(
        "--methylation-source",
        choices=["counts", "percentage"],
        default="counts",
        help="Use exact methylated/(methylated+unmethylated) counts or the percentage column.",
    )
    parser.add_argument("--min-coverage", type=int, default=1, help="Drop rows with total read coverage below this value.")
    parser.add_argument("--no-sort", action="store_true", help="Preserve input row order instead of sorting by chr/start/end.")
    args = parser.parse_args()

    out = parse_bismark_coverage(
        args.cov,
        args.output,
        methylation_source=args.methylation_source,
        min_coverage=args.min_coverage,
        sort=not args.no_sort,
    )
    print(f"wrote {len(out)} WGBS rows to {args.output}")


if __name__ == "__main__":
    main()
