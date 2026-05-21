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
    parser.add_argument(
        "--deduplicate-cpg-strands",
        action="store_true",
        help="Collapse adjacent opposite-strand rows that represent the same CpG dyad.",
    )
    parser.add_argument(
        "--deduplicate-method",
        choices=["first", "second", "mean", "coverage-weighted"],
        default="coverage-weighted",
        help=(
            "How to collapse paired CpG-strand rows. 'coverage-weighted' combines methylated/unmethylated "
            "counts before calculating WGBS; 'second' matches the original README convention."
        ),
    )
    args = parser.parse_args()

    out = parse_bismark_coverage(
        args.cov,
        args.output,
        methylation_source=args.methylation_source,
        min_coverage=args.min_coverage,
        sort=not args.no_sort,
        deduplicate_cpg_strands=args.deduplicate_cpg_strands,
        deduplicate_method=args.deduplicate_method,
    )
    print(f"wrote {len(out)} WGBS rows to {args.output}")


if __name__ == "__main__":
    main()
