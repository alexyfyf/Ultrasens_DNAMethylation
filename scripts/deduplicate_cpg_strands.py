#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.wgbs import deduplicate_adjacent_cpg_strands


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collapse adjacent opposite-strand WGBS rows that represent the same CpG dyad."
    )
    parser.add_argument("--input", required=True, help="4-column processed WGBS BED: chr start end WGBS.")
    parser.add_argument("--output", required=True, help="Deduplicated 4-column WGBS BED.")
    parser.add_argument(
        "--method",
        choices=["first", "second", "mean"],
        default="second",
        help="How to collapse paired rows. 'second' matches the original README's NR %% 2 == 0 behavior.",
    )
    parser.add_argument(
        "--max-methylation-diff",
        type=float,
        default=None,
        help="Only collapse adjacent rows when their methylation fractions differ by at most this value.",
    )
    parser.add_argument("--no-sort", action="store_true", help="Preserve input order instead of sorting first.")
    args = parser.parse_args()

    out = deduplicate_adjacent_cpg_strands(
        args.input,
        args.output,
        method=args.method,
        sort=not args.no_sort,
        max_methylation_diff=args.max_methylation_diff,
    )
    print(f"wrote {len(out)} deduplicated WGBS rows to {args.output}")


if __name__ == "__main__":
    main()
