#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.intersect import intersect_wgbs_with_density


def main() -> None:
    parser = argparse.ArgumentParser(description="Python replacement for WGBS_CpGIntersect_AllData.command")
    parser.add_argument("--wgbs-bed", required=True)
    parser.add_argument("--density-bed", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    out = intersect_wgbs_with_density(args.wgbs_bed, args.density_bed, args.output)
    print(f"wrote {len(out)} intersected CpG rows to {args.output}")


if __name__ == "__main__":
    main()

