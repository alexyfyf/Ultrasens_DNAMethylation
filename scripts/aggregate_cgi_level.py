#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cgi_analysis import aggregate_cgi_level


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate WGBS/CGI intersections to CGI-level means.")
    parser.add_argument("--cgi-intersection", required=True, help="Columns: chr,start,end,WGBS,CGIno,CpGNum,CGIlen")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    table = aggregate_cgi_level(args.cgi_intersection, args.output)
    print(f"wrote {len(table)} CGI-level rows to {args.output}")


if __name__ == "__main__":
    main()

