#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cpg_analysis import (
    plot_switch_parameter_panels,
    summarize_cpg_density_file,
    switch_record_from_summary,
)


REQUIRED_TABLE_COLUMNS = {"label", "HillK", "HillN", "DirectED50", "DirectSlope"}


def read_metric_table(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    missing = REQUIRED_TABLE_COLUMNS - set(table.columns)
    if missing:
        raise SystemExit(f"{path} is missing required columns: {', '.join(sorted(missing))}")
    if "cohort" not in table.columns:
        table["cohort"] = "All"
    if "ko_type" not in table.columns:
        table["ko_type"] = "None"
    return table


def records_from_inputs(args: argparse.Namespace) -> pd.DataFrame:
    paths = [Path(p) for p in args.inputs]
    labels = args.labels or [p.stem.replace("_CpGsOnly_Chr1", "") for p in paths]
    cohorts = args.cohorts or ["All"] * len(paths)
    ko_types = args.ko_types or ["None"] * len(paths)
    if not (len(paths) == len(labels) == len(cohorts) == len(ko_types)):
        raise SystemExit("--inputs, --labels, --cohorts, and --ko-types must have matching lengths")
    rows = []
    for path, label, cohort, ko_type in zip(paths, labels, cohorts, ko_types):
        summary = summarize_cpg_density_file(path, hypo_cutoff=args.hypo_cutoff)
        rows.append(switch_record_from_summary(summary, label=label, cohort=cohort, ko_type=ko_type))
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Fig. 5B/C-style switch parameter plots: Hill n vs K and direct slope vs ED50."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--table",
        help=(
            "CSV with columns label,HillK,HillN,DirectED50,DirectSlope and optional cohort,ko_type, "
            "plus optional symmetric/asymmetric error columns such as HillKErr or HillKErrLow/HillKErrHigh."
        ),
    )
    source.add_argument(
        "--inputs",
        nargs="+",
        help="Intersected WGBS/density files with columns: chr start end WGBS Density.",
    )
    parser.add_argument("--labels", nargs="+", default=None, help="Point labels, one per input.")
    parser.add_argument("--cohorts", nargs="+", default=None, help="Panel group labels, one per input.")
    parser.add_argument("--ko-types", nargs="+", default=None, help="KO color labels, one per input.")
    parser.add_argument("--cohort-order", nargs="+", default=None)
    parser.add_argument("--output", required=True, help="Output figure path.")
    parser.add_argument("--summary-csv", default=None, help="Optional CSV of plotted point values.")
    parser.add_argument("--hypo-cutoff", type=float, default=0.2)
    args = parser.parse_args()

    records = read_metric_table(args.table) if args.table else records_from_inputs(args)
    plot_switch_parameter_panels(records, args.output, cohort_order=args.cohort_order)
    if args.summary_csv:
        Path(args.summary_csv).parent.mkdir(parents=True, exist_ok=True)
        records.to_csv(args.summary_csv, index=False)
    print(f"wrote switch parameter plot to {args.output}")


if __name__ == "__main__":
    main()
