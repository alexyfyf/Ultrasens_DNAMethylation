#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.intersect import intersect_wgbs_with_density


DATASETS = [
    ("GSM1112841_HUES64WT_WGBS_proc.bed", "HUES64WT_CpGsOnly_Chr1.bed"),
    ("GSM1545002_DNMT3A_KO_Early_proc.bed", "HUES64_DNMT3Ako_early_CpGsOnly_Chr1.bed"),
    ("GSM1545003_DNMT3B_KO_Early_proc.bed", "HUES64_DNMT3Bko_early_CpGsOnly_Chr1.bed"),
    ("GSM1545005_DNMT3A_KO_Late_proc.bed", "HUES64_DNMT3Ako_late_CpGsOnly_Chr1.bed"),
    ("GSM1545006_DNMT3B_KO_Late_proc.bed", "HUES64_DNMT3Bko_late_CpGsOnly_Chr1.bed"),
    ("GSM1545004_DKO_Early_proc.bed", "HUES64_DNMT3_dko_early_CpGsOnly_Chr1.bed"),
    ("GSM1545007_DKO_Late_proc.bed", "HUES64_DNMT3_dko_late_CpGsOnly_Chr1.bed"),
    ("GSM3618718_HUES8_WT_WGBS_proc.bed", "HUES8WT_CpGsOnly_Chr1.bed"),
    ("GSM4458672_WGBS_HUES8_DKO_P6_proc.bed", "HUES8_DKO_P6_CpGsOnly_Chr1.bed"),
    ("GSM3618720_HUES8_TKO_WGBS_proc.bed", "HUES8_TKO_CpGsOnly_Chr1.bed"),
    ("GSM3618719_HUES8_QKO_WGBS_proc.bed", "HUES8_QKO_CpGsOnly_Chr1.bed"),
    ("GSM3618721_HUES8_PKO_WGBS_proc.bed", "HUES8_PKO_P0_CpGsOnly_Chr1.bed"),
    ("GSM3662266_HUES8_PKO_P6_WGBS_proc.bed", "HUES8_PKO_P6_CpGsOnly_Chr1.bed"),
    ("GSM4458671_WGBS_HUES8_PKO_P20_proc.bed", "HUES8_PKO_P20_CpGsOnly_Chr1.bed"),
    ("GSM432687_IMR90_WGBS_proc.bed", "IMR90WT_CpGsOnly_Chr1.bed"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Intersect the repository's known WGBS dataset filenames with CpG density.")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--density-bed", default=None, help="Defaults to <data-dir>/CpGDensities_W50.bed.")
    parser.add_argument("--output-dir", default=None, help="Defaults to <data-dir>.")
    parser.add_argument("--require-all", action="store_true", help="Fail if any original input file is missing.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir) if args.output_dir else data_dir
    density_bed = Path(args.density_bed) if args.density_bed else data_dir / "CpGDensities_W50.bed"
    if not density_bed.exists():
        raise FileNotFoundError(f"Missing density BED: {density_bed}")

    missing = []
    for input_name, output_name in DATASETS:
        input_path = data_dir / input_name
        output_path = output_dir / output_name
        if not input_path.exists():
            missing.append(input_name)
            print(f"missing {input_name}; skipped")
            continue
        out = intersect_wgbs_with_density(input_path, density_bed, output_path)
        print(f"wrote {len(out)} rows to {output_path}")

    if missing and args.require_all:
        missing_list = "\n".join(f"- {name}" for name in missing)
        raise SystemExit(f"Missing {len(missing)} required inputs:\n{missing_list}")


if __name__ == "__main__":
    main()
