#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="End-to-end Python pipeline with placeholder inputs.")
    parser.add_argument("--cpg-csv", required=True, help="PLACEHOLDER: chr1 CpG CSV from seqkit locate.")
    parser.add_argument("--wgbs-bed", required=True, help="PLACEHOLDER: processed WGBS BED with chr,start,end,methyl_fraction.")
    parser.add_argument("--cgi-intersection", required=True, help="PLACEHOLDER: WGBS/CGI intersection table.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--density-window", type=int, default=50)
    parser.add_argument("--n-cpg", type=int, default=27)
    parser.add_argument("--fit-iterations", type=int, default=60)
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    py = sys.executable
    density_bed = out / f"CpGDensities_W{args.density_window}.bed"
    cpg_intersect = out / "WGBS_CpGsOnly_Chr1.bed"
    island = out / "IslandLvl_agg.csv"
    data_prefix = out / "Save_WGBS_CpGsOnly_Chr1"
    fit_prefix = out / "Fit_WGBS_CpGsOnly_Chr1"
    temporal_prefix = out / "Temporal_WGBS_CpGsOnly_Chr1"

    run([py, str(here / "calculate_cpg_density.py"), "--cpg-csv", args.cpg_csv, "--window", str(args.density_window), "--chrom", args.chrom, "--output", str(density_bed)])
    run([py, str(here / "intersect_wgbs_cpg_density.py"), "--wgbs-bed", args.wgbs_bed, "--density-bed", str(density_bed), "--output", str(cpg_intersect)])
    run([py, str(here / "aggregate_cgi_level.py"), "--cgi-intersection", args.cgi_intersection, "--output", str(island)])
    run([py, str(here / "bivariate_histogram_human_wt.py"), "--island-files", str(island), "--output-json", str(out / "cgi_bivariate_summary.json")])
    run([py, str(here / "read_plot_data_wt.py"), "--input", str(cpg_intersect), "--output-prefix", str(data_prefix), "--summary-json", str(out / "cpg_density_summary.json")])
    run([py, str(here / "fit_cme_methylation.py"), "--data-struct", str(data_prefix.with_suffix(".npz")), "--output-prefix", str(fit_prefix), "--n-cpg", str(args.n_cpg), "--iterations", str(args.fit_iterations)])
    run([py, str(here / "run_temporal_cme_example.py"), "--fit-result", str(fit_prefix.with_suffix(".npz")), "--output-prefix", str(temporal_prefix), "--n-cpg", str(args.n_cpg)])


if __name__ == "__main__":
    main()

