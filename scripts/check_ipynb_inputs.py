#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate data/ inputs used by Analysis of Human ESCs notebook.")
    parser.add_argument("--wgbs-cgi", default="data/HUES64WT_WGBS_CGI_int.bed")
    parser.add_argument("--cgi-annotation", default="data/CpGIsAnn_hg19_dec14.bed")
    parser.add_argument("--output-dir", default="data/ipynb_run_check")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.wgbs_cgi, sep="\t", header=None)
    df.columns = ["Chrom", "CpGstart", "CpGEnd", "WGBS", "CpGNum", "CGIlen", "CGIno"]

    df_hol = pd.read_csv(args.cgi_annotation, sep="\t", header=None)
    df_hol.columns = ["Chrom", "CpGstart", "CpGEnd", "CpGNum", "CGIlen", "CGIno"]

    island = df.groupby("CGIno")[["WGBS", "CpGNum", "CGIlen"]].mean().reset_index()
    counts = df["CGIno"].value_counts().rename_axis("CGIno").reset_index(name="CpGCounts")
    coverage = island.join(counts.set_index("CGIno"), on="CGIno")
    coverage["CpGCoverage"] = coverage["CpGCounts"] / coverage["CpGNum"]

    # These are the intended exports referenced by the notebook's final cell.
    df.to_csv(out / "CGIs_WGBS_Raw.csv", index=False)
    island[["CGIno", "CpGNum", "WGBS"]].to_csv(out / "CpGvsWGBS_mean.csv", index=False)
    island[["CGIno", "CGIlen", "WGBS"]].to_csv(out / "CGIlengthvWGBS_mean.csv", index=False)
    island.to_csv(out / "IslandLvl_agg_HUES64WT.csv", index=False)
    coverage.to_csv(out / "CGI_coverage_HUES64WT.csv", index=False)

    x = df["CGIlen"].to_numpy(dtype=float)
    y = df["CpGNum"].to_numpy(dtype=float)
    coef, intercept = np.polyfit(x, y, 1)
    pearson = np.corrcoef(x, y)[0, 1]
    summary = {
        "wgbs_rows": int(len(df)),
        "annotation_rows": int(len(df_hol)),
        "number_CGIs_reported_from_UCSC": int(df_hol["CGIno"].max()),
        "number_CGIs_with_any_WGBS_values": int(df["CGIno"].nunique()),
        "island_level_rows": int(len(island)),
        "CpGNum_mean": float(df["CpGNum"].mean()),
        "CpGNum_std": float(df["CpGNum"].std()),
        "CGIlen_mean": float(df["CGIlen"].mean()),
        "CGIlen_std": float(df["CGIlen"].std()),
        "pearson_CGIlen_CpGNum": float(pearson),
        "linear_coef_CpGNum_vs_CGIlen": float(coef),
        "linear_intercept_CpGNum_vs_CGIlen": float(intercept),
    }
    text = "\n".join(f"{key}={value}" for key, value in summary.items()) + "\n"
    (out / "summary.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()

