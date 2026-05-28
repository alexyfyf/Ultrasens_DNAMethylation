#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cpg_analysis import CPG_WGBS_DENSITY_COLUMNS, summarize_cpg_density_file
from ultrasens.utils import ensure_dir, maybe_named_columns, optional_pyplot, write_json


def _bivariate_counts(input_path: str | Path, hypo_cutoff: float = 0.2):
    df = pd.read_csv(input_path, sep="\t", header=None)
    df = maybe_named_columns(df, CPG_WGBS_DENSITY_COLUMNS)
    x = df["Density"].to_numpy(float)
    y = df["WGBS"].to_numpy(float)
    if np.nanmax(y) > 0:
        y = y / np.nanmax(y)
    density_vals = np.sort(np.unique(x))
    ygrid = np.linspace(0, np.nanmax(y), 50)
    counts = np.zeros((len(density_vals), len(ygrid)), dtype=int)
    for xi, yi in zip(x, y):
        ix = int(np.searchsorted(density_vals, xi, side="left"))
        iy = min(int(np.searchsorted(ygrid, yi, side="left")), len(ygrid) - 1)
        counts[ix, iy] += 1
    hyper_cutoff = 1 - hypo_cutoff
    hypo_bins = ygrid < hypo_cutoff
    hyper_bins = ygrid > hyper_cutoff
    inter_bins = ~(hypo_bins | hyper_bins)
    totals = np.maximum(counts[:, hypo_bins].sum(axis=1) + counts[:, hyper_bins].sum(axis=1) + counts[:, inter_bins].sum(axis=1), 1)
    return {
        "densityvals": density_vals,
        "ygrid": ygrid,
        "probability": counts / max(len(x), 1),
        "Hypo": counts[:, hypo_bins].sum(axis=1) / totals,
        "Hyper": counts[:, hyper_bins].sum(axis=1) / totals,
        "Inter": counts[:, inter_bins].sum(axis=1) / totals,
    }


def _plot_combined(inputs, labels, summaries, output_path: str | Path, hypo_cutoff: float = 0.2) -> None:
    plt = optional_pyplot()
    if plt is None:
        return
    n = len(inputs)
    fig, axes = plt.subplots(2, n, figsize=(4.2 * n, 7), squeeze=False)
    for idx, (input_path, label) in enumerate(zip(inputs, labels)):
        bivar = _bivariate_counts(input_path, hypo_cutoff=hypo_cutoff)
        density = bivar["densityvals"]
        ygrid = bivar["ygrid"]
        heat = -np.log(np.where(bivar["probability"] > 0, bivar["probability"], np.nan)).T

        ax = axes[0, idx]
        mesh = ax.imshow(
            heat,
            origin="lower",
            aspect="auto",
            extent=[float(density.min()), float(density.max()), float(ygrid.min()), float(ygrid.max())],
        )
        ax.axhline(hypo_cutoff, linestyle="--", color="green", linewidth=1)
        ax.axhline(1 - hypo_cutoff, linestyle="--", color="red", linewidth=1)
        ax.set_title(label)
        ax.set_xlabel("Local CpG Density")
        ax.set_ylabel("Methyl. Frac.")
        fig.colorbar(mesh, ax=ax, label="-log(Probability)")

        summary = summaries[label]
        ax = axes[1, idx]
        ax.plot(summary["densityvals"], summary["Hypo"], "-o", linewidth=2, color="blue", label=f"Hypo < {hypo_cutoff}")
        ax.plot(summary["densityvals"], summary["Hyper"], "-o", linewidth=2, color="red", label=f"Hyper > {1 - hypo_cutoff}")
        ax.plot(summary["densityvals"], summary["Inter"], "-o", linewidth=2, color="black", label="Inter")
        ax.set_title(f"ED50 = {summary['ED50']:.3f}")
        ax.set_xlabel("Local CpG Density")
        ax.set_ylabel("Frac. of CpGs")
        ax.set_ylim(0, 1)
        ax.grid(True)
        if idx == n - 1:
            ax.legend()
    fig.tight_layout()
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch analysis of CpG methylation classes and ED50 versus local CpG density.")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--labels", nargs="+", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--plot", default=None)
    parser.add_argument("--hypo-cutoff", type=float, default=0.2)
    args = parser.parse_args()

    labels = args.labels or [Path(p).stem.replace("_CpGsOnly_Chr1", "") for p in args.inputs]
    if len(labels) != len(args.inputs):
        raise ValueError("--labels must have the same length as --inputs")
    output_dir = ensure_dir(args.output_dir)

    summaries = {}
    for input_path, label in zip(args.inputs, labels):
        summary = summarize_cpg_density_file(
            input_path,
            hypo_cutoff=args.hypo_cutoff,
            output_prefix=output_dir / f"{label}_CpGsOnly_Chr1",
        )
        summaries[label] = summary

    serial = {
        label: {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in summary.items()}
        for label, summary in summaries.items()
    }
    write_json(args.summary_json, serial)
    if args.plot:
        _plot_combined(args.inputs, labels, summaries, args.plot, hypo_cutoff=args.hypo_cutoff)
    for label, summary in summaries.items():
        print(f"{label}: ED50={summary['ED50']:.6g}, HillN={summary['HillN']:.6g}, density points={len(summary['densityvals'])}")


if __name__ == "__main__":
    main()
