#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cme import temporal_model
from ultrasens.cpg_analysis import fit_switch_values, fit_switch_values_temporal
from ultrasens.cme import model_curves
from ultrasens.utils import ensure_dir, load_data_struct, optional_pyplot, save_data_struct, write_json


def moving_average(y: np.ndarray, window: int = 5) -> np.ndarray:
    if len(y) < window:
        return y
    pad = window // 2
    padded = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(window) / window
    return np.convolve(padded, kernel, mode="valid")


def plot_temporal_outputs(model: dict[str, np.ndarray], base_vals: np.ndarray, all_vals: np.ndarray, plot_prefix: str | Path) -> None:
    plt = optional_pyplot()
    if plt is None:
        print("matplotlib is not installed; skipped temporal plots")
        return
    plot_prefix = Path(plot_prefix)
    ensure_dir(plot_prefix.parent)
    density = model["densityvals"]
    mean = model["Mean"]

    fig, ax = plt.subplots(figsize=(2.75, 2.0))
    for tind in range(max(0, mean.shape[0] - 1)):
        color = "black" if tind == 0 else (0.6, 0.6, 0.6)
        ax.plot(density, moving_average(mean[tind]), "-", linewidth=2, color=color)
    ax.set_xlabel("Local CpG Density")
    ax.set_ylabel("Mean CpG Methylation")
    ax.set_xlim(0, 0.5)
    fig.tight_layout()
    fig.savefig(plot_prefix.with_name(plot_prefix.name + "_mean_methylation.png"), dpi=300)
    plt.close(fig)

    km_list = all_vals[:, 4] - base_vals[4]
    nm_list = all_vals[:, 5] - base_vals[5]
    fig, ax = plt.subplots(figsize=(1.1, 1.1))
    if len(km_list) > 1:
        ax.scatter(km_list[:-1], nm_list[:-1], s=20, edgecolors=(0.6, 0.6, 0.6), facecolors=(0.6, 0.6, 0.6), alpha=0.35)
    ax.scatter([km_list[0]], [nm_list[0]], s=20, edgecolors="black", facecolors="black", alpha=0.35)
    ax.grid(True)
    ax.set_xlim(-0.1, 0.1)
    ax.set_ylim(-2, 2)
    ax.set_xlabel("K*")
    ax.set_ylabel("n*")
    fig.tight_layout()
    fig.savefig(plot_prefix.with_name(plot_prefix.name + "_inset.png"), dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Python rewrite of RunTemporalCMEExample/Plot_ExampleTemporal.m")
    parser.add_argument("--fit-result", required=True, help=".npz/.json from fit_cme_methylation.py")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--n-cpg", type=int, default=27)
    parser.add_argument("--tspan", default="0,20,40,60,100,200")
    parser.add_argument("--k3-scale", type=float, default=0.88)
    parser.add_argument("--plot-prefix", default=None, help="Optional output prefix for Plot_ExampleTemporal-style PNGs.")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    fit = load_data_struct(args.fit_result)
    params1 = fit["Parameters"]
    dist1 = fit["DistLength"]
    params2 = params1.copy()
    params2[2] = params1[2] * args.k3_scale
    tspan = np.asarray([float(x) for x in args.tspan.split(",")], dtype=float)
    model = temporal_model(params1, params2, dist1, dist1.copy(), tspan, n_cpg=args.n_cpg)
    base_model = model_curves(params1, dist1, n_cpg=args.n_cpg)
    base_vals = fit_switch_values(base_model)
    all_vals = fit_switch_values_temporal(model)
    km_list = all_vals[:, 4] - base_vals[4]
    nm_list = all_vals[:, 5] - base_vals[5]
    save_data_struct(
        args.output_prefix,
        **model,
        tspan=tspan,
        BaseSwitchValues=base_vals,
        TemporalSwitchValues=all_vals,
        KStar=km_list,
        NStar=nm_list,
    )
    write_json(
        Path(args.output_prefix).with_suffix(".switch_summary.json"),
        {
            "columns": ["ED50", "slope", "KFD", "nFD", "KFM", "nFM"],
            "BaseSwitchValues": base_vals.tolist(),
            "TemporalSwitchValues": all_vals.tolist(),
            "KStar": km_list.tolist(),
            "NStar": nm_list.tolist(),
        },
    )
    if not args.no_plots:
        plot_temporal_outputs(model, base_vals, all_vals, args.plot_prefix or args.output_prefix)
    print(f"wrote temporal CME output to {args.output_prefix}.npz/.json")


if __name__ == "__main__":
    main()
