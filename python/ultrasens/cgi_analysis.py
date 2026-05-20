from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import ensure_dir, maybe_named_columns, optional_pyplot


CGI_COLUMNS = ["chr", "start", "end", "WGBS", "CpGNum", "CGIlen", "CGIno"]


def aggregate_cgi_level(input_path: str | Path, output_path: str | Path) -> pd.DataFrame:
    """Python version of the notebook's island-level aggregation."""
    df = pd.read_csv(input_path, sep=None, engine="python")
    df = maybe_named_columns(df, CGI_COLUMNS)
    grouped = df.groupby("CGIno", as_index=False)[["WGBS", "CpGNum", "CGIlen"]].mean()
    ensure_dir(Path(output_path).parent)
    grouped.to_csv(output_path, index=False)
    return grouped


def cgi_bivariate_summary(
    island_files: list[str | Path],
    hypo_cutoff: float = 0.2,
    nbins: int = 50,
) -> dict[str, dict[str, float | list[float]]]:
    out: dict[str, dict[str, float | list[float]]] = {}
    hyper_cutoff = 1 - hypo_cutoff
    for file in island_files:
        df = pd.read_csv(file)
        df = maybe_named_columns(df, ["CGIno", "WGBS", "CpGNum", "CGIlen"])
        x = df["CGIlen"].to_numpy(float)
        y = df["WGBS"].to_numpy(float)
        bins = np.percentile(x, np.linspace(0, 100, nbins))
        ygrid = np.arange(0, 1.0001, 0.05)
        bivariate = np.zeros((len(bins), len(ygrid)), dtype=int)
        for xi, yi in zip(x, y):
            x_inds = np.where(xi < bins)[0]
            y_inds = np.where(yi < ygrid)[0]
            x_ind = int(x_inds[0]) if x_inds.size else bivariate.shape[0] - 1
            y_ind = int(y_inds[0]) if y_inds.size else bivariate.shape[1] - 1
            bivariate[x_ind, y_ind] += 1
        hypo_bins = np.where(ygrid < hypo_cutoff)[0]
        hyper_bins = np.where(ygrid > hyper_cutoff)[0]
        inter_bins = np.setdiff1d(np.arange(len(ygrid)), np.r_[hypo_bins, hyper_bins])
        hypo = bivariate[1:-1][:, hypo_bins].sum(axis=1)
        hyper = bivariate[1:-1][:, hyper_bins].sum(axis=1)
        inter = bivariate[1:-1][:, inter_bins].sum(axis=1)
        xax = bins[1:-1]
        arr = np.column_stack([xax, hypo, hyper, inter]).astype(float)
        cross_idx = int(np.argmin(np.abs(hypo - hyper)))
        crossover = float(xax[cross_idx])
        out[str(file)] = {
            "cgi_crossover_bp": crossover,
            "length_bins": arr[:, 0].tolist(),
            "hypo_counts": arr[:, 1].tolist(),
            "hyper_counts": arr[:, 2].tolist(),
            "intermediate_counts": arr[:, 3].tolist(),
        }
    return out


def compare_cgi_methylation(
    file_a: str | Path,
    file_b: str | Path,
    diff_bins=(-1, -0.6, -0.2, 0.2, 0.6, 1),
    percentile_step: int = 2,
) -> pd.DataFrame:
    """Match IndividualCGIMethylationChange: pair CGIs by CpGNum/CGIlen and nearest CGIno."""
    a = maybe_named_columns(pd.read_csv(file_a), ["CGIno", "WGBS", "CpGNum", "CGIlen"])
    b = maybe_named_columns(pd.read_csv(file_b), ["CGIno", "WGBS", "CpGNum", "CGIlen"])
    rows = []
    for _, row in a.iterrows():
        candidates = b[(b["CpGNum"] == row["CpGNum"]) & (b["CGIlen"] == row["CGIlen"])]
        if candidates.empty:
            continue
        best = candidates.iloc[(candidates["CGIno"] - row["CGIno"]).abs().argmin()]
        rows.append((row["CGIno"], row["CGIlen"], row["WGBS"] - best["WGBS"]))
    paired = pd.DataFrame(rows, columns=["CGIno", "CGIlen", "DeltaMeth"])
    if paired.empty:
        return paired
    perc = np.arange(0, 100 + percentile_step, percentile_step)
    len_bins = np.unique(np.percentile(paired["CGIlen"], perc))
    paired["LengthBin"] = pd.cut(paired["CGIlen"], len_bins, include_lowest=True, duplicates="drop")
    records = []
    for length_bin, group in paired.groupby("LengthBin", observed=True):
        counts, _ = np.histogram(group["DeltaMeth"], bins=diff_bins)
        denom = max(len(group), 1)
        record = {
            "LengthBin": str(length_bin),
            "MeanLength": group["CGIlen"].mean(),
            "MeanDiff": group["DeltaMeth"].mean(),
            "MedianDiff": group["DeltaMeth"].median(),
            "Prc25": group["DeltaMeth"].quantile(0.25),
            "Prc75": group["DeltaMeth"].quantile(0.75),
        }
        for i, count in enumerate(counts):
            record[f"Class{i + 1}_Pct"] = count / denom * 100
        records.append(record)
    return pd.DataFrame(records)


def plot_cgi_bivariate(summary: dict, output_path: str | Path) -> None:
    plt = optional_pyplot()
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    for name, vals in summary.items():
        ax.plot(vals["length_bins"], vals["hypo_counts"], label=f"{Path(name).stem} hypo")
        ax.plot(vals["length_bins"], vals["hyper_counts"], label=f"{Path(name).stem} hyper")
        ax.plot(vals["length_bins"], vals["intermediate_counts"], label=f"{Path(name).stem} inter")
    ax.set_xlabel("CGI Length")
    ax.set_ylabel("No. of Islands")
    ax.legend()
    fig.tight_layout()
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path)
    plt.close(fig)


def plot_cgi_bivariate_matlab_style(
    island_files: list[str | Path],
    output_path: str | Path,
    hypo_cutoff: float = 0.2,
) -> None:
    """Generate a compact version of BivariateHistogram_HumanWT_example.m."""
    plt = optional_pyplot()
    if plt is None:
        return
    hyper_cutoff = 1 - hypo_cutoff
    n = len(island_files)
    fig, axes = plt.subplots(2, n, figsize=(4.2 * n, 7), squeeze=False)
    for j, file in enumerate(island_files):
        df = pd.read_csv(file)
        df = maybe_named_columns(df, ["CGIno", "WGBS", "CpGNum", "CGIlen"])
        x = df["CGIlen"].to_numpy(float)
        y = df["WGBS"].to_numpy(float)
        xgrid = np.linspace(0, max(x), 500)
        ygrid = np.linspace(0, max(y), 100)
        hist, xedges, yedges = np.histogram2d(x, y, bins=[xgrid, ygrid])
        prob = hist / max(len(x), 1)
        heat = -np.log(np.where(prob > 0, prob, np.nan)).T
        ax = axes[0, j]
        mesh = ax.pcolormesh(xedges, yedges, heat, shading="auto")
        ax.axhline(hypo_cutoff, linestyle="--", color="green", linewidth=1)
        ax.axhline(hyper_cutoff, linestyle="--", color="red", linewidth=1)
        ax.set_xlim(205, 4500)
        ax.set_xlabel("CGI Length (bp)")
        ax.set_ylabel("Methyl. Frac.")
        ax.set_title(Path(file).stem)
        fig.colorbar(mesh, ax=ax, label="-log(Probability)")

        nbins = 50
        bins = np.percentile(x, np.linspace(0, 100, nbins))
        ygrid3 = np.arange(0, 1.0001, 0.05)
        bivariate = np.zeros((len(bins), len(ygrid3)), dtype=int)
        for xi, yi in zip(x, y):
            x_inds = np.where(xi < bins)[0]
            y_inds = np.where(yi < ygrid3)[0]
            x_ind = int(x_inds[0]) if x_inds.size else bivariate.shape[0] - 1
            y_ind = int(y_inds[0]) if y_inds.size else bivariate.shape[1] - 1
            bivariate[x_ind, y_ind] += 1
        hypo_bins = np.where(ygrid3 < hypo_cutoff)[0]
        hyper_bins = np.where(ygrid3 > hyper_cutoff)[0]
        inter_bins = np.setdiff1d(np.arange(len(ygrid3)), np.r_[hypo_bins, hyper_bins])
        arr = np.column_stack(
            [
                bins[1:-1],
                bivariate[1:-1][:, hypo_bins].sum(axis=1),
                bivariate[1:-1][:, hyper_bins].sum(axis=1),
                bivariate[1:-1][:, inter_bins].sum(axis=1),
            ]
        ).astype(float)
        ax = axes[1, j]
        if arr.size:
            ax.plot(arr[:, 0], arr[:, 1], "-o", linewidth=2, color="blue", label=f"Hypo < {hypo_cutoff}")
            ax.plot(arr[:, 0], arr[:, 2], "-o", linewidth=2, color="red", label=f"Hyper > {hyper_cutoff}")
            ax.plot(arr[:, 0], arr[:, 3], "-o", linewidth=2, color="black", label="Inter")
            cross_idx = int(np.argmin(np.abs(arr[:, 1] - arr[:, 2])))
            ax.set_title(f"CGI-Cross = {arr[cross_idx, 0]:.1f}bp")
        ax.set_xlim(205, 1500)
        ax.set_xlabel("CGI Length")
        ax.set_ylabel("No. of Islands")
        ax.legend()
    fig.tight_layout()
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_cgi_methylation_change(summary: pd.DataFrame, output_path: str | Path, title: str = "") -> None:
    plt = optional_pyplot()
    if plt is None or summary.empty:
        return
    class_labels = {
        "Class1_Pct": "Strong decrease (Delta <= -0.6)",
        "Class2_Pct": "Moderate decrease (-0.6 < Delta <= -0.2)",
        "Class3_Pct": "Minimal change (-0.2 < Delta <= 0.2)",
        "Class4_Pct": "Moderate increase (0.2 < Delta <= 0.6)",
        "Class5_Pct": "Strong increase (Delta > 0.6)",
    }
    class_cols = [c for c in summary.columns if c.startswith("Class") and c.endswith("_Pct")]
    if {"Class1_Pct", "Class2_Pct", "Class3_Pct", "Class4_Pct", "Class5_Pct"}.issubset(class_cols):
        class_cols = ["Class1_Pct", "Class2_Pct", "Class4_Pct", "Class5_Pct"]
    fig, ax = plt.subplots(figsize=(6, 4))
    for col in class_cols:
        ax.semilogx(summary["MeanLength"], summary[col], "-", linewidth=2, label=class_labels.get(col, col.replace("_Pct", "")))
    ax.set_xlabel("CGI Length (bp)")
    ax.set_ylabel("% of CGIs")
    ax.set_xlim(200, max(summary["MeanLength"]))
    ax.grid(True)
    if title:
        ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
