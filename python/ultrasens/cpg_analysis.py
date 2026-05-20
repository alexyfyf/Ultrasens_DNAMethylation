from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import ensure_dir, interp1, maybe_named_columns, optional_pyplot, save_data_struct


CPG_WGBS_DENSITY_COLUMNS = ["chr", "start", "end", "WGBS", "Density"]


def ed50_from_hypo_hyper(density: np.ndarray, hypo: np.ndarray, hyper: np.ndarray) -> float:
    fine_x = np.linspace(float(np.min(density)), float(np.max(density)), 300)
    diff = interp1(density, hypo, fine_x) - interp1(density, hyper, fine_x)
    return float(fine_x[int(np.argmin(np.abs(diff)))])


def get_ed50_and_slope(hypo: np.ndarray, hyper: np.ndarray, xax: np.ndarray) -> tuple[float, float, int]:
    """Match Plot_ExampleTemporal.m/GetED50."""
    fine_x = np.linspace(float(np.min(xax)), float(np.max(xax)), 300)
    hypo_fine = interp1(xax, hypo, fine_x)
    hyper_fine = interp1(xax, hyper, fine_x)
    diff = hypo_fine - hyper_fine
    keep = np.where((diff > -0.2) & (diff < 0.2))[0]
    if keep.size <= 3:
        return 0.0, 0.0, 1
    x_keep = fine_x[keep]
    d_keep = diff[keep]
    order = np.argsort(d_keep)
    ed50 = float(np.interp(0.0, d_keep[order], x_keep[order]))
    slope = float(np.polyfit(x_keep, d_keep, 1)[0])
    return ed50, slope, 0


def hill_decreasing(x: np.ndarray, k: float, n: float) -> np.ndarray:
    return (k**n) / (k**n + x**n)


def hill_increasing(x: np.ndarray, k: float, n: float) -> np.ndarray:
    return (x**n) / (k**n + x**n)


def fit_hill_one_param(ed50: float, x: np.ndarray, y: np.ndarray) -> float:
    grid = np.linspace(0.005, 50, 2000)
    errs = [np.sum((y - hill_decreasing(x, ed50, h)) ** 2) for h in grid]
    return float(grid[int(np.argmin(errs))])


def fit_hill_two_param(x: np.ndarray, y: np.ndarray, direction: int = 0) -> tuple[float, float, np.ndarray, float]:
    """Match RunTemporalCMEExample/FitHill2Param.m.

    Returns (n, K, model_on_original_x, sum_squared_error).
    """
    base_x = np.asarray(x, dtype=float)
    fine_x = np.linspace(float(np.min(base_x)), float(np.max(base_x)), 300)
    fine_y = interp1(base_x, np.asarray(y, dtype=float), fine_x)
    max_y = float(np.max(fine_y))
    if max_y <= 0:
        return 0.0, 0.0, np.zeros_like(base_x), 0.0
    data = fine_y / max_y

    def model(xx: np.ndarray, k: float, n: float) -> np.ndarray:
        return hill_increasing(xx, k, n) if direction else hill_decreasing(xx, k, n)

    k_grid = np.linspace(0.05, 20, 80)
    n_grid = np.linspace(0.05, 20, 80)
    best = (float("inf"), 1.1, 1.1)
    for k in k_grid:
        pred_cache = fine_x
        for n in n_grid:
            err = float(np.sum((data - model(pred_cache, k, n)) ** 2))
            if err < best[0]:
                best = (err, float(k), float(n))
    for span in (1.0, 0.2):
        _, best_k, best_n = best
        k_grid = np.linspace(max(1e-6, best_k - span), min(20.0, best_k + span), 60)
        n_grid = np.linspace(max(1e-6, best_n - span), min(20.0, best_n + span), 60)
        for k in k_grid:
            for n in n_grid:
                err = float(np.sum((data - model(fine_x, k, n)) ** 2))
                if err < best[0]:
                    best = (err, float(k), float(n))
    err, k, n = best
    return n, k, max_y * model(base_x, k, n), err


def fit_switch_values(model_struct: dict[str, np.ndarray]) -> np.ndarray:
    """Return [ED50, slope, KFD, nFD, KFM, nFM] like Plot_ExampleTemporal.m."""
    hypo = np.asarray(model_struct["Hypo"], dtype=float)
    hyper = np.asarray(model_struct["Hyper"], dtype=float)
    mean = np.asarray(model_struct["MeanMeth" if "MeanMeth" in model_struct else "Mean"], dtype=float)
    xax = np.asarray(model_struct["densityvals"], dtype=float)
    ed50, slope, _flag = get_ed50_and_slope(hypo, hyper, xax)
    diff = hyper - hypo
    cdiff = diff + 1
    n_fd, k_fd, _model, _err = fit_hill_two_param(xax, cdiff, direction=0)
    n_fm, k_fm, _model, _err = fit_hill_two_param(xax, mean, direction=0)
    return np.asarray([ed50, slope, k_fd, n_fd, k_fm, n_fm], dtype=float)


def fit_switch_values_temporal(model_struct: dict[str, np.ndarray]) -> np.ndarray:
    hyper = np.asarray(model_struct["Hyper"], dtype=float)
    rows = []
    for i in range(hyper.shape[0]):
        rows.append(
            fit_switch_values(
                {
                    "densityvals": model_struct["densityvals"],
                    "Hypo": np.asarray(model_struct["Hypo"])[i],
                    "Hyper": np.asarray(model_struct["Hyper"])[i],
                    "Mean": np.asarray(model_struct["Mean"])[i],
                }
            )
        )
    return np.vstack(rows)


def fit_hill_log_transform(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    y = np.clip(y, 1e-6, 1 - 1e-6)
    x = np.clip(x, 1e-12, None)
    midpoint = int(np.argmin(np.abs(y - 0.5)))
    start = max(0, midpoint - 3)
    end = min(len(x), midpoint + 4)
    lx = np.log10(x[start:end])
    ly = np.log10(y[start:end] / (1 - y[start:end]))
    if len(lx) < 2:
        return float(x[midpoint]), 0.0
    slope, intercept = np.polyfit(lx, ly, 1)
    n = -float(slope)
    k = float(10 ** (intercept / n)) if n != 0 else float(x[midpoint])
    return k, n


def summarize_cpg_density_file(
    input_path: str | Path,
    hypo_cutoff: float = 0.2,
    output_prefix: str | Path | None = None,
) -> dict[str, np.ndarray | float]:
    df = pd.read_csv(input_path, sep=None, engine="python")
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
    hypo_counts = counts[:, hypo_bins].sum(axis=1)
    hyper_counts = counts[:, hyper_bins].sum(axis=1)
    inter_counts = counts[:, inter_bins].sum(axis=1)
    totals = np.maximum(hypo_counts + hyper_counts + inter_counts, 1)
    hypo = hypo_counts / totals
    hyper = hyper_counts / totals
    inter = inter_counts / totals
    means = np.array([np.mean(y[x == d]) for d in density_vals])
    medians = np.array([np.median(y[x == d]) for d in density_vals])
    prc25 = np.array([np.percentile(y[x == d], 25) for d in density_vals])
    prc75 = np.array([np.percentile(y[x == d], 75) for d in density_vals])
    ed50 = ed50_from_hypo_hyper(density_vals, hypo, hyper)
    norm_mean = (means - np.min(means)) / max(np.max(means) - np.min(means), 1e-12)
    hill_n = fit_hill_one_param(ed50, density_vals, norm_mean)
    hill_k_log, hill_n_log = fit_hill_log_transform(density_vals, norm_mean)
    result = {
        "densityvals": density_vals,
        "Hypo": hypo,
        "Hyper": hyper,
        "Inter": inter,
        "MeanMeth": means,
        "MedMeth": medians,
        "Prc25": prc25,
        "Prc75": prc75,
        "ED50": ed50,
        "HillN": hill_n,
        "HillKLog": hill_k_log,
        "HillNLog": hill_n_log,
    }
    direct_ed50, direct_slope, direct_flag = get_ed50_and_slope(hypo, hyper, density_vals)
    result.update(
        {
            "DirectED50": direct_ed50,
            "DirectSlope": direct_slope,
            "DirectFlag": direct_flag,
        }
    )
    if output_prefix is not None:
        save_data_struct(
            output_prefix,
            densityvals=density_vals,
            Hypo=hypo,
            Hyper=hyper,
            Inter=inter,
            MeanMeth=means,
            MedMeth=medians,
            Prc25=prc25,
            Prc75=prc75,
        )
    return result


def plot_cpg_summary(summary: dict, output_path: str | Path) -> None:
    plt = optional_pyplot()
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    density = summary["densityvals"]
    ax.plot(density, summary["Hypo"], "-o", label="Hypo")
    ax.plot(density, summary["Hyper"], "-o", label="Hyper")
    ax.plot(density, summary["Inter"], "-o", label="Inter")
    ax.set_xlabel("Local CpG Density")
    ax.set_ylabel("Frac. of CpGs")
    ax.legend()
    fig.tight_layout()
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path)
    plt.close(fig)


def hill_fit_to_mean(summary: dict) -> np.ndarray:
    """Return the ReadPlotData_WT_example.m Hill fit on the raw mean-methylation scale."""
    density = np.asarray(summary["densityvals"], dtype=float)
    means = np.asarray(summary["MeanMeth"], dtype=float)
    min_mean = float(np.nanmin(means))
    max_mean = float(np.nanmax(means))
    span = max(max_mean - min_mean, 1e-12)
    norm_fit = hill_decreasing(density, float(summary["ED50"]), float(summary["HillN"]))
    return norm_fit * span + min_mean


def plot_figure5a_density_panels(
    summaries: list[dict],
    labels: list[str],
    output_path: str | Path,
    max_neighbors_per_100bp: float = 50.0,
) -> None:
    """Plot Fig. 5A-style mean/median/IQR methylation versus local CpG density."""
    plt = optional_pyplot()
    if plt is None:
        return
    n = len(summaries)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 3.2), squeeze=False)
    for ax, summary, label in zip(axes[0], summaries, labels):
        density = np.asarray(summary["densityvals"], dtype=float)
        mean_meth = np.asarray(summary["MeanMeth"], dtype=float)
        med_meth = np.asarray(summary["MedMeth"], dtype=float)
        prc25 = np.asarray(summary["Prc25"], dtype=float)
        prc75 = np.asarray(summary["Prc75"], dtype=float)
        hill = hill_fit_to_mean(summary)

        ax.fill_between(density, prc25, prc75, color="#bfc4ff", alpha=0.45, label="IQR")
        ax.plot(density, med_meth, color="#5963ff", linewidth=1.4, label="MedMeth")
        ax.plot(density, mean_meth, color="#e9531a", linewidth=3.0, label="MeanMeth")
        ax.plot(density, hill, color="#00d43a", linestyle=(0, (4, 4)), linewidth=2.7, label="Hill Model")

        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("Local CpG Density")
        ax.set_ylabel("Methyl. Frac")
        ax.set_xlim(0, max(0.5, float(np.nanmax(density))))
        ax.set_ylim(0, 1)
        ax.set_xticks([0, 0.5])
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(direction="out")
        ax.set_box_aspect(1)

        secax = ax.secondary_xaxis(
            "bottom",
            functions=(lambda x: x * max_neighbors_per_100bp, lambda x: x / max_neighbors_per_100bp),
        )
        secax.set_xlabel("# CpG Neighbors in 100 bp")
        secax.set_xticks([0, 10, 20])
        secax.spines["bottom"].set_position(("outward", 40))
        secax.tick_params(direction="out")

    handles, legend_labels = axes[0][-1].get_legend_handles_labels()
    unique = dict(zip(legend_labels, handles))
    axes[0][-1].legend(unique.values(), unique.keys(), loc="upper right", frameon=True)
    fig.tight_layout()
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


SWITCH_KO_COLORS = {
    "None": "#4d4d4d",
    "DNMT3AB": "#ff35e6",
    "DNMT3(single)": "#ff8a8a",
    "DNMT3 single": "#ff8a8a",
    "DNMT1": "#fff22e",
    "Combinatorial": "#2955ff",
    "TET1-3": "#47e6a0",
}


def switch_record_from_summary(
    summary: dict,
    label: str,
    cohort: str = "All",
    ko_type: str = "None",
) -> dict[str, str | float]:
    """Return the Figure 5B/C point values from one CpG-density summary."""
    return {
        "label": label,
        "cohort": cohort,
        "ko_type": ko_type,
        "HillK": float(summary["HillKLog"]),
        "HillN": float(summary["HillNLog"]),
        "DirectED50": float(summary.get("DirectED50", summary["ED50"])),
        "DirectSlope": float(summary.get("DirectSlope", 0.0)),
    }


def _symmetric_or_asymmetric_errors(row: pd.Series, x_base: str, y_base: str):
    def one_axis(base: str):
        low = row.get(f"{base}ErrLow")
        high = row.get(f"{base}ErrHigh")
        err = row.get(f"{base}Err")
        if pd.notna(low) and pd.notna(high):
            return [[float(low)], [float(high)]]
        if pd.notna(err):
            return float(err)
        return None

    return one_axis(x_base), one_axis(y_base)


def plot_switch_parameter_panels(
    records: pd.DataFrame,
    output_path: str | Path,
    cohort_order: list[str] | None = None,
    ko_colors: dict[str, str] | None = None,
) -> None:
    """Plot Figure 5B/C-style Hill K/n and direct ED50/slope switch summaries."""
    plt = optional_pyplot()
    if plt is None:
        return
    records = records.copy()
    if "cohort" not in records:
        records["cohort"] = "All"
    if "ko_type" not in records:
        records["ko_type"] = "None"
    colors = {**SWITCH_KO_COLORS, **(ko_colors or {})}
    cohorts = cohort_order or list(dict.fromkeys(records["cohort"].astype(str)))
    fig, axes = plt.subplots(2, len(cohorts), figsize=(4.0 * len(cohorts), 6.6), squeeze=False)

    specs = [
        ("HillK", "HillN", "K (threshold CpG dens.)", "n (steepness)", "Hill Fit Method"),
        ("DirectED50", "DirectSlope", "ED50 (threshold CpG dens.)", "slope (steepness)", "Direct Method"),
    ]
    for col, cohort in enumerate(cohorts):
        subset = records[records["cohort"].astype(str) == str(cohort)]
        for row_i, (x_col, y_col, x_label, y_label, method_label) in enumerate(specs):
            ax = axes[row_i][col]
            wt = subset[subset["ko_type"].astype(str) == "None"]
            if not wt.empty:
                ax.axvline(float(wt.iloc[0][x_col]), color="0.45", linestyle=(0, (7, 5)), linewidth=1.1)
                ax.axhline(float(wt.iloc[0][y_col]), color="0.45", linestyle=(0, (7, 5)), linewidth=1.1)
            for _, point in subset.iterrows():
                ko_type = str(point["ko_type"])
                color = colors.get(ko_type, "#888888")
                xerr, yerr = _symmetric_or_asymmetric_errors(point, x_col, y_col)
                ax.errorbar(
                    [float(point[x_col])],
                    [float(point[y_col])],
                    xerr=xerr,
                    yerr=yerr,
                    fmt="o",
                    markersize=6,
                    color=color,
                    ecolor=color,
                    markeredgecolor="0.25",
                    elinewidth=1.1,
                    capsize=3,
                    alpha=0.9,
                )
                ax.annotate(
                    str(point["label"]),
                    (float(point[x_col]), float(point[y_col])),
                    xytext=(4, 2),
                    textcoords="offset points",
                    fontsize=8,
                )
            ax.set_title(f"{cohort}, {method_label}", fontweight="bold")
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.grid(True, color="0.85", linewidth=0.8)
            ax.tick_params(direction="out")

    legend_items = []
    legend_labels = []
    for ko_type in dict.fromkeys(records["ko_type"].astype(str)):
        handle = plt.Line2D([0], [0], marker="o", linestyle="", color=colors.get(ko_type, "#888888"))
        legend_items.append(handle)
        legend_labels.append(ko_type)
    axes[0][-1].legend(legend_items, legend_labels, title="Enzyme KO", loc="best", frameon=True)
    fig.tight_layout()
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
