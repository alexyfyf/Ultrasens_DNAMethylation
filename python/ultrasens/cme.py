from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .cpg_density import compute_cpg_density
from .utils import ensure_dir, interp1, optional_pyplot


@dataclass
class CMEResult:
    pvec_msm: np.ndarray
    mbins: np.ndarray
    pvec: np.ndarray
    prob_ind: np.ndarray
    ind_cpg_p: np.ndarray
    matrix: np.ndarray
    types: np.ndarray
    meth: np.ndarray


def enumerate_types(n_cpg: int) -> np.ndarray:
    rows = []
    for nm in range(n_cpg + 1):
        for nh in range(n_cpg - nm + 1):
            nu = n_cpg - nm - nh
            rows.append((nu, nh, nm))
    return np.asarray(rows, dtype=int)


def stationary_distribution(generator: np.ndarray) -> np.ndarray:
    a = np.asarray(generator, dtype=float).copy()
    b = np.zeros(a.shape[0], dtype=float)
    a[-1, :] = 1.0
    b[-1] = 1.0
    try:
        p = np.linalg.solve(a, b)
    except np.linalg.LinAlgError:
        p = np.linalg.lstsq(a, b, rcond=None)[0]
    p = np.real(p)
    p[p < 0] = 0
    total = p.sum()
    return p / total if total > 0 else np.ones_like(p) / len(p)


def build_cme_matrix(n_cpg: int, parameters, cpg_positions, dist_length) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    k = np.asarray(parameters, dtype=float)
    dist_length = np.asarray(dist_length, dtype=float)
    positions = np.asarray(cpg_positions, dtype=float)
    dist_mat = np.abs(positions[:, None] - positions[None, :]) + np.eye(n_cpg) * 1e6
    types = enumerate_types(n_cpg)
    type_to_idx = {tuple(row): i for i, row in enumerate(types)}
    n_macro = len(types)
    collab = np.arange(4, 12)
    vals, inverse = np.unique(dist_length[collab], return_inverse=True)
    dist_assign = np.full(12, -1, dtype=int)
    dist_assign[collab] = inverse
    inv_lengths = 1 / vals
    fu = np.zeros((n_macro, len(vals)))
    fh = np.zeros_like(fu)
    fm = np.zeros_like(fu)
    meth = types[:, 1] * 0.5 + types[:, 2]
    for m, (nu, nh, nm) in enumerate(types):
        pu, ph, pm = np.asarray([nu, nh, nm], dtype=float) / n_cpg
        for i, inv_len in enumerate(inv_lengths):
            exp_sum = np.exp(-inv_len * dist_mat).sum()
            fu[m, i] = pu * exp_sum / n_cpg
            fh[m, i] = ph * exp_sum / n_cpg
            fm[m, i] = pm * exp_sum / n_cpg
    std = np.array(
        [
            [-k[0], k[1] + k[13], 0],
            [k[0], -(k[1] + k[13] + k[2]), k[3] + k[12]],
            [0, k[2], -(k[3] + k[12])],
        ],
        dtype=float,
    )
    prob_ind = stationary_distribution(std)
    msm = np.zeros((n_macro, n_macro), dtype=float)
    stoich = np.asarray([[-1, 1, 0], [1, -1, 0], [0, -1, 1], [0, 1, -1]], dtype=int)
    for m, (nu, nh, nm) in enumerate(types):
        fu9, fu11 = fu[m, dist_assign[8]], fu[m, dist_assign[10]]
        fh5, fh7 = fh[m, dist_assign[4]], fh[m, dist_assign[6]]
        fh10, fh12 = fh[m, dist_assign[9]], fh[m, dist_assign[11]]
        fm6, fm8 = fm[m, dist_assign[5]], fm[m, dist_assign[7]]
        collab_matrix = np.array(
            [
                [-k[6] * fh7 - k[7] * fm8, k[8] * fu9 + k[9] * fh10, 0],
                [
                    k[6] * fh7 + k[7] * fm8,
                    -k[8] * fu9 - k[9] * fh10 - k[4] * fh5 - k[5] * fm6,
                    k[11] * fh12 + k[10] * fu11,
                ],
                [0, k[4] * fh5 + k[5] * fm6, -k[11] * fh12 - k[10] * fu11],
            ],
            dtype=float,
        )
        rm = std + collab_matrix
        counts = [nu, nh, nh, nm]
        rates = [rm[1, 0], rm[0, 1], rm[2, 1], rm[1, 2]]
        for rx, change in enumerate(stoich):
            new_type = np.asarray([nu, nh, nm]) + change
            if np.all(new_type >= 0):
                idx = type_to_idx[tuple(new_type.tolist())]
                msm[idx, m] += counts[rx] * rates[rx]
    msm = msm - np.diag(msm.sum(axis=0))
    return msm, types, meth, prob_ind


def cme_model(n_cpg: int, parameters, cpg_positions, dist_length) -> CMEResult:
    msm, types, meth, prob_ind = build_cme_matrix(n_cpg, parameters, cpg_positions, dist_length)
    pvec = stationary_distribution(msm)
    mbins = np.arange(0, n_cpg + 0.5, 0.5)
    pvec_msm = np.asarray([pvec[meth == val].sum() for val in mbins])
    ind_cpg_p = np.array(
        [
            np.sum(types[:, 0] * pvec) / n_cpg,
            np.sum(types[:, 1] * pvec) / n_cpg,
            np.sum(types[:, 2] * pvec) / n_cpg,
        ]
    )
    return CMEResult(pvec_msm, mbins, pvec, prob_ind, ind_cpg_p, msm, types, meth)


def model_curves(parameters, dist_length, n_cpg: int = 27, ds=None, density_window: int = 50) -> dict[str, np.ndarray]:
    if ds is None:
        ds = [2, 4, 5, 6, 7, 8, 9, 11, 13, 17, 26, 110]
    ds = list(reversed(ds))
    hyper_val = 0.8
    hypo_val = 0.2
    density_vals = []
    hyper = []
    hypo = []
    mean = []
    for d in ds:
        positions = np.arange(1, n_cpg * d + 1, d)
        _, densities = compute_cpg_density(positions, density_window)
        result = cme_model(n_cpg, parameters, positions, dist_length)
        p = np.clip(result.pvec_msm, 0, None)
        p = p / p.sum()
        meth_ratio = result.mbins / n_cpg
        density_vals.append(float(np.mean(densities)))
        hyper.append(float(p[meth_ratio > hyper_val].sum()))
        hypo.append(float(p[meth_ratio < hypo_val].sum()))
        mean.append(float(np.sum(p * meth_ratio)))
    return {
        "densityvals": np.asarray(density_vals),
        "Hyper": np.asarray(hyper),
        "Hypo": np.asarray(hypo),
        "MeanMeth": np.asarray(mean),
    }


def propagate_cme(generator: np.ndarray, p0: np.ndarray, tspan) -> np.ndarray:
    vals, vecs = np.linalg.eig(generator)
    inv_vecs = np.linalg.pinv(vecs)
    coeff = inv_vecs @ p0
    rows = []
    for t in tspan:
        p = np.real(vecs @ (np.exp(vals * t) * coeff))
        p[p < 0] = 0
        rows.append(p / p.sum())
    return np.asarray(rows)


def temporal_model(parameters1, parameters2, dist_length1, dist_length2, tspan, n_cpg: int = 27, ds=None) -> dict[str, np.ndarray]:
    if ds is None:
        ds = [2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 17, 26, 110]
    ds = list(reversed(ds))
    hyper_val = 0.8
    hypo_val = 0.2
    density_vals = []
    mean_rows = []
    hyper_rows = []
    hypo_rows = []
    for d in ds:
        positions = np.arange(1, n_cpg * d + 1, d)
        _, densities = compute_cpg_density(positions, 50)
        start = cme_model(n_cpg, parameters1, positions, dist_length1)
        matrix2, _types, meth, _prob = build_cme_matrix(n_cpg, parameters2, positions, dist_length2)
        p_micro_t = propagate_cme(matrix2, start.pvec, tspan)
        steady2 = cme_model(n_cpg, parameters2, positions, dist_length2)
        density_vals.append(float(np.mean(densities)))
        means = []
        hypers = []
        hypos = []
        mbins = start.mbins
        meth_ratio = mbins / n_cpg
        for p_micro in list(p_micro_t) + [steady2.pvec]:
            p_macro = np.asarray([p_micro[meth == val].sum() for val in mbins])
            p_macro = np.clip(p_macro, 0, None)
            p_macro = p_macro / p_macro.sum()
            means.append(float(np.sum(p_macro * meth_ratio)))
            hypers.append(float(p_macro[meth_ratio > hyper_val].sum()))
            hypos.append(float(p_macro[meth_ratio < hypo_val].sum()))
        mean_rows.append(means)
        hyper_rows.append(hypers)
        hypo_rows.append(hypos)
    return {
        "densityvals": np.asarray(density_vals),
        "Mean": np.asarray(mean_rows).T,
        "Hyper": np.asarray(hyper_rows).T,
        "Hypo": np.asarray(hypo_rows).T,
    }


def standard_parameters(parameters) -> np.ndarray:
    """Return a CME parameter vector with collaborative reactions removed."""
    params = np.asarray(parameters, dtype=float).copy()
    params[4:12] = 0.0
    return params


def cme_fit_comparison(parameters, dist_length, data_struct: dict, n_cpg: int = 27, ds=None) -> dict[str, dict[str, np.ndarray]]:
    """Return data, fitted collaborative model, and standard-model curves."""
    if ds is None:
        ds = [2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 17, 26, 110]
    collab = model_curves(parameters, dist_length, n_cpg=n_cpg, ds=ds)
    standard = model_curves(standard_parameters(parameters), dist_length, n_cpg=n_cpg, ds=ds)
    data = {
        "densityvals": np.asarray(data_struct["densityvals"], dtype=float),
        "Hyper": np.asarray(data_struct["Hyper"], dtype=float),
        "Hypo": np.asarray(data_struct["Hypo"], dtype=float),
    }
    data["Inter"] = 1.0 - data["Hyper"] - data["Hypo"]
    for curves in (collab, standard):
        curves["Inter"] = 1.0 - curves["Hyper"] - curves["Hypo"]
    return {"data": data, "collaborative": collab, "standard": standard}


def cme_fit_ssd(model: dict[str, np.ndarray], data: dict[str, np.ndarray]) -> float:
    """SSD against data using the CME fitting fields and interpolation used by MATLAB."""
    total = 0.0
    for field in ("Hyper", "Hypo", "MeanMeth"):
        if field in model and field in data:
            data_y = interp1(np.asarray(data["densityvals"]), np.asarray(data[field]), np.asarray(model["densityvals"]))
            total += float(np.sum((np.asarray(model[field]) - data_y) ** 2))
    return total


def _arrow_width(value: float, max_value: float, base: float = 0.8, scale: float = 5.0) -> float:
    if max_value <= 0:
        return base
    return base + scale * max(float(value), 0.0) / max_value


def _draw_arrow(ax, start, end, color, width, rad: float = 0.0) -> None:
    from matplotlib.patches import FancyArrowPatch

    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=13,
        linewidth=width,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=4,
        shrinkB=4,
    )
    ax.add_patch(arrow)


def draw_cme_schematic(ax, parameters, collaborative: bool = True) -> None:
    """Draw u/h/m reaction schematic with arrow widths scaled by fitted rates."""
    k = np.asarray(parameters, dtype=float)
    x_u, x_h, x_m = 0.15, 0.50, 0.85
    y_mid = 0.48
    for x, label in ((x_u, "u"), (x_h, "h"), (x_m, "m")):
        ax.text(x, y_mid, label, ha="center", va="center", fontsize=14, fontweight="bold")

    standard_rates = {
        "u_to_h": k[0],
        "h_to_u": k[1] + k[13],
        "h_to_m": k[2],
        "m_to_h": k[3] + k[12],
    }
    collab_rates = {
        "u_to_h": k[6] + k[7],
        "h_to_u": k[8] + k[9],
        "h_to_m": k[4] + k[5],
        "m_to_h": k[10] + k[11],
    }
    all_rates = list(standard_rates.values()) + (list(collab_rates.values()) if collaborative else [])
    max_rate = max(max(all_rates), 1e-12)

    _draw_arrow(ax, (x_u + 0.05, y_mid + 0.08), (x_h - 0.05, y_mid + 0.08), "black", _arrow_width(standard_rates["u_to_h"], max_rate))
    _draw_arrow(ax, (x_h + 0.05, y_mid + 0.08), (x_m - 0.05, y_mid + 0.08), "black", _arrow_width(standard_rates["h_to_m"], max_rate))
    _draw_arrow(ax, (x_h - 0.05, y_mid - 0.08), (x_u + 0.05, y_mid - 0.08), "#1bb35e", _arrow_width(standard_rates["h_to_u"], max_rate))
    _draw_arrow(ax, (x_m - 0.05, y_mid - 0.08), (x_h + 0.05, y_mid - 0.08), "#1bb35e", _arrow_width(standard_rates["m_to_h"], max_rate))

    if collaborative:
        _draw_arrow(ax, (x_u + 0.03, y_mid + 0.22), (x_h - 0.06, y_mid + 0.17), "red", _arrow_width(collab_rates["u_to_h"], max_rate), rad=-0.35)
        _draw_arrow(ax, (x_h + 0.05, y_mid + 0.20), (x_m - 0.05, y_mid + 0.17), "red", _arrow_width(collab_rates["h_to_m"], max_rate), rad=-0.35)
        _draw_arrow(ax, (x_h - 0.05, y_mid - 0.20), (x_u + 0.05, y_mid - 0.17), "#1f58ff", _arrow_width(collab_rates["h_to_u"], max_rate), rad=0.35)
        _draw_arrow(ax, (x_m - 0.05, y_mid - 0.21), (x_h + 0.05, y_mid - 0.17), "#1f58ff", _arrow_width(collab_rates["m_to_h"], max_rate), rad=0.35)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def plot_cme_fit_comparison_panel(
    parameters,
    dist_length,
    data_struct: dict,
    output_path: str | Path,
    n_cpg: int = 27,
    label: str = "HUES8",
    ds=None,
) -> None:
    """Plot the Figure-style CME collaborative fit and no-collaboration standard model."""
    plt = optional_pyplot()
    if plt is None:
        return
    comparison = cme_fit_comparison(parameters, dist_length, data_struct, n_cpg=n_cpg, ds=ds)
    data = comparison["data"]
    fig = plt.figure(figsize=(6.7, 5.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 1.35], hspace=0.22, wspace=0.28)
    ax_scheme_collab = fig.add_subplot(gs[0, 0])
    ax_scheme_standard = fig.add_subplot(gs[0, 1])
    ax_collab = fig.add_subplot(gs[1, 0])
    ax_standard = fig.add_subplot(gs[1, 1])

    draw_cme_schematic(ax_scheme_collab, parameters, collaborative=True)
    draw_cme_schematic(ax_scheme_standard, parameters, collaborative=False)

    for ax, title, key in ((ax_collab, "Collaborative", "collaborative"), (ax_standard, "Standard", "standard")):
        model = comparison[key]
        ax.plot(model["densityvals"], model["Hyper"], color="red", linewidth=2.7, label="Model hyper")
        ax.plot(model["densityvals"], model["Hypo"], color="blue", linewidth=2.7, label="Model hypo")
        ax.plot(model["densityvals"], model["Inter"], color="black", linewidth=2.0, label="Model inter")
        ax.scatter(data["densityvals"], data["Hyper"], s=30, color="#ff7f7f", edgecolor="0.35", alpha=0.85, label=f"{label} hyper")
        ax.scatter(data["densityvals"], data["Hypo"], s=30, color="#8888ff", edgecolor="0.35", alpha=0.85, label=f"{label} hypo")
        ax.scatter(data["densityvals"], data["Inter"], s=24, color="0.65", edgecolor="0.35", alpha=0.85, label=f"{label} inter")
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Local CpG Density")
        ax.set_ylabel("Frac. of CpGs")
        ax.set_xlim(0, max(0.45, float(np.nanmax(data["densityvals"]))))
        ax.set_ylim(0, 1.02)
        ax.tick_params(direction="out")
        ax.grid(False)
    ax_collab.legend(loc="center right", fontsize=8, frameon=False)
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.10, top=0.94, hspace=0.25, wspace=0.28)
    ensure_dir(Path(output_path).parent)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
