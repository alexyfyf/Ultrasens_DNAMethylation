from __future__ import annotations

import numpy as np

from .cme import model_curves
from .utils import interp1


def update_parameters(p, parameters, dist_length):
    parameters = np.asarray(parameters, dtype=float).copy()
    dist_length = np.asarray(dist_length, dtype=float).copy()
    parameters[0] = p[0]
    parameters[1] = p[1]
    parameters[3] = p[2]
    parameters[5] = p[3]
    parameters[7] = p[4]
    parameters[8] = p[5]
    parameters[10] = p[6]
    dist_length[4:8] = p[7]
    dist_length[8:12] = p[8]
    return parameters, dist_length


def cme_fit_error(p, data_struct, parameters, dist_length, n_cpg: int = 27, ds=None) -> float:
    params, dists = update_parameters(p, parameters, dist_length)
    model = model_curves(params, dists, n_cpg=n_cpg, ds=ds)
    xq = model["densityvals"]
    fields = ["Hyper", "Hypo", "MeanMeth"]
    err = 0.0
    for field in fields:
        if field not in data_struct:
            continue
        data_y = interp1(data_struct["densityvals"], data_struct[field], xq)
        err += float(np.sum((model[field] - data_y) ** 2))
    return err


def fit_cme_parameters(parameters, dist_length, data_struct, n_cpg: int = 27, iterations: int = 60, seed: int = 1, ds=None):
    """Particleswarm-like bounded random search fallback; scipy users can replace with differential_evolution."""
    rng = np.random.default_rng(seed)
    p0 = np.asarray(
        [parameters[0], parameters[1], parameters[3], parameters[5], parameters[7], parameters[8], parameters[10], dist_length[4], dist_length[8]],
        dtype=float,
    )
    lo = np.ones_like(p0) * 1e-5
    hi = np.ones_like(p0) * 25
    hi[:3] = 1
    hi[7:9] = 100
    best = np.clip(p0, lo, hi)
    best_err = cme_fit_error(best, data_struct, parameters, dist_length, n_cpg=n_cpg, ds=ds)
    for i in range(iterations):
        scale = max(0.05, 1 - i / max(iterations, 1))
        if i < len(p0):
            cand = best.copy()
            cand[i] = rng.uniform(lo[i], hi[i])
        else:
            cand = best * np.exp(rng.normal(0, scale, size=best.shape))
        cand = np.clip(cand, lo, hi)
        err = cme_fit_error(cand, data_struct, parameters, dist_length, n_cpg=n_cpg, ds=ds)
        if err < best_err:
            best, best_err = cand, err
    return best, best_err, update_parameters(best, parameters, dist_length)


def cme_fit_bounds(parameters) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p0 = np.asarray(
        [parameters[0], parameters[1], parameters[3], parameters[5], parameters[7], parameters[8], parameters[10], 33.6308, 43.3724],
        dtype=float,
    )
    lo = np.ones_like(p0) * 1e-5
    hi = np.ones_like(p0) * 25
    hi[:3] = 1
    hi[7:9] = 100
    return p0, lo, hi


def fit_cme_parameters_scipy(
    parameters,
    dist_length,
    data_struct,
    n_cpg: int = 27,
    maxiter: int = 35,
    popsize: int = 8,
    seed: int = 1,
    ds=None,
    polish: bool = True,
):
    """Differential-evolution fit with L-BFGS-B polish, analogous to particleswarm + fmincon."""
    try:
        from scipy.optimize import differential_evolution
    except Exception as exc:
        raise RuntimeError("SciPy is required for --optimizer scipy") from exc

    p0, lo, hi = cme_fit_bounds(parameters)
    p0[7] = dist_length[4]
    p0[8] = dist_length[8]
    bounds = list(zip(lo, hi))

    def objective(p):
        return cme_fit_error(np.asarray(p), data_struct, parameters, dist_length, n_cpg=n_cpg, ds=ds)

    result = differential_evolution(
        objective,
        bounds,
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        polish=polish,
        updating="immediate",
        workers=1,
        x0=np.clip(p0, lo, hi),
    )
    best = np.asarray(result.x, dtype=float)
    best_err = float(result.fun)
    return best, best_err, update_parameters(best, parameters, dist_length), result


def fit_cme_parameters_lbfgsb(
    parameters,
    dist_length,
    data_struct,
    n_cpg: int = 27,
    maxiter: int = 300,
    maxfun: int = 2000,
    ds=None,
):
    """Bounded local fit from the MATLAB default starting point."""
    try:
        from scipy.optimize import minimize
    except Exception as exc:
        raise RuntimeError("SciPy is required for --optimizer lbfgsb") from exc

    p0, lo, hi = cme_fit_bounds(parameters)
    p0[7] = dist_length[4]
    p0[8] = dist_length[8]
    bounds = list(zip(lo, hi))

    def objective(p):
        return cme_fit_error(np.asarray(p), data_struct, parameters, dist_length, n_cpg=n_cpg, ds=ds)

    result = minimize(
        objective,
        np.clip(p0, lo, hi),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": maxiter, "maxfun": maxfun, "ftol": 1e-10},
    )
    best = np.asarray(result.x, dtype=float)
    best_err = float(result.fun)
    return best, best_err, update_parameters(best, parameters, dist_length), result


def update_example_parameters(p, parameters):
    parameters = np.asarray(parameters, dtype=float).copy()
    parameters[3] = p[0]
    parameters[5] = p[1]
    parameters[7] = p[2]
    parameters[8] = p[3]
    parameters[10] = p[4]
    return parameters


def example_fit_error(p, data_struct, parameters, dist_length, n_cpg: int = 27, ds=None) -> float:
    params = update_example_parameters(p, parameters)
    model = model_curves(params, dist_length, n_cpg=n_cpg, ds=ds)
    xq = model["densityvals"]
    err = 0.0
    for field in ("Hyper", "Hypo"):
        data_y = interp1(data_struct["densityvals"], data_struct[field], xq)
        err += float(np.sum((model[field] - data_y) ** 2))
    return err


def fit_example_cme_parameters(parameters, dist_length, data_struct, n_cpg: int = 27, iterations: int = 60, seed: int = 1, ds=None):
    """ExampleCMEFit-style bounded fit over k4, k6, k8, k9, and k11."""
    rng = np.random.default_rng(seed)
    p0 = np.asarray([parameters[3], parameters[5], parameters[7], parameters[8], parameters[10]], dtype=float)
    lo = np.ones_like(p0) * 1e-2
    hi = np.asarray([3, 25, 25, 25, 25], dtype=float)
    best = np.clip(p0, lo, hi)
    best_err = example_fit_error(best, data_struct, parameters, dist_length, n_cpg=n_cpg, ds=ds)
    for i in range(iterations):
        scale = max(0.05, 1 - i / max(iterations, 1))
        if i < len(p0):
            cand = best.copy()
            cand[i] = rng.uniform(lo[i], hi[i])
        else:
            cand = best * np.exp(rng.normal(0, scale, size=best.shape))
        cand = np.clip(cand, lo, hi)
        err = example_fit_error(cand, data_struct, parameters, dist_length, n_cpg=n_cpg, ds=ds)
        if err < best_err:
            best, best_err = cand, err
    return best, best_err, update_example_parameters(best, parameters)
