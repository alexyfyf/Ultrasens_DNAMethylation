#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.fit import fit_cme_parameters, fit_cme_parameters_lbfgsb, fit_cme_parameters_scipy
from ultrasens.utils import load_data_struct, save_data_struct, write_json


def parse_vector(value: str) -> list[float]:
    return [float(x) for x in value.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser(description="Python rewrite of CMEFitting/LoopFit.m + Fit_CME_Methylation_PS.m")
    parser.add_argument("--data-struct", required=True, help=".npz/.json exported by read_plot_data_wt.py, or .mat with scipy.")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--n-cpg", type=int, default=27)
    parser.add_argument("--iterations", type=int, default=60)
    parser.add_argument("--optimizer", choices=["random", "scipy", "lbfgsb"], default="random")
    parser.add_argument("--maxiter", type=int, default=35, help="SciPy differential-evolution generations.")
    parser.add_argument("--popsize", type=int, default=8, help="SciPy differential-evolution population multiplier.")
    parser.add_argument("--maxfun", type=int, default=2000, help="Maximum function evaluations for L-BFGS-B.")
    parser.add_argument("--no-polish", action="store_true", help="Disable SciPy L-BFGS-B polishing.")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--parameters",
        default="0.0999,0.1000,1,0.1431,0,3.5747,0,0.7281,0.0639,0,24.7399,0,0.0416666667,0.0208333333",
    )
    parser.add_argument(
        "--dist-length",
        default="0,0,0,0,33.6308,33.6308,33.6308,33.6308,43.3724,43.3724,43.3724,43.3724",
    )
    args = parser.parse_args()

    data = load_data_struct(args.data_struct)
    parameters = np.asarray(parse_vector(args.parameters), dtype=float)
    dist_length = np.asarray(parse_vector(args.dist_length), dtype=float)
    if args.optimizer == "scipy":
        best, ssd, (fit_params, fit_dist), opt_result = fit_cme_parameters_scipy(
            parameters,
            dist_length,
            data,
            n_cpg=args.n_cpg,
            maxiter=args.maxiter,
            popsize=args.popsize,
            seed=args.seed,
            polish=not args.no_polish,
        )
        optimizer_summary = {
            "optimizer": "scipy.differential_evolution",
            "maxiter": args.maxiter,
            "popsize": args.popsize,
            "polish": not args.no_polish,
            "message": str(opt_result.message),
            "nit": int(opt_result.nit),
            "nfev": int(opt_result.nfev),
            "success": bool(opt_result.success),
        }
    elif args.optimizer == "lbfgsb":
        best, ssd, (fit_params, fit_dist), opt_result = fit_cme_parameters_lbfgsb(
            parameters,
            dist_length,
            data,
            n_cpg=args.n_cpg,
            maxiter=args.maxiter,
            maxfun=args.maxfun,
        )
        optimizer_summary = {
            "optimizer": "scipy.minimize.L-BFGS-B",
            "maxiter": args.maxiter,
            "maxfun": args.maxfun,
            "message": str(opt_result.message),
            "nit": int(opt_result.nit),
            "nfev": int(opt_result.nfev),
            "success": bool(opt_result.success),
        }
    else:
        best, ssd, (fit_params, fit_dist) = fit_cme_parameters(
            parameters, dist_length, data, n_cpg=args.n_cpg, iterations=args.iterations, seed=args.seed
        )
        optimizer_summary = {"optimizer": "bounded_random_search", "iterations": args.iterations}
    save_data_struct(args.output_prefix, Parameters=fit_params, DistLength=fit_dist, BestFitVector=best, SSD=np.asarray([ssd]))
    write_json(
        Path(args.output_prefix).with_suffix(".summary.json"),
        {
            "SSD": ssd,
            "BestFitVector": best.tolist(),
            "Parameters": fit_params.tolist(),
            "DistLength": fit_dist.tolist(),
            **optimizer_summary,
        },
    )
    print(f"wrote CME fit result to {args.output_prefix}.npz/.json")


if __name__ == "__main__":
    main()
