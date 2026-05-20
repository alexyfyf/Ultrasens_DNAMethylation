from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .intersect import WGBS_COLUMNS
from .utils import ensure_dir


BISMARK_COV_COLUMNS = [
    "chr",
    "start_1based",
    "end_1based",
    "methylation_percentage",
    "count_methylated",
    "count_unmethylated",
]


def parse_bismark_coverage(
    input_path: str | Path,
    output_path: str | Path | None = None,
    methylation_source: str = "counts",
    min_coverage: int = 1,
    sort: bool = True,
) -> pd.DataFrame:
    """Parse Bismark .cov/.cov.gz into the pipeline's 4-column WGBS BED format.

    Bismark coverage files use 1-based genomic coordinates:
    chr, start, end, methylation percentage, count methylated, count unmethylated.

    This function returns/writes chr, start, end, WGBS where start/end are converted
    to 0-based half-open BED coordinates and WGBS is a methylation fraction in [0, 1].
    """
    if methylation_source not in {"counts", "percentage"}:
        raise ValueError("methylation_source must be 'counts' or 'percentage'")
    if min_coverage < 0:
        raise ValueError("min_coverage must be non-negative")

    df = pd.read_csv(
        input_path,
        sep="\t",
        header=None,
        names=BISMARK_COV_COLUMNS,
        compression="infer",
        comment="#",
    )
    if df.empty:
        out = pd.DataFrame(columns=WGBS_COLUMNS)
    else:
        numeric_cols = BISMARK_COV_COLUMNS[1:]
        df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="raise")
        if (df["start_1based"] < 1).any() or (df["end_1based"] < df["start_1based"]).any():
            raise ValueError("Bismark coordinates must be 1-based with end >= start")

        total = df["count_methylated"] + df["count_unmethylated"]
        keep = total >= min_coverage
        df = df.loc[keep].copy()
        total = total.loc[keep]
        if methylation_source == "counts":
            wgbs = np.divide(
                df["count_methylated"].to_numpy(float),
                total.to_numpy(float),
                out=np.full(len(df), np.nan, dtype=float),
                where=total.to_numpy(float) > 0,
            )
        else:
            wgbs = df["methylation_percentage"].to_numpy(float) / 100.0

        out = pd.DataFrame(
            {
                "chr": df["chr"].astype(str),
                "start": df["start_1based"].astype(np.int64) - 1,
                "end": df["end_1based"].astype(np.int64),
                "WGBS": wgbs,
            }
        )
        out = out.dropna(subset=["WGBS"])
        if ((out["WGBS"] < 0) | (out["WGBS"] > 1)).any():
            raise ValueError("Parsed methylation fractions must be within [0, 1]")
        if sort and not out.empty:
            out = out.sort_values(["chr", "start", "end"], kind="mergesort").reset_index(drop=True)

    if output_path is not None:
        output_path = Path(output_path)
        ensure_dir(output_path.parent)
        out.to_csv(output_path, sep="\t", header=False, index=False)
    return out
