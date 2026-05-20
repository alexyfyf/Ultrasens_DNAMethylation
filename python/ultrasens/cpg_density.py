from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import ensure_dir


def compute_cpg_density(positions, window: int) -> tuple[np.ndarray, np.ndarray]:
    """Match CpGDensities_Function.m: neighbors within +/- window, normalized by 2*ceil(W/2)."""
    pos = np.asarray(positions, dtype=int)
    order = np.argsort(pos)
    pos = pos[order]
    max_nn = 2 * int(np.ceil(window / 2))
    starts = np.searchsorted(pos, pos - window, side="left")
    ends = np.searchsorted(pos, pos + window, side="right")
    counts = ends - starts - 1
    return counts.astype(int), counts / max_nn


def compute_density_table(
    positions: np.ndarray,
    window: int,
    chrom: str,
    name: str = "CpG",
    trim_edges_like_matlab: bool = False,
) -> pd.DataFrame:
    pos = np.asarray(positions, dtype=int)
    counts, densities = compute_cpg_density(pos, window)
    if trim_edges_like_matlab:
        wcpg = int(np.ceil(window / 2))
        keep = np.arange(wcpg, max(wcpg, len(pos) - wcpg))
        pos = pos[keep]
        counts = counts[keep]
        densities = densities[keep]
    return pd.DataFrame(
        {
            "chr": chrom,
            "CpGstart": pos,
            "CpGend": pos + 1,
            "Name": name,
            "Density": densities,
            "NumNeighbors": counts,
        }
    )


def load_cpg_positions(csv_path: str | Path, start_column: str = "start") -> np.ndarray:
    table = pd.read_csv(csv_path)
    if start_column not in table.columns:
        raise KeyError(f"{csv_path} must contain a '{start_column}' column.")
    return table[start_column].to_numpy(dtype=int)


def write_density_bed(table: pd.DataFrame, output_path: str | Path) -> None:
    output_path = Path(output_path)
    ensure_dir(output_path.parent)
    table[["chr", "CpGstart", "CpGend", "Name", "Density"]].to_csv(
        output_path, sep="\t", header=False, index=False
    )

