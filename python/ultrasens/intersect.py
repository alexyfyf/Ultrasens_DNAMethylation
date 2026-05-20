from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import ensure_dir, maybe_named_columns


WGBS_COLUMNS = ["chr", "start", "end", "WGBS"]
DENSITY_COLUMNS = ["chr", "start", "end", "name", "density"]


def intersect_wgbs_with_density(wgbs_path: str | Path, density_bed_path: str | Path, output_path: str | Path) -> pd.DataFrame:
    """Pure-Python replacement for WGBS_CpGIntersect_AllData.command."""
    wgbs = pd.read_csv(wgbs_path, sep="\t", header=None)
    den = pd.read_csv(density_bed_path, sep="\t", header=None)
    wgbs = maybe_named_columns(wgbs, WGBS_COLUMNS)
    den = maybe_named_columns(den, DENSITY_COLUMNS)
    wgbs = wgbs.sort_values(["chr", "start", "end"], kind="mergesort").reset_index(drop=True)
    den = den.sort_values(["chr", "start", "end"], kind="mergesort").reset_index(drop=True)
    out_parts = []
    for chrom, wgbs_chrom in wgbs.groupby("chr", sort=False):
        den_chrom = den[den["chr"] == chrom]
        if den_chrom.empty:
            continue
        den_starts = den_chrom["start"].to_numpy()
        den_ends = den_chrom["end"].to_numpy()
        den_density = den_chrom["density"].to_numpy()
        starts = wgbs_chrom["start"].to_numpy()
        ends = wgbs_chrom["end"].to_numpy()
        left = np.searchsorted(den_starts, starts, side="left")
        right = np.searchsorted(den_starts, ends, side="left")
        rows = []
        for row_idx, den_left, den_right in zip(range(len(wgbs_chrom)), left, right):
            if den_left == den_right:
                continue
            row = wgbs_chrom.iloc[row_idx]
            overlaps = den_ends[den_left:den_right] > row["start"]
            for density in den_density[den_left:den_right][overlaps]:
                rows.append((row["chr"], row["start"], row["end"], row["WGBS"], density))
        if rows:
            out_parts.append(pd.DataFrame(rows, columns=["chr", "start", "end", "WGBS", "density"]))
    out = pd.concat(out_parts, ignore_index=True) if out_parts else pd.DataFrame(columns=["chr", "start", "end", "WGBS", "density"])
    ensure_dir(Path(output_path).parent)
    out.to_csv(output_path, sep="\t", header=False, index=False, float_format="%.6g")
    return out
