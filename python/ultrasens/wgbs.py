from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .intersect import WGBS_COLUMNS
from .utils import ensure_dir, maybe_named_columns


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
    deduplicate_cpg_strands: bool = False,
    deduplicate_method: str = "coverage-weighted",
    chunksize: int | None = None,
) -> pd.DataFrame:
    """Parse Bismark .cov/.cov.gz into the pipeline's 4-column WGBS BED format.

    Bismark coverage files use 1-based genomic coordinates:
    chr, start, end, methylation percentage, count methylated, count unmethylated.

    This function returns/writes chr, start, end, WGBS where start/end are converted
    to 0-based half-open BED coordinates and WGBS is a methylation fraction in [0, 1].
    """
    if methylation_source not in {"counts", "percentage"}:
        raise ValueError("methylation_source must be 'counts' or 'percentage'")
    if deduplicate_method not in {"first", "second", "mean", "coverage-weighted"}:
        raise ValueError("deduplicate_method must be 'first', 'second', 'mean', or 'coverage-weighted'")
    if methylation_source == "percentage" and deduplicate_method == "coverage-weighted":
        raise ValueError("coverage-weighted deduplication requires methylation_source='counts'")
    if min_coverage < 0:
        raise ValueError("min_coverage must be non-negative")
    if chunksize is not None:
        if chunksize <= 0:
            raise ValueError("chunksize must be positive")
        if sort:
            raise ValueError("chunked Bismark parsing requires sort=False; sort the output with bedtools")
        return parse_bismark_coverage_chunked(
            input_path,
            output_path,
            methylation_source=methylation_source,
            min_coverage=min_coverage,
            deduplicate_cpg_strands=deduplicate_cpg_strands,
            deduplicate_method=deduplicate_method,
            chunksize=chunksize,
        )

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
        if deduplicate_cpg_strands and deduplicate_method == "coverage-weighted":
            df = deduplicate_bismark_cpg_strands_counts(df, sort=sort)
            total = df["count_methylated"] + df["count_unmethylated"]
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
        if deduplicate_cpg_strands and deduplicate_method != "coverage-weighted":
            out = deduplicate_adjacent_cpg_strands(out, method=deduplicate_method, sort=False)

    if output_path is not None:
        output_path = Path(output_path)
        ensure_dir(output_path.parent)
        out.to_csv(output_path, sep="\t", header=False, index=False)
    return out


def parse_bismark_coverage_chunked(
    input_path: str | Path,
    output_path: str | Path | None,
    methylation_source: str = "counts",
    min_coverage: int = 1,
    deduplicate_cpg_strands: bool = False,
    deduplicate_method: str = "coverage-weighted",
    chunksize: int = 1_000_000,
) -> pd.DataFrame:
    """Stream large Bismark coverage files without materializing all rows.

    Chunked parsing preserves input order. If sorted output is needed for bedtools,
    write to disk and sort the result after this function completes.
    """
    if output_path is None:
        raise ValueError("chunked Bismark parsing requires an output_path")
    output_path = Path(output_path)
    ensure_dir(output_path.parent)
    total_rows = 0
    prev: tuple[str, int, int, float, float, float] | None = None
    opener = gzip.open if str(input_path).endswith(".gz") else open

    def parse_line(line: str) -> tuple[str, int, int, float, float, float] | None:
        if not line.strip() or line.startswith("#"):
            return None
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            raise ValueError("Bismark coverage rows must have at least 6 tab-delimited columns")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        if start < 1 or end < start:
            raise ValueError("Bismark coordinates must be 1-based with end >= start")
        pct = float(fields[3])
        meth = float(fields[4])
        unmeth = float(fields[5])
        return chrom, start, end, pct, meth, unmeth

    def fraction(row: tuple[str, int, int, float, float, float]) -> float:
        _chrom, _start, _end, pct, meth, unmeth = row
        total = meth + unmeth
        return meth / total if methylation_source == "counts" and total > 0 else pct / 100.0

    def write_row(handle, row: tuple[str, int, int, float, float, float], value: float | None = None) -> int:
        chrom, start, end, _pct, meth, unmeth = row
        if meth + unmeth < min_coverage:
            return 0
        wgbs = fraction(row) if value is None else value
        if not np.isfinite(wgbs):
            return 0
        if wgbs < 0 or wgbs > 1:
            raise ValueError("Parsed methylation fractions must be within [0, 1]")
        handle.write(f"{chrom}\t{start - 1}\t{end}\t{wgbs:.12g}\n")
        return 1

    def write_pair(
        handle,
        first: tuple[str, int, int, float, float, float],
        second: tuple[str, int, int, float, float, float],
    ) -> int:
        if deduplicate_method == "first":
            return write_row(handle, first)
        if deduplicate_method == "second":
            return write_row(handle, second)
        if deduplicate_method == "mean":
            return write_row(handle, second, (fraction(first) + fraction(second)) / 2.0)
        meth = first[4] + second[4]
        unmeth = first[5] + second[5]
        total = meth + unmeth
        return write_row(handle, (second[0], second[1], second[2], 100.0 * meth / total if total else np.nan, meth, unmeth))

    with opener(input_path, "rt") as source, output_path.open("w") as handle:
        for line_no, line in enumerate(source, start=1):
            row = parse_line(line)
            if row is None or row[4] + row[5] < min_coverage:
                continue
            if line_no % chunksize == 0:
                handle.flush()
                print(
                    f"[parse_bismark_cov] processed {line_no:,} input rows; wrote {total_rows:,} WGBS rows",
                    file=sys.stderr,
                    flush=True,
                )
            if not deduplicate_cpg_strands:
                total_rows += write_row(handle, row)
                continue
            if prev is None:
                prev = row
                continue
            adjacent = prev[0] == row[0] and prev[2] + 1 == row[1]
            if adjacent:
                total_rows += write_pair(handle, prev, row)
                prev = None
            else:
                total_rows += write_row(handle, prev)
                prev = row
        if prev is not None:
            total_rows += write_row(handle, prev)

    return pd.DataFrame({"rows_written": [total_rows]})


def deduplicate_bismark_cpg_strands_counts(df: pd.DataFrame, sort: bool = True) -> pd.DataFrame:
    """Coverage-weighted collapse of adjacent Bismark strand rows before fraction calculation."""
    work = df[BISMARK_COV_COLUMNS].copy()
    work["chr"] = work["chr"].astype(str)
    if sort:
        work = work.sort_values(["chr", "start_1based", "end_1based"], kind="mergesort").reset_index(drop=True)

    rows = []
    i = 0
    while i < len(work):
        current = work.iloc[i].copy()
        if i + 1 < len(work):
            nxt = work.iloc[i + 1]
            adjacent = current["chr"] == nxt["chr"] and int(current["end_1based"]) + 1 == int(nxt["start_1based"])
            if adjacent:
                meth = float(current["count_methylated"]) + float(nxt["count_methylated"])
                unmeth = float(current["count_unmethylated"]) + float(nxt["count_unmethylated"])
                total = meth + unmeth
                keep = nxt.copy()
                keep["count_methylated"] = meth
                keep["count_unmethylated"] = unmeth
                keep["methylation_percentage"] = 100.0 * meth / total if total > 0 else np.nan
                rows.append(keep)
                i += 2
                continue
        rows.append(current)
        i += 1
    return pd.DataFrame(rows, columns=BISMARK_COV_COLUMNS).reset_index(drop=True)


def deduplicate_adjacent_cpg_strands(
    wgbs: pd.DataFrame | str | Path,
    output_path: str | Path | None = None,
    method: str = "second",
    sort: bool = True,
    max_methylation_diff: float | None = None,
) -> pd.DataFrame:
    """Collapse adjacent opposite-strand rows that represent the same CpG dyad.

    Some WGBS exports report both strand observations for a CpG dinucleotide as
    adjacent 1-bp BED rows. For example, a CpG at positions 10469-10470 may be
    exported as 10468-10469 and 10469-10470 after BED conversion. Downstream
    scripts expect one row per CpG dyad, so this helper collapses such adjacent
    pairs.

    The output coordinate follows the kept row. The default `method="second"`
    matches the original repository's `awk 'NR % 2 == 0'` convention.
    """
    if method not in {"first", "second", "mean"}:
        raise ValueError("method must be 'first', 'second', or 'mean'")
    if isinstance(wgbs, (str, Path)):
        df = pd.read_csv(wgbs, sep="\t", header=None)
        df = maybe_named_columns(df, WGBS_COLUMNS)
    else:
        df = maybe_named_columns(wgbs.copy(), WGBS_COLUMNS)
    if df.empty:
        out = pd.DataFrame(columns=WGBS_COLUMNS)
    else:
        df = df[WGBS_COLUMNS].copy()
        df["chr"] = df["chr"].astype(str)
        df[["start", "end"]] = df[["start", "end"]].astype(np.int64)
        df["WGBS"] = pd.to_numeric(df["WGBS"], errors="raise")
        if sort:
            df = df.sort_values(["chr", "start", "end"], kind="mergesort").reset_index(drop=True)

        rows = []
        i = 0
        while i < len(df):
            current = df.iloc[i]
            if i + 1 < len(df):
                nxt = df.iloc[i + 1]
                adjacent = current["chr"] == nxt["chr"] and int(current["end"]) == int(nxt["start"])
                methyl_ok = True
                if max_methylation_diff is not None:
                    methyl_ok = abs(float(current["WGBS"]) - float(nxt["WGBS"])) <= max_methylation_diff
                if adjacent and methyl_ok:
                    if method == "first":
                        keep = current.copy()
                    else:
                        keep = nxt.copy()
                        if method == "mean":
                            keep["WGBS"] = (float(current["WGBS"]) + float(nxt["WGBS"])) / 2.0
                    rows.append(keep)
                    i += 2
                    continue
            rows.append(current.copy())
            i += 1
        out = pd.DataFrame(rows, columns=WGBS_COLUMNS).reset_index(drop=True)

    if output_path is not None:
        output_path = Path(output_path)
        ensure_dir(output_path.parent)
        out.to_csv(output_path, sep="\t", header=False, index=False)
    return out
