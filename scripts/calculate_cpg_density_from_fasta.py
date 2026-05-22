#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ultrasens.cpg_density import compute_density_table


def open_text(path: str | Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def write_chrom_density(
    handle,
    chrom: str,
    sequence_parts: list[str],
    window: int,
    start_offset: int,
) -> int:
    started = time.time()
    sequence = "".join(sequence_parts).upper()
    cpg_positions = []
    pos = sequence.find("CG")
    while pos != -1:
        cpg_positions.append(pos + start_offset)
        pos = sequence.find("CG", pos + 1)
    positions = np.asarray(cpg_positions, dtype=np.int64)
    if positions.size == 0:
        print(f"[cpg_density_from_fasta] {chrom}: no CpGs", file=sys.stderr, flush=True)
        return 0
    table = compute_density_table(positions, window, chrom)
    table[["chr", "CpGstart", "CpGend", "Name", "Density"]].to_csv(
        handle, sep="\t", header=False, index=False, mode="a"
    )
    elapsed = time.time() - started
    print(
        f"[cpg_density_from_fasta] {chrom}: wrote {len(table):,} density rows in {elapsed:.1f}s",
        file=sys.stderr,
        flush=True,
    )
    return int(len(table))


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate CpG-density BED directly from FASTA/FASTA.GZ.")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument(
        "--chroms",
        nargs="+",
        default=None,
        help="Optional chromosome/contig names to include. Defaults to every FASTA record.",
    )
    parser.add_argument(
        "--start-offset",
        type=int,
        default=0,
        help="Offset added to each zero-based CpG start. Use 1 to align with Bismark second-strand dyad rows.",
    )
    args = parser.parse_args()

    include = set(args.chroms) if args.chroms else None
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    current_chrom: str | None = None
    sequence_parts: list[str] = []
    with open_text(args.fasta) as source, out_path.open("w") as out:
        for line in source:
            if line.startswith(">"):
                if current_chrom is not None and (include is None or current_chrom in include):
                    total += write_chrom_density(out, current_chrom, sequence_parts, args.window, args.start_offset)
                current_chrom = line[1:].strip().split()[0]
                if include is None or current_chrom in include:
                    print(f"[cpg_density_from_fasta] reading {current_chrom}", file=sys.stderr, flush=True)
                sequence_parts = []
                continue
            if current_chrom is not None and (include is None or current_chrom in include):
                sequence_parts.append(line.strip())
        if current_chrom is not None and (include is None or current_chrom in include):
            total += write_chrom_density(out, current_chrom, sequence_parts, args.window, args.start_offset)

    print(f"wrote {total} CpG density rows to {args.output}")


if __name__ == "__main__":
    main()
