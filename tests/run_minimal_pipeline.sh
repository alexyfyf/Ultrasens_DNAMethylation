#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT_DIR/tests/output/minimal_pipeline"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

"$PYTHON_BIN" "$ROOT_DIR/scripts/run_pipeline.py" \
  --cpg-csv "$ROOT_DIR/tests/minimal_data/chr1_cpgs.csv" \
  --wgbs-bed "$ROOT_DIR/tests/minimal_data/wgbs_proc.bed" \
  --cgi-intersection "$ROOT_DIR/tests/minimal_data/cgi_intersection.csv" \
  --output-dir "$OUT_DIR" \
  --density-window 10 \
  --n-cpg 4 \
  --fit-iterations 3

test -s "$OUT_DIR/CpGDensities_W10.bed"
test -s "$OUT_DIR/WGBS_CpGsOnly_Chr1.bed"
test -s "$OUT_DIR/IslandLvl_agg.csv"
test -s "$OUT_DIR/Save_WGBS_CpGsOnly_Chr1.npz"
test -s "$OUT_DIR/Fit_WGBS_CpGsOnly_Chr1.npz"
test -s "$OUT_DIR/Temporal_WGBS_CpGsOnly_Chr1.npz"

"$PYTHON_BIN" "$ROOT_DIR/scripts/deduplicate_cpg_strands.py" \
  --input "$ROOT_DIR/tests/minimal_data/wgbs_all_strands.bed" \
  --output "$OUT_DIR/wgbs_dedup_second.bed" \
  --method second
grep -q $'chr1\t11\t12\t0.2' "$OUT_DIR/wgbs_dedup_second.bed"
grep -q $'chr1\t21\t22\t0.9' "$OUT_DIR/wgbs_dedup_second.bed"

"$PYTHON_BIN" "$ROOT_DIR/scripts/parse_bismark_cov.py" \
  --cov "$ROOT_DIR/tests/minimal_data/bismark_example.cov" \
  --output "$OUT_DIR/bismark_proc.bed" \
  --deduplicate-cpg-strands \
  --deduplicate-method mean
test -s "$OUT_DIR/bismark_proc.bed"

echo "minimal pipeline ok"
