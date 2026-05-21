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
  --deduplicate-cpg-strands
test -s "$OUT_DIR/bismark_proc.bed"

printf 'chr1\t10\t10\t75.0\t3\t1\nchr1\t11\t11\t0.0\t0\t6\n' > "$OUT_DIR/bismark_paired.cov"
"$PYTHON_BIN" "$ROOT_DIR/scripts/parse_bismark_cov.py" \
  --cov "$OUT_DIR/bismark_paired.cov" \
  --output "$OUT_DIR/bismark_paired_weighted.bed" \
  --deduplicate-cpg-strands
grep -q $'chr1\t10\t11\t0.3' "$OUT_DIR/bismark_paired_weighted.bed"

echo "minimal pipeline ok"
