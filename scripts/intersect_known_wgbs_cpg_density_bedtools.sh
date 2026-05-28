#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="${1:-data}"
OUTPUT_DIR="${2:-$DATA_DIR/bedtools_density_intersections}"
DENSITY_BED="${3:-$DATA_DIR/CpGDensities_W50.bed}"

mkdir -p "$OUTPUT_DIR"

run_one() {
  local input_name="$1"
  local output_name="$2"
  local input_path="$DATA_DIR/$input_name"
  local output_path="$OUTPUT_DIR/$output_name"
  if [[ ! -s "$input_path" ]]; then
    echo "missing $input_name; skipped"
    return
  fi
  bedtools intersect -a "$input_path" -b "$DENSITY_BED" -wa -wb -sorted \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$9}' \
    > "$output_path"
  echo "wrote $(wc -l < "$output_path") rows to $output_path"
}

run_one "GSM1112841_HUES64WT_WGBS_proc.bed" "HUES64WT_CpGsOnly_Chr1.bed"
run_one "GSM1545002_DNMT3A_KO_Early_proc.bed" "HUES64_DNMT3Ako_early_CpGsOnly_Chr1.bed"
run_one "GSM1545003_DNMT3B_KO_Early_proc.bed" "HUES64_DNMT3Bko_early_CpGsOnly_Chr1.bed"
run_one "GSM1545005_DNMT3A_KO_Late_proc.bed" "HUES64_DNMT3Ako_late_CpGsOnly_Chr1.bed"
run_one "GSM1545006_DNMT3B_KO_Late_proc.bed" "HUES64_DNMT3Bko_late_CpGsOnly_Chr1.bed"
run_one "GSM1545004_DKO_Early_proc.bed" "HUES64_DNMT3_dko_early_CpGsOnly_Chr1.bed"
run_one "GSM1545007_DKO_Late_proc.bed" "HUES64_DNMT3_dko_late_CpGsOnly_Chr1.bed"
run_one "GSM3618718_HUES8_WT_WGBS_proc.bed" "HUES8WT_CpGsOnly_Chr1.bed"
run_one "GSM4458672_WGBS_HUES8_DKO_P6_proc.bed" "HUES8_DKO_P6_CpGsOnly_Chr1.bed"
run_one "GSM3618720_HUES8_TKO_WGBS_proc.bed" "HUES8_TKO_CpGsOnly_Chr1.bed"
run_one "GSM3618719_HUES8_QKO_WGBS_proc.bed" "HUES8_QKO_CpGsOnly_Chr1.bed"
run_one "GSM3618721_HUES8_PKO_WGBS_proc.bed" "HUES8_PKO_P0_CpGsOnly_Chr1.bed"
run_one "GSM3662266_HUES8_PKO_P6_WGBS_proc.bed" "HUES8_PKO_P6_CpGsOnly_Chr1.bed"
run_one "GSM4458671_WGBS_HUES8_PKO_P20_proc.bed" "HUES8_PKO_P20_CpGsOnly_Chr1.bed"
run_one "GSM432687_IMR90_WGBS_proc.bed" "IMR90WT_CpGsOnly_Chr1.bed"
