#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
if [[ -z "${PYTHON_BIN:-}" && -x "$HOME/.conda/envs/data_env/bin/python" ]]; then
  PYTHON_BIN="$HOME/.conda/envs/data_env/bin/python"
else
  PYTHON_BIN="${PYTHON_BIN:-python}"
fi
MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mplconfig}"
export MPLCONFIGDIR

log() {
  printf '[%s] [bismark_test] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

finish_step() {
  local path="$1"
  if [[ -e "$path" ]]; then
    log "wrote $path ($(du -sh "$path" | awk '{print $1}'))"
  else
    log "expected output missing: $path"
  fi
}

COV="$SCRIPT_DIR/CB_A_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
SAMPLE="CB_A"
SCOPE="${SCOPE:-chr1}"
DENSITY_WINDOW="${DENSITY_WINDOW:-50}"
MIN_COVERAGE="${MIN_COVERAGE:-1}"
INCLUDE_NON_STANDARD_CONTIGS="${INCLUDE_NON_STANDARD_CONTIGS:-0}"
DATA_DIR="$SCRIPT_DIR/data"
RAW_DIR="$DATA_DIR/raw"
SAMPLE_DIR="$DATA_DIR/$SAMPLE"
WORK_DIR="$SAMPLE_DIR/work"
RESULTS_DIR="$SAMPLE_DIR/results"
FIGURES_DIR="$SAMPLE_DIR/figures"

mkdir -p "$RAW_DIR" "$WORK_DIR" "$RESULTS_DIR" "$FIGURES_DIR" "$MPLCONFIGDIR"

if ! [[ "$MIN_COVERAGE" =~ ^[0-9]+$ ]]; then
  echo "MIN_COVERAGE must be a non-negative integer; got '$MIN_COVERAGE'" >&2
  exit 2
fi

STANDARD_CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
STANDARD_CHROM_REGEX='^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$'
case "$INCLUDE_NON_STANDARD_CONTIGS" in
  1|true|TRUE|yes|YES)
    CGI_STANDARD_CHROMOSOMES_ONLY=0
    CGI_CONTIG_SCOPE="standard + random/alt/unplaced contigs"
    ;;
  *)
    CGI_STANDARD_CHROMOSOMES_ONLY=1
    CGI_CONTIG_SCOPE="standard chromosomes only"
    ;;
esac

if [[ "$SCOPE" == "chr1" ]]; then
  ANALYSIS_SUFFIX="Chr1"
  INDIVIDUAL_CPG_SCOPE="chr1 only"
  STANDARD_CHROMOSOMES_ONLY=1
  FASTA_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz"
  FASTA="$RAW_DIR/chr1.hg38.fa.gz"
elif [[ "$SCOPE" == "genome" ]]; then
  case "$INCLUDE_NON_STANDARD_CONTIGS" in
    1|true|TRUE|yes|YES)
      ANALYSIS_SUFFIX="hg38_all_contigs"
      INDIVIDUAL_CPG_SCOPE="standard + random/alt/unplaced contigs"
      STANDARD_CHROMOSOMES_ONLY=0
      ;;
    *)
      ANALYSIS_SUFFIX="hg38_standard_chromosomes"
      INDIVIDUAL_CPG_SCOPE="standard chromosomes only"
      STANDARD_CHROMOSOMES_ONLY=1
      ;;
  esac
  FASTA_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
  FASTA="$RAW_DIR/hg38.fa.gz"
else
  echo "SCOPE must be 'genome' or 'chr1'; got '$SCOPE'" >&2
  exit 2
fi

log "starting hg38 Bismark pipeline"
log "sample=$SAMPLE individual_cpg_scope=$SCOPE individual_cpg_contigs=$INDIVIDUAL_CPG_SCOPE cpg_island_contigs=$CGI_CONTIG_SCOPE density_window=$DENSITY_WINDOW min_coverage=$MIN_COVERAGE python=$PYTHON_BIN"
log "data_dir=$DATA_DIR"
log "sample_output_dir=$SAMPLE_DIR"
if [[ "$INCLUDE_NON_STANDARD_CONTIGS" != "1" && "$INCLUDE_NON_STANDARD_CONTIGS" != "true" && "$INCLUDE_NON_STANDARD_CONTIGS" != "TRUE" && "$INCLUDE_NON_STANDARD_CONTIGS" != "yes" && "$INCLUDE_NON_STANDARD_CONTIGS" != "YES" ]]; then
  log "default analysis excludes random/alt/unplaced contigs; set INCLUDE_NON_STANDARD_CONTIGS=1 to include them"
fi

if [[ ! -s "$RAW_DIR/cpgIslandExt_hg38.txt.gz" ]]; then
  log "downloading hg38 CpG island annotation"
  curl -L \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz \
    -o "$RAW_DIR/cpgIslandExt_hg38.txt.gz"
else
  log "using cached hg38 CpG island annotation: $RAW_DIR/cpgIslandExt_hg38.txt.gz"
fi

log "converting hg38 CpG island annotation to BED"
gzip -cd "$RAW_DIR/cpgIslandExt_hg38.txt.gz" \
  | awk -v standard_only="$CGI_STANDARD_CHROMOSOMES_ONLY" -v standard_regex="$STANDARD_CHROM_REGEX" \
      'BEGIN{OFS="\t"} !standard_only || $2 ~ standard_regex {print $2,$3,$4,$8,$7,++n}' \
  > "$WORK_DIR/hg38_cpg_islands.bed"
finish_step "$WORK_DIR/hg38_cpg_islands.bed"

log "sorting hg38 CpG island BED"
bedtools sort -i "$WORK_DIR/hg38_cpg_islands.bed" \
  > "$WORK_DIR/hg38_cpg_islands.sorted.bed"
finish_step "$WORK_DIR/hg38_cpg_islands.sorted.bed"

log "parsing Bismark coverage to processed WGBS BED; this is a long step"
"$PYTHON_BIN" "$ROOT_DIR/scripts/parse_bismark_cov.py" \
  --cov "$COV" \
  --output "$WORK_DIR/processed_methylation.unsorted.bed" \
  --min-coverage "$MIN_COVERAGE" \
  --deduplicate-cpg-strands \
  --no-sort \
  --chunksize 1000000
finish_step "$WORK_DIR/processed_methylation.unsorted.bed"

log "sorting processed WGBS BED; this may take a while for whole-genome data"
bedtools sort -i "$WORK_DIR/processed_methylation.unsorted.bed" \
  > "$WORK_DIR/processed_methylation.sorted.bed"
finish_step "$WORK_DIR/processed_methylation.sorted.bed"

log "intersecting processed WGBS with hg38 CpG islands"
bedtools intersect \
  -a "$WORK_DIR/processed_methylation.sorted.bed" \
  -b "$WORK_DIR/hg38_cpg_islands.sorted.bed" \
  -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$8,$9,$10}' \
  > "$WORK_DIR/processed_methylation_cpg_island_intersection.bed"
finish_step "$WORK_DIR/processed_methylation_cpg_island_intersection.bed"

log "aggregating CpG-island-level methylation"
"$PYTHON_BIN" "$ROOT_DIR/scripts/aggregate_cgi_level.py" \
  --cgi-intersection "$WORK_DIR/processed_methylation_cpg_island_intersection.bed" \
  --output "$RESULTS_DIR/cpg_island_methylation_aggregate.csv"
finish_step "$RESULTS_DIR/cpg_island_methylation_aggregate.csv"

log "generating CpG-island-level figures and summary"
"$PYTHON_BIN" "$ROOT_DIR/scripts/bivariate_histogram_human_wt.py" \
  --island-files "$RESULTS_DIR/cpg_island_methylation_aggregate.csv" \
  --output-json "$RESULTS_DIR/cpg_island_methylation_by_length_summary.json" \
  --output-plot "$FIGURES_DIR/cpg_island_methylation_classes_by_length.png" \
  --matlab-style-plot "$FIGURES_DIR/cpg_island_methylation_heatmap_and_classes_by_length.png"
finish_step "$RESULTS_DIR/cpg_island_methylation_by_length_summary.json"
finish_step "$FIGURES_DIR/cpg_island_methylation_classes_by_length.png"
finish_step "$FIGURES_DIR/cpg_island_methylation_heatmap_and_classes_by_length.png"

if [[ ! -s "$FASTA" ]]; then
  log "downloading hg38 FASTA for $SCOPE individual CpG analysis"
  curl -L "$FASTA_URL" -o "$FASTA"
else
  log "using cached FASTA: $FASTA"
fi

log "calculating CpG density from FASTA; progress reports one line per included chromosome/contig"
if [[ "$SCOPE" == "chr1" ]]; then
  "$PYTHON_BIN" "$ROOT_DIR/scripts/calculate_cpg_density_from_fasta.py" \
    --fasta "$FASTA" \
    --window "$DENSITY_WINDOW" \
    --start-offset 1 \
    --chroms chr1 \
    --output "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.unsorted.bed"
elif [[ "$STANDARD_CHROMOSOMES_ONLY" -eq 1 ]]; then
  "$PYTHON_BIN" "$ROOT_DIR/scripts/calculate_cpg_density_from_fasta.py" \
    --fasta "$FASTA" \
    --window "$DENSITY_WINDOW" \
    --start-offset 1 \
    --chroms "${STANDARD_CHROMS[@]}" \
    --output "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.unsorted.bed"
else
  "$PYTHON_BIN" "$ROOT_DIR/scripts/calculate_cpg_density_from_fasta.py" \
    --fasta "$FASTA" \
    --window "$DENSITY_WINDOW" \
    --start-offset 1 \
    --output "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.unsorted.bed"
fi
finish_step "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.unsorted.bed"

log "sorting CpG density BED"
bedtools sort -i "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.unsorted.bed" \
  > "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.bed"
finish_step "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.bed"

if [[ "$SCOPE" == "chr1" ]]; then
  log "filtering processed WGBS to chr1 for quick mode"
  awk 'BEGIN{OFS="\t"} $1=="chr1" {print $0}' \
    "$WORK_DIR/processed_methylation.sorted.bed" \
    > "$WORK_DIR/processed_methylation_${ANALYSIS_SUFFIX}.sorted.bed"
  finish_step "$WORK_DIR/processed_methylation_${ANALYSIS_SUFFIX}.sorted.bed"
  WGBS_FOR_CPG="$WORK_DIR/processed_methylation_${ANALYSIS_SUFFIX}.sorted.bed"
elif [[ "$STANDARD_CHROMOSOMES_ONLY" -eq 1 ]]; then
  log "filtering processed WGBS to standard chromosomes"
  awk -v standard_regex="$STANDARD_CHROM_REGEX" 'BEGIN{OFS="\t"} $1 ~ standard_regex {print $0}' \
    "$WORK_DIR/processed_methylation.sorted.bed" \
    > "$WORK_DIR/processed_methylation_${ANALYSIS_SUFFIX}.sorted.bed"
  finish_step "$WORK_DIR/processed_methylation_${ANALYSIS_SUFFIX}.sorted.bed"
  WGBS_FOR_CPG="$WORK_DIR/processed_methylation_${ANALYSIS_SUFFIX}.sorted.bed"
else
  WGBS_FOR_CPG="$WORK_DIR/processed_methylation.sorted.bed"
fi

log "intersecting processed WGBS with CpG density BED for individual CpG analysis"
bedtools intersect \
  -a "$WGBS_FOR_CPG" \
  -b "$WORK_DIR/cpg_density_${ANALYSIS_SUFFIX}_window_${DENSITY_WINDOW}.bed" \
  -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$9}' \
  > "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}.bed"
finish_step "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}.bed"

log "summarizing individual CpG density methylation classes/statistics"
"$PYTHON_BIN" "$ROOT_DIR/scripts/read_plot_data_wt.py" \
  --input "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}.bed" \
  --output-prefix "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}" \
  --summary-json "$RESULTS_DIR/individual_cpg_methylation_by_density_summary.json" \
  --plot "$FIGURES_DIR/individual_cpg_methylation_classes_by_density.png"
finish_step "$RESULTS_DIR/individual_cpg_methylation_by_density_summary.json"
finish_step "$FIGURES_DIR/individual_cpg_methylation_classes_by_density.png"

log "generating mean/median/Hill CpG-density plot"
"$PYTHON_BIN" "$ROOT_DIR/scripts/plot_figure5a_cpg_density.py" \
  --inputs "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}.bed" \
  --labels "$SAMPLE" \
  --output "$FIGURES_DIR/individual_cpg_mean_median_hill_by_density.png" \
  --summary-json "$RESULTS_DIR/individual_cpg_mean_median_hill_by_density_summary.json"
finish_step "$FIGURES_DIR/individual_cpg_mean_median_hill_by_density.png"

log "generating switch-parameter plot and CSV"
"$PYTHON_BIN" "$ROOT_DIR/scripts/plot_switch_parameters.py" \
  --inputs "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}.bed" \
  --labels "$SAMPLE" \
  --cohorts "$SAMPLE" \
  --ko-types None \
  --output "$FIGURES_DIR/individual_cpg_switch_parameters.png" \
  --summary-csv "$RESULTS_DIR/individual_cpg_switch_parameters.csv"
finish_step "$RESULTS_DIR/individual_cpg_switch_parameters.csv"
finish_step "$FIGURES_DIR/individual_cpg_switch_parameters.png"

log "fitting CME model"
"$PYTHON_BIN" "$ROOT_DIR/scripts/fit_cme_methylation.py" \
  --data-struct "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}.npz" \
  --output-prefix "$RESULTS_DIR/cme_fit_individual_cpg_${ANALYSIS_SUFFIX}" \
  --n-cpg 27 \
  --iterations 60 \
  --seed 20260522
finish_step "$RESULTS_DIR/cme_fit_individual_cpg_${ANALYSIS_SUFFIX}.summary.json"

log "generating CME fit comparison plot"
"$PYTHON_BIN" "$ROOT_DIR/scripts/plot_cme_fit_comparison.py" \
  --data-struct "$RESULTS_DIR/individual_cpg_methylation_density_${ANALYSIS_SUFFIX}.npz" \
  --fit-result "$RESULTS_DIR/cme_fit_individual_cpg_${ANALYSIS_SUFFIX}.npz" \
  --label "$SAMPLE" \
  --output "$FIGURES_DIR/cme_fit_model_comparison.png" \
  --summary-json "$RESULTS_DIR/cme_fit_model_comparison_summary.json"
finish_step "$FIGURES_DIR/cme_fit_model_comparison.png"

log "running temporal CME example and plots"
"$PYTHON_BIN" "$ROOT_DIR/scripts/run_temporal_cme_example.py" \
  --fit-result "$RESULTS_DIR/cme_fit_individual_cpg_${ANALYSIS_SUFFIX}.npz" \
  --output-prefix "$RESULTS_DIR/temporal_cme" \
  --plot-prefix "$FIGURES_DIR/temporal_cme" \
  --n-cpg 27
finish_step "$RESULTS_DIR/temporal_cme.switch_summary.json"
finish_step "$FIGURES_DIR/temporal_cme_mean_methylation.png"
finish_step "$FIGURES_DIR/temporal_cme_inset.png"

log "writing consolidated metrics summary"
"$PYTHON_BIN" - "$RESULTS_DIR" "$SAMPLE" "$ANALYSIS_SUFFIX" "$INDIVIDUAL_CPG_SCOPE" "$CGI_CONTIG_SCOPE" "$MIN_COVERAGE" <<'PY'
import csv
import json
import sys
from pathlib import Path

results = Path(sys.argv[1])
sample = sys.argv[2]
analysis_suffix = sys.argv[3]
individual_cpg_scope = sys.argv[4]
cgi_contig_scope = sys.argv[5]
min_coverage = int(sys.argv[6])
cgi = json.loads((results / "cpg_island_methylation_by_length_summary.json").read_text())
cpg = json.loads((results / "individual_cpg_methylation_by_density_summary.json").read_text())
switch = list(csv.DictReader((results / "individual_cpg_switch_parameters.csv").open()))[0]
fit = json.loads((results / f"cme_fit_individual_cpg_{analysis_suffix}.summary.json").read_text())
temporal = json.loads((results / "temporal_cme.switch_summary.json").read_text())
island = next(iter(cgi.values()))
summary = {
    "sample": sample,
    "genome": "hg38",
    "cpg_island_contig_scope": cgi_contig_scope,
    "min_coverage": min_coverage,
    "cgi_level": {
        "LD50_bp": island["cgi_crossover_bp"],
        "cgi_crossover_bp": island["cgi_crossover_bp"],
    },
    "individual_cpg": {
        "scope": analysis_suffix,
        "contig_scope": individual_cpg_scope,
        "ED50": cpg["ED50"],
        "DirectED50": cpg["DirectED50"],
        "DirectSlope": cpg["DirectSlope"],
        "HillK": cpg["HillKLog"],
        "HillN": cpg["HillNLog"],
        "HillN_one_param": cpg["HillN"],
    },
    "switch_plot_record": switch,
    "cme_fit": {
        "SSD": fit["SSD"],
        "Parameters": fit["Parameters"],
        "DistLength": fit["DistLength"],
    },
    "temporal_base_switch": {
        "columns": temporal["columns"],
        "values": temporal["BaseSwitchValues"],
    },
}
(results / "pipeline_metrics_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(results / "pipeline_metrics_summary.csv").write_text(
    "\n".join(
        [
            "metric,value",
            f"Individual_CpG_scope,{analysis_suffix}",
            f"Individual_CpG_contig_scope,{individual_cpg_scope}",
            f"CpG_island_contig_scope,{cgi_contig_scope}",
            f"Min_coverage,{min_coverage}",
            f"CGI_LD50_bp,{summary['cgi_level']['LD50_bp']}",
            f"CpG_ED50,{cpg['ED50']}",
            f"CpG_DirectED50,{cpg['DirectED50']}",
            f"CpG_DirectSlope,{cpg['DirectSlope']}",
            f"CpG_HillK,{cpg['HillKLog']}",
            f"CpG_HillN,{cpg['HillNLog']}",
            f"CME_SSD,{fit['SSD']}",
        ]
    )
    + "\n"
)
print("wrote metrics summary")
PY
finish_step "$RESULTS_DIR/pipeline_metrics_summary.csv"
finish_step "$RESULTS_DIR/pipeline_metrics_summary.json"
log "pipeline complete"
