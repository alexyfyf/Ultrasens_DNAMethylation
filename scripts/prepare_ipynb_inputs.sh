#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${DATA_DIR:-$ROOT_DIR/data}"
RAW_DIR="$DATA_DIR/raw"

mkdir -p "$DATA_DIR" "$RAW_DIR"

UCSC_CGI_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz"
GSM1112841_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1112841&file=GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.wig.gz&format=file"

CGI_RAW="$RAW_DIR/cpgIslandExt_hg19.txt.gz"
WIG_GZ="$DATA_DIR/GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.wig.gz"
WIG_BED="$DATA_DIR/GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.bed"
WGBS_PROC="$DATA_DIR/GSM1112841_HUES64WT_WGBS_proc.bed"
CGI_BED="$DATA_DIR/CpGIsAnn_hg19_dec14.bed"
CGI_CHR1_BED="$DATA_DIR/CpGIsAnn_hg19_chr1.bed"
IPYNB_INPUT="$DATA_DIR/HUES64WT_WGBS_CGI_int.bed"

command -v curl >/dev/null
command -v gzip >/dev/null
command -v awk >/dev/null
command -v sort >/dev/null
command -v wig2bed >/dev/null
command -v bedtools >/dev/null

if [[ ! -s "$CGI_RAW" ]]; then
  curl -L "$UCSC_CGI_URL" -o "$CGI_RAW"
fi

# Notebook format: Chrom, CpGstart, CpGEnd, CpGNum, CGIlen, CGIno
# UCSC cpgIslandExt's "name" field is "CpG: N", so whitespace splitting makes:
# $2 chrom, $3 chromStart, $4 chromEnd, $7 length, $8 cpgNum.
gzip -cd "$CGI_RAW" \
  | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$8,$7,++n}' \
  > "$CGI_BED"

awk 'BEGIN{OFS="\t"} $1=="chr1" {print $0}' "$CGI_BED" > "$CGI_CHR1_BED"

if [[ ! -s "$WIG_GZ" ]]; then
  curl -L -C - -o "$WIG_GZ" "$GSM1112841_URL"
fi

if [[ ! -s "$WIG_BED" ]]; then
  gzip -cd "$WIG_GZ" | wig2bed > "$WIG_BED"
fi

# Repo README format: keep every second WIG-derived row and write chr/start/end/methyl_fraction.
awk 'BEGIN{OFS="\t"} NR % 2 == 0 {print $1,$2,$3,$5}' "$WIG_BED" \
  | sort -k1,1 -k2,2n \
  > "$WGBS_PROC"

# Notebook input format: Chrom, CpGstart, CpGEnd, WGBS, CpGNum, CGIlen, CGIno
bedtools intersect -a "$WGBS_PROC" -b "$CGI_BED" -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$8,$9,$10}' \
  > "$IPYNB_INPUT"

echo "Wrote $CGI_BED"
echo "Wrote $CGI_CHR1_BED"
echo "Wrote $WGBS_PROC"
echo "Wrote $IPYNB_INPUT"

