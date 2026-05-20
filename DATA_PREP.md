# Input Preparation For The Notebook

This documents the commands used to prepare the inputs expected by `Analysis of Human ESCs (HUES64WT_example).ipynb`.

## Outputs

- `data/CpGIsAnn_hg19_dec14.bed`
  - Full hg19 UCSC CpG island annotation.
  - No header, tab-delimited.
  - Columns: `Chrom`, `CpGstart`, `CpGEnd`, `CpGNum`, `CGIlen`, `CGIno`.

- `data/HUES64WT_WGBS_CGI_int.bed`
  - Intersected HUES64 WT WGBS methylation proportions and CpG island annotation.
  - No header, tab-delimited.
  - Columns: `Chrom`, `CpGstart`, `CpGEnd`, `WGBS`, `CpGNum`, `CGIlen`, `CGIno`.

## Reproducible Script

Run:

```bash
bash scripts/prepare_ipynb_inputs.sh
```

The script downloads:

- UCSC hg19 `cpgIslandExt.txt.gz`
- GEO `GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.wig.gz`

It then runs:

```bash
gzip -cd data/GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.wig.gz \
  | wig2bed \
  > data/GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.bed

awk 'BEGIN{OFS="\t"} NR % 2 == 0 {print $1,$2,$3,$5}' \
  data/GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.bed \
  | sort -k1,1 -k2,2n \
  > data/GSM1112841_HUES64WT_WGBS_proc.bed

bedtools intersect \
  -a data/GSM1112841_HUES64WT_WGBS_proc.bed \
  -b data/CpGIsAnn_hg19_dec14.bed \
  -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$8,$9,$10}' \
  > data/HUES64WT_WGBS_CGI_int.bed
```

## Current Validation

- `data/CpGIsAnn_hg19_dec14.bed`: 30,344 CpG islands, matching the notebook's expected whole-genome annotation count.
- `data/HUES64WT_WGBS_CGI_int.bed`: 2,004,040 intersected WGBS rows.
- Unique CGIs with WGBS values: 27,265, matching the notebook output.

To repeat the notebook input check without the plotting/Jupyter dependencies:

```bash
python scripts/check_ipynb_inputs.py
```

This writes notebook-equivalent CSV exports to `data/ipynb_run_check/`.

## Additional WT Bivariate Inputs

The bivariate CGI-length analysis also needs HUES8 WT and IMR90 WT island-level files. I prepared:

- `data/IslandLvl_agg_HUES64WT`
- `data/IslandLvl_agg_HUES8WT`
- `data/IslandLvl_agg_IMR90WT`

Sources:

- HUES8 WT: `GSM3618718_HUES8_WT_WGBS.bed.gz`
- IMR90 WT, paper/repo-matching input: `GSM432687_UCSD.IMR90.Bisulfite-Seq.combined.wig.gz`
- IMR90 WT, alternate input tested: `GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.bed.gz`

The normalized WGBS BED files are:

- `data/GSM3618718_HUES8_WT_WGBS_proc.bed`
- `data/GSM432687_IMR90_WGBS_proc.bed`
- `data/GSM1204464_IMR90_WGBS_proc.bed`
- `data/GSM1204464_IMR90_WGBS_ratio_proc.bed`

The original repository command `WGBS_CpGIntersect_AllData.command` points IMR90 WT to
`GSM432687_IMR90_WGBS_proc.bed`. `GSM1204464` is another IMR90 WGBS sample and gives a
different crossover, so it should not be used for reproducing the paper's Figure 2
IMR90 WT panel.

Prepared island-level IMR90 files:

- `data/IslandLvl_agg_IMR90WT`: canonical paper-matching island file from `GSM432687`.
- `data/IslandLvl_agg_IMR90WT_GSM432687`: same paper-matching file, with accession in name.
- `data/IslandLvl_agg_IMR90WT_GSM1204464_score`: alternate file using the rounded BED score.
- `data/IslandLvl_agg_IMR90WT_GSM1204464_ratio`: alternate file using exact methylated/total ratio.

The MATLAB-exact bivariate results are:

- HUES64 WT: 351 bp
- HUES8 WT: 321 bp
- IMR90 WT, `GSM432687`: 214 bp
- IMR90 WT, alternate `GSM1204464`: 306 bp

The combined result files are in `data/figure_checks/`:

- `bivariate_human_wt_summary.json`
- `bivariate_human_wt_classes.png`
- `bivariate_human_wt_matlab_style.png`
