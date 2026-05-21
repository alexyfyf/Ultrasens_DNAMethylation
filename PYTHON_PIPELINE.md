# Python Pipeline Guide

This guide describes how to run the Python rewrite in the same order as the original MATLAB README:

1. CpG island-level analysis and plots.
2. Individual CpG-level density analysis and plots.
3. CME model fitting and model-derived plots.

Reusable code lives in `python/ultrasens/`. Command-line scripts live in `scripts/`.

## Environment

Activate the project environment before running the commands below:

```bash
conda activate data_env
```

If you are creating a fresh environment, install the Python packages plus the external genomics tools:

```bash
conda install -c conda-forge numpy pandas scipy matplotlib bedtools bedops seqkit
```

`matplotlib` may need a writable cache directory on some systems:

```bash
export MPLCONFIGDIR=/private/tmp/mplconfig
```

## Input Conventions

The Python scripts use these normalized file formats:

- Processed WGBS BED: tab-delimited `chr start end WGBS`, where `WGBS` is a methylation fraction from 0 to 1 and each CpG dyad appears once.
- CGI annotation BED: tab-delimited `Chrom CpGstart CpGEnd CpGNum CGIlen CGIno`.
- WGBS/CGI intersection: tab-delimited `Chrom CpGstart CpGEnd WGBS CpGNum CGIlen CGIno`.
- CpG-density BED: tab-delimited `chr start end Density`.
- WGBS/CpG-density intersection: tab-delimited `chr start end WGBS Density`.

See `DATA_PREP.md` for download and conversion details.

## One-To-One Script Map

- `CpGDensity_Calc.m` -> `scripts/calculate_cpg_density.py`
- `WGBS_CpGIntersect_AllData.command` -> `scripts/intersect_wgbs_cpg_density.py`, `scripts/intersect_wgbs_cpg_density_all.py`, or `scripts/intersect_wgbs_cpg_density_all_bedtools.sh`
- Opposite-strand CpG dyad deduplication for WIG/BED/COV inputs -> `scripts/deduplicate_cpg_strands.py` or `scripts/parse_bismark_cov.py --deduplicate-cpg-strands`
- `BivariateHistogram_HumanWT_example.m` -> `scripts/bivariate_histogram_human_wt.py`
- `IndividualCGIMethylationChange_HUES64_example.m` -> `scripts/individual_cgi_methylation_change.py`
- `ReadPlotData_WT_example.m` -> `scripts/read_plot_data_wt.py`, `scripts/read_plot_data_wt_all.py`, `scripts/plot_figure5a_cpg_density.py`, and `scripts/plot_switch_parameters.py`
- `CMEFitting/LoopFit.m`, `CMEFitting/Fit_CME_Methylation_PS.m`, `CMEFitting/CMEModel.m` -> `scripts/fit_cme_methylation.py` and `python/ultrasens/cme.py`
- `ExampleCMEFit/CallFitting.m` -> `scripts/fit_example_cme.py`
- `RunTemporalCMEExample/*.m` -> `scripts/run_temporal_cme_example.py`

## A. CpG Island-Level Analysis

### A1. Prepare WGBS/CGI Intersections

For the HUES64 WT notebook example:

```bash
bash scripts/prepare_ipynb_inputs.sh
```

This writes:

- `data/CpGIsAnn_hg19_dec14.bed`
- `data/CpGIsAnn_hg19_chr1.bed`
- `data/GSM1112841_HUES64WT_WGBS_proc.bed`
- `data/HUES64WT_WGBS_CGI_int.bed`

For any other sample, first make a processed WGBS BED and then intersect it with the CGI annotation:

```bash
bedtools intersect \
  -a data/SAMPLE_WGBS_proc.bed \
  -b data/CpGIsAnn_hg19_dec14.bed \
  -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$8,$9,$10}' \
  > data/SAMPLE_WGBS_CGI_int.bed
```

### A2. Aggregate WGBS Rows To Island-Level Values

This reproduces the notebook output used by the MATLAB island-level scripts:

```bash
python scripts/aggregate_cgi_level.py \
  --cgi-intersection data/HUES64WT_WGBS_CGI_int.bed \
  --output data/IslandLvl_agg_HUES64WT
```

Run the same command for each sample you want to compare. The output columns are:

```text
CGIno,CpGNum,CGIlen,WGBS
```

### A3. WT Bivariate Histograms And EL50/ED50-Style Crossovers

This replaces `BivariateHistogram_HumanWT_example.m`.

```bash
python scripts/bivariate_histogram_human_wt.py \
  --island-files \
    data/IslandLvl_agg_HUES64WT \
    data/IslandLvl_agg_HUES8WT \
    data/IslandLvl_agg_IMR90WT \
  --output-json data/figure_checks/bivariate_human_wt_summary.json \
  --output-plot data/figure_checks/bivariate_human_wt_classes.png \
  --matlab-style-plot data/figure_checks/bivariate_human_wt_matlab_style.png
```

Main outputs:

- `*_summary.json`: crossover values and binned fractions.
- `*_classes.png`: hypo/hyper/inter island count curves.
- `*_matlab_style.png`: MATLAB-style bivariate heatmap and class curves.

### A4. CGI Methylation Change Between WT And KO

This replaces `IndividualCGIMethylationChange_HUES64_example.m`.

```bash
python scripts/individual_cgi_methylation_change.py \
  --file-a data/IslandLvl_agg_HUES64_DNMT3_dko_early \
  --file-b data/IslandLvl_agg_HUES64WT \
  --output data/HUES64_DNMT3_dko_early_CGI_change.csv \
  --plot data/HUES64_DNMT3_dko_early_CGI_change.png \
  --title "HUES64 DNMT3 DKO early"
```

Here `DeltaMeth = file-a WGBS - file-b WGBS`, so the example above reports KO minus WT. The plot names the classes explicitly:

- Class 1: strong decrease, `DeltaMeth <= -0.6`
- Class 2: moderate decrease, `-0.6 < DeltaMeth <= -0.2`
- Class 3: minimal change, `-0.2 < DeltaMeth <= 0.2`
- Class 4: moderate increase, `0.2 < DeltaMeth <= 0.6`
- Class 5: strong increase, `DeltaMeth > 0.6`

## B. Individual CpG-Level Analysis

### B1. Calculate Local CpG Density

This replaces `CpGDensity_Calc.m`.

Input is a CpG coordinate CSV from `seqkit locate`, with at least a `start` column. For the paper-style individual CpG analysis, use a 50 bp window:

```bash
python scripts/calculate_cpg_density.py \
  --cpg-csv data/chr1_cpgs.csv \
  --window 50 \
  --chrom chr1 \
  --output data/CpGDensities_W50.bed
```

The output is a sorted BED-like table:

```text
chr1  start  end  Density
```

### B2. Intersect WGBS With CpG Density

This replaces `WGBS_CpGIntersect_AllData.command`.

Single sample, pure Python:

```bash
python scripts/intersect_wgbs_cpg_density.py \
  --wgbs-bed data/GSM3618718_HUES8_WT_WGBS_proc.bed \
  --density-bed data/CpGDensities_W50.bed \
  --output data/HUES8WT_CpGsOnly_Chr1.bed
```

All original repo dataset names, pure Python:

```bash
python scripts/intersect_wgbs_cpg_density_all.py \
  --data-dir data \
  --density-bed data/CpGDensities_W50.bed \
  --output-dir data/wgbs_density_intersections
```

All original repo dataset names, using installed `bedtools`:

```bash
bash scripts/intersect_wgbs_cpg_density_all_bedtools.sh \
  data \
  data/bedtools_density_intersections \
  data/CpGDensities_W50.bed
```

Each output has:

```text
chr  start  end  WGBS  Density
```

### B3. Classify CpGs And Fit Density Switches

This replaces the core of `ReadPlotData_WT_example.m`.

Single sample:

```bash
python scripts/read_plot_data_wt.py \
  --input data/bedtools_density_intersections/HUES8WT_CpGsOnly_Chr1.bed \
  --output-prefix data/read_plot_wt_checks/HUES8WT_CpGsOnly_Chr1 \
  --summary-json data/read_plot_wt_checks/HUES8WT_CpGsOnly_Chr1_summary.json \
  --plot data/read_plot_wt_checks/HUES8WT_CpGsOnly_Chr1_curves.png
```

Important output fields:

- `Hypo`, `Hyper`, `Inter`: fraction of CpGs by methylation class versus local density.
- `MeanMeth`, `MedMeth`, `Prc25`, `Prc75`: methylation summaries versus density.
- `ED50`: direct hypo/hyper crossover estimate.
- `HillKLog`, `HillNLog`: Hill threshold and steepness from the log-transform method.
- `DirectED50`, `DirectSlope`: direct crossover and slope values.

Multiple WT samples:

```bash
python scripts/read_plot_data_wt_all.py \
  --inputs \
    data/bedtools_density_intersections/HUES64WT_CpGsOnly_Chr1.bed \
    data/bedtools_density_intersections/HUES8WT_CpGsOnly_Chr1.bed \
    data/bedtools_density_intersections/IMR90WT_CpGsOnly_Chr1.bed \
  --labels HUES64WT HUES8WT IMR90WT \
  --output-dir data/read_plot_wt_combined \
  --summary-json data/read_plot_wt_combined/read_plot_data_wt_summary.json \
  --plot data/read_plot_wt_combined/read_plot_data_wt_combined.png
```

### B4. Figure 5A-Style Mean/Median/IQR Plot

This plots mean methylation, median/IQR methylation, and the Hill fit to mean methylation:

```bash
python scripts/plot_figure5a_cpg_density.py \
  --inputs \
    data/bedtools_density_intersections/HUES8WT_CpGsOnly_Chr1.bed \
    data/bedtools_density_intersections/HUES8_TKO_CpGsOnly_Chr1.bed \
    data/bedtools_density_intersections/HUES8_DKO_P6_CpGsOnly_Chr1.bed \
  --labels HUES8 "TET1-3 (TKO)" "DNMT3AB (DKO)" \
  --output data/figure5a_cpg_density.png \
  --summary-json data/figure5a_cpg_density.json
```

### B5. Figure 5B/C-Style Switch Parameter Plot

This plots Hill `n` versus `K` and direct slope versus `ED50`.

From CpG-density input files:

```bash
python scripts/plot_switch_parameters.py \
  --inputs \
    data/bedtools_density_intersections/HUES64WT_CpGsOnly_Chr1.bed \
    data/bedtools_density_intersections/HUES8WT_CpGsOnly_Chr1.bed \
  --labels HUES64 HUES8 \
  --cohorts HUES64 HUES8 \
  --ko-types None None \
  --output data/figure5bc_switch_parameters.png \
  --summary-csv data/figure5bc_switch_parameters.csv
```

Or from a curated CSV table with optional error bars:

```bash
python scripts/plot_switch_parameters.py \
  --table data/figure5bc_switch_parameters.csv \
  --cohort-order HUES64 HUES8 \
  --output data/figure5bc_switch_parameters.png
```

Required table columns:

```text
label,cohort,ko_type,HillK,HillN,DirectED50,DirectSlope
```

Optional error columns include `HillKErr`, `HillNErr`, `DirectED50Err`, `DirectSlopeErr`, or asymmetric `HillKErrLow/HillKErrHigh` style pairs.

## C. CME Fitting And Model Plots

### C1. Fit The CME Model

This replaces `CMEFitting/LoopFit.m` and `CMEFitting/Fit_CME_Methylation_PS.m`.

Fast smoke fit:

```bash
python scripts/fit_cme_methylation.py \
  --data-struct data/read_plot_wt_checks/HUES8WT_CpGsOnly_Chr1.npz \
  --output-prefix data/cme_fit_checks/Fit_HUES8WT_smoke \
  --iterations 10
```

Bounded local SciPy fit from the MATLAB default starting parameters:

```bash
python scripts/fit_cme_methylation.py \
  --optimizer lbfgsb \
  --maxiter 300 \
  --maxfun 2000 \
  --data-struct data/read_plot_wt_checks/HUES8WT_CpGsOnly_Chr1.npz \
  --output-prefix data/cme_fit_checks/Fit_HUES8WT_lbfgsb
```

Global differential-evolution fit:

```bash
python scripts/fit_cme_methylation.py \
  --optimizer scipy \
  --maxiter 35 \
  --popsize 8 \
  --data-struct data/read_plot_wt_checks/HUES8WT_CpGsOnly_Chr1.npz \
  --output-prefix data/cme_fit_checks/Fit_HUES8WT_full_scipy \
  --seed 20260520
```

Outputs:

- `<prefix>.npz` and `<prefix>.json`: `Parameters`, `DistLength`, `BestFitVector`, `SSD`.
- `<prefix>.summary.json`: optimizer settings and convergence metadata.

### C2. Plot Collaborative Versus Standard CME Fit

This generates the fitted collaborative model, the no-collaboration standard model, and a schematic with arrow widths scaled by fitted rates:

```bash
python scripts/plot_cme_fit_comparison.py \
  --data-struct data/read_plot_wt_checks/HUES8WT_CpGsOnly_Chr1.npz \
  --fit-result data/cme_fit_checks/Fit_HUES8WT_lbfgsb.npz \
  --label HUES8 \
  --output data/cme_fit_checks/Fit_HUES8WT_model_comparison.png \
  --summary-json data/cme_fit_checks/Fit_HUES8WT_model_comparison.summary.json
```

The standard model is created by setting collaborative rates `k5:k12` to zero and rerunning the CME model.

### C3. Run Temporal CME Example

This replaces `RunTemporalCMEExample/Plot_ExampleTemporal.m` and related model files:

```bash
python scripts/run_temporal_cme_example.py \
  --fit-result data/cme_fit_checks/Fit_HUES8WT_lbfgsb.npz \
  --output-prefix data/cme_fit_checks/Temporal_HUES8WT \
  --plot-prefix data/cme_fit_checks/Temporal_HUES8WT
```

Outputs include:

- `Temporal_HUES8WT.npz`
- `Temporal_HUES8WT.json`
- `Temporal_HUES8WT.switch_summary.json`
- `Temporal_HUES8WT_mean_methylation.png`
- `Temporal_HUES8WT_inset.png`

## D. End-To-End Placeholder Pipeline

For small tests or new data, `scripts/run_pipeline.py` chains the major steps:

```bash
python scripts/run_pipeline.py \
  --cpg-csv /path/to/chr1_cpgs.csv \
  --wgbs-bed /path/to/processed_wgbs_chr1.bed \
  --cgi-intersection /path/to/wgbs_cgi_intersection.bed \
  --output-dir /path/to/output \
  --chrom chr1 \
  --density-window 50 \
  --n-cpg 27 \
  --fit-iterations 60
```

Run the included minimal fixture pipeline:

```bash
PYTHON_BIN=$(which python) MPLCONFIGDIR=/private/tmp/mplconfig bash tests/run_minimal_pipeline.sh
```

Expected final line:

```text
minimal pipeline ok
```
