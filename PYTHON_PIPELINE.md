# Python Pipeline

This is the minimal Python rewrite of the necessary MATLAB analysis scripts. Reusable logic lives in `python/ultrasens/`; command-line entry points live in `scripts/`.

## One-To-One Rewrites

- `CpGDensity_Calc.m` -> `scripts/calculate_cpg_density.py`
- `WGBS_CpGIntersect_AllData.command` -> `scripts/intersect_wgbs_cpg_density.py`
- `BivariateHistogram_HumanWT_example.m` -> `scripts/bivariate_histogram_human_wt.py`
- `IndividualCGIMethylationChange_HUES64_example.m` -> `scripts/individual_cgi_methylation_change.py`
- `ReadPlotData_WT_example.m` -> `scripts/read_plot_data_wt.py`
- `CMEFitting/LoopFit.m`, `CMEFitting/Fit_CME_Methylation_PS.m`, `CMEFitting/CMEModel.m`, `CMEFitting/CpGDensities_Function.m` -> `scripts/fit_cme_methylation.py` plus `python/ultrasens/cme.py`
- `RunTemporalCMEExample/*.m` -> `scripts/run_temporal_cme_example.py` plus `python/ultrasens/cme.py`

`ExampleCMEFit/` is an older duplicate fitting example and is not rewritten separately.

## Full Pipeline With Placeholders

```bash
python scripts/run_pipeline.py \
  --cpg-csv /path/to/chr1_cpgs.csv \
  --wgbs-bed /path/to/processed_wgbs_chr1.bed \
  --cgi-intersection /path/to/wgbs_cgi_intersection.csv \
  --output-dir /path/to/output \
  --chrom chr1 \
  --density-window 50 \
  --n-cpg 27 \
  --fit-iterations 60
```

Expected placeholder input formats:

- `--cpg-csv`: CSV with a `start` column, as from `seqkit locate`.
- `--wgbs-bed`: tab-delimited `chr start end methyl_fraction`.
- `--cgi-intersection`: `chr,start,end,WGBS,CGIno,CpGNum,CGIlen`.

For Bismark coverage output, convert `.cov`/`.cov.gz` first:

```bash
python scripts/parse_bismark_cov.py \
  --cov /path/to/sample.bismark.cov.gz \
  --output /path/to/sample_wgbs_proc.bed
```

The Bismark input is interpreted as 1-based `chrom start end methylation_percentage count_methylated count_unmethylated`.
The output is the pipeline's 0-based half-open `chr start end WGBS` format. By default `WGBS`
is computed from the exact counts; use `--methylation-source percentage` to use the percentage column.

The pipeline writes CpG-density BED, WGBS-density intersections, CGI summaries, `DataStruct`-style `.npz/.json` files for CME fitting, fitted CME parameters, and temporal CME output.

`scripts/run_temporal_cme_example.py` now matches the MATLAB temporal plotting script more closely: it exports `BaseSwitchValues`, `TemporalSwitchValues`, `KStar`, and `NStar`, and writes `*.switch_summary.json`. When `matplotlib` is installed, it also writes:

- `<plot-prefix>_mean_methylation.png`
- `<plot-prefix>_inset.png`

Use `--plot-prefix /path/to/ParameterSweepsExample` to control those filenames, or `--no-plots` to skip plotting.

## Minimal Test

```bash
PYTHON_BIN=/Users/yan.a/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  bash tests/run_minimal_pipeline.sh
```
