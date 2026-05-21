# Data Download And Input Preparation

This document records how to prepare the files consumed by the Python pipeline. Large generated and downloaded files live under `data/`, which is intentionally ignored by Git.

## Required Tools

Install these before preparing data:

```bash
conda install -c conda-forge bedtools bedops seqkit
```

The scripts also require the Python packages in `requirements.txt`.

## Directory Layout

Recommended layout:

```text
data/
  raw/
  CpGIsAnn_hg19_dec14.bed
  CpGIsAnn_hg19_chr1.bed
  *_WGBS_proc.bed
  *_WGBS_CGI_int.bed
  IslandLvl_agg_*
  CpGDensities_W50.bed
  bedtools_density_intersections/
```

## 1. CpG Island Annotation

Download the hg19 UCSC CpG island annotation:

```bash
mkdir -p data/raw
curl -L \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz \
  -o data/raw/cpgIslandExt_hg19.txt.gz
```

Convert it to the notebook/pipeline format:

```bash
gzip -cd data/raw/cpgIslandExt_hg19.txt.gz \
  | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$8,$7,++n}' \
  > data/CpGIsAnn_hg19_dec14.bed

awk 'BEGIN{OFS="\t"} $1=="chr1" {print $0}' \
  data/CpGIsAnn_hg19_dec14.bed \
  > data/CpGIsAnn_hg19_chr1.bed
```

Output columns:

```text
Chrom  CpGstart  CpGEnd  CpGNum  CGIlen  CGIno
```

The helper script `scripts/prepare_ipynb_inputs.sh` performs this step automatically for the HUES64 WT notebook example.

## 2. WGBS Downloads

Download WGBS files from GEO using the accession pages listed in the paper. GEO file links can be downloaded with the URL pattern:

```bash
curl -L -C - \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=ACCESSION&file=FILENAME&format=file" \
  -o data/FILENAME
```

Known original-repo dataset names used by `scripts/intersect_wgbs_cpg_density_all.py`:

| Sample | Processed WGBS filename expected in `data/` | Output CpG-density filename |
|---|---|---|
| HUES64 WT | `GSM1112841_HUES64WT_WGBS_proc.bed` | `HUES64WT_CpGsOnly_Chr1.bed` |
| HUES64 DNMT3A KO early | `GSM1545002_DNMT3A_KO_Early_proc.bed` | `HUES64_DNMT3Ako_early_CpGsOnly_Chr1.bed` |
| HUES64 DNMT3B KO early | `GSM1545003_DNMT3B_KO_Early_proc.bed` | `HUES64_DNMT3Bko_early_CpGsOnly_Chr1.bed` |
| HUES64 DNMT3A KO late | `GSM1545005_DNMT3A_KO_Late_proc.bed` | `HUES64_DNMT3Ako_late_CpGsOnly_Chr1.bed` |
| HUES64 DNMT3B KO late | `GSM1545006_DNMT3B_KO_Late_proc.bed` | `HUES64_DNMT3Bko_late_CpGsOnly_Chr1.bed` |
| HUES64 DKO early | `GSM1545004_DKO_Early_proc.bed` | `HUES64_DNMT3_dko_early_CpGsOnly_Chr1.bed` |
| HUES64 DKO late | `GSM1545007_DKO_Late_proc.bed` | `HUES64_DNMT3_dko_late_CpGsOnly_Chr1.bed` |
| HUES8 WT | `GSM3618718_HUES8_WT_WGBS_proc.bed` | `HUES8WT_CpGsOnly_Chr1.bed` |
| HUES8 DNMT3AB DKO P6 | `GSM4458672_WGBS_HUES8_DKO_P6_proc.bed` | `HUES8_DKO_P6_CpGsOnly_Chr1.bed` |
| HUES8 TET1-3 TKO | `GSM3618720_HUES8_TKO_WGBS_proc.bed` | `HUES8_TKO_CpGsOnly_Chr1.bed` |
| HUES8 QKO | `GSM3618719_HUES8_QKO_WGBS_proc.bed` | `HUES8_QKO_CpGsOnly_Chr1.bed` |
| HUES8 PKO P0 | `GSM3618721_HUES8_PKO_WGBS_proc.bed` | `HUES8_PKO_P0_CpGsOnly_Chr1.bed` |
| HUES8 PKO P6 | `GSM3662266_HUES8_PKO_P6_WGBS_proc.bed` | `HUES8_PKO_P6_CpGsOnly_Chr1.bed` |
| HUES8 PKO P20 | `GSM4458671_WGBS_HUES8_PKO_P20_proc.bed` | `HUES8_PKO_P20_CpGsOnly_Chr1.bed` |
| IMR90 WT | `GSM432687_IMR90_WGBS_proc.bed` | `IMR90WT_CpGsOnly_Chr1.bed` |

For IMR90 WT, use `GSM432687` to reproduce the paper/repo value. `GSM1204464` is another IMR90 sample and gives a different crossover.

## 3. Normalize WGBS Files

The pipeline expects processed WGBS BED files:

```text
chr  start  end  WGBS
```

`WGBS` must be a methylation fraction from 0 to 1.

Some WGBS exports report the two opposite-strand observations for the same CpG
dinucleotide as adjacent genomic rows. The downstream analysis expects one row
per CpG dyad, so collapse those paired strand rows when they are present. The
helper used below keeps the second row by default, matching the original README
command `awk 'NR % 2 == 0'`.

### WIG Or WIG.GZ Files

Convert WIG to BED with `wig2bed`:

```bash
gzip -cd data/SAMPLE.wig.gz | wig2bed > data/SAMPLE.bed
```

Then normalize and deduplicate adjacent opposite-strand CpG rows:

```bash
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' \
  data/SAMPLE.bed \
  | sort -k1,1 -k2,2n \
  > data/SAMPLE_WGBS_all_strands.bed

python scripts/deduplicate_cpg_strands.py \
  --input data/SAMPLE_WGBS_all_strands.bed \
  --output data/SAMPLE_WGBS_proc.bed \
  --method second
```

If you have verified that the file is already one row per CpG dyad, you can skip
the deduplication step and write directly to `*_WGBS_proc.bed`.

### BED Or BED.GZ Files

If the downloaded BED already has `chr start end methyl_fraction`, sort it:

```bash
gzip -cd data/SAMPLE.bed.gz \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' \
  | sort -k1,1 -k2,2n \
  > data/SAMPLE_WGBS_all_strands.bed

python scripts/deduplicate_cpg_strands.py \
  --input data/SAMPLE_WGBS_all_strands.bed \
  --output data/SAMPLE_WGBS_proc.bed \
  --method second
```

If the BED is known to be one row per CpG dyad, sort directly:

```bash
gzip -cd data/SAMPLE.bed.gz \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' \
  | sort -k1,1 -k2,2n \
  > data/SAMPLE_WGBS_proc.bed
```

If the BED stores methylated and total counts, compute the fraction explicitly,
then deduplicate if paired strand rows are present:

```bash
gzip -cd data/SAMPLE.bed.gz \
  | awk 'BEGIN{OFS="\t"} $5>0 {print $1,$2,$3,$4/$5}' \
  | sort -k1,1 -k2,2n \
  > data/SAMPLE_WGBS_all_strands.bed

python scripts/deduplicate_cpg_strands.py \
  --input data/SAMPLE_WGBS_all_strands.bed \
  --output data/SAMPLE_WGBS_proc.bed \
  --method mean
```

Adjust `$4/$5` to the correct count columns for the downloaded file. Use
`--method mean` when the two strand rows can have different methylation values
and you want a dyad-level average; use `--method second` to reproduce the older
README convention.

You can require near-equal paired values before collapsing:

```bash
python scripts/deduplicate_cpg_strands.py \
  --input data/SAMPLE_WGBS_all_strands.bed \
  --output data/SAMPLE_WGBS_proc.bed \
  --method second \
  --max-methylation-diff 0.001
```

### Bismark Coverage Files

Bismark coverage files have:

```text
chromosome  start_position  end_position  methylation_percentage  count_methylated  count_unmethylated
```

Convert `.cov` or `.cov.gz` to processed WGBS BED:

```bash
python scripts/parse_bismark_cov.py \
  --cov data/SAMPLE.cov.gz \
  --output data/SAMPLE_WGBS_proc.bed
```

By default the methylation fraction is computed from exact counts. Use `--methylation-source percentage` to use the percentage column instead.

If the Bismark coverage file contains adjacent opposite-strand CpG rows, collapse
them during parsing. The default deduplication method for Bismark is
coverage-weighted, which combines methylated/unmethylated counts from the paired
strand rows before calculating the dyad methylation fraction:

```bash
python scripts/parse_bismark_cov.py \
  --cov data/SAMPLE.cov.gz \
  --output data/SAMPLE_WGBS_proc.bed \
  --deduplicate-cpg-strands
```

The coverage-weighted dyad value is:

```text
dyad_WGBS = (meth1 + meth2) / (meth1 + unmeth1 + meth2 + unmeth2)
```

Use `--deduplicate-method mean` only if you explicitly want an unweighted
average of the two strand-level methylation fractions.

## 4. Prepare Notebook And CGI-Level Inputs

The HUES64 WT notebook example can be prepared with:

```bash
bash scripts/prepare_ipynb_inputs.sh
```

This downloads:

- UCSC hg19 `cpgIslandExt.txt.gz`
- GEO `GSM1112841_BI.HUES64.Bisulfite-Seq.WGBS_Lib_39.wig.gz`

It writes:

- `data/CpGIsAnn_hg19_dec14.bed`
- `data/CpGIsAnn_hg19_chr1.bed`
- `data/GSM1112841_HUES64WT_WGBS_proc.bed`
- `data/HUES64WT_WGBS_CGI_int.bed`

Validate the notebook-compatible files:

```bash
python scripts/check_ipynb_inputs.py
```

For any processed WGBS file, create the CGI intersection with:

```bash
bedtools intersect \
  -a data/SAMPLE_WGBS_proc.bed \
  -b data/CpGIsAnn_hg19_dec14.bed \
  -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$8,$9,$10}' \
  > data/SAMPLE_WGBS_CGI_int.bed
```

Then aggregate to island level:

```bash
python scripts/aggregate_cgi_level.py \
  --cgi-intersection data/SAMPLE_WGBS_CGI_int.bed \
  --output data/IslandLvl_agg_SAMPLE
```

## 5. Prepare CpG Coordinates For Individual CpG Analysis

Download hg19 chromosome FASTA files from UCSC. For chr1 only:

```bash
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz data/
```

Find CpG positions with `seqkit`:

```bash
gzip -cd data/chr1.fa.gz \
  | seqkit locate -P -p cg \
  > data/chr1_cpgs.csv
```

The density script expects a CSV with a `start` column, which is the default `seqkit locate` output.

## 6. Prepare CpG-Density Intersections

Calculate local CpG density:

```bash
python scripts/calculate_cpg_density.py \
  --cpg-csv data/chr1_cpgs.csv \
  --window 50 \
  --chrom chr1 \
  --output data/CpGDensities_W50.bed
```

Intersect one processed WGBS file with CpG density:

```bash
python scripts/intersect_wgbs_cpg_density.py \
  --wgbs-bed data/SAMPLE_WGBS_proc.bed \
  --density-bed data/CpGDensities_W50.bed \
  --output data/SAMPLE_CpGsOnly_Chr1.bed
```

Or use `bedtools`:

```bash
bedtools intersect \
  -a data/SAMPLE_WGBS_proc.bed \
  -b data/CpGDensities_W50.bed \
  -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$9}' \
  > data/SAMPLE_CpGsOnly_Chr1.bed
```

Run all original-repo names with:

```bash
bash scripts/intersect_wgbs_cpg_density_all_bedtools.sh \
  data \
  data/bedtools_density_intersections \
  data/CpGDensities_W50.bed
```

## 7. Known Prepared Results

The following paper-matching WT checks were generated from the prepared local files:

- HUES64 WT: 351 bp CGI crossover.
- HUES8 WT: 321 bp CGI crossover.
- IMR90 WT, `GSM432687`: 214 bp CGI crossover.
- IMR90 WT, alternate `GSM1204464`: 306 bp CGI crossover.

The combined bivariate check outputs are written under `data/figure_checks/` when using the commands in `PYTHON_PIPELINE.md`.
