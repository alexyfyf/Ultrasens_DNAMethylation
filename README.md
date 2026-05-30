# Ultrasens_DNAMethylation
This repository is a collection of code for the associated publication, "A Tunable, Ultrasensitive Threshold in Enzymatic Activity governs the DNA Methylation Landscape". It contains various terminal/console commands, Python scripts written in Jupyter Notebook, and MATLAB files which manipulate the required genomic datasets, perform parameterization of ensuing models to experimental data, and create the visualizations shown in the publication.

In directory `paper_matlab_reference/CMEFitting/`
Contains scripts for fitting the density-dependent methylation curves extracted from WGBS data using a CME model. The data-derived curves from HUES8 are stored in Save_HUES8WT_CpGsOnly_Chr1.mat. Run LoopFit.m to call the function which loads the data and performs the fit.

In directory `paper_matlab_reference/RunTemporalCMEExample/`
Contains scripts to compute the methylation versus density curves over time, according to the CME model, from a specified initial condition. Run Plot_ExampleTemporal.m.

## Python rewrite

The MATLAB analysis workflow has been rewritten as Python command-line scripts under `scripts/`, with reusable functions in `python/ultrasens/`.
The original paper MATLAB/notebook reference files are retained unchanged under `paper_matlab_reference/`.

- See `PYTHON_PIPELINE.md` for the complete step-by-step Python workflow, following the original order: CpG island analysis and plots first, then individual CpG density analysis, switch plots, and CME fitting.
- See `DATA_PREP.md` for data downloads, WGBS normalization, CGI annotation preparation, Bismark coverage parsing, and CpG-density input preparation.

## Goal
Taking .bed file of WGBS (or other methylation) dataset(s), perform CGI-Level analysis of methylation macrostates, Individual CpG analysis of methyl fraction dependence on local CpG landscape, perform Hill function fitting of Individual CpG mean methylation, and finally fit landscape using Master Equation.

## Files Needed
1. CpG Island Annotation (.bed File)
	1. Available directly from the UCSC Genome/Table Browser under "Regulation". Note that this is specifically created for the hg19 Assembly, but [This Page](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=cpgIslandExt) describes the schema, along with how to manually perform the island annotation.
    
1. Whole-Genome Bisulfite Sequencing (WGBS) Datasets (.bed file(s) )
   1. Via their, GEO Accession numbers (shown in Table 1 within the Main Text), these files can be downloaded directly from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
   2. Note that some files may be in the ".wig" format, but can be converted within the terminal environment using the following command (available via the bedops package):
      
```
wig2bed < $WGBS_filename.wig > $WGBS_filename.bed
```
      
1. Entire Genome Sequence (hg19) (.fasta file)
   1. This may be downloaded also from the UCSC Genome Browser using the following command (Mac OS)

```
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ .
```
   2. Note that each chromosome is separated into its own file; the entire genome requires ~3.2gb of space.
   3. Using the "Seqkit" package, one can filter each chromosome for CpGs using the following command from the terminal environment; the example here performs the operation only on Chromosome 1:
```
cat chr1.fa | seqkit locate -P -p  cg > {OUTPUT_FILENAME}.csv
```

## Included Scripts
### CGI-Level Analysis
1. paper_matlab_reference/Analysis of Human ESCS (HUES64WT_example)
   1. Written for Python/Jupyter Notebook (.ipynb)
   2. Provides visualization of CpG Number vs. CGI Length, along with outputting required "IslandLvl_Agg_###", which are used in the ensuing MATLAB scripts
   3. The WGBS (.bed) files must be sorted and intersected with the genome assembly (using the bedtools package). This can be accomplished using the following code block from the terminal environment:
```
# Step 1: Sorting WGBS file
f=$WGBS_filename.bed
g=$WGBS_filename”_proc”.bed
touch $g 
out_fp=$WGBS_filename”.tmp”
awk 'NR % 2 == 0' $f | awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"$5}' | sort -k1,1 -k2,2n > $out_fp
mv $out_fp $g

#Step 2: Intersect WGBS and CGI Annotation Bed files
##Generic
bedtools intersect -a WGBS_BED_PATH -b  CGI_PATH -wa -wb -sorted | awk ‘{print $1”\t”$2”\t”$3”\t”$4”\t”$8"\t"$9"\t"$10}’ > TARGET_OUTPUT_PATH.bed
```

2. paper_matlab_reference/BivariateHistogram_HumanWT_example
   1. Written for MATLAB  2023b (.m)
   2. Provides bivariate histograms and CpG island classifications shown in the paper figures.
   3. Needs the "IslandLvl_Agg_###" files for HUES64, HUES8, and IMR90 WT cell lines for the ensuing example.

3. paper_matlab_reference/IndividualCGIMethylationChange_HUES64_example
   1. Written for MATLAB 2023b (.m)
   2. Calculates percent change of each CGI between the WT and its ensuing Knockout.
   3. Needs the "IslandLvl_Agg_###" files for HUES64WT, along with DKO (Early/Late) for the ensuing example.

### Individual CpG-Level Analysis
1. paper_matlab_reference/CpGDensity_Calc
   1. Written for MATLAB 2023 (.m)
   2. Using the .csv containing a list of CpG locations, this script calculates the local density of each CpG listed. Generates "CpGDensities_W##" which is used to intersect with the raw WGBS datasets.

2. paper_matlab_reference/WGBS_CpGIntersect_AllData
   1. Written to be run from the terminal environment; requires the sorted/processed file of WGBS data, the bedtools package and CpGDensities_W##.csv output file.
   2. In an coordinated manner, this .command file intersects the WGBS information with CpG density calculation (given a certain window size ##), then outputs .bed files with each CpG containing a value for its local density and associated methyl fraction.

3. paper_matlab_reference/ReadPlotData_WT_example
   1. Written for MATLAB2023 (.m)
   2. Peforms classification into hypo/hypermethylated, or intermediate states, calculates the Mean/Median and 25th/75th percentiles with respect to local CpG density, then finally performs the Simple/Log-transformed Hill function fittting to the mean of each dataset (see Figure ###/Supplemental Figure ###)
   3. Needs the file containing the WGBS/Local CpG Density intersection files for HUES64, HUES8, and IMR90 WT cell lines for the ensuing example.
  
4. paper_matlab_reference/ExampleCMEFit (Folder)
   1. Contains MATLAB scripts, initial parameters, example data (of HUES8WT, in a .mat file) which performs fitting on the derived experimental landscape, and returns a fitted parameter set.
   2. To fit ensuing data, run "CallFitting.m", and change "DataName" to file containing data-derived curve for hyper-/hypomethylated fractions as a function of CpG Density.
