#! /bin/bash
CpGFilename=CpGDensities_W50
bedtools intersect -a GSM1112841_HUES64WT_WGBS_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES64WT_CpGsOnly_Chr1.bed
bedtools intersect -a GSM1545002_DNMT3A_KO_Early_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES64_DNMT3Ako_early_CpGsOnly_Chr1.bed
bedtools intersect -a GSM1545003_DNMT3B_KO_Early_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES64_DNMT3Bko_early_CpGsOnly_Chr1.bed
bedtools intersect -a GSM1545005_DNMT3A_KO_Late_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES64_DNMT3Ako_late_CpGsOnly_Chr1.bed
bedtools intersect -a GSM1545006_DNMT3B_KO_Late_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES64_DNMT3Bko_late_CpGsOnly_Chr1.bed
bedtools intersect -a GSM1545004_DKO_Early_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES64_DNMT3_dko_early_CpGsOnly_Chr1.bed
bedtools intersect -a GSM1545007_DKO_Late_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES64_DNMT3_dko_late_CpGsOnly_Chr1.bed
bedtools intersect -a GSM3618718_HUES8_WT_WGBS_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES8WT_CpGsOnly_Chr1.bed
bedtools intersect -a GSM4458672_WGBS_HUES8_DKO_P6_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES8_DKO_P6_CpGsOnly_Chr1.bed
bedtools intersect -a GSM3618720_HUES8_TKO_WGBS_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES8_TKO_CpGsOnly_Chr1.bed
bedtools intersect -a GSM3618719_HUES8_QKO_WGBS_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES8_QKO_CpGsOnly_Chr1.bed
bedtools intersect -a GSM3618721_HUES8_PKO_WGBS_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES8_PKO_P0_CpGsOnly_Chr1.bed
bedtools intersect -a GSM3662266_HUES8_PKO_P6_WGBS_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES8_PKO_P6_CpGsOnly_Chr1.bed
bedtools intersect -a GSM4458671_WGBS_HUES8_PKO_P20_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > HUES8_PKO_P20_CpGsOnly_Chr1.bed
bedtools intersect -a GSM432687_IMR90_WGBS_proc.bed -b  $CpGFilename.bed -wa -wb -sorted | awk ' {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > IMR90WT_CpGsOnly_Chr1.bed




			

			

				

			


