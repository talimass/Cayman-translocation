## Fixation index (Fst) estimation by RNA-seq SNPs analysis of _P. astreoides_ samples 

This script is based on the pipeline for SNP analysis of [Coral_Stress_Phenome](https://github.com/hputnam/Coral_Stress_Phenome/tree/main/Genotype_Analysis/Pocillopora_acuta_PacBio_Assembly/RNAseq_short_variant_analysis) project and on the work of [Federica Scuccia](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs), with some modifications.

### Programs
GATK4 pipeline was installed on the remote cluster through conda using this [guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4). PLINK2 and bcftools were also installed via conda. rgsam was installed through [github](https://github.com/djhshih/rgsam), and then the Makefile was manually changed to ensure installation to the local home bin folder.

### Prepare BAM files: Convert reads to BAM, add read group info, merge alignments, filter 

[Fastqtobam.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/fastqtobam.sh) script converts paired fastq files to BAM files sorted by read name, [rgsam.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/rgsam.sh) adds read group info (RG) to aligned reads, [dictreference.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/dictreference.sh) indexes and creates a dictionary from the reference, [mergebam.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/mergebam.sh) merges unaligned BAM files with aligned BAMs and additionally filters the alinged read (e.g. removes secondary alignments). 

### Mark duplicates
[markdup.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/markdup.sh) marks potential PCR duplicates.

### Split reads that contain Ns in their CIGAR string 
[splitsigar.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/splitsigar.sh) splits reads with N in CIGAR - used for post-processing of the alignment.

### Call SNPs and indels
[haplotcall.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/haplotcall.sh) calls variants (SNPs and indels) simultaneously via local de-novo assembly of haplotypes in an active region. 

### Combine *.g.vcf.gz files and call genotypes
[combinegvcf.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/combinegvcf.sh) combines all gvcf files to a cohort file, [genotypegvcf.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/genotypegvcf.sh) performs joint genotyping on a cohort. 

### Filter SNPs and indels
[filtvar.sh]() firstly selects SNPs and indels and extracts variant quality scores for making diagnostics plots with [diagnostics.R](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/diagnostics.R), then performs 1st-pass filtering of variants with the parameters chosen according to the diagnostics plots, then performs 2nd and 3rd-pass filtering with the parameters chosen according to the diagnostics plots drawn by [plot.dp.scores.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/plot.dp.scores.sh) script.   

Number of variants passed:
```
# Check the number of PASSED after the first filtering
zcat "${OUT}.SNPs.vcf.gz" | grep -v '^#' | wc -l
# 4142116
zcat "${OUT}_SNPs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l
# 2791170
zcat "${OUT}.INDELs.vcf.gz" | grep -v '^#' | wc -l
# 838000
zcat "${OUT}_INDELs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l
# 622066
```
```
## Check number of variants that PASSED after the 2nd-pass filtering
#for F in *.GT.DP.txt; do awk 'NR>1' $F; done | wc -l
# 430694419
awk 'NR>1' "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 69361
awk 'NR>1' "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 4995
```
### Filter for linkage disequilibrium
[bcfplink.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/bcfplink.sh) performs pruning of SNPs that are strongly genetically linked and extracts these SNPs from the vcf file.

### Fst analysis
[R_run.sh](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/R_run.sh) runs the script [R_Fst.R](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/R_Fst.R) to install necessary libraries and then [Fst.R](https://github.com/talimass/Cayman-translocation/blob/main/RNAseq_analysis/Connectivity/Fst.R) to run the Fst analysis between both different depths and sites of the samples.

```
# The total number of samples: 31 
# The total number of SNPs: 48473
```
