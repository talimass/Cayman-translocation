#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=8 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=8 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=64G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=%J.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=%J.output ## stdout file name(The %J will print job ID number)
#SBATCH --mail-type=NONE ## Send your job status via e-mail: Valid type values are NONE, BEGIN, END, FAIL, REQUEUE, ALL
########### Job information #############
echo "================================"
echo "Start at `date`"
echo "Job id is $SLURM_JOBID"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NTASKS processors."
echo "================================"
#########################################

######## Load required modules ##########
#. /etc/profile.d/modules.sh # Required line for modules environment to work
#module load openmpi/1.8.4 python/2.7 # Load modules that are required by your program
#########################################
source /lustre1/home/mass/eskalon/miniconda/bin/activate gatk2
### Below you can enter your program job command ###


REF="./Pastreoides.genome.fasta"
OUT="GVCFall"

SNP_QD_MIN=20.00
SNP_MQ_MIN=50.00
SNP_FS_MAX=200.00
SNP_SOR_MAX=4.00

INDEL_QD_MIN=20.00
INDEL_MQ_MIN=45.00
INDEL_FS_MAX=60.00
INDEL_SOR_MAX=4.00

QUAL=30.00


# Extract variant quality scores and make diagnostic plots 
gatk java-options "-Xmx64G -XX:ParallelGCThreads=8" SelectVariants \
     -R ${REF} -V  cohort_genotypes.vcf.gz -O ${OUT}.SNPs.vcf.gz \
      -select-type SNP 1> ${OUT}_SNPs.vcf.gz.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" SelectVariants \
     -R ${REF} -V  cohort_genotypes.vcf.gz -O ${OUT}.INDELs.vcf.gz \
     --select-type-to-include INDEL 1> ${OUT}_INDELs.vcf.gz.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantsToTable\
    --variant ${OUT}.SNPs.vcf.gz \
    --output ${OUT}_SNPs.table \
    -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum \
    -F ReadPosRankSum 1> ${OUT}_SNPs.table.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantsToTable \
     --variant ${OUT}.INDELs.vcf.gz \
     --output ${OUT}_INDELs.table -F CHROM -F POS -F QUAL -F QD -F DP \
     -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
     1> ${OUT}_INDELs.table.log 2>&1

# 1st-pass filtering - values are chosed according to the diagnostic plots
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" \
        VariantFiltration --reference ${REF} --variant ${OUT}.SNPs.vcf.gz  \
        --output ${OUT}_SNPs_VarScores_filter.qd.vcf.gz \
        -filter "QUAL < $QUAL"  --filter-name "VarScores_filter_QUAL"  \
        -filter  "QD < $SNP_QD_MIN"  --filter-name "VarScores_filter_QD" \
        -filter  "MQ < $SNP_MQ_MIN" --filter-name "VarScores_filter_MQ"  \
        -filter "FS > $SNP_FS_MAX"  --filter-name "VarScores_filter_FS"  \
        -filter  "SOR > $SNP_SOR_MAX" --filter-name "VarScores_filter_SOR" \
        1> ${OUT}_SNPs_VarScores_filter.qd.vcf.gz.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" \
        VariantFiltration --reference ${REF} --variant ${OUT}.INDELs.vcf.gz \
        --output ${OUT}_INDELs_VarScores_filter.qd.vcf.gz \
        -filter "QUAL < $QUAL"  --filter-name "VarScores_filter_QUAL"  \
	-filter  "QD < $SNP_QD_MIN"  --filter-name "VarScores_filter_QD" \
	-filter  "MQ < $SNP_MQ_MIN" --filter-name "VarScores_filter_MQ"  \
	-filter "FS > $SNP_FS_MAX"  --filter-name "VarScores_filter_FS"  \
	-filter  "SOR > $SNP_SOR_MAX" --filter-name "VarScores_filter_SOR" \
 	1> ${OUT}_INDELs_VarScores_filter.qd.vcf.gz.log 2>&1

# Check the number of PASSED after the first filtering
zcat "${OUT}.SNPs.vcf.gz" | grep -v '^#' | wc -l
zcat "${OUT}_SNPs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l
zcat "${OUT}.INDELs.vcf.gz" | grep -v '^#' | wc -l
zcat "${OUT}_INDELs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l


# Extract only variants that PASSED filtering
zcat "${OUT}_SNPs_VarScores_filter.qd.vcf.gz" | grep -E '^#|PASS' > "${OUT}_SNPs_VarScores_filterPASSED.vcf"
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" IndexFeatureFile --input "${OUT}_SNPs_VarScores_filterPASSED.vcf" 1> "${OUT}_SNPs_VarScores_filterPASSED.vcf.log" 2>&1
zcat "${OUT}_INDELs_VarScores_filter.qd.vcf.gz" | grep -E '^#|PASS' > "${OUT}_INDELs_VarScores_filterPASSED.vcf"
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" IndexFeatureFile --input "${OUT}_INDELs_VarScores_filterPASSED.vcf" 1> "${OUT}_INDELs_VarScores_filterPASSED.vcf.log" 2>&1

# Extract and plot DP info for each sample (from all samples before previous filtering; shows us overall variant coverage, independent of other factors)
gatk VariantsToTable --reference "${REF}" --variant cohort_genotypes.vcf.gz --output "${OUT}.DP.table" -F CHROM -F POS -GF GT -GF DP 1> "${OUT}.DP.table.log" 2>&1

# Each sample has 2 columns (GT, DP); get the GT column index for each sample first, extract idx and idx+1, filter using these two columns.
while read i;
do
        GT=$(cut -f $i "${OUT}.DP.table" | head -n 1)
        cut -f $i,$((i+1)) "${OUT}.DP.table" | awk '$1 != "./." && $1 != ".|." {print $2}' > $GT.DP.txt
done < <(head -n 1 "${OUT}.DP.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')

# Make diagnostic plots and choose the apropriate DP value
DP_MIN=20.00

# 2nd-pass filtering
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantFiltration \
 --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED.vcf" \
 --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf" \
	        --genotype-filter-name "DP_filter" \
 --genotype-filter-expression "DP < $DP_MIN" \
	        1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantFiltration \
 --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED.vcf" \
 --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf" \
	        --genotype-filter-name "DP_filter" \
 --genotype-filter-expression "DP < $DP_MIN" \
	        1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1

# Set filtered sites to no call
gatk SelectVariants --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
	1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1
gatk SelectVariants --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
	1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1


## Check number of variants that PASSED after filtering
#for F in *.GT.DP.txt; do awk 'NR>1' $F; done | wc -l
# 430694419
awk 'NR>1' "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 69361
awk 'NR>1' "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 4995


# Extract and plot DP info for each sample (combine it all together and plot; we dont care about each sample as we just want to check that filtering worked)
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table" -F CHROM -F POS -GF GT -GF DP \
	1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table.log" 2>&1
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table" -F CHROM -F POS -GF GT -GF DP \
	1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table.log" 2>&1

while read i;
do
        cut -f $i,$((i+1)) "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk '$1 != "./." && $1 != ".|." {print $2}' > "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt"
done < <(head -n 1 "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')
while read i;
do
        cut -f $i,$((i+1)) "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk '$1 != "./." && $1 != ".|." {print $2}' > "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt"
done < <(head -n 1 "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')


gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" -F CHROM -F POS -GF GT -GF AD -GF DP \
	1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table.log" 2>&1

# Each sample has 3 columns (GT, AD, DP); get the GT column index for each sample first, extract idx - idx+2, filter using these two columns.
while read i;
do
        GT=$(cut -f $i "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | head -n 1)
        cut -f $i-$((i+2)) "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '$1 != "./." && $1 != ".|." { if(NR==1) {print $2} else {split($2,a,","); for (i in a) {print a[i]/$3} } }' > $GT.AD.txt
done < <(head -n 1 "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')


# Generate allelic depth plots using sum(AD) as denominator instead of DP
while read i;
do
        GT=$(cut -f $i "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | head -n 1)
        cut -f $i-$((i+2)) "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '$1 != "./." && $1 != ".|." { if(NR==1) {print $2} else {split($2,a,","); SUM=0; for (i in a) {SUM=SUM+a[i]}; if(SUM>0){for (i in a) {print a[i]/SUM}} } }' > $GT.AD.sumSDdenominator.txt
done < <(head -n 1 "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')


# VCF to Table (use for phylogeny)
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf"   --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.table"   -F CHROM -F POS -GF GT \
	1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.table.log" 2>&1
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.table" -F CHROM -F POS -GF GT \
	1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.table.log" 2>&1
