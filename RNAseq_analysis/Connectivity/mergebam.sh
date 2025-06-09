#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=16 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=16 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=128G ## memory to allocate for your job in MB
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

TAGs=( \
123-CC-DD-81 \
117-MF-SD-44 \
114-MF-SD-45 \
108-CC-DD-94 \
102-CC-DD-82 \
100-CC-SD-73 \
99-MF-SD-43 \
98-MF-SS-69 \
97-MF-SS-64 \
94-CC-SS-37 \
91-CC-DS-89 \
84-CC-DD-80 \
81-CC-DS-87 \
80-1-CC-DS-86 \
78-MF-DS-54 \
75-CC-DD-84 \
73-CC-SS-38 \
66-MF-DS-52 \
63-MF-SS-60 \
49-CC40-1b \
48-CC10-6b \
42-MF10-10b \
37-CC40-8b \
35-MF40-1b \
32-CC10-1b \
31-MF40-5b \
28-CC10-9b \
27-MF40-2b \
25-CC40-2b \
23-MF10-1a \
4-MF10-4a \
)

REF="./Pastreoides.genome.fasta"

for R1 in ${TAGs[@]}; do
 BASENAME=${R1}
 OUT="./$BASENAME"
  
 echo "aligned ${R1}, unaligned ${OUT}.FastqToSam.unmapped.rg.bam"
 
 gatk --java-options "-Xmx128G -XX:ParallelGCThreads=16" MergeBamAlignment \
                	--REFERENCE_SEQUENCE ${REF} \
                	--UNMAPPED_BAM ${OUT}.FastqToSam.unmapped.rg.bam \
                	--ALIGNED_BAM \
 /lustre1/home/mass/eskalon/Porites/analysis/mapping/bams2/${R1}_Aligned.sortedByCoord.out.bam \
                	--OUTPUT ${OUT}.MergeBamAlignment.merged.bam \
                	--INCLUDE_SECONDARY_ALIGNMENTS false \
                	--VALIDATION_STRINGENCY SILENT --TMP_DIR ./temp
done


