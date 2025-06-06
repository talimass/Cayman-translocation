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
source /lustre1/home/mass/eskalon/miniconda/bin/activate rnaseq
#########################################

### Below you can enter your program job command ###

for R1 in /lustre1/home/mass/eskalon/Porites/analysis/trimming/trimmed-data/*_trim_1P.gz; do
R2="${R1%_trim_1P.gz}_trim_2P.gz"
R3="${R1%_trim_1P.gz}_"

STAR --genomeDir index_Pastr.gtf \
 --readFilesIn "$R1" "$R2" \
 --readFilesCommand gunzip -c \
 --runThreadN 8 --outSAMtype BAM SortedByCoordinate \
 --outReadsUnmapped Fastx \
 --outFileNamePrefix "$R3" \
 --limitBAMsortRAM 64424509440  --outBAMsortingBinsN 100
done
