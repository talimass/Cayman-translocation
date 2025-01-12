#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=8 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=8 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=32G ## memory to allocate for your job in MB
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
for R1 in  /lustre1/home/mass/eskalon/Porites/raw_data/30-1082101046/00_fastq/*_R1_001.fastq.gz; do
R2="${R1%_R1_001.fastq.gz}_R2_001.fastq.gz"
R3="${R1%_R1_001.fastq.gz}_trim.gz"

trimmomatic PE -threads 8 "$R1" "$R2" -baseout "$R3" \
ILLUMINACLIP:/lustre1/home/mass/eskalon/Porites/analysis/trimming/Sequencing_adaptors.fasta:2:30:10 \
SLIDINGWINDOW:4:5 MAXINFO:50:0.6 MINLEN:25


done

fastqc --noextract -o /lustre1/home/mass/eskalon/Porites/analysis/trimming/fastqc-after2 \
 -t 8 /lustre1/home/mass/eskalon/Porites/raw_data/30-1082101046/00_fastq/*P*.gz
