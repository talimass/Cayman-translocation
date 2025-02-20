#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=16 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=16 ## number of cpu's to allocate per each node
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

TAGs=(\
27-MF40-2b \
31-MF40-5b \
35-MF40-1b \
37-CC40-8b \
49-CC40-1b \
25-CC40-2b \
102-CC-DD-82 \
108-CC-DD-94 \
123-CC-DD-81 \
75-CC-DD-84 \
84-CC-DD-80 \
66-MF-DS-52 \
78-MF-DS-54 \
80-1-CC-DS-86 \
81-CC-DS-87 \
91-CC-DS-89 \
23-MF10-1a \
28-CC10-9b \
32-CC10-1b \
4-MF10-4a \
42-MF10-10b \
48-CC10-6b \
100-CC-SD-73 \
114-MF-SD-45 \
117-MF-SD-44 \
99-MF-SD-43 \
63-MF-SS-60 \
73-CC-SS-38 \
94-CC-SS-37 \
97-MF-SS-64 \
98-MF-SS-69 \
)

for i in ${TAGs[@]}; do

rm -rf ./sortmerna/run/kvdb/

echo "/lustre1/home/mass/eskalon/Porites/analysis/trimming/trimmed-data/${i}_trim_1P.gz"

sortmerna --ref ./sortmerna/rRNA_databases_v4/smr_v4.3_default_db.fasta \
  --reads  /lustre1/home/mass/eskalon/Porites/analysis/trimming/trimmed-data/${i}_trim_1P.gz \
  --reads /lustre1/home/mass/eskalon/Porites/analysis/trimming/trimmed-data/${i}_trim_2P.gz \
  --workdir ./sortmerna/run \
  --aligned "${i}_rRNA" --other "${i}_clean" \
  --threads 16 --paired_in --fastx

done
