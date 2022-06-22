#!/bin/bash
#SBATCH -J brm_reg_sympt
#SBATCH -n 6               # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-24:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=24000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o jobmessages/job%j-%a.out
#SBATCH -e jobmessages/jobERR%j-%a.out
#SBATCH --array=1-6
echo $SLURM_ARRAY_TASK_ID

mkdir -p jobout/${SLURM_JOB_NAME}
mkdir -p jobmessages/${SLURM_JOB_NAME}

module load gcc/9.3.0-fasrc01 #Load gcc
module load R/4.0.2-fasrc01 #Load R module

export R_LIBS_USER=$HOME/apps/R_4.0.2

R CMD BATCH --quiet --no-restore --no-save cluster/brm_regression_cluster_sympt.R jobout/${SLURM_JOB_NAME}/run_${SLURM_ARRAY_TASK_ID}.Rout
