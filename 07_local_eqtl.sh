#!/bin/bash
#SBATCH --qos=batch
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=8G # memory pool for all cores
#SBATCH --time=3:00:00 # time (D-HH:MM)
#SBATCH --array=1-2

# --array=1-220

# 21929 genes

##### VARIABLES #####

CHUNK_SIZE=100

START=$(( ${SLURM_ARRAY_TASK_ID} * ${CHUNK_SIZE} - ${CHUNK_SIZE} + 1 ))

END=$(( ${SLURM_ARRAY_TASK_ID} * ${CHUNK_SIZE} ))

R=~/containers/r_qtl2.sif

##### MAIN #####

module load singularity

singularity exec ${R} Rscript 07_local_eqtl.R ${START} ${END}


