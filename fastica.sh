#!/bin/bash
#SBATCH --qos=dev
#SBATCH --partition=dev
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=30 # number of cores
#SBATCH --mem=180G # memory pool for all cores
#SBATCH --time=08:00:00 # time (D-HH:MM)

### VARIABLES ###

# Top level directory for project.
BASE_DIR=/projects/bolcun-filas-lab/DO_Superovulation/scripts/DO_superovulation

### PROGRAM ###

cd ${BASE_DIR}

module load singularity

# Run FastQC on all fastq files.
singularity exec ~/containers/fastica.sif Rscript fastica.R

