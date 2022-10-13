#!/bin/bash
#SBATCH --qos=dev
#SBATCH --partition dev
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=32G # memory pool for all cores
#SBATCH --time=1:00:00 # time (D-HH:MM)

##### VARIABLES #####

# R container
R=~/containers/bioconductor.sif

##### MAIN #####

module load singularity

singularity run ${R} Rscript --no-save 04_gather_gbrs_data.R

