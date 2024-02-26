#!/bin/bash
#SBATCH --qos=dev
#SBATCH --partition=dev
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=8G # memory pool for all cores
#SBATCH --time=2:00:00 # time (D-HH:MM)

##### VARIABLES #####

GENOME=grcm38

ENSEMBL=102

R=~/containers/r_qtl2.sif

##### MAIN #####

module load singularity

singularity exec ${R} Rscript local_eqtl.R ${GENOME} ${ENSEMBL}



