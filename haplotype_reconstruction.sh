#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 30 # number of cores
#SBATCH --mem 180G # memory pool for all cores
#SBATCH -t 0-4:00 # time (D-HH:MM)

GENOME=grcm39

CONTAINER=~/containers/r_qtl2.sif

RSCRIPT=/projects/bolcun-filas-lab/DO_Superovulation/scripts/DO_superovulation/haplotype_reconstruction.R

module load singularity

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${GENOME}


