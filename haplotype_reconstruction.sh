#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 30 # number of cores
#SBATCH --mem 180G # memory pool for all cores
#SBATCH -t 0-4:00 # time (D-HH:MM)

cd /projects/bolcun-filas-lab/DO_Superovulation/scripts/DO_superovulation

module load singularity

singularity exec ~/containers/r_qtl2.sif Rscript haplotype_reconstruction.R grcm39


