#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 64G
#SBATCH --time 0-1:00

##### VARIABLES #####

# Genome build: either grcm38 or grcm39.
GENOME=grcm39

# Output directory for qtl2 input files.
OUTDIR=/fastscratch/dgatti/qtl2

# Singularity container.
CONTAINER=/projects/compsci/vmp/USERS/dgatti/containers/bioconductor.sif

##### MAIN #####

module load singularity

singularity exec ${CONTAINER} Rscript ${GENOME} ${OUTDIR}
