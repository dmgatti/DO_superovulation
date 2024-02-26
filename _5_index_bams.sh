#!/bin/bash
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 24G # memory pool for all cores
#SBATCH -t 0-4:00:00 # time (D-HH:MM)

################################################################################
# Sort and index genome BAM files created by STAR & RSEM using samtools.
# https://deweylab.github.io/RSEM/
# http://www.htslib.org/
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-09-17
################################################################################


### VARIABLES ###

# Base directory for this user.
BASE_DIR=/fastscratch/education/dgatti

# Full path to directory where genome BAM files are stored.
BAM_DIR=${BASE_DIR}/rnaseq/results/rsem

# Temporary directory where intermediate files will be stored.
TEMP_DIR=${BASE_DIR}/temp

# Get a list of all of the BAM files in the BAM_DIR.
BAM_FILES=`ls ${BAM_DIR}/*genome.bam`

# Full path to the Samtools container.
CONTAINER=/projects/education/containers/samtools.sif


### PROGRAM ###

# Index each BAM file.
for F in ${BAM_FILES}
do

   echo Sorting ${F}
   singularity exec ${CONTAINER} samtools sort -@ 8 \
                                               -T ${TEMP_DIR} \
                                               -o ${F}_sorted.bam \
                                               ${F}

   echo Indexing ${F}
   singularity exec ${CONTAINER} samtools index ${F}_sorted.bam

   # Remove the unsorted BAM file.
   rm ${F}

done

