#!/bin/bash
#SBATCH --qos=batch
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=4G # memory pool for all cores
#SBATCH --time=04:00:00 # time (D-HH:MM)

################################################################################
# Copy DO FASTQ files over to /fastscratch
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-09-28
################################################################################


### VARIABLES ###

# Source directory where original FASTQ files are stored.
SRC_DIR1=/projects/bolcun-filas-lab/DO_Superovulation/data/rnaseq/fastq/run1
SRC_DIR2=/projects/bolcun-filas-lab/DO_Superovulation/data/rnaseq/fastq/run2

# Destination directory on /fastscratch.
DEST_DIR=/fastscratch/dgatti/fastq


### MAIN ###

mkdir -p ${DEST_DIR}

cp ${SRC_DIR1}/SODO*.fastq.gz ${DEST_DIR}

cp ${SRC_DIR2}/SODO*.fastq.gz ${DEST_DIR}


