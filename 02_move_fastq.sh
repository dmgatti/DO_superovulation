#!/bin/bash
#SBATCH --qos=batch
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=4G # memory pool for all cores
#SBATCH --time=04:00:00 # time (D-HH:MM)

################################################################################
# copy FASTQ files over to /fastscratch
#
# Daniel Gatti
# 2022-09-28
################################################################################


### VARIABLES ###

# Source directory where original FASTQ files are stored.
SRC_DIR=/projects/bolcun-filas-lab/DO_Superovulation/data/rnaseq/fastq

# Destination directory on /fastscratch.
DEST_DIR=/fastscratch/dgatti/fastq


### MAIN ###

mkdir -p ${DEST_DIR}

cp ${SRC_DIR}/*.gz ${DEST_DIR}

