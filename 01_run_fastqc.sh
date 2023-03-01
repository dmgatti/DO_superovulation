#!/bin/bash
#SBATCH --qos=dev
#SBATCH --partition=dev
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=30 # number of cores
#SBATCH --mem=32G # memory pool for all cores
#SBATCH --time=08:00:00 # time (D-HH:MM)

################################################################################
# Perform quality control on FASTQ files using FastQC.
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# 
# Daniel Gatti
# 2022-09-28
################################################################################


### VARIABLES ###

# Top level directory for project.
BASE_DIR=/projects/bolcun-filas-lab/DO_Superovulation

# Full path to the FASTQ file directories.
FASTQ_DIR1=${BASE_DIR}/data/rnaseq/fastq/run1
FASTQ_DIR2=${BASE_DIR}/data/rnaseq/fastq/run2
FASTQ_DIR3=${BASE_DIR}/data/rnaseq/fastq/run3

# Full path to the FastQC output directory.
OUTPUT_DIR=${BASE_DIR}/results/fastqc

# Full path to the FastQC and MultiQC containers.
FASTQC=~/containers/fastqc.sif
MULTIQC=~/containers/multiqc.sif

### PROGRAM ###

module load singularity

mkdir -p ${BASE_DIR}/results/fastqc

# Run FastQC on all fastq files.
echo "Processing run1"
singularity exec ${FASTQC} fastqc \
                              --threads 30 \
                              --outdir ${OUTPUT_DIR} \
                              ${FASTQ_DIR1}/*.gz

echo "Processing run2"
singularity exec ${FASTQC} fastqc \
                              --threads 30 \
                              --outdir ${OUTPUT_DIR} \
                              ${FASTQ_DIR2}/*.gz

echo "Processing run3"
singularity exec ${FASTQC} fastqc \
                              --threads 30 \
                              --outdir ${OUTPUT_DIR} \
                              ${FASTQ_DIR3}/*.gz
                              

# Run MultiQC to summarize FastQC results.
singularity exec ${MULTIQC} multiqc ${OUTPUT_DIR}

