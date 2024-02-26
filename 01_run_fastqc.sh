#!/bin/bash
#SBATCH --qos=dev
#SBATCH --partition=dev
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=30 # number of cores
#SBATCH --mem=16G # memory pool for all cores
#SBATCH --time=8:00:00 # time (D-HH:MM)
#SBATCH --array=1-3

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
FASTQ_DIR=${BASE_DIR}/data/rnaseq/fastq/run${SLURM_ARRAY_TASK_ID}

# Full path to the FastQC output directory.
OUTPUT_DIR=${BASE_DIR}/results/fastqc

# Full path to the FastQC and MultiQC containers.
FASTQC=~/containers/fastqc.sif
MULTIQC=~/containers/multiqc.sif

### PROGRAM ###

module load singularity

mkdir -p ${BASE_DIR}/results/fastqc

# Run FastQC on all fastq files.
echo "Processing run${SLURM_ARRAY_TASK_ID}"
singularity exec ${FASTQC} fastqc \
                              --threads 30 \
                              --outdir ${OUTPUT_DIR} \
                              ${FASTQ_DIR}/SO*.fastq.gz

# Run MultiQC to summarize FastQC results.
singularity exec ${MULTIQC} multiqc --outdir ${OUTPUT_DIR} \
                                    --filename multiqc_report${SLURM_ARRAY_TASK_ID}.html \
                                    ${OUTPUT_DIR}

