#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of tasks
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G # memory pool for all cores
#SBATCH --time 0-2:00 # time (D-HH:MM)

################################################################################
# Build RSEM/STAR genome index for GRCm39.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2022-12-06
################################################################################

##### VARIABLES #####

# Scratch directory.
SCRATCH_DIR=/fastscratch/dgatti

# Genome build
GENOME=grcm38
#GENOME=grcm39

# Ensembl version
ENSEMBL=102
#ENSEMBL=106

# Output directory for STAR/RSEM index.
GENOME_DIR=${SCRATCH_DIR}/${GENOME}_index/rsem_

# Reference annotation directory in omics_share
ANNOT_DIR=/projects/omics_share/mouse/${GENOME}

# fasta file directory for mouse genome
FASTA_DIR=${ANNOT_DIR}/genome/sequence/ensembl/v${ENSEMBL}

# Name of whole genome FASTA file (without gz suffix).
FASTA_FILE=${FASTA_DIR}/Mus_musculus.${GENOME}.dna.primary_assembly.fa

# Transcript GTF file from Ensembl.
GTF_FILE=${ANNOT_DIR}/transcriptome/annotation/ensembl/v${ENSEMBL}/Mus_musculus.${GENOME}.${ENSEMBL}.gtf

# Read length of FASTQ files.
READ_LENGTH=151

# RSEM/STAR Singularity container.
STAR_RSEM=~/containers/star_rsem.sif


##### MAIN #####

echo 'Building genome reference'

mkdir -p ${GENOME_DIR}

module load singularity

singularity exec ${STAR_RSEM} rsem-prepare-reference \
   --num-threads 8 \
   --star \
   --star-sjdboverhang $(($READ_LENGTH - 1)) \
   --gtf ${GTF_FILE} \
   ${FASTA_FILE} \
   ${GENOME_DIR}

