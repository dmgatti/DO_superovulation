#!/bin/bash
#SBATCH --qos batch
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of tasks
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G # memory pool for all cores
#SBATCH --time 0-12:00 # time (D-HH:MM)
#SBATCH --array=347

################################################################################
# Use RSEM & STAR to estimate gene counts aligned to  GRCm39.
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

# Base directory for results.
BASE_DIR=/projects/bolcun-filas-lab/DO_Superovulation

# Results directory.
RESULTS_DIR=${BASE_DIR}/results/star_rsem_${GENOME}

# FASTQ file directory.
FASTQ_DIR=${SCRATCH_DIR}/fastq

# Number of threads.
NUM_THREADS=8

# Sample ID.
FASTQ_FILES=`ls -1 ${FASTQ_DIR}/*.gz | grep -E "SODO(_)?${SLURM_ARRAY_TASK_ID}_"`
FASTQ_FILES=($FASTQ_FILES)

if [ -z $FASTQ_FILES ]
then
   echo No Sample ${SLURM_ARRAY_TASK_ID}
   exit 0
fi

SAMPLE=`basename ${FASTQ_FILES[0]}`
SAMPLE=`echo ${SAMPLE} | sed -E "s/_GT22-[0-9]+_[ACGT]+-[ACGT]+_S[0-9]+_L00[0-9]+_R[0-9]_001.fastq.gz$//"`

echo SAMPLE=${SAMPLE}

# Some samples have multiple files. Separate out the R1 and R2 files. 
R1_FILES=`ls -1 ${FASTQ_DIR}/*.gz | grep -E "SODO(_)?${SLURM_ARRAY_TASK_ID}_.+R1"`
R2_FILES=`ls -1 ${FASTQ_DIR}/*.gz | grep -E "SODO(_)?${SLURM_ARRAY_TASK_ID}_.+R2"`

R1_FILES=${R1_FILES[@]}
R2_FILES=${R2_FILES[@]}

R1_FILES=${R1_FILES/ /,}
R2_FILES=${R2_FILES/ /,}

# RSEM output directory.
OUT_DIR=${SCRATCH_DIR}/sodo/${GENOME}

# Temporary working directory.
TEMP_DIR=${SCRATCH_DIR}/temp/${SAMPLE}

# STAR/RSEM Singularity container.
STAR_RSEM=~/containers/star_rsem.sif

### PROGRAM ###

mkdir -p ${TEMP_DIR}
mkdir -p ${OUT_DIR}
mkdir -p ${RESULTS_DIR}

module load singularity

echo Aligning sample ${SAMPLE}

echo ${R1_FILES[0]}
echo ${R2_FILES[1]}

# Using STAR to align and  RSEM to quantify.
singularity exec ${STAR_RSEM} rsem-calculate-expression \
   --paired-end \
   --num-threads ${NUM_THREADS} \
   --star \
   --estimate-rspd \
   --append-names \
   --star-gzipped-read-file \
   --temporary-folder ${TEMP_DIR} \
   ${R1_FILES} ${R2_FILES} \
   ${GENOME_DIR} \
   ${OUT_DIR}/${SAMPLE}_

# Copy results to a stable location.
echo Copying ${OUT_DIR}/${SAMPLE}_*.results and ${OUT_DIR}/${SAMPLE}_*.log to ${RESULTS_DIR}
cp ${OUT_DIR}/${SAMPLE}_*.results ${RESULTS_DIR}
cp ${OUT_DIR}/${SAMPLE}_*.log ${RESULTS_DIR}

# NOTE: STAR cleans up its temp directory.

