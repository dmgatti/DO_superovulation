#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 1 # number of tasks
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G # memory pool for all cores
#SBATCH --time 0-2:00 # time (D-HH:MM)
#SBATCH --array=2-3

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

# Output directory for STAR index.
STAR_INDEX=${SCRATCH_DIR}/rsem_index

# Output directory for RSEM index.
RSEM_INDEX=${SCRATCH_DIR}/rsem_index/rsem_index

# Base directory for results.
BASE_DIR=/projects/bolcun-filas-lab/DO_Superovulation

# Results directory.
RESULTS_DIR=${BASE_DIR}/results/star_rsem

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

# RSEM output directory.
OUT_DIR=${SCRATCH_DIR}/star_rsem

# Temporary working directory.
TEMP_DIR=${SCRATCH_DIR}/temp/${SAMPLE}

# STAR/RSEM Singularity container.
STAR_RSEM=~/containers/star_rsem.sif

### PROGRAM ###

mkdir -p ${TEMP_DIR}
mkdir -p ${OUT_DIR}

module load singularity

echo Aligning sample ${SAMPLE}

echo ${FASTQ_FILES[0]}
echo ${FASTQ_FILES[1]}

# Use STAR to align first, then call RSEM on the STAR BAMs.
# Setting mismatch filter to 8 per read.
singularity exec ${STAR_RSEM} STAR --runThreadN 8 \
                                   --readFilesIn ${FASTQ_FILES[0]} ${FASTQ_FILES} \
                                   --readFilesCommand zcat \
                                   --genomeDir ${STAR_INDEX} \
                                   --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_ \
                                   --outFilterType BySJout \
                                   --outSAMtype BAM SortedByCoordinate \
                                   --outSAMattributes Standard \
                                   --quantMode TranscriptomeSAM

BAM_FILE=${OUT_DIR}/${SAMPLE}_Aligned.toTranscriptome.out.bam

singularity exec ${STAR_RSEM} rsem-calculate-expression \
                                --paired-end \
                                --num-threads 8 \
                                --alignments \
                                --estimate-rspd \
                                --append-names \
                                --output-genome-bam \
                                --bam ${BAM_FILE} \
                                ${RSEM_INDEX} \
                                ${OUT_DIR}/${SAMPLE}_


# Old combined START/RSEM command.
#singularity exec ${STAR_RSEM} rsem-calculate-expression \
#   --paired-end \
#   --num-threads ${NUM_THREADS} \
#   --star \
#   --estimate-rspd \
#   --append-names \
#   --star-gzipped-read-file \
#   --temporary-folder ${TEMP_DIR} \
#   ${FASTQ_FILES[0]} ${FASTQ_FILES[1]} \
#   ${GENOME_DIR} \
#   ${OUT_DIR}/${SAMPLE}

cp ${OUT_DIR}/${SAMPLE}*.results ${RESULTS_DIR}

# Remove BAMs.
#rm ${OUT_DIR}/${SAMPLE}*.bam

#rm -r ${TEMP_DIR}
