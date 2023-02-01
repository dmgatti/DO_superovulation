#!/bin/bash
#SBATCH --qos=batch
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=4G # memory pool for all cores
#SBATCH --time=10:00:00 # time (D-HH:MM)
#SBATCH --array=2-50

################################################################################
# Run GBRS on each sample.
# GRCm38.p6. Gencode M23. 
#
# Daniel Gatti
# 2022-09-28
################################################################################


##### VARIABLES #####

# Base directory for project.
BASE_DIR=/projects/bolcun-filas-lab/DO_Superovulation

# Base /fastscratch directory.
FASTSCRATCH=/fastscratch/dgatti

# FASTQ file directory.
FASTQ_DIR=${FASTSCRATCH}/fastq

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

# Output directory.
OUT_DIR=${FASTSCRATCH}/results/${SAMPLE}

# Final results directory (in backed up space)
DEST_DIR=${BASE_DIR}/results/gbrs/${SAMPLE}

# GBRS Nextflow pipeline
GBRS_PATH=${FASTSCRATCH}/gbrs/main.nf

# Number of threads.
NUM_THREADS=8

# Temporary working directory.
TMP_DIR=${FASTSCRATCH}/tmp/${SAMPLE}

# Singularity cache directory.
export NXF_SINGULARITY_CACHEDIR=${FASTSCRATCH}/singularity_cache

##### MAIN #####

mkdir -p ${OUT_DIR}
mkdir -p ${DEST_DIR}
mkdir -p ${TMP_DIR}
mkdir -p ${NXF_SINGULARITY_CACHEDIR}

module load singularity

cd ${TMP_DIR}

nextflow run ${GBRS_PATH} -profile singularity,sumner \
                          -resume \
                          --fastqR1 ${FASTQ_FILES[0]} \
                          --fastqR2 ${FASTQ_FILES[1]} \
                          --outputDir ${OUT_DIR} \
                          --generation G40 \
                          --sex F \
                          --threads ${NUM_THREADS} \
                          --tmpdir ${TMP_DIR}

# Copy output files from OUT_DIR to DEST_DIR.
SAMPLE_PATH=${OUT_DIR}/`ls ${OUT_DIR}`
#cp ${SAMPLE_PATH}/*counts ${DEST_DIR}
#cp ${SAMPLE_PATH}/*.tsv ${DEST_DIR}
#cp ${SAMPLE_PATH}/*.tpm ${DEST_DIR}
#cp ${SAMPLE_PATH}/*.pdf ${DEST_DIR}

#QC_FILE=`find . -name *summary_stats.txt`
#cp ./${QC_FILE} ${DEST_DIR}

# Cleanup
# Remove TMP_DIR
cd ../..
rm -rf ${TMP_DIR}
