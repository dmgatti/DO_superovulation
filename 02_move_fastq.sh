#!/bin/bash
#SBATCH --qos=batch
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=48G # memory pool for all cores
#SBATCH --time=0-12:00 # time (D-HH:MM)

################################################################################
# Copy DO FASTQ files over to /fastscratch.
# Concatenate samples with more than two files to two files.
# Downsample files with too many reads to 40 million reads.
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

# Copy the files to fastscratch.
echo Moving Run 1
cp ${SRC_DIR1}/SODO*.fastq.gz ${DEST_DIR}

echo Moving Run 2
cp ${SRC_DIR2}/SODO*.fastq.gz ${DEST_DIR}


# Loop through each sample ID and get the FASTQ files associated with it.
# If there are more than two files, combine them into two files.
# If the files are larger than 6 GB, subsample them.
echo "Checking for replicate and large files"
for i in {1..350}
do
  
  # Get FASTQ files for this sample.
  FASTQ_FILES=`ls -1 ${DEST_DIR}/*.gz | grep -E "SODO(_)?${i}_"`
  FASTQ_FILES=($FASTQ_FILES)

  if [ -z $FASTQ_FILES ]
  then
     echo No Sample ${i}
     continue
  fi
  
  SAMPLE=`basename ${FASTQ_FILES[0]}`
  SAMPLE=`echo ${SAMPLE} | sed -E "s/_GT22-[0-9]+_[ACGT]+-[ACGT]+_S[0-9]+_L00[0-9]+_R[0-9]_001.fastq.gz$//"`

  echo SAMPLE=${SAMPLE}

  if [ ${#FASTQ_FILES[@]} -eq 4 ]
  then
    echo Four files for sample ${SAMPLE}
  
    # Create new temporary filenames.
    NEW_FQ1=`dirname ${FASTQ_FILES[0]}`/${SAMPLE}_R1_001.fastq.gz
    NEW_FQ2=`dirname ${FASTQ_FILES[0]}`/${SAMPLE}_R2_001.fastq.gz
  
    # Concatenate the read 1 and read 2 files.
    cat ${FASTQ_FILES[0]} ${FASTQ_FILES[2]} > ${NEW_FQ1}
    cat ${FASTQ_FILES[1]} ${FASTQ_FILES[3]} > ${NEW_FQ2}
  
    # Remove the old files.
    rm ${FASTQ_FILES[@]}
  
    # Rename the new files to match the first two old files.
    mv ${NEW_FQ1} ${FASTQ_FILES[0]}
    mv ${NEW_FQ2} ${FASTQ_FILES[1]}

  fi

  # Get the current filesize. 
  FILESIZE=`stat -c %s ${FASTQ_FILES[0]}`

  if [ ${FILESIZE} -gt 6000000000 ]
  then

    echo Subsampling ${SAMPLE}

    # Seqtk container.
    SEQTK=~/containers/seqtk_latest.sif

    # New temporary fastq filenames.
    NEW_FQ1=`dirname ${FASTQ_FILES[0]}`/${SAMPLE}_R1_001.fastq
    NEW_FQ2=`dirname ${FASTQ_FILES[0]}`/${SAMPLE}_R2_001.fastq

    # Subsample to 40 million reads and zip up new fastq file.
    singularity exec ${SEQTK} seqtk sample -s1234567 ${FASTQ_FILES[0]} 40000000 | gzip > ${NEW_FQ1}.gz
    singularity exec ${SEQTK} seqtk sample -s1234567 ${FASTQ_FILES[1]} 40000000 | gzip > ${NEW_FQ2}.gz

    # Move subsampled file to original filename.
    mv ${NEW_FQ1}.gz ${FASTQ_FILES[0]}
    mv ${NEW_FQ2}.gz ${FASTQ_FILES[1]}

  fi
done
