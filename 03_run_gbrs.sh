#!/bin/bash
#SBATCH --qos=batch
#SBATCH --ntasks=1 # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=1G # memory pool for all cores
#SBATCH --time=48:00:00 # time (D-HH:MM)

################################################################################
# Run GBRS on each sample.
# GRCm39. Ensembl 105. 
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-09-28
################################################################################


##### VARIABLES #####

# Base directory for project.
BASE_DIR=/projects/bolcun-filas-lab/DO_Superovulation

# Sample metadata file.
SAMPLE_META_FILE=${BASE_DIR}/data/sodo_gbrs_metadata.csv

# Base /fastscratch directory.
FASTSCRATCH=/fastscratch/dgatti

# Temporary directory.
TMP_DIR=${FASTSCRATCH}/tmp

# Output directory.
OUT_DIR=${FASTSCRATCH}/results

# Final results directory (in backed up space)
DEST_DIR=${BASE_DIR}/results/gbrs/grcm39

# GBRS Reference Directory
GBRS_REF_DIR=/projects/churchill-lab/projects/GBRS_GRCm39

# GBRS Nextflow pipeline.
GBRS_PATH=${GBRS_REF_DIR}/github_tool_clones/ngs-ops-nf-pipelines/main.nf

# GBRS Transcripts.
GBRS_TRANSCRIPTS=${GBRS_REF_DIR}/eight_way_transcriptome/full_transcriptome/primary_assembly/emase/emase.fullTranscripts.info

# GBRS/EMASE gene to transcript file.
GBRS_GENE2TRANS=${GBRS_REF_DIR}/eight_way_transcriptome/full_transcriptome/primary_assembly/emase/emase.gene2transcripts.tsv

# GBRS Full Transcript.
GBRS_FULL_TRANSCRIPTS=${GBRS_REF_DIR}/eight_way_transcriptome/full_transcriptome/primary_assembly/emase/emase.pooled.fullTranscripts.info

# GBRS Emission Probs.
GBRS_EMIS_PROBS=${GBRS_REF_DIR}/emission_probabilities/gbrs_emissions_all_tissues.avecs.npz

# GBRS Transmission Probs.
GBRS_TRANS_PROBS=${GBRS_REF_DIR}/transition_probabilities/full_transcriptome/nextflow_output

# Ensembl 105 gene positions.
ENSEMBL_105=${GBRS_REF_DIR}/transition_probabilities/full_transcriptome/nextflow_output/ref.gene_pos.ordered_ensBuild_105.npz

# GBRS 69K Marker Grid.
MARKER_GRID=${GBRS_REF_DIR}/transition_probabilities/full_transcriptome/nextflow_output/ref.genome_grid.GRCm39.tsv

# Bowtie index for GBRS.
BOWTIE_INDEX=/projects/churchill-lab/projects/GBRS_GRCm39/eight_way_transcriptome/primary_assembly/bowtie/bowtie.transcripts

# Singularity cache directory.
export NXF_SINGULARITY_CACHEDIR=${FASTSCRATCH}/singularity_cache

##### MAIN #####

mkdir -p ${OUT_DIR}
mkdir -p ${DEST_DIR}
mkdir -p ${TMP_DIR}
mkdir -p ${NXF_SINGULARITY_CACHEDIR}

module load singularity

cd ${TMP_DIR}

nextflow ${GBRS_PATH} \
         -profile sumner \
         --workflow gbrs \
         --pubdir ${OUT_DIR} \
         -w ${TMP_DIR} \
         --bowtie_index ${BOWTIE_INDEX} \
         --csv_input ${SAMPLE_META_FILE} \
         --transcripts_info ${GBRS_TRANSCRIPTS} \
         --gene2transcript_csv ${GBRS_GENE2TRANS} \
         --full_transcript_info ${GBRS_FULL_TRANSCRIPTS} \
         --emission_prob_avecs ${GBRS_EMIS_PROBS} \
         --trans_prob_dir ${GBRS_TRANS_PROBS} \
         --gene_position_file ${ENSEMBL_105} \
         --genotype_grid ${MARKER_GRID} \
         -resume

# Copy output files from OUT_DIR to DEST_DIR.
DIRS=`ls ${OUT_DIR}`

for i in ${DIRS}
do

  CURR_DEST_DIR=${DEST_DIR}/$i

  mkdir -p ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/stats/* ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*_counts ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*.tsv ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*.pdf ${CURR_DEST_DIR}

done



