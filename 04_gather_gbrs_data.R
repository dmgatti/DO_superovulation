################################################################################
# Gather the GBRS data together.
# Create counts matrix.
# Create genotype-at-gene matrix.
# Create qtl2-style genoprobs object.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-13
################################################################################
options(stringsAsFactors = FALSE)

library(qtl2convert)
library(qtl2)
library(tidyverse)

##### VARIABLES #####

base_dir = '/projects/bolcun-filas-lab/DO_Superovulation'
gbrs_dir = file.path(base_dir, 'results', 'gbrs')
out_dir  = file.path(base_dir, 'results')

marker_dir  = '/projects/omics_share/mouse/GRCm38/supporting_files/gbrs/rel_1505_v5'
marker_file = file.path(marker_dir, 'ref.genome_grid.69k.noYnoMT_KBEdit.txt')

##### MAIN #####

# Read in markers.
markers = read_delim(marker_file, delim = '\t')

# Get sample directories.
sample_dirs = dir(gbrs_dir)

# Gather expected counts for genes & isoforms, genotypes at genes,
# and genoprobs.
sample_dir = file.path(gbrs_dir, sample_dirs[1])
count_file = dir(sample_dir, pattern = 'diploid.genes.expected_read_counts$',
                 full.names = TRUE)
iso_file   = dir(sample_dir, pattern = 'diploid.isoforms.expected_read_counts$',
                 full.names = TRUE) 
ct  = read_delim(count_file, delim = '\t')
iso = read_delim(iso_file,   delim = '\t')
genes       = ct$`#target_id`
transcripts = iso$`#target_id`

# Create output data structures.
gene_counts = array(0, dim = c(length(sample_dirs), 8, length(genes)), 
                    dimnames = list(sample_dirs, LETTERS[1:8], genes))
iso_counts  = array(0, dim = c(length(sample_dirs), 8, length(transcripts)), 
                    dimnames = list(sample_dirs, LETTERS[1:8], transcripts))
geno        = matrix('', nrow = length(genes), ncol = length(sample_dirs), 
                     dimnames = list(genes, sample_dirs))
probs       = array(0, dim = c(length(sample_dirs), 8, nrow(markers)), 
                    dimnames = list(sample_dirs, LETTERS[1:8], markers$marker))

for(d in sample_dirs) {

  print(d)

  sample_dir = file.path(gbrs_dir, d)

  # Expected gene counts file.
  count_file = dir(sample_dir, pattern = 'diploid.genes.expected_read_counts$',
                   full.names = TRUE)
  
  if(length(count_file) == 0) {
    next
  }
  
  ct = read_delim(count_file, delim = '\t')
  gene_counts[d,,] = t(as.matrix(ct[,LETTERS[1:8]]))

  # Genotype at gene.
  geno[,d] = ct$notes

  # Expected isoform counts file.
  iso_file = dir(sample_dir, pattern = 'diploid.isoforms.expected_read_counts$',
                 full.names = TRUE)
  
  iso = read_delim(iso_file, delim = '\t')
  iso_counts[d,,] = t(as.matrix(iso[,LETTERS[1:8]]))

  # Genoprobs.
  probs_file = dir(sample_dir, pattern = 'gbrs.interpolated.genoprobs.tsv$',
                 full.names = TRUE)
  pr = read_delim(probs_file, delim = '\t')
  probs[d,,] = t(as.matrix(pr))

} # for(d)

# Write out results files.
saveRDS(gene_counts, file = file.path(out_dir, 'sodo_gene_allele_counts.rds'))
saveRDS(iso_counts,  file = file.path(out_dir, 'sodo_transcript_allele_counts.rds'))
saveRDS(probs,       file = file.path(out_dir, 'sodo_gbrs_genoprobs_69K.rds'))
write_csv(as.data.frame(geno), file = file.path(out_dir, 'sodo_gene_genotypes.csv'))

# Gather total gene and transcript counts.
counts = apply(gene_counts, c(1, 3), sum, na.rm = TRUE)
counts = as.data.frame(t(counts))
counts = counts %>%
           rownames_to_column(var = 'ensembl')
write_csv(counts, file = file.path(out_dir, 'sodo_gene_counts.csv'))

counts = apply(iso_counts, c(1, 3), sum, na.rm = TRUE)
counts = as.data.frame(t(counts))
counts = counts %>%
           rownames_to_column(var = 'ensembl')
write_csv(counts, file = file.path(out_dir, 'sodo_transcript_counts.csv'))

