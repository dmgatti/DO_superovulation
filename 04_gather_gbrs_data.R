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

library(AnnotationHub)
library(ensembldb)
library(DESeq2)
library(qtl2convert)
library(qtl2)
library(readxl)
library(tidyverse)

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)

genome = args[1]

base_dir = '/projects/bolcun-filas-lab/DO_Superovulation'
gbrs_dir = file.path(base_dir, 'results', 'gbrs', genome)
out_dir  = file.path(base_dir, 'results')

marker_dir  = '/projects/omics_share/mouse/GRCm38/supporting_files/gbrs/rel_1505_v5'
marker_file = file.path(marker_dir, 'ref.genome_grid.69k.noYnoMT_KBEdit.txt')

if(genome == 'grcm39') {

  marker_dir  = gbrs_dir
  marker_file = file.path(marker_dir, 'ref.genome_grid.GRCm39.tsv')

} # if(genome == 'grcm39')

# GBRS docs say that they used GencodeM23 for transcript annotation.
# This corresponds to Ensembl 98 according to the header in:
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
hub   = AnnotationHub()
ensdb = NULL

if(genome == 'grcm38') {

  hub   = query(hub, c('ensdb', 'mus musculus', '98'))
  ensdb = hub[['AH75036']] # ensembl 98
  #ensdb = hub[['AH78811']] # ensembl 99

} else {

  # New GRCm39 GBRS pipeline uses Ensembl 105.
  hub   = query(hub, c('ensdb', 'mus musculus', '105'))
  ensdb = hub[['AH98078']] # ensembl 105

} # else

##### MAIN #####

# Get gene annotation.
annot = genes(ensdb)

# Read in phenotypes.
pheno = read_xlsx(file.path(base_dir, 'data', 'JDO9376 all 350.xlsx'),
                  sheet = 'Sheet1')

# Read in markers.
markers = read_delim(marker_file, delim = '\t', show_col_types = FALSE)

# Get sample directories.
sample_dirs = dir(gbrs_dir, pattern = 'SODO')

# Gather expected counts for genes & isoforms, genotypes at genes,
# and genoprobs.
sample_dir = file.path(gbrs_dir, sample_dirs[1])
count_file = dir(sample_dir, pattern = 'diploid.genes.expected_read_counts$',
                 full.names = TRUE)
iso_file   = dir(sample_dir, pattern = 'diploid.isoforms.expected_read_counts$',
                 full.names = TRUE) 
ct  = read_delim(count_file, delim = '\t', show_col_types = FALSE)
iso = read_delim(iso_file,   delim = '\t', show_col_types = FALSE)
genes       = ct$`locus`
transcripts = iso$`locus`

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
    print('DIRECTORY WAS EMPTY!!!')
    next
  }
  
  ct = read_delim(count_file, delim = '\t', show_col_types = FALSE)
  gene_counts[d,,] = t(as.matrix(ct[,LETTERS[1:8]]))

  # Genotype at gene.
  geno[,d] = ct$notes

  # Expected isoform counts file.
  iso_file = dir(sample_dir, pattern = 'diploid.isoforms.expected_read_counts$',
                 full.names = TRUE)
  
  iso = read_delim(iso_file, delim = '\t', show_col_types = FALSE)
  iso_counts[d,,] = t(as.matrix(iso[,LETTERS[1:8]]))

  # Genoprobs.
  probs_file = dir(sample_dir, pattern = 'gbrs.interpolated.genoprobs.tsv$',
                 full.names = TRUE)
  pr = read_delim(probs_file, delim = '\t', show_col_types = FALSE)
  probs[d,,] = t(as.matrix(pr))

} # for(d)

# Remove '_' from sample names.
rownames(gene_counts) = sub('_', '', rownames(gene_counts))
rownames(iso_counts)  = sub('_', '', rownames(iso_counts))
rownames(probs)       = sub('_', '', rownames(probs))

# If there are samples for which all counts equal zero, then set the counts 
# for that sample to NA.
wh = which(rowSums(gene_counts) == 0)

print(paste(paste(rownames(gene_counts)[wh], sep = ','), sep = ' '))

if(length(wh) > 0) {

  gene_counts = gene_counts[-wh,,]
  iso_counts  = iso_counts[-wh,,]
  probs       = probs[-wh,,]
  
} # if(length(wh) > 0)

# Write out results files.
saveRDS(gene_counts, file = file.path(out_dir, paste0('sodo_gene_allele_counts_', genome, '.rds')))
saveRDS(iso_counts,  file = file.path(out_dir, paste0('sodo_transcript_allele_counts_', genome, '.rds')))
saveRDS(probs,       file = file.path(out_dir, paste0('sodo_gbrs_genoprobs_69K_', genome, '.rds')))
geno = as.data.frame(geno) %>% 
         rownames_to_column(var = 'ensembl')
write_csv(geno, file = file.path(out_dir, paste0('sodo_gene_genotypes_', genome, '.csv')))

# Gather total gene counts.
counts = apply(gene_counts, c(1, 3), sum, na.rm = TRUE)
counts = as.data.frame(t(counts))
counts = counts %>%
           rownames_to_column(var = 'ensembl')

# Write gene annotation.
common_genes = intersect(names(annot), counts$ensembl)
annot  = subset(annot,  names(annot)   %in% common_genes)
counts = subset(counts, counts$ensembl %in% common_genes)
annot  = as.data.frame(annot)
annot  = annot[counts$ensembl,]
annot  = select(annot, chr = seqnames, start, end, strand, gene_id, gene_name,                                            gene_biotype, symbol, entrezid)

stopifnot(counts$ensembl == rownames(annot))

if(genome == 'grcm38') {

  write_csv(annot, file = file.path(out_dir, 'gene_annotation_ens98.csv'))

} else {

  write_csv(annot, file = file.path(out_dir, 'gene_annotation_ens105.csv'))

} # else

# Write gene counts.
stopifnot(nrow(counts) > 0)
write_csv(counts, file = file.path(out_dir, paste0('sodo_gene_counts_', genome, '.csv')))

# Get VST normalized data.
counts = data.frame(counts)
rownames(counts) = counts$ensembl
counts = as.matrix(counts[,-1])
counts = round(counts)

metadata = pheno[, c('Female number', 'FreshIVF_IVFDISHES::IVF Date')]
colnames(metadata) = c('mouse', 'date')
metadata$mouse = paste0('SODO', metadata$mouse)
metadata = as.data.frame(metadata)
rownames(metadata) = metadata$mouse
metadata$date = factor(metadata$date)
metadata = metadata[colnames(counts),]

# Normalize counts using DESeq's VST.
dds = DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ date)
dds = DESeq(dds)
vst = assay(vst(dds))
vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ensembl') %>%
  write_csv(file = file.path(out_dir, paste0('sodo_gene_counts_norm_', genome, '.csv')))

# Gather total transcript counts.
counts = apply(iso_counts, c(1, 3), sum, na.rm = TRUE)
counts = as.data.frame(t(counts))
counts = counts %>%
           rownames_to_column(var = 'ensembl')
write_csv(counts, file = file.path(out_dir, paste0('sodo_transcript_counts_', genome, '.csv')))

