################################################################################
# Gather the STAR/RSEM counts into one file.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-12-07
################################################################################
options(stringsAsfactors = FALSE)

library(AnnotationHub)
library(DESeq2)
library(ensembldb)
library(tidyverse)

##### VARIABLES #####

genome = 'GRCm39'
ensembl_version = 105

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results', paste0('star_rsem_', genome))

# Get Ensembl  genes.
hub = AnnotationHub()
hub = query(hub, c('mus musculus', 'ensdb', ensembl_version))
ensembl = hub[[1]]

##### MAIN #####

# Read in the gene and isoform files.
gene_files = dir(results_dir, pattern = 'genes\\.results$',    full.names = TRUE)
iso_files  = dir(results_dir, pattern = 'isoforms\\.results$', full.names = TRUE)

genes = lapply(gene_files, read_tsv, col_select = c('gene_id', 'expected_count'))
iso   = lapply(iso_files,  read_tsv, col_select = c('transcript_id', 'gene_id', 
               'expected_count'))

# Collapse genes into a matrix.
gene_names = genes[[1]]$gene_id
genes = lapply(genes, function(z) { dplyr::select(z, 'expected_count') })
genes = do.call(cbind, genes)
dimnames(genes) = list(gene_names,
                       str_replace(basename(gene_files), 
                       pattern = '\\.genes\\.results$', replacement = ''))
colnames(genes) = gsub('_', '', colnames(genes))

# Collapse transcripts into a matrix.
trans_names = iso[[1]][,c('transcript_id', 'gene_id')]
trans = lapply(iso, function(z) { dplyr::select(z, 'expected_count') })
trans = do.call(cbind, trans)
trans = cbind(trans_names, trans)
dimnames(trans) = list(pull(trans_names, transcript_id),
                       c('transcript_id', 'gene_id', 
                       str_replace(basename(iso_files), 
                       pattern = '\\.isoforms\\.results$', replacement = '')))
colnames(trans) = gsub('_', '', colnames(trans))

# Write out raw counts files.
#####
##### NOTE: Commented out so that I don't accidently overwrite.
#####
saveRDS(genes, file = file.path(results_dir, paste0('star_rsem_gene_expected_counts_', genome,'.rds')))
saveRDS(trans, file = file.path(results_dir, paste0('star_rsem_transcript_expected_counts_', genome, '.rds')))

# Normalize gene counts using DESeq2 VST.

# Read in sample metadata.
meta = read_csv(file.path(data_dir, 'do_superovulation_animal_info.csv')) %>%
         dplyr::select(id = `Animal ID`, gen = `DO generation`) %>%
         mutate(gen = as.factor(gen),
                id  = str_replace(id, '_T$', ''),
                id  = str_replace(id, '^SO', 'SODO')) %>%
         column_to_rownames(var = 'id')

meta = meta[colnames(genes),,drop = FALSE]
stopifnot(rownames(meta) == colnames(genes))

# Remove gene symbol from gene rownames.
rownames(genes) = sub('_[A-Z,a-z,0-9,_,(,),\\.,\\-]+$', '', rownames(genes))

# Filter out genes with zero counts in less than 20% of samples.
keep = which(rowMeans(genes > 0) > 0.2)
genes = genes[keep,]

print(str_c('Number of genes: ', nrow(genes)))

# Create DESeq object.
dds = DESeqDataSetFromMatrix(countData = round(genes), 
                             colData   = meta, 
                             design    = ~ gen)
dds = DESeq(dds)

norm = assays(vst(dds))[[1]]

saveRDS(norm, file = file.path(results_dir, paste0('star_rsem_gene_vst_counts_', genome, '.rds')))

# Write out gene annotation.
annot = ensembldb::select(x       = ensembl, 
               keys    = rownames(genes), 
               columns = c('SEQNAME', 'GENESEQSTART', 'GENESEQEND',  'SEQSTRAND',
                           'GENEID',  'GENENAME',     'GENEBIOTYPE', 'SYMBOL',
                           'ENTREZID'), 
               keytype = 'GENEID')
dup = unique(annot$GENEID[duplicated(annot$GENEID)])

for(g in dup) {

  rows = which(annot$GENEID == g)
  annot$ENTREZID[rows[1]] = str_c(annot$ENTREZID[rows], collapse = ';')
  annot = annot[-rows[-1],]

} # for(g)

rownames(annot) = annot$GENEID

stopifnot(all(rownames(genes) == rownames(annot)))

# Set colnames to what the QTL viewer may expect.
colnames(annot) = c('chr',       'start',         'end',     'strand', 'gene_id', 
                    'gene_name', 'gene_biotype', 'symbol', 'entrezid')
annot = annot %>%
          mutate(strand = as.character(strand),
                 strand = if_else(strand == '1', '+', '-'))
saveRDS(annot, file = file.path(results_dir, paste0('gene_annotation_ens',  
        ensembl_version, '.rds')))







