################################################################################
# Gather eQTL results once they have been run in parallel.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-14
################################################################################
options(stringsAsFactors = FALSE)

library(AnnotationHub)

##### VARIABLES #####

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
results_dir = file.path(base_dir, 'results', 'eqtl')
marker_dir  = '/projects/omics_share/mouse/GRCm38/supporting_files/muga_annotation/gigamuga'
marker_file = file.path(marker_dir, 'gm_uwisc_v1.csv')

# This assumes that the 07_local_eqtl.R script has been run and that
# the *rds output files exist.
eqtl_files = dir(path = results_dir, pattern = 'rds$', full.names = TRUE)

hub = AnnotationHub()
hub = query(hub, c('mus musculus', 'ensdb', '99'))
ensdb = hub[['AH78811']]

ens_genes = genes(ensdb)

##### MAIN #####

# Read in markers.
markers = read.csv(marker_file)
markers = markers[,1:3]
rownames(markers) = markers$marker
colnames(markers)[3] = 'pos'
markers$pos = markers$pos * 1e-6

# Read in the eqtl results files.
results = lapply(eqtl_files, readRDS)

# Bind the results into one object.
results = do.call(cbind, results)

# Write out eqtl results.
saveRDS(results, file.path(results_dir, '..', 'superovulation_eqtl_lod.rds'))

#results = readRDS(file.path(results_dir, '..', 'superovulation_eqtl_lod.rds'))

# Subset the markers to retain the ones used in mapping.
markers = markers[rownames(results),]

stopifnot(rownames(markers) == rownames(results))

# For each gene, get the maximum LOD and it's location.
max_lod     = apply(results, 2, max)
max_lod_pos = apply(results, 2, which.max)

# There are a few genes that are in the results,but not in the annotation.
m = match(colnames(results), names(ens_genes))
wh_na = which(is.na(m))
gene_chr = rep(NA_character_, ncol(results))
gene_chr[(1:ncol(results))[-wh_na]] = seqnames(ens_genes)[m[-wh_na]]

gene_pos = rep(NA_real_, ncol(results))
gene_pos[(1:ncol(results))[-wh_na]] = start(ens_genes)[m[-wh_na]] * 1e-6

symbol = rep(NA_character_, ncol(results))
symbol[(1:ncol(results))[-wh_na]] = ens_genes$gene_name[m[-wh_na]]

max_qtl = data.frame(ensembl  = colnames(results),
                     symbol   = symbol,
                     gene_chr = gene_chr,
                     gene_pos = gene_pos,
                     qtl_chr  = markers$chr[max_lod_pos],
                     qtl_pos  = markers$pos[max_lod_pos],
                     lod      = max_lod)

write.csv(max_qtl, file = file.path(results_dir, '..', 'sodo_max_eqtl_summary.csv'),
          row.names = FALSE, quote = FALSE)

