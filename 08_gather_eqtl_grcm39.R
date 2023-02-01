################################################################################
# Gather eQTL results once they have been run in parallel. GRCm39.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-14
################################################################################
options(stringsAsFactors = FALSE)

##### VARIABLES #####

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
star_rsem_dir = file.path(base_dir, 'results', 'star_rsem')
results_dir = file.path(base_dir, 'results', 'eqtl_grcm39')
marker_dir  = '/projects/omics_share/mouse/GRCm39/supporting_files/muga_annotation/gigamuga'
marker_file = file.path(marker_dir, 'gm_uwisc_v2.csv')

# This assumes that the 07_local_eqtl.R script has been run and that
# the *rds output files exist.
eqtl_files = dir(path = results_dir, pattern = 'rds$', full.names = TRUE)

annot = readRDS(file.path(star_rsem_dir, 'gene_annotation_ens105.rds'))

##### MAIN #####

# Read in markers.
markers = read.csv(marker_file)
markers = markers[,c(1:2, 4)]
rownames(markers) = markers$marker
colnames(markers)[3] = 'pos'
markers$pos = markers$pos * 1e-6

# Read in the eqtl results files.
results = lapply(eqtl_files, readRDS)

# Bind the results into one object.
results = do.call(cbind, results)

# Write out eqtl results.
saveRDS(results, file.path(results_dir, '..', 'superovulation_eqtl_lod_grcm39.rds'))

#results = readRDS(file.path(results_dir, '..', 'superovulation_eqtl_lod.rds'))

# Subset the markers to retain the ones used in mapping.
markers = markers[rownames(results),]

stopifnot(rownames(markers) == rownames(results))

# For each gene, get the maximum LOD and it's location.
max_lod     = apply(results, 2, max)
max_lod_pos = apply(results, 2, which.max)

# Synch up gene names between results and annotation.
annot = annot[colnames(results),]
stopifnot(colnames(results) == rownames(annot))

gene_chr = annot$chr
gene_pos = annot$start * 1e-6

symbol = annot$symbol

max_qtl = data.frame(ensembl  = colnames(results),
                     symbol   = symbol,
                     gene_chr = gene_chr,
                     gene_pos = gene_pos,
                     qtl_chr  = markers$chr[max_lod_pos],
                     qtl_pos  = markers$pos[max_lod_pos],
                     lod      = max_lod)

write.csv(max_qtl, file = file.path(results_dir, '..', 'sodo_max_eqtl_summary_grcm39.csv'),
          row.names = FALSE, quote = FALSE)

