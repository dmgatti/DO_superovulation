################################################################################
# Perform local eQTL mapping on all genes.
# Use the GBRS genoprobs so that we don't have to deal with sample mismatches.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-14
################################################################################
options(stringsAsFactors = FALSE)

library(qtl2convert)
library(qtl2)

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)
start = as.numeric(args[1])
end   = as.numeric(args[2])

print(paste('START =', start, ' END =', end))

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results', 'eqtl_grcm39')

# This assumes that the 06_gather_qtl_data.R script has been run and that
# 'superovulation_qtl2_grcm39_120K.Rdata' exists.
data_file = file.path(data_dir, 'superovulation_qtl2_grcm39_120K.Rdata')

##### MAIN #####

# Read in data.
load(data_file)

# Keep the genes that we are mapping based on the start and end.
end = min(end, nrow(norm_expr))

norm_expr = t(norm_expr[start:end,])

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

lod = scan1(genoprobs = probs, pheno = norm_expr, kinship = K, addcovar = addcovar)

print('Writing results file')

saveRDS(lod, file = file.path(results_dir, paste0('eqtl_lod_', start, '_', end, '.rds')))
