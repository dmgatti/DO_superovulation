################################################################################
# Haplotype Reconstruction on GRCm38 forDO Superovulation study.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-12-06
################################################################################

library(qtl2convert)
library(qtl2)

args = commandArgs(trailingOnly = TRUE)

genome = args[1]

base_dir = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir = file.path(base_dir, 'data')

scratch_dir = '/fastscratch/dgatti'
qtl2_dir    = file.path(scratch_dir, 'sodo_qtl2')

# Haplotype reconstruction.
print('Reading cross')
cross = read_cross2(file.path(qtl2_dir, paste0('do_superovulation_', genome, '.json')))

print(cross)

print('Writing cross object.')
saveRDS(cross, file = file.path(data_dir, paste0('sodo_cross_', genome, '.rds')))

print('Running HR')
probs = calc_genoprob(cross, map = cross$gmap, cores = 20, quiet = FALSE)
print('Writing 36 state genoprobs')
saveRDS(probs, file = file.path(data_dir, paste0('sodo_genoprobs_', genome, '_123K.rds')))

print('Converting 36 state genoprobs to 8 state allele probs')
aprobs = genoprob_to_alleleprob(probs, cores = 20, quiet = FALSE)
print('Writing allele probs')
saveRDS(aprobs, file = file.path(data_dir, paste0('sodo_alleleprobs_', genome, '_123K.rds')))

