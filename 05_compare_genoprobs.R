################################################################################
# Compare the GigaMUGA and GBRS haplotype reconstructions.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-13
################################################################################
options(stringsAsFactors = FALSE)

library(qtl2convert)
library(qtl2)

source('interpolate_genoprobs.R')

##### VARIABLES #####

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results')
out_dir     = file.path(base_dir, 'results')

# Gigamuga marker file (GRCm38)
gm_marker_file = '/projects/omics_share/mouse/GRCm39/supporting_files/muga_annotation/gigamuga/gm_uwisc_v4.csv'

# 69K grid marker file.
#gbrs_marker_file = file.path(results_dir, 'ref.genome_grid.69k.noYnoMT_KBEdit.txt')
gbrs_marker_file = file.path(results_dir, 'ref.genome_grid.GRCm39.tsv')

# Comparison function.
# a & b: numeric matrix, founders in rows, markers in columns.
cos_sim = function(a, b) {
  v = rep(1L, nrow(a))
  mean((v %*% (a * b)) / sqrt(v %*% a^2 * v %*% b^2))
} # cos_sim()

##### MAIN #####

# Read in the Gigamuga markers.
gm_markers     = read.csv(gm_marker_file)
gm_markers$pos = gm_markers$bp_grcm39
gm_markers     = subset(gm_markers, gm_markers$chr %in% c(1:19, 'X'))
gm_map         = qtl2convert::map_df_to_list(gm_markers, pos_column = 'pos')

# Read in the 69K grid markers.
gbrs_markers = read.delim(gbrs_marker_file)
gbrs_map     = qtl2convert::map_df_to_list(gbrs_markers, pos_column = 'pos')

# Read in the Gigamuga genoprobs and interpolate to the 69K grid.
gm_probs = readRDS(file.path(data_dir, 'sodo_alleleprobs_grcm39_137K.rds'))
for(chr in seq_along(gm_probs)) {

  gm_map[[chr]] = gm_map[[chr]][dimnames(gm_probs[[chr]])[[3]]]

} # for(chr)

gm_probs = interpolate_genoprobs(gm_probs, gm_map, gbrs_map)
gm_probs = qtl2convert::probs_qtl2_to_doqtl(gm_probs)

# Read in the 69K genoprobs.
gbrs_probs = readRDS(file.path(results_dir, 'sodo_gbrs_genoprobs_69K_grcm39.rds'))

# Synch up sample IDs between the genoprobs.
samples = intersect(rownames(gm_probs), rownames(gbrs_probs)) 

gm_probs   = gm_probs[samples,,]
gbrs_probs = gbrs_probs[samples,,]

# Verify that the Gigamuga and GBRS genoprobs have the same structure.
stopifnot(rownames(gm_probs) == rownames(gbrs_probs))
stopifnot(colnames(gm_probs) == colnames(gbrs_probs))
stopifnot(dimnames(gm_probs)[[3]] == dimnames(gbrs_probs)[[3]])

# Compare GM and GBRS probs for each sample.
prob_cor = data.frame(id = rownames(gm_probs), sim = 0)

for(i in 1:nrow(gm_probs)) {

  prob_cor$sim[i] = cos_sim(gm_probs[i,,], gbrs_probs[i,,])

} # for(i) 

# Plot the similarities.
png(file.path(results_dir, 'genoprobs_cor_grcm39.png'))
plot(prob_cor$sim)
wh = which(prob_cor$sim < 0.5)
text(wh, prob_cor$sim[wh], labels = prob_cor$id[wh])
dev.off()

# Compare SODO188 & SODO189.
mixup_comp = cos_sim(gm_probs['SODO188',,], gbrs_probs['SODO189',,]) 

# There are two samples that don't match.
# SODO188 matches SODO189.
print(mixup_comp)

# Sample 188 and 189 are swapped. 

# Fix the two mixed up samples.
# I'm going to change the RNAseq files and leave genoprobs and phenotypes 
# the same.

# Read in the counts data and write out a corrected counts file.
file_table = data.frame(input  = c('sodo_gene_counts_grcm39.csv', 
                                   'sodo_gene_counts_norm_grcm39.csv', 
                                   'sodo_transcript_counts_grcm39.csv'),
                        output = c('sodo_gene_counts_grcm39_corrected.csv',
                                   'sodo_gene_counts_norm_grcm39_corrected.csv',
                                   'sodo_transcript_counts_grcm39_corrected.csv'))

for(i in 1:nrow(file_table)) {

  counts = read.csv(file.path(results_dir, file_table$input[i]))

  tmp                = counts[,'SODO188']
  counts[,'SODO188'] = counts[,'SODO189']
  counts[,'SODO189'] = tmp

  print(nrow(counts))

  write.csv(counts, file = file.path(results_dir, file_table$output[i]),
            row.names = FALSE, quote = FALSE)

} # for(i)

# Read in the allele counts data and write out a corrected counts file.
file_table = data.frame(input  = c('sodo_gene_allele_counts_grcm39.rds', 
                                   'sodo_transcript_allele_counts_grcm39.rds'),
                        output = c('sodo_gene_allele_counts_grcm39_corrected.rds', 
                                   'sodo_transcript_allele_counts_grcm39_corrected.rds'))

for(i in 1:nrow(file_table)) {

  counts = readRDS(file.path(results_dir, file_table$input[i]))

  tmp                 = counts['SODO188',,]
  counts['SODO188',,] = counts['SODO189',,]
  counts['SODO189',,] = tmp

  saveRDS(counts, file = file.path(results_dir, file_table$output[i]))

} # for(i)



