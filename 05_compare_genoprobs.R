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

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation/'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results')
out_dir     = file.path(base_dir, 'results')

# Gigamuga marker file (GRCm38)
gm_marker_file = '/projects/omics_share/mouse/GRCm38/supporting_files/muga_annotation/gigamuga/gm_uwisc_v1.csv'

# 69K grid marker file.
gbrs_marker_file = file.path(results_dir, 'ref.genome_grid.69k.noYnoMT_KBEdit.txt')

##### MAIN #####

# Read in the Gigamuga markers.
gm_markers = read.csv(gm_marker_file)
gm_markers$pos = gm_markers$bp_mm10 * 1e-6
gm_markers = subset(gm_markers, gm_markers$chr %in% c(1:19, 'X'))
gm_map = qtl2convert::map_df_to_list(gm_markers, pos_column = 'pos')

# Read in the 69K grid markers.
gbrs_markers = read.delim(gbrs_marker_file)
gbrs_map = qtl2convert::map_df_to_list(gbrs_markers, pos_column = 'pos')

# Read in the Gigamuga genoprobs and interpolate to the 69K grid.
gm_probs = readRDS(file.path(data_dir, 'bolcun-filas1__GigaMUGA_genoprobs_8state.rds'))
for(chr in seq_along(gm_probs)) {

  gm_map[[chr]] = gm_map[[chr]][dimnames(gm_probs[[chr]])[[3]]]

} # for(chr)
gm_probs = interpolate_genoprobs(gm_probs, gm_map, gbrs_map)
gm_probs = qtl2convert::probs_qtl2_to_doqtl(gm_probs)
rownames(gm_probs) = gsub('^Jackson_Lab_Filas_MURGIGV01_20220929_|_T_[A-H][0-9]+$', '',
                          rownames(gm_probs))
rownames(gm_probs) = gsub('_[A-H][0-9]+$', '',      rownames(gm_probs))
rownames(gm_probs) = gsub('^DO_QTL_T ',    'SO', rownames(gm_probs))

# Read in the 69K genoprobs.
gbrs_probs = readRDS(file.path(results_dir, 'sodo_gbrs_genoprobs_69K.rds'))
rownames(gbrs_probs) = sub('DO(_)?', '', rownames(gbrs_probs)) 

# Synch up sample IDs between the genoprobs.
samples = intersect(rownames(gm_probs), rownames(gbrs_probs)) 

gm_probs   = gm_probs[samples,,]
gbrs_probs = gbrs_probs[samples,,]

# Verify that the Gigamuga and GBRS genoprobs have the same structure.
stopifnot(rownames(gm_probs) == rownames(gbrs_probs))
stopifnot(colnames(gm_probs) == colnames(gbrs_probs))
stopifnot(dimnames(gm_probs)[[3]] == dimnames(gbrs_probs)[[3]])

# Reformat the genoprobs to be two dimensional.
# Each sample will be one column.
n_samples  = nrow(gm_probs)
n_founders = ncol(gm_probs)
n_markers  = dim(gm_probs)[[3]]
gm_vec = matrix(0, nrow = n_founders * n_markers, ncol = n_samples,
                dimnames = list(NULL, samples))
gbrs_vec = matrix(0, nrow = n_founders * n_markers, ncol = n_samples,
                  dimnames = list(NULL, samples))
for(i in 1:nrow(gm_probs)) {

  gm_vec[,i]   = as.vector(gm_probs[i,,])
  gbrs_vec[,i] = as.vector(gbrs_probs[i,,])

} # for(i)

# Get the correlation of the genoprobs with each other.
prob_cor = cor(gm_vec, gbrs_vec)

png(file.path(results_dir, 'genoprobs_cor.png'))
x = diag(prob_cor)
plot(x)
wh = which(x < 0.5)
text(wh, x[wh], labels = rownames(prob_cor)[wh])
dev.off()

# There are 5 samples that don't match.
png(file.path(results_dir, 'genoprobs_cor_mismatches.png'))
image(1:length(wh), 1:length(wh), prob_cor[wh,wh], axes = F, ann = F)
axis(side = 1, at = 1:length(wh), labels = rownames(prob_cor)[wh])
axis(side = 2, at = 1:length(wh), labels = rownames(prob_cor)[wh])
dev.off()

# Sample 188 and 189 are swapped. I'm not sure about the other three.
# They don't match any other samples.

# Fix the two mixed up samples and remove the three that don't match
# any other samples.

# Read in the original Gigamuga genoprobs.
gm_probs = readRDS(file.path(data_dir, 'bolcun-filas1__GigaMUGA_genoprobs_8state.rds'))
for(chr in seq_along(gm_probs)) {

  gm_map[[chr]] = gm_map[[chr]][dimnames(gm_probs[[chr]])[[3]]]

} # for(chr)

# Fix the sample mismatches. 
for(chr in seq_along(gm_probs)) {

  rownames(gm_probs[[chr]]) = gsub('^Jackson_Lab_Filas_MURGIGV01_20220929_|_T_[A-H][0-9]+$', '',
                            rownames(gm_probs[[chr]]))
  rownames(gm_probs[[chr]]) = gsub('_[A-H][0-9]+$', '',      rownames(gm_probs[[chr]]))
  rownames(gm_probs[[chr]]) = gsub('^DO_QTL_T ',    'SO', rownames(gm_probs[[chr]]))

  tmp                        = gm_probs[[chr]]['SO188',,]
  gm_probs[[chr]]['SO188',,] = gm_probs[[chr]]['SO189',,]
  gm_probs[[chr]]['SO189',,] = tmp
 
  remove = which(rownames(gm_probs[[chr]]) %in% c('SO289', 'SO290', 'SO291'))
  gm_probs[[chr]] = gm_probs[[chr]][-remove,,]

} # for(chr)

# Write out the corrected genoprobs.
saveRDS(gm_probs, 
        file = file.path(results_dir, 'bolcun-filas1__GigaMUGA_genoprobs_8state_corrected.rds'))

# Read in the counts data and write out a corrected counts file.
file_table = data.frame(input  = c('sodo_gene_counts.csv', 
                                   'sodo_gene_counts_norm.csv', 
                                   'sodo_transcript_counts.csv'),
                        output = c('sodo_gene_counts_corrected.csv',
                                   'sodo_gene_counts_norm_corrected.csv',
                                   'sodo_transcript_counts_corrected.csv'))

for(i in 1:nrow(file_table)) {

  counts = read.csv(file.path(results_dir, file_table$input[i]))

  tmp                 = counts[,'SODO188']
  counts[,'SODO188'] = counts[,'SODO189']
  counts[,'SODO189'] = tmp

  remove = which(colnames(counts) %in% c('SODO289', 'SODO290', 'SODO291'))
  counts = counts[,-remove]

  write.csv(counts, file = file.path(results_dir, file_table$output[i]),
            row.names = FALSE, quote = FALSE)

} # for(i)

# Read in the allele counts data and write out a corrected counts file.
file_table = data.frame(input  = c('sodo_gene_allele_counts.rds', 
                                   'sodo_transcript_allele_counts.rds'),
                        output = c('sodo_gene_allele_counts_corrected.rds', 
                                   'sodo_transcript_allele_counts_corrected.rds'))

for(i in 1:nrow(file_table)) {

  counts = readRDS(file.path(results_dir, file_table$input[i]))

  tmp                 = counts['SODO188',,]
  counts['SODO188',,] = counts['SODO189',,]
  counts['SODO189',,] = tmp

  remove = which(colnames(counts) %in% c('SODO289', 'SODO290', 'SODO291'))
  counts = counts[-remove,,]

  saveRDS(counts, file = file.path(results_dir, file_table$output[i]))

} # for(i)



