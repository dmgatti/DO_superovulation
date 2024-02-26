options(stringsAsFactors = FALSE)

library(readxl)
library(qtl2)
library(qtl2convert)

base_dir = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir = file.path(base_dir, 'data')
fig_dir  = file.path(base_dir, 'figures')
results_dir = file.path(base_dir, 'results', 'qtl2')
marker_dir = '/projects/omics_share/mouse/GRCm38/supporting_files/muga_annotation/gigamuga'

qtl2_dir   = '/projects/compsci/vmp/USERS/dgatti/data/qtl2/'
ccsnp_file = file.path(qtl2_dir, 'cc_variants_v3.sqlite')
gene_file  = file.path(qtl2_dir, 'mouse_genes_mgi_v8.sqlite')

snp_fxn  = create_variant_query_func(ccsnp_file)
gene_fxn = create_gene_query_func(gene_file)

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

pheno = read_xlsx(file.path(data_dir, 'JDO_Superovulation_collect_1_2_3.xlsx'), 
                  sheet = 'Data')
covar = read.csv(file.path(data_dir, 'do_superovulation_animal_info.csv'))
probs = readRDS(file.path(data_dir, 'bolcun-filas1__GigaMUGA_genoprobs_8state.rds'))
markers = read.csv(file.path(marker_dir, 'gm_uwisc_v1.csv'))
markers$bp_mm10 = markers$bp_mm10 * 1e-6
markers = subset(markers, chr %in% c(1:19, 'X'))
map = qtl2convert::map_df_to_list(markers, pos_column = 'bp_mm10')

colnames(covar) = c('mouse', 'sex', 'gen')
covar$mouse = gsub('^SO|_T$' , '', covar$mouse)
covar$mouse = paste0('SODO', covar$mouse)
rownames(covar) = covar$mouse

pheno = pheno[,c('Female number', 'FreshIVF_IVFDISHES::IVF Date', 'Weight', 
                 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 'DAY1NumbFrag', 
                 'Total oocytes ovulated', 'DAY2Numb2cells', 'Fertilized (%)', 'DAY2NumbDead', 'DAY2NumbFrag')]
colnames(pheno)[1:2] = c('mouse', 'date')
colnames(pheno)[colnames(pheno) == 'Fertilized (%)']  = 'Fertilized'
colnames(pheno)[colnames(pheno) == 'Total oocytes ovulated']  = 'TotalOocytes'
pheno$Weight = as.numeric(sub('g$', '', pheno$Weight))
pheno = as.data.frame(pheno)
pheno$mouse = paste0('SODO', pheno$mouse) 
rownames(pheno) = pheno$mouse

covar = merge(covar, pheno[,c('mouse', 'date', 'TotalOocytes', 'NumbClutches', 'DAY1Numb1cells', 'DAY2Numb2cells')], by = 'mouse', sort = FALSE)
rownames(covar) = covar$mouse
covar$date = factor(covar$date)
pheno = as.matrix(pheno[,-(1:2)])

for(chr in seq_along(probs)) {

  rownames(probs[[chr]]) = gsub('^Jackson_Lab_Filas_MURGIGV01_20220929_SO|_T_[A-H][0-9]+$', '', rownames(probs[[chr]]))
  rownames(probs[[chr]]) = gsub('^Jackson_Lab_Filas_MURGIGV01_20220929_DO_QTL_T |_[A-H][0-9]+$', '', rownames(probs[[chr]]))
  rownames(probs[[chr]]) = paste0('SODO', rownames(probs[[chr]]))
  probs[[chr]] = probs[[chr]][!duplicated(rownames(probs[[chr]])),,]
  map[[chr]]   = map[[chr]][dimnames(probs[[chr]])[[3]]]

} # for(chr)

# Synch up samples.
samples = sort(intersect(rownames(pheno), rownames(probs[[1]])))
pheno = pheno[samples,]
covar = covar[samples,]
probs = probs[samples,]

K = calc_kinship(probs, type = 'loco', cores = 4)

# Save Rdata file with all QTL mapping data.
save(pheno, covar, probs, K, map,
     file = file.path(data_dir, 'superovulation_qtl2.Rdata'))

# Weight
addcovar = model.matrix(~date, data = covar)[,-1,drop = FALSE]
ph = pheno[,'Weight', drop = FALSE]
ph[,'Weight'] = log(ph[,'Weight'])
lod = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
plot(lod, map, main = 'Body Weight, covar: date')
ph[,1] = rankZ(ph[,1])
lod_rz = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
               kinship = K, cores = 4)
plot(lod, map, main = 'Body Weight, rankZ, covar: date')


# NumbClutches
pheno_name = 'NumbClutches' 
addcovar = model.matrix(~date, data = covar)[,-1,drop = FALSE]
ph = pheno[,pheno_name, drop = FALSE]
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
plot(qtl, map, main = 'Num Clutches, covar: date')
lod = cbind(lod, qtl)

peaks = find_peaks(lod, map, threshold = 8, prob = 0.95)

chr = '5'
feff = scan1blup(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
                 kinship = K[[chr]], cores = 4)
curr_peak = peaks[peaks$chr == chr,]
saveRDS(feff, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '.rds')))
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '.png'))),
    width = 1000, height = 800, res = 128)
plot_coefCC(feff, map, scan1_output = lod[,pheno_name, drop = FALSE],
            main = pheno_name, xlim = curr_peak$pos + c(-5, 5),
            ylim = c(-0.7, 0.2), legend = 'bottomleft')
dev.off()

chr = '15'
feff = scan1blup(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
                 kinship = K[[chr]], cores = 4)
curr_peak = peaks[peaks$chr == chr,]
saveRDS(feff, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '.rds')))
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '.png'))))
plot_coefCC(feff, map, scan1_output = lod[,pheno_name, drop = FALSE],
            main = pheno_name, xlim = c(100, 104))
dev.off()


# DAY1Numb1cells
pheno_name = 'DAY1Numb1cells' 
addcovar = model.matrix(~date + TotalOocytes, data = covar)[,-1,drop = FALSE]
ph = pheno[,pheno_name, drop = FALSE]
ph[,pheno_name] = log1p(ph[,pheno_name])
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
plot(qtl, map, main = 'Day 1 Num 1-cells, covar: date + TotalOocytes')
lod = cbind(lod, qtl)

ph[,1] = rankZ(ph[,1])
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
plot(qtl, map, main = 'Day 1 Num 1-cells, rankZ, covar: date + TotalOocytes')

peaks = find_peaks(lod, map, threshold = 7, prob = 0.95)

chr = '2'
curr_peak = peaks[peaks$lodcolumn == pheno_name & peaks$chr == chr,]
feff = scan1blup(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
                 kinship = K[[chr]], cores = 4)
saveRDS(feff, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '.rds')))
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '.png'))),
    width = 1000, height = 800, res = 128)
plot_coefCC(feff, map, scan1_output = lod[,pheno_name, drop = FALSE],
            main = pheno_name, xlim = curr_peak$pos + c(-5, 5),
            legend = 'bottomleft')
dev.off()

start = curr_peak$pos - 2
end   = curr_peak$pos + 2
assoc = scan1snps(genoprobs = probs[,chr], map = map[chr], pheno = ph, 
                  addcovar = addcovar, kinship = K[[chr]], query_func = snp_fxn, 
                  chr = chr, start = start, end = end, keep_all_snps = TRUE)

genes = gene_fxn(chr = chr, start = start, end = end)
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '_assoc.png'))),
    width = 1000, height = 800, res = 128)
plot_snpasso(scan1output = assoc$lod, snpinfo = assoc$snpinfo, genes = genes, 
             drop_hilit = 1, sdp_panel = TRUE, panel_prop = c(0.2, 0.2, 0.6),
             colors = 'black')
dev.off()

snp_map = qtl2:::snpinfo_to_map(assoc$snpinfo)
tmp     = qtl2:::expand_snp_results(snp_results = assoc$lod, snp_map, assoc$snpinfo)
tmp     = data.frame(assoc$snpinfo, lod = tmp$lod)
write.csv(tmp, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '_assoc.csv')), 
          row.names = FALSE, quote = FALSE)

chr = '5'
curr_peak = peaks[peaks$lodcolumn == pheno_name & peaks$chr == chr,]
feff = scan1blup(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
                 kinship = K[[chr]], cores = 4)
feff[abs(feff) > 1] = NA
saveRDS(feff, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '.rds')))
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '.png'))))
plot_coefCC(feff, map, scan1_output = lod[,pheno_name, drop = FALSE],
            main = pheno_name, xlim = curr_peak$pos + c(-3, 8), ylim = c(-0.6, 0.4),
            legend = 'bottomleft')
dev.off()

start = curr_peak$pos + 2
end   = curr_peak$pos + 6
assoc = scan1snps(genoprobs = probs[,chr], map = map[chr], pheno = ph, 
                  addcovar = addcovar, kinship = K[[chr]], query_func = snp_fxn, 
                  chr = chr, start = start, end = end, keep_all_snps = TRUE)

genes = gene_fxn(chr = chr, start = start, end = end)
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '_assoc.png'))),
    width = 1000, height = 800, res = 128)
plot_snpasso(scan1output = assoc$lod, snpinfo = assoc$snpinfo, genes = genes, 
             drop_hilit = 1, sdp_panel = TRUE, panel_prop = c(0.2, 0.2, 0.6),
             colors = 'black')
dev.off()

snp_map = qtl2:::snpinfo_to_map(assoc$snpinfo)
tmp     = qtl2:::expand_snp_results(snp_results = assoc$lod, snp_map, assoc$snpinfo)
tmp     = data.frame(assoc$snpinfo, lod = tmp$lod)
write.csv(tmp, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '_assoc.csv')), 
          row.names = FALSE, quote = FALSE)

chr = '19'
curr_peak = peaks[peaks$lodcolumn == pheno_name & peaks$chr == chr,]
feff = scan1blup(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
                 kinship = K[[chr]], cores = 4)
saveRDS(feff, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '.rds')))
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '.png'))))
plot_coefCC(feff, map, scan1_output = lod[,pheno_name, drop = FALSE],
            main = pheno_name, xlim = curr_peak$pos + c(-3, 5),
            legend = 'bottomright')
dev.off()

start = curr_peak$pos - 0.5
end   = curr_peak$pos + 1.5
assoc = scan1snps(genoprobs = probs[,chr], map = map[chr], pheno = ph, 
                  addcovar = addcovar, kinship = K[[chr]], query_func = snp_fxn, 
                  chr = chr, start = start, end = end, keep_all_snps = TRUE)

genes = gene_fxn(chr = chr, start = start, end = end)
png(file.path(fig_dir, paste0(paste0(pheno_name, '_chr', chr, '_assoc.png'))),
    width = 1000, height = 800, res = 128)
plot_snpasso(scan1output = assoc$lod, snpinfo = assoc$snpinfo, genes = genes, 
             drop_hilit = 1, sdp_panel = TRUE, panel_prop = c(0.2, 0.2, 0.6),
             colors = 'black')
dev.off()

snp_map = qtl2:::snpinfo_to_map(assoc$snpinfo)
tmp     = qtl2:::expand_snp_results(snp_results = assoc$lod, snp_map, assoc$snpinfo)
tmp     = data.frame(assoc$snpinfo, lod = tmp$lod)
write.csv(tmp, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '_assoc.csv')), 
          row.names = FALSE, quote = FALSE)

# Map DAY1Numb1cells with NumClutches as a covariate since they seem to share
# a QTL.
addcovar = model.matrix(~date + TotalOocytes + NumbClutches, data = covar)[,-1,drop = FALSE]
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
plot(qtl, map)

chr = '13'
peaks = find_peaks(qtl, map, threshold = 7) 
curr_peak = peaks[peaks$lodcolumn == pheno_name & peaks$chr == chr,]
feff = scan1blup(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
                 kinship = K[[chr]], cores = 4)
plot_coefCC(feff, map, scan1_output = qtl, main = pheno_name, 
            legend = 'bottomright')


# DAY1NumbDead
pheno_name = 'DAY1NumbDead' 
addcovar = model.matrix(~date + TotalOocytes, data = covar)[,-1,drop = FALSE]
ph = pheno[,pheno_name, drop = FALSE]
ph[,pheno_name] = log1p(ph[,pheno_name])
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
lod = cbind(lod, qtl)

chr = '1'
feff = scan1blup(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
                 kinship = K[[chr]], cores = 4)
saveRDS(feff, file = file.path(results_dir, paste0(pheno_name, '_chr', chr, '.rds')))


# DAY1NumbFrag
addcovar = model.matrix(~date + TotalOocytes, data = covar)[,-1,drop = FALSE]
ph = pheno[,'DAY1NumbFrag', drop = FALSE]
ph[,'DAY1NumbFrag'] = log1p(ph[,'DAY1NumbFrag'])
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
lod = cbind(lod, qtl)

# TotalOocytes
addcovar = model.matrix(~date + NumbClutches, data = covar)[,-1,drop = FALSE]
ph = pheno[,'TotalOocytes', drop = FALSE]
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
lod = cbind(lod, qtl)

# DAY2Numb2cells
addcovar = model.matrix(~date + DAY1Numb1cells, data = covar)[,-1,drop = FALSE]
ph = pheno[,'DAY2Numb2cells', drop = FALSE]
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
lod = cbind(lod, qtl)

# Fertilized
addcovar = model.matrix(~date , data = covar)[,-1,drop = FALSE]
ph = pheno[,'Fertilized', drop = FALSE]
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
lod = cbind(lod, qtl)

# DAY2NumbDead
addcovar = model.matrix(~date + DAY2Numb2cells, data = covar)[,-1,drop = FALSE]
ph = pheno[,'DAY2NumbDead', drop = FALSE]
ph[,'DAY2NumbDead'] = log1p(ph[,'DAY2NumbDead'])
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
lod = cbind(lod, qtl)

# DAY2NumbFrag
addcovar = model.matrix(~date + DAY2Numb2cells, data = covar)[,-1,drop = FALSE]
ph = pheno[,'DAY2NumbFrag', drop = FALSE]
ph[,'DAY2NumbFrag'] = log1p(ph[,'DAY2NumbFrag'])
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)
lod = cbind(lod, qtl)


saveRDS(lod, file = file.path(base_dir, 'results', 'qtl2', 'lod.rds'))

# Plot QTL peaks
for(i in 1:ncol(lod)) {
  png(file.path(fig_dir, paste0(colnames(lod)[i], '_lod.png')),
      width = 1000, height = 800, res = 128)
  plot_scan1(lod, map, lodcolumn = i, main = colnames(lod)[i])
  dev.off()
} # for(i)

# Find Peaks
peaks = find_peaks(lod, map, threshold = 7)
peaks = subset(peaks, (lodcolumn == 'NumbClutches' & lod > 10.0) |
               lodcolumn %in% c('DAY1Numb1cells', 'DAY1NumbDead', 'DAY2Numb2cells', 
               'Fertilized'))
write.csv(peaks, file = file.path(results_dir, 'sodo_qlt2_peaks.csv'))

################################################################################
# Mediation analysis using RNAseq data.

options(stringsAsFactors = FALSE)

library(AnnotationHub)
library(DESeq2)
library(readxl)
library(qtl2)
library(qtl2convert)

base_dir = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir = file.path(base_dir, 'data')
fig_dir  = file.path(base_dir, 'figures')
results_dir = file.path(base_dir, 'results', 'qtl2')
marker_dir = '/projects/omics_share/mouse/GRCm38/supporting_files/muga_annotation/gigamuga'

# Load in qtl2 data.
load(file = file.path(data_dir, 'superovulation_qtl2.Rdata'))

annot = read.csv(file.path(results_dir, '..', 'gene_annotation_ens98.csv'))
annot = subset(annot, !is.na(annot$gene_id))
rownames(annot) = annot$gene_id

counts = read.csv(file.path(results_dir, '..', 'sodo_gene_counts_correted.csv'))
colnames(counts) = sub('_', '', colnames(counts))
rownames(counts) = counts$ensembl
counts = as.matrix(round(counts[,-1]))
counts = counts[,colMeans(counts) > 0]

# NOTE: There are 708 genes in counts that are NOT in Ensembl 98 annotation.
sum(!rownames(counts) %in% annot$gene_id)

# Synch up genes.
common_genes = intersect(rownames(annot), rownames(counts))
counts = counts[common_genes,]
annot = annot[common_genes,]
stopifnot(rownames(annot) == rownames(counts))

# Synch up samples.
samples = intersect(rownames(pheno), colnames(counts))

pheno    = pheno[samples,]
covar    = covar[samples,]
probs    = probs[samples,]
K        = lapply(K, function(z) { z[samples, samples] })
counts   = round(counts[,samples])
metadata = covar[samples,'date', drop = FALSE]

stopifnot(colnames(counts) == rownames(metadata))

# Filter low counts genes. Keep genes expressed in at least half of the samples.
rm_gt0 = rowMeans(counts > 0)
dim(counts)
counts = counts[rm_gt0 > 0.5,]
dim(counts)
annot = annot[rownames(counts),]

dds = DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ date)
dds = DESeq(dds)
vst = assay(vst(dds))

write.csv(vst, file = file.path(results_dir, '..', 'sodo_gene_counts_correted_norm.csv'),
          row.names = FALSE, quote = FALSE)

# Mediate on Chr 13 for DAY1DAY1Numb1cells.
pheno_name = 'DAY1Numb1cells' 
chr = '13'
addcovar = model.matrix(~date + TotalOocytes, data = covar)[,-1,drop = FALSE]
ph = pheno[,pheno_name, drop = FALSE]
ph[,pheno_name] = log1p(ph[,pheno_name])
ph[,pheno_name] = rankZ(ph[,pheno_name])
lod = scan1(genoprobs = probs[,chr], pheno = ph, addcovar = addcovar, 
            kinship = K[[chr]], cores = 4)

max_pos  = max(lod, map)
pr       = pull_genoprobpos(probs, map, max_pos$chr, max_pos$pos)
base_mod = fit1(genoprobs = pr, pheno = ph, kinship = K[[chr]], addcovar = addcovar)
base_lod = base_mod$lod

# Get Chr 5 genes.
chr5_annot = annot[annot$chr == chr,]
chr5_vst   = vst[rownames(chr5_annot),]

med_out = data.frame(chr      = chr5_annot$chr,
                     pos      = chr5_annot$start * 1e-6,
                     gene_id  = chr5_annot$gene_id,
                     symbol   = chr5_annot$symbol,
                     base_lod = base_lod,
                     med_lod  = 0)

for(i in 1:nrow(chr5_vst)) {

  curr_covar = cbind(addcovar, chr5_vst[i,])
  med_mod    = fit1(genoprobs = pr, pheno = ph, kinship = K[[chr]], 
                    addcovar = curr_covar)
  med_out$med_lod[i] = med_mod$lod
  
} # for(i)

med_out$lod_drop = med_out$med_lod - med_out$base_lod

plot(lod_drop ~ pos, data = med_out, main = pheno_name)
abline(v = max_pos$pos, col = 2)

# Mediate on Chr 9 for DAY2Numb2cells
pheno_name = 'DAY2Numb2cells' 
chr = '9'
addcovar = model.matrix(~date + DAY1Numb1cells, data = covar)[,-1,drop = FALSE]
ph = pheno[,'DAY2Numb2cells', drop = FALSE]
qtl = scan1(genoprobs = probs, pheno = ph, addcovar = addcovar, 
            kinship = K, cores = 4)

max_pos  = max(qtl, map)
pr       = pull_genoprobpos(probs, map, max_pos$chr, max_pos$pos)
base_mod = fit1(genoprobs = pr, pheno = ph, kinship = K[[chr]], addcovar = addcovar)
base_lod = base_mod$lod

# Get Chr 9 genes.
chr9_annot = annot[annot$chr == chr,]
chr9_vst   = vst[rownames(chr9_annot),]

med_out = data.frame(chr      = chr9_annot$chr,
                     pos      = chr9_annot$start * 1e-6,
                     gene_id  = chr9_annot$gene_id,
                     symbol   = chr9_annot$symbol,
                     base_lod = base_lod,
                     med_lod  = 0)

for(i in 1:nrow(chr9_vst)) {

  curr_covar = cbind(addcovar, chr9_vst[i,])
  med_mod    = fit1(genoprobs = pr, pheno = ph, kinship = K[[chr]], 
                    addcovar = curr_covar)
  med_out$med_lod[i] = med_mod$lod
  
} # for(i)

med_out$lod_drop = med_out$med_lod - med_out$base_lod

plot(lod_drop ~ pos, data = med_out, main = pheno_name)
abline(v = max_pos$pos, col = 2)

