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

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)

genome          = args[1]
ensembl_version = args[2]

print(paste(genome, ensembl_version))

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results', paste0('eqtl_', genome))

# This assumes that the 06_gather_qtl_data.R script has been run and that
# 'superovulation_qtl2.Rdata' exists.
#data_file = file.path(data_dir, paste0('superovulation_qtl2_', genome, '_137K.Rdata'))
pheno_file = file.path(data_dir, 'superovulation_all_pheno.csv')

# All genoprobs file.
all_probs_file = file.path(data_dir, paste0('superovulation_all_probs_', genome, '.rds'))

# Gene annotation file.
annot_file   = file.path(results_dir, paste0('../gene_annotation_ens', ensembl_version, '.csv'))

# All expression file.
all_expr_file = file.path(results_dir, paste0('../soall_gene_counts_norm.csv'))

# Gigamuga markers file.
marker_file = 'https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/gm_uwisc_v2.csv'

##### MAIN #####

# Read in phenotypes.
pheno = readr::read_csv(pheno_file, show_col_types = FALSE)
pheno = as.data.frame(pheno)
rownames(pheno) = pheno$id 
pheno$ivf_date = format(pheno$ivf_date, format = '%Y-%m-%d')
pheno$ivf_date = factor(pheno$ivf_date) 

# Read in  allele probs.
probs = readRDS(all_probs_file)
probs_attr = attributes(probs)

# Read in markers.
markers     = readr::read_csv(marker_file, show_col_types = FALSE)
markers     = as.data.frame(markers)
markers$pos = markers$bp_mm10 * 1e-6

# Synch up samples between pheno & probs.
common_samples = intersect(pheno$id, rownames(probs[[1]]))
pheno = subset(pheno, id %in% common_samples)
probs = lapply(probs, function(z) { z[common_samples,,] })
stopifnot(pheno$id == rownames(probs[[1]]))

# Synch up markers between map and probs.
map = map_df_to_list(markers, pos_column = 'pos')
map = map[c(1:19, 'X')]

for(i in seq_along(probs)) {

  common_markers = intersect(dimnames(probs[[i]])[[3]], names(map[[i]]))
  probs[[i]]     = probs[[i]][,,common_markers] 
  map[[i]]       = map[[i]][common_markers]
  stopifnot(dimnames(probs[[i]])[[3]] == names(map[[i]]))
  
} # for(i)

# Create kinship matrices.
attributes(probs) = probs_attr
K = calc_kinship(probs, type = 'loco')

# Read in gene expression & annotation.
annot = read.csv(annot_file)
expr  = read.csv(all_expr_file)

rownames(annot) = annot$gene_id
rownames(expr)  = expr$ensembl

# Synch up genes between expression & annotation.
common_genes = intersect(rownames(annot), rownames(expr))
annot = annot[common_genes,]
expr  = expr[common_genes,]
stopifnot(rownames(annot) == rownames(expr)) 

# Convert expression to matrix and rankZ transform. 
expr = as.matrix(expr[,-1])
expr = t(expr)
expr = apply(expr, 2, rankZ)

# Only keep genes with genomic positions (i.e. not scaffolds).
# We are mapping local eQTL, so Y and MT don't make sense here.
annot = subset(annot, chr %in% c(1:19, 'X'))
expr  = expr[,rownames(annot)]

# Add covariates.
addcovar = model.matrix(~ivf_date, data = pheno)[,-1,drop = FALSE]

eqtl = data.frame(gene_id    = annot$gene_id,
                  symbol     = annot$gene_name,
                  gene_chr   = annot$chr,
                  gene_start = annot$start,
                  gene_end   = annot$end,
                  strand     = annot$strand,
                  marker     = '',
                  lod        = 0,
                  matrix(0, nrow = nrow(annot), ncol = 8, dimnames = list(NULL, LETTERS[1:8])),
                  matrix(0, nrow = nrow(annot), ncol = 8, dimnames = list(NULL, paste0(LETTERS[1:8], '_se'))))

# Loop through each gene, find the marker nearest the start, and map.
for(i in 1:nrow(annot)) {

  if(i %% 100 == 0) print(i)
  
  mkr = find_marker(map = map, chr = annot$chr[i], pos = annot$start[i] * 1e-6)
  pr  = probs[[annot$chr[i]]][,,mkr]
  
  mod = fit1(genoprobs = pr, pheno = expr[,i,drop = FALSE],
             kinship = K[[annot$chr[i]]], addcovar = addcovar,
             blup = FALSE)

  eqtl$lod[i]    = mod$lod
  eqtl$marker[i] = mkr
  
  mod = fit1(genoprobs = pr, pheno = expr[,i,drop = FALSE],
             kinship = K[[annot$chr[i]]], addcovar = addcovar,
             blup = TRUE)
  
  eqtl[i, LETTERS[1:8]] = mod$coef[1:8]
  eqtl[i, paste0(LETTERS[1:8], '_se')] = mod$SE[1:8]
  
} # for(i)

print('Writing results file')

saveRDS(eqtl, file = file.path(results_dir, paste0('local_eqtl_all_', genome, '.rds')))


#############################
# Subset to only use OD mice.

pheno = subset(pheno, grepl('^SODO', pheno$id))
probs = probs[pheno$id,]
for(i in seq_along(K)) {

  K[[i]] = K[[i]][rownames(probs[[1]]), rownames(probs[[1]])]

} # for(i)

# Add covariates.
addcovar = model.matrix(~ivf_date, data = pheno)[,-1,drop = FALSE]

eqtl = data.frame(gene_id    = annot$gene_id,
                  symbol     = annot$gene_name,
                  gene_chr   = annot$chr,
                  gene_start = annot$start,
                  gene_end   = annot$end,
                  strand     = annot$strand,
                  marker     = '',
                  lod        = 0,
                  matrix(0, nrow = nrow(annot), ncol = 8, dimnames = list(NULL, LETTERS[1:8])),
                  matrix(0, nrow = nrow(annot), ncol = 8, dimnames = list(NULL, paste0(LETTERS[1:8], '_se'))))

# Loop through each gene, find the marker nearest the start, and map.
for(i in 1:nrow(annot)) {

  if(i %% 100 == 0) print(i)
  
  mkr = find_marker(map = map, chr = annot$chr[i], pos = annot$start[i] * 1e-6)
  pr  = probs[[annot$chr[i]]][,,mkr]
  
  mod = fit1(genoprobs = pr, pheno = expr[,i,drop = FALSE],
             kinship = K[[annot$chr[i]]], addcovar = addcovar,
             blup = FALSE)

  eqtl$lod[i]    = mod$lod
  eqtl$marker[i] = mkr
  
  mod = fit1(genoprobs = pr, pheno = expr[,i,drop = FALSE],
             kinship = K[[annot$chr[i]]], addcovar = addcovar,
             blup = TRUE)
  
  eqtl[i, LETTERS[1:8]] = mod$coef[1:8]
  eqtl[i, paste0(LETTERS[1:8], '_se')] = mod$SE[1:8]
  
} # for(i)

print('Writing results file')

saveRDS(eqtl, file = file.path(results_dir, paste0('local_eqtl_do_', genome, '.rds')))



