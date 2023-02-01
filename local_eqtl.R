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
rna_dir     = file.path(base_dir, 'results', paste0('star_rsem_', genome))
results_dir = file.path(base_dir, 'results', paste0('eqtl_', genome))

# This assumes that the 06_gather_qtl_data.R script has been run and that
# 'superovulation_qtl2.Rdata' exists.
data_file = file.path(data_dir, paste0('superovulation_qtl2_', genome, '.Rdata'))

# Gene annotation file.
annot_file = file.path(rna_dir, paste0('gene_annotation_ens', ensembl_version, '.rds'))

# STAR expression file.
expr_file = file.path(rna_dir, paste0('star_rsem_gene_vst_counts_', genome, '.rds'))

##### MAIN #####

# Read in data.
load(data_file)

# Read in gene annotation.
annot = readRDS(annot_file)
expr  = readRDS(expr_file)
norm_expr = t(expr)
norm_expr = apply(norm_expr, 2, rankZ)

common_genes = intersect(rownames(annot), colnames(norm_expr))
annot = annot[common_genes,]
norm_expr = norm_expr[,common_genes]
stopifnot(rownames(annot) == colnames(norm_expr))

# Only keep genes with genomic positions (i.e. not scaffolds).
annot = subset(annot, chr %in% c(1:19, 'X'))
norm_expr = norm_expr[,rownames(annot)]

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

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
  pr  = pull_genoprobpos(genoprobs = probs, marker = mkr)
  
  mod = fit1(genoprobs = pr, pheno = norm_expr[,i,drop = FALSE],
             kinship = K[[annot$chr[i]]], addcovar = addcovar,
             blup = FALSE)

  eqtl$lod[i]    = mod$lod
  eqtl$marker[i] = mkr
  
  mod = fit1(genoprobs = pr, pheno = norm_expr[,i,drop = FALSE],
             kinship = K[[annot$chr[i]]], addcovar = addcovar,
             blup = TRUE)
  
  eqtl[i, LETTERS[1:8]] = mod$coef[1:8]
  eqtl[i, paste0(LETTERS[1:8], '_se')] = mod$SE[1:8]
  
} # for(i)

print('Writing results file')

saveRDS(eqtl, file = file.path(results_dir, paste0('local_eqtl_', genome, '.rds')))


