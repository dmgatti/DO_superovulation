################################################################################
# Gather the data needed for QTL mapping.
# phenotypes
# covariates
# genoprobs
# kinship
# marker map
# normalized expression
# gene annotation
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-19
################################################################################

##### LIBRARIES #####

library(AnnotationHub)
library(ensembldb)
library(readxl)
library(qtl2)
library(qtl2convert)

##### VARIABLES #####

base_dir = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir = file.path(base_dir, 'data')
fig_dir  = file.path(base_dir, 'figures')
results_dir = file.path(base_dir, 'results')
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

##### MAIN #####

# Read in phenotypes.
pheno = read_xlsx(file.path(data_dir, 'JDO_Superovulation_collect_1_2_3.xlsx'), 
                  sheet = 'Data')
# Read in covariates.
covar = read.csv(file.path(data_dir, 'do_superovulation_animal_info.csv'))
# Read in genoprobs.
probs = readRDS(file.path(data_dir, 'bolcun-filas1__GigaMUGA_genoprobs_8state.rds'))
# Read in markers.
markers = read.csv(file.path(marker_dir, 'gm_uwisc_v1.csv'))
markers$bp_mm10 = markers$bp_mm10 * 1e-6
markers = subset(markers, chr %in% c(1:19, 'X'))
map = qtl2convert::map_df_to_list(markers, pos_column = 'bp_mm10')

# NOTE: I'm removing the '_' from samples IDs and adding 'SODO' to the
# beginning of each sample ID. 'SODO' stands for 'superovulation diversity outbred'.

# Wrangle covariates.
colnames(covar) = c('mouse', 'sex', 'gen')
covar$mouse = gsub('^SO|_T$' , '', covar$mouse)
covar$mouse = paste0('SODO', covar$mouse)
rownames(covar) = covar$mouse

# Wrangle phenotypes.
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

# Select covariates that we need for mapping phenotypes.
covar = merge(covar, pheno[,c('mouse', 'date', 'TotalOocytes', 'NumbClutches', 'DAY1Numb1cells', 'DAY2Numb2cells')], by = 'mouse', sort = FALSE)
rownames(covar) = covar$mouse
covar$date = factor(covar$date)
pheno = as.matrix(pheno[,-(1:2)])

# Subset probs to contain the correct markers and samples.
# Also remove markers with allele frequencies less than 1 for any one founder.
for(chr in seq_along(probs)) {

  rownames(probs[[chr]]) = gsub('^Jackson_Lab_Filas_MURGIGV01_20220929_SO|_T_[A-H][0-9]+$', '', rownames(probs[[chr]]))
  rownames(probs[[chr]]) = gsub('^Jackson_Lab_Filas_MURGIGV01_20220929_DO_QTL_T |_[A-H][0-9]+$', '', rownames(probs[[chr]]))
  rownames(probs[[chr]]) = paste0('SODO', rownames(probs[[chr]]))
  probs[[chr]] = probs[[chr]][!duplicated(rownames(probs[[chr]])),,]
  
  # Check allele frequencies.
  freq = apply(probs[[chr]], 2:3, sum)
  min_freq = apply(freq, 2, min)
  
  # Remove markers with minimum founder allele frequencies <= 1.
  probs[[chr]] = probs[[chr]][,,min_freq > 1]
  map[[chr]]   = map[[chr]][dimnames(probs[[chr]])[[3]]]
  
} # for(chr)

# Synch up samples.
samples = sort(intersect(rownames(pheno), rownames(probs[[1]])))
pheno = pheno[samples,]
covar = covar[samples,]
probs = probs[samples,]

# Create rankZ phenotypes.
pheno_rz = apply(pheno, 2, rankZ)

# Create LOCO kinship matrices.
K = calc_kinship(probs, type = 'loco', cores = 4)

# Read in the expression data.
counts    = read.csv(file.path(results_dir, 'sodo_gene_counts_corrected.csv'))
rownames(counts) = counts$ensembl
counts = as.matrix(counts[,-1])
norm_expr = read.csv(file.path(results_dir, 'sodo_gene_counts_norm_corrected.csv'))
rownames(norm_expr) = norm_expr$ensembl
norm_expr = as.matrix(norm_expr[,-1])
annot     = read.csv(file.path(results_dir, 'gene_annotation_ens98.csv'))
annot     = annot[!is.na(annot$gene_id),]
rownames(annot) = annot$gene_id

# Synch up genes.
common_genes = intersect(rownames(annot), rownames(counts))
counts    = counts[common_genes,]
norm_expr = norm_expr[common_genes,]
annot     = annot[common_genes,]


# Save Rdata file with all QTL mapping data.
save(pheno, pheno_rz, covar, probs, K, map, counts, norm_expr, annot,
     file = file.path(data_dir, 'superovulation_qtl2.Rdata'))



