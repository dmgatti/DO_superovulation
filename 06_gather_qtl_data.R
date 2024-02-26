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

# Get command line arguments.
args = commandArgs(trailingOnly = TRUE)

# Genome build: either 'grcm38' or 'grcm39'.
genome = args[1]


# Base directories.
base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
fig_dir     = file.path(base_dir, 'figures')
results_dir = file.path(base_dir, 'results')
expr_dir    = results_dir

scratch_dir = '/flashscratch/dgatti'

# Minimum allele probability for each founder below which we will
# set the probs to zero. 
probs_threshold = 1e-6

# Phenotypes.
do_pheno_file = file.path(data_dir, 'JDO9376 all 350.xlsx')
cc_pheno_file = file.path(data_dir, 'CC lines export.xlsx')
founder_pheno_file = file.path(data_dir, 'Founder strains export.xlsx')

do_probs_file = file.path(data_dir, paste0('sodo_alleleprobs_', genome, '_137K.rds'))

# RNAseq: GBRS build 38
expr_file  = file.path(expr_dir, paste0('sodo_gene_counts_norm_', genome, '_corrected.csv'))
annot_file = dir(expr_dir, pattern = '^gene_annotation_', full.names = TRUE)

if(genome == 'grcm38') {

  # Ensembl 98.
  annot_file = annot_file[grep('98', annot_file)]

} else {

  # Ensembl 105.
  annot_file = annot_file[grep('105', annot_file)]

} # else

# CC genorpobs.
#cc_dir        = '/projects/compsci/vmp/USERS/dgatti/data/cc/genoprobs'
#cc_probs_file = file.path(cc_dir, 'cc_consensus_genoprobs_137K.rds')
cc_probs_url  = 'https://thejacksonlaboratory.box.com/shared/static/nofljklyykez5i4gu4g8wzyr1j8b9b68.rds'
cc_probs_file = file.path(scratch_dir, 'cc_alleleprobs.rds')

# Markers.
marker_dir    = '/projects/omics_share/mouse/GRCm38/supporting_files/muga_annotation/gigamuga'
marker_file   = file.path(marker_dir, 'gm_uwisc_v1.csv')

if(genome == 'grcm39') {

  marker_dir    = '/projects/omics_share/mouse/GRCm39/supporting_files/muga_annotation/gigamuga'
  marker_file   = file.path(marker_dir, 'gm_uwisc_v4.csv')  

} # if(genome == 'grcm39')

# qtl2 sqlite files.
qtl2_dir   = '/projects/omics_share/mouse/GRCm38/supporting_files/qtl2'
ccsnp_file = file.path(qtl2_dir, 'cc_variants_v3.sqlite')
gene_file  = file.path(qtl2_dir, 'mouse_genes_mgi_v8.sqlite')

if(genome =='grcm39') {

  qtl2_dir   = '/projects/omics_share/mouse/GRCm39/supporting_files/qtl2'
  ccsnp_file = file.path(qtl2_dir, 'fv.2021.snps.db3')
  gene_file  = file.path(qtl2_dir, 'fv.2021.snps.db3')

} # if(genome = ='grcm39')

# Create qtl2 gene and SNP query functions.
snp_fxn  = create_variant_query_func(ccsnp_file)
gene_fxn = create_gene_query_func(gene_file)

if(genome == 'grcm39') {

  snp_fxn  = create_variant_query_func(ccsnp_file, id_field = 'variants_id')
  gene_fxn = create_gene_query_func(gene_file, 
                                    chr_field   = 'chromosome', 
                                    name_field  = 'symbol',
                                    start_field = 'start_position', 
                                    stop_field  = 'end_position')

} # if(genome == 'grcm39')

##### FUNCTIONS #####

# RankZ.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

##### MAIN #####

dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)

# Read in phenotypes.
pheno    = read_xlsx(do_pheno_file, sheet = 'Sheet1')
cc_pheno = read_xlsx(cc_pheno_file, sheet = 'Sheet1')
cc_table = read.csv('/projects/compsci/vmp/USERS/dgatti/data/cc/cc_strain_jr_table.csv')
founder_pheno = read_xlsx(founder_pheno_file, sheet = 'Sheet1')
founder_table = read_xlsx(founder_pheno_file, sheet = 'Sample labels info') 

# Read in covariates.
covar = read.csv(file.path(data_dir, 'genotypes', 'do_superovulation_animal_info.csv'))

# Read in DO genoprobs.
do_probs = readRDS(do_probs_file)

# I had trouble with downloading files in the Singularity container. So I'm downloading
# the file and then reading it in locally.
if(!file.exists(cc_probs_file)) {
  download.file(url = cc_probs_url, destfile = cc_probs_file)
} # if(!file.exists(cc_probs_file))
cc_probs = readRDS(cc_probs_file)

# Read in markers.
markers = read.csv(marker_file)
markers$bp_mm10 = markers$bp_mm10 * 1e-6
markers = subset(markers, chr %in% c(1:19, 'X'))
map = NULL
if(genome == 'grcm38') {
  map = qtl2convert::map_df_to_list(markers, pos_column = 'bp_mm10')
} else {
  map = qtl2convert::map_df_to_list(markers, pos_column = 'bp_grcm39')
} # else

# NOTE: I'm removing the '_' from samples IDs and adding 'SODO' to the
# beginning of each sample ID. 'SODO' stands for 'superovulation diversity outbred'.

# Wrangle covariates.
colnames(covar) = c('mouse', 'sex', 'gen')
covar$sex = 'F'
covar$mouse = gsub('^SO|_T$' , '', covar$mouse)
covar  = subset(covar, !grepl('^DO_QTL_T ', covar$mouse))
covar$mouse = paste0('SODO', covar$mouse)
rownames(covar) = covar$mouse

# Wrangle phenotypes.

# DO phenotypes.
pheno = pheno[,c('Female number', 'FreshIVF_IVFDISHES::IVF Date', 'Weight', 
                 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 'DAY1NumbFrag', 
                 'Total oocytes ovulated', 'DAY2Numb2cells', 'Fertilized (%)', 
                 'DAY2NumbDead', 'DAY2NumbFrag')]
colnames(pheno) = c('mouse', 'date', 'Weight', 
                 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 'DAY1NumbFrag', 
                 'TotalOocytes', 'DAY2Numb2cells', 'Fertilized', 
                 'DAY2NumbDead', 'DAY2NumbFrag')
pheno$mouse[!grepl('^SODO', pheno$mouse)] = paste0('SODO', pheno$mouse[!grepl('^SODO', pheno$mouse)])
pheno$Weight = as.numeric(sub('g$', '', pheno$Weight))
pheno        = data.frame(mouse  = pheno$mouse,
                    strain = 'DO',
                    pheno[,c('date', 'Weight', 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead',
                             'DAY1NumbFrag', 'TotalOocytes', 'DAY2Numb2cells', 'Fertilized', 
                             'DAY2NumbDead', 'DAY2NumbFrag')])
rownames(pheno) = pheno$mouse
pheno = subset(pheno, NumbClutches > 0)

# CC phenotypes.
cc_pheno = cc_pheno[,c('Female number', 'Line number', 'FreshIVF_IVFDISHES::IVF Date', 
                       'Weight', 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 
                       'DAY1NumbFrag', 'Total oocytes ovulated', 'DAY2Numb2cells',
                       'Fertilized(%)', 'DAY2NumbDead', 'DAY2NumbFrag')]
colnames(cc_pheno) = c('mouse', 'jr', 'date', 
                       'Weight', 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 
                       'DAY1NumbFrag', 'TotalOocytes', 'DAY2Numb2cells',
                       'Fertilized', 'DAY2NumbDead', 'DAY2NumbFrag')
cc_pheno = merge(cc_pheno, cc_table, by = 'jr', all.x = TRUE, sort = FALSE)
cc_pheno = cc_pheno[,c('mouse', 'cc_strain', 'date', 
                       'Weight', 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 
                       'DAY1NumbFrag', 'TotalOocytes', 'DAY2Numb2cells',
                       'Fertilized', 'DAY2NumbDead', 'DAY2NumbFrag')]
colnames(cc_pheno)[colnames(cc_pheno) == 'cc_strain'] = 'strain'
cc_pheno$Weight = as.numeric(sub('g$', '', cc_pheno$Weight))
cc_pheno = subset(cc_pheno, NumbClutches > 0)
cc_pheno$date = as.character(cc_pheno$date)
# Condense each strain to strain means.
cc_pheno = data.frame(mouse  = unique(gsub('/(Geni|Tau|Unc)+J$', '', cc_pheno$strain)),
                      strain = unique(cc_pheno$strain),
                      date   = sapply(split(cc_pheno$date, cc_pheno$strain), unique),
                      aggregate(cc_pheno[,-(1:3)], list(cc_pheno$strain), mean, na.rm = TRUE))
cc_pheno = cc_pheno[,c('mouse', 'strain', 'date', 'Weight', 'NumbClutches', 'DAY1Numb1cells',     
                       'DAY1NumbDead', 'DAY1NumbFrag', 'TotalOocytes', 'DAY2Numb2cells',
                       'Fertilized', 'DAY2NumbDead', 'DAY2NumbFrag')]

# Founder phenotypes.
founder_pheno = founder_pheno[,c('Female number', 'JR', 'FreshIVF_IVFDISHES::IVF Date', 
                                 'Weight', 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 
                                 'DAY1NumbFrag', 'Total oocytes ovulated', 'DAY2Numb2cells',
                                 'Fertilizaed (%)', 'DAY2NumbDead', 'DAY2NumbFrag')]
colnames(founder_pheno) = c('mouse', 'jr', 'date', 
                                 'Weight', 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 
                                 'DAY1NumbFrag', 'TotalOocytes', 'DAY2Numb2cells',
                                 'Fertilized', 'DAY2NumbDead', 'DAY2NumbFrag')
founder_table = founder_table[,1:2]
colnames(founder_table) = c('strain', 'jr')
founder_table$strain = sub(' [0-9]+_[0-9]$', '', founder_table$strain)
founder_table$jr     = sub('-[0-9]$', '', founder_table$jr)
founder_table        = unique(founder_table)
founder_pheno = merge(founder_pheno, founder_table, by = 'jr', all.x = TRUE, sort = FALSE)
founder_pheno = founder_pheno[,c('mouse', 'strain', 'date', 
                                 'Weight', 'NumbClutches', 'DAY1Numb1cells', 'DAY1NumbDead', 
                                 'DAY1NumbFrag', 'TotalOocytes', 'DAY2Numb2cells',
                                 'Fertilized', 'DAY2NumbDead', 'DAY2NumbFrag')]
founder_pheno$Weight = as.numeric(sub('g$', '', founder_pheno$Weight))
founder_pheno = subset(founder_pheno, NumbClutches > 0)
founder_pheno$date = as.character(founder_pheno$date)
# Condense each strain to strain means.
founder_pheno = data.frame(mouse  = unique(founder_pheno$strain),
                      strain = unique(founder_pheno$strain),
                      date   = sapply(split(founder_pheno$date, founder_pheno$strain), unique),
                      aggregate(founder_pheno[,-(1:3)], list(founder_pheno$strain), mean, na.rm = TRUE))
founder_pheno = founder_pheno[,c('mouse', 'strain', 'date', 'Weight', 'NumbClutches', 
                                 'DAY1Numb1cells', 'DAY1NumbDead', 'DAY1NumbFrag', 'TotalOocytes',
                                 'DAY2Numb2cells', 'Fertilized', 'DAY2NumbDead', 'DAY2NumbFrag')]


stopifnot(colnames(pheno) == colnames(cc_pheno))
stopifnot(colnames(pheno) == colnames(founder_pheno))

# Combine phenotypes.
pheno = rbind(pheno, cc_pheno)
pheno = rbind(pheno, founder_pheno)
rownames(pheno) = pheno$mouse

# Select covariates that we need for mapping phenotypes.
non_do = which(pheno$strain != 'DO')
covar = rbind(covar, data.frame(mouse = pheno$mouse[non_do], 
                                sex   = 'F',
                                gen   = ifelse(grepl('^CC', pheno$strain[non_do]), '4', '1'),
                                row.names = pheno$mouse[non_do]))

covar = merge(covar, pheno[,c('mouse', 'strain', 'date', 'TotalOocytes', 'NumbClutches', 
                              'DAY1Numb1cells', 'DAY2Numb2cells')], by = 'mouse', all = TRUE,
                              sort = FALSE)
rownames(covar) = covar$mouse

stopifnot(pheno$mouse %in% covar$mouse)

covar      = covar[rownames(pheno),]
covar$date = factor(covar$date)
covar$gen  = factor(covar$gen)
pheno      = as.matrix(pheno[,-(1:3)])

stopifnot(rownames(pheno) == rownames(covar))


# Combine the DO and CC genoprobs.
probs = do_probs

# Subset probs to contain the correct markers and samples.
# Also remove markers with allele frequencies less than 1 for any one founder.
for(chr in seq_along(probs)) {

  print(chr)

  # Keep unique DO samples.
  probs[[chr]] = probs[[chr]][!duplicated(rownames(probs[[chr]])),,]
  
  samples = c(rownames(probs[[chr]]), rownames(cc_probs[[chr]]))
  
  common_markers = intersect(dimnames(probs[[chr]])[[3]], dimnames(cc_probs[[chr]])[[3]])
  probs[[chr]]    = probs[[chr]][,,common_markers]
  cc_probs[[chr]] = cc_probs[[chr]][,,common_markers]
  stopifnot(dimnames(probs[[chr]])[[3]] == dimnames(cc_probs[[chr]])[[3]])
  
  curr_pr = array(0, dim = c(length(samples), ncol(probs[[chr]]), dim(probs[[chr]])[3]), 
                  dimnames = list(samples, colnames(probs[[chr]]), dimnames(probs[[chr]])[[3]]))
  
  curr_pr[1:nrow(probs[[chr]]),,] = probs[[chr]]
  curr_pr[(nrow(probs[[chr]]) + 1):nrow(curr_pr),,] = cc_probs[[chr]]
  
  probs[[chr]] = curr_pr
  
  # Set low allele frequencies equal to zero.
  probs[[chr]][probs[[chr]] < probs_threshold] = 0
  
  # Rescale the probs to sum to 1 on each row. 
  for(i in 1:dim(probs[[chr]])[[3]]) {
  
    probs[[chr]][,,i] = probs[[chr]][,,i] / rowSums(probs[[chr]][,,i]) 
  
  } # for(i)
  
  # Synch up markers between map & probs.
  map[[chr]] = map[[chr]][dimnames(probs[[chr]])[[3]]]
  
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
norm_expr = readr::read_csv(expr_file) %>%
              as.data.frame()
annot     = readr::read_csv(annot_file) %>%
              as.data.frame()

rownames(norm_expr) = norm_expr$ensembl
norm_expr           = as.matrix(norm_expr[,-1]) 
rownames(annot)     = annot$gene_id 

stopifnot(rownames(norm_expr) == rownames(annot))

# Filter expression to retain expressed genes.
rm = rowMeans(norm_expr)
norm_expr = norm_expr[rm > 7,]

# Create RankZ transformed expression.
rz_expr = apply(norm_expr, 2, rankZ)

# Save Rdata file with all QTL mapping data.
save(pheno, pheno_rz, covar, probs, K, map, norm_expr, rz_expr, annot,
     file = file.path(data_dir, paste0('superovulation_qtl2_', genome, '_120K.Rdata')))



