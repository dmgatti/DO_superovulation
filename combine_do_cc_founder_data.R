################################################################################
# Combine the DO, CC, and founder phenotype and genoprobs data.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-06-09
################################################################################
options(stringsAsFactors = FALSE)

##### LIBRARIES #####

library(qtl2convert)
library(qtl2)
library(readxl)
library(DESeq2)
library(tidyverse)
library(lubridate)


##### FUNCTIONS #####

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

##### VARIABLES #####

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results')

# DO genoprob file.
do_probs_file = file.path(data_dir, 'sodo_alleleprobs_grcm39_137K.rds') 

# CC genoprobs file.
cc_probs_file = '/projects/compsci/vmp/USERS/dgatti/data/cc/genoprobs/cc_consensus_genoprobs_137K.rds'

# DO phenotype file.
do_pheno_file = file.path(data_dir, 'JDO9376 all 350.xlsx')

# CC phenotype file.
cc_pheno_file = file.path(data_dir, 'CC lines export.xlsx')

# CC metadata file.
cc_meta_file  = '/projects/compsci/vmp/USERS/dgatti/data/cc/cc_strain_jr_table.csv'

# Founder phenotype file.
fs_pheno_file = file.path(data_dir,   'Founder strains export.xlsx')

# Gene annotation file.
annot_file   = file.path(results_dir, 'gene_annotation_ens105.csv')

# GBRS DO expression file.
do_expr_file = file.path(results_dir, 'sodo_gene_counts_grcm39_corrected.csv')

# GBRS CC expression file.
cc_expr_file = file.path(results_dir, 'socc_gene_counts_grcm39.csv')

# CC expression metadata file.
cc_meta_file = file.path(results_dir, 'socc_metadata.csv')

# Founder expression file.
fs_expr_file = file.path(results_dir, 'sofs_gene_counts.csv')

# Founder expression metadata file.
fs_meta_file = file.path(results_dir, 'sofs_metadata.csv')

##### MAIN #####

##############################
# Read in the phenotype files.

# DO phenotypes.
# NOTE: Removing samples with no data and setting Day 2 phenotypes to NA where 
# Dish Notes indicate trouble.  
do_pheno = read_xlsx(do_pheno_file, sheet = 'Sheet1') %>%
             mutate(`Dish Notes` = if_else(is.na(`Dish Notes`), '', `Dish Notes`),
                    `IASE injection time` = str_replace(`IASE injection time`, '^IASe():? ', ''),
                    `IASE injection time` = parse_time(`IASE injection time`, format = '%H:%M')) %>%
             filter(!str_detect(`Dish Notes`, 'female injured during injections and was euthanized')) %>%
             filter(!is.na(NumbClutches)) %>%
             mutate(DAY2Numb2cells = if_else(`Dish Notes` == ' cumulus cells still attached at changeover, only partial sperm drop entered IVF drop', NA, DAY2Numb2cells),
                    DAY2NumbDead  = if_else(`Dish Notes` == ' cumulus cells still attached at changeover, only partial sperm drop entered IVF drop', NA, DAY2NumbDead),
                    DAY2NumbFrag  = if_else(`Dish Notes` == ' cumulus cells still attached at changeover, only partial sperm drop entered IVF drop', NA, DAY2NumbFrag),
                    `Fertilized (%)` = if_else(`Dish Notes` == ' cumulus cells still attached at changeover, only partial sperm drop entered IVF drop', NA, `Fertilized (%)`),
                    DAY2Numb2cells = if_else(`Dish Notes` == ' female injured during injections and was euthanized', NA, DAY2Numb2cells),
                    DAY2NumbDead  = if_else(`Dish Notes` == ' female injured during injections and was euthanized', NA, DAY2NumbDead),
                    DAY2NumbFrag  = if_else(`Dish Notes` == ' female injured during injections and was euthanized', NA, DAY2NumbFrag),
                    `Fertilized (%)` = if_else(`Dish Notes` == ' female injured during injections and was euthanized', NA, `Fertilized (%)`),
                    DAY2Numb2cells   = if_else(`Total oocytes ovulated` == 0, NA, DAY2Numb2cells),
                    DAY2NumbDead     = if_else(`Total oocytes ovulated` == 0, NA, DAY2NumbDead),
                    DAY2NumbFrag     = if_else(`Total oocytes ovulated` == 0, NA, DAY2NumbFrag),
                    `Fertilized (%)` = if_else(`Total oocytes ovulated` == 0, NA, `Fertilized (%)`),
                    strain = 'DO',
                    `Female number` = str_c('SODO', `Female number`)) %>%
             select(-`IASE injection time`) %>%
             rename_with(str_to_lower)

print(str_c('There are ', nrow(do_pheno), ' DO mice.'))

## CC phenotypes.
cc_pheno = read_xlsx(cc_pheno_file, sheet = 'Sheet1') %>%
             filter(!is.na(DAY1Numb1cells)) %>%
             mutate(`Female number` = str_c('SOCC', as.character(`Female number`))) %>%
             rename(id = `Female number`, jr = `Line number`)
cc_meta  = readr::read_csv(cc_meta_file, show_col_types = FALSE) %>%
             mutate(id = str_replace(id, '_', ''))
cc_pheno = left_join(cc_pheno, cc_meta, by = 'id') %>%
             dplyr::rename(strain = cc_strain,
                           'FreshIVF_IVFDISHES::female born date' = `Female born date`) %>%
             select(-`IASE injection time`) %>%
             mutate('female number' = id,
                    strain = str_replace_all(strain, '/(Geni|Tau|Unc)+J$', '')) %>%
             rename_with(str_to_lower)

print(str_c('There are ', nrow(cc_pheno), ' CC mice.'))

# Founder phenotypes.
fs_meta  = read_xlsx(fs_pheno_file, sheet = 'Sample labels info') %>%
             rename('Female number' = `Female number RS`) %>%
             mutate(jr        = str_replace(`Female number`, '-[0-9]$', ''),
                    `GT Code` = str_replace(`GT Code`,       '_',       '')) %>%
             separate(`Sample identity`, into = c('strain', 'id'), sep = ' ') %>%
             rename_with(str_to_lower)

fs_pheno = read_xlsx(fs_pheno_file, sheet = 'Sheet1') %>%
             filter(!is.na(DAY1Numb1cells)) %>%
             rename('FreshIVF_IVFDISHES::female born date' = `female born date`) %>%
             select(-`IASE Injection time`) %>%
             rename_with(str_to_lower) %>%
             mutate(jr = as.character(jr)) %>%
             left_join(select(fs_meta, `female number`, `gt code`, strain) %>% distinct(), 
                       by = 'female number') %>%
             select(-`female number`) %>%
             rename(`female number` = `gt code`)

print(str_c('There are ', nrow(fs_pheno), ' Founder mice.'))

# Combine all of the phenotype data.
common_colnames = intersect(colnames(cc_pheno), colnames(do_pheno))
common_colnames = intersect(common_colnames,    colnames(fs_pheno))

do_pheno = do_pheno[,common_colnames]
cc_pheno = cc_pheno[,common_colnames]
fs_pheno = fs_pheno[,common_colnames]

pheno = bind_rows(do_pheno, cc_pheno, fs_pheno) %>%
          select(id       = `female number`,
                 coat     = `coat color`,
                 strain,
                 birth_date = `freshivf_ivfdishes::female born date`,
                 ivf_date   = `freshivf_ivfdishes::ivf date`,
                 ivf_time,
                 n_clutches    = numbclutches,
                 n_1cells      = day1numb1cells,
                 n_1cells_dead = day1numbdead,
                 n_1cells_frag = day1numbfrag,
                 n_oocytes     = `total oocytes ovulated`,
                 n_2cells      = day2numb2cells,
                 n_2cells_dead = day2numbdead,
                 n_2cells_frag = day2numbfrag) %>%
          mutate(ivf_date = parse_date_time2(str_c(year(ivf_date), '-', month(ivf_date), '-',
                            day(ivf_date), ' ' , hour(ivf_time), ':',  minute(ivf_time), ':',
                            second(ivf_time)), orders = 'Ymd HMS')) %>%
          select(-ivf_time)

print(str_c('There are ', nrow(pheno), ' total mice.'))

# Write out phenotypes.
readr::write_csv(pheno, file = file.path(data_dir, 'superovulation_all_pheno.csv'))


####################
# Read in genoprobs.
do_probs = readRDS(do_probs_file)
cc_probs = readRDS(cc_probs_file)

# Get rownames for new probs.
new_rownames = c(pheno$id[grep('^(CC|DO)', pheno$strain)])

# Get match between CC strains and CC samples ids.
cc_map = cc_pheno %>%
           select(id = `female number`, strain) %>%
           mutate(mapping = match(strain, rownames(cc_probs[[1]])))

n_do     = nrow(do_probs[[1]])
n_cc     = length(grep('SOCC', pheno$id))

new_probs = NULL
new_rn   = c(rownames(do_probs[[1]]), pheno$id[grep('SOCC', pheno$id)])

for(i in seq_along(do_probs)) {

  common_markers = intersect(dimnames(do_probs[[i]])[[3]], 
                             dimnames(cc_probs[[i]])[[3]])

  new_probs[[i]] = array(0, dim = c(length(new_rn), ncol(do_probs[[i]]),  length(common_markers)), 
                         dimnames = list(new_rn, colnames(do_probs[[i]]), common_markers))
  
  new_probs[[i]][rownames(do_probs[[i]]),,] = do_probs[[i]][,,common_markers]
  new_probs[[i]][cc_map$id,,] = cc_probs[[i]][cc_map$mapping,,common_markers]
  

} # for(i)

attributes(new_probs) = attributes(do_probs)

saveRDS(new_probs, file = file.path(data_dir, 'superovulation_all_probs_grcm39.rds'))

rm(do_probs, cc_probs, new_probs)
gc()

#########################
# Gather expression data.

# Read expression data.
do_expr = readr::read_csv(do_expr_file)
cc_expr = readr::read_csv(cc_expr_file)
fs_expr = readr::read_csv(fs_expr_file)

# Get common genes.
common_genes = intersect(do_expr$ensembl, cc_expr$ensembl)
common_genes = intersect(common_genes,    fs_expr$ensembl)

# Subset expression to contain the same genes.
do_expr = do_expr %>%
            filter(ensembl %in% common_genes)
cc_expr = cc_expr %>%
            filter(ensembl %in% common_genes)
fs_expr = fs_expr %>%
            filter(ensembl %in% common_genes)

all_expr = left_join(do_expr, cc_expr, by = 'ensembl') %>%
             left_join(fs_expr, by = 'ensembl')

# Read in metadata files.
cc_meta = read_csv(cc_meta_file, show_col_types = FALSE) %>%
            rename(strain = cc_strain) %>%
            mutate(id = str_replace(id, '_', ''))
fs_meta = read_csv(fs_meta_file, show_col_types = FALSE)

metadata = data.frame(id     = colnames(all_expr)[-1]) %>%
             left_join(bind_rows(cc_meta, fs_meta)) %>%
             mutate(strain = if_else(str_detect(id, 'SODO'), 'DO', strain)) %>%
             select(-jr, -num, -new_id)

# Normalize counts using DESeq's VST.
counts = all_expr %>%
             column_to_rownames(var = 'ensembl') %>%
             as.matrix()

# Subset counts to retain most highly expressed genes.
counts = counts[rowMeans(log1p(counts), na.rm = T) > 3,]

# There are NA values for the C57BL/6J samples (427 - 431) for 90 genes.
# I'm not sure why right now, but I need to revisit this in GRCm39.
# Remove these genes for now.
counts = counts[rowSums(is.na(counts)) == 0,]

print(str_c('There are ', nrow(counts), ' genes.'))

counts = round(counts)
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData   = metadata, 
                             design    = ~ strain)
dds = DESeq(dds)
vst = assay(vst(dds))
vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ensembl') %>%
  write_csv(file = file.path(results_dir, 'soall_gene_counts_norm_grcm39.csv'))



