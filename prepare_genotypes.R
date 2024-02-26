################################################################################
# Gather genotype data for haplotype construction using GRCm38.
# 
# Daniel Gatti
# dan.dgatti@jax.org
# 2022-12-05
################################################################################
options(stringsAsFactors = FALSE)

library(data.table)
library(tidyverse)
library(qtl2convert)
library(qtl2)

##### VARIABLES #####

genome = 'grcm39'

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
geno_dir    = file.path(data_dir, 'genotypes')
neogen_dirs = dir(geno_dir, pattern = 'Jackson')

meta_file = file.path(geno_dir, 'do_superovulation_animal_info.csv')

scratch_dir = '/fastscratch/dgatti'
qtl2_dir    = file.path(scratch_dir, 'sodo_qtl2')

marker_dir  = '/projects/omics_share/mouse/GRCm38/supporting_files/muga_annotation/gigamuga'
marker_file = file.path(marker_dir, 'gm_uwisc_v1.csv')

if(genome == 'grcm39') {

  marker_dir  = '/projects/omics_share/mouse/GRCm39/supporting_files/muga_annotation/gigamuga'
  marker_file = file.path(marker_dir, 'gm_uwisc_v4.csv')

} # if(genome == 'grcm39')

# OLd file with all markers. 
#founder_geno_file = '/projects/compsci/vmp/USERS/dgatti/data/muga/gigamuga_consensus_geno.rds'

# New file with selected markers.
founder_geno_file = '/projects/compsci/vmp/USERS/dgatti/data/muga/GigaMUGA_founder_consensus_genotypes.csv'

acgt = setNames(c('A', 'C', 'G', 'T'), c('T', 'G', 'C', 'A'))

##### MAIN #####

dir.create(qtl2_dir, showWarnings = FALSE)

# Read in markers.
markers = readr::read_csv(marker_file)
alleles = NULL

if(genome == 'grcm38') {

  markers = markers[,c('marker', 'chr', 'bp_mm10', 'cM_cox')]
  markers[,'pos'] = markers[,'bp_mm10'] * 1e-6

} else {

  markers = markers[,c('marker', 'chr', 'bp_grcm39', 'cM_cox', 'snp')]
  markers[,'pos'] = markers[,'bp_grcm39'] * 1e-6
  alleles = strsplit(markers$snp, '')
  alleles = matrix(unlist(alleles), nrow = length(alleles), ncol = 2,
                   byrow = TRUE, dimnames = list(markers$marker, c('A', 'B'))) 

} # else

# Read in the FinalReport files.
geno  = NULL
inten = NULL

for(i in 1:length(neogen_dirs)) {

  print(neogen_dirs[i])

  fr_file = dir(file.path(geno_dir, neogen_dirs[i]), pattern = 'FinalReport',
                full.names = TRUE)

  fr = fread(cmd = paste0('gunzip -cq ', fr_file), skip = 9,
             sep = '\t', header = TRUE, colClasses = rep(c('character', 'numeric'), c(8, 3))) %>%
         select(`SNP Name`:`Allele1 - Forward`, `SNP Name`:`Allele2 - Forward`, X:Y) %>%
         rename(marker = `SNP Name`, id = `Sample ID`) %>%
         unite('geno', `Allele1 - Forward`, `Allele2 - Forward`, sep = '')

  g = fr %>%
        select(marker, id, geno) %>%
        pivot_wider(names_from = id, values_from = geno)

  if(is.null(geno)) {
    geno = g
  } else {
    geno = full_join(geno, g, by = 'marker')
  } # else

  ii = fr %>%
         select(marker, id, X, Y) %>%
         left_join(markers, by = 'marker') %>%
         filter(chr %in% c('X', 'Y')) %>%
         pivot_longer(cols = X:Y, names_to = 'channel', values_to = 'intensity') %>%
         group_by(chr, id) %>%
         summarize(mean_x = mean(intensity, na.rm = TRUE),
                   .groups = 'drop') %>%
         pivot_wider(names_from = chr, values_from = mean_x) %>%
         mutate(id = as.character(id))

  rm(g, fr)
  gc()

  inten = bind_rows(inten, ii)

} # for(i)

# Change sample names to start with "SODO" and remove samples that are not part
# of this dataset.
geno = geno %>%
         select(-starts_with('DO_QTL_T')) %>%
         rename_with(.fn = str_replace, pattern = '_T$', replacement = '') %>%
         rename_with(.fn = str_replace, pattern = '^SO', replacement = 'SODO')

inten = inten %>%
          filter(!startsWith(id, 'DO_QTL_T')) %>%
          mutate(id = str_replace(id, pattern = '_T$', replacement = ''),
                 id = str_replace(id, pattern = '^SO', replacement = 'SODO'))

inten = inten[match(colnames(geno)[-1], inten$id),]
stopifnot(colnames(geno)[-1] == inten$id)

geno = as.data.frame(geno)

saveRDS(geno,  file = file.path(data_dir, 'genotypes.rds'))
saveRDS(inten, file = file.path(data_dir, 'xy_intensities.rds'))

######################
# Prepare qtl2 files.

# Founder genotypes.
fgeno = read.csv(founder_geno_file)
#fgeno = readRDS(founder_geno_file)
fgeno[is.na(fgeno)] = 'N'
#fgeno = data.frame(marker = rownames(fgeno), fgeno)

# Genetic and physical maps.

# Filter markers to retain the markers in founder geno.
markers = markers[markers$marker %in% fgeno$marker,]
fgeno   = fgeno[match(markers$marker, fgeno$marker),]
stopifnot(markers$marker == fgeno$marker)

# pmap
markers %>%
  select(marker, chr, pos) %>%
  write_csv(file = file.path(qtl2_dir, 'pmap.csv'))

# gmap
markers %>%
  select(marker, chr, cM_cox) %>%
  rename(pos = cM_cox) %>%
  write_csv(file = file.path(qtl2_dir, 'gmap.csv'))

# Sample genotypes.
geno = readRDS(file = file.path(data_dir, 'genotypes.rds'))
geno = geno %>%
          filter(marker %in% markers$marker)

tmp = full_join(fgeno, geno, by = 'marker') %>%
        as.data.frame()
rownames(tmp) = tmp$marker
tmp = as.matrix(tmp[,-1])

# Build 39 file contains the correct alleles.
if(genome == 'grcm38') {

  b39_markers = read.csv('/projects/omics_share/mouse/GRCm39/supporting_files/muga_annotation/gigamuga/gm_uwisc_v2.csv')
  
  b39_markers = b39_markers %>%
                  dplyr::select(marker, snp) %>%
                  right_join(dplyr::select(markers, marker)) %>%
                  separate(snp, into = c('junk', 'A', 'B'), sep = '') %>%
                  select(-junk)

  alleles = matrix(c(b39_markers$A, b39_markers$B), ncol = 2,
                   dimnames = list(b39_markers$marker, c('A', 'B')))
  alleles = alleles[markers$marker,]

  stopifnot(rownames(alleles) == markers$marker)
  
} # if(genome == 'grcm38') 

alleles = alleles[rownames(tmp),]
stopifnot(nrow(alleles) == nrow(tmp))

tmp = encode_geno(geno = tmp, allele_codes = alleles)

fgeno = tmp[,1:8]
fgeno = data.frame(marker = rownames(fgeno), fgeno)
geno  = tmp[,-(1:8)]
geno  = data.frame(marker = rownames(geno), geno)
colnames(geno) = sub('^X', '', colnames(geno))

readr::write_csv(fgeno, file = file.path(qtl2_dir, 'founder_geno.csv'))
readr::write_csv(geno,  file = file.path(qtl2_dir, 'sample_geno.csv'))

rm(tmp, alleles)

# Phenotypes and covariates.
meta = readr::read_csv(meta_file) %>%
         select(id = `Animal ID`, sex = Sex, gen = `DO generation`) %>%
         mutate(sex = 'female',
                id  = str_replace(id, pattern = '_T$', replacement = ''),
                id  = str_replace(id, pattern = '^SO', replacement = 'SODO'))

covar = meta
readr::write_csv(covar, file = file.path(qtl2_dir, 'covar.csv'))

pheno = data.frame(id = covar$id, pheno = rnorm(nrow(covar)))
readr::write_csv(pheno, file = file.path(qtl2_dir, 'pheno.csv'))

# JSON control file

json = '{
  "description": "DO Superovulation Study",
  "crosstype": "do",
  "sep": ",",
  "na.strings": ["-", "NA"],
  "comment.char": "#",
  "geno": "sample_geno.csv",
  "founder_geno": "founder_geno.csv",
  "gmap": "gmap.csv",
  "pmap": "pmap.csv",
  "pheno": "pheno.csv",
  "covar": "covar.csv",
  "alleles": ["A", "B", "C", "D", "E", "F", "G", "H"],
  "x_chr": "X",
  "genotypes": {
    "A": 1,
    "H": 2,
    "B": 3
  },
  "geno_transposed": true,
  "founder_geno_transposed": true,
  "sex": {
    "covar": "sex",
    "female": "female",
    "male": "male"
  },
  "cross_info": {
    "covar": "gen"
  }
}'

writeLines(text = json, con = file.path(qtl2_dir, paste0('do_superovulation_', genome, '.json')))


