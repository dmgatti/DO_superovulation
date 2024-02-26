################################################################################
# Gather the data required to create a qtlviewer and format it as 
# described at: https://qtlviewer.jax.org/docs/RDataObject.pdf.
# DO superovulation phenotypes and ovary RNAseq.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-11-08
################################################################################
options(stringsAsFactors = FALSE)

library(qtl2api)
library(qtl2)
library(tidyverse)

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

###### VARIABLES #####

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
qtl2_file   = file.path(data_dir, 'superovulation_qtl2_grcm39.Rdata')
marker_dir  = '/projects/omics_share/mouse/GRCm39/supporting_files/muga_annotation/gigamuga'
marker_file = file.path(marker_dir, 'gm_uwisc_v4.csv')

# GBRS files say that they use Ensembl 105.
# https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/config/gbrs.config 
ensembl_release = 105

markers = readr::read_csv(marker_file) %>%
            select(marker:bp_mm10) %>%
            rename(marker_id = marker,
                   pos       = bp_mm10) %>%
            mutate(pos = pos * 1e-6)

# Fertility phenotypes
dataset.oocytes = NULL

# Ovary RNAseq dataset.
dataset.ovary_rnaseq = NULL


##### MAIN #####

# Read in the data from qtl mapping.
load(qtl2_file)

common_markers = unlist(sapply(map, names))
markers = markers[markers$marker_id %in% common_markers,]

genoprobs = probs
rm(probs)

# Oocyte phenotypes

pheno = data.frame(sample_id = rownames(pheno), pheno)

annot_phenotype = tibble(data_name  = colnames(pheno),
                         short_name = c('Mouse ID', 'Body Weight', 'Num Clutches', 
                                        'Healthy Oocytes: Day 1', 'Dead Oocytes: Day 1',
                                        'Fragmented Oocytes: Day 1',
                                        'Total Oocytes', 'Fertilized Oocytes: Day 2',
                                        'Percent Oocytes Fertilized', 'Dead Oocytes: Day 2', 
                                        'Fragmented Oocytes: Day 2'),
                         category   = rep('Fertility', ncol(pheno)),
                         description = c('Mouse ID', 'Body weight at euthanasia', 
                                        'Number of oocyte clutches',
                                        'Number of healthy oocytes on day 1',
                                        'Number of dead oocytes on day 1',
                                        'Number of fragmented oocytes on day 1',
                                        'Total oocytes (sum of all oocytes on day 1)',
                                        'Number of fertilized oocytes on day 2',
                                        'Percent of oocytes on day 1 which were fertilized on day 2',
                                        'Number of dead oocytes on day 2',
                                        'Number of fragmented oocytes on day 2'),
                         is_id      = c(TRUE,  rep(FALSE, ncol(pheno) - 1)),
                         is_pheno   = c(FALSE, rep(TRUE,  ncol(pheno) - 1)),
                         is_numeric = c(FALSE, rep(TRUE,  ncol(pheno) - 1)),
                         omit       = c(TRUE,  rep(FALSE, ncol(pheno) - 1)),
                         is_date    = rep(FALSE, ncol(pheno)),
                         is_factor  = rep(FALSE, ncol(pheno)),
                         factor_levels = NA,
                         use_covar  = c(NA, 'date', 'date', 'date:TotalOocytes', 
                                        'date:TotalOocytes:NumbClutches', 'date:DAY1Numb1cells',
                                        'date:DAY1Numb1cells', 'date:NumbClutches',
                                        'date', 'date:DAY2Numb2cells', 'date:DAY2Numb2cells'))

annot_samples = as_tibble(covar) %>%
                  rename(sample_id = mouse)

covar_info = tibble(sample_column = 'date',
                    display_name  = 'Date',
                    interactive   = FALSE,
                    primary       = TRUE,
                    lod_peaks     = NA)

data = list(raw = as.matrix(pheno[,-1]),
            log = as.matrix(log1p(pheno[,-1])),
            rz  = pheno_rz)
rownames(data$raw) = pheno$sample_id
rownames(data$log) = pheno$sample_id

lod_peaks = read_csv(file.path(base_dir, 'results', 'qtl2', 'sodo_qtl2_peaks.csv')) %>%
              select(lodcolumn:lod) %>%
              rename(data_name = lodcolumn) %>%
              mutate(marker_id = find_marker(map = map, chr = chr, pos = pos),
                     additive = TRUE) %>%
              select(-chr, -pos)


dataset.oocytes = list(annot_phenotype = annot_phenotype,
                       annot_samples   = annot_samples,
                       covar_info      = covar_info,
                       data            = data,
                       datatype        = 'phenotype',
                       display_name    = 'Oocyte Phenotypes',
                       lod_peaks       = lod_peaks)

# Ovary RNAseq

annot_mrna = tibble(gene_id = annot$gene_id,
                    symbol  = annot$gene_name,
                    chr     = annot$chr,
                    start   = annot$start * 1e-6,
                    end     = annot$end   * 1e-6)

annot_samples = annot_samples[annot_samples$sample_id %in% colnames(norm_expr),]

data = list(norm = t(norm_expr),
            rz   = apply(t(norm_expr), 2, rankZ))

lod_peaks = read_csv(file.path(base_dir, 'results', 'eqtl', 'sodo_max_eqtl_summary.csv')) %>%
              select(ensembl, qtl_chr:lod) %>%
              rename(gene_id   = ensembl) %>%
              mutate(marker_id = find_marker(map = map, chr = qtl_chr, pos = qtl_pos),
                     additive  = rep(TRUE, nrow(.))) %>%
              select(-qtl_chr, -qtl_pos)

dataset.ovary_rnaseq = list(annot_mrna    = annot_mrna,
                            annot_samples = annot_samples,
                            covar_info    = covar_info,
                            data          = data,
                            datatype      = 'mrna',
                            display_name  = 'Ovary RNAseq',
                            lod_peaks     = lod_peaks)

qtl2api::validate_environment()
                            
save(ensembl_version, genoprobs, K, map, markers, 
     file = file.path(base_dir, 'core.ovary_qtlviewer.Rdata'))

saveRDS(dataset.oocytes,      file = file.path(base_dir, 'dataset.phenotype.DO_oocytes.rds'))
saveRDS(dataset.ovary_rnaseq, file = file.path(base_dir, 'dataset.mrna.DO_oocytes.rds'))



