################################################################################
# Map the DO superovulation phenotypes using miqtl.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-05-23
################################################################################
options(stringsAsFactors = FALSE)

library(readxl)
library(tidyverse)
library(qtl2)
library(miqtl)

##### VARIABLES #####

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results')

# Allele probs file from qtl2.
probs_file = file.path(data_dir, 'sodo_genoprobs_grcm39_123K.rds')

# Cross object from qtl2.
cross_file = file.path(data_dir, 'sodo_cross_grcm39.rds')

# Phenotype file.
pheno_file  = file.path(data_dir, 'JDO9376 all 350.xlsx')

# Path to HAPPY genome cache.
happy_tmp_cache = '/fastscratch/dgatti/happy_tmp_cache'

# Path to HAPPY collapsed cache.
happy_cache = '/fastscratch/dgatti/happy_cache'

##### MAIN #####


# Read in the genoprobs object.
probs = readRDS(probs_file)

# If the HAPPY cache has been removed, create it again.
if(!file.exists(happy_cache)) {

  # Create HAPPY cache directories.
  dir.create(happy_tmp_cache, showWarnings = FALSE)
  dir.create(happy_cache,     showWarnings = FALSE)

  # Read in the cross object.
  cross = readRDS(cross_file)

  # Convert the qtl2-style data to HAPPY format.
  convert.qtl2.to.HAPPY(qtl2.object       = probs, 
                        cross.object      = cross, 
                        HAPPY.output.path = happy_cache, 
                        diplotype.order   = 'qtl2')

  # Reduce the size of the genome cache (if possible).
  collapse.genomecache(original.cache = happy_cache, 
                       new.cache      = happy_collapsed, 
                       criterion      = 'l2.norm',
                       model          = 'additive')
} # if(!file.exists(happy_cache))
        
# Read in phenotypes.             
pheno = read_xlsx(pheno_file) %>%
          select(id = `Female number`, ivf_date = `FreshIVF_IVFDISHES::IVF Date`,
                 weight = Weight, NumbClutches:DAY2NumbFrag) %>%
          rename_with(str_to_lower) %>%
          rename(total_oocytes = `total oocytes ovulated`,
                 pct_fertilized = `fertilized (%)`) %>%
          filter(id != 'Female 50') %>%
          mutate(weight = as.numeric(str_replace(weight, 'g', '')),
                 ivf_date = factor(ivf_date))




miqtl.rop.scan = scan.h2lmm(genomecache = happy_cache, 
                             data       = pheno, 
                             pheno.id   = 'id', 
                             geno.id    = 'id', 
                             formula    = total_oocytes ~ ivf_date, 
                             use.multi.impute      = FALSE, 
                             return.allele.effects = TRUE)

genome.plotter.whole(scan.list = list(ROP = miqtl.rop.scan))





