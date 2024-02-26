################################################################################
# Estimate the heritability of a full set of samples.
# Take subsets of samples and re-estimate heritability, recording the number
# of samples and the generation.
# Using 350 female DO mice form DO superovulation study.
#
# Danitl Gatti
# dan.gatti@jax.org
# 2023-03-11
################################################################################
options(stringsAsFActors = FALSE)

library(readxl)
library(qtl2)

##### VARIABLES #####

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir    = file.path(base_dir, 'data')
results_dir = file.path(base_dir, 'results', 'herit_sim')

pheno_file  = file.path(data_dir, 'JDO9376 all 350.xlsx')
probs_file = file.path(data_dir, 'sodo_alleleprobs_grcm38_137K.rds')

##### MAIN #####

# Read in sample metadata.
pheno = read_xlsx(pheno_file, sheet = 'Sheet1')
pheno = pheno[,c('Female number', 'FreshIVF_IVFDISHES::female born date', 'Weight',
               'NumbClutches',  'DAY1Numb1cells',  'DAY1NumbDead', 'DAY1NumbFrag', 
               'Total oocytes ovulated', 'DAY2Numb2cells', 'Fertilized (%)', 
               'DAY2NumbDead', 'DAY2NumbFrag')]
colnames(pheno) = c('id', 'birth_date', 'weight', 'n_clutches', 'healthy_1cells',
                    'dead_1cells', 'frag_1cells', 'oocytes', 'health_2cells',
                    'pct_fertilized', 'dead_2cells', 'frag_2cells')
pheno$id = paste0('SODO', pheno$id)
pheno = as.data.frame(pheno)
rownames(pheno) = pheno$id
pheno$birth_date = factor(pheno$birth_date)
pheno$weight = as.numeric(sub('g', '', pheno$weight)) 
pheno$gen = '46'
pheno$gen[pheno$birth_date == '6/21/2022'] = '47'
pheno$gen[pheno$birth_date %in%  c('7/12/2022', '7/21/2022')] = '48'

# Keep only 2 clutch samples.
pheno = subset(pheno, n_clutches == 2)

# Read in probs.
probs = readRDS(probs_file)

# Create kinship matrix.
K = calc_kinship(probs, type = 'overall')
rm(probs)

# Synch up samples.
common_samples = intersect(rownames(K), pheno$id)
pheno = pheno[common_samples,]
K     = K[common_samples, common_samples]

# Set up covariates.
addcovar = model.matrix(~gen, data = pheno)[,-1]

# Make pheno a matrix.
pheno = as.matrix(pheno[,c('weight', 'healthy_1cells', 'dead_1cells', 
                  'frag_1cells', 'oocytes', 'health_2cells',
                  'pct_fertilized', 'dead_2cells', 'frag_2cells')])

# Get overall kinship estimate.
h2_overall = est_herit(pheno = pheno, kinship = K, addcovar = addcovar)


# Set up simulations by sample size.
sample_sizes = c(50, 100, 150, 200, 250, 300)
# Number of trials per sample_size.
n_sim = 100

sim_record = data.frame(sample_size = nrow(pheno), 
                        sim         = 1,
                        samples     = rownames(pheno))
herit      = data.frame(sample_size = nrow(pheno), 
                        sim         = 1,
                        t(h2_overall))

for(s in sample_sizes) {

  print(paste('Sample size:', s))
  
  for(i in 1:n_sim) {
  
    print(paste('   Sim:', i))

    # Select samples.
    samples = sample(rownames(pheno), size = s)
    
    sim_record = rbind(sim_record,
                       data.frame(sample_size = s, 
                                  sim         = i,
                                  samples     = samples))
    
    sim_h2 = est_herit(pheno    = pheno[samples,],
                       kinship  = K[samples, samples],
                       addcovar = addcovar[samples,])

    herit = rbind(herit, data.frame(sample_size = s, 
                                    sim         = i,
                                    t(sim_h2)))
    
  
  } # for(i)
  
  saveRDS(sim_record, file = file.path(results_dir, 'herit_sim_samples.rds'))
  saveRDS(herit,      file = file.path(results_dir, 'herit_sim_results.rds'))

} # for(s)


# Plot results.
library(tidyverse)


png(file.path(results_dir, 'herit_by_sample_size.png'),
    width = 1000, height = 800, res = 128)
herit %>%
  pivot_longer(cols = weight:frag_2cells) %>%
  mutate(sample_size = factor(sample_size)) %>%
  ggplot(aes(sample_size, value)) +
    geom_boxplot() + 
    geom_point() +
    facet_wrap(~name) +
    labs(title = 'Heritability Estimates by Sample Size',
         x     = 'Sample Size', y = 'Heritability')
dev.off()         






