################################################################################
# Create the GBRS sample metadata file for the DO samples.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-07-11
################################################################################

base_dir = '/projects/bolcun-filas-lab/DO_Superovulation'
data_dir = file.path(base_dir, 'data')

fastq_dir = '/fastscratch/dgatti/fastq'

# Craete metadata data frame.
meta = data.frame(sampleID = 1:350,
                  sex      = 'F',
                  generation = rep(46:48, c(150, 150, 50)),
                  lane     = '',
                  fastq_1  = '',
                  fastq_2  = '')

for(i in 1:nrow(meta)) {

  files = dir(fastq_dir, pattern = paste0('SODO(_)?', i, '_'),
              full.names = TRUE)

  if(length(files) > 0) {
  
    meta$fastq_1[i] = files[1]
    meta$fastq_2[i] = files[2]
  
    lane         = strsplit(files[1], '_')[[1]]
    meta$lane[i] = lane[grep('^L', lane)]
    
  } # if(length(files) > 0)

} # for(i()

meta = meta[meta$lane != '',]

write.csv(meta, file = file.path(data_dir, 'sodo_gbrs_metadata.csv'),
          quote = FALSE, row.names = FALSE)
