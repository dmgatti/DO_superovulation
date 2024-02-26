################################################################################
# Gather the RSEM counts files into a single file.
# Daniel Gatti
# dan.gatti@jax.org
# 2021-10-29 
################################################################################

# Arguments:
# 1. rsem directory: full path to rsem output directory containing the gene counts.
# 2. output file:    full path of the output file. Should end in 'csv'
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
  stop('usage: gather_rsem_counts.R <rsem directory> <output_file>')
}

RSEM_DIR = args[1]
OUT_FILE = args[2]

options(stringsAsFactors = FALSE)

gene_files = dir(RSEM_DIR, pattern = '\\.genes\\.results', full.names = TRUE)
samples = gsub(paste0('^', RSEM_DIR, '/|\\.genes\\.results$'), '', gene_files)

if(length(gene_files) == 0) {
  stop(paste('No gene counts file found in', RSEM_DIR))
}

# Read in all of the gene count files.
data = lapply(gene_files, read.delim)

# Get the gene IDs from the first file. (all files should have the same genes in the same order)
genes = data[[1]]$gene_id

# Keep the 'expected_count' column in each file.
data = lapply(data, function(z) { z[,'expected_count',drop = FALSE] })

# Gather the counts into a matrix.
data = do.call(cbind, data)
dimnames(data) = list(genes, samples)
data = data.frame(gene_id = rownames(data), data)

# Write out the results.
write.csv(data, file = OUT_FILE, quote = FALSE, row.names = FALSE)

print(paste('Created file:', OUT_FILE))
