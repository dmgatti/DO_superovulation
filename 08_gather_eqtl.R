################################################################################
# Gather eQTL results once they have been run in parallel.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-14
################################################################################
options(stringsAsFactors = FALSE)

library(AnnotationHub)

##### VARIABLES #####

genome = 'grcm39'

base_dir    = '/projects/bolcun-filas-lab/DO_Superovulation'
results_dir = file.path(base_dir, 'results', paste0('eqtl_', genome))
marker_dir  = '/projects/omics_share/mouse/GRCm38/supporting_files/muga_annotation/gigamuga'
marker_file = file.path(marker_dir, 'gm_uwisc_v1.csv')

if(genome == 'grcm39') {

  marker_dir  = '/projects/omics_share/mouse/GRCm39/supporting_files/muga_annotation/gigamuga'
  marker_file = file.path(marker_dir, 'gm_uwisc_v4.csv')

} # if(genome == 'grcm39')

# This assumes that the 07_local_eqtl.R script has been run and that
# the *rds output files exist.
eqtl_files = dir(path = results_dir, pattern = '[0-9]\\.rds$', full.names = TRUE)

hub = AnnotationHub()
ensdb = NULL
if(genome == 'grcm38') {

  hub = query(hub, c('mus musculus', 'ensdb', '98'))
  ensdb = hub[['AH75036']]

} else {

  hub = query(hub, c('mus musculus', 'ensdb', '105'))
  ensdb = hub[['AH98078']]

} # else

annot = genes(ensdb)
annot = as.data.frame(annot)

##### MAIN #####

# Read in markers.
markers = read.csv(marker_file)
markers = markers[,1:5]
rownames(markers) = markers$marker
markers$pos = markers$bp_grcm39 * 1e-6

# Read in the eqtl results files.
results = lapply(eqtl_files, readRDS)

# Bind the results into one object.
results = do.call(cbind, results)

# Write out eqtl results.
#saveRDS(results, file.path(results_dir, '..', paste0('superovulation_eqtl_', genome,'_lod.rds')))

results = readRDS(file.path(results_dir, '..', paste0('superovulation_eqtl_', genome, '_lod.rds')))

# Subset the markers to retain the ones used in mapping.
markers = markers[rownames(results),]

stopifnot(rownames(markers) == rownames(results))

# For each gene, get the maximum LOD and it's location.
#max_lod     = apply(results, 2, max)
#max_lod_pos = apply(results, 2, which.max)

# Changing to a slower way to catch genes with multiple QTL and 
# to screen for false spikes.
# Make a data.frame with marker and chr to use when splitting results.
mkr           = markers[,c('marker', 'chr')]
rownames(mkr) = mkr$marker
mkr$chr       = factor(mkr$chr, levels = c(1:19, 'X'))
mkr = mkr[rownames(results),]
stopifnot(rownames(mkr) == rownames(results))

### Using base R because tidyverse threw an error. May revisit this.

# Got through each gene and harvest eQTL over LOD = 6.
max_qtl = NULL

for(i in 1:ncol(results)) {

  if(i %% 100 == 0) print(i)

  # Join the markers & lOD.
  lod = cbind(mkr, results[,i])
  colnames(lod)[3] = 'lod'
  
  # Split lod by chormosome.
  lod_max = split(lod, lod$chr)
  
  # Get the maximum LOD on each chromosome.
  lod_max = lapply(lod_max, function(z) { z[which.max(z$lod),,drop = FALSE] })
  lod_max = do.call(rbind, lod_max)
  
  # Retain peaks with LOD >= 6.
  lod_max = lod_max[lod_max$lod >= 6,]
  
  # Get the markers and LODs adjacent to the max peaks.
  num_loci = nrow(lod_max)
  
  if(num_loci > 0) {
  
    m = match(lod_max$marker, lod$marker)
  
    # Get the number of markers around each peak at which LOD >= 6.
    lod_max$peakwidth = NA
    for(j in seq_along(m)) {
  
      lod_max$peakwidth[j] = sum(lod$lod[max(1, (m[j]-20)):(m[j]+20)] >= 6)
  
    } # for(j)
  
     lod_max             = cbind(colnames(results)[i], lod_max)
     colnames(lod_max)[1] = 'gene_id'
 
     max_qtl = rbind(max_qtl, lod_max)
 
  } # if(num_loci > 0)
 
} # for(i)

# Get the median and the median absolute deviation of hte width of each QTL hotspot.
hotspots = sort(table(max_qtl$marker), decreasing = TRUE)
res = data.frame(marker = names(hotspots)[1:length(hotspots)], 
                 n      = as.vector(hotspots[1:length(hotspots)]),
                 lod_mean  = 0,
                 width_med = 0,
                 width_mad = 0)
for(i in seq_along(hotspots)) {

  x = subset(max_qtl, marker == names(hotspots)[i])
  res$lod_mean[i]  = mean(x$lod)
  res$width_med[i] = median(x$peakwidth)
  res$width_mad[i] = mad(x$peakwidth)

} # for(i)

# Retain hotspots with more than 5 genes.
res = subset(res, n > 5)

# Retain hotspots with a small width and zero mad width.
res = subset(res, width_med <= 5 & width_mad == 0)

# Get the genes in each of these hotspots. Remove the marker with
# the highest LOD & get the next highest.
for(i in 1:nrow(res)) {

  curr_mkr = res$marker[i]
  curr_chr = markers[curr_mkr,'chr'] 
  
  # Get this hotspot's genes.
  g = subset(max_qtl, marker == curr_mkr)

  # Get the LOD profiles for this hotspot's genes.
  lod = results[,g$gene_id]
  
  # Subset LODs to retain the ones on the chromosome where the current 
  # marker lies.
  chr_mkr = subset(markers, chr == curr_chr)
  lod     = lod[chr_mkr$marker,]
  
  # Remove the hotspot marker from the LODs.
  lod = lod[rownames(lod) != curr_mkr,]
  
  # Obtain the maximum LOD for each gene on this chr.
  curr_max_qtl = data.frame(gene_id = colnames(lod),
                            marker  = rownames(lod)[apply(lod, 2, which.max)],
                            chr     = curr_chr,
                            lod     = apply(lod, 2, max),
                            peakwidth = 0) 

  # Subset to retain peaks with LOD > 6.
  curr_max_qtl  = subset(curr_max_qtl, lod >= 6)
  lod = lod[,curr_max_qtl$gene_id,drop = FALSE]
  stopifnot(colnames(lod) == curr_max_qtl$gene_id)

  # j indexes genes in curr_max_lod & lod.
  for(j in 1:nrow(curr_max_qtl)) {
  
    m = curr_max_qtl$marker[j]
    m = which(rownames(lod) == m)
    curr_max_qtl$peakwidth[j] = sum(lod[max(1, (m-20)):(m+20), j] >= 6)
  
  } # for(j)

  # Remove the current marker from the max qtl results.
  max_qtl = subset(max_qtl, marker != curr_mkr)
  
  # Add the new results into max_qtl.
  max_qtl = rbind(max_qtl, curr_max_qtl)

} # for(i)

# Join the lod to the annotation.
colnames(max_qtl)[colnames(max_qtl) == 'chr'] = 'qtl_chr'
max_qtl = merge(max_qtl, 
                annot[,c('gene_id', 'seqnames', 'start', 'end', 'strand', 'symbol')], 
                by = 'gene_id', all.x = TRUE, sort = FALSE)
max_qtl = merge(max_qtl,
                markers[,c('marker', 'chr', 'pos')],
                by = 'marker', all.x = TRUE, sort = FALSE)

max_qtl = data.frame(ensembl    = max_qtl$gene_id, 
                     symbol     = max_qtl$symbol,  
                     gene_chr   = max_qtl$seqnames, 
                     gene_start = max_qtl$start, 
                     gene_end   = max_qtl$end, 
                     strand     = max_qtl$strand,
                     marker     = max_qtl$marker, 
                     qtl_chr    = max_qtl$qtl_chr,
                     qtl_pos    = max_qtl$pos,
                     lod        = max_qtl$lod)
                     
max_qtl$gene_start = max_qtl$gene_start * 1e-6
max_qtl$gene_end   = max_qtl$gene_end   * 1e-6 

write.csv(max_qtl, file = file.path(results_dir, '..', paste0('sodo_max_eqtl_summary_', genome, '.csv')),
          row.names = FALSE, quote = FALSE)

 