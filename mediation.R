################################################################################
# Mediation analysis function.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2022-10-19
################################################################################

library(qtl2)

# Mediation analysis using RNAseq data from the same samples.
# Arguments:
# pheno: Numeric matrix containing phenotype. All dimensions are named.
#        Samples in rows, phenotypes in columns. Should contain one phenotype.
# covar: Numeric matrix containing covariates. All dimensions are named.
#        Samples in rows, covariates in columns. Should be created with
#        model.matrix() with the intercept removed and be ready for qtl2
# probs: qtl2-style list of 3D genoprobs arrays. All dimensions are named.
#        In each slice, samples in rows, founders in columns, markers in slices.
# K:     qtl2-style list of LOCO kinship matrices. All dimensions are named.
#        In each matrix, samples in rows, samples in columns.
# map:   qtl2-style list of named numeric vectors containing marker positions.
#        Marker names in names.
# expr:  Numeric matrix containing normalized expression data. Genes in 
#        rownames, samples in colnames.
# annot: Data.frame containing gene annotation. Genes in rows.
# chr:   Character. Chromosome where mediation should be done.
mediation = function(pheno, covar, probs, K, map, expr, annot, chr) {
  
  tmp   = synch_samples(pheno, covar, probs, K, expr)
  pheno = tmp$pheno
  covar = tmp$covar
  probs = tmp$probs
  K     = tmp$K
  expr  = tmp$expr

  tmp   = synch_markers(probs, map)
  probs = tmp$probs
  map   = tmp$map
  
  tmp   = synch_genes(expr, annot)
  expr  = tmp$expr
  annot = tmp$annot
  rm(tmp)

  # Map the base phenotype first.
  qtl = scan1(genoprobs = probs, pheno = pheno, kinship = K, 
              addcovar = covar)
  
  # Get the peak on the requested chromosome.
  peaks = find_peaks(qtl, map)
  peaks = peaks[peaks$chr == chr,]
  
  # Fit the model just at the marker with the highest LOD.
  max_pos  = max(qtl, map)
  pr       = pull_genoprobpos(probs, map, max_pos$chr, max_pos$pos)
  base_mod = fit1(genoprobs = pr, pheno = ph, kinship = K[[chr]], addcovar = addcovar)
  
  # Get the genes on the requested chromosome.
  chr_annot = annot[annot$chr == chr,]
  chr_expr  = expr[rownames(chr_annot),]

  # Mediate the phenotype with each gene.
  retval = data.frame(gene_id  = rownames(chr_annot),
                      symbol   = chr_annot$symbol,
                      chr      = chr_annot$chr,
                      pos      = chr_annot$start * 1e-6,
                      base_lod = base_mod$lod,
                      med_lod  = 0)
  for(i in 1:nrow(chr_expr)) {
    
    curr_covar = cbind(covar, chr_expr[i,])
    med_mod    = fit1(genoprobs = pr, pheno = ph, kinship = K[[chr]],
                      addcovar = curr_covar)
    retval$med_lod[i] = med_mod$lod
    
  } # for(i)
  
  retval$lod_drop = retval$med_lod - retval$base_lod

  return(retval)

} # mediation()


# Synch samples between objects.
synch_samples = function(pheno, covar, probs, K, expr) {
  
  common_samples = intersect(rownames(pheno), rownames(probs[[1]]))
  common_samples = intersect(common_samples,  colnames(expr))
  
  pheno = pheno[common_samples,]
  covar = covar[common_samples,]
  probs = probs[common_samples,]
  K     = lapply(K, function(z) { z[common_samples, common_samples] })
  expr  = expr[,common_samples]
  
  return(list(pheno = pheno, covar = covar, probs = probs, K = K, expr = expr))

} # synch_samples()


synch_markers = function(probs, map) {
  
  stopifnot(length(probs) == length(map))
  stopifnot(names(probs)  == names(map))
  
  for(i in seq_along(probs)) {

    common_markers = intersect(dimnames(probs[[i]])[[3]], names(map[[i]]))
    probs[[i]] = probs[[i]][,,common_markers]
    map[[i]]   = map[[i]][common_markers]

  } # for(i)
  
  return(list(probs = probs, map = map))
  
} # synch_markers()


synch_genes = function(expr, annot) {
  
  common_genes = intersect(rownames(expr), rownames(annot))
  expr         = expr[common_genes,]
  annot        = annot[common_genes,]
  
  return(list(expr = expr, annot = annot))
  
} # synch_genes()

