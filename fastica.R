library(fastICA)
setwd('/projects/bolcun-filas-lab/DO_Superovulation/')
load('data/superovulation_qtl2_grcm38_137K.Rdata')
rm(probs, K)

norm_expr = t(norm_expr)

# These samples wree outliers when I ran fastICA the first time. 
# I can't see any mean or variance reason why.
norm_expr = norm_expr[!rownames(norm_expr) %in% c('SODO293', 'SODO3', ' SODO53'),]

res = fastICA(X       = norm_expr,
              n.comp  = 20,
              alg.typ = 'parallel', 
              method  = 'R', 
              verbose = TRUE)
saveRDS(res, file = 'results/fastica.rds')

