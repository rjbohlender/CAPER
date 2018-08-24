#!/usr/bin/env Rscript
require(SKAT)

args = commandArgs(trailingOnly = T)

matrix_file = args[1]
covar_file = args[2]

df <- read.table(matrix_file, header=T)
covars <- read.table(covar_file, header=F)

run_analysis <- function(df, covars, nperm) {
  obj <- SKAT_Null_Model(formula("V1 ~ 1 + V2 + V3"), data=covars, out_type = "D", n.Resampling = nperm, type.Resampling = "permutation", Adjustment = T)
  
  for(gene in levels(df$Gene)) {
    s = df[df$Gene == gene,]
    s$Transcript = factor(s$Transcript)
    
    for(ts in levels(s$Transcript)) {
      Z = t(as.matrix(s[s$Transcript == ts,4:ncol(s)]))
      Z[Z > 2] = 0
      res <- SKAT(Z, obj, kernel="linear", method="optimal.adj")
      msg <- sprintf("%s %s %6.5f %6.5f", gene, ts, res$p.value, sum(res$p.value.resampling <= res$p.value) / nperm)
      print(msg)
    }
  }
}

system.time(run_analysis(df, covars, 1))