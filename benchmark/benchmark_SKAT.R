#!/usr/bin/env Rscript
require(SKAT)
require(data.table)

eprintf <- function(...) {
  cat(paste(sprintf(...), "\n"), sep='', file=stderr())
  flush.console()
}

args = commandArgs(trailingOnly = T)
# setwd("~/CLionProjects/Permute_Associate/")
# args = c("test.sim/test_data.txt", "test.sim/test_data.cov", "test.sim/test_data.ped")

run_analysis <- function(args) {

  matrix_file = args[1]
  covar_file = args[2]
  ped_file = args[3]

  eprintf("Loading matrix.")
  df <- fread(matrix_file, header=T)
  eprintf("Loaded matrix.")
  cov <- fread(covar_file, header=F)
  eprintf("Loaded covariates.")
  ped <- fread(ped_file, header=T)
  eprintf("Loaded ped.")

  df <- as.data.frame(df)
  colnames(df)[1] = "Gene"
  colnames(df)[2] = "Transcript"
  colnames(df)[3] = "Location"
  df$Gene = as.factor(df$Gene)
  cov <- as.data.frame(cov)
  ped <- as.data.frame(ped)

  case_control = as.vector(ped[,6])
  for(i in 1:length(case_control)) {
    if(case_control[i] == 2) {
      case_control[i] = 1
    } else {
      case_control[i] = 0
    }
  }
  
  nperm <- 10000

  cov$V1 = case_control
  null_formula <- paste(c("V1 ~ 1", colnames(cov)[3:ncol(cov) - 1]), sep = " + ", collapse = " + ")
  eprintf(paste("formula: ", null_formula))
  obj <- SKAT_Null_Model(formula(null_formula), data=cov, out_type = "D", n.Resampling = nperm, type.Resampling = "permutation", Adjustment = F)
  
  for(gene in levels(df$Gene)) {
    s = df[df$Gene == gene,]
    s$Transcript = factor(s$Transcript)
    
    for(ts in levels(s$Transcript)) {
      Z = t(as.matrix(s[s$Transcript == ts,4:ncol(s)]))
      Z[Z > 2] = 0
      res <- SKAT(Z, obj, kernel="linear", method="davies")
      eprintf("%s %s %6.5f", gene, ts, res$p.value)
    }
  }
}

system.time(run_analysis(args))
