#!/usr/bin/env Rscript
require(SKATR)
require(CompQuadForm)
require(data.table)

eprintf <- function(...) {
  cat(paste(sprintf(...), "\n"), sep='', file=stderr())
  flush.console()
}

# args = commandArgs(trailingOnly = T)
setwd("~/CLionProjects/Permute_Associate/")
# args = c("test.sim/test_data.txt", "test.sim/test_data.cov", "test.sim/test_data.ped")
args = c("~/Downloads/MDA_Ovarian/case_control.report.matrix.sort.txt", "~/Downloads/MDA_Ovarian/samples_pca.matrix.txt", "~/Downloads/MDA_Ovarian/xqc.new.ped")

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
  
  names <- colnames(df)[4:ncol(df)]
  name_order <- NULL
  for(n in names) {
    name_order <- c(name_order, which(ped[,2] == n))
  }
  
  case_control = ped[name_order, 6]
  for(i in 1:length(case_control)) {
    case_control[i] = ifelse(case_control[i] == 2, 1, 0)
  }
  
  name_order <- NULL
  for(n in names) {
    name_order <- c(name_order, which(cov$V1 == n))
  }
  
  cov$V1 = case_control
  
  for(i in 2:ncol(cov)) {
    cov[,i] = cov[name_order,i]
  }

  obj <- KAT.null(as.vector(cov$V1), as.matrix(cov[,2:ncol(cov)]))
  
  for(gene in levels(df$Gene)) {
    if(gene != "BRCA1") {
      next;
    }
    s = df[df$Gene == gene,]
    s$Transcript = factor(s$Transcript)
    
    for(ts in levels(s$Transcript)) {
      Z = t(as.matrix(s[s$Transcript == ts,4:ncol(s)]))
      Z[Z > 2] = 0
      res <- SKATh(obj, Z)
      eprintf("%s %s %6.5f", gene, ts, res)
    }
  }
}

run_analysis(args)

system.time(run_analysis(args))

