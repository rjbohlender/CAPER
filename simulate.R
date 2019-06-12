# Quick and dirty genetic simulation

library(Matrix)

d <- "./test.sim/"

print("Writing matrix.")

sink(paste0(d, "test.matrix"))
total <- 100000

# Log odds 
lodds <- log(10)

ngenes <- 1000
nsnps <- rpois(ngenes, 20)
ncausal <- 10

freq <- list()
for(i in 1:ngenes) {
  freq[[i]] <- runif(nsnps[i], 0.001, 0.01)
}

data = Matrix(0L, nrow = total, ncol = ngenes, sparse = T)
location = 0
for(i in 1:ngenes) {
  gene = Matrix(rbinom(total * nsnps[i], 2, freq[[i]]), nrow = nsnps[i], ncol = total, sparse = T)
  for(j in 1:nsnps[i]) {
    location = location + 1
    cat(
      paste(paste0("gene_", i), 
            paste0("transcript_", i), 
            paste("chr1", location, location, "SNV", sep="-"),
            paste(gene[j,], sep = "", collapse = '\t'),
            sep = "\t"))
    cat("\n")
  }
  data[, i] = colSums(gene)
}

sink()

print("Finished writing matrix.")

# Find best val for noncausal genes
func <- function(x) {
  ret <- NULL
  for(i in x) {
    b <- c(rep(lodds, ncausal), rep(log(i), ngenes - ncausal))
    ret <- c(ret, mean(binomial()$linkinv(as.matrix(data %*% b))) - 0.5)
  }
  ret
}

odds_nc <- uniroot(func, c(0.001, 2))$root

betas <- c(rep(lodds, ncausal), rep(log(odds_nc), ngenes - ncausal))
case_control <- rbinom(total, 1, binomial()$linkinv(as.matrix(data %*% betas)))

case_count = 0;
control_count = 0;

labels = NULL
for(v in case_control) {
  if(v == 1) {
    case_count = case_count + 1
    labels <- c(labels, sprintf("case_%g", case_count))
  } else {
    control_count = control_count + 1
    labels <- c(labels, sprintf("control_%g", control_count))
  }
}
header <- paste("#Gene\tTranscript\tLocation", paste(labels, collapse = "\t"), sep="\t")
header <- paste0(header, "\n")

# Output the matrix header

sink(paste0(d, "header.txt"))
cat(header)
sink()

# Write a ped file

sink(paste0(d, "test_data.ped"))
cat("#fid\tiid\tpid\tmid\tsex\taff\n")
for(i in 1:length(labels)) {
  if(case_control[i] == 1) {
    cat(0, labels[i], 0, 0, 0, 2, sep = "\t")
  } else {
    cat(0, labels[i], 0, 0, 0, 1, sep = "\t")
  }
  cat("\n")
}
sink()

# Write a covariate matrix

sink(paste0(d, "test_data.cov"))
for(i in 1:length(labels)) {
  cat(paste(labels[i], paste(rnorm(10), collapse="\t"), sep="\t"))
  cat("\n")
}
sink()


### TEST
if(sys.nframe() == 0L) {
  library(lme4)
  print("Running tests.")
  
  setwd("~/CLionProjects/Permute_Associate/")
  
  df <- read.table(paste0(d, "test_data.txt"), header=T, comment.char = "~")
  cov <- read.table(paste0(d, "test_data.cov"), header=F)
  ped <- read.table(paste0(d, "test_data.ped"), header=F)
  
  case_control = ped[,6]
  for(i in 1:length(case_control)) {
    if(case_control[i] == 2) {
      case_control[i] = 1
    } else {
      case_control[i] = 0
    }
  }
  
  cov$V1 = rep(1, length(case_control))
  
  for(gene in levels(df$X.Gene)) {
    s = df[df$X.Gene == gene,]
    s$Transcript = factor(s$Transcript)
    
    for(ts in levels(s$Transcript)) {
      Z = t(as.matrix(s[s$Transcript == ts,4:ncol(s)]))
      Z[Z > 2] = 0
      Z <- Z[, colSums(Z) > 0]
      p <- colMeans(Z) / 2
      Z.sw <- sweep(Z, 2, 2 * p)
      Z.sc <- sweep(Z.sw, 2, 2 * p * (1 - p), FUN = "/")
      grm <- Z.sc %*% t(Z.sw) / ncol(Z)
      new_df <- cbind(case_control, cov, Z, grm)
      # This doesn't work. You have to use the GRM to specify the variance structure. Apparently not
      # possible in lmer. 
      f1 <- formula(paste("case_control ~ . + ( 0 |", paste(colnames(grm), collapse = "+"), ")"))
      res <- glmer(f1, data = new_df, family = binomial)
    }
  }
}
