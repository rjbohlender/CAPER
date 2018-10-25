# Quick and dirty genetic simulation

sink("~/CLionProjects/Permute_Associate/test.sim/test.matrix")
total <- 500000

# Log odds 
lodds <- log(10)

ngenes <- 1000
nsnps <- rpois(ngenes, 20)
ncausal <- 10

freq <- list()
for(i in 1:ngenes) {
  freq[[i]] <- runif(nsnps[i], 5 / total, 65 / total)
}


data = matrix(0L, nrow = total, ncol = ngenes)
location = 0
for(i in 1:ngenes) {
  gene = matrix(rbinom(total * nsnps[i], 2, freq[[i]]), nrow = nsnps[i], ncol = total)
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

betas <- c(rep(lodds, ncausal), rep(log(0.1) / (ngenes - ncausal), ngenes - ncausal))
case_control <- rbinom(total, 1, binomial()$linkinv(data %*% betas))

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

sink("~/CLionProjects/Permute_Associate/test.sim/header.txt")
cat(header)
sink()

# Write a ped file

sink("~/CLionProjects/Permute_Associate/test.sim/test_data.ped")
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

sink("~/CLionProjects/Permute_Associate/test.sim/test_data.cov")
for(i in 1:length(labels)) {
  cat(paste(labels[i], paste(rnorm(10), collapse="\t"), sep="\t"))
  cat("\n")
}
sink()
