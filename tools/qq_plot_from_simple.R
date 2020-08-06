#!/usr/bin/env Rscript
# Title     : Create a QQ plot from the simple file output of a run
# Created by: rjbohlender
# Created on: 2019-01-30

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(Haplin))

spec <- matrix(c(
  'simple', 's', 1, 'character',
  'output', 'o', 1, 'character',
  'test_statistic', 't', 0, 'logical',
  'help', 'h', 0, 'logical'
), byrow = T, nrow = 4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=T))
  q(status = 1)
}

df <- read.table(opt$simple, header=T)

pdf(opt$output, 6.5, 6.5)

if (opt$test_statistic) {
   pvals <- df$Test_Statistic
} else {
  pvals <- df$Empirical_MidP
}
names(pvals) <- df$Gene

pQQ(pvals[!duplicated(df$Gene)], nlabs=5)
dev.off()
