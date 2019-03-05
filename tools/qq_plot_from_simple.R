#!/usr/bin/env Rscript
# Title     : Create a QQ plot from the simple file output of a run
# Created by: rjbohlender
# Created on: 2019-01-30

library(Haplin)

args <- commandArgs(trailingOnly=T)

if(length(args) != 2) {
    stop("Usage: qq_plot_from_simple.R <simple_file> <pdf_output>")
}

df <- read.table(args[1], header=T)

pdf(args[2], 6.5, 6.5)

pvals <- df$Empirical_MidP
names(pvals) <- df$Gene

pQQ(pvals[!duplicated(df$Gene)])
dev.off()
