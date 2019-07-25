# Title     : TODO
# Objective : TODO
# Created by: rjbohlender
# Created on: 2019-07-22

library(Haplin)

args <- commandArgs(trailingOnly=T)

if(length(args) != 2) {
    stop("Usage: qq_plot_from_simple.R <simple_file> <pdf_output>")
}

df <- read.table(args[1], header=T)

pdf(args[2], 6.5, 6.5)

pvals <- df$Test_Statistic
names(pvals) <- df$Gene

pQQ(pvals[!duplicated(df$Gene)])
dev.off()
