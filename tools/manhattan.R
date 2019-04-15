# Title     : TODO
# Objective : TODO
# Created by: rjbohlender
# Created on: 2019-02-15

library(qqman)

args <- commandArgs(T)

if(length(args) != 3) {
    stop("Usage: manhattan.R <qqman_file> <pdf_output> <title>")
}

df <- read.table(args[1], header=T)

pdf(args[2], 6.5, 6.5)
manhattan(df, main = args[3], annotatePval = 0.001)
dev.off()
