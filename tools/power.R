# Title     : Plot Power Analysis Results
# Objective : Produce a usable plot of the power analysis results
# Created by: rjbohlender
# Created on: 2019-04-01

library(ggplot2)
library(gridExtra)
library(gtable)

args <- commandArgs(T)

if(length(args) != 2) {
    stop("Usage: power.R <power_file> <pdf_output>")
}

df <- read.table(args[1], header=T)

g <- NULL
pdf(args[2], 7, 5)
for(gene in levels(df$Gene)) {
    subs <- subset(df, Gene == gene)
    p <- ggplot(data=subs, aes(x=Alpha, y=Ratio, col=interaction(Ncases, Ncontrols))) +
        geom_line() +
        geom_point() +
        facet_wrap(Transcript ~ Gene)
    print(p)
}
dev.off()


