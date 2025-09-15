# Title     : TODO
# Objective : TODO
# Created by: rjbohlender
# Created on: 8/28/20

suppressPackageStartupMessages(library(getopt))

spec <- matrix(c(
  'matrix', 'm', 1, 'character',
  'output', 'o', 1, 'character',
  'transcript', 't', 1, 'character',
  'help', 'h', 0, 'logical'
), byrow = T, nrow = 4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=T))
  q(status = 1)
}

# df <- read.table(opt$matrix, header=T)
df <- read.table('/Volumes/epi/huff/jchen15/Jointcall0520/xpat_run_temps/BBA_data/final_MDA_PAN/MDA_PAN.case_control.new.sort2.txt', header=T)

df[df$Gene == 'ATM',]



