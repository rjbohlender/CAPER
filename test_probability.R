
setwd("~/ClionProjects/Permute_Associate/")

s1pset <- read.table("db/stage1_pset.txt")
prob <- read.table("db/stage1_pset.txt.lr")
s2pset <- read.table("db/stage2_pset.txt")
gMpset <- read.table("db/gMean_pset.txt")
rmspset <- read.table("db/rms_pset.txt")

s1pset_prob <- as.matrix(colMeans(s1pset))
s2pset_prob <- as.matrix(colMeans(s2pset))
prob <- as.matrix(prob)
gMean_prob <- colMeans(gMpset)
rms_prob <- colMeans(rmspset)

idx = which(s2pset_prob != 1)
idx = intersect(idx, which(s2pset_prob != 0))

paste(s1pset_prob[idx], s2pset_prob[idx], prob[idx])

sqrt(mean((s2pset_prob[idx] - prob[idx])^2))
sqrt(mean((s1pset_prob[idx] - prob[idx])^2))

cd63idx = c(1, 53, 56, 228, 248, 265, 292, 457, 477, 694, 707, 775, 781, 830)
repIdx = c(1, 56, 228, 248, 265, 292, 457, 477, 694, 775, 781, 1134, 1788, 1911)

paste(s1pset_prob[cd63idx], s2pset_prob[cd63idx], prob[cd63idx])
paste(s1pset_prob[repIdx], s2pset_prob[repIdx], prob[repIdx])

sqrt(mean((s2pset_prob[cd63idx] - prob[cd63idx])^2))
sqrt(mean((s1pset_prob[cd63idx] - prob[cd63idx])^2))

sqrt(mean((s2pset_prob[repIdx] - prob[repIdx])^2))
sqrt(mean((s1pset_prob[repIdx] - prob[repIdx])^2))

# Compare P-Value distributions
s1 <- read.table('test_s1/VAAST.simple', header=T)
s2 <- read.table('test_s2/VAAST.simple', header=T)
