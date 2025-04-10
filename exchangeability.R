library(AssotesteR)
# Test exchangeability
df <- read.table('/Users/rjbohlender/Downloads/MDA_Ovarian/brca1.txt', header=T, check.names = F)
brca1.dat <- t(as.matrix(df[,4:length(df)]))
brca1.vars <- dim(brca1.dat)[2]
brca1.dat[brca1.dat == 9] = 0
# brca1.dat <- rowSums(brca1.dat)

dat <- brca1.dat
dat <- dat[order(rownames(dat)), ]

ped <- read.table('/Users/rjbohlender/Downloads/MDA_Ovarian/xqc.new.ped', header=F)
ped$V2 <- make.names(ped$V2)
ped <- ped[order(ped$V2), ]

ped$V6[ped$V6 == 1] = 0
ped$V6[ped$V6 == 2] = 1

Y = ped$V6

covars <- read.table('/Users/rjbohlender/Downloads/MDA_Ovarian/samples_pca.matrix.txt', header=F)
covars <- as.matrix(covars[order(covars[,1]),2:11])

perm <- sample.int(length(ped$V6), replace = F)

glm(Y ~ covars, data=df, family=binomial())

# Get minor and major allele carriers
mac <- which(rowSums(dat) > 0)
maj <- which(rowSums(dat) == 0)

altY <- Y
altY[maj] <- sample(Y[maj])

WSS(Y, dat)
# wss.stat 1291618
WSS(altY, dat)
# wss.stat 1291618 -- Shuffling maj allele carriers has no effect
CMC(Y, dat, maf = 0.01)
# cmc.stat 9.985264
CMC(altY, dat, maf = 0.01)
# cmc.stat 9.985264
CALPHA(Y, dat)
# calpha.stat 33.053950
CALPHA(altY, dat)
# calpha.stat 33.053950
RVT1(Y, dat)
# rvt1.stat 3.829717
RVT1(altY, dat)
# rvt1.stat 3.829717
RVT2(Y, dat)
# rvt2.stat 4.837249
RVT2(altY, dat)
# rvt2.stat 4.837249
VT(Y, dat)
# vt.stat 2.215397
VT(altY, dat)
# vt.stat 2.215397
