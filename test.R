library(Haplin)
library(statmod)

df <- read.table('~/Downloads/MDA_Ovarian/results/RVT1.simple', header=T)
pvals <- as.vector(df$Empirical_P)
pQQ(pvals)
pvals <- as.vector(df$Empirical_MidP)
pQQ(pvals)
pvals <- as.vector(df$MGIT)
pQQ(pvals)
pvals <- as.vector(df$Original)
pQQ(pvals)

df$Gene[pvals < 0.0005]

bw <- bw.SJ(df$Original)
plot(density(df$Original, bw=bw, kernel="gaussian"))

bw <- bw.SJ(df$Empirical_P)
plot(density(df$Empirical_P, bw=bw, kernel="gaussian"))

########
library(AssotesteR)

df <- read.table('~/PycharmProjects/Permute_Associate/test/input.select.sort.matrix', header=T)

bard1 = t(as.matrix(df[df[,1] == 'BARD1', 4:ncol(df)]))
bard1[bard1 > 2] = 0
cov <- read.table('~/PycharmProjects/Permute_Associate/test/vaast-covinput.txt', header=F)
phenotypes = as.vector(cov$V1)

CMC(phenotypes, bard1, maf=0.001)
SKAT(phenotypes, bard1, kernel = "linear")

brca1 = df[df[,1] == 'BRCA1',]
brca1 = t(as.matrix(brca1[brca1[,2] == 'NM_007294', 4:ncol(df)]))
brca1[brca1 > 2] = 0
CMC(phenotypes, brca1, maf=0.001)
SKAT(phenotypes, brca1, kernel = "linear")

#########
# Test Cases
#########
library(AssotesteR)

X1 = matrix(c(0, 1, 0, 1, 1, 1, 1, 1, 2, 0, 0, 0, 2, 0, 0, 1), nrow=4)
X2 = matrix(c(0, 1, 0, 1, 1, 1, 0, 1, 2, 0, 0, 0, 2, 0, 0, 1), nrow=4)
Y = c(1, 1, 0, 0)

CALPHA(Y, X1)
CALPHA(Y, X2)

CMC(Y, X1)
CMC(Y, X2)

SKAT(Y, X1)
SKAT(Y, X2)

VT(Y, X1)
VT(Y, X2)

WSS(Y, X1)
WSS(Y, X2)

########
# Test SKAT
########

library(SKAT)

# df <- read.table("~/PycharmProjects/Permute_Associate/test/input.select.sort.matrix", header=T)
# covars <- read.table("~/PycharmProjects/Permute_Associate/test/vaast-covinput.txt", header=F)
# df <- read.table("~/Downloads/tcga-ov/test3.matrix.txt", header=T)
# covars <- read.table("~/Downloads/tcga-ov/ov.pca.cov", header=F)

df <- read.table('~/Downloads/tcga-ov/1000genes.matrix.txt', header=T)
covars <- read.table('~/Downloads/tcga-ov/ov.pca.cov', header=F)

obj <- SKAT_Null_Model(formula("V1 ~ 1 + ."), data=covars, out_type = "D", n.Resampling = 1, type.Resampling = "permutation", Adjustment = T)

for(gene in levels(df$Gene)) {
  s = df[df$Gene == gene,]
  s$Transcript = factor(s$Transcript)
  
  for(ts in levels(s$Transcript)) {
    Z = t(as.matrix(s[s$Transcript == ts,4:ncol(s)]))
    Z[Z > 2] = 0
    res <- SKAT(Z, obj, kernel="linear", method="liu")
    msg <- sprintf("%s %s %6.5f", gene, ts, res$p.value)
    print(msg)
  }
}

######
# More SKAT Test
######
library(SKAT)

df <- read.table('~/Downloads/tcga-ov/test3.matrix.txt', header=T)
covars <- read.table('~/Downloads/tcga-ov/ov.pca.cov', header=F)

obj <- SKAT_Null_Model(formula("V1 ~ 1 + V2 + V3"), data=covars, out_type = "D", n.Resampling = 0, type.Resampling = "permutation", Adjustment = T)

Z = t(as.matrix(df[df$Gene == "BRCA2",4:ncol(df)]))
Z[Z > 2] = 0

res <- SKAT(Z, obj, kernel="linear", method="optimal.adj")


### From biostatistics and ajhg papers
library(expm)

df <- read.table('~/Downloads/tcga-ov/test3.matrix.txt', header=T)
covars <- read.table('~/Downloads/tcga-ov/ov.pca.cov', header=F)

G = t(as.matrix(df[df$Gene == "BRCA2",4:ncol(df)]))
G[G > 2] = 0

n = dim(G)[1] # Samples
p = dim(G)[2] # Variants

X1 = model.matrix(formula("V1 ~ 1 + V2 + V3"), data = covars)
glmfit = glm(formula("V1 ~ 1 + V2 + V3"), data = covars, family=binomial(link=logit))
mu = glmfit$fitted.values
eta = glmfit$linear.predictors
delta = diag(eta)
pi_1 = mu * (1 - mu) # Variance of a bernnouli random variable

res = covars$V1 - mu

V = diag(pi_1)

r.all = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 0.999)
n.r = length(r.all)
W = list()
R = list()
Kr = list()
Q = list()
lambda = list()

invV = diag(1 / pi_1) # Inverse of a diagonal matrix is just the inverse of the elements
P = invV - invV %*% X1 %*% solve(t(X1) %*% invV %*% X1) %*% t(X1) %*% invV
rtV = sqrt(V) # Square root of diagonal matrix is just the sqrt of the elements
rtP = sqrtm(P)
for(i in 1:n.r) {
  W[[i]] = diag(p)
  R[[i]] = (1 - r.all[i]) * diag(p) + rsolve.all[i] * matrix(rep(1, p * p), nrow=p)
  Kr[[i]] = G %*% (W[[i]] %*% (R[[i]] %*% (W[[i]] %*% t(G))))
  Q[[i]] = res %*% delta %*% invV %*% Kr[[i]] %*% invV %*% delta %*% res
  lambda[[i]] = eigen(rtP %*% Kr[[i]] %*% rtP)$values
}

out = list()
for(i in 1:n.r) {
  r = length(lambda[[i]])
  nc1 = rep(0, r)
  df1 = rep(1, r)
  
  out[[i]] = .C("qfc",
                lambdas=as.double(lambda[[i]]),
                noncentral=as.double(nc1),
                df=as.integer(df1),
                r=as.integer(r),
                sigma=as.double(0),
                q=as.double(Q[[i]]),
                lim=as.integer(10000),
                acc=as.double(1e-6),
                trace=as.double(rep(0,7)),
                ifault=as.integer(0),
                res=as.double(0),
                PACKAGE="SKAT")
  print(1-out[[i]]$res)
}


# Following Kai Wang (2016)
# Calculating SKAT

df <- read.table('~/Downloads/tcga-ov/test3.matrix.txt', header=T)
covars <- read.table('~/Downloads/tcga-ov/ov.pca.cov', header=F)

ff = formula("V1 ~ 1 + V2 + V3")
X1 = model.matrix(ff, data=covars)
glmfit = glm(ff, data=covars, family=binomial(link=logit))

pi_1 = glmfit$fitted.values * (1 - glmfit$fitted.values)

G = t(df[,4:ncol(df)])
G[G > 2] = 0

V = diag(pi_1)

P = V - V %*% X1 %*% solve(t(X1) %*% V %*% X1) %*% t(X1) %*% V

n = dim(G)[1]
p = dim(G)[2]

sigma = 1 / (n - p - 1) * (t(G) %*% P %*% G)

### SKATR

library(SKATR)
library(CompQuadForm)

df <- read.table('~/Downloads/tcga-ov/1000genes.matrix.txt', header=T)
covars <- read.table('~/Downloads/tcga-ov/ov.pca.cov', header=F)

ff = formula('V1 ~ 1 + .')
design = model.matrix(ff, data=covars)
obj = KAT.null(covars$V1, design)

for(gene in levels(df$Gene)) {
  s = df[df$Gene == gene,]
  s$Transcript = factor(s$Transcript)
  
  for(ts in levels(s$Transcript)) {
    Z = t(as.matrix(s[s$Transcript == ts,4:ncol(s)]))
    Z[Z > 2] = 0
    res <- SKATh(obj = obj, G = Z)
    msg <- sprintf("%s %s %6.5f", gene, ts, res)
    print(msg)
  }
}


######
# Logistic Regression
######

library(arm)
df <- read.table('~/Downloads/MDA_Ovarian/case_control.report.matrix.sort.txt', header=T, nrows=10000)
covars <- read.table('~/Downloads/MDA_Ovarian/samples_pca.matrix.txt', header=F, row.names = 1)

ped <- read.table('~/Downloads/MDA_Ovarian/xqc.new.ped', header=F)

case_control = sort(ped[,6], T)
for(i in 1:length(case_control)) {
  if(case_control[i] == 2) {
    case_control[i] = 1
  } else {
    case_control[i] = 0
  }
}

ff <- case_control ~ 1 + .

frame <- model.frame(ff, covars)

for(gene in levels(df$Gene)) {
  if(gene != "A2ML1") {
    next
  }
  s = df[df$Gene == gene,]
  s$Transcript = factor(s$Transcript)
  
  for(ts in levels(s$Transcript)) {
    Z = t(s[s$Transcript == ts,4:ncol(s)])
    Z[Z > 2] = 0
  }
  break
}

scale.frame = scale(frame[,2:ncol(frame)], center = T, scale = T)
merged = data.frame(cbind(case_control, rowSums(Z), scale.frame))

fit.1 = bayesglm(case_control ~ ., data=merged, family = binomial(link="logit"))
exp(coef(fit.1))
display(fit.1)
coef(fit.1)

### RVT test

#!/usr/bin/env Rscript
require(data.table)

eprintf <- function(...) {
  cat(paste(sprintf(...), "\n"), sep='', file=stderr())
  flush.console()
}

# args = commandArgs(trailingOnly = T)
setwd("~/CLionProjects/Permute_Associate/")
# args = c("test.sim/test_data.txt", "test.sim/test_data.cov", "test.sim/test_data.ped")
args = c("~/Downloads/MDA_Ovarian/case_control.report.matrix.sort.txt", "~/Downloads/MDA_Ovarian/samples_pca.matrix.txt", "~/Downloads/MDA_Ovarian/xqc.new.ped")

matrix_file = args[1]
covar_file = args[2]
ped_file = args[3]

eprintf("Loading matrix.")
df <- fread(matrix_file, header=T)
eprintf("Loaded matrix.")
cov <- fread(covar_file, header=F)
eprintf("Loaded covariates.")
ped <- fread(ped_file, header=T)
eprintf("Loaded ped.")

df <- as.data.frame(df)
colnames(df)[1] = "Gene"
colnames(df)[2] = "Transcript"
colnames(df)[3] = "Location"
df$Gene = as.factor(df$Gene)
cov <- as.data.frame(cov)
ped <- as.data.frame(ped)

names <- colnames(df)[4:ncol(df)]
name_order <- NULL
for(n in names) {
  name_order <- c(name_order, which(ped[,2] == n))
}

case_control = ped[name_order, 6]
for(i in 1:length(case_control)) {
  case_control[i] = ifelse(case_control[i] == 2, 1, 0)
}

name_order <- NULL
for(n in names) {
  name_order <- c(name_order, which(cov$V1 == n))
}

cov$V1 = case_control

for(i in 2:ncol(cov)) {
  cov[,i] = cov[name_order,i]
}

obj <- KAT.null(as.vector(cov$V1), as.matrix(cov[,2:ncol(cov)]))

for(gene in levels(df$Gene)) {
  if(gene != "BRCA1") {
    next;
  }
  s = df[df$Gene == gene,]
  s$Transcript = factor(s$Transcript)
  
  for(ts in levels(s$Transcript)) {
    Z = t(as.matrix(s[s$Transcript == ts,4:ncol(s)]))
    Z[Z > 2] = 0
    
    design <- cbind(1, as.matrix(cov[,2:ncol(cov)]))
    
    fit.1 <- glm(case_control ~ design, family = binomial("logit"))
    
    r <- rowSums(Z) / dim(Z)[2]
    
    fit.2 <- glm(case_control ~ design + r, family = binomial("logit"))
  }
}


