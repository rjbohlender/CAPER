library(Hotelling)

set.seed(1)

# Sample allele frequencies from the following
shape1 = 1
shape2 = 7
nvars = 10
nsamples = 1000

curve(dbeta(x, shape1, shape2))
hist(rbeta(10000, shape1, shape2))

freq = rbeta(nvars, shape1, shape2)

test_data = NULL
for(i in 1:nvars) {
  test_data = cbind(test_data, rbinom(nsamples, 2, freq[i]))
}

# Check frequency distribution
maf = colMeans(test_data) / 2
maf
freq

phen = rbinom(nsamples, 1, 0.5)

rare = maf < 0.005


s1 = test_data[which(phen == 1),]
s2 = test_data[which(phen == 0),]

fit = hotelling.test(s1, s2)
fit

# Print data in format to include in C++ test code
for(i in 1:nsamples) {
  cat(paste('{', paste(paste(test_data[i,], sep='', collapse = ","), '},\n')))
}

cat(paste(phen, sep='', collapse = ", "))

