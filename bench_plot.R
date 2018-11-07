library(ggplot2)

df <- data.frame(
  tool = c(rep("PA", 8), rep("RVTests", 8)),
  method = rep(c("VAAST", "BURDEN", "RVT1", "RVT2", "SKAT", "SKATO", "CMC", "VT"), 2),
  time = c(338, 1513, 75, 76, NA, NA, 449, 399, NA, NA, 1756, NA, 7130, 13849, 1876, 15227),
  cputime = c(7549, 40692, 553, 553, NA, NA, 10254, 9050, NA, NA, 1562, NA, 30108, 47467, 1706, 15091)
)

ggplot(df, aes(x=time, y=method, col=tool)) + geom_point() + ggtitle("Runtime in seconds of PA vs. RVTest for 2000 samples")
ggplot(df, aes(x=cputime, y=method, col=tool)) + geom_point() + ggtitle("CPU time in seconds of PA vs. RVTests for 2000 samples")
