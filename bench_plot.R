library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

df <- data.frame(
  Tool = c(rep("CARVA", 7), rep("RVTests", 7)),
  Method = rep(c("VAAST", "WSS", "RVT1", "SKAT", "SKATO", "CMC", "VT"), 2),
  Time = c(338, 1822, 75, 212, 12770, 449, 399, NA, 68613, 1756, 7130, 13849, 1876, 15227)
  #cputime = length(c(7549, 40692, 553, 553, 2798, 339677, 10254, 9050, NA, NA, 1562, NA, 30108, 47467, 1706, 15091))
)

right_label <- df %>%
  group_by(Method) %>%
  arrange(desc(Time)) %>%
  top_n(1)

left_label <- df %>%
  group_by(Method) %>%
  arrange(desc(Time)) %>%
  slice(2)

ggplot(df, aes(x=Time, y=Method)) +
  geom_point(size = 2.5, aes(color=Tool)) + 
  geom_line(aes(group=Method)) + 
  geom_text(data=right_label, aes(color=Tool, label = round(Time, 0)), size = 3, hjust = -.5) +
  geom_text(data=left_label, aes(color=Tool, label = round(Time, 0)), size = 3, hjust = 1.5) +
  scale_x_continuous(limits = c(0, max(df$Time, na.rm=T) + 0.1 * max(df$Time, na.rm=T))) +
  theme_minimal() +
  ggtitle("Runtime in seconds of CARVA vs. RVTest for 2000 samples")
# ggplot(df, aes(x=cputime, y=method, col=Tool)) + geom_point() + ggtitle("CPU time in seconds of PA vs. RVTests for 2000 samples")

#######
# 1 Million Samples
#######

df <- data.frame(
  Tool = c(rep("CARVA", 7), rep("RVTests", 7)),
  Method = rep(c("VAAST", "WSS", "RVT1", "SKAT", "SKATO", "CMC", "VT"), 2),
  Time = c(115398, 525474, 18337, 279385, 43653, 322916, 98245, NA, 864006, 21868, NA, 189022, 89743, 50572) / 60 / 60
  #cputime = length(c(7549, 40692, 553, 553, 2798, 339677, 10254, 9050, NA, NA, 1562, NA, 30108, 47467, 1706, 15091))
)

right_label <- df %>%
  group_by(Method) %>%
  arrange(desc(Time)) %>%
  top_n(1)

left_label <- df %>%
  group_by(Method) %>%
  arrange(desc(Time)) %>%
  slice(2)

ggplot(df, aes(x=Time, y=Method)) +
  geom_point(size = 2.5, aes(color=Tool)) + 
  geom_line(aes(group=Method)) + 
  geom_text(data=right_label, aes(color=Tool, label = round(Time, 0)), size = 3, hjust = -.5) +
  geom_text(data=left_label, aes(color=Tool, label = round(Time, 0)), size = 3, hjust = 1.5) +
  scale_x_continuous(limits = c(0, max(df$Time, na.rm=T) + 0.1 * max(df$Time, na.rm=T))) +
  theme_minimal() +
  ggtitle("Runtime in hours of CARVA vs. RVTest for 1,000,000 samples")
