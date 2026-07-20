###############################
#  Correlation: SDW vs log2_comp
#  File: coorelation_sdw_comp_jake.xlsx
###############################

## 0. Install required packages (run once if needed)
# install.packages("readxl")
# install.packages("ggplot2")

library(readxl)
library(ggplot2)

###############################
# 1. Read the Excel file
###############################

# Set working directory to the folder where your Excel file is
setwd("/Users/taniagupta/Desktop/genomes_paper") # <-- EDIT THIS PATH

dat <- read_excel("coorelation_sdw_comp_jake.xlsx")

# Quick check
head(dat)
str(dat)

# At this point you should see columns:
# Strain, log2_comp, sdw

###############################
# 2. Normality checks
###############################

## 2A. Histograms + Q–Q plots for log2_comp
par(mfrow = c(1, 2))  # two plots in one row

hist(dat$log2_comp,
     main = "Histogram of log2_comp",
     xlab  = "log2_comp")

qqnorm(dat$log2_comp,
       main = "Q-Q plot of log2_comp")
qqline(dat$log2_comp, col = "orange")

par(mfrow = c(1, 1))  # reset layout

## 2B. Histograms + Q–Q plots for sdw
par(mfrow = c(1, 2))

hist(dat$sdw,
     main = "Histogram of SDW",
     xlab  = "SDW")

qqnorm(dat$sdw,
       main = "Q-Q plot of SDW")
qqline(dat$sdw, col = "orange")

par(mfrow = c(1, 1))

## 2C. Shapiro–Wilk normality tests
shapiro_log2 <- shapiro.test(dat$log2_comp)
shapiro_sdw  <- shapiro.test(dat$sdw)

shapiro_log2
shapiro_sdw

# Look at the p-values:
# p < 0.05  -> not normal
# p >= 0.05 -> approximately normal

###############################
# 3. Correlation analyses
###############################

## 3A. Spearman (correct test for non-normal data)
spearman_res <- cor.test(dat$log2_comp, dat$sdw,
                         method = "spearman")

spearman_res
# spearman_res$estimate  # rho value
# spearman_res$p.value   # p-value

## 3B. Pearson (for comparison only; assumptions not met)
pearson_res <- cor.test(dat$log2_comp, dat$sdw,
                        method = "pearson")

pearson_res

###############################
# 4. Scatter plot with regression line
###############################

# Base R version
plot(dat$log2_comp, dat$sdw,
     xlab = "log2 competitiveness",
     ylab = "Shoot dry weight (SDW)",
     main = "SDW vs log2 competitiveness")

abline(lm(sdw ~ log2_comp, data = dat),
       col = "orange", lwd = 2)

# ggplot2 version (nicer)
ggplot(dat, aes(x = log2_comp, y = sdw)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "log2 competitiveness",
       y = "Shoot dry weight (SDW)",
       title = "SDW vs log2 competitiveness") +
  theme_bw()

###############################
# 5. Print a compact summary
###############################

cat("\nShapiro–Wilk normality tests:\n")
print(shapiro_log2)
print(shapiro_sdw)

cat("\nSpearman correlation (main result):\n")
print(spearman_res)

cat("\nPearson correlation (for reference only):\n")
print(pearson_res)
pdf("normality_plots.pdf", width = 8, height = 6)

# Plot 1
hist(dat$log2_comp, main="Histogram log2_comp")

# Plot 2
qqnorm(dat$log2_comp, main="Q-Q log2_comp")
qqline(dat$log2_comp, col="orange")

# Plot 3
hist(dat$sdw, main="Histogram SDW")

# Plot 4
qqnorm(dat$sdw, main="Q-Q SDW")
qqline(dat$sdw, col="orange")

dev.off()

png("scatter_sdw_vs_comp.png", width = 1200, height = 900, res = 150)
plot(dat$log2_comp, dat$sdw,
     xlab="log2_comp",
     ylab="SDW",
     main="SDW vs log2_comp")
abline(lm(sdw ~ log2_comp, data = dat), col="orange", lwd=2)
dev.off()


