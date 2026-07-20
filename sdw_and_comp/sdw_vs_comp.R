# ================================
# Competition vs SDW Scatter Plot and corelation between the two traits
# ================================

# Install if needed
# install.packages(c("readxl", "ggplot2", "ggrepel"))

library(readxl)
library(ggplot2)

# ---- Read Excel file ----
file_path <- "/Users/taniagupta/Desktop/genomes_paper/sdw_vs_comp/sdw_vs_comp.xlsx"
dat <- read_excel(file_path)

# ---- Clean column names ----
colnames(dat) <- make.names(colnames(dat))

# ---- Spearman correlation ----
cor_res <- cor.test(dat$log2_comp_jake_analysis, dat$sdw, method = "spearman")
print(cor_res)

# ---- Plot (no labels + confidence band) ----
p <- ggplot(dat, aes(x = log2_comp_jake_analysis, y = sdw)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.25) +  # CI band on
  labs(
    x = "Competition (log2 transformed)",
    y = "Shoot Dry Weight (SDW)",
    title = "Competition vs Symbiotic Effectiveness",
    subtitle = paste0("Spearman rho = ", round(cor_res$estimate, 3),
                      ", p = ", signif(cor_res$p.value, 3))
  ) +
  theme_classic(base_size = 14)

print(p)

# ---- Save figure ----
ggsave("Competition_vs_SDW_scatter.pdf", plot = p, width = 7, height = 6)
ggsave("Competition_vs_SDW_scatter.png", plot = p, width = 7, height = 6, dpi = 300)



#--------------------------------------------------------------------------------------------------
## to make the quadrant plot based on the median value of each trait. 

# install.packages(c("readxl","dplyr","ggplot2"))
library(readxl)
library(dplyr)
library(ggplot2)

file_path <- "/Users/taniagupta/Desktop/genomes_paper/sdw_vs_comp/sdw_vs_comp.xlsx"   # <-- change
dat <- read_excel(file_path)

# standardize names
colnames(dat) <- make.names(colnames(dat))

# keep only needed columns & drop missing
dat <- dat %>%
  transmute(
    Strain = as.character(Strain),
    comp = as.numeric(log2_comp_jake_analysis),
    sdw  = as.numeric(sdw)
  ) %>%
  filter(!is.na(comp) & !is.na(sdw))


x_med <- median(dat$comp, na.rm = TRUE)
y_med <- median(dat$sdw, na.rm = TRUE)

dat_quad <- dat %>%
  mutate(
    comp_cat = ifelse(comp >= x_med, "High competition", "Low competition"),
    sdw_cat  = ifelse(sdw  >= y_med, "High SDW", "Low SDW"),
    quadrant = factor(paste(comp_cat, sdw_cat, sep = " / "),
                      levels = c("Low competition / Low SDW",
                                 "Low competition / High SDW",
                                 "High competition / Low SDW",
                                 "High competition / High SDW"))
  )

p2 <- ggplot(dat_quad, aes(x = comp, y = sdw)) +
  geom_point(aes(color = quadrant), size = 2.8, alpha = 0.85) +
  geom_vline(xintercept = x_med, linetype = "dashed") +
  geom_hline(yintercept = y_med, linetype = "dashed") +
  labs(
    x = "Competition (log2 transformed)",
    y = "Shoot Dry Weight (SDW)",
    color = "Quadrant",
    title = "Competition vs Symbiotic Effectiveness"
  ) +
  theme_classic(base_size = 14)

print(p2)

ggsave("Competition_SDW_quadrants.pdf", p2, width = 7, height = 5.5)
ggsave("Competition_SDW_quadrants.png", p2, width = 7, height = 5.5, dpi = 300)