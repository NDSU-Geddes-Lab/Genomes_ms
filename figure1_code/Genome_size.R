
#For generating the dot plot of size of each replicon
# install.packages(c("readxl","tidyverse"))   # if needed

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


# 1) Let the user pick the file
file_path <- file.choose()

# 2) Read in your Excel (or switch to read.csv if CSV)
df <- read_excel(file_path)

# 3) Pivot to “long” form, dropping any NA sizes
df_long <- df %>%
  pivot_longer(
    cols         = -`Sample ID`,
    names_to     = "Replicon",
    values_to    = "Size",
    values_drop_na = TRUE
  )

# 4) Define the exact order of your replicons
rep_levels <- c(
  "Chromosome",
  "pSymA",
  "pSymB",
  "Chromosome_psymB_integrated",
  "pSymA_pSymB_integrated",
  "Accessory_plasmid"
)
pal <- c(
  "Chromosome"                  = "#FF6F61",
  "pSymA"                       = "#7FB069",
  "pSymB"                       = "#4CAEA8",
  "Chromosome_psymB_integrated" = "#FFA500",
  "pSymA_pSymB_integrated"      = "#FF69B4",
  "Accessory_plasmid"        = "#C77CFF")
  
# 5) And the 2‑line wrapped labels you want on the x‑axis
wrapped_labels <- c(
  "Chromosome\n",
  "pSymA\n",
  "pSymB\n",
  "Chromosome\npsymB integrated",
  "pSymA\npsymB integrated",
  "Accessory\nplasmid"
)

# 6) Turn Replicon into a factor with that order
df_long$Replicon <- factor(df_long$Replicon, levels = rep_levels)

p <- ggplot(df_long, aes(x = Replicon, y = Size, color = Replicon)) +
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.7) +
  scale_color_manual(values = pal, guide = "none") +   # <- change here
  scale_x_discrete(labels = wrapped_labels) +
  scale_y_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) +
  labs(title = "Replicon Sizes Across Samples", x = "Replicon", y = "Size (kbp)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 8), size = 12, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 8), size = 12, face = "bold")
  )


print(p)

#For creating bar graphh of strains contianing replicon

# Required packages
# install.packages(c("readxl","tidyr","dplyr","ggplot2","scales"))
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)

# 1) Pick the Excel file interactively
xlsx_path <- file.choose()   # opens a file picker

# 2) Read the sheet
sheet_name <- "No_of_replicon_2"
df <- read_excel(xlsx_path, sheet = sheet_name)

# If your sheet is a single row of counts with headers as categories,
# pivot it to long format:
df_long <- df %>%
  pivot_longer(everything(),
               names_to = "Replicon",
               values_to = "n") %>%
  mutate(Replicon = as.character(Replicon))

# Keep bars in the same order as the columns in Excel
replicon_order <- colnames(df)
df_long$Replicon <- factor(df_long$Replicon, levels = replicon_order)

# 3) Color mapping (exactly as you specified)
pal <- c(
  "Chromosome"                  = "#FF6F61",
  "pSymA"                       = "#7FB069",
  "pSymB"                       = "#4CAEA8",
  "Chromosome_psymB_integrated" = "#FFA500",
  "pSymA_pSymB_integrated"      = "#FF69B4",
  "0_Accessory_plasmids"        = "#C77CFF",
  "1_Accessory_plasmid"         = "#C77CFF",
  "2_Accessory_plasmids"        = "#C77CFF",
  "3_Accessory_plasmids"        = "#C77CFF",
  "4_Accessory_plasmids"        = "#C77CFF"
)

# 5) And the 2‑line wrapped labels you want on the x‑axis
wrapped_labels <- c(
  "Chromosome\n",
  "pSymA\n",
  "pSymB\n",
  "Chromosome\npsymB integrated",
  "pSymA\npsymB integrated",
  "0_Accessory\nplasmid",
  "1_Accessory\nplasmid",
  "2_Accessory\nplasmid",
  "3_Accessory\nplasmid",
  "4_Accessory\nplasmid"
  )

# 4) Plot
p <- ggplot(df_long, aes(x = Replicon, y = n, fill = Replicon)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.6) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = wrapped_labels) +
  scale_y_continuous(
    limits = c(0, 200),           # set y-axis min and max
    expand = expansion(mult = c(0, 0.05)), 
    breaks = seq(0, 200, 50)      # nice breaks at 0, 50, 100, 150, 200
  ) +
  labs(x = "Replicon", y = "No. of strains") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 8), size = 12, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 8), size = 12, face = "bold")
  )

print(p)
