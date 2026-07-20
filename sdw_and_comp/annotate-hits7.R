library(tidyverse)

## NOTE: This script assumes you have your Roary/Scoary results and
##       .gff files in the current working directory.

## Files
library(tidyverse)

scoary_file <- "result.tsv"
roary_file <- "gene_presence_absence.csv"

# Change this if your KH46 column has a slightly different name
strain <- "KH46"

# -----------------------------
# Read Scoary results
# -----------------------------
scoary_results <- read_tsv(scoary_file, show_col_types = FALSE) %>%
  mutate(
    odds_ratio = as.double(odds_ratio),
    effect = if_else(odds_ratio > 1, "positive", "negative")
  )

# -----------------------------
# Read Roary gene presence/absence table
# -----------------------------
roary_results <- read_csv(roary_file, show_col_types = FALSE)

# Check KH46 column
if (!strain %in% colnames(roary_results)) {
  cat("Could not find column:", strain, "\n")
  cat("Possible KH46 columns are:\n")
  print(grep("KH46", colnames(roary_results), value = TRUE))
  stop("Please update the strain name.")
}

# Check Annotation column
if (!"Annotation" %in% colnames(roary_results)) {
  stop("The gene_presence_absence.csv file does not contain an Annotation column.")
}

# -----------------------------
# Join Scoary hits with Roary annotation
# -----------------------------
kh46_hits_annotated <- scoary_results %>%
  left_join(roary_results, by = "Gene") %>%
  mutate(
    KH46_gene_hit = .data[[strain]]
  ) %>%
  filter(
    !is.na(KH46_gene_hit),
    KH46_gene_hit != ""
  ) %>%
  separate_rows(KH46_gene_hit, sep = "\\t|;|,") %>%
  mutate(
    KH46_gene_hit = str_trim(KH46_gene_hit)
  ) %>%
  filter(KH46_gene_hit != "") %>%
  select(
    Gene,
    KH46_gene_hit,
    Annotation,
    effect,
    odds_ratio,
    everything()
  )

# -----------------------------
# Write output files
# -----------------------------
write_tsv(
  kh46_hits_annotated,
  "KH46_Scoary_hits_with_Roary_annotation.tsv"
)

write_tsv(
  kh46_hits_annotated %>% filter(effect == "positive"),
  "KH46_positive_hits_with_Roary_annotation.tsv"
)

write_tsv(
  kh46_hits_annotated %>% filter(effect == "negative"),
  "KH46_negative_hits_with_Roary_annotation.tsv"
)

# Print summary
cat("Total KH46 annotated Scoary hits:", nrow(kh46_hits_annotated), "\n")
cat("Positive hits:", sum(kh46_hits_annotated$effect == "positive"), "\n")
cat("Negative hits:", sum(kh46_hits_annotated$effect == "negative"), "\n")
