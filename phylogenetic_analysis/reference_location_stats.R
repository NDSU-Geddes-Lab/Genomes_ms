# ============================================================
# Mantel test:
# Existing region classification vs phylogenetic distance
#
# Metadata:
# /Users/taniagupta/Desktop/genomes_paper/reference_stats/location_of_strains.xlsx
#
# Tree:
# /Users/taniagupta/Desktop/genomes_paper/reference_stats/pSymA_core_gene_alignment.aln.treefile
# ============================================================

library(ape)
library(vegan)
library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)

# -----------------------------
# Input paths
# -----------------------------

meta_xlsx <- "/Users/taniagupta/Desktop/genomes_paper/reference_stats/location_of_strains.xlsx"
tree_file <- "/Users/taniagupta/Desktop/genomes_paper/reference_stats/chr_core_gene_alignment.aln.treefile"

outdir <- dirname(meta_xlsx)

analysis_name <- tools::file_path_sans_ext(basename(tree_file))

# -----------------------------
# Read tree and metadata
# -----------------------------

tree <- read.tree(tree_file)
meta <- read_excel(meta_xlsx, sheet = 1)

names(meta) <- str_replace_all(names(meta), "\\s+", "_")

if (!"Strain" %in% names(meta)) {
  stop("Metadata must contain a column named Strain")
}

if (!"Location" %in% names(meta)) {
  stop("Metadata must contain a column named Location")
}

# -----------------------------
# Clean tree tip names
# -----------------------------
# Examples:
# TG102_chromosome -> TG102
# TG102_pSymA -> TG102
# TG102_pSymB -> TG102
# TG102_accessory1 -> TG102

original_tip_labels <- tree$tip.label

tree$tip.label <- tree$tip.label %>%
  str_trim() %>%
  str_replace("_chromosome.*$", "") %>%
  str_replace("_pSymA.*$", "") %>%
  str_replace("_pSymB.*$", "") %>%
  str_replace("_accessory.*$", "")

tree_tip_table <- tibble(
  original_tip = original_tip_labels,
  cleaned_tip = tree$tip.label,
  match_key = str_to_lower(tree$tip.label)
)

write.csv(
  tree_tip_table,
  file.path(outdir, paste0(analysis_name, "_cleaned_tree_tip_names.csv")),
  row.names = FALSE
)

# -----------------------------
# Use existing Location column as Region
# -----------------------------
# Your Location column should contain:
# Mediterranean
# Temperate Europe
# Central Asia
# North America

meta2 <- meta %>%
  transmute(
    Strain = as.character(Strain),
    Region_raw = as.character(Location)
  ) %>%
  mutate(
    Strain = str_trim(Strain),
    Region_raw = str_squish(Region_raw),
    Region_clean = str_to_lower(Region_raw),
    match_key = str_to_lower(Strain),
    
    Region = case_when(
      Region_clean == "mediterranean" ~ "Mediterranean",
      Region_clean == "temperate europe" ~ "Temperate Europe",
      Region_clean == "central asia" ~ "Central Asia",
      Region_clean == "north america" ~ "North America",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Strain), Strain != "", !is.na(Region_raw), Region_raw != "")

if (any(is.na(meta2$Region))) {
  cat("\nThese values in Location were not recognized:\n")
  print(unique(meta2$Region_raw[is.na(meta2$Region)]))
  stop("Fix spelling in the Location column or add the value to case_when().")
}

cat("\n================ Region counts ================\n")
print(table(meta2$Region))

write.csv(
  meta2,
  file.path(outdir, paste0(analysis_name, "_strain_region_classification.csv")),
  row.names = FALSE
)

# -----------------------------
# Match tree tips and metadata
# -----------------------------

tree_keys <- str_to_lower(tree$tip.label)
meta_keys <- meta2$match_key

common_keys <- intersect(tree_keys, meta_keys)

tree_not_in_metadata <- setdiff(tree_keys, meta_keys)
metadata_not_in_tree <- setdiff(meta_keys, tree_keys)

cat("\n================ Matching diagnostics ================\n")

cat("\nTree tips after cleaning:", length(tree$tip.label), "\n")
cat("Metadata rows:", nrow(meta2), "\n")
cat("Common strains used:", length(common_keys), "\n")

cat("\nTree tips NOT found in metadata:", length(tree_not_in_metadata), "\n")
print(tree_tip_table %>% filter(match_key %in% tree_not_in_metadata))

cat("\nMetadata strains NOT found in tree:", length(metadata_not_in_tree), "\n")
print(meta2 %>% filter(match_key %in% metadata_not_in_tree) %>% select(Strain, Region_raw, Region))

write.csv(
  tree_tip_table %>% filter(match_key %in% tree_not_in_metadata),
  file.path(outdir, paste0(analysis_name, "_tree_tips_not_in_metadata.csv")),
  row.names = FALSE
)

write.csv(
  meta2 %>% filter(match_key %in% metadata_not_in_tree) %>% select(Strain, Region_raw, Region),
  file.path(outdir, paste0(analysis_name, "_metadata_strains_not_in_tree.csv")),
  row.names = FALSE
)

if (length(common_keys) < 5) {
  cat("\nToo few strains matched.\n")
  cat("\nExample cleaned tree tip labels:\n")
  print(head(tree$tip.label, 30))
  cat("\nExample metadata strain names:\n")
  print(head(meta2$Strain, 30))
  stop("Fix strain name mismatches and rerun.")
}

# -----------------------------
# Keep matched strains only
# -----------------------------

tips_to_keep <- tree$tip.label[tree_keys %in% common_keys]
tree2 <- keep.tip(tree, tips_to_keep)

meta_matched <- meta2 %>%
  filter(match_key %in% str_to_lower(tree2$tip.label)) %>%
  distinct(match_key, .keep_all = TRUE)

meta_matched <- meta_matched[match(str_to_lower(tree2$tip.label), meta_matched$match_key), ]

stopifnot(all(str_to_lower(tree2$tip.label) == meta_matched$match_key))

write.csv(
  meta_matched,
  file.path(outdir, paste0(analysis_name, "_matched_strains_used_for_region_mantel.csv")),
  row.names = FALSE
)

# -----------------------------
# Phylogenetic distance matrix
# -----------------------------

phylo_dist <- cophenetic(tree2)

# -----------------------------
# Region distance matrix
# -----------------------------
# Same region = 0
# Different region = 1

region_vector <- meta_matched$Region
names(region_vector) <- meta_matched$Strain

region_dist <- outer(
  region_vector,
  region_vector,
  FUN = function(x, y) ifelse(x == y, 0, 1)
)

rownames(region_dist) <- meta_matched$Strain
colnames(region_dist) <- meta_matched$Strain

# -----------------------------
# Mantel test
# -----------------------------

set.seed(1)

mantel_res <- mantel(
  as.dist(phylo_dist),
  as.dist(region_dist),
  method = "pearson",
  permutations = 9999
)

cat("\n================ Mantel test ================\n")
print(mantel_res)

sink(file.path(outdir, paste0(analysis_name, "_mantel_region_classification_result.txt")))
cat("================ Mantel test: region classification ================\n")
print(mantel_res)
sink()

# -----------------------------
# PERMANOVA / adonis2
# Useful because Region is categorical
# -----------------------------

set.seed(1)

adonis_res <- adonis2(
  as.dist(phylo_dist) ~ Region,
  data = meta_matched,
  permutations = 9999
)

cat("\n================ PERMANOVA / adonis2 ================\n")
print(adonis_res)

sink(file.path(outdir, paste0(analysis_name, "_adonis_region_classification_result.txt")))
cat("================ PERMANOVA/adonis2: phylogenetic distance by region ================\n")
print(adonis_res)
sink()

# -----------------------------
# Plot: phylogenetic distance within vs between regions
# -----------------------------

plot_df <- tibble(
  phylogenetic_distance = as.vector(as.dist(phylo_dist)),
  region_distance = as.vector(as.dist(region_dist))
) %>%
  mutate(
    comparison = ifelse(region_distance == 0, "Same region", "Different region")
  )

p <- ggplot(plot_df, aes(x = comparison, y = phylogenetic_distance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  theme_bw() +
  labs(
    x = "",
    y = "Phylogenetic distance",
    title = "Phylogenetic distance within vs between geographic regions",
    subtitle = paste0(
      "Mantel r = ",
      round(mantel_res$statistic, 4),
      ", p = ",
      signif(mantel_res$signif, 4)
    )
  )

print(p)

ggsave(
  filename = file.path(outdir, paste0(analysis_name, "_region_phylogenetic_distance_boxplot.png")),
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)

# -----------------------------
# Final output summary
# -----------------------------

cat("\n================ Saved outputs ================\n")
cat(file.path(outdir, paste0(analysis_name, "_cleaned_tree_tip_names.csv")), "\n")
cat(file.path(outdir, paste0(analysis_name, "_strain_region_classification.csv")), "\n")
cat(file.path(outdir, paste0(analysis_name, "_tree_tips_not_in_metadata.csv")), "\n")
cat(file.path(outdir, paste0(analysis_name, "_metadata_strains_not_in_tree.csv")), "\n")
cat(file.path(outdir, paste0(analysis_name, "_matched_strains_used_for_region_mantel.csv")), "\n")
cat(file.path(outdir, paste0(analysis_name, "_mantel_region_classification_result.txt")), "\n")
cat(file.path(outdir, paste0(analysis_name, "_adonis_region_classification_result.txt")), "\n")
cat(file.path(outdir, paste0(analysis_name, "_region_phylogenetic_distance_boxplot.png")), "\n")