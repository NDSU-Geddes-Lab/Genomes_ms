#!/usr/bin/env Rscript

library(tidyverse)

#### ================== PARAMETERS ================== ####

# Odds ratio threshold for candidate genes
odds_threshold <- 1   # <-- as requested

# Base directory with per-strain Prokka output (each strain has its own folder + .gff)
prokka_base <- "/mmfs1/projects/barney.geddes/tania.gupta/TG_comparitive_genomics/whole_genome/A_renamed_whole_genome_prokka_roary_results/prokka_output_whole_renamned"

# Directory with Scoary + Roary outputs (also where we’ll write the final/intermediate CSVs) create this directory- gene_mapping and add result.tsv and gene_presence_absence.csv in it
out_dir <- "/mmfs1/projects/barney.geddes/tania.gupta/TG_comparitive_genomics/whole_genome/A_renamed_whole_genome_prokka_roary_results/roary_output_whole_renamned/85/scoary2_changed_script/scoary2_log2_transform_comp/gene_mapping"

scoary_file <- file.path(out_dir, "result.tsv")
roary_file  <- file.path(out_dir, "gene_presence_absence.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### ================== READ SCOARY + ROARY ================== ####


message("Reading Scoary results from: ", scoary_file)
scoary_results <- read_table(scoary_file) %>%
  mutate(
    odds_ratio = as.double(odds_ratio),
    effect     = if_else(odds_ratio > 1, "positive", "negative")
  )

message("Reading Roary gene_presence_absence from: ", roary_file)
roary_results <- read_csv(roary_file)

# Standard Roary metadata columns (rest are candidates to be strain columns)
meta_cols <- c(
  "Gene", "Non-unique Gene name", "Annotation",
  "No. isolates", "No. sequences", "Avg sequences per isolate",
  "Genome Fragment", "Order within Fragment",
  "Accessory Fragment", "Accessory Order with Fragment",
  "QC"
)

# Keep only those meta_cols that actually exist in this file
meta_cols <- intersect(meta_cols, colnames(roary_results))

# First, get all non-metadata columns
candidate_cols <- setdiff(colnames(roary_results), meta_cols)

# Among those, keep only character/factor columns (usually the strain locus_tag columns)
strain_names <- candidate_cols[sapply(
  roary_results[candidate_cols],
  function(x) is.character(x) || is.factor(x)
)]

message("Strains detected from Roary (character columns only): ",
        paste(strain_names, collapse = ", "))

if (length(strain_names) == 0) {
  stop("No strain columns detected. Check that gene_presence_absence.csv has character columns for strains.")
}

#### ================== SELECT CANDIDATE GENES (odds_ratio < 1) ================== ####

message("Filtering Scoary hits for odds_ratio < ", odds_threshold)
cand_hits <- scoary_results %>%
  filter(odds_ratio < odds_threshold) %>%   # <-- CHANGED (> to <)
  left_join(roary_results, by = "Gene")

if (nrow(cand_hits) == 0) {
  stop("No candidate genes with odds_ratio < ", odds_threshold, " found in Scoary results.")
}

# Keep only Gene, odds_ratio, effect, and strain locus_tag columns
cand_hits <- cand_hits %>%
  select(Gene, odds_ratio, effect, all_of(strain_names))

# ==== INTERMEDIATE FILE 1: candidates with per-strain locus tags ====
cand_hits_file <- file.path(
  out_dir,
  paste0("candidate_genes_gene_presence_absence_OR_lt", odds_threshold, ".csv")  # <-- CHANGED (gt to lt)
)
write_csv(cand_hits, cand_hits_file)
message("Wrote candidate genes × strains (locus_tag) table to: ", cand_hits_file)

#### ================== LONG FORMAT: Gene × strain × locus_tag ================== ####

cand_long <- cand_hits %>%
  pivot_longer(
    cols      = all_of(strain_names),
    names_to  = "strain",
    values_to = "locus_tag"
  ) %>%
  filter(!is.na(locus_tag)) %>%
  # IMPORTANT: split cells that contain multiple locus tags like "TG10_05199;TG10_06646"
  separate_rows(locus_tag, sep = ";")

message("Total candidate Gene × strain × locus_tag entries after splitting: ",
        nrow(cand_long))

# ==== INTERMEDIATE FILE 2: long-format Gene × strain × locus_tag before GFF mapping ====
cand_long_file <- file.path(
  out_dir,
  paste0("candidate_genes_strain_locus_long_OR_lt", odds_threshold, ".csv")  # <-- CHANGED (gt to lt)
)
write_csv(cand_long, cand_long_file)
message("Wrote long-format candidate Gene × strain × locus_tag table to: ",
        cand_long_file)

#### ================== BUILD locus_tag → replicon MAP FOR ALL STRAINS ================== ####

message("Building locus_tag → seqid/replicon map from GFFs...")

gff_map_all <- map_dfr(strain_names, function(strain) {

  gff_file <- file.path(prokka_base, strain, paste0(strain, ".gff"))

  if (!file.exists(gff_file)) {
    warning("GFF file not found for strain ", strain, ": ", gff_file)
    return(tibble(strain = character(), locus_tag = character(),
                  seqid = character(), replicon = character()))
  }

  gff <- read_delim(
    gff_file, "\t", comment = "#",
    col_names = c("seqid", "source", "type", "start", "end",
                  "score", "strand", "phase", "attributes"),
    show_col_types = FALSE
  ) %>%
    filter(!is.na(source), type == "gene")

  # Extract locus_tag from attributes (Prokka style: locus_tag=XXXX)
  gff_lut <- gff %>%
    mutate(
      locus_tag = str_extract(attributes, "locus_tag=[^;]+"),
      locus_tag = str_remove(locus_tag, "locus_tag="),
      # Classify replicon based on seqid; tweak patterns if needed
      replicon = case_when(
        str_detect(seqid, regex("chrom", ignore_case = TRUE)) ~ "chromosome",
        str_detect(seqid, regex("psymA", ignore_case = TRUE)) ~ "pSymA",
        str_detect(seqid, regex("psymB", ignore_case = TRUE)) ~ "pSymB",

        # accessory replicons
        str_detect(seqid, regex("accessory1", ignore_case = TRUE)) ~ "accessory1",
        str_detect(seqid, regex("accessory2", ignore_case = TRUE)) ~ "accessory2",
        str_detect(seqid, regex("accessory3", ignore_case = TRUE)) ~ "accessory3",
        str_detect(seqid, regex("accessory4", ignore_case = TRUE)) ~ "accessory4",
        str_detect(seqid, regex("accessory",  ignore_case = TRUE)) ~ "accessory",

        TRUE ~ "other"  # if nothing matches any of the above
      )
    ) %>%
    filter(!is.na(locus_tag)) %>%
    select(locus_tag, seqid, replicon) %>%
    mutate(strain = strain)

  message("  Processed GFF for strain: ", strain, " (",
          nrow(gff_lut), " genes with locus_tag)")

  return(gff_lut)
})

if (nrow(gff_map_all) == 0) {
  stop("No GFF mapping information could be built. Check prokka_base path and GFF naming.")
}

#### ================== JOIN: ADD REPLICON CLASS ================== ####

cand_long_rep <- cand_long %>%
  left_join(gff_map_all, by = c("strain", "locus_tag"))

# Check unmapped entries
unmapped <- cand_long_rep %>% filter(is.na(replicon))
if (nrow(unmapped) > 0) {
  warning("There are ", nrow(unmapped),
          " Gene × strain × locus_tag entries where locus_tag was not found in GFFs.")

  # ==== INTERMEDIATE FILE 3: list of unmapped locus_tags ====
  unmapped_file <- file.path(out_dir, "unmapped_gene_locus_pairs.csv")
  write_csv(
    unmapped %>% select(Gene, strain, locus_tag) %>% arrange(strain),
    unmapped_file
  )
  message("Wrote unmapped Gene × strain × locus_tag entries to: ",
          unmapped_file)
}

#### ================== WIDE MATRIX: Gene × strain → replicon ================== ####

# If a Gene has multiple copies in a strain on different replicons,
# collapse them as "chromosome;pSymA" etc.
cand_matrix <- cand_long_rep %>%
  group_by(Gene, strain) %>%
  summarise(
    replicon = paste(sort(unique(replicon)), collapse = ";"),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = strain,
    values_from = replicon
  )

# Attach odds_ratio & effect from Scoary for each Gene
cand_matrix <- cand_hits %>%
  distinct(Gene, odds_ratio, effect) %>%
  left_join(cand_matrix, by = "Gene")

#### ================== WRITE FINAL OUTPUT ================== ####

out_csv <- file.path(
  out_dir,
  paste0("candidate_genes_replicon_matrix_OR_lt", odds_threshold, ".csv")  # <-- CHANGED (gt to lt)
)
write_csv(cand_matrix, out_csv)

message("==============================================")
message("Finished! Wrote candidate gene × strain → replicon sheet to:")
message("  ", out_csv)
message("Rows (genes): ", nrow(cand_matrix))
message("Columns (including metadata): ", ncol(cand_matrix))
message("==============================================")
