# ---------------------------
# Clade vs Rest: overlap + Venn + Excel outputs
# ---------------------------

# If needed (run once):
install.packages(c("ggVennDiagram"))

library(readxl)
library(openxlsx)
library(ggVennDiagram)
library(ggplot2)

# ---- paths ----
base_dir <- "/Users/taniagupta/Desktop/genomes_paper/accessory_different_clade"
gpa_file <- file.path(base_dir, "gene_presence_absence.csv")
clade_xlsx <- file.path(base_dir, "acc_different_clade.xlsx")

# outputs
out_xlsx <- file.path(base_dir, "clade_vs_rest_outputs.xlsx")
out_venn_png <- file.path(base_dir, "clade_vs_rest_venn.png")

# ---- helpers ----
clean_name <- function(x) {
  x <- gsub("\\.gff$|\\.gbk$|\\.fa$|\\.fasta$|\\.ffn$|\\.faa$", "", x, ignore.case = TRUE)
  trimws(x)
}

is_present_vec <- function(x) {
  !is.na(x) & trimws(as.character(x)) != ""
}

# ---- load Roary ----
gpa <- read.csv(gpa_file, stringsAsFactors = FALSE, check.names = FALSE)

# ---- load Excel strain list ----
clade_df <- read_excel(clade_xlsx)

# Assumption: strain names are in the FIRST column of the excel file
clade_strains <- clade_df[[1]]
clade_strains <- clean_name(as.character(clade_strains))
clade_strains <- unique(clade_strains[!is.na(clade_strains) & clade_strains != ""])

# ---- identify strain columns in Roary ----
meta_cols <- c(
  "Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
  "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
  "Accessory Fragment", "Accessory Order with Fragment", "QC",
  "Min group size nuc", "Max group size nuc", "Avg group size nuc"
)
meta_cols <- meta_cols[meta_cols %in% colnames(gpa)]

strain_cols_all <- setdiff(colnames(gpa), meta_cols)

# map cleaned strain colnames -> original
strain_cols_clean <- clean_name(strain_cols_all)
clean_to_original <- setNames(strain_cols_all, strain_cols_clean)

# match clade strains to Roary columns
clade_in_gpa_clean <- intersect(clade_strains, strain_cols_clean)
missing_clean <- setdiff(clade_strains, strain_cols_clean)

if (length(clade_in_gpa_clean) == 0) {
  stop("ERROR: None of the clade strains matched Roary column names.\n",
       "Fix: check naming in Excel vs gene_presence_absence.csv strain headers.")
}

clade_cols <- unname(clean_to_original[clade_in_gpa_clean])
rest_cols  <- setdiff(strain_cols_all, clade_cols)

# ---- build presence/absence boolean matrix ----
pres <- as.data.frame(lapply(gpa[strain_cols_all], is_present_vec), check.names = FALSE)

# ---- gene ID column ----
gene_id_col <- if ("Gene" %in% colnames(gpa)) "Gene" else colnames(gpa)[1]

# ---- define sets (ANY presence) ----
in_clade_any <- if (length(clade_cols) > 0) apply(pres[clade_cols], 1, any) else rep(FALSE, nrow(gpa))
in_rest_any  <- if (length(rest_cols)  > 0) apply(pres[rest_cols],  1, any) else rep(FALSE, nrow(gpa))

# Categories
clade_only_any <- in_clade_any & !in_rest_any
shared_any     <- in_clade_any &  in_rest_any
rest_only_any  <- !in_clade_any & in_rest_any

# Also (optional) strict: present in ALL clade strains
in_clade_all <- if (length(clade_cols) > 0) apply(pres[clade_cols], 1, all) else rep(FALSE, nrow(gpa))
clade_all_only_anyrest <- in_clade_all & !in_rest_any
clade_all_shared <- in_clade_all & in_rest_any

# ---- create data frames for export ----
meta_keep <- unique(c(meta_cols, gene_id_col))
meta_keep <- meta_keep[meta_keep %in% colnames(gpa)]

df_clade_only <- gpa[clade_only_any, meta_keep]
df_shared     <- gpa[shared_any,     meta_keep]
df_rest_only  <- gpa[rest_only_any,  meta_keep]

df_allclade_only <- gpa[clade_all_only_anyrest, meta_keep]
df_allclade_shared <- gpa[clade_all_shared, meta_keep]

counts <- data.frame(
  Metric = c(
    "Total genes in Roary table",
    "Clade strains requested (Excel)",
    "Clade strains matched in Roary",
    "Clade strains missing from Roary",
    "Rest strains (all others)",
    "Genes present in >=1 clade strain (ANY)",
    "Genes present in >=1 rest strain (ANY)",
    "Shared genes (ANY) [clade ∩ rest]",
    "Clade-only genes (ANY) [clade \\ rest]",
    "Rest-only genes (ANY) [rest \\ clade]",
    "Genes present in ALL clade strains",
    "ALL-clade shared with rest",
    "ALL-clade only (absent in rest)"
  ),
  Value = c(
    nrow(gpa),
    length(clade_strains),
    length(clade_cols),
    length(missing_clean),
    length(rest_cols),
    sum(in_clade_any),
    sum(in_rest_any),
    nrow(df_shared),
    nrow(df_clade_only),
    nrow(df_rest_only),
    sum(in_clade_all),
    nrow(df_allclade_shared),
    nrow(df_allclade_only)
  )
)

missing_df <- data.frame(Missing_Clade_Strain = missing_clean)

# ---- Venn diagram ----
# Use gene IDs for the sets (Venn needs unique items)
genes_clade_any <- unique(gpa[[gene_id_col]][in_clade_any])
genes_rest_any  <- unique(gpa[[gene_id_col]][in_rest_any])

venn_list <- list(
  Clade = genes_clade_any,
  Rest  = genes_rest_any
)

p <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "lightblue", high = "lightblue") +
  guides(fill = "none") +
  ggtitle("Clade vs Rest gene overlap (ANY presence)") +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(out_venn_png, p, width = 7, height = 6, dpi = 300)



# ---- Write Excel workbook ----
wb <- createWorkbook()

addWorksheet(wb, "Summary_counts")
writeData(wb, "Summary_counts", counts)

addWorksheet(wb, "Shared_ANY")
writeData(wb, "Shared_ANY", df_shared)

addWorksheet(wb, "Clade_only_ANY")
writeData(wb, "Clade_only_ANY", df_clade_only)

addWorksheet(wb, "Rest_only_ANY")
writeData(wb, "Rest_only_ANY", df_rest_only)

addWorksheet(wb, "Present_in_ALL_clade")
writeData(wb, "Present_in_ALL_clade", gpa[in_clade_all, meta_keep])

addWorksheet(wb, "ALLclade_shared_with_rest")
writeData(wb, "ALLclade_shared_with_rest", df_allclade_shared)

addWorksheet(wb, "ALLclade_only_absent_rest")
writeData(wb, "ALLclade_only_absent_rest", df_allclade_only)

addWorksheet(wb, "Clade_strains_used")
writeData(wb, "Clade_strains_used", data.frame(Clade_strain = clade_cols))

addWorksheet(wb, "Missing_clade_strains")
writeData(wb, "Missing_clade_strains", missing_df)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("DONE ✅")
message("Excel written: ", out_xlsx)
message("Venn PNG saved: ", out_venn_png)
if (length(missing_clean) > 0) message("Some strains did not match (see Excel sheet: Missing_clade_strains)")

