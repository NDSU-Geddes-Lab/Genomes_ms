# Combine 5 candidate-gene CSVs into ONE table (RStudio-friendly)
#To assign colors to respective replicon
#to create the gene presence absence matrix

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
})

# ---- Set your working directory ----
setwd("/Users/taniagupta/Desktop/genomes_paper/scoary2_figure")

# ---- Inputs ----
tree_file <- "core_gene_whole_85.treefile"

file_map <- c(
  "candidate_genes_replicon_comp_negative.csv" = "comp_neg",
  "candidate_genes_replicon_comp_positive.csv" = "comp_pos",
  "candidate_genes_replicon_positive_sdw.csv"  = "sdw_pos",
  "paper_gene_hits.csv"                        = "paper",
  "candidate_genes_sdw_negative.csv"           = "sdw_neg"
)

# ---- Output ----
out_file <- "combined_gene_hits_all_sources.csv"

# ---- Read strains from tree ----
tree <- ape::read.tree(tree_file)
strain_order <- tree$tip.label

# ---- Helper to read one file and force all tree strains as columns ----
read_and_align <- function(csv_path, source_label, strain_order) {
  df <- readr::read_csv(csv_path, show_col_types = FALSE)
  
  if (!("Gene" %in% colnames(df))) {
    stop("File missing required column 'Gene': ", csv_path)
  }
  
  # Drop common metadata cols if present (keep ONLY Gene + strain columns)
  meta_drop <- intersect(
    c("odds_ratio","effect","pvalue","p_value","sensitivity","specificity"),
    colnames(df)
  )
  df <- df %>% select(-all_of(meta_drop))
  
  # Ensure all tree strains exist as columns; create missing columns as NA
  missing_cols <- setdiff(strain_order, colnames(df))
  if (length(missing_cols) > 0) {
    df[missing_cols] <- NA
  }
  
  # Keep columns in consistent order: source, Gene, (tree strain order)
  df %>%
    mutate(source = source_label, .before = 1) %>%
    select(source, Gene, all_of(strain_order))
}

# ---- Build combined table ----
combined <- purrr::imap_dfr(file_map, function(label, file) {
  if (!file.exists(file)) stop("Missing file: ", file)
  read_and_align(file, label, strain_order)
})

# (Optional) sort rows by source then Gene for readability
combined <- combined %>% arrange(source, Gene)

# ---- Write combined output ----
readr::write_csv(combined, out_file)
message("Saved combined table to: ", file.path(getwd(), out_file))



#change the order you want it in and also add spaces between the gene hits if you want
# also add sdw and competiton binary data in the exel sheet over here




#To assign colors to respective replicon

library(readr)
library(openxlsx)
library(stringr)

setwd("/Users/taniagupta/Desktop/genomes_paper/scoary2_figure")

in_csv   <- "combined_gene_hits_all_sources.csv"
out_xlsx <- "combined_gene_hits_all_sources_COLORCODES.xlsx"

# ---- colors (your legend + TRUE/FALSE) ----
col_chrom   <- "#1F77B4"
col_psyma   <- "#FF7F0E"
col_psymb   <- "#2CA02C"
col_acc     <- "#9467BD"
col_contig  <- "#8C564B"
col_other   <- "#7F7F7F"
col_multi   <- "#E377C2"
col_absent  <- "#D9D9D9"

col_true    <- "#8B0000"  # red
col_false   <- "#FFF44F"  # yellow

df <- read_csv(in_csv, show_col_types = FALSE)
stopifnot(all(c("source","Gene") %in% colnames(df)))
strain_cols <- setdiff(colnames(df), c("source","Gene"))

# classify each cell -> hex code
to_hex <- function(x) {
  # TRUE/FALSE handling (case-insensitive)
  if (!is.na(x)) {
    x_chr <- as.character(x)
    xl0 <- tolower(trimws(x_chr))
    if (xl0 == "true")  return(col_true)
    if (xl0 == "false") return(col_false)
  }
  
  # NA/blank -> Absent
  if (is.na(x) || trimws(as.character(x)) == "" || as.character(x) %in% c("NA",".","NaN")) return(col_absent)
  
  x_chr <- as.character(x)
  xl <- tolower(gsub("\\s+", "", x_chr))
  
  if (str_detect(xl, "chrom"))     return(col_chrom)
  if (str_detect(xl, "psyma"))     return(col_psyma)
  if (str_detect(xl, "psymb"))     return(col_psymb)
  if (str_detect(xl, "accessory")) return(col_acc)     # accessory1-4 included
  if (str_detect(xl, "contig"))    return(col_contig)  # contig1-4 included
  if (str_detect(xl, "^other$"))   return(col_other)
  if (str_detect(x_chr, ";"))      return(col_multi)
  
  # anything unknown -> Other (change to col_absent if you prefer)
  col_other
}

hex_df <- df
hex_df[strain_cols] <- lapply(hex_df[strain_cols], function(v) vapply(v, to_hex, character(1)))

wb <- createWorkbook()
addWorksheet(wb, "Colorcodes")
writeData(wb, "Colorcodes", hex_df)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Wrote:", out_xlsx, "\n")





#to create the gene presence absence matrix


suppressPackageStartupMessages({
  library(ape)
  library(openxlsx)  # for .xlsx reading
  library(readr)     # for .csv reading
})

setwd("/Users/taniagupta/Desktop/genomes_paper/scoary2_figure")

# ---- INPUTS ----
tree_file <- "core_gene_whole_85.treefile"

# Your colorcode table (works if it's .xlsx OR .csv)
# If your file is named exactly "combined_gene_hits_all_sources_COLORCODES.xlsx", keep this:
heatmap_file <- "combined_gene_hits_all_sources_COLORCODES.xlsx"
# heatmap_file <- "combined_gene_hits_all_sources_COLORCODES.csv"  # if you saved as CSV instead

sheet_name <- "Colorcodes"   # the sheet name used when you wrote the xlsx (change if different)

# ---- OUTPUT ----
out_png <- "tree_plus_heatmap_colorcodes.png"

# ---- OPTIONS ----
draw_cell_grid <- TRUE     # TRUE = outline every cell with thin grid lines
grid_lwd <- 0.2
grid_col <- "black"

# ----------------------------
# Read tree
# ----------------------------
tree <- ape::read.tree(tree_file)
strain_order <- tree$tip.label

# ----------------------------
# Read heatmap table (hex codes)
# Expected columns: source, Gene, TG1, TG2, ...
# ----------------------------
if (grepl("\\.xlsx$", heatmap_file, ignore.case = TRUE)) {
  df <- openxlsx::read.xlsx(heatmap_file, sheet = sheet_name)
} else {
  df <- readr::read_csv(heatmap_file, show_col_types = FALSE)
}

stopifnot(all(c("source", "Gene") %in% colnames(df)))

# Strain columns = everything except source/Gene
strain_cols <- setdiff(colnames(df), c("source", "Gene"))

# Keep only strains that are in the tree (and in tree order)
common_strains <- intersect(strain_order, strain_cols)
if (length(common_strains) < 2) {
  stop("Too few matching strains between tree tips and heatmap columns.\n",
       "Example tree tips: ", paste(head(strain_order, 10), collapse = ", "), "\n",
       "Example heatmap cols: ", paste(head(strain_cols, 10), collapse = ", "))
}
df <- df[, c("source", "Gene", common_strains)]
common_strains <- strain_order[strain_order %in% common_strains]  # enforce tree order
df <- df[, c("source", "Gene", common_strains)]

# Replace blanks/NA with white (optional; your Absent is already #D9D9D9)
for (cn in common_strains) {
  v <- as.character(df[[cn]])
  v[is.na(v) | trimws(v) == "" | v %in% c("NA", ".", "NaN")] <- "#FFFFFF"
  df[[cn]] <- v
}

# Build matrix: strains (rows) x genes (cols)
gene_labels <- paste(df$source, df$Gene, sep = "|")  # keeps separators (blank rows) distinct-ish
hex_mat <- t(as.matrix(df[, common_strains, drop = FALSE]))  # now rows=strains, cols=genes
colnames(hex_mat) <- gene_labels
rownames(hex_mat) <- common_strains

# Map hex colors -> integers for image()
hex_vals <- unique(as.vector(hex_mat))
hex_vals <- hex_vals[!is.na(hex_vals)]
hex_vals <- unique(toupper(hex_vals))
# ensure all are like "#RRGGBB"
hex_vals <- hex_vals[grepl("^#[0-9A-F]{6}$", hex_vals)]
if (length(hex_vals) == 0) stop("No valid hex colors like #RRGGBB found in heatmap table.")

hex_to_int <- setNames(seq_along(hex_vals), hex_vals)

int_mat <- matrix(
  hex_to_int[toupper(hex_mat)],
  nrow = nrow(hex_mat),
  ncol = ncol(hex_mat),
  dimnames = dimnames(hex_mat)
)

# ----------------------------
# Plot (tree left, heatmap right)
# ----------------------------
pdf("tree_plus_heatmap_colorcodes.pdf", width = 13, height = 7.5)  # inches

#png(out_png, width = 2600, height = 1500, res = 200)

layout(matrix(c(1, 2), nrow = 1), widths = c(0.6, 2.9))

# ---- PANEL 1: TREE ----
par(mar = c(4, 0.2, 2, 1), xpd = NA)  # small right margin, allow drawing outside
plot(tree, show.tip.label = TRUE, cex = 0.25, label.offset = 0.0005)
par(xpd = FALSE)

# ---- PANEL 2: HEATMAP (variable-width columns) ----
par(mar = c(4, 1, 2, 1))

# Column widths: sdw + comp are 2x wider
base_w <- 1
wide_w <- 2

gene_source <- as.character(df$source)  # one per gene column
col_w <- ifelse(gene_source %in% c("sdw", "comp"), wide_w, base_w)

# keep your global scaling
box_width_factor <- 3.0
col_w <- col_w * box_width_factor

# x boundaries for each column
x_left  <- c(0, cumsum(col_w))[1:length(col_w)]
x_right <- x_left + col_w

# y boundaries for strains (rows)
nR <- nrow(hex_mat)     # strains
nC <- ncol(hex_mat)     # gene columns
y_bottom <- 0:(nR - 1)
y_top    <- y_bottom + 1

# set plotting region for heatmap
plot.new()
plot.window(
  xlim = c(min(x_left), max(x_right)),
  ylim = c(0, nR),
  xaxs = "i", yaxs = "i"
)

# draw heatmap cells
for (j in seq_len(nC)) {
  for (i in seq_len(nR)) {
    col_hex <- toupper(hex_mat[i, j])
    if (is.na(col_hex) || !grepl("^#[0-9A-F]{6}$", col_hex)) col_hex <- "#FFFFFF"
    
    rect(
      xleft   = x_left[j],
      xright  = x_right[j],
      ybottom = y_bottom[i],
      ytop    = y_top[i],
      col     = col_hex,
      border  = NA
    )
  }
}

# Optional: grid lines
if (draw_cell_grid) {
  # vertical boundaries (variable widths)
  abline(v = c(x_left, tail(x_right, 1)), col = grid_col, lwd = grid_lwd)
  # horizontal boundaries
  abline(h = 0:nR, col = grid_col, lwd = grid_lwd)
}

box()

dev.off()
cat("Saved:", out_png, "\n")













