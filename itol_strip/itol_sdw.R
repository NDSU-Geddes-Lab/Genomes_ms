library(readxl)
setwd("/Users/taniagupta/Desktop/genomes_paper/itol_work/competition")   # change path to your folder

library(readxl)

excel_file <- "itol_strip_competition.xlsx"
sheet <- 1      

# Read everything as text (works in all versions of readxl)
df <- read_excel(excel_file, sheet = sheet, col_types = "text")

# Use the first column as ID
id_col_name <- colnames(df)[1]
cat("Using ID column:", id_col_name, "\n")

ids <- df[[id_col_name]]

# Generate one file per column (except ID)
for (col_name in colnames(df)) {
  if (col_name == id_col_name) next
  
  colors <- df[[col_name]]
  out_file <- paste0("itol_", col_name, "_colorstrip_comp.txt")
  cat("Writing", out_file, "\n")
  
  con <- file(out_file, open = "w", encoding = "UTF-8")
  
  # iTOL HEADER
  writeLines("DATASET_COLORSTRIP", con)
  writeLines("SEPARATOR\tTAB", con)
  writeLines(paste0("DATASET_LABEL\t", col_name), con)
  writeLines("COLOR\t#000000", con)
  writeLines("", con)
  
  # DATA SECTION
  writeLines("DATA", con)
  
  for (i in seq_along(ids)) {
    node <- as.character(ids[i])
    col  <- as.character(colors[i])
    
    if (is.na(col) || col == "") next
    
    writeLines(paste(node, col, sep = "\t"), con)
  }
  
  close(con)
}

cat("Done. All files written.\n")

