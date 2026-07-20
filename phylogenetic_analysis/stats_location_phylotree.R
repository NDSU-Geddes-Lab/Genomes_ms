# ============================================================
# Mantel test: Geographic distance (lat/lon) vs Phylogenetic distance
# Inputs:
#   Metadata (Excel): /Users/taniagupta/Desktop/genomes_paper/maps/Stats_location.xlsx
#   Tree (Newick):    /Users/taniagupta/Downloads/psymb_core_genes.treefile
# Output:
#   - Mantel test result printed
#   - Scatter plot (geo vs phylo) shown + saved as PNG
# ============================================================

# ---- Packages ----
#install.packages(c("ape","vegan","geosphere","readxl","dplyr","stringr"))
library(ape)
library(vegan)
library(geosphere)
library(readxl)
library(dplyr)
library(stringr)

# ---- Paths (yours) ----
meta_xlsx <- "/Users/taniagupta/Desktop/genomes_paper/maps/Stats_location.xlsx"
tree_file <- "/Users/taniagupta/Downloads/chr_core_gene_alignment.aln.treefile"

# ---- Read tree ----
tree <- read.tree(tree_file)

# ---- Read metadata (Excel) ----
# If your coordinates are not in the FIRST sheet, set sheet = "SheetName"
meta <- read_excel(meta_xlsx, sheet = 1)

# ---- Standardize column names ----
# EXPECTED columns (case-insensitive): Strain (or Sample ID), lat, lon
# If your file uses different names (e.g., Latitude/Longitude), this tries to catch them.
names(meta) <- str_replace_all(names(meta), "\\s+", "_")

# Try to detect columns robustly
strain_col <- names(meta)[tolower(names(meta)) %in% c("strain","sample_id","sampleid","sample","isolate","id")]
lat_col    <- names(meta)[tolower(names(meta)) %in% c("lat","latitude","y","gps_lat")]
lon_col    <- names(meta)[tolower(names(meta)) %in% c("lon","long","longitude","x","gps_lon")]

if (length(strain_col) != 1) stop("Could not uniquely detect the Strain column. Rename it to 'Strain' or 'Sample ID'.")
if (length(lat_col)    != 1) stop("Could not uniquely detect latitude column. Rename it to 'lat' or 'Latitude'.")
if (length(lon_col)    != 1) stop("Could not uniquely detect longitude column. Rename it to 'lon' or 'Longitude'.")

coords <- meta %>%
  transmute(
    Strain = as.character(.data[[strain_col]]),
    lat    = as.numeric(.data[[lat_col]]),
    lon    = as.numeric(.data[[lon_col]])
  ) %>%
  mutate(
    Strain = str_trim(Strain)
  ) %>%
  filter(!is.na(Strain) & Strain != "" & !is.na(lat) & !is.na(lon))

# ---- Match strains between tree tips and metadata ----
tree$tip.label <- str_trim(tree$tip.label)

common <- intersect(tree$tip.label, coords$Strain)

cat("Tree tips:", length(tree$tip.label), "\n")
cat("Metadata rows with coords:", nrow(coords), "\n")
cat("Common strains used:", length(common), "\n")

if (length(common) < 5) {
  # Helpful diagnostics if matching fails
  cat("\n⚠️ Too few strains matched.\n")
  cat("Example tree tip labels:\n")
  print(head(tree$tip.label, 20))
  cat("\nExample metadata Strain values:\n")
  print(head(coords$Strain, 20))
  stop("Fix strain name mismatches (case/spaces/prefixes) and re-run.")
}

# Drop tips not in metadata
tree2 <- drop.tip(tree, setdiff(tree$tip.label, common))

# Keep only matched rows and order to tree tip order
coords2 <- coords %>%
  filter(Strain %in% common) %>%
  distinct(Strain, .keep_all = TRUE)

coords2 <- coords2[match(tree2$tip.label, coords2$Strain), ]

stopifnot(all(tree2$tip.label == coords2$Strain))

# ---- Build phylogenetic distance matrix ----
phylo_dist <- cophenetic(tree2)  # strains x strains

# ---- Build geographic distance matrix (km) ----
# IMPORTANT: distm expects (lon, lat) order
geo_mat_m <- distm(coords2[, c("lon","lat")], fun = distHaversine)
geo_dist <- geo_mat_m / 1000  # km
rownames(geo_dist) <- coords2$Strain
colnames(geo_dist) <- coords2$Strain

# ---- Mantel test ----
set.seed(1)
mantel_res <- mantel(as.dist(phylo_dist), as.dist(geo_dist),
                     method = "pearson", permutations = 9999)

cat("\n================ Mantel test ================\n")
print(mantel_res)

# ---- Plot: geographic vs phylogenetic distances ----
x <- as.vector(as.dist(geo_dist))
y <- as.vector(as.dist(phylo_dist))

plot(x, y,
     xlab = "Geographic distance (km)",
     ylab = "Phylogenetic distance",
     pch = 16)

fit <- lm(y ~ x)
abline(fit)

mtext(sprintf("Mantel r = %.4f, p = %.4g", mantel_res$statistic, mantel_res$signif),
      side = 3, line = 0.2)

# ---- Save plot next to your metadata file directory ----
out_png <- file.path(dirname(meta_xlsx), "mantel_geo_vs_phylo_psymb.png")
png(out_png, width = 1800, height = 1400, res = 300)
plot(x, y,
     xlab = "Geographic distance (km)",
     ylab = "Phylogenetic distance",
     pch = 16)
abline(fit)
mtext(sprintf("Mantel r = %.4f, p = %.4g", mantel_res$statistic, mantel_res$signif),
      side = 3, line = 0.2)
dev.off()

cat("\nSaved plot to:\n", out_png, "\n")
