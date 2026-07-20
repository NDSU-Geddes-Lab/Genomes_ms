# =========================
# Pagel's lambda + Blomberg's K for ACCESSORY_PLASMID + plot on tree + Mantel test
# =========================

# Install packages (run once)
# install.packages(c("ape", "phytools", "readxl"))

library(ape)
library(phytools)
library(readxl)

# ---- User inputs ----
treefile <- "/Users/taniagupta/Downloads/chr_core_gene_alignment.aln.treefile"
excel_file <- "/Users/taniagupta/Desktop/genomes_paper/stats_phylotree_trait/sdw_stats.xlsx"  # keep same file if column is in here
excel_sheet <- 1

# ---- Read data ----
tree <- read.tree(treefile)

dat <- read_excel(excel_file, sheet = excel_sheet)
dat <- as.data.frame(dat)

# Clean
dat$Strain <- trimws(as.character(dat$Strain))

# ---- CHANGE THIS LINE to your accessory plasmid column name exactly as in Excel ----
# Examples: dat$accessory_plasmid, dat$Accessory_plasmid, dat$accessory, dat$acc_plasmid
dat$accessory_plasmid <- as.numeric(dat$accessory_plasmid)

# ---- Match strains between tree and data ----
keep <- intersect(tree$tip.label, dat$Strain)

cat("Tips in tree:", length(tree$tip.label), "\n")
cat("Rows in data:", nrow(dat), "\n")
cat("Shared strains:", length(keep), "\n")

missing_in_tree <- setdiff(dat$Strain, tree$tip.label)
missing_in_data <- setdiff(tree$tip.label, dat$Strain)

if (length(missing_in_tree) > 0) {
  cat("\nStrains in Excel but NOT in tree:\n")
  print(missing_in_tree)
}
if (length(missing_in_data) > 0) {
  cat("\nTips in tree but NOT in Excel:\n")
  print(missing_in_data)
}

# Prune tree and subset data
tree2 <- drop.tip(tree, setdiff(tree$tip.label, keep))
dat2 <- dat[match(tree2$tip.label, dat$Strain), ]

# Make named trait vector in the same order as tree tips
acc <- dat2$accessory_plasmid
names(acc) <- dat2$Strain

# Safety checks
stopifnot(all(names(acc) == tree2$tip.label))
stopifnot(!any(is.na(acc)))

# =========================
# Pagel's lambda test
# =========================
lambda_res <- phylosig(tree2, acc, method = "lambda", test = TRUE)

cat("\n=== Pagel's lambda result ===\n")
print(lambda_res)

# =========================
# Blomberg's K test
# =========================
K_res <- phylosig(tree2, acc, method = "K", test = TRUE)

cat("\n=== Blomberg's K result ===\n")
print(K_res)

# =========================
# Save both results to a text file (nice for records)
# =========================
sink("Accessory_plasmid_phylo_signal_results.txt")
cat("Treefile:", treefile, "\n")
cat("Excel:", excel_file, "\n")
cat("N tips (original):", length(tree$tip.label), "\n")
cat("N tips (used):", length(tree2$tip.label), "\n\n")

cat("=== Pagel's lambda ===\n")
print(lambda_res)
cat("\n")

cat("=== Blomberg's K ===\n")
print(K_res)
cat("\n")
sink()

cat("\nSaved results: Accessory_plasmid_phylo_signal_results.txt\n")

# =========================
# Plot ACCESSORY_PLASMID on the tree
# =========================
pdf("Accessory_plasmid_trait_on_tree.pdf", width = 8, height = 10)
plotTree(tree2, ftype = "i", fsize = 0.6)

# contMap returns an object; plot=TRUE draws it
contMap(tree2, acc, plot = TRUE, ftype = "i", fsize = 0.6, legend = TRUE)
dev.off()

cat("\nSaved plot: Accessory_plasmid_trait_on_tree.pdf\n")

# =========================
# Mantel test: phylogenetic distance vs ACCESSORY_PLASMID distance
# =========================
# install.packages("vegan")  # run once if needed
library(vegan)
library(ape)

# 1) Pairwise phylogenetic distances between tips (patristic distances)
Dphy <- cophenetic.phylo(tree2)   # matrix

# 2) Pairwise ACCESSORY_PLASMID distances
Dacc <- dist(acc, method = "euclidean")  # returns 'dist' object

# 3) Mantel test
set.seed(123)
mantel_res <- mantel(as.dist(Dphy), Dacc,
                     method = "spearman",
                     permutations = 9999)

cat("\n=== Mantel test result (phylogenetic distance vs ACCESSORY_PLASMID distance) ===\n")
print(mantel_res)

# Optional: scatterplot of pairwise distances
pdf("Mantel_scatter_phyloDist_vs_AccessoryPlasmidDist.pdf", width = 6, height = 5)
plot(as.vector(as.dist(Dphy)),
     as.vector(Dacc),
     xlab = "Phylogenetic distance (patristic)",
     ylab = "Accessory plasmid distance (|Δ|)",
     pch = 16, cex = 0.6)
mtext(paste0("Mantel Spearman r = ", round(mantel_res$statistic, 3),
             ", p = ", signif(mantel_res$signif, 3)))
dev.off()

cat("\nSaved plot: Mantel_scatter_phyloDist_vs_AccessoryPlasmidDist.pdf\n")

