library(ggplot2)

df <- read.csv("roary_total_gene_clusters_by_blastp.csv")

p <- ggplot(df, aes(x = blastp_identity, y = total_gene_clusters)) +
  geom_line(linewidth = 1.2, color = "steelblue4") +
  geom_point(size = 3, color = "dodgerblue3") +
  geom_vline(xintercept = 85, linetype = "dashed", color = "red3", linewidth = 1) +
  theme_classic(base_size = 14) +
  labs(
    x = "Roary BLAST-P identity threshold (%)",
    y = "Total gene clusters",
    title = "Roary sensitivity analysis"
  )

ggsave("roary_blastp_gene_cluster_sensitivity_colored.pdf", p, width = 6, height = 4)
ggsave("roary_blastp_gene_cluster_sensitivity_colored.png", p, width = 6, height = 4, dpi = 300)
