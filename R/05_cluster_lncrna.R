# 05_cluster_lncrna.R
# Repeat the clustering analysis from 04_cluster.R using lncRNAs instead
# of protein-coding genes, and compare the resulting clustering both to
# the known PAM50 subtypes and to the protein-coding clustering.
#
# Outputs:
#   data/processed/pca_tumor_lncrna.rds
#   data/processed/hclust_assignments_lncrna.rds
#   data/processed/top_variable_lncrnas.rds
#   figures/pca_tumors_by_pam50_lncrna.png
#   figures/hclust_confusion_heatmap_lncrna.png

library(SummarizedExperiment)
library(edgeR)
library(limma)
library(dplyr)
library(ggplot2)
library(mclust)


# --- 1. Load data --------------------------------------------------------

brca             <- readRDS("data/processed/brca_clean_SE.rds")
pc_assignments   <- readRDS("data/processed/hclust_assignments.rds")


# --- 2. Subset to lncRNAs and the four main PAM50 tumor subtypes --------

is_lnc <- rowData(brca)$gene_type == "lncRNA"
cat("lncRNAs in dataset:", sum(is_lnc), "\n")

keep_samples <- brca$sample_type == "Primary Tumor" &
  brca$pam50 %in% c("LumA", "LumB", "Her2", "Basal")

brca_lnc <- brca[is_lnc, keep_samples]
labels   <- factor(brca_lnc$pam50, levels = c("LumA", "LumB", "Her2", "Basal"))

cat("Tumors retained:", ncol(brca_lnc), "\n")


# --- 3. Filter + voom log-CPM -------------------------------------------

# lncRNAs are generally lower-expressed than protein-coding genes, so
# filterByExpr's adaptive threshold handles them appropriately.
group  <- factor(rep("tumor", ncol(brca_lnc)))
counts <- assay(brca_lnc)
colnames(counts) <- colnames(brca_lnc)

dge  <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge, group = group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")

cat("lncRNAs after filter:", nrow(dge), "\n")

v    <- voom(dge)
expr <- v$E


# --- 4. Select and scale top 2,000 most variable lncRNAs ----------------

# If fewer than 2,000 lncRNAs survive filtering, take all of them.
n_top     <- min(2000, nrow(expr))
gene_vars <- apply(expr, 1, var)
top_ids   <- names(sort(gene_vars, decreasing = TRUE))[1:n_top]
expr_top  <- expr[top_ids, ]

expr_scaled <- t(scale(t(expr_top)))

cat("Genes used for clustering:", nrow(expr_scaled), "\n")


# --- 5. PCA --------------------------------------------------------------

pca <- prcomp(t(expr_scaled), center = FALSE, scale. = FALSE)

pc_df <- data.frame(
  PC1   = pca$x[, 1],
  PC2   = pca$x[, 2],
  pam50 = labels
)

var_exp <- summary(pca)$importance[2, 1:2] * 100


# --- 6. Hierarchical clustering ------------------------------------------

d    <- dist(t(expr_scaled), method = "euclidean")
hc   <- hclust(d, method = "ward.D2")
clus <- cutree(hc, k = 4)


# --- 7. Compare to PAM50 and to the protein-coding clustering -----------

# Align the protein-coding cluster assignments to the same samples.
pc_aligned <- pc_assignments[colnames(expr_scaled)]

ari_vs_pam50   <- adjustedRandIndex(clus, labels)
ari_vs_pc      <- adjustedRandIndex(clus, pc_aligned)

cat("\n=== lncRNA clustering results ===\n")
cat("ARI vs. PAM50:                    ", round(ari_vs_pam50, 3), "\n")
cat("ARI vs. protein-coding clustering:", round(ari_vs_pc,    3), "\n\n")

cat("Confusion matrix (lncRNA cluster vs. PAM50):\n")
print(table(cluster = clus, pam50 = labels))

cat("\nCross-clustering (lncRNA cluster vs. protein-coding cluster):\n")
print(table(lncrna = clus, protein_coding = pc_aligned))


# --- 8. Plots ------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

p_pca <- ggplot(pc_df, aes(PC1, PC2, color = pam50)) +
  geom_point(alpha = 0.7, size = 1.5) +
  labs(
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color    = "PAM50",
    title    = "TCGA-BRCA tumors by PAM50 subtype (lncRNAs)",
    subtitle = sprintf("PCA on top %d variable lncRNAs", n_top)
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/pca_tumors_by_pam50_lncrna.png", p_pca,
       width = 7, height = 5, dpi = 300)

conf <- as.data.frame(table(cluster = clus, pam50 = labels))

p_conf <- ggplot(conf, aes(factor(cluster), pam50, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white") +
  scale_fill_viridis_c() +
  labs(
    x        = "Unsupervised cluster (k = 4)",
    y        = "PAM50 subtype",
    fill     = "Samples",
    title    = "lncRNA clustering vs. PAM50 labels",
    subtitle = sprintf("Adjusted Rand Index = %.3f", ari_vs_pam50)
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/hclust_confusion_heatmap_lncrna.png", p_conf,
       width = 6, height = 5, dpi = 300)


# --- 9. Save objects -----------------------------------------------------

saveRDS(pca,     "data/processed/pca_tumor_lncrna.rds")
saveRDS(clus,    "data/processed/hclust_assignments_lncrna.rds")
saveRDS(top_ids, "data/processed/top_variable_lncrnas.rds")