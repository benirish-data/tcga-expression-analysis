# 04_cluster.R
# Cluster TCGA-BRCA primary tumors on the top 2,000 most variable
# protein-coding genes and compare the unsupervised clustering to
# PAM50 subtype labels. PCA for visualization, hierarchical clustering
# (Euclidean, Ward.D2) for quantitative comparison via adjusted Rand
# index.
#
# Outputs:
#   data/processed/pca_tumor_pc.rds
#   data/processed/hclust_assignments.rds
#   data/processed/top_variable_genes.rds
#   figures/pca_tumors_by_pam50.png
#   figures/hclust_confusion_heatmap.png

library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(mclust)


# --- 1. Load data --------------------------------------------------------

brca    <- readRDS("data/processed/brca_clean_SE.rds")
vst_mat <- readRDS("data/processed/vst_matrix.rds")


# --- 2. Subset to tumors in the four main PAM50 subtypes -----------------

keep_samples <- brca$sample_type == "Primary Tumor" &
  brca$pam50 %in% c("LumA", "LumB", "Her2", "Basal")

brca_t <- brca[, keep_samples]
expr   <- vst_mat[, colnames(brca_t)]
labels <- factor(brca_t$pam50, levels = c("LumA", "LumB", "Her2", "Basal"))

cat("Tumors retained:", ncol(expr), "\n")
print(table(labels))


# --- 3. Select and scale top 2,000 most variable genes -------------------

gene_vars <- apply(expr, 1, var)
top_ids   <- names(sort(gene_vars, decreasing = TRUE))[1:2000]
expr_top  <- expr[top_ids, ]

# Z-score each gene across samples so no single highly-expressed gene
# dominates the distance metric.
expr_scaled <- t(scale(t(expr_top)))


# --- 4. PCA --------------------------------------------------------------

# prcomp expects samples-as-rows.
pca <- prcomp(t(expr_scaled), center = FALSE, scale. = FALSE)

pc_df <- data.frame(
  PC1   = pca$x[, 1],
  PC2   = pca$x[, 2],
  PC3   = pca$x[, 3],
  pam50 = labels
)

var_exp <- summary(pca)$importance[2, 1:3] * 100


# --- 5. Hierarchical clustering ------------------------------------------

d    <- dist(t(expr_scaled), method = "euclidean")
hc   <- hclust(d, method = "ward.D2")
clus <- cutree(hc, k = 4)

ari <- adjustedRandIndex(clus, labels)

cat("\nAdjusted Rand Index (hclust vs. PAM50):", round(ari, 3), "\n\n")
cat("Confusion matrix (rows = hclust cluster, cols = PAM50):\n")
print(table(cluster = clus, pam50 = labels))


# --- 6. Plots ------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

p_pca <- ggplot(pc_df, aes(PC1, PC2, color = pam50)) +
  geom_point(alpha = 0.7, size = 1.5) +
  labs(
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color    = "PAM50",
    title    = "TCGA-BRCA tumors by PAM50 subtype",
    subtitle = "PCA on top 2,000 variable protein-coding genes"
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/pca_tumors_by_pam50.png", p_pca,
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
    title    = "Hierarchical clustering vs. PAM50 labels",
    subtitle = sprintf("Adjusted Rand Index = %.3f", ari)
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/hclust_confusion_heatmap.png", p_conf,
       width = 6, height = 5, dpi = 300)


# --- 7. Save objects for downstream use ----------------------------------

saveRDS(pca,     "data/processed/pca_tumor_pc.rds")
saveRDS(clus,    "data/processed/hclust_assignments.rds")
saveRDS(top_ids, "data/processed/top_variable_genes.rds")
