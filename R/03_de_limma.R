# 03_de_limma.R
# Differential expression of TCGA-BRCA primary tumor vs. solid tissue
# normal, protein-coding genes only, using limma-voom.
#
# limma-voom chosen over DESeq2 because DESeq2's dispersion estimation
# exceeds available memory on an 8 GB machine at this cohort size
# (n = 1,208). Agreement between the two methods on well-powered bulk
# RNA-seq is well documented (Law et al. 2014; Costa-Silva et al. 2017).
#
# Outputs:
#   data/processed/de_results_tumor_vs_normal.rds
#   data/processed/vst_matrix.rds   (log2-CPM from voom; used by 04_cluster.R)

mem.maxVSize(24000)

library(SummarizedExperiment)
library(edgeR)
library(limma)
library(dplyr)


# --- 1. Load cleaned data ------------------------------------------------

brca <- readRDS("data/processed/brca_clean_SE.rds")
cat("Loaded", ncol(brca), "samples x", nrow(brca), "genes\n")
print(table(brca$sample_type))


# --- 2. Subset to protein-coding genes -----------------------------------

# Act 1 reproduces the canonical PAM50 subtyping, defined on protein-coding
# genes. lncRNA analysis is handled separately in a later script.
is_pc <- rowData(brca)$gene_type == "protein_coding"
brca  <- brca[is_pc, ]
cat("Protein-coding genes:", nrow(brca), "\n")


# --- 3. Build DGEList and filter low-count genes -------------------------

group <- factor(
  brca$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)

counts_mat <- assay(brca)
colnames(counts_mat) <- colnames(brca)    # preserve TCGA barcodes
dge <- DGEList(counts = counts_mat, group = group)
keep <- filterByExpr(dge, group = group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes after filter:", nrow(dge), "\n")


# --- 4. Normalize, voom, fit ---------------------------------------------

dge    <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ group)

v   <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)


# --- 5. Extract results --------------------------------------------------

res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# Attach gene symbols and biotypes; rename columns to match DESeq2
# conventions so downstream scripts remain method-agnostic.
row_meta <- rowData(brca)[rownames(res), ]
res_df <- res |>
  tibble::rownames_to_column("gene_id") |>
  mutate(
    gene_name = row_meta$gene_name,
    gene_type = row_meta$gene_type
  ) |>
  rename(
    log2FoldChange = logFC,
    pvalue         = P.Value,
    padj           = adj.P.Val
  )


# --- 6. Summary ----------------------------------------------------------

cat("\n=== DE summary ===\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("padj < 0.05:",
    sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("padj < 0.05 & |log2FC| > 1:",
    sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
        na.rm = TRUE), "\n\n")

cat("Top 10 by padj:\n")
print(head(res_df[, c("gene_name", "gene_type",
                      "log2FoldChange", "padj")], 10))


# --- 7. Save -------------------------------------------------------------

saveRDS(v$E,    "data/processed/vst_matrix.rds")
saveRDS(res_df, "data/processed/de_results_tumor_vs_normal.rds")