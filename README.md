# TCGA Breast Cancer Expression Analysis

Differential gene expression analysis of TCGA-BRCA RNA-seq data, reproducing the canonical finding that breast tumors cluster into molecular subtypes (Luminal A/B, HER2+, Basal-like). The analysis is extended by asking whether the same subtype structure emerges from long non-coding RNA (lncRNA) expression alone.

## Background

Breast cancer is not a single disease. It comprises at least four molecularly distinct subtypes with different treatment responses and prognoses, originally defined by protein-coding gene expression patterns. This project reproduces that canonical subtyping using public TCGA data, then asks whether lncRNAs — a largely non-protein-coding class of transcripts — capture the same biological structure or reveal different axes of variation.

## Data

- **Source:** TCGA-BRCA cohort via the NCI Genomic Data Commons (GDC)
- **Samples:** 1,224 total — 1,111 primary tumor, 113 solid tissue normal
- **Features:** 60,660 genes (Ensembl annotations, includes protein-coding and lncRNAs)
- **Assay:** RNA-seq gene expression (STAR aligner, raw counts)
- **Access:** programmatic via the `TCGAbiolinks` R package
- **Date accessed:** April 21, 2026

## Methods

- Differential expression (tumor vs. normal) with `DESeq2`
- Variance-stabilizing transformation for downstream visualization
- Principal component analysis and hierarchical clustering for subtype discovery
- Separate parallel pipelines for protein-coding genes and lncRNAs
- Quantitative comparison of the two clusterings (adjusted Rand index)
- Pathway enrichment with `clusterProfiler`

## Key Findings

*(To be filled in after analysis is complete.)*

## Reproducibility

1. Clone this repository
2. Open `tcga-expression-analysis.Rproj` in RStudio
3. Run `renv::restore()` to install the exact package versions used
4. Source scripts in `R/` in numerical order (`01_download.R` → `02_clean.R` → ...)

The download step requires ~5 GB of disk space and 30 min – 2 hr depending on connection speed. Subsequent steps load cached data and run in minutes.

## Tech Stack

**R packages:** TCGAbiolinks, DESeq2, SummarizedExperiment, edgeR, limma, ComplexHeatmap, EnhancedVolcano, clusterProfiler, org.Hs.eg.db, tidyverse, pheatmap, factoextra, survival, survminer

**Environment:** managed via `renv` for reproducibility.

## Project Structure

```
tcga-expression-analysis/
├── R/                     # analysis scripts, numbered in execution order
├── data/
│   ├── raw/              # downloaded TCGA data (gitignored)
│   └── processed/        # intermediate results (gitignored)
├── figures/              # exported plots
├── analysis.qmd          # main Quarto report (narrative)
├── analysis.html         # rendered output
└── renv.lock             # pinned package versions
```