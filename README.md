# TCGA Breast Cancer Expression Analysis

I built this project to teach myself how to handle real bulk RNA sequencing data, reproducing the established finding that breast tumors cluster into four molecular subtypes (Luminal A/B, HER2+, Basal-like). I then extended the analysis by asking whether the same subtype structure emerges from long non-coding RNA (lncRNA) expression alone.

## Background

In 2000, Dr. Perou and his colleagues at Stanford discovered that breast cancer is not actually a single disease. It comprises at least four distinct subtypes with different treatment needs. A 2009 follow-up study by Parker et al. produced the PAM50 classifier: a 50-gene signature that reliably assigns a new tumor to one of these subtypes, and the source of the labels I use here. In this project, I reproduced that subtyping using public TCGA data. I then asked whether a largely non-protein-coding class of transcripts (lncRNAs) captures the same pattern or a different one.

## Data

- **Source:** TCGA-BRCA cohort via the NCI Genomic Data Commons (GDC)
- **Samples:** 1,224 downloaded; 1,208 retained after patient-level deduplication (1,095 primary tumor, 113 solid tissue normal)
- **Features:** 60,660 genes (Ensembl annotations, including protein-coding genes and lncRNAs)
- **Assay:** RNA-seq gene expression (STAR-aligned raw counts)
- **Subtype labels:** PAM50 calls from curated TCGA subtype annotations, merged via `TCGAquery_subtype()`
- **Download method:** `TCGAbiolinks` R package
- **Date accessed:** April 21, 2026

PAM50 distribution across tumors: LumA 562, LumB 209, Basal 190, HER2 82, Normal-like 40, unclassified 12. These proportions match published TCGA-BRCA breakdowns, which served as the first sanity check to ensure the data came through correctly.

## Methods

A few choices worth flagging:

- **DE method: limma-voom over DESeq2.** DESeq2's dispersion estimation on a cohort this size (n = 1,208) crashed R on the 8 GB machine I was working on. limma-voom is well-documented to agree closely with DESeq2 on well-powered bulk RNA-seq, so this was a practical substitution.

- **Patient-level deduplication.** TCGA contains patients with multiple samples of the same type. I kept one per patient × sample-type, breaking ties on library size. Without this, downstream statistical tests would treat non-independent samples as independent.

- **Gene filtering deferred to the analysis step.** The cleaning script keeps all genes. Filtering happens per-analysis because DESeq2, limma-voom, PCA, and clustering each benefit from different filters.

- **Protein-coding genes only for the main analysis.** Both the original subtype discovery (Perou 2000) and the PAM50 classifier (Parker 2009) are defined on protein-coding genes. The lncRNA analysis runs the same pipeline on the non-coding subset separately.

**Parameters and thresholds:**

- Low-count filter: `edgeR::filterByExpr()` with default settings (grouped by tumor/normal)
- DE significance: adjusted p-value < 0.05 (Benjamini-Hochberg) and |log2 fold change| > 1
- Clustering input: top 2,000 most variable protein-coding genes after log-CPM transformation, z-scored per gene
- Clustering method: hierarchical, Euclidean distance, Ward.D2 linkage, k = 4
- Comparison to PAM50: adjusted Rand index and confusion matrix

## Key Findings

Unsupervised clustering on variable-gene expression partially recovers PAM50 subtypes (ARI = 0.35).

Basal tumors separate cleanly — 185 of 190 land in a single cluster. Basal is known to be the most distinct subtype, so this is expected.

HER2 is messier. 52 of 82 form their own cluster, but about a third get pulled into the luminal group. That's consistent with HER2+ tumors often co-expressing hormone receptors.

The hard case is Luminal A vs. Luminal B. Unsupervised clustering doesn't separate them — they end up mixed together in one cluster, with a second cluster holding what looks like pure LumA. This isn't surprising. The LumA/LumB distinction depends on proliferation genes that PAM50's trained 50-gene signature picks up, but variance-based clustering on ~2,000 genes doesn't. To reliably tell them apart, you need a supervised method.

The lncRNA analysis was a bit more interesting. Despite using only the ~2,000 lncRNAs that survived expression filtering (versus a much larger pool of ~17,000 protein-coding genes), lncRNAs recover PAM50 subtypes essentially as well as protein-coding genes do (ARI 0.353 vs. 0.347). Basal isolates cleanly in both. HER2 partially splits in both. LumA/LumB mix in both.

The two clusterings moderately agree with each other (ARI 0.29): both recover Basal the same way, but they disagree on which specific luminal tumors belong with which. That suggests lncRNAs and protein-coding genes see partially different aspects of the biology, even when the clusters they produce look similar overall.

## Project Structure
```
tcga-expression-analysis/
├── R/
│   ├── 01_download.R          # Download raw TCGA-BRCA data from GDC
│   ├── 02_clean.R             # Attach subtype labels, deduplicate patients
│   ├── 03_de_limma.R          # Differential expression (tumor vs normal) via limma-voom
│   ├── 04_cluster.R           # PCA and hierarchical clustering by PAM50 subtype
│   └── 05_cluster_lncrna.R    # Repeat clustering on lncRNAs; compare to protein-coding
├── data/
│   ├── raw/                   # Raw downloaded data (gitignored)
│   └── processed/             # Cleaned and analyzed objects (gitignored)
├── figures/                   # Exported plots
├── analysis.qmd               # Main Quarto report (narrative + figures)
├── renv.lock                  # Pinned package versions
└── README.md
```

## Reproducibility

1. Clone this repository
2. Open `tcga-expression-analysis.Rproj` in RStudio
3. Run `renv::restore()` in the Console to install the exact package versions used
4. Source scripts in `R/` in numerical order

The download step (`01_download.R`) requires ~5 GB of disk space and 30 min – 2 hr depending on connection speed. All downstream scripts load cached data and complete within a few minutes each.

Memory note: assembling the full SummarizedExperiment (script 01) can exceed R's default 16 GB vector memory limit on macOS; the scripts raise it at the top.

## Tech Stack

**Core:** R 4.5, managed via [`renv`](https://rstudio.github.io/renv/)

**Bioconductor:** TCGAbiolinks, SummarizedExperiment, edgeR, limma

**CRAN:** dplyr, ggplot2, mclust

## A note on scope

This reproduces known findings. I'm not claiming a new scientific discovery. The point is to demonstrate the full bioinformatics workflow on real public data with an original twist (the lncRNA clustering comparison).