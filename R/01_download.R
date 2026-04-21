# installing bioconductor
install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "edgeR", "limma", 
                       "SummarizedExperiment", "ComplexHeatmap", 
                       "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("tidyverse", "survminer", "survival", "pheatmap", "factoextra"))

# snapshotting environment
install.packages("renv")
renv::init()
file.exists("renv.lock")
list.files("renv")

# building query
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# downloading packages
library(TCGAbiolinks)
library(SummarizedExperiment)
.libPaths()
renv::status()
.libPaths()
renv::hydrate()
renv::dependencies()
list.files("R")
