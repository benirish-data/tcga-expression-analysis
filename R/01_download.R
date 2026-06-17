# 01_download.R
# Download TCGA-BRCA RNA-seq expression data and clinical annotations
# from the GDC using TCGAbiolinks. Writes a raw SummarizedExperiment
# to data/raw/brca_raw_SE.rds for use by downstream scripts.

# Raise R's memory ceiling; the full SE exceeds the default on macOS.
mem.maxVSize(32000)

library(TCGAbiolinks)
library(SummarizedExperiment)


# --- 1. Build the query --------------------------------------------------

query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)


# --- 2. Download the files -----------------------------------------------

# Smaller chunks (50 files vs. default ~250) are more reliable on
# home connections and avoid a known retry bug in TCGAbiolinks.
GDCdownload(query, directory = "data/raw/GDCdata", files.per.chunk = 50)


# --- 3. Assemble into a SummarizedExperiment -----------------------------

brca_data <- GDCprepare(
  query     = query,
  directory = "data/raw/GDCdata"
)


# --- 4. Save for downstream scripts --------------------------------------

saveRDS(brca_data, file = "data/raw/brca_raw_SE.rds")


# --- 5. Check -----------------------------------------------------

cat("Samples:", ncol(brca_data), "\n")
cat("Genes:",   nrow(brca_data), "\n")
cat("Sample types:\n")
print(table(brca_data$sample_type))