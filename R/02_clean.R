# 02_clean.R
# Load the raw TCGA-BRCA SummarizedExperiment produced by 01_download.R,
# attach PAM50 subtype labels, deduplicate patients, and save a cleaned
# object to data/processed/brca_clean_SE.rds for downstream analysis.
#
# Low-expression gene filtering is deferred to the analysis scripts,
# since different methods (DESeq2, clustering, PCA) benefit from
# different filters.

mem.maxVSize(32000)

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)


# --- 1. Load raw data ----------------------------------------------------

brca_raw <- readRDS("data/raw/brca_raw_SE.rds")

n_samples_start <- ncol(brca_raw)
n_genes_start   <- nrow(brca_raw)


# --- 2. Pull PAM50 molecular subtype calls -------------------------------

# TCGAquery_subtype returns the curated subtype table from the 2012 TCGA
# BRCA paper; PAM50 calls are not in the default TCGAbiolinks output.
subtypes <- TCGAquery_subtype(tumor = "BRCA")

subtypes_small <- subtypes |>
  select(patient = patient, pam50 = BRCA_Subtype_PAM50)


# --- 3. Attach subtype labels to sample metadata -------------------------

sample_meta <- as.data.frame(colData(brca_raw)) |>
  mutate(patient = substr(barcode, 1, 12)) |>
  left_join(subtypes_small, by = "patient") |>
  # PAM50 is only defined on tumors — blank it out for normal samples
  # (otherwise normals inherit their patient's tumor subtype via the join).
  mutate(pam50 = if_else(sample_type == "Primary Tumor", pam50, NA_character_))

cat("Samples with PAM50 call:", sum(!is.na(sample_meta$pam50)), "\n")
print(table(sample_meta$pam50, sample_meta$sample_type, useNA = "ifany"))


# --- 4. Deduplicate patients ---------------------------------------------

# If a patient has multiple samples of the same type, keep the one
# with the highest library size. Breaks ties on data quality.
lib_sizes <- colSums(assay(brca_raw))

sample_meta <- sample_meta |>
  mutate(lib_size = lib_sizes[barcode]) |>
  group_by(patient, sample_type) |>
  slice_max(lib_size, n = 1, with_ties = FALSE) |>
  ungroup()

brca_clean <- brca_raw[, sample_meta$barcode]

# Reattach metadata; match() guarantees row order aligns with colnames.
colData(brca_clean) <- DataFrame(
  sample_meta[match(colnames(brca_clean), sample_meta$barcode), ]
)


# --- 5. Save -------------------------------------------------------------

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(brca_clean, file = "data/processed/brca_clean_SE.rds")


# --- 6. Summary ----------------------------------------------------------

cat("\n=== Cleaning summary ===\n")
cat("Samples before:", n_samples_start, "\n")
cat("Samples after: ", ncol(brca_clean), "\n")
cat("Genes:         ", nrow(brca_clean), "(filtering deferred)\n\n")

cat("Sample types after cleaning:\n")
print(table(brca_clean$sample_type))

cat("\nPAM50 subtype distribution (tumors only):\n")
print(table(brca_clean$pam50[brca_clean$sample_type == "Primary Tumor"],
            useNA = "ifany"))