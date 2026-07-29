# =====================================
# UCEC Project
# 01 - Data Exploration
# =====================================

library(tidyverse)
library(readr)

counts_path <- paste0(
  "TCGA Data/subtype_tumor_only/",
  "tcga_ucec_subtype_counts_geneSymbol.csv"
)

metadata_path <- paste0(
  "TCGA Data/subtype_tumor_only/",
  "tcga_ucec_subtype_metadata.csv"
)

if (!file.exists(counts_path)) {
  stop("TCGA count file not found: ", counts_path)
}

if (!file.exists(metadata_path)) {
  stop("TCGA metadata file not found: ", metadata_path)
}

tcga_counts <- readr::read_csv(
  counts_path,
  show_col_types = FALSE
)

tcga_meta <- readr::read_csv(
  metadata_path,
  show_col_types = FALSE
)

if (!"barcode" %in% names(tcga_meta)) {
  stop("The TCGA metadata does not contain a 'barcode' column.")
}

if (!"subtype3" %in% names(tcga_meta)) {
  stop("The TCGA metadata does not contain a 'subtype3' column.")
}

expression_barcodes <- colnames(tcga_counts)[-1]

if (!setequal(expression_barcodes, tcga_meta$barcode)) {
  stop("TCGA expression sample identifiers do not match metadata barcodes.")
}

tcga_meta <- tcga_meta[
  match(expression_barcodes, tcga_meta$barcode),
  ,
  drop = FALSE
]

stopifnot(
  identical(expression_barcodes, tcga_meta$barcode),
  nrow(tcga_meta) == 232
)

cat("\nTCGA count table dimensions:\n")
print(dim(tcga_counts))

cat("\nTCGA metadata dimensions:\n")
print(dim(tcga_meta))

cat("\nThree-class subtype distribution:\n")
print(table(tcga_meta$subtype3, useNA = "ifany"))
