# ============================================================
# Install packages required for the UCEC analysis workflow
# Student ID: 3173204
#
# Run this file once before running the R Markdown workflow:
# source("install_packages.R")
# ============================================================

cran_packages <- c(
  "BiocManager",
  "caret",
  "flextable",
  "ggplot2",
  "ggrepel",
  "glmnet",
  "knitr",
  "officer",
  "patchwork",
  "pheatmap",
  "ranger",
  "readr",
  "rmarkdown",
  "stringr",
  "tibble",
  "tidyr",
  "tidyverse"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

bioconductor_packages <- c(
  "AnnotationDbi",
  "clusterProfiler",
  "edgeR",
  "enrichplot",
  "limma",
  "org.Hs.eg.db"
)

missing_bioc <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

cat("\nPackage installation check completed.\n")
