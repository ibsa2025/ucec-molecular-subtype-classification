# =====================================
# UCEC Project
# 05 - Differential Expression Analysis
# limma-voom
# =====================================

library(tidyverse)
library(edgeR)
library(limma)
library(tibble)

dir.create(
  "Results",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/DE",
  recursive = TRUE,
  showWarnings = FALSE
)

tcga_counts <- readr::read_csv(
  paste0(
    "TCGA Data/subtype_tumor_only/",
    "tcga_ucec_subtype_counts_geneSymbol.csv"
  ),
  show_col_types = FALSE
)

tcga_meta <- readr::read_csv(
  paste0(
    "TCGA Data/subtype_tumor_only/",
    "tcga_ucec_subtype_metadata.csv"
  ),
  show_col_types = FALSE
)

gene_names <- tcga_counts[[1]]

counts_matrix <- tcga_counts[, -1] |>
  as.data.frame() |>
  as.matrix()

rownames(counts_matrix) <- gene_names

# Align metadata and expression samples.
if (!setequal(colnames(counts_matrix), tcga_meta$barcode)) {
  stop("TCGA count-matrix columns and metadata barcodes do not match.")
}

tcga_meta <- tcga_meta[
  match(colnames(counts_matrix), tcga_meta$barcode),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(counts_matrix), tcga_meta$barcode)
)

# Resolve duplicated gene symbols using the highest total count.
if (anyDuplicated(rownames(counts_matrix)) > 0) {
  total_expression <- rowSums(counts_matrix, na.rm = TRUE)
  keep_order <- order(total_expression, decreasing = TRUE)
  counts_matrix <- counts_matrix[keep_order, , drop = FALSE]
  counts_matrix <- counts_matrix[
    !duplicated(rownames(counts_matrix)),
    ,
    drop = FALSE
  ]
}

# Remove mitochondrial genes.
counts_matrix <- counts_matrix[
  !grepl("^MT-", rownames(counts_matrix)),
  ,
  drop = FALSE
]

group <- factor(
  tcga_meta$subtype3,
  levels = c(
    "CN_high",
    "CN_low",
    "Hypermutated"
  )
)

dge <- edgeR::DGEList(
  counts = counts_matrix,
  group = group
)

keep_expression <- edgeR::filterByExpr(
  dge,
  group = group
)

dge <- dge[
  keep_expression,
  ,
  keep.lib.sizes = FALSE
]

dge <- edgeR::calcNormFactors(dge)

design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- levels(group)

v <- limma::voom(
  dge,
  design,
  plot = FALSE
)

fit <- limma::lmFit(
  v,
  design
)

contrast_matrix <- limma::makeContrasts(
  CN_high_vs_CN_low = CN_high - CN_low,
  CN_high_vs_Hypermutated = CN_high - Hypermutated,
  CN_low_vs_Hypermutated = CN_low - Hypermutated,
  levels = design
)

fit2 <- limma::contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- limma::eBayes(fit2)

extract_de_table <- function(
    fitted_model,
    coefficient_name
) {
  limma::topTable(
    fitted_model,
    coef = coefficient_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  ) |>
    tibble::rownames_to_column("Gene")
}

de_CNhigh_CNlow <- extract_de_table(
  fit2,
  "CN_high_vs_CN_low"
)

de_CNhigh_Hyper <- extract_de_table(
  fit2,
  "CN_high_vs_Hypermutated"
)

de_CNlow_Hyper <- extract_de_table(
  fit2,
  "CN_low_vs_Hypermutated"
)

sig_CNhigh_CNlow <- de_CNhigh_CNlow |>
  dplyr::filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )

sig_CNhigh_Hyper <- de_CNhigh_Hyper |>
  dplyr::filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )

sig_CNlow_Hyper <- de_CNlow_Hyper |>
  dplyr::filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )

readr::write_csv(
  de_CNhigh_CNlow,
  "Results/DE_CNhigh_vs_CNlow_limma_voom.csv"
)

readr::write_csv(
  de_CNhigh_Hyper,
  "Results/DE_CNhigh_vs_Hypermutated_limma_voom.csv"
)

readr::write_csv(
  de_CNlow_Hyper,
  "Results/DE_CNlow_vs_Hypermutated_limma_voom.csv"
)

readr::write_csv(
  sig_CNhigh_CNlow,
  "Results/Sig_DE_CNhigh_vs_CNlow_limma_voom.csv"
)

readr::write_csv(
  sig_CNhigh_Hyper,
  "Results/Sig_DE_CNhigh_vs_Hypermutated_limma_voom.csv"
)

readr::write_csv(
  sig_CNlow_Hyper,
  "Results/Sig_DE_CNlow_vs_Hypermutated_limma_voom.csv"
)

de_summary <- tibble::tibble(
  Comparison = c(
    "CN-high vs CN-low",
    "CN-high vs Hypermutated",
    "CN-low vs Hypermutated"
  ),
  Significant_Genes = c(
    nrow(sig_CNhigh_CNlow),
    nrow(sig_CNhigh_Hyper),
    nrow(sig_CNlow_Hyper)
  )
)

readr::write_csv(
  de_summary,
  "Results/DE_Significant_Gene_Counts.csv"
)

cat("\nGenes retained after filterByExpr:", nrow(dge), "\n")
print(de_summary)
