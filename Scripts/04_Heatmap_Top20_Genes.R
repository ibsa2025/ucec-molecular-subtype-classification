# =====================================
# UCEC Project
# 04 - Heatmap of Top 20 Variable Genes
# =====================================

library(tidyverse)
library(edgeR)
library(pheatmap)

dir.create(
  "Figures/Exploratory",
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

tcga_meta <- tcga_meta[
  match(colnames(counts_matrix), tcga_meta$barcode),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(counts_matrix), tcga_meta$barcode)
)

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

counts_matrix <- counts_matrix[
  !grepl("^MT-", rownames(counts_matrix)),
  ,
  drop = FALSE
]

dge_heatmap <- edgeR::DGEList(counts = counts_matrix)
dge_heatmap <- edgeR::calcNormFactors(dge_heatmap)

log_cpm_heatmap <- edgeR::cpm(
  dge_heatmap,
  log = TRUE,
  prior.count = 1
)

gene_variance <- apply(log_cpm_heatmap, 1, var)

top20_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[1:20]

heatmap_matrix <- log_cpm_heatmap[
  top20_genes,
  ,
  drop = FALSE
]

# Standardise each gene across samples.
heatmap_scaled <- t(
  scale(
    t(heatmap_matrix)
  )
)

annotation_col <- data.frame(
  Subtype = factor(
    tcga_meta$subtype3,
    levels = c(
      "CN_high",
      "CN_low",
      "Hypermutated"
    )
  )
)

rownames(annotation_col) <- tcga_meta$barcode

heatmap_file <- "Figures/Exploratory/Heatmap_Top20_Variable_Genes.png"

heatmap_object <- pheatmap::pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Top 20 Variable Genes Across TCGA-UCEC Samples",
  silent = TRUE
)

png(
  filename = heatmap_file,
  width = 12,
  height = 7,
  units = "in",
  res = 300
)

grid::grid.newpage()
grid::grid.draw(heatmap_object$gtable)
dev.off()

file.info(heatmap_file)[, c("size", "mtime")]
