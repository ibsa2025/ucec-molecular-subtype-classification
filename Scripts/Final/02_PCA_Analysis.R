# =====================================
# UCEC Project
# 02 - Principal Component Analysis
# =====================================

library(tidyverse)
library(edgeR)
library(patchwork)

dir.create(
  "Figures/PCA",
  recursive = TRUE,
  showWarnings = FALSE
)

counts_path <- paste0(
  "TCGA Data/subtype_tumor_only/",
  "tcga_ucec_subtype_counts_geneSymbol.csv"
)

metadata_path <- paste0(
  "TCGA Data/subtype_tumor_only/",
  "tcga_ucec_subtype_metadata.csv"
)

tcga_counts <- readr::read_csv(
  counts_path,
  show_col_types = FALSE
)

tcga_meta <- readr::read_csv(
  metadata_path,
  show_col_types = FALSE
)

gene_names <- tcga_counts[[1]]

counts_matrix <- tcga_counts[, -1] |>
  as.data.frame() |>
  as.matrix()

rownames(counts_matrix) <- gene_names

# Align metadata with expression columns.
if (!setequal(colnames(counts_matrix), tcga_meta$barcode)) {
  stop("TCGA expression samples and metadata barcodes do not match.")
}

tcga_meta <- tcga_meta[
  match(colnames(counts_matrix), tcga_meta$barcode),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(counts_matrix), tcga_meta$barcode)
)

# Remove duplicated gene symbols by retaining the row with the
# highest total expression.
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

subtype <- factor(
  tcga_meta$subtype3,
  levels = c("CN_high", "CN_low", "Hypermutated")
)

# Filter and normalise counts before PCA.
dge_pca <- edgeR::DGEList(
  counts = counts_matrix,
  group = subtype
)

keep_expression <- edgeR::filterByExpr(
  dge_pca,
  group = subtype
)

dge_pca <- dge_pca[
  keep_expression,
  ,
  keep.lib.sizes = FALSE
]

dge_pca <- edgeR::calcNormFactors(dge_pca)

log_cpm_pca <- edgeR::cpm(
  dge_pca,
  log = TRUE,
  prior.count = 1
)

# Select the 1,000 most variable genes.
gene_variance <- apply(log_cpm_pca, 1, var)

top_gene_count <- min(1000, length(gene_variance))

top1000_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[seq_len(top_gene_count)]

pca_input <- t(
  log_cpm_pca[
    top1000_genes,
    ,
    drop = FALSE
  ]
)

pca_result <- stats::prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE
)

# Percentage of total variance explained by each principal component.
pca_var <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

pca_df <- tibble::tibble(
  Barcode = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4],
  Subtype = subtype
)

subtype_labels <- c(
  "CN_high" = "CN-high",
  "CN_low" = "CN-low",
  "Hypermutated" = "Hypermutated"
)

p_pc12 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = Subtype
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    title = "A",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
    colour = "Molecular subtype"
  ) +
  scale_colour_discrete(
    labels = subtype_labels
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_pc34 <- ggplot(
  pca_df,
  aes(
    x = PC3,
    y = PC4,
    colour = Subtype
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    title = "B",
    x = paste0("PC3 (", round(pca_var[3], 1), "%)"),
    y = paste0("PC4 (", round(pca_var[4], 1), "%)"),
    colour = "Molecular subtype"
  ) +
  scale_colour_discrete(
    labels = subtype_labels
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

pca_combined <- (
  p_pc12 +
    p_pc34 +
    patchwork::plot_layout(guides = "collect")
) &
  theme(legend.position = "bottom")

print(pca_combined)

ggsave(
  filename = "Figures/PCA/PCA_TCGA_Subtypes_Top1000_Combined.png",
  plot = pca_combined,
  width = 11,
  height = 5.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Figures/PCA/PCA_TCGA_Subtypes_Top1000_Combined.pdf",
  plot = pca_combined,
  width = 11,
  height = 5.8,
  bg = "white"
)

readr::write_csv(
  tibble::tibble(
    Principal_Component = paste0("PC", seq_along(pca_var)),
    Variance_Percent = pca_var
  ),
  "Figures/PCA/PCA_Variance_Explained.csv"
)

cat("\nVariance explained by PC1-PC4 (%):\n")
print(round(pca_var[1:4], 1))
