# =====================================
# UCEC Project
# 03 - Top 20 Most Variable Genes
# =====================================

library(tidyverse)
library(edgeR)

dir.create(
  "Results/Exploratory",
  recursive = TRUE,
  showWarnings = FALSE
)

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

gene_names <- tcga_counts[[1]]

counts_matrix <- tcga_counts[, -1] |>
  as.data.frame() |>
  as.matrix()

rownames(counts_matrix) <- gene_names

# Retain the highest-total-expression row for duplicated symbols.
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

dge_variable <- edgeR::DGEList(counts = counts_matrix)
dge_variable <- edgeR::calcNormFactors(dge_variable)

log_cpm_variable <- edgeR::cpm(
  dge_variable,
  log = TRUE,
  prior.count = 1
)

gene_variance <- apply(log_cpm_variable, 1, var)

top20_variable_genes <- tibble::tibble(
  Gene = names(
    sort(
      gene_variance,
      decreasing = TRUE
    )
  )[1:20],
  Variance = sort(
    gene_variance,
    decreasing = TRUE
  )[1:20]
)

print(top20_variable_genes)

readr::write_csv(
  top20_variable_genes,
  "Results/Exploratory/Top20_Most_Variable_NonMito_Genes.csv"
)

p_top20_variable <- ggplot(
  top20_variable_genes,
  aes(
    x = reorder(Gene, Variance),
    y = Variance
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Top 20 Most Variable Non-Mitochondrial Genes",
    x = "Gene",
    y = "Variance of log2 CPM"
  )

print(p_top20_variable)

ggsave(
  "Figures/Exploratory/Top20_Most_Variable_Genes.png",
  plot = p_top20_variable,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
