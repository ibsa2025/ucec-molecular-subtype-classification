# =====================================
# UCEC Project
# 07 - Top Differentially Expressed Genes
# =====================================

library(tidyverse)

de_CNhigh_CNlow <- readr::read_csv(
  "Results/DE_CNhigh_vs_CNlow_limma_voom.csv",
  show_col_types = FALSE
)

de_CNhigh_Hyper <- readr::read_csv(
  "Results/DE_CNhigh_vs_Hypermutated_limma_voom.csv",
  show_col_types = FALSE
)

de_CNlow_Hyper <- readr::read_csv(
  "Results/DE_CNlow_vs_Hypermutated_limma_voom.csv",
  show_col_types = FALSE
)

extract_top_genes <- function(de_table, n = 20) {
  de_table |>
    dplyr::arrange(
      adj.P.Val,
      dplyr::desc(abs(logFC))
    ) |>
    dplyr::slice_head(n = n)
}

top20_CNhigh_CNlow <- extract_top_genes(
  de_CNhigh_CNlow
)

top20_CNhigh_Hyper <- extract_top_genes(
  de_CNhigh_Hyper
)

top20_CNlow_Hyper <- extract_top_genes(
  de_CNlow_Hyper
)

readr::write_csv(
  top20_CNhigh_CNlow,
  "Results/Top20_DE_CNhigh_vs_CNlow.csv"
)

readr::write_csv(
  top20_CNhigh_Hyper,
  "Results/Top20_DE_CNhigh_vs_Hypermutated.csv"
)

readr::write_csv(
  top20_CNlow_Hyper,
  "Results/Top20_DE_CNlow_vs_Hypermutated.csv"
)

cat("\nTop genes: CN-high vs CN-low\n")
print(
  top20_CNhigh_CNlow |>
    dplyr::select(
      Gene,
      logFC,
      adj.P.Val
    ) |>
    head(10)
)

cat("\nTop genes: CN-high vs Hypermutated\n")
print(
  top20_CNhigh_Hyper |>
    dplyr::select(
      Gene,
      logFC,
      adj.P.Val
    ) |>
    head(10)
)

cat("\nTop genes: CN-low vs Hypermutated\n")
print(
  top20_CNlow_Hyper |>
    dplyr::select(
      Gene,
      logFC,
      adj.P.Val
    ) |>
    head(10)
)
