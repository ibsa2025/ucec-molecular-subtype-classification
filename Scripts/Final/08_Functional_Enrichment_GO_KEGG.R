# =====================================
# UCEC Project
# 08 - Functional Enrichment Analysis
# GO Biological Process and KEGG
# =====================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(patchwork)

dir.create(
  "Results/Enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/Enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

sig_CNhigh_CNlow <- readr::read_csv(
  "Results/Sig_DE_CNhigh_vs_CNlow_limma_voom.csv",
  show_col_types = FALSE
)

sig_CNhigh_Hyper <- readr::read_csv(
  "Results/Sig_DE_CNhigh_vs_Hypermutated_limma_voom.csv",
  show_col_types = FALSE
)

sig_CNlow_Hyper <- readr::read_csv(
  "Results/Sig_DE_CNlow_vs_Hypermutated_limma_voom.csv",
  show_col_types = FALSE
)

convert_symbols_to_entrez <- function(gene_symbols) {
  clusterProfiler::bitr(
    unique(gene_symbols),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) |>
    dplyr::distinct(
      ENTREZID,
      .keep_all = TRUE
    )
}

CNhigh_CNlow_entrez <- convert_symbols_to_entrez(
  sig_CNhigh_CNlow$Gene
)

CNhigh_Hyper_entrez <- convert_symbols_to_entrez(
  sig_CNhigh_Hyper$Gene
)

CNlow_Hyper_entrez <- convert_symbols_to_entrez(
  sig_CNlow_Hyper$Gene
)

run_go_bp <- function(entrez_table) {
  clusterProfiler::enrichGO(
    gene = entrez_table$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

run_kegg <- function(entrez_table) {
  clusterProfiler::enrichKEGG(
    gene = entrez_table$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
}

GO_CNhigh_CNlow <- run_go_bp(
  CNhigh_CNlow_entrez
)

GO_CNhigh_Hyper <- run_go_bp(
  CNhigh_Hyper_entrez
)

GO_CNlow_Hyper <- run_go_bp(
  CNlow_Hyper_entrez
)

KEGG_CNhigh_CNlow <- run_kegg(
  CNhigh_CNlow_entrez
)

KEGG_CNhigh_Hyper <- run_kegg(
  CNhigh_Hyper_entrez
)

KEGG_CNlow_Hyper <- run_kegg(
  CNlow_Hyper_entrez
)

readr::write_csv(
  as.data.frame(GO_CNhigh_CNlow),
  "Results/Enrichment/GO_BP_CNhigh_vs_CNlow.csv"
)

readr::write_csv(
  as.data.frame(GO_CNhigh_Hyper),
  "Results/Enrichment/GO_BP_CNhigh_vs_Hypermutated.csv"
)

readr::write_csv(
  as.data.frame(GO_CNlow_Hyper),
  "Results/Enrichment/GO_BP_CNlow_vs_Hypermutated.csv"
)

readr::write_csv(
  as.data.frame(KEGG_CNhigh_CNlow),
  "Results/Enrichment/KEGG_CNhigh_vs_CNlow.csv"
)

readr::write_csv(
  as.data.frame(KEGG_CNhigh_Hyper),
  "Results/Enrichment/KEGG_CNhigh_vs_Hypermutated.csv"
)

readr::write_csv(
  as.data.frame(KEGG_CNlow_Hyper),
  "Results/Enrichment/KEGG_CNlow_vs_Hypermutated.csv"
)

make_enrichment_plot <- function(
    enrichment_object,
    title_text,
    categories = 12
) {
  enrichplot::dotplot(
    enrichment_object,
    showCategory = categories,
    title = title_text
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      axis.text.y = element_text(size = 8)
    )
}

p_go_1 <- make_enrichment_plot(
  GO_CNhigh_CNlow,
  "A  CN-high vs CN-low"
)

p_go_2 <- make_enrichment_plot(
  GO_CNhigh_Hyper,
  "B  CN-high vs Hypermutated"
)

p_go_3 <- make_enrichment_plot(
  GO_CNlow_Hyper,
  "C  CN-low vs Hypermutated"
)

go_combined <- (
  p_go_1 |
    p_go_2
) /
  p_go_3

print(go_combined)

ggsave(
  "Figures/Enrichment/GO_BP_Combined.png",
  plot = go_combined,
  width = 13,
  height = 11,
  dpi = 600,
  bg = "white"
)

p_kegg_1 <- make_enrichment_plot(
  KEGG_CNhigh_CNlow,
  "A  CN-high vs CN-low"
)

p_kegg_2 <- make_enrichment_plot(
  KEGG_CNhigh_Hyper,
  "B  CN-high vs Hypermutated"
)

p_kegg_3 <- make_enrichment_plot(
  KEGG_CNlow_Hyper,
  "C  CN-low vs Hypermutated"
)

kegg_combined <- (
  p_kegg_1 |
    p_kegg_2
) /
  p_kegg_3

print(kegg_combined)

ggsave(
  "Figures/Enrichment/KEGG_Combined.png",
  plot = kegg_combined,
  width = 13,
  height = 11,
  dpi = 600,
  bg = "white"
)

enrichment_summary <- tibble::tibble(
  Comparison = c(
    "CN-high vs CN-low",
    "CN-high vs Hypermutated",
    "CN-low vs Hypermutated"
  ),
  GO_BP_Terms = c(
    nrow(as.data.frame(GO_CNhigh_CNlow)),
    nrow(as.data.frame(GO_CNhigh_Hyper)),
    nrow(as.data.frame(GO_CNlow_Hyper))
  ),
  KEGG_Pathways = c(
    nrow(as.data.frame(KEGG_CNhigh_CNlow)),
    nrow(as.data.frame(KEGG_CNhigh_Hyper)),
    nrow(as.data.frame(KEGG_CNlow_Hyper))
  )
)

readr::write_csv(
  enrichment_summary,
  "Results/Enrichment/Enrichment_Summary.csv"
)

print(enrichment_summary)
