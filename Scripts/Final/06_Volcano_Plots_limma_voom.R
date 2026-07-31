# =====================================
# UCEC Project
# 06 - Volcano Plots and DEG Counts
# =====================================

library(tidyverse)
library(patchwork)

dir.create(
  "Figures/DE",
  recursive = TRUE,
  showWarnings = FALSE
)

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

make_volcano <- function(
    de_table,
    panel_label,
    comparison_label
) {
  plot_data <- de_table |>
    dplyr::mutate(
      Regulation = dplyr::case_when(
        adj.P.Val < 0.05 & logFC > 1 ~ "Positive log2FC",
        adj.P.Val < 0.05 & logFC < -1 ~ "Negative log2FC",
        TRUE ~ "Not significant"
      ),
      Negative_Log10_Adjusted_P =
        -log10(
          pmax(
            adj.P.Val,
            .Machine$double.xmin
          )
        )
    )

  label_data <- plot_data |>
    dplyr::filter(
      Regulation != "Not significant"
    ) |>
    dplyr::arrange(adj.P.Val) |>
    dplyr::slice_head(n = 4)

  ggplot(
    plot_data,
    aes(
      x = logFC,
      y = Negative_Log10_Adjusted_P,
      colour = Regulation
    )
  ) +
    geom_point(
      alpha = 0.55,
      size = 0.9
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_text(
      data = label_data,
      aes(label = Gene),
      colour = "black",
      size = 2.7,
      check_overlap = TRUE,
      vjust = -0.4,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = c(
        "Negative log2FC" = "#377EB8",
        "Not significant" = "grey75",
        "Positive log2FC" = "#E41A1C"
      )
    ) +
    theme_classic() +
    labs(
      title = paste0(panel_label, "  ", comparison_label),
      x = "log2 fold change",
      y = "-log10 adjusted P",
      colour = "Gene status"
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      legend.position = "bottom"
    )
}

p_volcano_1 <- make_volcano(
  de_CNhigh_CNlow,
  "A",
  "CN-high vs CN-low"
)

p_volcano_2 <- make_volcano(
  de_CNhigh_Hyper,
  "B",
  "CN-high vs Hypermutated"
)

p_volcano_3 <- make_volcano(
  de_CNlow_Hyper,
  "C",
  "CN-low vs Hypermutated"
)

volcano_combined <- (
  p_volcano_1 +
    p_volcano_2 +
    p_volcano_3 +
    patchwork::plot_layout(
      ncol = 3,
      guides = "collect"
    )
) &
  theme(legend.position = "bottom")

print(volcano_combined)

ggsave(
  "Figures/DE/Volcano_Plots_Combined.png",
  plot = volcano_combined,
  width = 14,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figures/DE/Volcano_Plots_Combined.pdf",
  plot = volcano_combined,
  width = 14,
  height = 5,
  bg = "white"
)

create_direction_counts <- function(
    de_table,
    comparison
) {
  de_table |>
    dplyr::filter(
      adj.P.Val < 0.05,
      abs(logFC) > 1
    ) |>
    dplyr::mutate(
      Direction = ifelse(
        logFC > 1,
        "Positive log2FC",
        "Negative log2FC"
      )
    ) |>
    dplyr::count(
      Direction,
      name = "Gene_Count"
    ) |>
    dplyr::mutate(
      Comparison = comparison
    )
}

deg_counts <- dplyr::bind_rows(
  create_direction_counts(
    de_CNhigh_CNlow,
    "CN-high vs CN-low"
  ),
  create_direction_counts(
    de_CNhigh_Hyper,
    "CN-high vs Hypermutated"
  ),
  create_direction_counts(
    de_CNlow_Hyper,
    "CN-low vs Hypermutated"
  )
)

readr::write_csv(
  deg_counts,
  "Results/DE_Direction_Counts.csv"
)

p_deg_counts <- ggplot(
  deg_counts,
  aes(
    x = Comparison,
    y = Gene_Count,
    fill = Direction
  )
) +
  geom_col(
    position = "dodge"
  ) +
  geom_text(
    aes(label = Gene_Count),
    position = position_dodge(width = 0.9),
    vjust = -0.25,
    size = 3.5
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of significant genes",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    legend.position = "top"
  )

print(p_deg_counts)

ggsave(
  "Figures/DE/DEG_Counts_by_Comparison.png",
  plot = p_deg_counts,
  width = 9,
  height = 6,
  dpi = 600,
  bg = "white"
)
