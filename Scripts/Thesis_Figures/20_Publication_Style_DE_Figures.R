# ============================================================
# UCEC Project
# Publication-style Differential Expression Figures
#
# Outputs:
# 1. Three-panel volcano figure
# 2. Significant DEG count summary
# 3. Three-panel MA figure
# 4. DEG summary tables
# ============================================================


# ============================================================
# 1. Project-relative execution
# ============================================================

# Run this script from the root of the UCEC RStudio project.
# No setwd() call is used so the workflow remains portable.


# ============================================================
# 2. Load packages
# ============================================================

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)


# Optional package for better gene-label placement.
# The script will still work if ggrepel is unavailable.

has_ggrepel <- requireNamespace(
  "ggrepel",
  quietly = TRUE
)


# ============================================================
# 3. Create output folders
# ============================================================

dir.create(
  "Figures/Differential_Expression",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Differential_Expression",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. Define differential-expression result files
# ============================================================

de_files <- c(
  CNhigh_vs_CNlow =
    "Results/DE_CNhigh_vs_CNlow_limma_voom.csv",
  
  CNhigh_vs_Hypermutated =
    "Results/DE_CNhigh_vs_Hypermutated_limma_voom.csv",
  
  CNlow_vs_Hypermutated =
    "Results/DE_CNlow_vs_Hypermutated_limma_voom.csv"
)


# Check that all files exist

file_check <- file.exists(de_files)

file_check

if (!all(file_check)) {
  
  stop(
    paste0(
      "The following DE result files were not found:\n",
      paste(
        names(de_files)[!file_check],
        collapse = "\n"
      )
    )
  )
}


# ============================================================
# 5. Set statistical thresholds
# ============================================================

fdr_threshold <- 0.05

logfc_threshold <- 1

genes_to_label_per_direction <- 4


# ============================================================
# 6. Function to load and prepare each DE comparison
# ============================================================

prepare_de_results <- function(
    file_path,
    comparison_name,
    comparison_display
) {
  
  de_data <- readr::read_csv(
    file_path,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "Gene",
    "logFC",
    "adj.P.Val"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(de_data)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        comparison_name,
        " is missing the following columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
  
  # Replace exact zero FDR values with a very small value
  # so that -log10 can be calculated.
  
  smallest_positive_fdr <- suppressWarnings(
    min(
      de_data$adj.P.Val[
        de_data$adj.P.Val > 0 &
          !is.na(de_data$adj.P.Val)
      ],
      na.rm = TRUE
    )
  )
  
  if (
    is.infinite(smallest_positive_fdr) |
    is.na(smallest_positive_fdr)
  ) {
    
    smallest_positive_fdr <- 1e-300
  }
  
  de_data %>%
    dplyr::mutate(
      Comparison = comparison_name,
      
      Comparison_Display = comparison_display,
      
      Adjusted_P_for_Plot = dplyr::case_when(
        is.na(adj.P.Val) ~ NA_real_,
        adj.P.Val <= 0 ~ smallest_positive_fdr / 10,
        TRUE ~ adj.P.Val
      ),
      
      Negative_Log10_FDR =
        -log10(Adjusted_P_for_Plot),
      
      Direction = dplyr::case_when(
        adj.P.Val < fdr_threshold &
          logFC > logfc_threshold ~
          "Positive log2FC",
        
        adj.P.Val < fdr_threshold &
          logFC < -logfc_threshold ~
          "Negative log2FC",
        
        TRUE ~
          "Not significant"
      ),
      
      Direction = factor(
        Direction,
        levels = c(
          "Not significant",
          "Negative log2FC",
          "Positive log2FC"
        )
      )
    ) %>%
    dplyr::filter(
      !is.na(Gene),
      !is.na(logFC),
      !is.na(Negative_Log10_FDR),
      is.finite(logFC),
      is.finite(Negative_Log10_FDR)
    )
}


# ============================================================
# 7. Load all three comparisons
# ============================================================

de_CNhigh_CNlow <- prepare_de_results(
  file_path =
    de_files["CNhigh_vs_CNlow"],
  
  comparison_name =
    "CNhigh_vs_CNlow",
  
  comparison_display =
    "CN-high vs CN-low"
)

de_CNhigh_Hyper <- prepare_de_results(
  file_path =
    de_files["CNhigh_vs_Hypermutated"],
  
  comparison_name =
    "CNhigh_vs_Hypermutated",
  
  comparison_display =
    "CN-high vs Hypermutated"
)

de_CNlow_Hyper <- prepare_de_results(
  file_path =
    de_files["CNlow_vs_Hypermutated"],
  
  comparison_name =
    "CNlow_vs_Hypermutated",
  
  comparison_display =
    "CN-low vs Hypermutated"
)


# Combine all results

de_all <- dplyr::bind_rows(
  de_CNhigh_CNlow,
  de_CNhigh_Hyper,
  de_CNlow_Hyper
)


de_all <- de_all %>%
  dplyr::mutate(
    Comparison_Display = factor(
      Comparison_Display,
      levels = c(
        "CN-high vs CN-low",
        "CN-high vs Hypermutated",
        "CN-low vs Hypermutated"
      )
    )
  )


# ============================================================
# 8. Check the imported results
# ============================================================

dim(de_CNhigh_CNlow)

dim(de_CNhigh_Hyper)

dim(de_CNlow_Hyper)

table(
  de_all$Comparison_Display,
  de_all$Direction,
  useNA = "ifany"
)


# ============================================================
# 9. Create DEG count summary
# ============================================================

deg_count_summary <- de_all %>%
  dplyr::filter(
    Direction != "Not significant"
  ) %>%
  dplyr::count(
    Comparison_Display,
    Direction,
    name = "Number_of_Genes"
  ) %>%
  dplyr::arrange(
    Comparison_Display,
    Direction
  )


deg_total_summary <- deg_count_summary %>%
  dplyr::group_by(
    Comparison_Display
  ) %>%
  dplyr::summarise(
    Total_Significant_DEGs =
      sum(Number_of_Genes),
    
    .groups = "drop"
  )


deg_count_summary

deg_total_summary


# Save DEG counts

readr::write_csv(
  deg_count_summary,
  paste0(
    "Results/Differential_Expression/",
    "Significant_DEG_Counts_by_Comparison.csv"
  )
)

readr::write_csv(
  deg_total_summary,
  paste0(
    "Results/Differential_Expression/",
    "Total_Significant_DEGs_by_Comparison.csv"
  )
)


# ============================================================
# 10. Select genes for volcano-plot labels
# ============================================================

select_volcano_labels <- function(de_data) {
  
  de_data %>%
    dplyr::filter(
      Direction != "Not significant"
    ) %>%
    dplyr::group_by(
      Direction
    ) %>%
    dplyr::arrange(
      adj.P.Val,
      dplyr::desc(abs(logFC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(
      n = genes_to_label_per_direction
    ) %>%
    dplyr::ungroup()
}


label_CNhigh_CNlow <- select_volcano_labels(
  de_CNhigh_CNlow
)

label_CNhigh_Hyper <- select_volcano_labels(
  de_CNhigh_Hyper
)

label_CNlow_Hyper <- select_volcano_labels(
  de_CNlow_Hyper
)


labelled_genes_all <- dplyr::bind_rows(
  label_CNhigh_CNlow,
  label_CNhigh_Hyper,
  label_CNlow_Hyper
) %>%
  dplyr::select(
    Comparison_Display,
    Gene,
    logFC,
    adj.P.Val,
    Direction
  )


readr::write_csv(
  labelled_genes_all,
  paste0(
    "Results/Differential_Expression/",
    "Volcano_Plot_Labelled_Genes.csv"
  )
)


# ============================================================
# 11. Define volcano colours
# ============================================================

volcano_colours <- c(
  "Not significant" = "#BDBDBD",
  "Negative log2FC" = "#3B82F6",
  "Positive log2FC" = "#D55E00"
)


# ============================================================
# 12. Function to create one volcano plot
# ============================================================

create_volcano_plot <- function(
    de_data,
    label_data,
    panel_title
) {
  
  volcano_plot <- ggplot(
    de_data,
    aes(
      x = logFC,
      y = Negative_Log10_FDR,
      colour = Direction
    )
  ) +
    
    geom_point(
      size = 1.25,
      alpha = 0.58
    ) +
    
    geom_vline(
      xintercept = c(
        -logfc_threshold,
        logfc_threshold
      ),
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey35"
    ) +
    
    geom_hline(
      yintercept =
        -log10(fdr_threshold),
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey35"
    ) +
    
    scale_colour_manual(
      values = volcano_colours,
      breaks = c(
        "Positive log2FC",
        "Negative log2FC",
        "Not significant"
      ),
      labels = c(
        "Positive log2FC",
        "Negative log2FC",
        "Not significant"
      ),
      drop = FALSE
    ) +
    
    labs(
      title = panel_title,
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~adjusted~italic(P)),
      colour = NULL
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12.5,
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 10.5
      ),
      
      axis.text = element_text(
        size = 9.5
      ),
      
      legend.position = "bottom",
      
      legend.text = element_text(
        size = 9.5
      ),
      
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
  
  # Use ggrepel when available.
  # Otherwise use ordinary labels with overlap checking.
  
  if (has_ggrepel) {
    
    volcano_plot <- volcano_plot +
      ggrepel::geom_text_repel(
        data = label_data,
        aes(
          label = Gene
        ),
        size = 3,
        colour = "black",
        box.padding = 0.35,
        point.padding = 0.2,
        min.segment.length = 0,
        segment.colour = "grey45",
        segment.linewidth = 0.35,
        max.overlaps = Inf,
        show.legend = FALSE,
        seed = 123
      )
    
  } else {
    
    volcano_plot <- volcano_plot +
      geom_text(
        data = label_data,
        aes(
          label = Gene
        ),
        size = 2.8,
        colour = "black",
        check_overlap = TRUE,
        nudge_y = 0.5,
        show.legend = FALSE
      )
  }
  
  volcano_plot
}


# ============================================================
# 13. Create the three volcano panels
# ============================================================

p_volcano_A <- create_volcano_plot(
  de_data =
    de_CNhigh_CNlow,
  
  label_data =
    label_CNhigh_CNlow,
  
  panel_title =
    "CN-high vs CN-low"
)


p_volcano_B <- create_volcano_plot(
  de_data =
    de_CNhigh_Hyper,
  
  label_data =
    label_CNhigh_Hyper,
  
  panel_title =
    "CN-high vs Hypermutated"
)


p_volcano_C <- create_volcano_plot(
  de_data =
    de_CNlow_Hyper,
  
  label_data =
    label_CNlow_Hyper,
  
  panel_title =
    "CN-low vs Hypermutated"
)


# ============================================================
# 14. Combine volcano plots into one journal-style figure
# ============================================================

volcano_combined <- (
  p_volcano_A +
    p_volcano_B +
    p_volcano_C
) +
  
  patchwork::plot_layout(
    ncol = 3,
    guides = "collect",
    widths = c(
      1,
      1,
      1
    )
  ) +
  
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    legend.position = "bottom",
    
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
  )


volcano_combined


# ============================================================
# 15. Save combined volcano figure
# ============================================================

ggsave(
  filename = paste0(
    "Figures/Differential_Expression/",
    "Combined_Volcano_Plots_600dpi.png"
  ),
  plot = volcano_combined,
  width = 15,
  height = 5.6,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = paste0(
    "Figures/Differential_Expression/",
    "Combined_Volcano_Plots.pdf"
  ),
  plot = volcano_combined,
  width = 15,
  height = 5.6,
  units = "in",
  bg = "white"
)


# ============================================================
# 16. Prepare DEG-count plot data
# ============================================================

deg_count_plot_data <- deg_count_summary %>%
  dplyr::mutate(
    Direction = factor(
      Direction,
      levels = c(
        "Positive log2FC",
        "Negative log2FC"
      )
    )
  )


# ============================================================
# 17. Create DEG count summary figure
# ============================================================

p_deg_counts <- ggplot(
  deg_count_plot_data,
  aes(
    x = Comparison_Display,
    y = Number_of_Genes,
    fill = Direction
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.68
  ) +
  
  geom_text(
    aes(
      label = Number_of_Genes
    ),
    position = position_dodge(
      width = 0.75
    ),
    vjust = -0.35,
    size = 3.5
  ) +
  
  scale_fill_manual(
    values = c(
      "Positive log2FC" = "#D55E00",
      "Negative log2FC" = "#3B82F6"
    ),
    labels = c(
      "Positive log2FC",
      "Negative log2FC"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Number of significant genes",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 15,
      hjust = 1,
      size = 10
    ),
    
    axis.title.y = element_text(
      size = 10.5
    ),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      8,
      10,
      8,
      8
    )
  )


p_deg_counts


# ============================================================
# 18. Save DEG count summary figure
# ============================================================

ggsave(
  filename = paste0(
    "Figures/Differential_Expression/",
    "Significant_DEG_Counts_600dpi.png"
  ),
  plot = p_deg_counts,
  width = 9,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = paste0(
    "Figures/Differential_Expression/",
    "Significant_DEG_Counts.pdf"
  ),
  plot = p_deg_counts,
  width = 9,
  height = 5.5,
  units = "in",
  bg = "white"
)


# ============================================================
# 19. Check whether MA plots can be created
# ============================================================

ma_plot_available <- all(
  c(
    "AveExpr"
  ) %in%
    names(de_CNhigh_CNlow)
) &&
  all(
    c(
      "AveExpr"
    ) %in%
      names(de_CNhigh_Hyper)
  ) &&
  all(
    c(
      "AveExpr"
    ) %in%
      names(de_CNlow_Hyper)
  )


# ============================================================
# 20. Function to create an MA plot
# ============================================================

create_ma_plot <- function(
    de_data,
    panel_title
) {
  
  ggplot(
    de_data,
    aes(
      x = AveExpr,
      y = logFC,
      colour = Direction
    )
  ) +
    
    geom_point(
      size = 1.15,
      alpha = 0.55
    ) +
    
    geom_hline(
      yintercept = 0,
      linewidth = 0.45,
      colour = "grey25"
    ) +
    
    geom_hline(
      yintercept = c(
        -logfc_threshold,
        logfc_threshold
      ),
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey40"
    ) +
    
    scale_colour_manual(
      values = volcano_colours,
      breaks = c(
        "Positive log2FC",
        "Negative log2FC",
        "Not significant"
      ),
      drop = FALSE
    ) +
    
    labs(
      title = panel_title,
      x = "Average log-expression",
      y = expression(log[2]~fold~change),
      colour = NULL
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12.5,
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 10.5
      ),
      
      axis.text = element_text(
        size = 9.5
      ),
      
      legend.position = "bottom",
      
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
}


# ============================================================
# 21. Create and save combined MA plots
# ============================================================

if (ma_plot_available) {
  
  p_ma_A <- create_ma_plot(
    de_CNhigh_CNlow,
    "CN-high vs CN-low"
  )
  
  p_ma_B <- create_ma_plot(
    de_CNhigh_Hyper,
    "CN-high vs Hypermutated"
  )
  
  p_ma_C <- create_ma_plot(
    de_CNlow_Hyper,
    "CN-low vs Hypermutated"
  )
  
  
  ma_combined <- (
    p_ma_A +
      p_ma_B +
      p_ma_C
  ) +
    
    patchwork::plot_layout(
      ncol = 3,
      guides = "collect"
    ) +
    
    patchwork::plot_annotation(
      tag_levels = "A"
    ) &
    
    theme(
      legend.position = "bottom",
      
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  
  
  ma_combined
  
  
  ggsave(
    filename = paste0(
      "Figures/Differential_Expression/",
      "Combined_MA_Plots_600dpi.png"
    ),
    plot = ma_combined,
    width = 15,
    height = 5.6,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  
  ggsave(
    filename = paste0(
      "Figures/Differential_Expression/",
      "Combined_MA_Plots.pdf"
    ),
    plot = ma_combined,
    width = 15,
    height = 5.6,
    units = "in",
    bg = "white"
  )
  
} else {
  
  message(
    paste0(
      "MA plots were not produced because the DE files ",
      "do not contain an AveExpr column."
    )
  )
}


# ============================================================
# 22. Create top-10 DEG summary for each comparison
# ============================================================

top10_degs_each_comparison <- de_all %>%
  dplyr::filter(
    Direction != "Not significant"
  ) %>%
  dplyr::group_by(
    Comparison_Display
  ) %>%
  dplyr::arrange(
    adj.P.Val,
    dplyr::desc(abs(logFC)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = 10
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Comparison_Display,
    Gene,
    logFC,
    adj.P.Val,
    Direction
  )


readr::write_csv(
  top10_degs_each_comparison,
  paste0(
    "Results/Differential_Expression/",
    "Top10_Significant_DEGs_Each_Comparison.csv"
  )
)


# ============================================================
# 23. Confirm saved outputs
# ============================================================

list.files(
  "Figures/Differential_Expression"
)

list.files(
  "Results/Differential_Expression"
)


file.info(
  list.files(
    "Figures/Differential_Expression",
    full.names = TRUE
  )
)[
  ,
  c(
    "size",
    "mtime"
  )
]
