# ============================================================
# FINAL IMPROVED FEATURE-SELECTION WORKFLOW
# Fixes:
# - wider canvas
# - wider side boxes
# - smaller text where needed
# - taller bottom boxes
# - 11-gene list placed in its own box
# ============================================================

library(ggplot2)
library(grid)

dir.create(
  "Figures/Feature_Selection",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# COLOURS
# ============================================================

blue_fill     <- "#DDEBFA"
blue_border   <- "#72A6D3"

orange_fill   <- "#FCE8D2"
orange_border <- "#E59A4A"

purple_fill   <- "#EEE3F8"
purple_border <- "#AE80D2"

green_fill    <- "#DFF1E7"
green_border  <- "#69B284"

grey_fill     <- "#F3F5F7"
grey_border   <- "#AEB8C2"

light_fill    <- "#EEF3F8"
light_border  <- "#B7C7D8"

dark_text     <- "#263238"
medium_text   <- "#566573"
arrow_colour  <- "#68747E"
white         <- "#FFFFFF"

# ============================================================
# HELPER FUNCTION: BOX
# ============================================================

add_workflow_box <- function(plot,
                             xmin,
                             xmax,
                             ymin,
                             ymax,
                             fill,
                             border,
                             title,
                             subtitle = NULL,
                             title_size = 5,
                             subtitle_size = 3.5) {
  
  centre_x <- (xmin + xmax) / 2
  centre_y <- (ymin + ymax) / 2
  box_height <- ymax - ymin
  
  title_y <- if (is.null(subtitle)) {
    centre_y
  } else {
    centre_y + box_height * 0.17
  }
  
  subtitle_y <- centre_y - box_height * 0.16
  
  plot <- plot +
    annotate(
      "rect",
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill,
      colour = border,
      linewidth = 1
    ) +
    annotate(
      "text",
      x = centre_x,
      y = title_y,
      label = title,
      fontface = "bold",
      size = title_size,
      colour = dark_text,
      lineheight = 0.95
    )
  
  if (!is.null(subtitle)) {
    plot <- plot +
      annotate(
        "text",
        x = centre_x,
        y = subtitle_y,
        label = subtitle,
        size = subtitle_size,
        colour = medium_text,
        lineheight = 1.05
      )
  }
  
  plot
}

# ============================================================
# HELPER FUNCTION: ARROW
# ============================================================

add_arrow <- function(plot, x, y, xend, yend) {
  plot +
    annotate(
      "segment",
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = arrow_colour,
      linewidth = 0.9,
      lineend = "round",
      arrow = arrow(
        length = unit(0.18, "cm"),
        type = "closed"
      )
    )
}

# ============================================================
# EMPTY CANVAS
# ============================================================

workflow_plot <- ggplot() +
  coord_cartesian(
    xlim = c(0, 14),
    ylim = c(0, 16),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = white, colour = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# ============================================================
# TITLE
# ============================================================

workflow_plot <- workflow_plot +
  annotate(
    "text",
    x = 7,
    y = 15.35,
    label = "Feature-Selection Strategy",
    fontface = "bold",
    size = 7.2,
    colour = "#111111"
  ) +
  annotate(
    "text",
    x = 7,
    y = 14.8,
    label = "Prioritisation of shared candidate biomarkers using complementary machine-learning approaches",
    size = 3.9,
    colour = medium_text
  ) +
  annotate(
    "segment",
    x = 1.1,
    xend = 12.9,
    y = 14.45,
    yend = 14.45,
    colour = "#D5DADD",
    linewidth = 0.7
  )

# ============================================================
# STEP 0: INPUT
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 2.2,
  xmax = 11.8,
  ymin = 12.55,
  ymax = 13.85,
  fill = light_fill,
  border = light_border,
  title = "1,000 Most Variable Training-Set Genes",
  subtitle = "Common predictor set used independently by both models",
  title_size = 5.3,
  subtitle_size = 3.7
)

workflow_plot <- add_arrow(workflow_plot, 7, 12.55, 4.1, 11.6)
workflow_plot <- add_arrow(workflow_plot, 7, 12.55, 9.9, 11.6)

# ============================================================
# STEP 1A: LASSO
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 0.9,
  xmax = 6.6,
  ymin = 10.0,
  ymax = 11.4,
  fill = blue_fill,
  border = blue_border,
  title = "Multinomial LASSO at lambda.min",
  subtitle = "56 genes retained with\nnon-zero coefficients",
  title_size = 4.8,
  subtitle_size = 3.45
)

# ============================================================
# STEP 1B: RANDOM FOREST
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 7.4,
  xmax = 13.1,
  ymin = 10.0,
  ymax = 11.4,
  fill = orange_fill,
  border = orange_border,
  title = "Random Forest Feature Ranking",
  subtitle = "Top 30 genes retained using\npermutation importance",
  title_size = 4.7,
  subtitle_size = 3.35
)

workflow_plot <- add_arrow(workflow_plot, 4.1, 10.0, 4.1, 9.2)
workflow_plot <- add_arrow(workflow_plot, 9.9, 10.0, 9.9, 9.2)

workflow_plot <- add_arrow(workflow_plot, 4.1, 9.2, 5.8, 8.5)
workflow_plot <- add_arrow(workflow_plot, 9.9, 9.2, 8.2, 8.5)

# ============================================================
# STEP 2: INTERSECTION
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 2.3,
  xmax = 11.7,
  ymin = 7.3,
  ymax = 8.65,
  fill = purple_fill,
  border = purple_border,
  title = "Feature-Set Intersection",
  subtitle = "Comparison of independently selected LASSO and Random Forest genes",
  title_size = 5.0,
  subtitle_size = 3.45
)

workflow_plot <- add_arrow(workflow_plot, 7, 7.3, 7, 6.45)

# ============================================================
# STEP 3: SHARED BIOMARKERS
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 2.3,
  xmax = 11.7,
  ymin = 5.0,
  ymax = 6.35,
  fill = green_fill,
  border = green_border,
  title = "11 Shared Candidate Biomarkers",
  subtitle = "Genes supported by both machine-learning approaches",
  title_size = 5.1,
  subtitle_size = 3.55
)

workflow_plot <- add_arrow(workflow_plot, 7, 5.0, 7, 4.2)

# ============================================================
# STEP 4: EVIDENCE-BASED PRIORITISATION
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 2.3,
  xmax = 11.7,
  ymin = 2.95,
  ymax = 4.15,
  fill = grey_fill,
  border = grey_border,
  title = "Evidence-Based Prioritisation",
  subtitle = "Integrated with differential-expression and biological evidence",
  title_size = 4.9,
  subtitle_size = 3.45
)

workflow_plot <- add_arrow(workflow_plot, 7, 2.95, 7, 2.2)

# ============================================================
# STEP 5: BIOLOGICAL INTERPRETATION
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 2.3,
  xmax = 11.7,
  ymin = 1.05,
  ymax = 2.45,
  fill = grey_fill,
  border = grey_border,
  title = "Biological Interpretation",
  subtitle = paste(
    "Subtype-associated expression",
    "Functional and pathway interpretation",
    "Published evidence and potential clinical relevance",
    sep = "\n"
  ),
  title_size = 4.8,
  subtitle_size = 3.15
)

# ============================================================
# STEP 6: GENE LIST BOX
# ============================================================

workflow_plot <- add_workflow_box(
  workflow_plot,
  xmin = 2.3,
  xmax = 11.7,
  ymin = 0.15,
  ymax = 0.9,
  fill = white,
  border = white,
  title = NULL,
  subtitle = NULL
)

workflow_plot <- workflow_plot +
  annotate(
    "text",
    x = 7,
    y = 0.52,
    label = paste(
      "L1CAM • AC009237.3 • ACTL8 • CLDN9 • AC092969.1 • HIF3A",
      "SDR42E1 • AC013724.1 • SLFN13 • RGS22 • SMIM31",
      sep = "\n"
    ),
    size = 3.25,
    colour = medium_text,
    lineheight = 1.05,
    fontface = "plain"
  )

# ============================================================
# DISPLAY
# ============================================================

workflow_plot

# ============================================================
# SAVE
# ============================================================

ggsave(
  filename = "Figures/Feature_Selection/Feature_Selection_Workflow_Final.png",
  plot = workflow_plot,
  width = 9.2,
  height = 10.4,
  dpi = 600,
  bg = "white"
)

if (capabilities("cairo")) {
  ggsave(
    filename = "Figures/Feature_Selection/Feature_Selection_Workflow_Final.pdf",
    plot = workflow_plot,
    width = 9.2,
    height = 10.4,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
} else {
  ggsave(
    filename = "Figures/Feature_Selection/Feature_Selection_Workflow_Final.pdf",
    plot = workflow_plot,
    width = 9.2,
    height = 10.4,
    device = "pdf",
    bg = "white"
  )
}
