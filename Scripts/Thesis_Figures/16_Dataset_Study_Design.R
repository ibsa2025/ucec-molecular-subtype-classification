# ============================================================
# DATASETS, STUDY DESIGN AND MOLECULAR SUBTYPES
# Final improved dissertation figure
# ============================================================

library(ggplot2)
library(grid)

dir.create(
  "Figures/Dataset_Slide",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# COLOUR PALETTE
# ============================================================

navy          <- "#1F5588"
blue_fill     <- "#DDEBFA"
blue_border   <- "#72A6D3"

orange_header <- "#CE751B"
orange_fill   <- "#FCE8D2"
orange_border <- "#E59A4A"

green_fill    <- "#DFF1E7"
green_border  <- "#69B284"

purple_fill   <- "#EEE3F8"
purple_border <- "#AE80D2"

red_row       <- "#FBE6E6"
blue_row      <- "#E7F0F9"
green_row     <- "#E7F5EC"

dark_text     <- "#263238"
medium_text   <- "#566573"
line_colour   <- "#67727C"
grid_colour   <- "#C7D0D8"
white         <- "#FFFFFF"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

add_box <- function(plot,
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                    fill,
                    border,
                    title,
                    subtitle = NULL,
                    title_size = 4.4,
                    subtitle_size = 3.3,
                    title_colour = dark_text,
                    subtitle_colour = medium_text,
                    border_size = 0.9,
                    title_adjust = 0.16) {
  
  centre_x <- (xmin + xmax) / 2
  centre_y <- (ymin + ymax) / 2
  box_height <- ymax - ymin
  
  title_y <- if (is.null(subtitle)) {
    centre_y
  } else {
    centre_y + box_height * title_adjust
  }
  
  subtitle_y <- centre_y - box_height * 0.18
  
  plot <- plot +
    annotate(
      "rect",
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill,
      colour = border,
      linewidth = border_size
    ) +
    annotate(
      "text",
      x = centre_x,
      y = title_y,
      label = title,
      fontface = "bold",
      size = title_size,
      colour = title_colour,
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
        colour = subtitle_colour,
        lineheight = 0.95
      )
  }
  
  plot
}

add_arrow <- function(plot,
                      x,
                      y,
                      xend,
                      yend,
                      colour = line_colour,
                      linewidth = 0.8) {
  
  plot +
    annotate(
      "segment",
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = colour,
      linewidth = linewidth,
      lineend = "round",
      arrow = arrow(
        length = unit(0.17, "cm"),
        type = "closed"
      )
    )
}

add_line <- function(plot,
                     x,
                     y,
                     xend,
                     yend,
                     colour = line_colour,
                     linewidth = 0.8) {
  
  plot +
    annotate(
      "segment",
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = colour,
      linewidth = linewidth,
      lineend = "round"
    )
}

# ============================================================
# EMPTY CANVAS
# ============================================================

dataset_plot <- ggplot() +
  coord_cartesian(
    xlim = c(0, 16),
    ylim = c(0, 10),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(
      fill = white,
      colour = NA
    )
  )

# ============================================================
# MAIN TITLE
# ============================================================

dataset_plot <- dataset_plot +
  annotate(
    "text",
    x = 0.7,
    y = 9.62,
    label = "Datasets, Study Design and Molecular Subtypes",
    hjust = 0,
    fontface = "bold",
    size = 7.6,
    colour = "#111111"
  ) +
  annotate(
    "segment",
    x = 0.7,
    xend = 15.3,
    y = 9.22,
    yend = 9.22,
    colour = "#D5DADD",
    linewidth = 0.7
  )

# ============================================================
# LEFT PANEL: TCGA-UCEC COHORT
# ============================================================

dataset_plot <- dataset_plot +
  annotate(
    "rect",
    xmin = 0.7,
    xmax = 7.45,
    ymin = 5.15,
    ymax = 8.82,
    fill = "#F2F7FC",
    colour = blue_border,
    linewidth = 1
  ) +
  annotate(
    "rect",
    xmin = 0.7,
    xmax = 7.45,
    ymin = 8.08,
    ymax = 8.82,
    fill = navy,
    colour = navy,
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 1.02,
    y = 8.45,
    label = "TCGA-UCEC COHORT",
    hjust = 0,
    fontface = "bold",
    size = 5.3,
    colour = white
  )

dataset_plot <- add_box(
  dataset_plot,
  xmin = 1.02, xmax = 2.90,
  ymin = 7.05, ymax = 7.78,
  fill = blue_fill,
  border = blue_border,
  title = "232",
  subtitle = "RNA-seq tumour samples",
  title_size = 5.8,
  subtitle_size = 2.95
)

dataset_plot <- add_box(
  dataset_plot,
  xmin = 3.15, xmax = 5.05,
  ymin = 7.05, ymax = 7.78,
  fill = blue_fill,
  border = blue_border,
  title = "~59,000",
  subtitle = "genes in raw counts",
  title_size = 5.4,
  subtitle_size = 3.0
)

dataset_plot <- add_box(
  dataset_plot,
  xmin = 5.30, xmax = 7.13,
  ymin = 7.05, ymax = 7.78,
  fill = blue_fill,
  border = blue_border,
  title = "3",
  subtitle = "analysed subtypes",
  title_size = 5.8,
  subtitle_size = 3.1
)

dataset_plot <- dataset_plot +
  annotate(
    "text",
    x = 1.02,
    y = 6.60,
    label = "Clinical metadata included molecular subtype annotations",
    hjust = 0,
    fontface = "bold",
    size = 3.55,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 1.02,
    y = 6.20,
    label = "Exploratory analyses retained the 1,000 most variable genes.",
    hjust = 0,
    size = 3.15,
    colour = medium_text
  ) +
  annotate(
    "text",
    x = 1.02,
    y = 5.82,
    label = "The TCGA cohort was used for:",
    hjust = 0,
    fontface = "bold",
    size = 3.45,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 1.12,
    y = 5.55,
    label = "•  Differential expression",
    hjust = 0,
    size = 3.30,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 1.12,
    y = 5.27,
    label = "•  Exploratory analysis",
    hjust = 0,
    size = 3.30,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 4.05,
    y = 5.55,
    label = "•  Machine-learning development",
    hjust = 0,
    size = 3.10,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 4.05,
    y = 5.27,
    label = "•  GO and KEGG enrichment",
    hjust = 0,
    size = 3.30,
    colour = dark_text
  )

# ============================================================
# RIGHT PANEL: MACHINE LEARNING WORKFLOW
# ============================================================

dataset_plot <- dataset_plot +
  annotate(
    "rect",
    xmin = 7.75,
    xmax = 15.30,
    ymin = 4.55,
    ymax = 8.82,
    fill = "#FEFBF8",
    colour = orange_border,
    linewidth = 1
  ) +
  annotate(
    "rect",
    xmin = 7.75,
    xmax = 15.30,
    ymin = 8.08,
    ymax = 8.82,
    fill = orange_header,
    colour = orange_header,
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 8.08,
    y = 8.45,
    label = "MACHINE LEARNING WORKFLOW",
    hjust = 0,
    fontface = "bold",
    size = 5.0,
    colour = white
  )

dataset_plot <- add_box(
  dataset_plot,
  xmin = 9.65, xmax = 13.40,
  ymin = 7.15, ymax = 7.78,
  fill = blue_fill,
  border = blue_border,
  title = "TCGA-UCEC",
  subtitle = "n = 232 RNA-seq tumour samples",
  title_size = 4.4,
  subtitle_size = 3.0
)

dataset_plot <- add_box(
  dataset_plot,
  xmin = 8.20, xmax = 11.00,
  ymin = 6.05, ymax = 6.76,
  fill = orange_fill,
  border = orange_border,
  title = "Training cohort",
  subtitle = "n = 186 (80%)",
  title_size = 4.15,
  subtitle_size = 3.15
)

dataset_plot <- add_box(
  dataset_plot,
  xmin = 12.05, xmax = 14.85,
  ymin = 6.05, ymax = 6.76,
  fill = green_fill,
  border = green_border,
  title = "Held-out test cohort",
  subtitle = "n = 46 (20%)",
  title_size = 3.80,
  subtitle_size = 3.15
)

dataset_plot <- add_box(
  dataset_plot,
  xmin = 8.20, xmax = 11.00,
  ymin = 5.00, ymax = 5.72,
  fill = "#FFF6EB",
  border = orange_border,
  title = "Feature selection and scaling\n5-fold cross-validation\nLASSO + Random Forest",
  title_size = 3.10
)

dataset_plot <- add_box(
  dataset_plot,
  xmin = 12.05, xmax = 14.85,
  ymin = 5.00, ymax = 5.72,
  fill = "#EDF8F1",
  border = green_border,
  title = "Internal validation",
  subtitle = "Evaluation on unseen test data",
  title_size = 3.95,
  subtitle_size = 2.80
)

# Workflow arrows
dataset_plot <- add_arrow(
  dataset_plot,
  x = 11.52, y = 7.15,
  xend = 9.60, yend = 6.80
)

dataset_plot <- add_arrow(
  dataset_plot,
  x = 11.52, y = 7.15,
  xend = 13.45, yend = 6.80
)

dataset_plot <- add_arrow(
  dataset_plot,
  x = 9.60, y = 6.05,
  xend = 9.60, yend = 5.76
)

dataset_plot <- add_arrow(
  dataset_plot,
  x = 13.45, y = 6.05,
  xend = 13.45, yend = 5.76
)

dataset_plot <- add_line(
  dataset_plot,
  x = 9.60, y = 5.00,
  xend = 9.60, yend = 4.90
)

dataset_plot <- add_line(
  dataset_plot,
  x = 13.45, y = 5.00,
  xend = 13.45, yend = 4.90
)

dataset_plot <- add_line(
  dataset_plot,
  x = 9.60, y = 4.90,
  xend = 13.45, yend = 4.90
)

dataset_plot <- add_arrow(
  dataset_plot,
  x = 11.52, y = 4.90,
  xend = 11.52, yend = 4.73
)

# Enlarged CPTAC box
dataset_plot <- add_box(
  dataset_plot,
  xmin = 9.55, xmax = 13.50,
  ymin = 4.10, ymax = 4.73,
  fill = purple_fill,
  border = purple_border,
  title = "CPTAC-UCEC  •  n = 101",
  subtitle = "Independent external validation\nNCI CPTAC programme (United States)",
  title_size = 3.65,
  subtitle_size = 2.45,
  title_adjust = 0.18
)

# ============================================================
# MOLECULAR SUBTYPES TABLE
# Lowered and vertically compressed
# ============================================================

dataset_plot <- dataset_plot +
  annotate(
    "text",
    x = 0.7,
    y = 3.92,
    label = "Molecular Subtypes",
    hjust = 0,
    fontface = "bold",
    size = 5.5,
    colour = "#111111"
  )

# Table boundaries
table_top    <- 3.62
header_low   <- 3.03
row1_low     <- 2.22
row2_low     <- 1.41
table_bottom <- 0.60

dataset_plot <- dataset_plot +
  annotate(
    "rect",
    xmin = 0.7,
    xmax = 15.30,
    ymin = table_bottom,
    ymax = table_top,
    fill = white,
    colour = grid_colour,
    linewidth = 0.9
  ) +
  annotate(
    "rect",
    xmin = 0.7,
    xmax = 15.30,
    ymin = header_low,
    ymax = table_top,
    fill = navy,
    colour = navy
  )

# Header labels
dataset_plot <- dataset_plot +
  annotate(
    "text",
    x = 1.05,
    y = 3.325,
    label = "Molecular subtype",
    hjust = 0,
    fontface = "bold",
    size = 4.05,
    colour = white
  ) +
  annotate(
    "text",
    x = 4.70,
    y = 3.325,
    label = "Key molecular and clinical characteristics",
    hjust = 0,
    fontface = "bold",
    size = 4.05,
    colour = white
  ) +
  annotate(
    "text",
    x = 13.92,
    y = 3.325,
    label = "TCGA, n (%)",
    fontface = "bold",
    size = 3.95,
    colour = white
  )

# Row backgrounds
dataset_plot <- dataset_plot +
  annotate(
    "rect",
    xmin = 0.7, xmax = 15.30,
    ymin = row1_low, ymax = header_low,
    fill = red_row,
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = 0.7, xmax = 15.30,
    ymin = row2_low, ymax = row1_low,
    fill = blue_row,
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = 0.7, xmax = 15.30,
    ymin = table_bottom, ymax = row2_low,
    fill = green_row,
    colour = NA
  )

# Dividers
dataset_plot <- dataset_plot +
  annotate(
    "segment",
    x = 4.35, xend = 4.35,
    y = table_bottom, yend = table_top,
    colour = grid_colour,
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = 12.55, xend = 12.55,
    y = table_bottom, yend = table_top,
    colour = grid_colour,
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = 0.7, xend = 15.30,
    y = row1_low, yend = row1_low,
    colour = grid_colour,
    linewidth = 0.65
  ) +
  annotate(
    "segment",
    x = 0.7, xend = 15.30,
    y = row2_low, yend = row2_low,
    colour = grid_colour,
    linewidth = 0.65
  )

# Row centres
row1_c <- (header_low + row1_low) / 2
row2_c <- (row1_low + row2_low) / 2
row3_c <- (row2_low + table_bottom) / 2

# Subtype labels
dataset_plot <- dataset_plot +
  annotate(
    "text",
    x = 1.05,
    y = row1_c + 0.12,
    label = "CN-high",
    hjust = 0,
    fontface = "bold",
    size = 4.45,
    colour = "#A12C2C"
  ) +
  annotate(
    "text",
    x = 1.05,
    y = row1_c - 0.14,
    label = "Copy-number high",
    hjust = 0,
    size = 3.30,
    colour = medium_text
  ) +
  annotate(
    "text",
    x = 1.05,
    y = row2_c + 0.12,
    label = "CN-low",
    hjust = 0,
    fontface = "bold",
    size = 4.45,
    colour = "#28699B"
  ) +
  annotate(
    "text",
    x = 1.05,
    y = row2_c - 0.14,
    label = "Copy-number low",
    hjust = 0,
    size = 3.30,
    colour = medium_text
  ) +
  annotate(
    "text",
    x = 1.05,
    y = row3_c + 0.12,
    label = "Hypermutated (MSI+POLE)", 
    hjust = 0,
    fontface = "bold",
    size = 4.45,
    colour = "#277B49"
  ) +
  annotate(
    "text",
    x = 1.05,
    y = row3_c - 0.14,
    label = "High mutation burden",
    hjust = 0,
    size = 3.30,
    colour = medium_text
  )

# Characteristics
dataset_plot <- dataset_plot +
  annotate(
    "text",
    x = 4.70,
    y = row1_c + 0.12,
    label = "Aggressive molecular subtype",
    hjust = 0,
    size = 3.95,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 4.70,
    y = row1_c - 0.14,
    label = "Associated with poorer prognosis",
    hjust = 0,
    size = 3.65,
    colour = medium_text
  ) +
  annotate(
    "text",
    x = 4.70,
    y = row2_c + 0.12,
    label = "Predominantly endometrioid tumours",
    hjust = 0,
    size = 3.95,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 4.70,
    y = row2_c - 0.14,
    label = "Intermediate prognosis",
    hjust = 0,
    size = 3.65,
    colour = medium_text
  ) +
  annotate(
    "text",
    x = 4.70,
    y = row3_c + 0.12,
    label = "High somatic mutation burden",
    hjust = 0,
    size = 3.95,
    colour = dark_text
  ) +
  annotate(
    "text",
    x = 4.70,
    y = row3_c - 0.14,
    label = "Heterogeneous prognosis",
    hjust = 0,
    size = 3.65,
    colour = medium_text
  )

# Counts
dataset_plot <- dataset_plot +
  annotate(
    "text",
    x = 13.92,
    y = row1_c,
    label = "60 (25.9%)",
    fontface = "bold",
    size = 4.05,
    colour = "#A12C2C"
  ) +
  annotate(
    "text",
    x = 13.92,
    y = row2_c,
    label = "90 (38.8%)",
    fontface = "bold",
    size = 4.05,
    colour = "#28699B"
  ) +
  annotate(
    "text",
    x = 13.92,
    y = row3_c,
    label = "82 (35.3%)",
    fontface = "bold",
    size = 4.05,
    colour = "#277B49"
  )

# ============================================================
# SAVE OUTPUTS
# ============================================================

ggsave(
  filename = "Figures/Dataset_Slide/Datasets_StudyDesign_MolecularSubtypes_Final.png",
  plot = dataset_plot,
  width = 13.33,
  height = 7.5,
  dpi = 600,
  bg = "white"
)

if (capabilities("cairo")) {
  ggsave(
    filename = "Figures/Dataset_Slide/Datasets_StudyDesign_MolecularSubtypes_Final.pdf",
    plot = dataset_plot,
    width = 13.33,
    height = 7.5,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
} else {
  ggsave(
    filename = "Figures/Dataset_Slide/Datasets_StudyDesign_MolecularSubtypes_Final.pdf",
    plot = dataset_plot,
    width = 13.33,
    height = 7.5,
    device = "pdf",
    bg = "white"
  )
}

dataset_plot
