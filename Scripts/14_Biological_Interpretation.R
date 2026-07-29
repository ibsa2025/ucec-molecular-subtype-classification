# =====================================
# UCEC Project
# 14 - Biological Interpretation of Candidate Biomarkers
# =====================================

# Aim:
# To identify and prioritise biologically relevant candidate biomarkers
# supported by the machine-learning and differential-expression analyses.
#
# This script:
# 1. Loads genes selected by LASSO
# 2. Loads Random Forest gene importance
# 3. Identifies genes supported by both models
# 4. Adds differential-expression evidence
# 5. Adds basic human gene annotations
# 6. Produces a consensus biomarker ranking
# 7. Visualises expression across molecular subtypes
# 8. Creates a literature-review template for the thesis
#
# This script is for biological interpretation and visualisation.
# It does not retrain or retest the classifiers.


# =====================================
# 1. Load packages
# =====================================

library(tidyverse)
library(edgeR)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)


# =====================================
# 2. Create output folders
# =====================================

dir.create(
  "Results/Biological_Interpretation",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/Biological_Interpretation",
  recursive = TRUE,
  showWarnings = FALSE
)


# =====================================
# 3. Load machine-learning results
# =====================================

lasso_genes <- read_csv(
  "Results/ML/LASSO_Selected_Genes_lambda_min.csv",
  show_col_types = FALSE
)

rf_importance <- read_csv(
  "Results/ML/RandomForest_All_Gene_Importance.csv",
  show_col_types = FALSE
)

shared_genes_saved <- read_csv(
  "Results/ML/LASSO_RandomForest_Shared_Genes.csv",
  show_col_types = FALSE
)


# Inspect imported files

dim(lasso_genes)
dim(rf_importance)

head(lasso_genes)
head(rf_importance)

colnames(lasso_genes)
colnames(rf_importance)


# =====================================
# 4. Summarise LASSO evidence by gene
# =====================================

# A gene can have a non-zero coefficient for more than one subtype.
#
# For each gene, retain:
# - largest absolute coefficient
# - subtype associated with that coefficient
# - signed coefficient
# - number of subtype models in which it appeared

lasso_gene_summary <- lasso_genes %>%
  group_by(Gene) %>%
  arrange(
    desc(Absolute_Coefficient),
    .by_group = TRUE
  ) %>%
  summarise(
    LASSO_Max_Absolute_Coefficient =
      dplyr::first(Absolute_Coefficient),
    
    LASSO_Coefficient =
      dplyr::first(Coefficient),
    
    LASSO_Associated_Subtype =
      dplyr::first(Class),
    
    LASSO_Number_of_Classes =
      n_distinct(Class),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(LASSO_Max_Absolute_Coefficient)
  )
# =====================================
# 5. Rank Random Forest genes
# =====================================

rf_gene_summary <- rf_importance %>%
  dplyr::select(
    Gene,
    RF_Importance = Importance
  ) %>%
  dplyr::arrange(
    dplyr::desc(RF_Importance)
  ) %>%
  dplyr::mutate(
    RF_Rank = dplyr::row_number()
  )

head(rf_gene_summary, 20)

# =====================================
# 6. Combine LASSO and Random Forest evidence
# =====================================

ml_gene_evidence <- full_join(
  lasso_gene_summary,
  rf_gene_summary,
  by = "Gene"
) %>%
  dplyr::mutate(
    Selected_by_LASSO =
      !is.na(LASSO_Max_Absolute_Coefficient),
    
    Selected_by_RF_Top30 =
      !is.na(RF_Rank) & RF_Rank <= 30,
    
    Selected_by_Both =
      Selected_by_LASSO &
      Selected_by_RF_Top30,
    
    ML_Method_Count =
      as.integer(Selected_by_LASSO) +
      as.integer(Selected_by_RF_Top30)
  )

table(
  ml_gene_evidence$Selected_by_Both,
  useNA = "ifany"
)


# =====================================
# 7. Load differential-expression results
# =====================================

de_CNhigh_CNlow <- read_csv(
  "Results/DE_CNhigh_vs_CNlow_limma_voom.csv",
  show_col_types = FALSE
)

de_CNhigh_Hyper <- read_csv(
  "Results/DE_CNhigh_vs_Hypermutated_limma_voom.csv",
  show_col_types = FALSE
)

de_CNlow_Hyper <- read_csv(
  "Results/DE_CNlow_vs_Hypermutated_limma_voom.csv",
  show_col_types = FALSE
)


# =====================================
# 8. Prepare DE evidence from each comparison
# =====================================

de_1 <- de_CNhigh_CNlow %>%
  transmute(
    Gene,
    CNhigh_vs_CNlow_logFC = logFC,
    CNhigh_vs_CNlow_FDR = adj.P.Val,
    Significant_CNhigh_vs_CNlow =
      adj.P.Val < 0.05 & abs(logFC) > 1
  )

de_2 <- de_CNhigh_Hyper %>%
  transmute(
    Gene,
    CNhigh_vs_Hyper_logFC = logFC,
    CNhigh_vs_Hyper_FDR = adj.P.Val,
    Significant_CNhigh_vs_Hyper =
      adj.P.Val < 0.05 & abs(logFC) > 1
  )

de_3 <- de_CNlow_Hyper %>%
  transmute(
    Gene,
    CNlow_vs_Hyper_logFC = logFC,
    CNlow_vs_Hyper_FDR = adj.P.Val,
    Significant_CNlow_vs_Hyper =
      adj.P.Val < 0.05 & abs(logFC) > 1
  )


# =====================================
# 9. Add DE evidence to ML genes
# =====================================

combined_gene_evidence <- ml_gene_evidence %>%
  left_join(
    de_1,
    by = "Gene"
  ) %>%
  left_join(
    de_2,
    by = "Gene"
  ) %>%
  left_join(
    de_3,
    by = "Gene"
  ) %>%
  mutate(
    across(
      starts_with("Significant_"),
      ~ replace_na(.x, FALSE)
    ),
    
    Significant_DE_Comparison_Count =
      as.integer(Significant_CNhigh_vs_CNlow) +
      as.integer(Significant_CNhigh_vs_Hyper) +
      as.integer(Significant_CNlow_vs_Hyper),
    
    Supported_by_DE =
      Significant_DE_Comparison_Count > 0
  )

head(combined_gene_evidence)

# =====================================
# 10. Add gene annotations
# =====================================

# org.Hs.eg.db supplies human gene symbols, Entrez IDs
# and descriptive gene names.
#
# Some non-coding or recently named genes may not map.

annotation_table <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(combined_gene_evidence$Gene),
  keytype = "SYMBOL",
  columns = c(
    "SYMBOL",
    "ENTREZID",
    "GENENAME"
  )
)

# Check the returned column names
colnames(annotation_table)

head(annotation_table)

# Convert to tibble and explicitly use dplyr functions
# to avoid package-masking conflicts.

annotation_table <- annotation_table %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    Gene = SYMBOL,
    Entrez_ID = ENTREZID,
    Gene_Name = GENENAME
  ) %>%
  dplyr::distinct(
    Gene,
    .keep_all = TRUE
  )

head(annotation_table)

dim(annotation_table)


# Add annotation information to the combined evidence table.

biomarker_table <- combined_gene_evidence %>%
  dplyr::left_join(
    annotation_table,
    by = "Gene"
  )

head(biomarker_table)

dim(biomarker_table)


# =====================================
# 11. Calculate evidence-based priority score
# =====================================

# Points:
# +2 if selected by both LASSO and top-30 Random Forest
# +1 if selected by LASSO
# +1 if present in the top-30 Random Forest genes
# +1 for each significant differential-expression comparison
#
# This score is only used to prioritise genes for interpretation.
# It does not prove that a gene is a clinical biomarker.

biomarker_table <- biomarker_table %>%
  dplyr::mutate(
    Evidence_Score =
      2 * as.integer(Selected_by_Both) +
      as.integer(Selected_by_LASSO) +
      as.integer(Selected_by_RF_Top30) +
      Significant_DE_Comparison_Count
  ) %>%
  dplyr::arrange(
    dplyr::desc(Evidence_Score),
    dplyr::desc(Selected_by_Both),
    dplyr::desc(LASSO_Max_Absolute_Coefficient),
    RF_Rank
  )

head(
  biomarker_table,
  30
)


# =====================================
# 12. Select priority biomarker genes
# =====================================

# Include:
# - genes selected by both models, plus
# - highest-ranked remaining genes
#
# Retain a maximum of 20 genes for readable figures.

priority_genes <- biomarker_table %>%
  filter(
    Selected_by_Both |
      Evidence_Score >= 3
  ) %>%
  slice_head(
    n = 20
  )

# Ensure at least 20 genes are retained if the score threshold
# produces fewer than 20.

if (nrow(priority_genes) < 20) {
  
  priority_genes <- biomarker_table %>%
    slice_head(
      n = 20
    )
}

priority_genes

priority_gene_names <- priority_genes$Gene

length(priority_gene_names)


# =====================================
# 13. Save combined biomarker tables
# =====================================

write_csv(
  biomarker_table,
  paste0(
    "Results/Biological_Interpretation/",
    "All_ML_DE_Biomarker_Evidence.csv"
  )
)

write_csv(
  priority_genes,
  paste0(
    "Results/Biological_Interpretation/",
    "Priority_Biomarker_Genes.csv"
  )
)


# =====================================
# 14. Load original TCGA counts and metadata
# =====================================

tcga_counts <- read_csv(
  paste0(
    "TCGA Data/subtype_tumor_only/",
    "tcga_ucec_subtype_counts_geneSymbol.csv"
  ),
  show_col_types = FALSE
)

tcga_meta <- read_csv(
  paste0(
    "TCGA Data/subtype_tumor_only/",
    "tcga_ucec_subtype_metadata.csv"
  ),
  show_col_types = FALSE
)


# =====================================
# 15. Create TCGA count matrix
# =====================================

gene_names <- tcga_counts[[1]]

tcga_count_matrix <- tcga_counts[, -1] %>%
  as.data.frame() %>%
  as.matrix()

rownames(tcga_count_matrix) <- gene_names

all(
  colnames(tcga_count_matrix) ==
    tcga_meta$barcode
)


# =====================================
# 16. Remove duplicated gene symbols
# =====================================

# Retain the row with the highest total count where duplicated symbols exist.

if (any(duplicated(rownames(tcga_count_matrix)))) {
  
  tcga_count_df <- as.data.frame(
    tcga_count_matrix
  ) %>%
    rownames_to_column("Gene") %>%
    mutate(
      Total_Count = rowSums(
        across(where(is.numeric)),
        na.rm = TRUE
      )
    ) %>%
    arrange(
      Gene,
      desc(Total_Count)
    ) %>%
    distinct(
      Gene,
      .keep_all = TRUE
    ) %>%
    select(
      -Total_Count
    )
  
  tcga_count_matrix <- tcga_count_df %>%
    column_to_rownames("Gene") %>%
    as.matrix()
}


# =====================================
# 17. Check priority-gene availability
# =====================================

priority_genes_present <- intersect(
  priority_gene_names,
  rownames(tcga_count_matrix)
)

priority_genes_missing <- setdiff(
  priority_gene_names,
  rownames(tcga_count_matrix)
)

length(priority_genes_present)
priority_genes_missing


# =====================================
# 18. Generate log2 CPM expression values
# =====================================

# This transformation is used for biological visualisation.
# It is separate from model fitting.

dge_visualisation <- DGEList(
  counts = tcga_count_matrix
)

dge_visualisation <- calcNormFactors(
  dge_visualisation
)

tcga_log_cpm <- cpm(
  dge_visualisation,
  log = TRUE,
  prior.count = 1
)

dim(tcga_log_cpm)


# =====================================
# 19. Prepare priority-gene heatmap matrix
# =====================================

heatmap_matrix <- tcga_log_cpm[
  priority_genes_present,
  ,
  drop = FALSE
]

# Scale each gene across samples:
# rows = genes
# columns = samples

heatmap_scaled <- t(
  scale(
    t(heatmap_matrix)
  )
)

sum(is.na(heatmap_scaled))
sum(is.infinite(heatmap_scaled))


# =====================================
# 20. Create subtype annotation
# =====================================

heatmap_annotation <- data.frame(
  Subtype = factor(
    tcga_meta$subtype3,
    levels = c(
      "CN_high",
      "CN_low",
      "Hypermutated"
    )
  )
)

rownames(heatmap_annotation) <- tcga_meta$barcode



# Define consistent colours for the three harmonised subtypes.
annotation_colors <- list(
  Subtype = c(
    "CN_high" = "#D95F02",
    "CN_low" = "#1B9E77",
    "Hypermutated" = "#7570B3"
  )
)

# =====================================
# 21. Order samples by subtype
# =====================================

sample_order <- order(
  heatmap_annotation$Subtype
)

heatmap_scaled_ordered <- heatmap_scaled[
  ,
  sample_order,
  drop = FALSE
]

heatmap_annotation_ordered <- heatmap_annotation[
  sample_order,
  ,
  drop = FALSE
]

png(
  filename = "Figures/Biological_Interpretation/Priority_Biomarker_Heatmap.png",
  width = 12,
  height = 8,
  units = "in",
  res = 600
)

pheatmap(
  heatmap_scaled_ordered,
  
  annotation_col = heatmap_annotation_ordered,
  
  annotation_colors = annotation_colors,
  
  color = colorRampPalette(
    c("#6BAED6", "#FFF7BC", "#FD8D3C")
  )(100),
  
  show_colnames = FALSE,
  
  show_rownames = TRUE,
  
  cluster_rows = TRUE,
  
  cluster_cols = FALSE,
  
  gaps_col = c(
    sum(heatmap_annotation_ordered$Subtype == "CN_high"),
    sum(heatmap_annotation_ordered$Subtype %in% c("CN_high", "CN_low"))
  ),
  
  main = "Priority Biomarker Expression Across TCGA-UCEC Molecular Subtypes",
  
  fontsize_row = 10,
  
  border_color = NA
)

dev.off()


# =====================================
# 22. Create and save 600 dpi heatmap
# =====================================

library(grid)
library(pheatmap)

heatmap_file <- paste0(
  "Figures/Biological_Interpretation/",
  "Priority_Biomarker_Heatmap_600dpi.png"
)

# Create the heatmap object without writing directly to a file
heatmap_object <- pheatmap::pheatmap(
  heatmap_scaled_ordered,
  
  annotation_col = heatmap_annotation_ordered,
  
  annotation_colors = annotation_colors,
  
  color = colorRampPalette(
    c("#6BAED6", "#FFF7BC", "#FD8D3C")
  )(100),
  
  show_colnames = FALSE,
  
  show_rownames = TRUE,
  
  cluster_rows = TRUE,
  
  cluster_cols = FALSE,
  
  gaps_col = c(
    sum(
      heatmap_annotation_ordered$Subtype == "CN_high"
    ),
    
    sum(
      heatmap_annotation_ordered$Subtype %in%
        c("CN_high", "CN_low")
    )
  ),
  
  main =
    "Priority Biomarker Expression Across TCGA-UCEC Molecular Subtypes",
  
  fontsize_row = 10,
  
  border_color = NA,
  
  silent = TRUE
)

# Open a high-resolution PNG device
png(
  filename = heatmap_file,
  width = 12,
  height = 8,
  units = "in",
  res = 600
)

# Explicitly draw the pheatmap object
grid::grid.newpage()
grid::grid.draw(heatmap_object$gtable)

# Close the graphics device
dev.off()

# Check that the file was created properly
file.info(heatmap_file)[, c("size", "mtime")]



# =====================================
# 23. Create long-format expression dataset
# =====================================

priority_expression_long <- heatmap_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(
    cols = -Gene,
    names_to = "Barcode",
    values_to = "Log2_CPM"
  ) %>%
  dplyr::left_join(
    tcga_meta %>%
      dplyr::select(
        barcode,
        subtype3
      ),
    by = c(
      "Barcode" = "barcode"
    )
  ) %>%
  dplyr::rename(
    Subtype = subtype3
  ) %>%
  dplyr::mutate(
    Subtype = factor(
      Subtype,
      levels = c(
        "CN_high",
        "CN_low",
        "Hypermutated"
      )
    )
  )

head(priority_expression_long)






# =====================================
# 24. Calculate mean expression by subtype
# =====================================

mean_expression_by_subtype <- priority_expression_long %>%
  group_by(
    Gene,
    Subtype
  ) %>%
  summarise(
    Mean_Log2_CPM =
      mean(
        Log2_CPM,
        na.rm = TRUE
      ),
    
    SD_Log2_CPM =
      sd(
        Log2_CPM,
        na.rm = TRUE
      ),
    
    Sample_Count =
      n(),
    
    .groups = "drop"
  )

mean_expression_by_subtype


write_csv(
  mean_expression_by_subtype,
  paste0(
    "Results/Biological_Interpretation/",
    "Priority_Gene_Mean_Expression_by_Subtype.csv"
  )
)


# =====================================
# 25. Create faceted expression boxplots
# =====================================

# Display the top 12 genes to keep the figure readable.

boxplot_genes <- priority_gene_names[
  seq_len(
    min(
      12,
      length(priority_gene_names)
    )
  )
]

p_expression_boxplots <- priority_expression_long %>%
  filter(
    Gene %in% boxplot_genes
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = boxplot_genes
    )
  ) %>%
  ggplot(
    aes(
      x = Subtype,
      y = Log2_CPM,
      fill = Subtype
    )
  ) +
  geom_boxplot(
    outlier.alpha = 0.25
  ) +
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    ncol = 3
  ) +
  theme_minimal() +
  labs(
    title =
      "Expression of Priority Genes by TCGA-UCEC Molecular Subtype",
    
    x = "Molecular subtype",
    
    y = "Log2 CPM expression",
    
    fill = "Subtype"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    
    legend.position = "top"
  )

p_expression_boxplots

ggsave(
  paste0(
    "Figures/Biological_Interpretation/",
    "Priority_Gene_Expression_Boxplots.png"
  ),
  plot = p_expression_boxplots,
  width = 13,
  height = 11,
  dpi = 600
)


# =====================================
# 26. Plot evidence scores
# =====================================

p_evidence_score <- priority_genes %>%
  mutate(
    Gene = reorder(
      Gene,
      Evidence_Score
    )
  ) %>%
  ggplot(
    aes(
      x = Gene,
      y = Evidence_Score,
      fill = Selected_by_Both
    )
  ) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title =
      "Evidence Supporting Priority Candidate Biomarkers",
    
    x = "Gene",
    
    y = "Evidence score",
    
    fill = "Selected by\nboth models"
  )

p_evidence_score

ggsave(
  paste0(
    "Figures/Biological_Interpretation/",
    "Priority_Biomarker_Evidence_Scores.png"
  ),
  plot = p_evidence_score,
  width = 9,
  height = 8,
  dpi = 600
)


# =====================================
# 27. Create model-overlap summary
# =====================================

model_overlap_summary <- tibble(
  Category = c(
    "LASSO-selected genes",
    "Top 30 Random Forest genes",
    "Selected by both models"
  ),
  
  Number_of_Genes = c(
    n_distinct(
      lasso_gene_summary$Gene
    ),
    
    sum(
      rf_gene_summary$RF_Rank <= 30
    ),
    
    sum(
      ml_gene_evidence$Selected_by_Both,
      na.rm = TRUE
    )
  )
)

model_overlap_summary


write_csv(
  model_overlap_summary,
  paste0(
    "Results/Biological_Interpretation/",
    "ML_Gene_Overlap_Summary.csv"
  )
)


# =====================================
# 28. Create literature-review template
# =====================================

# R can add database descriptions, but the endometrial-cancer evidence
# must be checked manually using peer-reviewed literature.
#
# Complete these columns while reading relevant papers.

literature_review_template <- priority_genes %>%
  dplyr::select(
    Gene,
    Gene_Name,
    Entrez_ID,
    Evidence_Score,
    Selected_by_LASSO,
    Selected_by_RF_Top30,
    Selected_by_Both,
    LASSO_Associated_Subtype,
    LASSO_Coefficient,
    RF_Importance,
    Significant_DE_Comparison_Count
  ) %>%
  dplyr::mutate(
    Main_Biological_Function = NA_character_,
    
    Evidence_in_Endometrial_Cancer =
      NA_character_,
    
    Reported_Clinical_Association =
      NA_character_,
    
    Agreement_with_Current_Results =
      NA_character_,
    
    Key_Reference =
      NA_character_,
    
    DOI_or_PMID =
      NA_character_
  )

readr::write_csv(
  literature_review_template,
  paste0(
    "Results/Biological_Interpretation/",
    "Priority_Gene_Literature_Review_Template.csv"
  )
)

head(literature_review_template)

dim(literature_review_template)

# =====================================
# 29. Produce concise presentation table
# =====================================

presentation_biomarker_table <- priority_genes %>%
  slice_head(
    n = 10
  ) %>%
  transmute(
    Gene,
    
    Gene_Name,
    
    Model_Support = case_when(
      Selected_by_Both ~
        "LASSO and Random Forest",
      
      Selected_by_LASSO ~
        "LASSO",
      
      Selected_by_RF_Top30 ~
        "Random Forest",
      
      TRUE ~
        "Other"
    ),
    
    Associated_Subtype =
      LASSO_Associated_Subtype,
    
    DE_Comparisons =
      Significant_DE_Comparison_Count,
    
    Evidence_Score
  )

presentation_biomarker_table


write_csv(
  presentation_biomarker_table,
  paste0(
    "Results/Biological_Interpretation/",
    "Top10_Biomarkers_for_Presentation.csv"
  )
)


# =====================================
# 30. Final summary
# =====================================

cat(
  "\nUnique genes selected by LASSO:",
  n_distinct(lasso_gene_summary$Gene),
  "\n"
)

cat(
  "Top Random Forest genes considered:",
  sum(rf_gene_summary$RF_Rank <= 30),
  "\n"
)

cat(
  "Genes selected by both models:",
  sum(
    ml_gene_evidence$Selected_by_Both,
    na.rm = TRUE
  ),
  "\n"
)




cat(
  "Priority biomarkers retained:",
  nrow(priority_genes),
  "\n"
)

cat(
  "Priority genes present in TCGA expression matrix:",
  length(priority_genes_present),
  "\n"
)


# =====================================
# 31. Confirm saved outputs
# =====================================

list.files(
  "Results/Biological_Interpretation"
)

list.files(
  "Figures/Biological_Interpretation"
)
