# =====================================
# UCEC Project
# 15 - Final Figures and Summary Tables
# =====================================

# Aim:
# To combine the main results from the TCGA transcriptomic analysis,
# machine-learning classification, external CPTAC validation and
# biomarker prioritisation into concise figures and tables suitable
# for the presentation and final thesis.
#
# This script does not perform new model training.
# It summarises results generated in Scripts 01-14.


# =====================================
# 1. Load packages
# =====================================

library(tidyverse)
library(patchwork)


# =====================================
# 2. Create final output folders
# =====================================

dir.create(
  "Figures/Final",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Final",
  recursive = TRUE,
  showWarnings = FALSE
)


# =====================================
# 3. Load internal model metrics
# =====================================

internal_metrics <- readr::read_csv(
  "Results/ML/Model_Comparison_Primary.csv",
  show_col_types = FALSE
)

internal_metrics

# =====================================
# 4. Load external CPTAC metrics
# =====================================

external_metrics <- readr::read_csv(
  "Results/External_Validation/CPTAC_Model_Metrics.csv",
  show_col_types = FALSE
)

external_metrics


# =====================================
# 5. Combine internal and external accuracy
# =====================================

internal_accuracy <- internal_metrics %>%
  dplyr::select(
    Model,
    Accuracy
  ) %>%
  dplyr::rename(
    Internal_Accuracy = Accuracy
  )

external_accuracy <- external_metrics %>%
  dplyr::select(
    Model,
    External_Accuracy
  )

internal_external_accuracy <- internal_accuracy %>%
  dplyr::left_join(
    external_accuracy,
    by = "Model"
  )

internal_external_accuracy


# Save numerical comparison

readr::write_csv(
  internal_external_accuracy,
  "Results/Final/Internal_External_Model_Accuracy.csv"
)

# =====================================
# 6. Prepare internal/external accuracy plot
# =====================================

accuracy_long <- internal_external_accuracy %>%
  tidyr::pivot_longer(
    cols = c(
      Internal_Accuracy,
      External_Accuracy
    ),
    names_to = "Validation",
    values_to = "Accuracy"
  ) %>%
  dplyr::mutate(
    Validation = dplyr::recode(
      Validation,
      "Internal_Accuracy" = "Internal TCGA test",
      "External_Accuracy" = "External CPTAC"
    ),
    
    Accuracy_Percent = Accuracy * 100
  )

accuracy_long


p_accuracy_summary <- ggplot(
  accuracy_long,
  aes(
    x = Model,
    y = Accuracy_Percent,
    fill = Validation
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        round(
          Accuracy_Percent,
          1
        ),
        "%"
      )
    ),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.4,
    size = 4
  ) +
  theme_minimal() +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(
      0,
      100,
      by = 10
    )
  ) +
  labs(
    title = "Internal and External Classification Performance",
    subtitle = "TCGA internal test set compared with independent CPTAC validation",
    x = NULL,
    y = "Accuracy (%)",
    fill = "Validation dataset"
  ) +
  theme(
    legend.position = "top"
  )

p_accuracy_summary


ggsave(
  "Figures/Final/Internal_vs_External_Model_Accuracy.png",
  plot = p_accuracy_summary,
  width = 9,
  height = 7,
  dpi = 300
)

# =====================================
# 7. Prepare full model performance plot
# =====================================

internal_performance_long <- internal_metrics %>%
  dplyr::select(
    Model,
    Accuracy,
    Macro_Sensitivity,
    Macro_Specificity,
    Macro_F1,
    Macro_Balanced_Accuracy
  ) %>%
  tidyr::pivot_longer(
    cols = -Model,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      "Accuracy" = "Accuracy",
      "Macro_Sensitivity" = "Sensitivity",
      "Macro_Specificity" = "Specificity",
      "Macro_F1" = "F1 score",
      "Macro_Balanced_Accuracy" = "Balanced accuracy"
    ),
    
    Value_Percent = Value * 100
  )


p_internal_performance <- ggplot(
  internal_performance_long,
  aes(
    x = Metric,
    y = Value_Percent,
    fill = Model
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        round(
          Value_Percent,
          1
        ),
        "%"
      )
    ),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.4,
    size = 3.3
  ) +
  theme_minimal() +
  scale_y_continuous(
    limits = c(0, 105)
  ) +
  labs(
    title = "Held-Out TCGA Test Performance",
    x = "Performance measure",
    y = "Performance (%)",
    fill = "Model"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "top"
  )

p_internal_performance


ggsave(
  "Figures/Final/Internal_Model_Performance_Comparison.png",
  plot = p_internal_performance,
  width = 10,
  height = 7,
  dpi = 300
)


# =====================================
# 8. Load priority biomarker results
# =====================================

priority_biomarkers <- readr::read_csv(
  paste0(
    "Results/Biological_Interpretation/",
    "Priority_Biomarker_Genes.csv"
  ),
  show_col_types = FALSE
)

priority_biomarkers


# Keep the top 10 genes for presentation

top10_biomarkers <- priority_biomarkers %>%
  dplyr::slice_head(
    n = 10
  )

top10_biomarkers


# =====================================
# 9. Plot top biomarker evidence scores
# =====================================

p_top_biomarkers <- top10_biomarkers %>%
  dplyr::mutate(
    Gene = reorder(
      Gene,
      Evidence_Score
    ),
    
    Model_Support = dplyr::case_when(
      Selected_by_Both ~
        "LASSO and Random Forest",
      
      Selected_by_LASSO ~
        "LASSO only",
      
      Selected_by_RF_Top30 ~
        "Random Forest only",
      
      TRUE ~
        "Other"
    )
  ) %>%
  ggplot(
    aes(
      x = Gene,
      y = Evidence_Score,
      fill = Model_Support
    )
  ) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top Candidate Biomarkers",
    subtitle = "Prioritised using machine-learning and differential-expression evidence",
    x = "Gene",
    y = "Evidence score",
    fill = "Model support"
  )

p_top_biomarkers


ggsave(
  "Figures/Final/Top10_Candidate_Biomarkers.png",
  plot = p_top_biomarkers,
  width = 9,
  height = 7,
  dpi = 300
)

# =====================================
# 10. Create final biomarker summary table
# =====================================

final_biomarker_table <- top10_biomarkers %>%
  dplyr::transmute(
    Gene,
    
    Gene_Name,
    
    Associated_Subtype =
      LASSO_Associated_Subtype,
    
    Model_Support =
      dplyr::case_when(
        Selected_by_Both ~
          "LASSO and Random Forest",
        
        Selected_by_LASSO ~
          "LASSO",
        
        Selected_by_RF_Top30 ~
          "Random Forest",
        
        TRUE ~
          "Other"
      ),
    
    LASSO_Coefficient,
    
    RF_Importance,
    
    Significant_DE_Comparisons =
      Significant_DE_Comparison_Count,
    
    Evidence_Score
  )

final_biomarker_table


readr::write_csv(
  final_biomarker_table,
  "Results/Final/Top10_Candidate_Biomarker_Summary.csv"
)



# =====================================
# 11. Load class-specific metrics
# =====================================

lasso_class_metrics <- readr::read_csv(
  "Results/ML/LASSO_Class_Specific_Metrics.csv",
  show_col_types = FALSE
)

rf_class_metrics <- readr::read_csv(
  "Results/ML/RandomForest_Class_Specific_Metrics.csv",
  show_col_types = FALSE
)


# Retain lambda.min for LASSO

lasso_class_main <- lasso_class_metrics %>%
  dplyr::filter(
    Model == "LASSO_lambda_min"
  ) %>%
  dplyr::mutate(
    Model = "LASSO"
  )

rf_class_main <- rf_class_metrics %>%
  dplyr::mutate(
    Model = "Random Forest"
  )


class_performance <- dplyr::bind_rows(
  lasso_class_main,
  rf_class_main
) %>%
  dplyr::mutate(
    Class = stringr::str_remove(
      Class,
      "^Class: "
    )
  )

class_performance


# =====================================
# 12. Plot subtype-specific sensitivity
# =====================================

p_subtype_sensitivity <- class_performance %>%
  dplyr::mutate(
    Sensitivity_Percent =
      Sensitivity * 100
  ) %>%
  ggplot(
    aes(
      x = Class,
      y = Sensitivity_Percent,
      fill = Model
    )
  ) +
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        round(
          Sensitivity_Percent,
          1
        ),
        "%"
      )
    ),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.4,
    size = 3.6
  ) +
  theme_minimal() +
  scale_y_continuous(
    limits = c(0, 105)
  ) +
  labs(
    title = "Subtype-Specific Sensitivity",
    x = "Molecular subtype",
    y = "Sensitivity (%)",
    fill = "Model"
  ) +
  theme(
    legend.position = "top"
  )

p_subtype_sensitivity


ggsave(
  "Figures/Final/Subtype_Specific_Sensitivity.png",
  plot = p_subtype_sensitivity,
  width = 9,
  height = 7,
  dpi = 300
)

# =====================================
# 13. Load CPTAC confusion matrices
# =====================================

cptac_lasso_cm <- readr::read_csv(
  paste0(
    "Results/External_Validation/",
    "CPTAC_LASSO_ConfusionMatrix.csv"
  ),
  show_col_types = FALSE
)

cptac_rf_cm <- readr::read_csv(
  paste0(
    "Results/External_Validation/",
    "CPTAC_RandomForest_ConfusionMatrix.csv"
  ),
  show_col_types = FALSE
)

cptac_lasso_cm
cptac_rf_cm


# =====================================
# 14. Summarise CPTAC predicted classes
# =====================================

cptac_predictions <- readr::read_csv(
  paste0(
    "Results/External_Validation/",
    "CPTAC_External_Predictions.csv"
  ),
  show_col_types = FALSE
)


lasso_prediction_counts <- cptac_predictions %>%
  dplyr::count(
    LASSO_Predicted,
    name = "Count"
  ) %>%
  dplyr::rename(
    Predicted_Subtype =
      LASSO_Predicted
  ) %>%
  dplyr::mutate(
    Model = "LASSO"
  )


rf_prediction_counts <- cptac_predictions %>%
  dplyr::count(
    Random_Forest_Predicted,
    name = "Count"
  ) %>%
  dplyr::rename(
    Predicted_Subtype =
      Random_Forest_Predicted
  ) %>%
  dplyr::mutate(
    Model = "Random Forest"
  )


external_prediction_distribution <- dplyr::bind_rows(
  lasso_prediction_counts,
  rf_prediction_counts
)

external_prediction_distribution


p_external_predictions <- ggplot(
  external_prediction_distribution,
  aes(
    x = Predicted_Subtype,
    y = Count,
    fill = Model
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.8
    )
  ) +
  geom_text(
    aes(label = Count),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.4
  ) +
  theme_minimal() +
  labs(
    title = "CPTAC Predicted Subtype Distribution",
    x = "Predicted subtype",
    y = "Number of samples",
    fill = "Model"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "top"
  )

p_external_predictions


ggsave(
  "Figures/Final/CPTAC_Predicted_Subtype_Distribution.png",
  plot = p_external_predictions,
  width = 9,
  height = 7,
  dpi = 300
)

# =====================================
# 15. Combine key ML figures
# =====================================

combined_ml_summary <- (
  p_internal_performance |
    p_accuracy_summary
) / (
  p_subtype_sensitivity |
    p_top_biomarkers
) +
  patchwork::plot_annotation(
    title =
      "Summary of Machine-Learning and Biomarker Results",
    
    subtitle =
      "TCGA internal performance, CPTAC external validation and candidate biomarkers"
  )

combined_ml_summary


ggsave(
  "Figures/Final/Combined_ML_and_Biomarker_Summary.png",
  plot = combined_ml_summary,
  width = 16,
  height = 12,
  dpi = 300
)

# =====================================
# 16. Create overall project summary
# =====================================

lasso_internal_accuracy <- internal_metrics %>%
  dplyr::filter(
    Model == "LASSO"
  ) %>%
  dplyr::pull(
    Accuracy
  )

rf_internal_accuracy <- internal_metrics %>%
  dplyr::filter(
    Model == "Random Forest"
  ) %>%
  dplyr::pull(
    Accuracy
  )

lasso_external_accuracy <- external_metrics %>%
  dplyr::filter(
    Model == "LASSO"
  ) %>%
  dplyr::pull(
    External_Accuracy
  )

rf_external_accuracy <- external_metrics %>%
  dplyr::filter(
    Model == "Random Forest"
  ) %>%
  dplyr::pull(
    External_Accuracy
  )


gene_overlap_summary <- readr::read_csv(
  paste0(
    "Results/Biological_Interpretation/",
    "ML_Gene_Overlap_Summary.csv"
  ),
  show_col_types = FALSE
)


shared_gene_count <- gene_overlap_summary %>%
  dplyr::filter(
    Category == "Selected by both models"
  ) %>%
  dplyr::pull(
    Number_of_Genes
  )


final_project_summary <- tibble::tibble(
  Measure = c(
    "TCGA samples",
    "TCGA training samples",
    "TCGA held-out test samples",
    "CPTAC external samples",
    "Input genes used for modelling",
    "LASSO-selected genes",
    "Top Random Forest genes",
    "Genes supported by both models",
    "LASSO internal test accuracy",
    "Random Forest internal test accuracy",
    "LASSO CPTAC external accuracy",
    "Random Forest CPTAC external accuracy"
  ),
  
  Result = c(
    "232",
    "186",
    "46",
    "101",
    "1000",
    "56",
    "30",
    as.character(shared_gene_count),
    paste0(
      round(
        lasso_internal_accuracy * 100,
        1
      ),
      "%"
    ),
    paste0(
      round(
        rf_internal_accuracy * 100,
        1
      ),
      "%"
    ),
    paste0(
      round(
        lasso_external_accuracy * 100,
        1
      ),
      "%"
    ),
    paste0(
      round(
        rf_external_accuracy * 100,
        1
      ),
      "%"
    )
  )
)

final_project_summary


readr::write_csv(
  final_project_summary,
  "Results/Final/Final_Project_Summary.csv"
)

# =====================================
# 17. Create key findings table
# =====================================

key_findings <- tibble::tibble(
  Analysis = c(
    "Differential expression",
    "Functional enrichment",
    "Multinomial LASSO",
    "Random Forest",
    "Model comparison",
    "External CPTAC validation",
    "Biomarker prioritisation"
  ),
  
  Main_Finding = c(
    paste(
      "Thousands of genes differed between the",
      "CN_high, CN_low and Hypermutated subtypes."
    ),
    
    paste(
      "GO and KEGG results highlighted subtype-associated",
      "ciliary, microtubule, immune and signalling pathways."
    ),
    
    paste(
      "LASSO achieved 95.7% held-out TCGA accuracy",
      "and selected 56 unique genes."
    ),
    
    paste(
      "Random Forest achieved 78.3% held-out TCGA accuracy",
      "and showed reduced sensitivity for Hypermutated tumours."
    ),
    
    paste(
      "LASSO outperformed Random Forest using the same",
      "training and test partitions."
    ),
    
    paste(
      "Accuracy fell to 56.4% for LASSO and 46.5%",
      "for Random Forest in the independent CPTAC cohort."
    ),
    
    paste(
      "Eleven genes were supported by both machine-learning methods,",
      "including L1CAM, CLDN9, ACTL8 and HIF3A."
    )
  )
)

key_findings


readr::write_csv(
  key_findings,
  "Results/Final/Key_Findings_for_Presentation.csv"
)

# =====================================
# 18. Final checks
# =====================================

list.files(
  "Results/Final"
)

list.files(
  "Figures/Final"
)


cat(
  "\nFinal figures and tables successfully created.\n"
)

cat(
  "LASSO internal accuracy:",
  round(
    lasso_internal_accuracy * 100,
    1
  ),
  "%\n"
)

cat(
  "Random Forest internal accuracy:",
  round(
    rf_internal_accuracy * 100,
    1
  ),
  "%\n"
)

cat(
  "LASSO external accuracy:",
  round(
    lasso_external_accuracy * 100,
    1
  ),
  "%\n"
)

cat(
  "Random Forest external accuracy:",
  round(
    rf_external_accuracy * 100,
    1
  ),
  "%\n"
)

cat(
  "Genes supported by both models:",
  shared_gene_count,
  "\n"
)
