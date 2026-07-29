# =====================================
# UCEC Project
# 12 - LASSO vs Random Forest Comparison
# =====================================

# Aim:
# To compare the predictive performance of multinomial LASSO and
# Random Forest for classifying the three TCGA-UCEC molecular subtypes:
#
# CN_high
# CN_low
# Hypermutated
#
# Both models used:
# - The same 186 training samples
# - The same 46 held-out test samples
# - The same 1000 training-selected variable genes
# - The same five cross-validation folds
#
# This allows a fair comparison between the two algorithms.


# =====================================
# 1. Load packages
# =====================================

library(tidyverse)
library(caret)


# =====================================
# 2. Load overall model metrics
# =====================================

lasso_metrics <- read_csv(
  "Results/ML/LASSO_Overall_Metrics.csv",
  show_col_types = FALSE
)

rf_metrics <- read_csv(
  "Results/ML/RandomForest_Overall_Metrics.csv",
  show_col_types = FALSE
)


# Inspect the imported tables

lasso_metrics
rf_metrics


# =====================================
# 3. Select the main LASSO model
# =====================================

# lambda.min is used as the main LASSO model because it achieved
# the lowest cross-validated classification error.
#
# lambda.1se is retained as a secondary, more strongly regularised model.

lasso_min_metrics <- lasso_metrics %>%
  filter(
    Model == "LASSO_lambda_min"
  )

lasso_1se_metrics <- lasso_metrics %>%
  filter(
    Model == "LASSO_lambda_1se"
  )


# =====================================
# 4. Create complete comparison table
# =====================================

model_comparison_all <- bind_rows(
  lasso_min_metrics,
  lasso_1se_metrics,
  rf_metrics
) %>%
  mutate(
    Model = recode(
      Model,
      "LASSO_lambda_min" = "LASSO: lambda.min",
      "LASSO_lambda_1se" = "LASSO: lambda.1se",
      "Random_Forest" = "Random Forest"
    )
  )

model_comparison_all


# Save complete comparison table

write_csv(
  model_comparison_all,
  "Results/ML/Model_Comparison_All.csv"
)


# =====================================
# 5. Create primary model comparison
# =====================================

# The primary comparison includes:
# - LASSO lambda.min
# - Random Forest
#
# lambda.1se remains a secondary sensitivity analysis.

model_comparison_primary <- bind_rows(
  lasso_min_metrics,
  rf_metrics
) %>%
  mutate(
    Model = recode(
      Model,
      "LASSO_lambda_min" = "LASSO",
      "Random_Forest" = "Random Forest"
    )
  )

model_comparison_primary


write_csv(
  model_comparison_primary,
  "Results/ML/Model_Comparison_Primary.csv"
)


# =====================================
# 6. Convert performance metrics to percentages
# =====================================

model_comparison_percent <- model_comparison_primary %>%
  mutate(
    Accuracy_Percent = Accuracy * 100,
    Macro_Sensitivity_Percent = Macro_Sensitivity * 100,
    Macro_Specificity_Percent = Macro_Specificity * 100,
    Macro_F1_Percent = Macro_F1 * 100,
    Macro_Balanced_Accuracy_Percent =
      Macro_Balanced_Accuracy * 100
  )

model_comparison_percent


# =====================================
# 7. Prepare metrics for plotting
# =====================================

comparison_long <- model_comparison_primary %>%
  select(
    Model,
    Accuracy,
    Macro_Sensitivity,
    Macro_Specificity,
    Macro_F1,
    Macro_Balanced_Accuracy
  ) %>%
  pivot_longer(
    cols = -Model,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Accuracy" = "Accuracy",
      "Macro_Sensitivity" = "Macro sensitivity",
      "Macro_Specificity" = "Macro specificity",
      "Macro_F1" = "Macro F1",
      "Macro_Balanced_Accuracy" = "Macro balanced accuracy"
    ),
    Percentage = Value * 100
  )

comparison_long


# =====================================
# 8. Plot model performance comparison
# =====================================

dir.create(
  "Figures/ML",
  recursive = TRUE,
  showWarnings = FALSE
)

p_model_comparison <- ggplot(
  comparison_long,
  aes(
    x = Metric,
    y = Percentage,
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
        round(Percentage, 1),
        "%"
      )
    ),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.4,
    size = 3.5
  ) +
  theme_minimal() +
  labs(
    title = "Test-Set Performance: LASSO vs Random Forest",
    x = "Performance metric",
    y = "Performance (%)",
    fill = "Model"
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(
      0,
      100,
      by = 10
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "top"
  )

p_model_comparison


ggsave(
  "Figures/ML/LASSO_vs_RandomForest_Performance.png",
  plot = p_model_comparison,
  width = 10,
  height = 7,
  dpi = 300
)


# =====================================
# 9. Create simpler accuracy comparison
# =====================================

p_accuracy <- model_comparison_primary %>%
  mutate(
    Accuracy_Percent = Accuracy * 100
  ) %>%
  ggplot(
    aes(
      x = Model,
      y = Accuracy_Percent,
      fill = Model
    )
  ) +
  geom_col(
    width = 0.6
  ) +
  geom_text(
    aes(
      label = paste0(
        round(Accuracy_Percent, 1),
        "%"
      )
    ),
    vjust = -0.5,
    size = 5
  ) +
  theme_minimal() +
  labs(
    title = "Held-Out Test Accuracy by Model",
    x = NULL,
    y = "Accuracy (%)"
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(
      0,
      100,
      by = 10
    )
  ) +
  theme(
    legend.position = "none"
  )

p_accuracy


ggsave(
  "Figures/ML/LASSO_vs_RandomForest_Accuracy.png",
  plot = p_accuracy,
  width = 8,
  height = 6,
  dpi = 300
)


# =====================================
# 10. Load class-specific metrics
# =====================================

lasso_class_metrics <- read_csv(
  "Results/ML/LASSO_Class_Specific_Metrics.csv",
  show_col_types = FALSE
)

rf_class_metrics <- read_csv(
  "Results/ML/RandomForest_Class_Specific_Metrics.csv",
  show_col_types = FALSE
)


# Keep lambda.min only for the main LASSO comparison

lasso_min_class_metrics <- lasso_class_metrics %>%
  filter(
    Model == "LASSO_lambda_min"
  ) %>%
  mutate(
    Model = "LASSO"
  )

rf_class_metrics <- rf_class_metrics %>%
  mutate(
    Model = "Random Forest"
  )


class_comparison <- bind_rows(
  lasso_min_class_metrics,
  rf_class_metrics
) %>%
  mutate(
    Class = str_remove(
      Class,
      "^Class: "
    )
  )

class_comparison


write_csv(
  class_comparison,
  "Results/ML/Model_Comparison_Class_Specific.csv"
)


# =====================================
# 11. Plot class-specific sensitivity
# =====================================

p_class_sensitivity <- class_comparison %>%
  mutate(
    Sensitivity_Percent = Sensitivity * 100
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
        round(Sensitivity_Percent, 1),
        "%"
      )
    ),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.4,
    size = 3.5
  ) +
  theme_minimal() +
  labs(
    title = "Subtype-Specific Sensitivity by Model",
    x = "Molecular subtype",
    y = "Sensitivity (%)",
    fill = "Model"
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(
      0,
      100,
      by = 10
    )
  ) +
  theme(
    legend.position = "top"
  )

p_class_sensitivity


ggsave(
  "Figures/ML/LASSO_vs_RandomForest_Class_Sensitivity.png",
  plot = p_class_sensitivity,
  width = 9,
  height = 7,
  dpi = 300
)


# =====================================
# 12. Plot class-specific F1 scores
# =====================================

p_class_f1 <- class_comparison %>%
  mutate(
    F1_Percent = F1 * 100
  ) %>%
  ggplot(
    aes(
      x = Class,
      y = F1_Percent,
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
        round(F1_Percent, 1),
        "%"
      )
    ),
    position = position_dodge(
      width = 0.8
    ),
    vjust = -0.4,
    size = 3.5
  ) +
  theme_minimal() +
  labs(
    title = "Subtype-Specific F1 Score by Model",
    x = "Molecular subtype",
    y = "F1 score (%)",
    fill = "Model"
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(
      0,
      100,
      by = 10
    )
  ) +
  theme(
    legend.position = "top"
  )

p_class_f1


ggsave(
  "Figures/ML/LASSO_vs_RandomForest_Class_F1.png",
  plot = p_class_f1,
  width = 9,
  height = 7,
  dpi = 300
)


# =====================================
# 13. Load test-set predictions
# =====================================

lasso_predictions <- read_csv(
  "Results/ML/LASSO_Test_Predictions.csv",
  show_col_types = FALSE
)

rf_predictions <- read_csv(
  "Results/ML/RandomForest_Test_Predictions.csv",
  show_col_types = FALSE
)


# Keep the primary LASSO lambda.min prediction.

lasso_predictions_main <- lasso_predictions %>%
  select(
    Barcode,
    Actual,
    LASSO_Predicted = Predicted_lambda_min,
    LASSO_Correct = Correct_lambda_min
  )


rf_predictions_main <- rf_predictions %>%
  select(
    Barcode,
    Actual,
    RF_Predicted = Predicted,
    RF_Correct = Correct
  )


# Join predictions using sample barcode.

paired_predictions <- lasso_predictions_main %>%
  inner_join(
    rf_predictions_main,
    by = "Barcode",
    suffix = c(
      "_LASSO",
      "_RF"
    )
  )


# Confirm actual labels match in both files.

all(
  paired_predictions$Actual_LASSO ==
    paired_predictions$Actual_RF
)


paired_predictions <- paired_predictions %>%
  mutate(
    Actual = Actual_LASSO,
    Both_Correct =
      LASSO_Correct & RF_Correct,
    LASSO_Only_Correct =
      LASSO_Correct & !RF_Correct,
    RF_Only_Correct =
      !LASSO_Correct & RF_Correct,
    Both_Incorrect =
      !LASSO_Correct & !RF_Correct
  ) %>%
  select(
    Barcode,
    Actual,
    LASSO_Predicted,
    RF_Predicted,
    LASSO_Correct,
    RF_Correct,
    Both_Correct,
    LASSO_Only_Correct,
    RF_Only_Correct,
    Both_Incorrect
  )

paired_predictions


write_csv(
  paired_predictions,
  "Results/ML/LASSO_RandomForest_Paired_Predictions.csv"
)


# =====================================
# 14. Summarise agreement between models
# =====================================

prediction_agreement_summary <- tibble(
  Category = c(
    "Both correct",
    "LASSO only correct",
    "Random Forest only correct",
    "Both incorrect"
  ),
  Count = c(
    sum(
      paired_predictions$Both_Correct
    ),
    sum(
      paired_predictions$LASSO_Only_Correct
    ),
    sum(
      paired_predictions$RF_Only_Correct
    ),
    sum(
      paired_predictions$Both_Incorrect
    )
  )
)

prediction_agreement_summary


write_csv(
  prediction_agreement_summary,
  "Results/ML/Model_Prediction_Agreement_Summary.csv"
)


# =====================================
# 15. Paired McNemar test
# =====================================

# Both models predicted the same test samples.
# McNemar's test compares their paired correct/incorrect outcomes.
#
# Interpretation should be cautious because the test set has only 46 samples
# and there may be relatively few discordant predictions.

paired_correctness_table <- table(
  LASSO_Correct = paired_predictions$LASSO_Correct,
  RF_Correct = paired_predictions$RF_Correct
)

paired_correctness_table


mcnemar_result <- mcnemar.test(
  paired_correctness_table,
  correct = TRUE
)

mcnemar_result


mcnemar_summary <- tibble(
  Statistic = unname(
    mcnemar_result$statistic
  ),
  Degrees_of_Freedom = unname(
    mcnemar_result$parameter
  ),
  P_Value = mcnemar_result$p.value,
  Method = mcnemar_result$method
)

mcnemar_summary


write_csv(
  mcnemar_summary,
  "Results/ML/LASSO_vs_RandomForest_McNemar_Test.csv"
)


# =====================================
# 16. Load shared important genes
# =====================================

shared_genes <- read_csv(
  "Results/ML/LASSO_RandomForest_Shared_Genes.csv",
  show_col_types = FALSE
)

shared_genes

nrow(shared_genes)


# =====================================
# 17. Create main numerical summary
# =====================================

lasso_accuracy <- model_comparison_primary %>%
  filter(
    Model == "LASSO"
  ) %>%
  pull(
    Accuracy
  )

rf_accuracy <- model_comparison_primary %>%
  filter(
    Model == "Random Forest"
  ) %>%
  pull(
    Accuracy
  )

accuracy_difference <- lasso_accuracy - rf_accuracy


model_summary <- tibble(
  LASSO_Test_Accuracy = lasso_accuracy,
  Random_Forest_Test_Accuracy = rf_accuracy,
  Accuracy_Difference = accuracy_difference,
  LASSO_Correct_Samples = sum(
    paired_predictions$LASSO_Correct
  ),
  Random_Forest_Correct_Samples = sum(
    paired_predictions$RF_Correct
  ),
  Test_Sample_Size = nrow(
    paired_predictions
  ),
  Shared_Important_Genes = nrow(
    shared_genes
  ),
  McNemar_P_Value = mcnemar_result$p.value
)

model_summary


write_csv(
  model_summary,
  "Results/ML/Final_Model_Comparison_Summary.csv"
)


# =====================================
# 18. Print final results
# =====================================

cat(
  "\nLASSO test accuracy:",
  round(
    lasso_accuracy * 100,
    1
  ),
  "%\n"
)

cat(
  "Random Forest test accuracy:",
  round(
    rf_accuracy * 100,
    1
  ),
  "%\n"
)

cat(
  "Difference in test accuracy:",
  round(
    accuracy_difference * 100,
    1
  ),
  "percentage points\n"
)

cat(
  "LASSO correctly classified:",
  sum(
    paired_predictions$LASSO_Correct
  ),
  "of",
  nrow(
    paired_predictions
  ),
  "samples\n"
)

cat(
  "Random Forest correctly classified:",
  sum(
    paired_predictions$RF_Correct
  ),
  "of",
  nrow(
    paired_predictions
  ),
  "samples\n"
)

cat(
  "Shared important genes:",
  nrow(
    shared_genes
  ),
  "\n"
)

cat(
  "McNemar test p-value:",
  round(
    mcnemar_result$p.value,
    4
  ),
  "\n"
)


# =====================================
# 19. Confirm saved outputs
# =====================================

list.files(
  "Results/ML"
)

list.files(
  "Figures/ML"
  
  
  
  
  #Results 
  #Using the same training samples, test samples, selected genes 
  #and cross-validation folds, multinomial LASSO achieved higher
  #held-out performance than Random Forest. LASSO correctly classified 
  #44 of 46 samples, compared with 36 of 46 for Random Forest. 
  #The largest Random Forest limitation was reduced sensitivity for 
  #the Hypermutated subtype. This may indicate that subtype-associated 
  #expression differences are captured more effectively by a sparse linear 
  #combination of genes than by the nonlinear tree-based model in this 
  #relatively small, high-dimensional cohort.
  
  #Important caution
  
  #The McNemar test is useful because both models predicted exactly the same 46 samples. However, the test set is small, so the p-value should be treated as supportive rather than definitive. The performance comparison should focus on:
    
    #held-out accuracy;
 # macro F1;
  #subtype-specific sensitivity;
  #confusion matrices;
  #external validation on CPTAC.
  
  
  
  
  
  
  
)
