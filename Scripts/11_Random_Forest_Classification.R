# =====================================
# UCEC Project
# 11 - Random Forest Classification
# =====================================

# Aim:
# To train and evaluate a multiclass Random Forest classifier for predicting
# the three TCGA-UCEC molecular subtypes:
#
# CN_high
# CN_low
# Hypermutated
#
# The model uses:
# - The same training and test sets used for LASSO
# - The same top 1000 training-selected variable genes
# - The same five cross-validation folds used for LASSO
# - Hyperparameter tuning using training data only
# - The independent test set for final evaluation
#
# This allows a fair comparison between LASSO and Random Forest.


# =====================================
# 1. Install package once if necessary
# =====================================

# Only run this if library(ranger) gives an installation error:
#




# =====================================
# 2. Load packages
# =====================================

library(ranger)
library(caret)
library(tidyverse)


dir.create(
  "Results/ML",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/ML",
  recursive = TRUE,
  showWarnings = FALSE
)

# =====================================
# 3. Set seed for reproducibility
# =====================================

set.seed(123)


# =====================================
# 4. Load prepared machine-learning data
# =====================================

x_train <- readRDS(
  "Results/ML/x_train_scaled.rds"
)

x_test <- readRDS(
  "Results/ML/x_test_scaled.rds"
)

y_train <- readRDS(
  "Results/ML/y_train.rds"
)

y_test <- readRDS(
  "Results/ML/y_test.rds"
)


# Load the same cross-validation fold assignments used for LASSO.

foldid <- readRDS(
  "Results/ML/LASSO_5fold_foldid.rds"
)


# =====================================
# 5. Check loaded objects
# =====================================

dim(x_train)
dim(x_test)

length(y_train)
length(y_test)

table(y_train)
table(y_test)

table(foldid)

identical(
  colnames(x_train),
  colnames(x_test)
)

sum(is.na(x_train))
sum(is.na(x_test))

sum(is.infinite(x_train))
sum(is.infinite(x_test))


# Expected:
#
# x_train = 186 samples x 1000 genes
# x_test  = 46 samples x 1000 genes
#
# identical columns = TRUE
# missing and infinite values = 0


# =====================================
# 6. Ensure correct formats
# =====================================

x_train <- as.data.frame(x_train)
x_test <- as.data.frame(x_test)

y_train <- factor(
  y_train,
  levels = c(
    "CN_high",
    "CN_low",
    "Hypermutated"
  )
)

y_test <- factor(
  y_test,
  levels = levels(y_train)
)


# Add the subtype outcome to the training dataset.

rf_train_data <- x_train %>%
  mutate(
    Subtype = y_train,
    .before = 1
  )

dim(rf_train_data)

head(rf_train_data[, 1:6])


# =====================================
# 7. Recreate LASSO cross-validation folds
# =====================================

# caret requires the training indices for each resampling fold.
#
# foldid identifies the validation fold for each sample.
# Therefore, for each fold, the training indices are all samples not
# belonging to that validation fold.

cv_training_indices <- lapply(
  sort(unique(foldid)),
  function(current_fold) {
    which(foldid != current_fold)
  }
)

names(cv_training_indices) <- paste0(
  "Fold",
  sort(unique(foldid))
)


# The corresponding validation samples can also be saved explicitly.

cv_validation_indices <- lapply(
  sort(unique(foldid)),
  function(current_fold) {
    which(foldid == current_fold)
  }
)

names(cv_validation_indices) <- paste0(
  "Fold",
  sort(unique(foldid))
)


# Check fold sizes.

sapply(
  cv_training_indices,
  length
)

sapply(
  cv_validation_indices,
  length
)


# Check subtype distribution in each validation fold.

for (i in seq_along(cv_validation_indices)) {
  
  cat(
    "\nValidation",
    names(cv_validation_indices)[i],
    "\n"
  )
  
  print(
    table(
      y_train[
        cv_validation_indices[[i]]
      ]
    )
  )
}


# =====================================
# 8. Define five-fold cross-validation
# =====================================

rf_control <- trainControl(
  method = "cv",
  number = 5,
  index = cv_training_indices,
  indexOut = cv_validation_indices,
  classProbs = TRUE,
  savePredictions = "final",
  allowParallel = TRUE,
  verboseIter = TRUE
)


# =====================================
# 9. Create Random Forest tuning grid
# =====================================

# mtry:
# Number of genes randomly considered at each tree split.
#
# splitrule:
# Rule used to determine the best split.
#
# min.node.size:
# Minimum number of observations allowed in a terminal node.
#
# With 1000 predictors, sqrt(1000) is approximately 32.
# Values below and above this are also tested.

rf_tune_grid <- expand.grid(
  mtry = c(
    10,
    32,
    100
  ),
  splitrule = "gini",
  min.node.size = c(
    1,
    5,
    10
  )
)

rf_tune_grid


# This creates nine candidate parameter combinations.

# =====================================
# 10. Train and tune Random Forest
# =====================================

# metric = "Accuracy" selects the model with the highest
# mean cross-validated accuracy.
#
# num.trees = 500 gives each candidate model 500 trees.
#
# importance = "permutation" allows genes to be ranked
# according to their contribution to prediction.
#
# Do not add probability = TRUE here because classProbs = TRUE
# in trainControl already requests probability predictions.

set.seed(123)

rf_model <- train(
  Subtype ~ .,
  data = rf_train_data,
  method = "ranger",
  metric = "Accuracy",
  trControl = rf_control,
  tuneGrid = rf_tune_grid,
  num.trees = 500,
  importance = "permutation"
)


# =====================================
# 11. Inspect tuned model
# =====================================

rf_model

rf_model$bestTune

rf_model$results


# Display cross-validation results ordered by accuracy.

rf_cv_results <- rf_model$results %>%
  arrange(
    desc(Accuracy)
  )

rf_cv_results


# =====================================
# 12. Plot tuning results
# =====================================

dir.create(
  "Figures/ML",
  recursive = TRUE,
  showWarnings = FALSE
)

p_rf_tuning <- ggplot(
  rf_cv_results,
  aes(
    x = factor(mtry),
    y = Accuracy,
    group = factor(min.node.size),
    linetype = factor(min.node.size)
  )
) +
  geom_point(
    size = 3
  ) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = Accuracy - AccuracySD,
      ymax = Accuracy + AccuracySD
    ),
    width = 0.1
  ) +
  theme_minimal() +
  labs(
    title = "Five-Fold Cross-Validation for Random Forest",
    x = "Number of genes considered at each split (mtry)",
    y = "Cross-validated accuracy",
    linetype = "Minimum node size"
  )

p_rf_tuning

ggsave(
  "Figures/ML/RandomForest_CrossValidation.png",
  plot = p_rf_tuning,
  width = 9,
  height = 7,
  dpi = 300
)


# =====================================
# 13. Predict held-out test samples
# =====================================

rf_test_predictions <- predict(
  rf_model,
  newdata = x_test,
  type = "raw"
)

rf_test_predictions <- factor(
  rf_test_predictions,
  levels = levels(y_test)
)

table(rf_test_predictions)


# =====================================
# 14. Test-set confusion matrix
# =====================================

rf_confusion <- confusionMatrix(
  data = rf_test_predictions,
  reference = y_test
)

rf_confusion


# =====================================
# 15. Extract overall performance metrics
# =====================================

rf_overall <- rf_confusion$overall

rf_class_metrics_matrix <- as.matrix(
  rf_confusion$byClass
)

rf_metrics <- tibble(
  Model = "Random_Forest",
  Accuracy = unname(
    rf_overall["Accuracy"]
  ),
  Kappa = unname(
    rf_overall["Kappa"]
  ),
  Macro_Sensitivity = mean(
    rf_class_metrics_matrix[, "Sensitivity"],
    na.rm = TRUE
  ),
  Macro_Specificity = mean(
    rf_class_metrics_matrix[, "Specificity"],
    na.rm = TRUE
  ),
  Macro_F1 = mean(
    rf_class_metrics_matrix[, "F1"],
    na.rm = TRUE
  ),
  Macro_Balanced_Accuracy = mean(
    rf_class_metrics_matrix[, "Balanced Accuracy"],
    na.rm = TRUE
  )
)

rf_metrics


# =====================================
# 16. Extract class-specific metrics
# =====================================

rf_class_metrics <- as.data.frame(
  rf_confusion$byClass
) %>%
  rownames_to_column(
    "Class"
  ) %>%
  mutate(
    Model = "Random_Forest",
    .before = 1
  )

rf_class_metrics


# =====================================
# 17. Create confusion-matrix table
# =====================================

rf_confusion_df <- as.data.frame(
  rf_confusion$table
) %>%
  rename(
    Predicted = Prediction,
    Actual = Reference,
    Count = Freq
  )

rf_confusion_df


# =====================================
# 18. Plot confusion matrix
# =====================================

p_rf_confusion <- ggplot(
  rf_confusion_df,
  aes(
    x = Actual,
    y = Predicted,
    fill = Count
  )
) +
  geom_tile() +
  geom_text(
    aes(label = Count),
    size = 6
  ) +
  theme_minimal() +
  labs(
    title = "Random Forest Test Confusion Matrix",
    x = "Actual subtype",
    y = "Predicted subtype"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

p_rf_confusion

ggsave(
  "Figures/ML/RandomForest_ConfusionMatrix.png",
  plot = p_rf_confusion,
  width = 8,
  height = 6,
  dpi = 300
)


# =====================================
# 19. Obtain test-set class probabilities
# =====================================

rf_test_probabilities <- predict(
  rf_model,
  newdata = x_test,
  type = "prob"
)

head(rf_test_probabilities)


# Create probability results table.

rf_probability_results <- tibble(
  Barcode = rownames(x_test),
  Actual = as.character(y_test),
  Predicted = as.character(
    rf_test_predictions
  )
) %>%
  bind_cols(
    as.data.frame(
      rf_test_probabilities
    )
  )

head(rf_probability_results)


# =====================================
# 20. Create test-prediction table
# =====================================

rf_test_results <- tibble(
  Barcode = rownames(x_test),
  Actual = as.character(y_test),
  Predicted = as.character(
    rf_test_predictions
  ),
  Correct = Actual == Predicted
)

head(rf_test_results)

table(rf_test_results$Correct)


# =====================================
# 21. Extract permutation importance
# =====================================

# The final ranger model is stored inside finalModel.
#
# Permutation importance measures the reduction in predictive performance
# when the values of a gene are randomly disrupted.

rf_importance_vector <- rf_model$finalModel$variable.importance

length(rf_importance_vector)

head(
  sort(
    rf_importance_vector,
    decreasing = TRUE
  )
)


# Convert importance values to a table.

rf_importance <- tibble(
  Gene = names(
    rf_importance_vector
  ),
  Importance = as.numeric(
    rf_importance_vector
  )
) %>%
  arrange(
    desc(Importance)
  )

head(rf_importance, 20)


# =====================================
# 22. Select top 30 important genes
# =====================================

top_rf_genes <- rf_importance %>%
  slice_head(
    n = 30
  )

top_rf_genes


# =====================================
# 23. Plot top Random Forest genes
# =====================================

p_rf_importance <- ggplot(
  top_rf_genes,
  aes(
    x = reorder(
      Gene,
      Importance
    ),
    y = Importance
  )
) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top Random Forest Genes",
    x = "Gene",
    y = "Permutation importance"
  )

p_rf_importance

ggsave(
  "Figures/ML/RandomForest_Top_Genes.png",
  plot = p_rf_importance,
  width = 9,
  height = 8,
  dpi = 300
)


# =====================================
# 24. Compare Random Forest genes with LASSO
# =====================================

lasso_selected_genes <- read_csv(
  "Results/ML/LASSO_Selected_Genes_lambda_min.csv",
  show_col_types = FALSE
)

lasso_gene_names <- unique(
  lasso_selected_genes$Gene
)

rf_top_gene_names <- top_rf_genes$Gene


# Find genes selected by both methods.

shared_lasso_rf_genes <- intersect(
  lasso_gene_names,
  rf_top_gene_names
)

shared_lasso_rf_genes

length(shared_lasso_rf_genes)


# Save the overlap as a table.

shared_gene_table <- tibble(
  Gene = shared_lasso_rf_genes
)

shared_gene_table


# =====================================
# 25. Save Random Forest outputs
# =====================================

write_csv(
  rf_cv_results,
  "Results/ML/RandomForest_CV_Results.csv"
)

write_csv(
  as_tibble(rf_model$bestTune),
  "Results/ML/RandomForest_Best_Tuning_Parameters.csv"
)

write_csv(
  rf_metrics,
  "Results/ML/RandomForest_Overall_Metrics.csv"
)

write_csv(
  rf_class_metrics,
  "Results/ML/RandomForest_Class_Specific_Metrics.csv"
)

write_csv(
  rf_confusion_df,
  "Results/ML/RandomForest_ConfusionMatrix.csv"
)

write_csv(
  rf_test_results,
  "Results/ML/RandomForest_Test_Predictions.csv"
)

write_csv(
  rf_probability_results,
  "Results/ML/RandomForest_Test_Probabilities.csv"
)

write_csv(
  rf_importance,
  "Results/ML/RandomForest_All_Gene_Importance.csv"
)

write_csv(
  top_rf_genes,
  "Results/ML/RandomForest_Top30_Genes.csv"
)

write_csv(
  shared_gene_table,
  "Results/ML/LASSO_RandomForest_Shared_Genes.csv"
)


# =====================================
# 26. Save fitted Random Forest model
# =====================================

saveRDS(
  rf_model,
  "Results/ML/RandomForest_Model.rds"
)


# =====================================
# 27. Final model summary
# =====================================

cat(
  "\nBest Random Forest tuning parameters:\n"
)

print(
  rf_model$bestTune
)

cat(
  "\nRandom Forest test accuracy:",
  round(
    unname(
      rf_confusion$overall["Accuracy"]
    ),
    3
  ),
  "\n"
)

cat(
  "Random Forest test Kappa:",
  round(
    unname(
      rf_confusion$overall["Kappa"]
    ),
    3
  ),
  "\n"
)

cat(
  "Number of shared genes between LASSO and top 30 Random Forest genes:",
  length(shared_lasso_rf_genes),
  "\n"
)


# =====================================
# 28. Confirm saved outputs
# =====================================

list.files("Results/ML")

list.files("Figures/ML")
