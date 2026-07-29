# =====================================
# UCEC Project
# 10 - Multinomial LASSO Classification
# =====================================

# Aim:
# To train and evaluate a multinomial LASSO classifier that predicts
# three TCGA-UCEC molecular subtypes:
# CN_high, CN_low and Hypermutated.
#
# The model uses:
# - The training and test datasets prepared in Script 09
# - Five-fold stratified cross-validation
# - LASSO penalty: alpha = 1
# - A completely held-out test set for final evaluation
#
# Two lambda values are examined:
# lambda.min  = lambda with the lowest cross-validation error
# lambda.1se  = a more regularised and usually simpler model within
#               one standard error of the minimum


# =====================================
# 1. Load packages
# =====================================

library(glmnet)
library(caret)
library(tidyverse)
library(Matrix)


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
# 2. Set seed for reproducibility
# =====================================

set.seed(123)


# =====================================
# 3. Load prepared ML datasets
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


# =====================================
# 4. Check loaded objects
# =====================================

dim(x_train)
dim(x_test)

length(y_train)
length(y_test)

table(y_train)
table(y_test)

identical(
  colnames(x_train),
  colnames(x_test)
)

sum(is.na(x_train))
sum(is.na(x_test))

sum(is.infinite(x_train))
sum(is.infinite(x_test))


# Expected:
# x_train = 186 samples x 1000 genes
# x_test  = 46 samples x 1000 genes
#
# identical columns = TRUE
# missing/infinite values = 0


# =====================================
# 5. Ensure correct matrix and factor formats
# =====================================

x_train <- as.matrix(x_train)
x_test <- as.matrix(x_test)

y_train <- factor(
  y_train,
  levels = c("CN_high", "CN_low", "Hypermutated")
)

y_test <- factor(
  y_test,
  levels = levels(y_train)
)


# =====================================
# 6. Create stratified five-fold CV IDs
# =====================================

# createFolds aims to preserve the class distribution within each fold.
# foldid is supplied directly to cv.glmnet so the same folds are used
# throughout model tuning.

foldid <- createFolds(
  y_train,
  k = 5,
  list = FALSE,
  returnTrain = FALSE
)

table(foldid)

# Check subtype balance within each fold
table(
  Fold = foldid,
  Subtype = y_train
)


# =====================================
# 7. Fit multinomial LASSO using CV
# =====================================

# alpha = 1 specifies pure LASSO regularisation.
#
# family = "multinomial" is required because there are three classes.
#
# type.measure = "class" selects lambda using classification error.
#
# standardize = FALSE because Script 09 already standardised the
# predictors using training-set means and standard deviations.
#
# keep = TRUE retains cross-validation predictions if needed later.

cv_lasso <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "multinomial",
  alpha = 1,
  foldid = foldid,
  type.measure = "class",
  standardize = FALSE,
  keep = TRUE
)


# =====================================
# 8. Inspect selected lambda values
# =====================================

cv_lasso$lambda.min
cv_lasso$lambda.1se

# lambda.min:
# Lowest cross-validated classification error.
#
# lambda.1se:
# More strongly regularised model within one standard error of the minimum.

# =====================================
# 9. Save cross-validation plot
# =====================================

dir.create(
  "Figures/ML",
  recursive = TRUE,
  showWarnings = FALSE
)

png(
  filename = "Figures/ML/LASSO_CrossValidation.png",
  width = 2400,
  height = 1800,
  res = 300
)

par(
  mar = c(5, 5, 7, 2)
)

plot(
  cv_lasso,
  main = ""
)

title(
  main = "Five-Fold Cross-Validation for Multinomial LASSO",
  line = 2,
  cex.main = 1.2,
  font.main = 2
)

dev.off()
# =====================================
# 10. Generate test-set predictions
# lambda.min
# =====================================

pred_min_raw <- predict(
  cv_lasso,
  newx = x_test,
  s = "lambda.min",
  type = "class"
)

# Convert prediction output into a factor
pred_min <- factor(
  as.vector(pred_min_raw),
  levels = levels(y_test)
)

table(pred_min)


# =====================================
# 11. Test-set confusion matrix
# lambda.min
# =====================================

cm_min <- confusionMatrix(
  data = pred_min,
  reference = y_test
)

cm_min


# =====================================
# 12. Generate test predictions
# lambda.1se
# =====================================

pred_1se_raw <- predict(
  cv_lasso,
  newx = x_test,
  s = "lambda.1se",
  type = "class"
)

pred_1se <- factor(
  as.vector(pred_1se_raw),
  levels = levels(y_test)
)

cm_1se <- confusionMatrix(
  data = pred_1se,
  reference = y_test
)

cm_1se


# =====================================
# 13. Extract overall performance metrics
# =====================================

# Helper function to calculate macro-average metrics.
# For multiclass classification, caret returns one row per class.

extract_lasso_metrics <- function(
    confusion_object,
    model_name
) {
  
  overall <- confusion_object$overall
  class_metrics <- confusion_object$byClass
  
  # Ensure byClass is treated as a matrix
  class_metrics <- as.matrix(class_metrics)
  
  tibble(
    Model = model_name,
    Accuracy = unname(overall["Accuracy"]),
    Kappa = unname(overall["Kappa"]),
    Macro_Sensitivity = mean(
      class_metrics[, "Sensitivity"],
      na.rm = TRUE
    ),
    Macro_Specificity = mean(
      class_metrics[, "Specificity"],
      na.rm = TRUE
    ),
    Macro_F1 = mean(
      class_metrics[, "F1"],
      na.rm = TRUE
    ),
    Macro_Balanced_Accuracy = mean(
      class_metrics[, "Balanced Accuracy"],
      na.rm = TRUE
    )
  )
}

metrics_min <- extract_lasso_metrics(
  cm_min,
  "LASSO_lambda_min"
)

metrics_1se <- extract_lasso_metrics(
  cm_1se,
  "LASSO_lambda_1se"
)

lasso_metrics <- bind_rows(
  metrics_min,
  metrics_1se
)

lasso_metrics


# =====================================
# 14. Extract class-specific metrics
# =====================================

class_metrics_min <- as.data.frame(
  cm_min$byClass
) %>%
  rownames_to_column("Class") %>%
  mutate(
    Model = "LASSO_lambda_min",
    .before = 1
  )

class_metrics_1se <- as.data.frame(
  cm_1se$byClass
) %>%
  rownames_to_column("Class") %>%
  mutate(
    Model = "LASSO_lambda_1se",
    .before = 1
  )

lasso_class_metrics <- bind_rows(
  class_metrics_min,
  class_metrics_1se
)

lasso_class_metrics


# =====================================
# 15. Create confusion-matrix data frames
# =====================================

confusion_min_df <- as.data.frame(
  cm_min$table
) %>%
  rename(
    Predicted = Prediction,
    Actual = Reference,
    Count = Freq
  )

confusion_1se_df <- as.data.frame(
  cm_1se$table
) %>%
  rename(
    Predicted = Prediction,
    Actual = Reference,
    Count = Freq
  )

confusion_min_df
confusion_1se_df


# =====================================
# 16. Plot confusion matrix
# lambda.min
# =====================================

p_cm_min <- ggplot(
  confusion_min_df,
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
    title = "LASSO Test Confusion Matrix: lambda.min",
    x = "Actual subtype",
    y = "Predicted subtype"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

p_cm_min

ggsave(
  "Figures/ML/LASSO_ConfusionMatrix_lambda_min.png",
  plot = p_cm_min,
  width = 8,
  height = 6,
  dpi = 300
)


# =====================================
# 17. Plot confusion matrix
# lambda.1se
# =====================================

p_cm_1se <- ggplot(
  confusion_1se_df,
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
    title = "LASSO Test Confusion Matrix: lambda.1se",
    x = "Actual subtype",
    y = "Predicted subtype"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

p_cm_1se

ggsave(
  "Figures/ML/LASSO_ConfusionMatrix_lambda_1se.png",
  plot = p_cm_1se,
  width = 8,
  height = 6,
  dpi = 300
)


# =====================================
# 18. Extract predicted probabilities
# =====================================

prob_min_raw <- predict(
  cv_lasso,
  newx = x_test,
  s = "lambda.min",
  type = "response"
)

# Multinomial predictions may be returned as a three-dimensional array.
# Drop the final singleton dimension.

prob_min <- drop(prob_min_raw)

prob_min <- as.data.frame(prob_min)

colnames(prob_min) <- levels(y_train)

probability_results <- tibble(
  Barcode = rownames(x_test),
  Actual = as.character(y_test),
  Predicted = as.character(pred_min)
) %>%
  bind_cols(prob_min)

head(probability_results)


# =====================================
# 19. Extract non-zero LASSO coefficients
# =====================================

# For multinomial glmnet models, coef() returns one coefficient matrix
# for each subtype.

extract_selected_genes <- function(
    fitted_model,
    lambda_value,
    lambda_name
) {
  
  coefficient_list <- coef(
    fitted_model,
    s = lambda_value
  )
  
  selected_tables <- lapply(
    names(coefficient_list),
    function(class_name) {
      
      coefficient_matrix <- coefficient_list[[class_name]]
      
      coefficient_df <- tibble(
        Gene = rownames(coefficient_matrix),
        Coefficient = as.numeric(coefficient_matrix)
      ) %>%
        filter(
          Gene != "(Intercept)",
          Coefficient != 0
        ) %>%
        mutate(
          Class = class_name,
          Lambda = lambda_name,
          Absolute_Coefficient = abs(Coefficient)
        ) %>%
        arrange(
          desc(Absolute_Coefficient)
        )
      
      coefficient_df
    }
  )
  
  bind_rows(selected_tables)
}

selected_genes_min <- extract_selected_genes(
  cv_lasso,
  cv_lasso$lambda.min,
  "lambda.min"
)

selected_genes_1se <- extract_selected_genes(
  cv_lasso,
  cv_lasso$lambda.1se,
  "lambda.1se"
)

head(selected_genes_min)
head(selected_genes_1se)

nrow(selected_genes_min)
nrow(selected_genes_1se)

n_distinct(selected_genes_min$Gene)
n_distinct(selected_genes_1se$Gene)


# =====================================
# 20. Select top coefficients for plotting
# =====================================

top_lasso_genes_min <- selected_genes_min %>%
  group_by(Class) %>%
  slice_max(
    order_by = Absolute_Coefficient,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

top_lasso_genes_min


# =====================================
# 21. Plot top selected genes
# =====================================

p_lasso_genes <- ggplot(
  top_lasso_genes_min,
  aes(
    x = reorder(
      paste(Gene, Class, sep = " — "),
      Absolute_Coefficient
    ),
    y = Absolute_Coefficient,
    fill = Class
  )
) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top LASSO-Selected Genes by Molecular Subtype",
    x = "Gene and subtype",
    y = "Absolute LASSO coefficient",
    fill = "Subtype"
  )

p_lasso_genes

ggsave(
  "Figures/ML/LASSO_Top_Selected_Genes.png",
  plot = p_lasso_genes,
  width = 10,
  height = 9,
  dpi = 300
)


# =====================================
# 22. Create test prediction table
# =====================================

test_predictions <- tibble(
  Barcode = rownames(x_test),
  Actual = as.character(y_test),
  Predicted_lambda_min = as.character(pred_min),
  Correct_lambda_min = Actual == Predicted_lambda_min,
  Predicted_lambda_1se = as.character(pred_1se),
  Correct_lambda_1se = Actual == Predicted_lambda_1se
)

head(test_predictions)


# =====================================
# 23. Save result tables
# =====================================

write_csv(
  lasso_metrics,
  "Results/ML/LASSO_Overall_Metrics.csv"
)

write_csv(
  lasso_class_metrics,
  "Results/ML/LASSO_Class_Specific_Metrics.csv"
)

write_csv(
  confusion_min_df,
  "Results/ML/LASSO_ConfusionMatrix_lambda_min.csv"
)

write_csv(
  confusion_1se_df,
  "Results/ML/LASSO_ConfusionMatrix_lambda_1se.csv"
)

write_csv(
  test_predictions,
  "Results/ML/LASSO_Test_Predictions.csv"
)

write_csv(
  probability_results,
  "Results/ML/LASSO_Test_Probabilities_lambda_min.csv"
)

write_csv(
  selected_genes_min,
  "Results/ML/LASSO_Selected_Genes_lambda_min.csv"
)

write_csv(
  selected_genes_1se,
  "Results/ML/LASSO_Selected_Genes_lambda_1se.csv"
)

write_csv(
  top_lasso_genes_min,
  "Results/ML/LASSO_Top_Selected_Genes_lambda_min.csv"
)


# =====================================
# 24. Save fitted model and CV settings
# =====================================

saveRDS(
  cv_lasso,
  "Results/ML/LASSO_cv_model.rds"
)

saveRDS(
  foldid,
  "Results/ML/LASSO_5fold_foldid.rds"
)


# =====================================
# 25. Final summary
# =====================================

cat(
  "\nLASSO lambda.min:",
  cv_lasso$lambda.min,
  "\n"
)

cat(
  "LASSO lambda.1se:",
  cv_lasso$lambda.1se,
  "\n"
)

cat(
  "Test accuracy at lambda.min:",
  round(
    unname(cm_min$overall["Accuracy"]),
    3
  ),
  "\n"
)

cat(
  "Test accuracy at lambda.1se:",
  round(
    unname(cm_1se$overall["Accuracy"]),
    3
  ),
  "\n"
)

cat(
  "Unique genes selected at lambda.min:",
  n_distinct(selected_genes_min$Gene),
  "\n"
)

cat(
  "Unique genes selected at lambda.1se:",
  n_distinct(selected_genes_1se$Gene),
  "\n"
)


# =====================================
# 26. Confirm saved outputs
# =====================================

list.files("Results/ML")
list.files("Figures/ML")


#The multinomial LASSO model correctly classified 44 of 46 held-out TCGA samples,
#giving an accuracy of 95.7% and a 95% confidence interval of 85.2%–99.5%. 
#Preprocessing and feature selection were restricted to the training data to 
#minimise direct test-set leakage. Nevertheless, the estimate should be
#interpreted cautiously because the held-out set was small and derived from the 
#same TCGA cohort. Robustness was therefore assessed using label permutation, 
#comparison with Random Forest and planned 
#independent CPTAC validation

 #95.7% necessarily not overfitting?
  
  
  
#the test set was separated before gene selection;
#low-expression filtering used training data only;
#the top 1,000 genes were selected from training data only;
#training means and standard deviations were applied to the test set;
#lambda was chosen by cross-validation within training data;
#the test set was held out from fitting.

#However, the result still needs cautious interpretation because:
  
#the test set contains only 46 samples;
#training and testing both come from TCGA;
#there are many more predictors than samples;
#variable-gene selection occurred before the inner cross-validation folds, although only within the training cohort.
