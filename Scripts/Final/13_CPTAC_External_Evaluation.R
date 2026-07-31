# =====================================
# UCEC Project
# 13 - CPTAC External Evaluation
# =====================================

# Aim:
# To evaluate whether the TCGA-trained multinomial LASSO and Random Forest
# models can classify molecular subtypes in an independent CPTAC-UCEC cohort.
#
# The CPTAC RNA file contains genes as rows and samples as columns.
# The expression values appear to be previously normalised and transformed,
# so they are not passed through edgeR CPM or log2 transformation again.
#
# CPTAC samples are aligned with their metadata, the same 1000 genes used
# during TCGA model development are extracted, and TCGA training means and
# standard deviations are applied before external prediction.
#
# Important:
# TCGA and CPTAC may differ in platform and preprocessing. Therefore, this is
# an external transfer analysis and performance should be interpreted alongside
# gene coverage and expression-scale compatibility.


# =====================================
# 1. Load packages
# =====================================

library(tidyverse)
library(glmnet)
library(caret)
library(ranger)


# =====================================
# 2. Set seed for reproducibility
# =====================================

set.seed(123)


# =====================================
# 3. Define CPTAC file paths
# =====================================

cptac_rna_path <- paste0(
  "CPTAC Data/subtype_tumor_only/",
  "cptac_ucec_subtype_rna.csv"
)

cptac_meta_path <- paste0(
  "CPTAC Data/subtype_tumor_only/",
  "cptac_ucec_subtype_metadata.csv"
)


# =====================================
# 4. Load CPTAC RNA and metadata
# =====================================

cptac_rna <- read_csv(
  cptac_rna_path,
  show_col_types = FALSE
)

cptac_meta <- read_csv(
  cptac_meta_path,
  show_col_types = FALSE
)


# =====================================
# 5. Inspect imported CPTAC data
# =====================================

dim(cptac_rna)
dim(cptac_meta)

head(cptac_rna[, 1:6])
head(cptac_meta)

names(cptac_meta)

# Expected:
# CPTAC RNA = 28,057 genes x 102 columns
#             one feature_id column and 101 sample columns
#
# CPTAC metadata = 101 samples x 22 variables


# =====================================
# 6. Confirm CPTAC expression type
# =====================================

rna_numeric_values <- cptac_rna %>%
  select(where(is.numeric)) %>%
  unlist(use.names = FALSE)

summary(rna_numeric_values)

range(
  rna_numeric_values,
  na.rm = TRUE
)

all_whole_numbers <- all(
  rna_numeric_values == floor(rna_numeric_values),
  na.rm = TRUE
)

contains_negative_values <- any(
  rna_numeric_values < 0,
  na.rm = TRUE
)

all_whole_numbers
contains_negative_values

# Your observed values were:
# range = 0 to 22.61
# whole numbers = FALSE
# negative values = FALSE
#
# This is consistent with already transformed/normalised expression data.
# Therefore, the CPTAC data will be used as supplied.
#
# Do not run CPM, voom or log2(x + 1) again here.


# =====================================
# 7. Create CPTAC expression matrix
# =====================================

# The first column, feature_id, contains gene symbols.

cptac_gene_names <- cptac_rna$feature_id

cptac_expression <- cptac_rna %>%
  select(-feature_id) %>%
  as.data.frame() %>%
  as.matrix()

rownames(cptac_expression) <- cptac_gene_names

dim(cptac_expression)

head(
  rownames(cptac_expression)
)

head(
  colnames(cptac_expression)
)


# =====================================
# 8. Deal with duplicated gene symbols
# =====================================

duplicated_gene_count <- sum(
  duplicated(
    rownames(cptac_expression)
  )
)

duplicated_gene_count

# If duplicated gene symbols exist, retain the row with the highest
# mean expression across CPTAC samples.

if (duplicated_gene_count > 0) {
  
  cptac_expression_df <- as.data.frame(
    cptac_expression
  ) %>%
    rownames_to_column("Gene") %>%
    mutate(
      Mean_Expression = rowMeans(
        across(
          where(is.numeric)
        ),
        na.rm = TRUE
      )
    ) %>%
    arrange(
      Gene,
      desc(Mean_Expression)
    ) %>%
    distinct(
      Gene,
      .keep_all = TRUE
    ) %>%
    select(
      -Mean_Expression
    )
  
  cptac_expression <- cptac_expression_df %>%
    column_to_rownames("Gene") %>%
    as.matrix()
}

dim(cptac_expression)


# =====================================
# 9. Prepare CPTAC subtype metadata
# =====================================

# subtype3_tcga already contains the TCGA-compatible subtype names:
#
# CN_high
# CN_low
# Hypermutated

table(
  cptac_meta$subtype3_tcga,
  useNA = "ifany"
)

cptac_meta_model <- cptac_meta %>%
  filter(
    subtype3_tcga %in%
      c(
        "CN_high",
        "CN_low",
        "Hypermutated"
      )
  ) %>%
  mutate(
    External_Subtype = factor(
      subtype3_tcga,
      levels = c(
        "CN_high",
        "CN_low",
        "Hypermutated"
      )
    )
  )

table(
  cptac_meta_model$External_Subtype
)


# =====================================
# 10. Match CPTAC RNA samples to metadata
# =====================================

shared_samples <- intersect(
  colnames(cptac_expression),
  cptac_meta_model$sample_id
)

length(shared_samples)

head(shared_samples)

if (length(shared_samples) == 0) {
  
  stop(
    paste(
      "No CPTAC RNA sample IDs matched the metadata.",
      "Check sample_id formatting."
    )
  )
}

# Reorder metadata and RNA into the same sample order.

cptac_meta_model <- cptac_meta_model[
  match(
    shared_samples,
    cptac_meta_model$sample_id
  ),
]

cptac_expression <- cptac_expression[
  ,
  shared_samples,
  drop = FALSE
]

all(
  colnames(cptac_expression) ==
    cptac_meta_model$sample_id
)

table(
  cptac_meta_model$External_Subtype
)


# =====================================
# 11. Remove mitochondrial genes
# =====================================

keep_non_mito <- !grepl(
  "^MT-",
  rownames(cptac_expression)
)

cptac_expression <- cptac_expression[
  keep_non_mito,
  ,
  drop = FALSE
]

dim(cptac_expression)


# =====================================
# 12. Load TCGA-trained models
# =====================================

lasso_model <- readRDS(
  "Results/ML/LASSO_cv_model.rds"
)

rf_model <- readRDS(
  "Results/ML/RandomForest_Model.rds"
)


# =====================================
# 13. Load TCGA preprocessing objects
# =====================================

tcga_top_genes <- readRDS(
  "Results/ML/top1000_training_variable_genes.rds"
)

tcga_training_means <- readRDS(
  "Results/ML/training_feature_means.rds"
)

tcga_training_sds <- readRDS(
  "Results/ML/training_feature_sds.rds"
)

length(tcga_top_genes)

head(tcga_top_genes)


# =====================================
# 14. Check gene overlap
# =====================================

genes_present <- intersect(
  tcga_top_genes,
  rownames(cptac_expression)
)

genes_missing <- setdiff(
  tcga_top_genes,
  rownames(cptac_expression)
)

gene_coverage_percent <- (
  length(genes_present) /
    length(tcga_top_genes)
) * 100

length(genes_present)
length(genes_missing)
gene_coverage_percent

head(
  genes_missing,
  30
)

# Higher coverage is preferable.
# Missing predictors are assigned the TCGA training mean below.
# Following scaling, these missing genes therefore receive a neutral value of 0.


# =====================================
# 15. Transpose CPTAC matrix
# =====================================

# Machine-learning models require:
# rows = samples
# columns = genes

cptac_x_available <- t(
  cptac_expression
)

dim(cptac_x_available)


# =====================================
# 16. Construct complete 1000-gene matrix
# =====================================

cptac_complete <- matrix(
  NA_real_,
  nrow = nrow(cptac_x_available),
  ncol = length(tcga_top_genes),
  dimnames = list(
    rownames(cptac_x_available),
    tcga_top_genes
  )
)

# Insert genes available in CPTAC.

cptac_complete[
  ,
  genes_present
] <- cptac_x_available[
  ,
  genes_present,
  drop = FALSE
]

# Fill missing CPTAC genes with their TCGA training means.
# This becomes zero after standardisation.

if (length(genes_missing) > 0) {
  
  for (missing_gene in genes_missing) {
    
    cptac_complete[
      ,
      missing_gene
    ] <- tcga_training_means[
      missing_gene
    ]
  }
}

dim(cptac_complete)

sum(
  is.na(cptac_complete)
)


# =====================================
# 17. Apply TCGA training scaling
# =====================================

# Do not calculate new means or SDs from CPTAC for the primary validation.
# The original TCGA training values are applied to simulate prediction
# of genuinely unseen samples.

cptac_scaled <- sweep(
  cptac_complete,
  2,
  tcga_training_means[
    tcga_top_genes
  ],
  FUN = "-"
)

cptac_scaled <- sweep(
  cptac_scaled,
  2,
  tcga_training_sds[
    tcga_top_genes
  ],
  FUN = "/"
)

dim(cptac_scaled)

identical(
  colnames(cptac_scaled),
  tcga_top_genes
)

sum(
  is.na(cptac_scaled)
)

sum(
  is.infinite(cptac_scaled)
)


# =====================================
# 18. Inspect cross-cohort scale compatibility
# =====================================

# Because TCGA and CPTAC may have undergone different expression preprocessing,
# inspect how extreme the CPTAC standardised values are.

cptac_scaled_summary <- summary(
  as.vector(cptac_scaled)
)

cptac_scaled_summary

extreme_value_percentage <- mean(
  abs(cptac_scaled) > 5
) * 100

extreme_value_percentage

# A large percentage of values above |5| would suggest considerable
# preprocessing or platform mismatch between TCGA and CPTAC.
#
# External accuracy should then be interpreted as an exploratory transfer
# result rather than definitive clinical validation.


# =====================================
# 19. External LASSO predictions
# =====================================

cptac_lasso_raw <- predict(
  lasso_model,
  newx = as.matrix(cptac_scaled),
  s = "lambda.min",
  type = "class"
)

cptac_lasso_predictions <- factor(
  as.vector(cptac_lasso_raw),
  levels = levels(
    cptac_meta_model$External_Subtype
  )
)

table(
  cptac_lasso_predictions
)


# =====================================
# 20. LASSO external confusion matrix
# =====================================

cptac_lasso_confusion <- confusionMatrix(
  data = cptac_lasso_predictions,
  reference = cptac_meta_model$External_Subtype
)

cptac_lasso_confusion


# =====================================
# 21. LASSO external probabilities
# =====================================

cptac_lasso_prob_raw <- predict(
  lasso_model,
  newx = as.matrix(cptac_scaled),
  s = "lambda.min",
  type = "response"
)

cptac_lasso_prob <- drop(
  cptac_lasso_prob_raw
) %>%
  as.data.frame()

colnames(cptac_lasso_prob) <- levels(
  cptac_meta_model$External_Subtype
)

head(cptac_lasso_prob)


# =====================================
# 22. External Random Forest predictions
# =====================================

cptac_rf_predictions <- predict(
  rf_model,
  newdata = as.data.frame(cptac_scaled),
  type = "raw"
)

cptac_rf_predictions <- factor(
  cptac_rf_predictions,
  levels = levels(
    cptac_meta_model$External_Subtype
  )
)

table(
  cptac_rf_predictions
)


# =====================================
# 23. Random Forest external confusion matrix
# =====================================

cptac_rf_confusion <- confusionMatrix(
  data = cptac_rf_predictions,
  reference = cptac_meta_model$External_Subtype
)

cptac_rf_confusion


# =====================================
# 24. Random Forest external probabilities
# =====================================

cptac_rf_prob <- predict(
  rf_model,
  newdata = as.data.frame(cptac_scaled),
  type = "prob"
)

head(cptac_rf_prob)


# =====================================
# 25. Function to extract model metrics
# =====================================

extract_external_metrics <- function(
    confusion_object,
    model_name
) {
  
  overall <- confusion_object$overall
  
  class_metrics <- as.matrix(
    confusion_object$byClass
  )
  
  tibble(
    Model = model_name,
    
    External_Accuracy = unname(
      overall["Accuracy"]
    ),
    
    External_Kappa = unname(
      overall["Kappa"]
    ),
    
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


# =====================================
# 26. Extract external model metrics
# =====================================

external_lasso_metrics <- extract_external_metrics(
  cptac_lasso_confusion,
  "LASSO"
)

external_rf_metrics <- extract_external_metrics(
  cptac_rf_confusion,
  "Random Forest"
)

external_model_metrics <- bind_rows(
  external_lasso_metrics,
  external_rf_metrics
)

external_model_metrics


# =====================================
# 27. Create external prediction table
# =====================================

external_predictions <- tibble(
  Sample_ID = rownames(cptac_scaled),
  
  Actual = as.character(
    cptac_meta_model$External_Subtype
  ),
  
  LASSO_Predicted = as.character(
    cptac_lasso_predictions
  ),
  
  Random_Forest_Predicted = as.character(
    cptac_rf_predictions
  ),
  
  LASSO_Correct =
    Actual == LASSO_Predicted,
  
  Random_Forest_Correct =
    Actual == Random_Forest_Predicted
)

external_predictions <- external_predictions %>%
  bind_cols(
    cptac_lasso_prob %>%
      rename_with(
        ~ paste0(
          "LASSO_Probability_",
          .
        )
      )
  ) %>%
  bind_cols(
    as.data.frame(cptac_rf_prob) %>%
      rename_with(
        ~ paste0(
          "RF_Probability_",
          .
        )
      )
  )

head(external_predictions)


# =====================================
# 28. Create output folders
# =====================================

dir.create(
  "Results/External_Validation",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/External_Validation",
  recursive = TRUE,
  showWarnings = FALSE
)


# =====================================
# 29. Prepare LASSO confusion-matrix table
# =====================================

lasso_external_cm_df <- as.data.frame(
  cptac_lasso_confusion$table
) %>%
  rename(
    Predicted = Prediction,
    Actual = Reference,
    Count = Freq
  )

lasso_external_cm_df


# =====================================
# 30. Prepare RF confusion-matrix table
# =====================================

rf_external_cm_df <- as.data.frame(
  cptac_rf_confusion$table
) %>%
  rename(
    Predicted = Prediction,
    Actual = Reference,
    Count = Freq
  )

rf_external_cm_df


# =====================================
# 31. Plot external LASSO confusion matrix
# =====================================

p_external_lasso <- ggplot(
  lasso_external_cm_df,
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
    title = "CPTAC External Validation: LASSO",
    x = "Actual subtype",
    y = "Predicted subtype"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

p_external_lasso

ggsave(
  "Figures/External_Validation/CPTAC_LASSO_ConfusionMatrix.png",
  plot = p_external_lasso,
  width = 8,
  height = 6,
  dpi = 300
)


# =====================================
# 32. Plot external RF confusion matrix
# =====================================

p_external_rf <- ggplot(
  rf_external_cm_df,
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
    title = "CPTAC External Validation: Random Forest",
    x = "Actual subtype",
    y = "Predicted subtype"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

p_external_rf

ggsave(
  "Figures/External_Validation/CPTAC_RandomForest_ConfusionMatrix.png",
  plot = p_external_rf,
  width = 8,
  height = 6,
  dpi = 300
)


# =====================================
# 33. Compare internal and external accuracy
# =====================================

internal_metrics <- read_csv(
  "Results/ML/Model_Comparison_Primary.csv",
  show_col_types = FALSE
)

internal_external_comparison <- internal_metrics %>%
  select(
    Model,
    Internal_Accuracy = Accuracy
  ) %>%
  left_join(
    external_model_metrics %>%
      select(
        Model,
        External_Accuracy
      ),
    by = "Model"
  ) %>%
  pivot_longer(
    cols = c(
      Internal_Accuracy,
      External_Accuracy
    ),
    names_to = "Validation",
    values_to = "Accuracy"
  ) %>%
  mutate(
    Validation = recode(
      Validation,
      
      "Internal_Accuracy" =
        "Internal TCGA test",
      
      "External_Accuracy" =
        "External CPTAC"
    ),
    
    Accuracy_Percent =
      Accuracy * 100
  )

internal_external_comparison


# =====================================
# 34. Plot internal vs external accuracy
# =====================================

p_internal_external <- ggplot(
  internal_external_comparison,
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
    limits = c(
      0,
      105
    )
  ) +
  labs(
    title = "Internal and External Model Performance",
    x = "Model",
    y = "Accuracy (%)",
    fill = "Validation dataset"
  )

p_internal_external

ggsave(
  "Figures/External_Validation/Internal_vs_External_Accuracy.png",
  plot = p_internal_external,
  width = 9,
  height = 7,
  dpi = 300
)


# =====================================
# 35. Save external evaluation results
# =====================================

write_csv(
  external_model_metrics,
  "Results/External_Validation/CPTAC_Model_Metrics.csv"
)

write_csv(
  external_predictions,
  "Results/External_Validation/CPTAC_External_Predictions.csv"
)

write_csv(
  lasso_external_cm_df,
  "Results/External_Validation/CPTAC_LASSO_ConfusionMatrix.csv"
)

write_csv(
  rf_external_cm_df,
  "Results/External_Validation/CPTAC_RandomForest_ConfusionMatrix.csv"
)

write_csv(
  internal_external_comparison,
  "Results/External_Validation/Internal_vs_External_Accuracy.csv"
)

write_csv(
  tibble(
    Total_TCGA_Model_Genes =
      length(tcga_top_genes),
    
    Genes_Present_in_CPTAC =
      length(genes_present),
    
    Genes_Missing_from_CPTAC =
      length(genes_missing),
    
    Gene_Coverage_Percent =
      gene_coverage_percent,
    
    Extreme_Scaled_Value_Percent =
      extreme_value_percentage
  ),
  "Results/External_Validation/CPTAC_Gene_and_Scale_Compatibility.csv"
)


# =====================================
# 36. Final external-validation summary
# =====================================

cat(
  "\nCPTAC samples included:",
  nrow(cptac_scaled),
  "\n"
)

cat(
  "CPTAC subtype distribution:\n"
)

print(
  table(
    cptac_meta_model$External_Subtype
  )
)

cat(
  "\nTCGA model genes present in CPTAC:",
  length(genes_present),
  "of",
  length(tcga_top_genes),
  "\n"
)

cat(
  "Gene coverage:",
  round(
    gene_coverage_percent,
    1
  ),
  "%\n"
)

cat(
  "CPTAC values outside ±5 TCGA SDs:",
  round(
    extreme_value_percentage,
    1
  ),
  "%\n"
)

cat(
  "CPTAC LASSO accuracy:",
  round(
    unname(
      cptac_lasso_confusion$overall[
        "Accuracy"
      ]
    ) * 100,
    1
  ),
  "%\n"
)

cat(
  "CPTAC Random Forest accuracy:",
  round(
    unname(
      cptac_rf_confusion$overall[
        "Accuracy"
      ]
    ) * 100,
    1
  ),
  "%\n"
)


# =====================================
# 37. Confirm saved outputs
# =====================================

list.files(
  "Results/External_Validation"
)

list.files(
  "Figures/External_Validation"
)
