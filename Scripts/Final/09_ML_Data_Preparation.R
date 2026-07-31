# =====================================
# UCEC Project
# 09 - Machine Learning Data Preparation
# =====================================

library(tidyverse)
library(edgeR)
library(caret)

set.seed(123)

dir.create(
  "Results/ML",
  recursive = TRUE,
  showWarnings = FALSE
)

tcga_counts <- readr::read_csv(
  paste0(
    "TCGA Data/subtype_tumor_only/",
    "tcga_ucec_subtype_counts_geneSymbol.csv"
  ),
  show_col_types = FALSE
)

tcga_meta <- readr::read_csv(
  paste0(
    "TCGA Data/subtype_tumor_only/",
    "tcga_ucec_subtype_metadata.csv"
  ),
  show_col_types = FALSE
)

gene_names <- tcga_counts[[1]]

counts_matrix <- tcga_counts[, -1] |>
  as.data.frame() |>
  as.matrix()

rownames(counts_matrix) <- gene_names

if (!setequal(colnames(counts_matrix), tcga_meta$barcode)) {
  stop("TCGA expression samples and metadata barcodes do not match.")
}

tcga_meta <- tcga_meta[
  match(colnames(counts_matrix), tcga_meta$barcode),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(counts_matrix), tcga_meta$barcode)
)

# Resolve duplicated gene symbols using the highest total expression row.
if (anyDuplicated(rownames(counts_matrix)) > 0) {
  total_expression <- rowSums(counts_matrix, na.rm = TRUE)
  keep_order <- order(total_expression, decreasing = TRUE)
  counts_matrix <- counts_matrix[keep_order, , drop = FALSE]
  counts_matrix <- counts_matrix[
    !duplicated(rownames(counts_matrix)),
    ,
    drop = FALSE
  ]
}

outcome <- factor(
  tcga_meta$subtype3,
  levels = c(
    "CN_high",
    "CN_low",
    "Hypermutated"
  )
)

counts_matrix <- counts_matrix[
  !grepl("^MT-", rownames(counts_matrix)),
  ,
  drop = FALSE
]

train_index <- caret::createDataPartition(
  outcome,
  p = 0.80,
  list = FALSE
)

train_counts <- counts_matrix[
  ,
  train_index,
  drop = FALSE
]

test_counts <- counts_matrix[
  ,
  -train_index,
  drop = FALSE
]

y_train <- outcome[train_index]
y_test <- outcome[-train_index]

dge_train <- edgeR::DGEList(
  counts = train_counts,
  group = y_train
)

keep_expression <- edgeR::filterByExpr(
  dge_train,
  group = y_train
)

train_counts_filtered <- train_counts[
  keep_expression,
  ,
  drop = FALSE
]

test_counts_filtered <- test_counts[
  keep_expression,
  ,
  drop = FALSE
]

# The original machine-learning analysis used library-size-adjusted log2 CPM
# without TMM normalisation. This is retained so the code reproduces the
# reported model outputs.
train_cpm <- edgeR::cpm(
  train_counts_filtered,
  log = TRUE,
  prior.count = 1,
  normalized.lib.sizes = FALSE
)

test_cpm <- edgeR::cpm(
  test_counts_filtered,
  log = TRUE,
  prior.count = 1,
  normalized.lib.sizes = FALSE
)

x_train_all <- t(train_cpm)
x_test_all <- t(test_cpm)

training_gene_variance <- apply(
  x_train_all,
  2,
  var
)

training_gene_variance <- sort(
  training_gene_variance,
  decreasing = TRUE
)

number_of_features <- min(
  1000,
  length(training_gene_variance)
)

top_genes <- names(
  training_gene_variance[
    seq_len(number_of_features)
  ]
)

x_train_selected <- x_train_all[
  ,
  top_genes,
  drop = FALSE
]

x_test_selected <- x_test_all[
  ,
  top_genes,
  drop = FALSE
]

train_means <- apply(
  x_train_selected,
  2,
  mean
)

train_sds <- apply(
  x_train_selected,
  2,
  sd
)

train_sds[
  is.na(train_sds) |
    train_sds == 0
] <- 1

x_train_scaled <- sweep(
  x_train_selected,
  2,
  train_means,
  FUN = "-"
)

x_train_scaled <- sweep(
  x_train_scaled,
  2,
  train_sds,
  FUN = "/"
)

x_test_scaled <- sweep(
  x_test_selected,
  2,
  train_means,
  FUN = "-"
)

x_test_scaled <- sweep(
  x_test_scaled,
  2,
  train_sds,
  FUN = "/"
)

stopifnot(
  identical(
    colnames(x_train_scaled),
    colnames(x_test_scaled)
  ),
  sum(is.na(x_train_scaled)) == 0,
  sum(is.na(x_test_scaled)) == 0,
  sum(is.infinite(x_train_scaled)) == 0,
  sum(is.infinite(x_test_scaled)) == 0
)

saveRDS(
  x_train_scaled,
  "Results/ML/x_train_scaled.rds"
)

saveRDS(
  x_test_scaled,
  "Results/ML/x_test_scaled.rds"
)

saveRDS(
  y_train,
  "Results/ML/y_train.rds"
)

saveRDS(
  y_test,
  "Results/ML/y_test.rds"
)

saveRDS(
  top_genes,
  "Results/ML/top1000_training_variable_genes.rds"
)

saveRDS(
  train_index,
  "Results/ML/train_sample_indices.rds"
)

saveRDS(
  train_means,
  "Results/ML/training_feature_means.rds"
)

saveRDS(
  train_sds,
  "Results/ML/training_feature_sds.rds"
)

readr::write_csv(
  tibble::tibble(Gene = top_genes),
  "Results/ML/Top1000_Training_Variable_Genes.csv"
)

readr::write_csv(
  tibble::tibble(
    Barcode = rownames(x_train_scaled),
    Subtype = as.character(y_train),
    Dataset = "Training"
  ),
  "Results/ML/Training_Samples.csv"
)

readr::write_csv(
  tibble::tibble(
    Barcode = rownames(x_test_scaled),
    Subtype = as.character(y_test),
    Dataset = "Test"
  ),
  "Results/ML/Test_Samples.csv"
)

preprocessing_summary <- tibble::tibble(
  Measure = c(
    "TCGA samples",
    "Training samples",
    "Held-out test samples",
    "Genes after training-only expression filtering",
    "Selected training-variable predictors"
  ),
  Value = c(
    length(outcome),
    length(y_train),
    length(y_test),
    sum(keep_expression),
    length(top_genes)
  )
)

readr::write_csv(
  preprocessing_summary,
  "Results/ML/ML_Preprocessing_Summary.csv"
)

cat("\nTraining subtype distribution:\n")
print(table(y_train))

cat("\nHeld-out test subtype distribution:\n")
print(table(y_test))

print(preprocessing_summary)
