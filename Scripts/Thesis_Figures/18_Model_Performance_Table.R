# ============================================================
# UCEC Project
# 18 - Internal and External Model Performance Table
# ============================================================

library(tidyverse)

dir.create(
  "Results/Tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Read internal TCGA metrics
# ------------------------------------------------------------

lasso_raw <- readr::read_csv(
  "Results/ML/LASSO_Overall_Metrics.csv",
  show_col_types = FALSE
)

rf_raw <- readr::read_csv(
  "Results/ML/RandomForest_Overall_Metrics.csv",
  show_col_types = FALSE
)

# Keep lambda.min as the primary LASSO model.
if ("Model" %in% names(lasso_raw) &&
    any(lasso_raw$Model == "LASSO_lambda_min")) {
  lasso_internal <- lasso_raw |>
    dplyr::filter(
      Model == "LASSO_lambda_min"
    )
} else {
  lasso_internal <- lasso_raw |>
    dplyr::slice_head(n = 1)
}

lasso_internal <- lasso_internal |>
  dplyr::transmute(
    Cohort = "TCGA held-out test",
    Model = "LASSO (lambda.min)",
    Accuracy,
    Kappa,
    Macro_F1,
    Balanced_Accuracy = Macro_Balanced_Accuracy,
    Sensitivity = Macro_Sensitivity,
    Specificity = Macro_Specificity
  )

rf_internal <- rf_raw |>
  dplyr::slice_head(n = 1) |>
  dplyr::transmute(
    Cohort = "TCGA held-out test",
    Model = "Random Forest",
    Accuracy,
    Kappa,
    Macro_F1,
    Balanced_Accuracy = Macro_Balanced_Accuracy,
    Sensitivity = Macro_Sensitivity,
    Specificity = Macro_Specificity
  )

# ------------------------------------------------------------
# 2. Read independent CPTAC metrics
# ------------------------------------------------------------

external_raw <- readr::read_csv(
  "Results/External_Validation/CPTAC_Model_Metrics.csv",
  show_col_types = FALSE
)

external_table <- external_raw |>
  dplyr::mutate(
    Model = dplyr::recode(
      Model,
      "LASSO" = "LASSO (lambda.min)",
      "Random_Forest" = "Random Forest",
      .default = Model
    )
  ) |>
  dplyr::transmute(
    Cohort = "CPTAC external evaluation",
    Model,
    Accuracy = External_Accuracy,
    Kappa = External_Kappa,
    Macro_F1,
    Balanced_Accuracy = Macro_Balanced_Accuracy,
    Sensitivity = Macro_Sensitivity,
    Specificity = Macro_Specificity
  )

# ------------------------------------------------------------
# 3. Combine and save a machine-readable table
# ------------------------------------------------------------

performance_table <- dplyr::bind_rows(
  lasso_internal,
  rf_internal,
  external_table
) |>
  dplyr::mutate(
    Cohort = factor(
      Cohort,
      levels = c(
        "TCGA held-out test",
        "CPTAC external evaluation"
      )
    ),
    Model = factor(
      Model,
      levels = c(
        "LASSO (lambda.min)",
        "Random Forest"
      )
    )
  ) |>
  dplyr::arrange(
    Cohort,
    Model
  ) |>
  dplyr::mutate(
    Cohort = as.character(Cohort),
    Model = as.character(Model)
  )

readr::write_csv(
  performance_table,
  "Results/Tables/Model_Performance_Table.csv"
)

print(performance_table)

# ------------------------------------------------------------
# 4. Create an editable Word version where packages are present
# ------------------------------------------------------------

if (
  requireNamespace("flextable", quietly = TRUE) &&
  requireNamespace("officer", quietly = TRUE)
) {
  metric_columns <- c(
    "Accuracy",
    "Kappa",
    "Macro_F1",
    "Balanced_Accuracy",
    "Sensitivity",
    "Specificity"
  )

  performance_ft <- flextable::flextable(
    performance_table
  )

  performance_ft <- flextable::set_header_labels(
    performance_ft,
    Cohort = "Validation cohort",
    Model = "Model",
    Accuracy = "Accuracy",
    Kappa = "Cohen's kappa",
    Macro_F1 = "Macro F1-score",
    Balanced_Accuracy = "Macro balanced accuracy",
    Sensitivity = "Macro sensitivity",
    Specificity = "Macro specificity"
  )

  performance_ft <- flextable::colformat_num(
    performance_ft,
    j = metric_columns,
    digits = 3
  )

  performance_ft <- flextable::merge_v(
    performance_ft,
    j = "Cohort"
  )

  performance_ft <- flextable::valign(
    performance_ft,
    j = "Cohort",
    valign = "center",
    part = "body"
  )

  performance_ft <- flextable::theme_booktabs(
    performance_ft
  )

  performance_ft <- flextable::bold(
    performance_ft,
    part = "header"
  )

  # Bold the better-performing LASSO row in each cohort.
  performance_ft <- flextable::bold(
    performance_ft,
    i = performance_table$Model == "LASSO (lambda.min)",
    j = c("Model", metric_columns),
    bold = TRUE,
    part = "body"
  )

  performance_ft <- flextable::align(
    performance_ft,
    j = c("Cohort", "Model"),
    align = "left",
    part = "body"
  )

  performance_ft <- flextable::align(
    performance_ft,
    j = metric_columns,
    align = "center",
    part = "all"
  )

  performance_ft <- flextable::align(
    performance_ft,
    align = "center",
    part = "header"
  )

  performance_ft <- flextable::font(
    performance_ft,
    fontname = "Arial",
    part = "all"
  )

  performance_ft <- flextable::fontsize(
    performance_ft,
    size = 8.5,
    part = "all"
  )

  performance_ft <- flextable::padding(
    performance_ft,
    padding.top = 5,
    padding.bottom = 5,
    padding.left = 4,
    padding.right = 4,
    part = "all"
  )

  performance_ft <- flextable::hline(
    performance_ft,
    i = 2,
    border = officer::fp_border(
      color = "black",
      width = 0.75
    ),
    part = "body"
  )

  performance_ft <- flextable::add_footer_lines(
    performance_ft,
    values = paste(
      "Note: Bold values indicate the better-performing model within each",
      "evaluation cohort. Macro-averaged metrics are the unweighted mean",
      "across the three UCEC molecular subtypes."
    )
  )

  performance_ft <- flextable::italic(
    performance_ft,
    part = "footer"
  )

  performance_ft <- flextable::autofit(
    performance_ft
  )

  document <- officer::read_docx()

  document <- officer::body_add_par(
    document,
    paste0(
      "Table. Internal and external performance of LASSO and ",
      "Random Forest models for UCEC molecular subtype classification."
    )
  )

  document <- flextable::body_add_flextable(
    document,
    value = performance_ft
  )

  print(
    document,
    target = file.path(
      "Results/Tables",
      "Model_Performance_Table_Editable.docx"
    )
  )
} else {
  message(
    paste(
      "CSV table was created.",
      "Install flextable and officer to also create the editable Word table."
    )
  )
}
