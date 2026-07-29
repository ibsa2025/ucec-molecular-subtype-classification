# ============================================================
# UCEC Project
# 19 - Shared Candidate Biomarker Table
# ============================================================

library(tidyverse)

dir.create(
  "Results/ML",
  recursive = TRUE,
  showWarnings = FALSE
)

lasso_selected <- readr::read_csv(
  "Results/ML/LASSO_Selected_Genes_lambda_min.csv",
  show_col_types = FALSE
)

rf_importance <- readr::read_csv(
  "Results/ML/RandomForest_All_Gene_Importance.csv",
  show_col_types = FALSE
)

lasso_summary <- lasso_selected |>
  dplyr::group_by(Gene) |>
  dplyr::arrange(
    dplyr::desc(Absolute_Coefficient),
    .by_group = TRUE
  ) |>
  dplyr::summarise(
    LASSO_Max_Absolute_Coefficient =
      dplyr::first(Absolute_Coefficient),
    LASSO_Coefficient =
      dplyr::first(Coefficient),
    LASSO_Associated_Subtype =
      dplyr::first(Class),
    LASSO_Number_of_Classes =
      dplyr::n_distinct(Class),
    .groups = "drop"
  )

rf_top30 <- rf_importance |>
  dplyr::arrange(
    dplyr::desc(Importance)
  ) |>
  dplyr::slice_head(n = 30) |>
  dplyr::mutate(
    RF_Rank = dplyr::row_number()
  ) |>
  dplyr::select(
    Gene,
    RF_Importance = Importance,
    RF_Rank
  )

shared_candidate_biomarkers_detailed <- dplyr::inner_join(
  lasso_summary,
  rf_top30,
  by = "Gene"
) |>
  dplyr::arrange(
    RF_Rank,
    dplyr::desc(
      LASSO_Max_Absolute_Coefficient
    )
  )

shared_candidate_biomarkers <- shared_candidate_biomarkers_detailed |>
  dplyr::select(Gene)

readr::write_csv(
  shared_candidate_biomarkers,
  "Results/ML/Shared_Candidate_Biomarkers.csv"
)

readr::write_csv(
  shared_candidate_biomarkers_detailed,
  "Results/ML/Shared_Candidate_Biomarkers_Detailed.csv"
)

stopifnot(
  nrow(shared_candidate_biomarkers) == 11
)

print(shared_candidate_biomarkers_detailed)
