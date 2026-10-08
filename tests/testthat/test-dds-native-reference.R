test_that("DDS branches and BIC scores match Joinpoint 6.1.0 exports", {
  reference_dir <- system.file(
    "extdata",
    "joinpoint-6.1.0-dds-disagreement",
    package = "crtrend"
  )
  if (!nzchar(reference_dir)) {
    reference_dir <- test_path(
      "..",
      "..",
      "inst",
      "extdata",
      "joinpoint-6.1.0-dds-disagreement"
    )
  }
  candidate_ref <- read.csv(file.path(reference_dir, "native-reference.csv"))
  selection_ref <- read.csv(file.path(reference_dir, "native-dds-results.csv"))
  cases <- c("bic", "bic3")
  expected_method <- c("BIC", "BIC3")
  for (index in seq_along(cases)) {
    case <- cases[index]
    input <- read.csv(file.path(
      reference_dir,
      paste0("dds-selects-", case, ".csv")
    ))
    fit <- joinpoint(
      input,
      "year",
      "rate",
      max_jp = 3,
      selection_method = "DDS",
      n_eq = 2
    )
    result <- calc_model(fit, ci_method = "Param")
    candidate <- candidate_ref[candidate_ref$case == case, ]
    native <- selection_ref[selection_ref$case == case, ]

    expect_equal(result$candidates$n_joinpoints, candidate$n_joinpoints)
    expect_equal(round(result$candidates$sse, 7), candidate$sse)
    expect_equal(round(result$candidates$bic, 7), candidate$bic)
    expect_equal(round(result$candidates$bic3, 7), candidate$bic3)
    expect_equal(
      result$selection$dds_method,
      expected_method[index]
    )
    expect_equal(
      result$selection$selected_n_joinpoints,
      native$final_selected_k
    )
    expect_lt(
      abs(result$selection$dds_bic_delta - native$delta_bic),
      6e-6
    )
    expect_lt(
      abs(result$selection$dds_bic3_delta - native$delta_bic3),
      6e-6
    )
  }
})
