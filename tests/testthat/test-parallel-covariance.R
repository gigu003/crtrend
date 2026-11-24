test_that("parallel covariance fits match serial fits for shuffled groups", {
  skip_on_cran()
  data <- subset(jp_data, sex %in% c(0, 1))
  data <- data[order(-data$year, data$sex), ]
  rownames(data) <- NULL
  index <- seq_len(nrow(data))
  log_covariance <- 0.0025 * 0.4^abs(outer(index, index, `-`))

  for (scale in c("log", "outcome")) {
    covariance <- if (scale == "log") {
      log_covariance
    } else {
      log_covariance * tcrossprod(data$rate)
    }
    args <- list(
      data = data,
      time = "year",
      outcome = "rate",
      group_vars = "sex",
      error = if (scale == "log") "GLS" else "Survey",
      cov = covariance,
      cov_scale = scale,
      max_jp = 1L,
      selection_method = "BIC",
      n_eq = 9L
    )
    serial <- do.call(joinpoint, args)
    concurrent <- do.call(
      joinpoint,
      c(args, list(parallel = TRUE, workers = 2L))
    )

    expect_identical(concurrent$model_info, serial$model_info)
    expect_identical(fitted(concurrent), fitted(serial))
    expect_identical(concurrent$selection_tests, serial$selection_tests)
    expect_identical(
      lapply(concurrent$models, function(model) model$resamples$draws),
      lapply(serial$models, function(model) model$resamples$draws)
    )
  }
})
