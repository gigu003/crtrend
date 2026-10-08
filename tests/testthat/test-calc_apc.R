test_that("APC statistics retain signed zero-variance and failed-group conventions", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year_index",
    "rate",
    se = "se",
    min_jp = 1,
    max_jp = 1,
    n_eq = 2
  )
  fit$models[[1L]]$inference$slope_estimates <- c(0, -0.01)
  fit$models[[1L]]$inference$slope_se <- c(0, 0)
  param <- calc_apc(fit, ci_method = "Param")
  expect_equal(param$statistic, c(0, -Inf))
  expect_equal(param$p_value, c(1, 0))

  fit$models[[2L]] <- fit$models[[1L]]
  fit$models[[2L]]$status <- "inference_failed"
  mixed <- calc_apc(fit, ci_method = "Param")
  expect_equal(mixed$statistic, c(0, -Inf, NA_real_))
  expect_identical(mixed$statistic_type, c("t", "t", NA_character_))
})
