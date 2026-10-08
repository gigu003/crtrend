test_that("AAPC statistics distinguish a single segment of a multi-segment model", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year_index",
    "rate",
    se = "se",
    min_jp = 1,
    max_jp = 1,
    n_eq = 2
  )
  apc <- calc_apc(fit, ci_method = "Param")
  for (segment in 1:2) {
    aapc <- calc_aapc(
      fit,
      apc$start[segment],
      apc$end[segment],
      ci_method = "Param"
    )
    expect_equal(aapc$statistic, apc$statistic[segment])
    expect_equal(aapc$statistic_type, "t")
    expect_equal(aapc$p_value, apc$p_value[segment])
  }

  outside <- calc_aapc(fit, -1, 21, ci_method = "Param")
  expect_equal(outside$status, "interval_out_of_range")
  expect_identical(outside$statistic, NA_real_)
  expect_identical(outside$statistic_type, NA_character_)
  fit$models[[1L]]$status <- "inference_failed"
  failed <- calc_aapc(fit, 0, 21, ci_method = "Param")
  expect_identical(failed$statistic, NA_real_)
  expect_identical(failed$statistic_type, NA_character_)
})

test_that("AAPC zero-variance statistics are consistent with their p-values", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year_index",
    "rate",
    se = "se",
    min_jp = 0,
    max_jp = 0,
    n_eq = 2
  )
  fit$models[[1L]]$inference$slope_covariance[,] <- 0
  for (slope in c(0, -0.01, 0.01)) {
    fit$models[[1L]]$inference$slope_estimates[] <- slope
    aapc <- calc_aapc(fit, 0, 21, ci_method = "Param")
    expect_equal(aapc$statistic, if (slope == 0) 0 else sign(slope) * Inf)
    expect_equal(aapc$p_value, if (slope == 0) 1 else 0)
  }
})
