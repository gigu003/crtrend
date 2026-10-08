test_that("reserved grouping names are rejected before labels are lost", {
  d <- data.frame(year = 2000:2005, rate = 1:6, status = "cohortA")
  expect_snapshot(
    error = TRUE,
    joinpoint(d, "year", "rate", "status", n_eq = 2)
  )
  for (name in c(
    "group_id",
    "segment",
    "trend_time",
    "apc_p_value",
    "inference_status"
  )) {
    message <- tryCatch(
      validate_joinpoint_group_names(name),
      error = conditionMessage
    )
    expect_match(message, "Reserved grouping column")
  }
})

test_that("standardization requires complete ages and common normalized weights", {
  d <- data.frame(
    year = c(2000, 2000, 2001, 2001),
    age = c("young", "old", "young", "old"),
    count = c(10, 100, 10, 100),
    population = 1000,
    standard_population = 50
  )
  expect_equal(
    calc_age_adjusted_rate(d, group_vars = "year")$rate,
    c(5500, 5500)
  )
  expect_snapshot(
    error = TRUE,
    calc_age_adjusted_rate(d[-4, ], group_vars = "year")
  )
  expect_snapshot(
    error = TRUE,
    calc_age_adjusted_rate(
      d,
      group_vars = "year",
      expected_ages = c("young", "old", "other")
    )
  )
  d$standard_population[4] <- 100
  expect_snapshot(error = TRUE, calc_age_adjusted_rate(d, group_vars = "year"))
  d$standard_population[3] <- 100
  expect_equal(
    calc_age_adjusted_rate(d[c(4, 3, 2, 1), ], group_vars = "year")$rate,
    c(5500, 5500)
  )
})

test_that("zero count groups do not interrupt positive groups", {
  d <- data.frame(
    year = rep(2000:2005, 3),
    cohort = rep(c("zero", "mixed", "valid"), each = 6),
    count = c(rep(0, 6), 0, 2:6, 1:6),
    population = 1000
  )
  rates <- calc_crude_rate(d, group_vars = c("year", "cohort"))
  fit <- joinpoint(
    rates,
    "year",
    "rate",
    "cohort",
    se = "se",
    max_jp = 0,
    n_eq = 19
  )
  status <- setNames(fit$model_info$status, fit$model_info$cohort)
  expect_equal(
    unname(status[c("zero", "mixed", "valid")]),
    c("zero_risk", "nonpositive_outcome", "ok")
  )
  expect_equal(calc_apc(fit)$status, fit$model_info$status)
  expect_equal(calc_model(fit)$statistics$status, fit$model_info$status)
  rates$se[rates$cohort == "valid"] <- 0
  expect_snapshot(
    error = TRUE,
    joinpoint(rates, "year", "rate", "cohort", se = "se", max_jp = 0, n_eq = 19)
  )
})

test_that("flat EQ distributions do not report significance", {
  d <- data.frame(year = 2000:2005, rate = 1, se = .1)
  fit <- joinpoint(d, "year", "rate", se = "se", max_jp = 0, n_eq = 19)
  apc <- calc_apc(fit)
  aapc <- calc_aapc(fit, 2000, 2005)
  expect_equal(apc$apc, 0)
  expect_equal(c(apc$apc_lower, apc$apc_upper), c(0, 0))
  expect_identical(apc$p_value, NA_real_)
  expect_identical(aapc$p_value, NA_real_)
  expect_equal(empirical_two_sided_p(c(-1, 0, 1, 2)), .5)
})

test_that("inference uses cached counts and does not reset the RNG", {
  d <- data.frame(year = 2000:2005, rate = c(10, 11, 13, 14, 15, 18), se = .1)
  fit <- joinpoint(d, "year", "rate", se = "se", max_jp = 0, n_eq = 19)
  set.seed(123)
  before <- .Random.seed
  expect_equal(calc_apc(fit), calc_apc(fit, n_resamples = 19))
  expect_equal(
    calc_aapc(fit, 2000, 2005),
    calc_aapc(fit, 2000, 2005, n_resamples = 19)
  )
  for (method in c("Param", "EQ")) {
    calc_apc(fit, ci_method = method, n_resamples = 23)
    calc_aapc(fit, 2000, 2005, ci_method = method, n_resamples = 23)
    expect_identical(.Random.seed, before)
    calc_apc(fit, ci_method = method, n_resamples = 23, seed = 1234)
    calc_aapc(
      fit,
      2000,
      2005,
      ci_method = method,
      n_resamples = 23,
      seed = 1234
    )
    expect_identical(.Random.seed, before)
  }
  fit$control$eq_seed <- NULL
  calc_apc(fit, n_resamples = 23)
  expect_identical(identical(.Random.seed, before), FALSE)
})
