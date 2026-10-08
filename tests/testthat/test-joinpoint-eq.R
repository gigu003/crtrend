test_that("the EQ seed produces the native uniform sequence without R RNG changes", {
  set.seed(42)
  state <- .Random.seed
  expect_equal(
    crtrend:::nci_mt_uniforms(10L, 10000L),
    c(
      0.72428673645481467,
      0.55590562429279089,
      0.26564709306694567,
      0.61694240500219166,
      0.041540267411619425,
      0.0031904426869004965,
      0.099032033002004027,
      0.38299144804477692,
      0.52762327436357737,
      0.63866342860274017
    ),
    tolerance = 1e-15
  )
  expect_identical(.Random.seed, state)
})

test_that("EQ matches all fourteen independent Joinpoint 6.1.0 exports", {
  audit_dir <- system.file(
    "extdata",
    "joinpoint-6.1.0-eq-audit",
    package = "crtrend"
  )
  reference <- read.csv(file.path(audit_dir, "nci-6.1-eq-reference.csv"))
  for (case in unique(reference$case)) {
    data <- read.delim(
      file.path(audit_dir, paste0(case, ".txt")),
      check.names = FALSE
    )
    fit <- joinpoint(
      data,
      "Year",
      "Rate",
      se = "Standard Error",
      min_jp = 0,
      max_jp = 1,
      selection_method = "WBIC",
      n_eq = 5001,
      eq_seed = 10000
    )
    apc <- calc_apc(fit)
    aapc <- calc_aapc(fit, min(data$Year), max(data$Year))
    target_apc <- reference[
      reference$case == case & reference$statistic == "APC",
    ]
    target_aapc <- reference[
      reference$case == case & reference$statistic == "AAPC",
    ]
    expect_equal(
      fit$models[[1L]]$n_joinpoints,
      nrow(target_apc) - 1L,
      info = case
    )
    expect_equal(apc$start, target_apc$start, info = case)
    expect_equal(apc$end, target_apc$end, info = case)
    expect_equal(round(apc$apc, 4), target_apc$estimate, info = case)
    expect_equal(round(apc$apc_lower, 4), target_apc$lower, info = case)
    expect_equal(round(apc$apc_upper, 4), target_apc$upper, info = case)
    expect_equal(round(apc$p_value, 6), target_apc$p_value, info = case)
    expect_equal(round(aapc$aapc, 4), target_aapc$estimate, info = case)
    expect_equal(round(aapc$aapc_lower, 4), target_aapc$lower, info = case)
    expect_equal(round(aapc$aapc_upper, 4), target_aapc$upper, info = case)
    expect_equal(round(aapc$p_value, 6), target_aapc$p_value, info = case)
    if (nrow(target_apc) > 1L) {
      tau <- fit$models[[1L]]$resamples$draws$joinpoints[, 1L]
      expect_equal(
        quantile(tau, c(.025, .975), names = FALSE),
        c(target_apc$tau_lower[1L], target_apc$tau_upper[1L]),
        info = case
      )
    }
  }
})

test_that("EQ includes the original estimate in a two-joinpoint DDS model", {
  reference_dir <- system.file(
    "extdata",
    "joinpoint-6.1.0-dds-disagreement",
    package = "crtrend"
  )
  data <- read.csv(file.path(reference_dir, "dds-selects-bic.csv"))
  fit <- joinpoint(
    data,
    "year",
    "rate",
    max_jp = 3,
    selection_method = "DDS",
    n_eq = 5001,
    eq_seed = 10000
  )
  selected <- fit$models[[1L]]
  draws <- selected$resamples$draws
  expect_equal(selected$n_joinpoints, 2L)
  expect_equal(nrow(draws$slopes), 5001L)
  expect_equal(draws$slopes[1L, ], unname(selected$inference$slope_estimates))
  expect_equal(draws$joinpoints[1L, ], selected$joinpoints)

  apc <- calc_apc(fit)
  expect_equal(round(apc$apc, 4), c(-1.9830, -0.0146, 0.2790))
  expect_equal(round(apc$apc_lower, 4), c(-2.7918, -2.6906, -0.4398))
  expect_equal(round(apc$apc_upper, 4), c(-1.3040, 0.8221, 1.1700))
  expect_equal(round(apc$p_value, 6), c(0.030394, 0.482304, 0.061988))

  full <- calc_aapc(fit, 0, 30)
  last_four <- calc_aapc(fit, 26, 30)
  last_segment <- calc_aapc(fit, 21, 30)
  expect_equal(round(full$aapc, 4), -0.5877)
  expect_equal(round(full$aapc_lower, 4), -0.6875)
  expect_equal(round(full$aapc_upper, 4), -0.5029)
  expect_equal(round(full$p_value, 6), 0)
  expect_equal(round(last_four$aapc_lower, 4), -0.3156)
  expect_equal(round(last_four$aapc_upper, 4), 0.8599)
  expect_equal(round(last_four$p_value, 6), 0.061988)
  expect_equal(round(last_segment$aapc_lower, 4), -0.0804)
  expect_equal(round(last_segment$aapc_upper, 4), 0.4781)
  expect_equal(round(last_segment$p_value, 6), 0.060788)
})

test_that("Sample.Data cohort zero matches the supplied 6.1.0 screenshots", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year_index",
    "rate",
    se = "se",
    min_jp = 0,
    max_jp = 4,
    n_eq = 5001,
    eq_seed = 10000
  )
  expect_equal(fit$models[[1L]]$joinpoints, 12)
  eq <- calc_apc(fit)
  eq_average <- calc_aapc(fit, 0, 21)
  expect_equal(round(eq$apc, 4), c(-2.7884, -1.4439))
  expect_equal(round(eq$apc_lower, 4), c(-4.5714, -2.1001))
  expect_equal(round(eq$apc_upper, 4), c(-2.3644, 1.2953))
  expect_equal(round(eq$p_value, 6), c(0.008398, 0.151970))
  expect_equal(round(eq_average$aapc, 4), -2.2144)
  expect_equal(
    round(c(eq_average$aapc_lower, eq_average$aapc_upper), 4),
    c(-2.4710, -1.9807)
  )
  expect_equal(eq_average$p_value, 0)
  expect_identical(eq$statistic, rep(NA_real_, 2L))
  expect_identical(eq$statistic_type, rep(NA_character_, 2L))
  expect_identical(eq_average$statistic, NA_real_)
  expect_identical(eq_average$statistic_type, NA_character_)
  expect_equal(
    quantile(
      fit$models[[1L]]$resamples$draws$joinpoints,
      c(.025, .975),
      names = FALSE
    ),
    c(3, 19)
  )
  param <- calc_apc(fit, ci_method = "Param")
  param_average <- calc_aapc(fit, 0, 21, ci_method = "Param")
  expect_equal(round(param$apc_lower, 4), c(-3.1766, -2.0678))
  expect_equal(round(param$apc_upper, 4), c(-2.3987, -0.8160))
  expect_lt(param$p_value[1L], 1e-6)
  expect_equal(round(param$p_value[2L], 6), 0.000156)
  expect_equal(round(param$statistic, 4), c(-14.9131, -4.8318))
  expect_equal(param$statistic_type, c("t", "t"))
  expect_equal(round(param_average$statistic, 4), -13.2917)
  expect_equal(param_average$statistic_type, "z")
  expect_equal(
    param$p_value,
    2 * pt(-abs(param$statistic), df = fit$models[[1L]]$inference$df)
  )
  expect_equal(param_average$p_value, 2 * pnorm(-abs(param_average$statistic)))
  expect_equal(
    round(c(param_average$aapc_lower, param_average$aapc_upper), 4),
    c(-2.5368, -1.8910)
  )
  expect_lt(param_average$p_value, 1e-6)
})
