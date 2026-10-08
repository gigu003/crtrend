test_that("opt-in weights reproduce the NCI 6.1 screenshot", {
  data <- subset(jp_data, sex == 0)
  fit <- joinpoint(
    data,
    "year",
    "rate",
    se = "se",
    n_eq = 9,
    wbic_weight = "joinpoint_6.1"
  )
  stable <- joinpoint(data, "year", "rate", se = "se", n_eq = 9)
  explicit_qr <- joinpoint(
    data,
    "year",
    "rate",
    se = "se",
    n_eq = 9,
    wbic_weight = "qr"
  )
  result <- calc_model(fit, ci_method = "Param")
  baseline <- calc_model(stable, ci_method = "Param")
  expect_equal(
    round(result$candidates$weight, 7),
    c(NA, .5184867, .6844505, .9146089, .9116191)
  )
  expect_equal(
    round(result$candidates$weighted_bic, 7),
    c(4.2970652, 3.9200959, 4.2041850, 4.6405566, 5.0180119)
  )
  expect_equal(result$candidates$weight_method, rep("joinpoint_6.1", 5))
  expect_equal(baseline$candidates$weight_method, rep("qr", 5))
  expect_equal(fit$control$wbic_weight, "joinpoint_6.1")
  expect_equal(fit$models[[1]]$control$wbic_weight, "joinpoint_6.1")
  expect_equal(stable$model_info, explicit_qr$model_info, tolerance = 0)
  expect_equal(
    stable$models[[1]]$selection$candidates,
    explicit_qr$models[[1]]$selection$candidates,
    tolerance = 0
  )
  expect_equal(
    result$candidates[c("sse", "bic", "bic3", "modified_bic", "joinpoints")],
    baseline$candidates[c("sse", "bic", "bic3", "modified_bic", "joinpoints")],
    tolerance = 0
  )
  expect_equal(result$joinpoints, stable$joinpoints)
  expect_equal(result$coefficients, baseline$coefficients, tolerance = 0)
  expect_equal(fitted(fit), fitted(stable), tolerance = 0)
  expect_equal(
    calc_apc(fit, ci_method = "Param"),
    calc_apc(stable, ci_method = "Param"),
    tolerance = 0
  )
  expect_equal(
    calc_aapc(fit, start = 2000, end = 2021, ci_method = "Param"),
    calc_aapc(stable, start = 2000, end = 2021, ci_method = "Param"),
    tolerance = 0
  )
  expect_equal(calc_apc(fit), calc_apc(stable), tolerance = 0)
  expect_equal(
    calc_aapc(fit, start = 2000, end = 2021),
    calc_aapc(stable, start = 2000, end = 2021),
    tolerance = 0
  )
  expect_equal(result$candidates$selected, c(FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(result$candidates$weighted_bic[2], fit$model_info$weighted_bic)
})

test_that("compatibility propagates through grouped and sorted data", {
  data <- jp_data[nrow(jp_data):1, ]
  fit <- joinpoint(
    data,
    "year",
    "rate",
    group_vars = "sex",
    se = "se",
    n_eq = 2,
    wbic_weight = "joinpoint_6.1"
  )
  result <- calc_model(fit, ci_method = "Param")
  expect_equal(unique(result$candidates$weight_method), "joinpoint_6.1")
  expect_equal(
    round(subset(result$candidates, sex == 0)$weight, 7),
    c(NA, .5184867, .6844505, .9146089, .9116191)
  )
  expect_equal(
    vapply(fit$models, \(x) x$control$wbic_weight, character(1)) |> unname(),
    rep("joinpoint_6.1", 3)
  )
})

test_that("compatibility restrictions are explicit", {
  data <- subset(jp_data, sex == 0)
  expect_snapshot(error = TRUE, {
    joinpoint(data, "year", "rate", wbic_weight = "unknown")
  })
  expect_snapshot(error = TRUE, {
    joinpoint(data, "year", "rate", wbic_weight = "joinpoint_6.1")
  })
  for (error_structure in c("AR1", "GLS", "Survey")) {
    expect_snapshot(error = TRUE, {
      joinpoint(
        data,
        "year",
        "rate",
        se = "se",
        error = error_structure,
        wbic_weight = "joinpoint_6.1"
      )
    })
  }
})

test_that("native-order kernels and numerical failure guards work", {
  a <- matrix(c(3, 1, 1, 2), 2)
  expect_equal(crtrend:::joinpoint_61_inverse(a), solve(a), tolerance = 1e-15)
  expect_equal(crtrend:::joinpoint_61_multiply(a, a), a %*% a, tolerance = 0)
  expect_snapshot(error = TRUE, {
    crtrend:::joinpoint_61_inverse(matrix(1, 2, 2))
  })
  expect_snapshot(error = TRUE, {
    crtrend:::joinpoint_61_inverse(matrix(Inf, 2, 2))
  })
  data <- data.frame(
    trend_time = 0:6,
    trend_outcome = rep(2, 7),
    trend_se = rep(1e-200, 7),
    trend_log = rep(log(2), 7)
  )
  expect_snapshot(error = TRUE, {
    crtrend:::joinpoint_61_partial_weight(data, pmax(data$trend_time - 3, 0))
  })
  data$trend_se <- rep(.1, 7)
  data$trend_log <- rep(0, 7)
  expect_snapshot(error = TRUE, {
    crtrend:::joinpoint_61_partial_weight(data, pmax(data$trend_time - 3, 0))
  })
})

test_that("zero-joinpoint compatibility uses BIC without an undefined weight", {
  data <- data.frame(year = 2000:2010, rate = 1, se = .1)
  fit <- joinpoint(
    data,
    "year",
    "rate",
    se = "se",
    max_jp = 0,
    n_eq = 2,
    wbic_weight = "joinpoint_6.1"
  )
  candidates <- calc_model(fit, ci_method = "Param")$candidates
  expect_equal(candidates$weight, NA_real_)
  expect_equal(candidates$weighted_bic, candidates$bic)
  expect_equal(candidates$weighted_bic, -Inf)
  expect_equal(candidates$weight_method, "joinpoint_6.1")
})

test_that("centered compatibility agrees with stable QR and legacy fits extract", {
  data <- subset(jp_data, sex == 0)
  data$year <- data$year - 2000
  native <- joinpoint(
    data,
    "year",
    "rate",
    se = "se",
    n_eq = 2,
    wbic_weight = "joinpoint_6.1"
  )
  stable <- joinpoint(data, "year", "rate", se = "se", n_eq = 2)
  expect_equal(
    native$models[[1]]$selection$candidates$weight,
    stable$models[[1]]$selection$candidates$weight,
    tolerance = 1e-9
  )
  stable$control$wbic_weight <- NULL
  stable$models[[1]]$control$wbic_weight <- NULL
  expect_equal(
    calc_model(stable, ci_method = "Param")$candidates$weight_method,
    rep("qr", 5)
  )
})
