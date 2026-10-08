test_that("calc_model reproduces the selected NCI 6.1 model outputs", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    se = "se",
    max_jp = 4,
    n_eq = 9
  )
  result <- calc_model(fit, ci_method = "Param")
  expect_named(
    result,
    c(
      "data",
      "statistics",
      "coefficients",
      "vcov",
      "correlation",
      "candidates",
      "selection",
      "selection_tests",
      "joinpoints"
    )
  )
  expect_equal(result$statistics$n_parameters, 4L)
  expect_equal(result$statistics$model_df, 18L)
  expect_equal(result$statistics$fixed_grid_df, 19L)
  expect_equal(result$statistics$inference_df, 17L)
  expect_equal(result$statistics$effective_n, 21L)
  expect_equal(round(result$statistics$sse, 4), 587.7386)
  expect_equal(round(result$statistics$mse, 4), 32.6521)
  expect_equal(result$statistics$inference_status, "ok")
  expect_equal(result$statistics$candidates_available, TRUE)
  expect_equal(result$joinpoints, fit$joinpoints)

  general <- subset(result$coefficients, parameterization == "general")
  standard <- subset(result$coefficients, parameterization == "standard")
  expect_equal(
    round(general$estimate, 6),
    c(60.600790, -.028280, 32.963267, -.014544)
  )
  expect_equal(round(general$se, 6), c(3.802859, .001896, 6.071286, .003010))
  expect_equal(
    round(general$statistic, 6),
    c(15.935587, -14.913115, 5.429371, -4.831790)
  )
  expect_equal(round(general$p_value, 6), c(0, 0, .000045, .000156))
  expect_equal(round(standard$estimate, 6), c(60.600790, -.028280, .013736))
  expect_equal(round(standard$se, 6), c(3.802859, .001896, .003558))
  expect_equal(round(standard$statistic, 6), c(15.935587, -14.913115, 3.861130))
  expect_equal(round(standard$p_value, 6), c(0, 0, .001253))
  expect_equal(unique(result$coefficients$statistic_type), "t")
  expect_equal(unique(result$coefficients$df), 17L)
  covariance <- result$vcov[["1"]]
  expect_equal(dimnames(covariance), list(general$term, general$term))
  reference <- matrix(
    c(
      14.461736,
      -.007211,
      0,
      0,
      -.007211,
      .000004,
      0,
      0,
      0,
      0,
      36.860516,
      -.018275,
      0,
      0,
      -.018275,
      .000009
    ),
    4
  )
  expect_equal(unname(round(covariance, 6)), reference)
  expect_equal(result$correlation[["1"]], cov2cor(covariance))

  modeled <- c(
    56.84,
    55.26,
    53.71,
    52.22,
    50.76,
    49.35,
    47.97,
    46.63,
    45.33,
    44.07,
    42.84,
    41.64,
    40.48,
    39.90,
    39.32,
    38.75,
    38.19,
    37.64,
    37.10,
    36.56,
    36.04,
    35.52
  )
  expect_equal(round(result$data$fitted, 2), modeled)
  expect_equal(result$data$fitted, as.numeric(fitted(fit)))
  expect_equal(result$data$residual, as.numeric(residuals(fit)))
  expect_equal(result$data$observed, subset(jp_data, sex == 0)$rate)
  expect_equal(result$data$se, subset(jp_data, sex == 0)$se)
  expect_equal(result$data$joinpoint[13], 1L)
  expect_equal(result$data$apc[13], NA_real_)
  expect_equal(result$data$segment[13], NA_integer_)
  expect_equal(result$data$apc_significant[-13], rep(TRUE, 21))
  expect_equal(
    round(result$data$apc[-13], 4),
    c(rep(-2.7884, 12), rep(-1.4439, 9))
  )
  expect_equal(result$candidates$n_joinpoints, 0:4)
  expect_equal(result$candidates$model_df, c(20L, 18L, 16L, 14L, 12L))
  expect_equal(result$candidates$selected, c(FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(result$candidates$weight[1], NA_real_)
  expect_equal(result$candidates$weighted_bic[2], fit$model_info$weighted_bic)
  expect_equal(result$selection$selected_n_joinpoints, 1L)
  expect_equal(result$selection$selection_method, "WBIC")
  expect_equal(nrow(result$selection_tests), 0L)
  expect_true(all(is.na(result$candidates$dds_delta)))
})

test_that("extraction reuses caches and leaves the model and RNG unchanged", {
  set.seed(123)
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    se = "se",
    max_jp = 1,
    n_eq = 9
  )
  original <- fit
  rng <- .Random.seed
  local_mocked_bindings(
    fit_candidate = function(...) stop("must not search"),
    nci_joinpoint_inference = function(...) stop("must reuse covariance"),
    nci_empirical_model_draws = function(...) stop("must reuse EQ")
  )
  eq <- calc_model(fit)
  param <- calc_model(fit, ci_method = "Param")
  expect_identical(fit, original)
  expect_identical(.Random.seed, rng)
  expect_equal(eq$coefficients, param$coefficients)
  expect_equal(eq$candidates, param$candidates)
  apc <- calc_apc(fit, n_resamples = 9)
  rows <- !is.na(eq$data$segment)
  expect_equal(eq$data$apc_p_value[rows], apc$p_value[eq$data$segment[rows]])
  expect_equal(unique(eq$data$apc_method), "EQ")
})

test_that("grouped outputs preserve original row order and failed groups", {
  data <- data.frame(
    time = rep(1:8, 3),
    sex = rep(c("b", "a", "zero"), each = 8)
  )
  data$rate <- exp(.1 * data$time + sin(data$time) / 20)
  data$rate[data$sex == "zero"] <- 0
  data <- data[rev(seq_len(nrow(data))), ]
  fit <- joinpoint(
    data,
    "time",
    "rate",
    group_vars = "sex",
    max_jp = 0,
    n_eq = 2
  )
  result <- calc_model(fit, ci_method = "Param")
  expect_equal(result$data$row_id, seq_len(nrow(data)))
  expect_equal(result$data$sex, data$sex)
  expect_equal(result$data$time, data$time)
  expect_equal(result$data$observed, data$rate)
  expect_equal(result$data$se, rep(NA_real_, nrow(data)))
  expect_equal(result$data$fitted, as.numeric(fitted(fit)))
  failed <- subset(result$statistics, sex == "zero")
  expect_equal(failed$status, "zero_risk")
  expect_equal(failed$model_df, NA_integer_)
  expect_equal(failed$sse, NA_real_)
  expect_equal(failed$candidates_available, FALSE)
  expect_null(result$vcov[[as.character(failed$group_id)]])
  expect_equal(nrow(subset(result$coefficients, sex == "zero")), 0L)
  expect_equal(nrow(subset(result$candidates, sex == "zero")), 0L)
  expect_equal(subset(result$data, sex == "zero")$apc, rep(NA_real_, 8))

  all_failed <- joinpoint(subset(data, sex == "zero"), "time", "rate", n_eq = 2)
  empty <- calc_model(all_failed, ci_method = "Param")
  expect_equal(nrow(empty$coefficients), 0L)
  expect_equal(nrow(empty$candidates), 0L)
  expect_equal(length(empty$vcov), 1L)
  expect_null(empty$vcov[[1]])
})

test_that("general coefficients reconstruct multi-joinpoint predictions", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year_index",
    "rate",
    se = "se",
    min_jp = 3,
    max_jp = 3,
    n_eq = 2
  )
  result <- calc_model(fit, ci_method = "Param")
  general <- subset(result$coefficients, parameterization == "general")
  intercepts <- general$estimate[c(TRUE, FALSE)]
  slopes <- general$estimate[c(FALSE, TRUE)]
  segment <- findInterval(result$data$time, fit$models[[1]]$joinpoints) + 1L
  expected <- exp(intercepts[segment] + slopes[segment] * result$data$time)
  expect_equal(expected, result$data$fitted, tolerance = 1e-12)
  expect_equal(dim(result$vcov[[1]]), c(8L, 8L))
  expect_equal(result$candidates$n_joinpoints, 3L)
  expect_equal(result$candidates$selected, TRUE)
  apc <- calc_apc(fit, ci_method = "Param")
  expect_equal(general$statistic[c(FALSE, TRUE)], apc$statistic)
  expect_equal(general$p_value[c(FALSE, TRUE)], apc$p_value)
})

test_that("all selectors retain candidates and their actual selected flag", {
  for (method in c("WBIC", "BIC", "BIC3", "MBIC", "DDS", "Perm")) {
    fit <- joinpoint(
      subset(jp_data, sex == 0),
      "year_index",
      "rate",
      se = "se",
      max_jp = 1,
      selection_method = method,
      n_perm = 9,
      n_eq = 2
    )
    result <- calc_model(fit, ci_method = "Param")
    expect_equal(result$candidates$n_joinpoints, 0:1)
    expect_equal(
      result$candidates$n_joinpoints[result$candidates$selected],
      fit$models[[1]]$n_joinpoints
    )
    expect_equal(unique(result$candidates$selection_method), method)
    expect_equal(
      result$selection$selected_n_joinpoints,
      fit$models[[1]]$n_joinpoints
    )
    expect_equal(nrow(result$selection_tests), if (method == "Perm") 1L else 0L)
  }
})

test_that("calc_model exposes DDS branch diagnostics and per-candidate Delta", {
  x <- 0:30
  d <- data.frame(
    x = x,
    y = exp(
      4 -
        .02 * x +
        .02 * pmax(x - 10, 0) +
        .003 * pmax(x - 22, 0) +
        .01 * sin(2 * x)
    )
  )
  fit <- joinpoint(d, "x", "y", max_jp = 3, selection_method = "DDS", n_eq = 2)
  result <- calc_model(fit, ci_method = "Param")
  expect_equal(
    result$candidates$dds_delta,
    fit$models[[1]]$selection$candidates$dds_delta
  )
  expect_equal(
    result$candidates$dds_score,
    fit$models[[1]]$selection$candidates$dds_score
  )
  expect_equal(
    result$selection$dds_bic_delta,
    fit$models[[1]]$selection$dds_bic_delta
  )
  expect_equal(
    result$selection$dds_bic3_delta,
    fit$models[[1]]$selection$dds_bic3_delta
  )
  expect_equal(
    result$selection$dds_delta_max,
    fit$models[[1]]$selection$dds_delta_max
  )
  expect_equal(
    result$selection$dds_delta_min,
    fit$models[[1]]$selection$dds_delta_min
  )
  expect_equal(result$selection$dds_variance_df, 23L)
  expect_equal(result$selection$dds_decision_rule, "delta_range_gt_200")
  expect_equal(result$selection$selected_joinpoints, "10;21")
  expect_false(anyNA(result$candidates$dds_delta))
  expect_equal(
    result$candidates$dds_delta[result$candidates$selected],
    result$selection$dds_delta_min
  )
})

test_that("coefficient covariance respects correlated and known errors", {
  data <- data.frame(time = 0:15)
  data$rate <- exp(
    2 - .03 * data$time + .04 * pmax(data$time - 7, 0) + sin(data$time) / 20
  )
  covariance <- .0025 * .4^abs(outer(data$time, data$time, `-`))
  for (method in c("AR1", "GLS", "Survey")) {
    args <- list(
      data,
      "time",
      "rate",
      min_jp = 1,
      max_jp = 1,
      n_eq = 2,
      error = method
    )
    if (method == "AR1") {
      args$rho <- .4
    } else {
      args$cov <- covariance
      args$cov_scale <- "log"
    }
    fit <- do.call(joinpoint, args)
    result <- calc_model(fit, ci_method = "Param")
    model <- fit$models[[1]]
    tau <- model$joinpoints
    keep <- data$time != tau
    time <- data$time[keep]
    x <- cbind(time < tau, time * (time < tau), time > tau, time * (time > tau))
    v <- if (method == "AR1") covariance / .0025 else covariance
    precision <- solve(v[keep, keep])
    expected <- solve(crossprod(x, precision %*% x))
    if (method == "AR1") {
      expected <- expected * model$inference$residual_variance
    }
    expect_equal(unname(result$vcov[[1]]), unname(expected), tolerance = 1e-9)
    general <- subset(result$coefficients, parameterization == "general")
    expect_equal(general$se, sqrt(diag(expected)))
    expect_equal(result$statistics$error, method)
    apc <- calc_apc(fit, ci_method = "Param")
    expect_equal(general$statistic[c(2, 4)], apc$statistic)
  }
})

test_that("older fits rebuild coefficient inference but never candidate searches", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    se = "se",
    max_jp = 1,
    n_eq = 2
  )
  reference <- calc_model(fit, ci_method = "Param")
  fit$models[[1]]$selection$candidates <- NULL
  fit$models[[1]]$inference$coefficient_covariance <- NULL
  local_mocked_bindings(fit_candidate = function(...) stop("must not search"))
  expect_snapshot(result <- calc_model(fit, ci_method = "Param"))
  expect_equal(result$coefficients, reference$coefficients)
  expect_equal(result$vcov, reference$vcov)
  expect_equal(nrow(result$candidates), 0L)
  expect_equal(result$statistics$candidates_available, FALSE)
})

test_that("unavailable and zero-variance inference has explicit missing values", {
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year_index",
    "rate",
    min_jp = 1,
    max_jp = 1,
    n_eq = 2
  )
  zero <- fit
  zero$models[[1]]$coefficients[] <- c(0, -.01, 0)
  zero$models[[1]]$inference$coefficient_covariance[] <- 0
  result <- calc_model(zero, ci_method = "Param")
  expect_equal(result$coefficients$statistic[1:3], c(0, -Inf, 0))
  expect_equal(result$coefficients$p_value[1:3], c(1, 0, 1))
  expect_equal(unname(result$correlation[[1]]), matrix(NA_real_, 4, 4))

  fit$models[[1]]$inference <- NULL
  local_mocked_bindings(nci_joinpoint_inference = function(...) NULL)
  missing <- calc_model(fit, ci_method = "Param")
  expect_equal(unique(missing$coefficients$status), "inference_failed")
  expect_equal(missing$coefficients$se, rep(NA_real_, 7))
  expect_equal(missing$coefficients$statistic_type, rep(NA_character_, 7))
  expect_null(missing$vcov[[1]])
  expect_equal(missing$statistics$inference_status, "inference_failed")
  expect_equal(missing$data$apc, rep(NA_real_, 22))
})

test_that("calc_model validates inputs and supports method aliases", {
  expect_snapshot(error = TRUE, calc_model(data.frame()))
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    max_jp = 0,
    n_eq = 2
  )
  expect_snapshot(error = TRUE, calc_model(fit, level = 1))
  expect_snapshot(error = TRUE, calc_model(fit, n_resamples = 1))
  expect_snapshot(error = TRUE, calc_model(fit, ci_method = "invalid"))
  expect_equal(
    calc_model(fit, ci_method = "parametric"),
    calc_model(fit, ci_method = "Param")
  )
  expect_equal(
    calc_model(fit, ci_method = "empirical_quantile"),
    calc_model(fit)
  )
})
