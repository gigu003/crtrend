test_that("joinpoint returns grouped model fits and calculated results", {
  years <- 1990:2023
  data <- expand.grid(
    year = years,
    sex = c(1, 2),
    cancer = c("A", "B")
  )
  data$risk <- with(
    data,
    exp(-3 + 0.02 * (year - 1990) + 0.01 * (sex == 2))
  )

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = c("sex", "cancer"),
    max_jp = 1,
    n_perm = 9,
    selection_method = "Perm",
    seed = 1
  )

  expect_s3_class(result, "crtrend_joinpoint")
  expect_equal(nrow(result$model_info), 4L)
  expect_equal(length(result$models), 4L)
  expect_equal(result$control$selection_method, "Perm")
  expect_true(all(result$model_info$status == "ok"))
  expect_true(all(
    c(
      "mc_p_value",
      "mc_p_value_lower",
      "mc_p_value_upper",
      "mc_n_permutations",
      "mc_stopped_early"
    ) %in%
      names(result$model_info)
  ))
  expect_equal(
    result$model_info$mc_p_value,
    result$model_info$selection_p_value
  )
  expect_true(all(
    result$model_info$mc_p_value_lower <= result$model_info$mc_p_value
  ))
  expect_true(all(
    result$model_info$mc_p_value <= result$model_info$mc_p_value_upper
  ))
  apc <- calc_apc(result)
  aapc <- calc_aapc(result, start = 1990, end = 2023)
  expect_true(all(is.finite(apc$apc)))
  expect_true(all(is.finite(aapc$aapc)))
  expect_true(all(is.finite(aapc$p_value)))
  expect_true(all(is.finite(apc$p_value)))
  expect_true(all(c("sex", "cancer") %in% names(apc)))
  expect_true(all(
    c("period", "apc_lower", "apc_upper", "p_value") %in% names(apc)
  ))
  expect_true(all(
    c("aapc", "aapc_se", "aapc_lower", "aapc_upper", "p_value") %in%
      names(aapc)
  ))
  expect_true(all(is.finite(aapc$aapc_se)))
  expect_type(apc$start, "integer")
  expect_type(apc$end, "integer")
  expect_type(aapc$start, "integer")
  expect_type(aapc$end, "integer")
  expect_equal(result$control$perm_method, "fixed")
})

test_that("fixed candidate fits agree with the Joinpoint V5.3 reference", {
  reference <- new.env(parent = emptyenv())
  load(
    test_path("fixtures", "joinpoint-v5.3-fitted.rda"),
    envir = reference
  )
  expected <- reference$v5_3_fitted

  for (n_joinpoints in 0:4) {
    fit <- joinpoint(
      jp_data,
      time = "year",
      outcome = "rate",
      group_vars = "sex",
      se = "se",
      min_jp = n_joinpoints,
      max_jp = n_joinpoints,
      n_eq = 19L
    )
    observed <- data.frame(
      model = n_joinpoints,
      sex = jp_data$sex,
      year = jp_data$year,
      fitted = as.numeric(fitted(fit))
    )
    expected_model <- expected[expected$model == n_joinpoints, ]
    expected_model <- expected_model[
      order(expected_model$sex, expected_model$year),
    ]
    row.names(expected_model) <- NULL
    expect_equal(
      observed[, c("model", "sex", "year", "fitted")],
      expected_model[, c("model", "sex", "year", "fitted")],
      tolerance = 5e-10
    )
  }
})

test_that("standard Jump sample fits agree with Joinpoint V5.3", {
  reference <- new.env(parent = emptyenv())
  load(
    test_path("fixtures", "joinpoint-v5.3-fitted.rda"),
    envir = reference
  )
  expected <- reference$v5_3_jump_standard

  for (n_joinpoints in 0:5) {
    fit <- joinpoint(
      jp_jump,
      time = "year",
      outcome = "rate",
      se = "se",
      min_jp = n_joinpoints,
      max_jp = n_joinpoints,
      n_eq = 19L
    )
    observed <- data.frame(
      model = n_joinpoints,
      year = jp_jump$year,
      fitted = as.numeric(fitted(fit))
    )
    expected_model <- expected[expected$model == n_joinpoints, ]
    row.names(expected_model) <- NULL
    expect_equal(
      observed,
      expected_model[, c("model", "year", "fitted")],
      tolerance = 5e-10
    )
  }
})

test_that("joinpoint follows standard grouped model-object conventions", {
  data <- data.frame(
    year = rep(2000:2007, 2),
    sex = rep(c("male", "female"), each = 8)
  )
  data$risk <- exp(
    -2 +
      0.03 * (data$year - 2000) +
      0.2 * (data$sex == "female") +
      rep(c(0.01, -0.01), 8)
  )
  data <- data[c(9:16, 1:8), , drop = FALSE]

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = "sex",
    max_jp = 0,
    n_eq = 19L
  )

  expect_equal(nobs(result), nrow(data))
  expect_equal(deparse(formula(result)), "log(risk) ~ year")
  expect_equal(model.frame(result), data[c("sex", "year", "risk")])
  expect_equal(
    fitted(result) + residuals(result),
    unname(data$risk),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    fitted(result, type = "link") + residuals(result, type = "link"),
    unname(log(data$risk)),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    predict(result, newdata = data),
    unname(fitted(result)),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_named(
    coef(result),
    c(
      "group_id",
      "sex",
      "term",
      "estimate",
      "component",
      "joinpoint"
    )
  )
  expect_equal(unique(coef(result)$term), c("(Intercept)", "year"))
  expect_s3_class(summary(result), "summary.crtrend_joinpoint")
  expect_equal(summary(result)$nobs, nrow(data))
  expect_equal(is.null(result$fits), TRUE)
  expect_s3_class(result$models[[1L]], "crtrend_joinpoint_model")
  expect_equal(
    fitted(result$models[[1L]]) + residuals(result$models[[1L]]),
    result$models[[1L]]$model$risk,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(nobs(result$models[[1L]]), 8L)
  expect_equal(
    predict(
      result$models[[1L]],
      newdata = result$models[[1L]]$model
    ),
    fitted(result$models[[1L]]),
    tolerance = 1e-12
  )
  expect_s3_class(
    summary(result$models[[1L]]),
    "summary.crtrend_joinpoint_model"
  )
  expect_named(
    result$control,
    c(
      "group_vars",
      "time",
      "outcome",
      "se",
      "error",
      "rho",
      "cov_scale",
      "min_jp",
      "max_jp",
      "min_obs",
      "selection_method",
      "wbic_weight",
      "n_perm",
      "alpha",
      "perm_method",
      "perm_alpha0",
      "perm_beta0",
      "perm_p1",
      "perm_conf",
      "n_eq",
      "eq_seed"
    )
  )
  expect_named(
    result$segments,
    c("group_id", "sex", "segment", "start", "end", "slope")
  )
  expect_named(
    result$joinpoints,
    c("group_id", "sex", "index", "joinpoint")
  )
  expect_equal("joinpoints" %in% names(result$model_info), FALSE)
})

test_that("joinpoint normalizes data-frame subclasses before setting row names", {
  data <- data.frame(
    year = rep(2000:2007, 2),
    sex = rep(c("male", "female"), each = 8),
    risk = exp(rep(seq(0, 0.3, length.out = 8), 2))
  )
  class(data) <- c("tbl_df", "tbl", "data.frame")

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = "sex",
    max_jp = 0,
    n_eq = 19L
  )

  expect_identical(class(result$model), "data.frame")
  expect_identical(class(result$models[[1L]]$model), "data.frame")
})

test_that("predict returns NA for groups without a successful fit", {
  data <- data.frame(
    year = rep(2000:2005, 2),
    sex = rep(c("zero", "positive"), each = 6),
    risk = c(rep(0, 6), exp(seq(0, 0.5, length.out = 6)))
  )
  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = "sex",
    max_jp = 0,
    n_eq = 19L
  )

  prediction <- predict(result, newdata = data)
  expect_equal(
    unname(is.na(prediction[data$sex == "zero"])),
    rep(TRUE, 6)
  )
  expect_equal(
    unname(is.finite(prediction[data$sex == "positive"])),
    rep(TRUE, 6)
  )
})

test_that("WBIC is the default model selection label", {
  data <- data.frame(year = 2000:2010, risk = exp(seq(0, 0.5, length.out = 11)))

  result <- joinpoint(data, time = "year", outcome = "risk", max_jp = 1)

  expect_equal(result$control$selection_method, "WBIC")
  expect_true(is.na(result$model_info$selection_p_value))
})

test_that("joinpoint caches empirical quantile resamples", {
  set.seed(21)
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(0.02 * (data$year - 2000) + rnorm(11, sd = 0.02))

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    n_eq = 19L,
    selection_method = "BIC"
  )

  expect_equal(result$control$n_eq, 19L)
  expect_equal(result$control$eq_seed, 10000L)
  expect_equal(result$models[[1L]]$resamples$n, 19L)
  expect_equal(is.null(result$models[[1L]]$resamples$draws), FALSE)
  expect_equal(is.null(result$models[[1L]]$inference), FALSE)
  apc <- calc_apc(result, n_resamples = 19L)
  expect_true(all(is.finite(apc$apc_lower)))
})

test_that("joinpoint isolates the empirical quantile seed", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(0.02 * (data$year - 2000) + sin(data$year) / 100)

  set.seed(1)
  first <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    n_eq = 19L,
    eq_seed = 123
  )
  set.seed(999)
  second <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    n_eq = 19L,
    eq_seed = 123
  )
  third <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    n_eq = 19L,
    eq_seed = 124
  )

  expect_equal(first$control$eq_seed, 123)
  expect_equal(
    first$models[[1L]]$resamples$draws,
    second$models[[1L]]$resamples$draws
  )
  expect_false(
    isTRUE(all.equal(
      first$models[[1L]]$resamples$draws,
      third$models[[1L]]$resamples$draws
    ))
  )
})

test_that("partial R-squared matches weighted residualization", {
  data <- data.frame(
    trend_time = 2000:2010,
    trend_log = c(0.1, 0.3, 0.2, 0.5, 0.7, 0.8, 1.2, 1.1, 1.4, 1.7, 1.9),
    trend_weight = seq(0.8, 1.8, length.out = 11)
  )
  reference <- function(weighted) {
    keep <- data$trend_time > 1999 & data$trend_time <= 2010
    x <- data$trend_time[keep]
    y <- data$trend_log[keep]
    z <- pmax(x - 2005, 0)
    root_weight <- if (weighted) sqrt(data$trend_weight[keep]) else rep(1, 11)
    design <- cbind(1, x) * root_weight
    solver <- qr(design)
    y_residual <- y * root_weight - qr.fitted(solver, y * root_weight)
    z_residual <- z * root_weight - qr.fitted(solver, z * root_weight)
    sum(y_residual * z_residual)^2 /
      (sum(y_residual^2) * sum(z_residual^2))
  }

  expect_equal(
    crtrend:::partial_joinpoint_r_squared(1L, 2005, data, FALSE),
    reference(FALSE),
    tolerance = 1e-12
  )
  expect_equal(
    crtrend:::partial_joinpoint_r_squared(1L, 2005, data, TRUE),
    reference(TRUE),
    tolerance = 1e-10
  )
})

test_that("joinpoint uses integer year locations by default", {
  data <- data.frame(year = 1990:2023)
  data$risk <- exp(ifelse(
    data$year <= 2003,
    -3 + 0.02 * (data$year - 1990),
    -3 + 0.02 * 13 + 0.08 * (data$year - 2003)
  ))

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    selection_method = "weighted_bic"
  )

  expect_equal(result$model_info$n_joinpoints, 1L)
  expect_equal(result$joinpoints$joinpoint, 2003)
  apc <- calc_apc(result)
  expect_equal(apc$period, c("1990-2003", "2003-2023"))
  expect_true(all(apc$start == round(apc$start)))
  expect_true(all(apc$end == round(apc$end)))
})

test_that("integer grid enforces observations between joinpoints", {
  data <- data.frame(year = 2000:2003, risk = exp(seq(0, 0.3, length.out = 4)))

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    min_obs = 2,
    selection_method = "weighted_bic"
  )

  expect_equal(result$model_info$n_joinpoints, 0L)
})

test_that("joinpoint rejects non-integer years", {
  data <- data.frame(year = 1990:2023)
  data$risk <- exp(
    -3 +
      0.02 * (data$year - 1990) +
      0.06 * pmax(data$year - 2003.5, 0) +
      0.001 * sin(data$year)
  )
  data$year <- data$year + 0.5

  error <- tryCatch(
    joinpoint(
      data,
      time = "year",
      outcome = "risk",
      max_jp = 1,
      selection_method = "weighted_bic"
    ),
    error = identity
  )

  expect_s3_class(error, "error")
  expect_match(error$message, "integer-valued years")
})

test_that("joinpoint supports standard-error weighted regression", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(
    1 +
      0.03 * (data$year - 2000) +
      c(0, 0.06, -0.03, 0.04, -0.02, 0.05, -0.04, 0.03, -0.02, 0.01, -0.03)
  )
  data$se <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 0,
    se = "se"
  )
  model <- result$models[[1L]]

  expect_equal(result$control$se, "se")
  expect_equal(model$weights, (data$risk / data$se)^2)
  expect_true(is.finite(calc_apc(result)$apc))
})

test_that("AR1 with zero correlation agrees with independent regression", {
  data <- data.frame(
    year = 2000:2010,
    risk = exp(1 + 0.03 * (0:10) + c(0.02, -0.01, rep(0, 9)))
  )
  independent <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 0,
    n_eq = 9L
  )
  ar1 <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 0,
    error = "AR1",
    rho = 0,
    n_eq = 9L
  )

  expect_equal(coef(ar1), coef(independent), tolerance = 1e-12)
  expect_equal(ar1$model_info$error, "AR1")
  expect_equal(ar1$models[[1L]]$control$rho, 0)
})

test_that("GLS uses the supplied full log-scale covariance", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(0.5 + 0.04 * (0:10) + sin(0:10) / 50)
  covariance <- 0.0025 * 0.45^abs(outer(0:10, 0:10, `-`))
  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 0,
    error = "GLS",
    cov = covariance,
    cov_scale = "log",
    n_eq = 9L
  )
  design <- cbind(`(Intercept)` = 1, trend_time = data$year)
  expected <- solve(
    crossprod(design, solve(covariance, design)),
    crossprod(design, solve(covariance, log(data$risk)))
  )

  expect_equal(
    unname(coef(result$models[[1L]])),
    as.numeric(expected),
    tolerance = 1e-10
  )
  expect_equal(result$control$error, "GLS")
})

test_that("Survey accepts outcome-scale design covariance", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(0.2 + 0.02 * (0:10))
  log_covariance <- 0.0016 * 0.3^abs(outer(0:10, 0:10, `-`))
  outcome_covariance <- log_covariance * tcrossprod(data$risk)
  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 0,
    error = "Survey",
    cov = outcome_covariance,
    cov_scale = "outcome",
    n_eq = 9L
  )
  design <- cbind(1, data$year)
  expected_covariance <- solve(crossprod(design, solve(log_covariance, design)))

  expect_equal(
    result$models[[1L]]$inference$slope_covariance,
    expected_covariance[2, 2, drop = FALSE],
    tolerance = 1e-8
  )
  expect_equal(result$model_info$error, "Survey")
})

test_that("general covariance inputs are validated", {
  data <- data.frame(year = 2000:2005, risk = exp(1 + 0.02 * (0:5)))
  expect_error(
    joinpoint(data, "year", "risk", error = "GLS"),
    "`cov` is required"
  )
  expect_error(
    joinpoint(
      data,
      "year",
      "risk",
      error = "GLS",
      cov = matrix(1, 6, 6),
      n_eq = 9L
    ),
    "positive definite"
  )
})

test_that("AAPC uses the APC t interval within one segment", {
  data <- data.frame(
    year = 2000:2010,
    risk = exp(0.02 * (2000:2010 - 2000))
  )
  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 0,
    selection_method = "weighted_bic"
  )
  apc <- calc_apc(result, level = 0.95, ci_method = "Param")
  aapc <- calc_aapc(
    result,
    start = 2002,
    end = 2008,
    level = 0.95,
    ci_method = "Param"
  )

  expect_equal(aapc$aapc_lower, apc$apc_lower)
  expect_equal(aapc$aapc_upper, apc$apc_upper)
  expect_equal(aapc$p_value, apc$p_value)
  expect_equal(aapc$statistic, apc$statistic)
  expect_equal(aapc$statistic_type, "t")
})

test_that("NCI inference deletes joinpoint observations for segment errors", {
  set.seed(11)
  data <- data.frame(year = 1990:2023)
  data$risk <- exp(
    ifelse(
      data$year <= 2003,
      -3 + 0.02 * (data$year - 1990),
      -3 + 0.02 * 13 + 0.08 * (data$year - 2003)
    ) +
      rnorm(nrow(data), sd = 0.05)
  )
  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    selection_method = "BIC"
  )
  expect_equal(result$model_info$n_joinpoints, 1L)
  inference <- crtrend:::nci_joinpoint_inference(result$models[[1L]])
  expect_equal(inference$effective_n, nrow(data) - 1L)
  expect_equal(inference$df, nrow(data) - 3L - 2L)
  expect_length(inference$slope_estimates, 2L)
  expect_true(all(is.finite(inference$slope_se)))
  apc <- calc_apc(result, ci_method = "Param")
  expect_equal(
    apc$apc,
    unname(100 * (exp(inference$slope_estimates) - 1)),
    tolerance = 1e-12
  )
})

test_that("Empirical quantile intervals refit the selected joinpoint count", {
  set.seed(12)
  data <- data.frame(year = 1990:2023)
  data$risk <- exp(
    ifelse(
      data$year <= 2003,
      -3 + 0.02 * (data$year - 1990),
      -3 + 0.02 * 13 + 0.08 * (data$year - 2003)
    ) +
      rnorm(nrow(data), sd = 0.05)
  )
  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    selection_method = "BIC"
  )
  apc <- calc_apc(result, n_resamples = 99, seed = 1)
  aapc <- calc_aapc(
    result,
    start = 1990,
    end = 2023,
    n_resamples = 99,
    seed = 1
  )
  expect_true(all(is.finite(apc$apc_lower)))
  expect_true(all(is.finite(apc$apc_upper)))
  expect_true(all(is.finite(apc$p_value)))
  expect_true(is.finite(aapc$aapc_lower))
  expect_true(is.finite(aapc$aapc_upper))
  expect_true(is.finite(aapc$p_value))
})

test_that("joinpoint can use BIC model selection", {
  data <- data.frame(year = 2000:2010, risk = exp(seq(0, 0.5, length.out = 11)))

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    selection_method = "weighted_bic"
  )

  expect_equal(result$control$selection_method, "WBIC")
  expect_true(is.na(result$model_info$selection_p_value))
})

test_that("joinpoint honors the minimum number of joinpoints", {
  data <- data.frame(year = 2000:2010, risk = exp(seq(0, 0.5, length.out = 11)))

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    min_jp = 1,
    max_jp = 1,
    selection_method = "weighted_bic"
  )

  expect_equal(result$control$min_jp, 1L)
  expect_equal(result$model_info$n_joinpoints, 1L)
})

test_that("joinpoint exposes separate BIC-family scores", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(ifelse(
    data$year <= 2005,
    0.02 * (data$year - 2000),
    0.1 + 0.08 * (data$year - 2005)
  ))

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    selection_method = "weighted_bic"
  )

  expect_true(all(
    c("bic", "bic3", "weighted_bic") %in%
      names(result$model_info)
  ))
  expect_true(is.finite(result$model_info$weighted_bic))
  expect_true(is.finite(result$model_info$partial_r_squared))
})

test_that("WBIC scores the least-SSE model for each joinpoint count", {
  data <- data.frame(
    trend_time = 2000:2020,
    trend_log = c(
      -0.00113875629863196,
      -0.0491486900977402,
      -0.0465104152391951,
      0.16262167050454,
      0.0297607713878162,
      -0.0974979665109614,
      0.136641391811953,
      0.208966200924805,
      -0.03350693677544,
      0.147294128369625,
      0.262939459069975,
      0.0764724896046638,
      -0.0163645958421821,
      0.134928758820908,
      0.186343452286635,
      0.287693476721379,
      0.203795410466745,
      0.196011937276681,
      0.252529974727508,
      0.296150049329679,
      0.222403340228658
    )
  )
  data$trend_outcome <- exp(data$trend_log)
  grid <- crtrend:::integer_joinpoint_grid(data$trend_time, 2L, 2L)
  scores <- lapply(
    grid,
    crtrend:::fixed_joinpoint_scores,
    data = data,
    weighted = FALSE
  )
  sse_best <- grid[[which.min(vapply(scores, `[[`, numeric(1), "sse"))]]
  wbic_best <- grid[[which.min(vapply(
    scores,
    `[[`,
    numeric(1),
    "weighted_bic"
  ))]]

  candidate <- crtrend:::fit_candidate(
    2L,
    data,
    weighted = FALSE,
    min_observations = 2L,
    selection_method = "weighted_bic"
  )

  expect_equal(candidate$joinpoints, sse_best)
  expect_equal(sse_best, c(2009, 2012))
  expect_equal(wbic_best, c(2007, 2011))
})

test_that("permutation SSE refits joinpoint locations over the full grid", {
  year <- 2000:2014
  data <- data.frame(trend_time = year)
  data$trend_log <- 0.01 *
    (year - 2000) +
    0.12 * pmax(year - 2005, 0) +
    0.003 * sin(year)
  data$trend_outcome <- exp(data$trend_log)
  candidate <- crtrend:::fit_candidate(
    1L,
    data,
    weighted = FALSE,
    min_observations = 2L,
    selection_method = "permutation"
  )
  permuted_log <- 0.01 *
    (year - 2000) +
    0.15 * pmax(year - 2009, 0) +
    0.004 * cos(year)
  permuted_data <- data
  permuted_data$trend_log <- permuted_log
  grid <- crtrend:::integer_joinpoint_grid(year, 1L, 2L)
  grid_sse <- vapply(
    grid,
    function(joinpoints) {
      crtrend:::fixed_joinpoint_scores(
        joinpoints,
        permuted_data,
        weighted = FALSE
      )$sse
    },
    numeric(1)
  )
  fixed_sse <- crtrend:::fixed_joinpoint_scores(
    candidate$joinpoints,
    permuted_data,
    weighted = FALSE
  )$sse

  expect_equal(candidate$joinpoints, 2005)
  expect_equal(grid[[which.min(grid_sse)]], 2009)
  expect_equal(
    crtrend:::candidate_sse(
      candidate,
      permuted_log,
      data,
      weighted = FALSE
    ),
    min(grid_sse),
    tolerance = 1e-10
  )
  expect_equal(
    crtrend:::candidate_sse_batch(
      candidate,
      cbind(permuted_log, rev(permuted_log)),
      data,
      weighted = FALSE
    ),
    c(
      crtrend:::candidate_sse(
        candidate,
        permuted_log,
        data,
        weighted = FALSE
      ),
      crtrend:::candidate_sse(
        candidate,
        rev(permuted_log),
        data,
        weighted = FALSE
      )
    ),
    tolerance = 1e-10
  )
  expect_gt(fixed_sse, min(grid_sse))
})

test_that("permutation selection uses the configured minimum joinpoints", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(seq(0, 0.5, length.out = 11))

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    min_jp = 1,
    max_jp = 2,
    selection_method = "Perm",
    n_perm = 9,
    seed = 1
  )

  expect_equal(result$model_info$n_joinpoints, 1L)
  expect_true(is.finite(result$model_info$selection_p_value))
})

test_that("permutation selection records every sequential comparison", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(
    0.02 *
      (data$year - 2000) +
      0.08 * pmax(data$year - 2005, 0) +
      0.002 * sin(data$year)
  )

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 4,
    selection_method = "Perm",
    n_perm = 19,
    n_eq = 19,
    seed = 11
  )

  expect_equal(result$model_info$effective_max_jp, 3L)
  expect_equal(nrow(result$selection_tests), 3L)
  expect_equal(result$selection_tests$step, 1:3)
  expect_equal(result$selection_tests$k_null[1L], 0L)
  expect_equal(result$selection_tests$k_alt[1L], 3L)
  expect_equal(result$selection_tests$alpha[1L], 0.05 / 3)
  expect_equal(
    tail(result$selection_tests$p_value, 1L),
    result$model_info$selection_p_value
  )
  expect_equal(
    result$models[[1L]]$selection$tests,
    result$selection_tests[names(result$models[[1L]]$selection$tests)]
  )
})

test_that("Fay tSPRT can stop before the maximum and fixed MC is available", {
  data <- data.frame(year = 2000:2010)
  data$risk <- exp(ifelse(
    data$year <= 2005,
    0.02 * (data$year - 2000),
    0.1 + 0.08 * (data$year - 2005)
  ))

  tsprt <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    n_perm = 99,
    selection_method = "Perm",
    perm_method = "tsprt",
    seed = 1
  )
  fixed <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    max_jp = 1,
    n_perm = 99,
    selection_method = "Perm",
    perm_method = "fixed",
    seed = 1
  )

  expect_true(tsprt$model_info$selection_n_permutations <= 99L)
  expect_true(is.logical(tsprt$model_info$selection_stopped_early))
  expect_true(is.finite(tsprt$model_info$selection_p_value))
  expect_true(is.finite(tsprt$model_info$mc_p_value_lower))
  expect_true(is.finite(tsprt$model_info$mc_p_value_upper))
  expect_equal(fixed$model_info$selection_n_permutations, 99L)
  expect_false(fixed$model_info$selection_stopped_early)
  expect_equal(fixed$model_info$mc_n_permutations, 99L)
  expect_true(
    fixed$model_info$mc_p_value_lower <= fixed$model_info$mc_p_value &&
      fixed$model_info$mc_p_value <= fixed$model_info$mc_p_value_upper
  )
})

test_that("joinpoint preserves zero-risk group status", {
  data <- data.frame(
    year = rep(2000:2005, 2),
    sex = rep(c(1, 2), each = 6),
    cancer = "A",
    risk = c(rep(0, 6), exp(seq(0, 0.5, length.out = 6)))
  )

  result <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = "sex",
    max_jp = 0
  )

  expect_equal(
    result$model_info$status[match(1, result$model_info$sex)],
    "zero_risk"
  )
  expect_equal(
    result$model_info$status[match(2, result$model_info$sex)],
    "ok"
  )
  aapc <- calc_aapc(result, start = 2000, end = 2005)
  expect_true(is.na(aapc$aapc[aapc$sex == 1]))
})

test_that("asr_inci is available as a demonstration dataset", {
  data(asr_inci, package = "crtrend")

  expect_s3_class(asr_inci, "data.frame")
  expect_equal(dim(asr_inci), c(2822L, 14L))
  expect_true(all(c("year", "sex", "cancer", "asr_wld85") %in% names(asr_inci)))
})

test_that("NCI unconstrained inference uses one pooled residual variance", {
  year <- 2000:2020
  outcome <- exp(
    1 +
      0.01 * (year - 2000) +
      0.04 * pmax(year - 2010, 0) +
      c(rep(0.01, 10), 0, rep(-0.01, 10))
  )
  result <- joinpoint(
    data.frame(year, outcome),
    "year",
    "outcome",
    min_jp = 1,
    max_jp = 1,
    selection_method = "BIC",
    n_eq = 9
  )
  fit <- result$models[[1L]]
  inference <- fit$inference
  analysis <- crtrend:::joinpoint_analysis_data(fit)
  keep <- analysis$trend_time != fit$joinpoints
  analysis <- analysis[keep, , drop = FALSE]
  segment <- findInterval(analysis$trend_time, fit$joinpoints) + 1L
  design <- matrix(0, nrow(analysis), 4L)
  for (index in 1:2) {
    rows <- segment == index
    design[rows, 2L * index - 1L] <- 1
    design[rows, 2L * index] <- analysis$trend_time[rows]
  }
  residual <- qr.resid(qr(design), analysis$trend_log)
  expected_variance <- sum(residual^2) / (nrow(analysis) - ncol(design))

  expect_equal(inference$df, nrow(analysis) - ncol(design))
  expect_equal(inference$residual_variance, expected_variance)
  expect_equal(
    inference$slope_covariance,
    expected_variance * solve(crossprod(design))[c(2L, 4L), c(2L, 4L)]
  )
})

test_that("AAPC parametric variance uses the full slope covariance", {
  year <- 2000:2020
  outcome <- exp(1 + 0.01 * (year - 2000) + 0.04 * pmax(year - 2010, 0))
  result <- joinpoint(
    data.frame(year, outcome),
    "year",
    "outcome",
    min_jp = 1,
    max_jp = 1,
    selection_method = "BIC",
    n_eq = 9
  )
  covariance <- matrix(c(0.01, 0.004, 0.004, 0.04), 2L, 2L)
  result$models[[1L]]$inference$slope_covariance <- covariance
  result$models[[1L]]$inference$slope_se <- sqrt(diag(covariance))
  aapc <- calc_aapc(result, 2000, 2020, ci_method = "Param")
  weights <- c(0.5, 0.5)
  log_estimate <- sum(weights * result$models[[1L]]$inference$slope_estimates)
  expected_se <- 100 *
    exp(log_estimate) *
    sqrt(drop(crossprod(weights, covariance %*% weights)))

  expect_equal(aapc$aapc_se, expected_se)
})

test_that("EmpQ2 endpoints cover both MAD and observed-extreme branches", {
  expect_equal(crtrend:::nci_empirical_endpoints(-2:2), c(-5.189, 5.189))
  expect_equal(
    crtrend:::nci_empirical_endpoints(c(-100, -1, 0, 1, 100)),
    c(-100.8, 100.8)
  )
  expect_equal(crtrend:::nci_empirical_residuals(rep(0, 5), 20L), rep(0, 20))
})

test_that("Empirical quantile p-values use the documented two-sided sign rule", {
  draws <- matrix(
    c(-2, -1, 0, 1, 2, 3, -1, 0, 1, 2),
    nrow = 5L,
    ncol = 2L
  )

  expect_equal(
    crtrend:::empirical_two_sided_p(draws),
    c(0.8, 0.4)
  )
  expect_equal(
    crtrend:::empirical_two_sided_p(c(-1, -0.5, 0, 1, 2)),
    0.8
  )
  expect_equal(
    crtrend:::empirical_two_sided_p(c(-1, NA_real_, 1)),
    NA_real_
  )
})

test_that("NCI 6.1 audit exports are identified as parametric inference", {
  cases <- list(
    list(
      file = "02-heteroscedastic-whitening.txt",
      apc_lower = c(-3.7818, 1.0543),
      apc_upper = c(-1.2795, 1.8383),
      apc_p = c(0.001423, 0.000015),
      aapc_lower = -0.5789,
      aapc_upper = 0.3697,
      aapc_p = 0.662293
    ),
    list(
      file = "03-endpoint-outlier.txt",
      apc_lower = c(-6.3812, 2.6094),
      apc_upper = c(-2.7060, 3.1388),
      apc_p = c(0.001555, 0.000001),
      aapc_lower = 0.8093,
      aapc_upper = 1.5380,
      aapc_p = 0.000000
    )
  )
  audit_dir <- system.file(
    "extdata",
    "joinpoint-6.1.0-eq-audit",
    package = "crtrend"
  )

  for (case in cases) {
    data <- utils::read.delim(
      file.path(audit_dir, case$file),
      check.names = FALSE
    )
    fit <- joinpoint(
      data,
      "Year",
      "Rate",
      se = "Standard Error",
      min_jp = 1,
      max_jp = 1,
      selection_method = "WBIC",
      n_eq = 9
    )
    apc <- calc_apc(fit, ci_method = "Param")
    aapc <- calc_aapc(
      fit,
      min(data$Year),
      max(data$Year),
      ci_method = "Param"
    )

    expect_equal(round(apc$apc_lower, 4L), case$apc_lower)
    expect_equal(round(apc$apc_upper, 4L), case$apc_upper)
    expect_equal(round(apc$p_value, 6L), case$apc_p)
    expect_equal(round(aapc$aapc_lower, 4L), case$aapc_lower)
    expect_equal(round(aapc$aapc_upper, 4L), case$aapc_upper)
    expect_equal(round(aapc$p_value, 6L), case$aapc_p)
  }
})

test_that("MBIC and DDS expose NCI model-selection diagnostics", {
  year <- 2000:2020
  outcome <- exp(1 + 0.01 * (year - 2000) + 0.04 * pmax(year - 2010, 0))
  data <- data.frame(year, outcome)
  mbic <- joinpoint(
    data,
    "year",
    "outcome",
    max_jp = 2,
    selection_method = "MBIC",
    n_eq = 9
  )
  dds <- joinpoint(
    data,
    "year",
    "outcome",
    max_jp = 2,
    selection_method = "DDS",
    n_eq = 9
  )

  expect_equal(mbic$control$selection_method, "MBIC")
  expect_equal(dds$control$selection_method, "DDS")
  expect_equal(is.finite(mbic$model_info$modified_bic), TRUE)
  expect_equal(is.finite(dds$model_info$dds_delta), TRUE)
  expect_match(dds$model_info$dds_method, "BIC")
})

test_that("permutation comparisons use the F-type statistic", {
  statistic <- crtrend:::permutation_f_statistic(100, 70, 0L, 2L, 30L)
  expected <- ((100 - 70) / 4) / (70 / 24)

  expect_equal(statistic, expected)
})
