test_that("DDS disagreement follows both official thresholds and their boundaries", {
  candidates <- list(
    list(n_joinpoints = 1L, bic = 0, bic3 = 2, sse = 20),
    list(n_joinpoints = 2L, bic = 1, bic3 = 0, sse = 14)
  )
  cases <- list(
    list(delta = c(3, 10), method = "BIC"),
    list(delta = c(3, 10 + 1e-6), method = "BIC3"),
    list(delta = c(3, 203), method = "BIC3"),
    list(delta = c(3, 203 + 1e-6), method = "BIC"),
    list(delta = c(11, 12), method = "BIC3"),
    list(delta = c(Inf, 3), method = "BIC")
  )
  for (case in cases) {
    delta <- case$delta
    local_mocked_bindings(
      dds_minimum_delta = function(candidate, ...) delta[candidate$n_joinpoints]
    )
    result <- select_dds_candidate(candidates, data.frame(x = 1:20), FALSE)
    expect_identical(result$dds_method, case$method)
    expected_k <- if (case$method == "BIC") 1L else 2L
    expect_identical(result$best$n_joinpoints, expected_k)
    expect_equal(result$best$dds_delta, delta[expected_k])
    expect_equal(result$best$dds_score, 0)
  }
  local_mocked_bindings(dds_minimum_delta = function(...) Inf)
  expect_error(
    select_dds_candidate(candidates, data.frame(x = 1:20), FALSE),
    "DDS cannot resolve BIC/BIC3 disagreement"
  )
  candidates[[1]]$bic3 <- -1
  expect_identical(
    select_dds_candidate(candidates, data.frame(x = 1:20), FALSE)$dds_method,
    "BIC=BIC3"
  )
})

test_that("DDS real disagreement uses maximum-model variance and local projection", {
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
  scores <- calc_model(fit, ci_method = "Param")$candidates
  expect_equal(scores$n_joinpoints[which.min(scores$bic)], 2L)
  expect_equal(scores$n_joinpoints[which.min(scores$bic3)], 1L)
  expect_identical(fit$model_info$dds_method, "BIC")
  expect_equal(fit$model_info$n_joinpoints, 2L)
  expect_equal(sum(scores$selected), 1)
  sigma2 <- scores$sse[scores$n_joinpoints == 3] / (31 - 8)
  # Independently calculate each adjacent-window projection with normal equations.
  data <- joinpoint_analysis_data(fit$models[[1]])
  candidate <- fit_candidate(2L, data, FALSE, 2L, "DDS")
  knots <- candidate$joinpoints
  changes <- coef(candidate$model)[3:4]
  endpoints <- c(min(data$trend_time) - 1, knots, max(data$trend_time))
  delta <- vapply(
    1:2,
    function(i) {
      keep <- data$trend_time > endpoints[i] &
        data$trend_time <= endpoints[i + 2]
      time <- data$trend_time[keep]
      design <- cbind(1, time)
      hinge <- pmax(time - knots[i], 0)
      residual <- hinge -
        design %*%
          solve(crossprod(design), crossprod(design, hinge))
      changes[i]^2 * sum(residual^2) / sigma2
    },
    numeric(1)
  )
  expect_equal(fit$model_info$dds_delta, min(delta), tolerance = 1e-8)
})

test_that("DDS can retain BIC3 despite BIC selecting an additional joinpoint", {
  x <- 0:30
  d <- data.frame(
    x = x,
    y = exp(
      4 -
        .02 * x +
        .005 * pmax(x - 10, 0) +
        .00325 * pmax(x - 22, 0) +
        .01 * sin(2 * x)
    )
  )
  fit <- joinpoint(d, "x", "y", max_jp = 3, selection_method = "DDS", n_eq = 2)
  scores <- calc_model(fit, ci_method = "Param")$candidates
  expect_equal(scores$n_joinpoints[which.min(scores$bic)], 2L)
  expect_equal(scores$n_joinpoints[which.min(scores$bic3)], 1L)
  expect_identical(fit$model_info$dds_method, "BIC3")
  expect_equal(fit$model_info$n_joinpoints, 1L)
  expect_equal(sum(scores$selected), 1)
})
