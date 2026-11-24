test_that("undefined permutation statistics cannot reject the null", {
  d <- data.frame(year = 1:12, rate = 1)
  for (method in c("fixed", "tsprt")) {
    fit <- joinpoint(
      d,
      "year",
      "rate",
      max_jp = 1,
      selection_method = "Perm",
      perm_method = method,
      n_perm = 99,
      n_eq = 2
    )
    expect_equal(fit$model_info$n_joinpoints, 0L)
    expect_identical(fit$selection_tests$p_value, NA_real_)
    expect_identical(fit$selection_tests$reject, NA)
    expect_equal(fit$selection_tests$n_permutations, 0L)
    expect_equal(fit$selection_tests$test_status, "unavailable")
  }
})

test_that("scalar and batch grid SSE retain near-zero residuals", {
  for (weighted in c(FALSE, TRUE)) {
    d <- data.frame(trend_time = 1:12, trend_log = log(10) + 1e-8 * sin(1:12))
    if (weighted) {
      d$trend_weight <- seq(1, 2, length.out = 12)
    }
    candidate <- fit_candidate(1L, d, weighted, 2L, "permutation")
    response <- cbind(d$trend_log, d$trend_log + .01 * cos(1:12))
    expected <- apply(response, 2, function(y) {
      min(vapply(
        integer_joinpoint_grid(d$trend_time, 1L),
        function(tau) {
          x <- joinpoint_whiten(d, joinpoint_design(tau, d))
          sum(lm.fit(x, as.numeric(joinpoint_whiten(d, y)))$residuals^2)
        },
        numeric(1)
      ))
    })
    actual <- candidate_sse_batch(candidate, response, d, weighted)
    expect_equal(actual[1] / expected[1], 1, tolerance = 1e-5)
    expect_equal(actual[2], expected[2], tolerance = 1e-9)
    expect_equal(
      candidate_sse(candidate, response[, 1], d, weighted),
      actual[1]
    )
  }
})

test_that("batched grid projections preserve explicit QR minima", {
  data <- data.frame(
    trend_time = 2000:2011,
    trend_outcome = exp(
      seq(-2, -1.5, length.out = 12) + sin(seq_len(12)) / 100
    ),
    trend_se = rep(0.01, 12),
    trend_row_id = seq_len(12)
  )
  data <- prepare_group_error_model(data, "Ind", 0, NULL, "outcome")
  data$trend_log <- log(data$trend_outcome)
  grid <- integer_joinpoint_grid(data$trend_time, 2L)
  solver <- compile_joinpoint_grid_solver(grid, data, TRUE)
  responses <- cbind(
    data$trend_log,
    data$trend_log + sin(seq_len(nrow(data))) / 100
  )
  transformed <- joinpoint_whiten(data, responses)
  actual <- grid_projection_minima(solver, transformed)
  reference <- vapply(
    grid,
    function(joinpoints) {
      design <- joinpoint_whiten(data, joinpoint_design(joinpoints, data))
      fit <- qr(design)
      colSums((transformed - qr.fitted(fit, transformed))^2)
    },
    numeric(ncol(transformed))
  )
  reference <- t(reference)
  expected_index <- apply(reference, 2L, which.min)

  expect_equal(actual$sse, reference[cbind(expected_index, 1:2)])
  expect_equal(solver$model_indices[actual$model_index], expected_index)
})

test_that("within-group whitening matrices are cached by row subset", {
  data <- data.frame(
    trend_time = 2000:2005,
    trend_outcome = exp(seq(-2, -1.8, length.out = 6)),
    trend_se = rep(0.01, 6),
    trend_row_id = seq_len(6)
  )
  data <- prepare_group_error_model(data, "Ind", 0, NULL, "outcome")
  cache <- attr(data, "crtrend_whitener_cache", exact = TRUE)

  expect_length(ls(cache), 0L)
  full <- joinpoint_whitener(data)
  expect_length(ls(cache), 1L)
  expect_identical(joinpoint_whitener(data), full)
  joinpoint_whitener(data[1:4, , drop = FALSE])
  expect_length(ls(cache), 2L)
})

test_that("invalid resamples abort both permutation modes without discarding draws", {
  d <- data.frame(trend_time = 1:12, trend_log = sin(1:12))
  reduced <- fit_candidate(0L, d, FALSE, 2L, "permutation")
  complex <- fit_candidate(1L, d, FALSE, 2L, "permutation")
  local_mocked_bindings(
    candidate_sse = function(...) NA_real_,
    candidate_sse_batch = function(candidate, response, ...) {
      rep(NA_real_, ncol(response))
    }
  )
  for (method in c("fixed", "tsprt")) {
    result <- monte_carlo_p_value(
      reduced,
      complex,
      d,
      FALSE,
      19L,
      .05,
      method,
      .01,
      .01,
      NULL
    )
    expect_identical(result$p_value, NA_real_)
    expect_identical(result$reject, NA)
  }
})
