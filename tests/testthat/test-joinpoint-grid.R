test_that("grid SSE agrees with full scoring across covariance structures", {
  time <- c(2000:2009, 2011:2018)
  raw <- data.frame(
    year = time,
    rate = exp(
      2 + 0.03 * (time - 2000) - 0.08 * pmax(time - 2008, 0) + 0.01 * sin(time)
    ),
    se = seq(0.1, 0.2, length.out = length(time))
  )
  for (error in c("Ind", "AR1", "GLS", "Survey")) {
    data <- crtrend:::prepare_joinpoint_data(raw, "year", "rate", NULL, "se")
    covariance <- 0.002 * 0.3^abs(outer(seq_along(time), seq_along(time), `-`))
    data <- crtrend:::prepare_group_error_model(
      data,
      error,
      0.4,
      covariance,
      "log"
    )
    data$trend_log <- log(data$trend_outcome)
    for (k in 0:3) {
      grid <- if (k == 0L) {
        list(numeric())
      } else {
        crtrend:::integer_joinpoint_grid(time, k, 2L)
      }
      expected <- vapply(
        grid,
        function(points) {
          crtrend:::fixed_joinpoint_scores(
            points,
            data,
            TRUE,
            calculate_partial = FALSE
          )$sse
        },
        numeric(1)
      )
      actual <- crtrend:::joinpoint_grid_sse(grid, data)
      expect_equal(actual, expected, tolerance = 1e-12, info = error)
      expect_identical(which.min(actual), which.min(expected))

      solver <- crtrend:::compile_joinpoint_grid_solver(grid, data, TRUE)
      response <- crtrend:::joinpoint_whiten(data, data$trend_log)
      projected <- crtrend:::grid_projection_minima(solver, response)
      expect_equal(projected$sse, min(expected), tolerance = 1e-9)
      expect_equal(
        solver$model_indices[projected$model_index],
        which.min(expected)
      )
    }
  }
})

test_that("grid scoring handles singular and almost exact fits", {
  data <- data.frame(trend_time = 2000:2015, trend_log = 0)
  grid <- list(2000, 2007, 2015, c(2005, 2005))
  expect_equal(crtrend:::joinpoint_grid_sse(grid, data), c(Inf, 0, Inf, Inf))
  data$trend_log <- 2 +
    0.03 * (data$trend_time - 2000) -
    0.1 * pmax(data$trend_time - 2007, 0) +
    1e-10 * sin(data$trend_time)
  grid <- crtrend:::integer_joinpoint_grid(data$trend_time, 1L, 2L)
  expected <- vapply(
    grid,
    function(points) {
      crtrend:::fixed_joinpoint_scores(points, data, FALSE, FALSE)$sse
    },
    numeric(1)
  )
  actual <- crtrend:::joinpoint_grid_sse(grid, data)
  expect_equal(actual, expected, tolerance = 1e-14)
  expect_identical(which.min(actual), which.min(expected))
})
