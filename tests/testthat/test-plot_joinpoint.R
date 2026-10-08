test_that("plot_joinpoint returns editable segmented trend layers", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(
    year = rep(2000:2011, 2),
    sex = rep(c("Male", "Female"), each = 12)
  )
  data$risk <- with(
    data,
    exp(
      -2 +
        0.02 * (year - 2000) +
        0.06 * pmax(year - 2005, 0) +
        0.15 * (sex == "Female")
    )
  )
  fit <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = "sex",
    min_jp = 1,
    max_jp = 1,
    n_eq = 19L
  )

  plot <- plot_joinpoint(
    fit,
    show_apc = TRUE,
    show_aapc = TRUE,
    n_resamples = 19L
  )

  expect_s3_class(plot, "ggplot")
  expect_gte(length(plot$layers), 5L)
  labels <- unlist(lapply(plot$layers, function(layer) layer$data$.label))
  expect_match(labels, "APC =", all = FALSE)
  expect_match(labels, "AAPC =", all = FALSE)
  expect_match(labels, "p [=<]", all = FALSE)
  built <- ggplot2::ggplot_build(plot)
  expect_length(unique(built$data[[2L]]$colour), 2L)
  expect_equal(plot$layers[[1L]]$geom_params$pointsize, NULL)
  expect_equal(plot$layers[[1L]]$aes_params$size, 1.2)
  expect_equal(plot$layers[[2L]]$aes_params$linewidth, 0.6)
  expect_equal(plot$scales$get_scales("colour")$guide, "none")
  apc_layer <- Filter(
    function(layer) ".vjust" %in% names(layer$data),
    plot$layers
  )
  expect_length(apc_layer, 1L)
  expect_setequal(unique(apc_layer[[1L]]$data$.vjust), c(-0.6, 1.5))
  expect_equal(apc_layer[[1L]]$aes_params$size, 3.0)
  aapc_layer <- Filter(
    function(layer) {
      ".label" %in% names(layer$data) && any(grepl("^AAPC", layer$data$.label))
    },
    plot$layers
  )
  expect_length(aapc_layer, 1L)
  expect_equal(aapc_layer[[1L]]$aes_params$size, 3.0)
})

test_that("plot_joinpoint filters groups and accepts a custom AAPC range", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(
    year = as.character(rep(2000:2010, 2)),
    sex = rep(c("Male", "Female"), each = 11),
    risk = as.character(exp(rep(seq(0, 0.4, length.out = 11), 2)))
  )
  fit <- joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = "sex",
    max_jp = 0,
    n_eq = 19L
  )

  plot <- plot_joinpoint(
    fit,
    group_id = 2L,
    show_aapc = TRUE,
    aapc_range = c(2002L, 2008L),
    ci_method = "Param"
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(unique(plot$data$group_id), 2L)
  expect_type(plot$data$.plot_time, "double")
  expect_type(plot$data$.plot_outcome, "double")
  label_data <- Filter(
    function(layer) ".label" %in% names(layer$data),
    plot$layers
  )
  expect_length(label_data, 1L)
  expect_match(label_data[[1L]]$data$.label, "AAPC =")
})

test_that("plot annotation formatting is stable", {
  expect_equal(
    crtrend:::format_joinpoint_p_value(c(NA, 0.0001, 0.1254), 3L),
    c("p = NA", "p < 0.001", "p = 0.125")
  )
  expect_equal(
    crtrend:::format_joinpoint_estimate(c(-1.234, 2), 2L),
    c("-1.23", "2.00")
  )
})

test_that("plotted NCI sample curves pass through predictions at every observed year", {
  skip_if_not_installed("ggplot2")
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    se = "se",
    max_jp = 4,
    n_eq = 9
  )
  original <- fit
  apc <- calc_apc(fit, ci_method = "Param")
  aapc <- calc_aapc(fit, 2000, 2021, ci_method = "Param")
  local_mocked_bindings(fit_candidate = function(...) stop("must not refit"))
  plot <- plot_joinpoint(
    fit,
    show_apc = TRUE,
    show_aapc = TRUE,
    ci_method = "Param"
  )
  curve <- plot$layers[[2L]]$data
  years <- fit$model$year
  annual <- curve$.plot_outcome[match(years, curve$.plot_time)]
  expect_equal(annual, as.numeric(predict(fit)), tolerance = 1e-12)
  expect_equal(round(annual[years == 2006], 4), 47.9692)
  expect_equal(round(annual[years == 2006], 2), 47.97)
  expect_gte(nrow(curve), 202L)
  expect_equal(
    curve$.plot_outcome[curve$.plot_time == 2012],
    rep(annual[13], 2)
  )
  expect_equal(curve$segment[curve$.plot_time == 2012], 1:2)
  built <- ggplot2::ggplot_build(plot)$data[[2L]]
  expect_equal(built$y[match(years, built$x)], annual, tolerance = 1e-12)
  expect_equal(length(unique(built$group)), 2L)
  expect_identical(fit, original)
  expect_equal(calc_apc(fit, ci_method = "Param"), apc)
  expect_equal(calc_aapc(fit, 2000, 2021, ci_method = "Param"), aapc)
})

test_that("dense curves handle irregular character times and zero joinpoints", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(time = as.character(c(2000, 2001, 2004, 2007, 2013, 2021)))
  data$rate <- exp(.1 * (as.numeric(data$time) - 2000))
  fit <- joinpoint(data, "time", "rate", max_jp = 0, n_eq = 2)
  curve <- plot_joinpoint(fit, show_points = FALSE)$layers[[1L]]$data
  newdata <- data.frame(time = curve$.plot_time)
  expect_equal(curve$.plot_outcome, as.numeric(predict(fit, newdata = newdata)))
  expect_equal(all(as.numeric(data$time) %in% curve$.plot_time), TRUE)
  expect_equal(anyDuplicated(curve$.plot_time), 0L)
  expect_equal(range(curve$.plot_time), c(2000, 2021))
  expect_equal(unique(curve$segment), 1L)
  expect_equal(is.unsorted(curve$.plot_time), FALSE)
  middle <- which.min(abs(curve$.plot_time - 2010.5))
  endpoint_chord <- mean(curve$.plot_outcome[c(1L, nrow(curve))])
  expect_gt(endpoint_chord - curve$.plot_outcome[middle], 1)
})

test_that("multi-joinpoint curves keep groups and shared endpoints separate", {
  skip_if_not_installed("ggplot2")
  fit <- joinpoint(
    jp_data,
    "year_index",
    "rate",
    group_vars = "sex",
    se = "se",
    min_jp = 2,
    max_jp = 2,
    n_eq = 2
  )
  selected <- fit$models[[2L]]
  plot <- plot_joinpoint(fit, group_id = 2L)
  curve <- plot$layers[[2L]]$data
  expect_equal(unique(curve$group_id), 2L)
  expect_equal(unique(curve$sex), selected$group$sex)
  expect_equal(unique(curve$segment), 1:3)
  expect_equal(
    curve$.plot_outcome,
    as.numeric(predict(
      selected,
      newdata = data.frame(year_index = curve$.plot_time)
    ))
  )
  for (tau in selected$joinpoints) {
    rows <- curve$.plot_time == tau
    expect_equal(sum(rows), 2L)
    expect_equal(diff(curve$.plot_outcome[rows]), 0)
  }
  for (segment in 1:3) {
    times <- curve$.plot_time[curve$segment == segment]
    expect_equal(
      range(times),
      unlist(selected$segments[segment, c("start", "end")]),
      ignore_attr = TRUE
    )
    expect_equal(anyDuplicated(times), 0L)
  }
})
