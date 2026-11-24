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
  expect_match(labels, "%\\*", all = FALSE)
  expect_no_match(labels, "APC =|p [=<]")
  expect_match(labels, "AAPC:", all = FALSE)
  built <- ggplot2::ggplot_build(plot)
  expect_length(unique(built$data[[2L]]$colour), 2L)
  expect_equal(plot$layers[[1L]]$geom_params$pointsize, NULL)
  expect_equal(plot$layers[[1L]]$aes_params$size, 1.2)
  expect_equal(plot$layers[[2L]]$aes_params$linewidth, 0.6)
  expect_equal(plot$scales$get_scales("colour")$guide, "none")
  expect_equal(unique(built$data[[3L]]$colour), "#D0D0D0")
  apc_layer <- Filter(
    function(layer) {
      ".vjust" %in% names(layer$data) && "segment" %in% names(layer$data)
    },
    plot$layers
  )
  expect_length(apc_layer, 1L)
  expect_setequal(unique(apc_layer[[1L]]$data$.vjust), c(-0.6, 1.5))
  expect_equal(
    rlang::as_label(apc_layer[[1L]]$mapping$size),
    ".apc_text_size"
  )
  expect_equal(
    rlang::as_label(apc_layer[[1L]]$mapping$angle),
    ".line_angle"
  )
  expect_equal(apc_layer[[1L]]$aes_params$family, "sans")
  expect_equal(apc_layer[[1L]]$data$.apc_text_size, rep(3, 4L))
  expect_equal(
    unique(built$data[[4L]]$colour),
    unique(built$data[[2L]]$colour)
  )
  aapc_layer <- Filter(
    function(layer) {
      ".label" %in% names(layer$data) && any(grepl("^AAPC:", layer$data$.label))
    },
    plot$layers
  )
  expect_length(aapc_layer, 1L)
  expect_equal(
    rlang::as_label(aapc_layer[[1L]]$mapping$size),
    ".aapc_text_size"
  )
  expect_equal(aapc_layer[[1L]]$aes_params$family, "sans")
  expect_equal(aapc_layer[[1L]]$aes_params$colour, "black")
  expect_equal(plot$theme$text$family, "sans")
  expect_equal(plot$theme$plot.caption$colour, "#999999")
  expect_equal(plot$theme$plot.caption$size, 3 * 2.845276 - 1)
  expect_equal(plot$theme$legend.position, "top")
  expect_equal(plot$theme$aspect.ratio, 0.58)
  expect_equal(plot$theme$panel.grid.minor, ggplot2::element_blank())
})

test_that("plot_joinpoint applies one font family to every text element", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(
    年份 = 2000:2010,
    发病率 = exp(seq(0, 0.3, length.out = 11))
  )
  fit <- joinpoint(data, "年份", "发病率", max_jp = 0, n_eq = 2)
  plot <- plot_joinpoint(
    fit,
    show_apc = TRUE,
    show_aapc = TRUE,
    ci_method = "Param",
    font_family = "Chinese Test Font"
  )

  annotation_layers <- Filter(
    function(layer) ".label" %in% names(layer$data),
    plot$layers
  )
  expect_equal(plot$labels$x, "年份")
  expect_equal(plot$labels$y, "发病率")
  expect_equal(plot$theme$text$family, "Chinese Test Font")
  expect_equal(
    unname(vapply(
      annotation_layers,
      \(layer) layer$aes_params$family,
      character(1L)
    )),
    rep("Chinese Test Font", 2L)
  )
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
  expect_identical(
    label_data[[1L]]$data$.label,
    "AAPC 2002-2008: 4.08%***"
  )
})

test_that("AAPC labels support multiple named ranges and lastN shortcuts", {
  skip_if_not_installed("ggplot2")
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    max_jp = 0,
    n_eq = 2
  )
  plot <- plot_joinpoint(
    fit,
    show_aapc = TRUE,
    aapc_range = list(
      overall = NULL,
      last10 = "last10",
      early = c(2000L, 2010L)
    ),
    ci_method = "Param"
  )
  labels <- plot$layers[[3L]]$data
  expect_identical(
    labels$.aapc_caption,
    c("AAPC:", "AAPC (last10):", "AAPC (early):")
  )
  expect_equal(
    labels[, c("start", "end")],
    data.frame(
      start = c(2000L, 2012L, 2000L),
      end = c(2021L, 2021L, 2010L)
    )
  )
  expect_equal(labels$.vjust, c(1.2, 3.4, 5.6))
  unnamed <- plot_joinpoint(
    fit,
    show_aapc = TRUE,
    aapc_range = list(c(2000L, 2010L), c(2014L, 2021L)),
    ci_method = "Param"
  )$layers[[3L]]$data
  expect_identical(
    unnamed$.aapc_caption,
    c("AAPC 2000-2010:", "AAPC 2014-2021:")
  )
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
  expect_equal(
    crtrend:::format_joinpoint_significance(
      c(NA, 0.05, 0.049, 0.01, 0.009, 0.001, 0.0009)
    ),
    c("", "", "*", "*", "**", "**", "***")
  )
})

test_that("APC labels use stars, follow their lines, and have a separate size", {
  skip_if_not_installed("ggplot2")
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    max_jp = 0,
    n_eq = 2
  )
  plot <- plot_joinpoint(
    fit,
    show_apc = TRUE,
    apc_text_size = 4.5,
    ci_method = "Param"
  )
  label_layer <- plot$layers[[3L]]
  expect_match(label_layer$data$.label, "^-2\\.26%\\*\\*\\*$")
  expect_equal(label_layer$data$.apc_text_size, 4.5)
  expect_equal(
    sign(label_layer$data$.line_angle),
    sign(label_layer$data$apc)
  )
  expect_gt(
    abs(label_layer$data$.line_angle),
    abs(atan(log1p(label_layer$data$apc / 100) / log(10)) * 180 / pi)
  )
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[3L]]$size, 4.5)
  expect_equal(built$data[[3L]]$colour, built$data[[2L]]$colour[1L])
  expect_equal(built$data[[3L]]$angle, label_layer$data$.line_angle)
})

test_that("AAPC labels use a configurable colour", {
  skip_if_not_installed("ggplot2")
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    max_jp = 0,
    n_eq = 2
  )
  plot <- plot_joinpoint(
    fit,
    show_aapc = TRUE,
    aapc_text_color = "#6A3D9A",
    ci_method = "Param"
  )
  expect_equal(plot$layers[[3L]]$aes_params$colour, "#6A3D9A")
  expect_equal(
    unique(ggplot2::ggplot_build(plot)$data[[3L]]$colour),
    "#6A3D9A"
  )
})

test_that("AAPC labels compact sex groups and append significance stars", {
  skip_if_not_installed("ggplot2")
  fit <- joinpoint(
    jp_data,
    "year",
    "rate",
    group_vars = "sex",
    max_jp = 0,
    n_eq = 2
  )
  plot <- plot_joinpoint(
    fit,
    overlay_by = "sex",
    show_aapc = TRUE,
    ci_method = "Param"
  )
  labels <- plot$layers[[3L]]$data$.label
  expect_identical(
    labels,
    c(
      "AAPC:",
      "合计 -2.26%***",
      "男性 -2.47%***",
      "女性 -2.15%***"
    )
  )
  column <- plot$layers[[3L]]$data
  expect_equal(column$.vjust, c(1.2, 3.4, 5.6, 7.8))
  row <- plot_joinpoint(
    fit,
    overlay_by = "sex",
    show_aapc = TRUE,
    aapc_layout = "row",
    ci_method = "Param"
  )$layers[[3L]]$data
  expect_equal(row$.vjust, rep(1.2, 4L))
  expect_equal(is.unsorted(row$.plot_time), FALSE)
  expect_equal(row$.hjust, rep(0, 4L))
  sized <- plot_joinpoint(
    fit,
    overlay_by = "sex",
    show_aapc = TRUE,
    aapc_text_size = 4.5,
    ci_method = "Param"
  )$layers[[3L]]$data
  expect_equal(sized$.aapc_text_size, rep(4.5, 4L))
  multiple <- plot_joinpoint(
    fit,
    overlay_by = "sex",
    show_aapc = TRUE,
    aapc_range = list(overall = NULL, last10 = "last10"),
    ci_method = "Param"
  )$layers[[3L]]$data
  expect_identical(
    multiple$.label[multiple$.aapc_caption],
    c("AAPC:", "AAPC (last10):")
  )
  expect_equal(multiple$.vjust, c(1.2, 3.4, 5.6, 7.8, 10, 12.2, 14.4, 16.6))
  multiple_row <- plot_joinpoint(
    fit,
    overlay_by = "sex",
    show_aapc = TRUE,
    aapc_range = list(overall = NULL, last10 = "last10"),
    aapc_layout = "row",
    ci_method = "Param"
  )$layers[[3L]]$data
  expect_equal(multiple_row$.vjust, rep(c(1.2, 3.4), each = 4L))
})

test_that("plotted NCI lines are straight on the APC log scale", {
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
  expect_equal(built$y[match(years, built$x)], log10(annual), tolerance = 1e-12)
  for (segment in unique(built$group)) {
    selected <- built$group == segment
    slopes <- diff(built$y[selected]) / diff(built$x[selected])
    expect_equal(slopes, rep(slopes[1L], length(slopes)), tolerance = 1e-10)
  }
  expect_equal(length(unique(built$group)), 2L)
  label_data <- plot$layers[[4L]]$data
  for (segment in unique(label_data$segment)) {
    curve_data <- curve[curve$segment == segment, , drop = FALSE]
    line_slope <- (log10(curve_data$.plot_outcome[nrow(curve_data)]) -
      log10(curve_data$.plot_outcome[1L])) /
      (curve_data$.plot_time[nrow(curve_data)] - curve_data$.plot_time[1L])
    label_angle <- label_data$.line_angle[label_data$segment == segment]
    expect_equal(sign(label_angle), sign(line_slope))
    expect_gt(abs(label_angle), abs(atan(line_slope) * 180 / pi))
  }
  expect_identical(fit, original)
  expect_equal(calc_apc(fit, ci_method = "Param"), apc)
  expect_equal(calc_aapc(fit, 2000, 2021, ci_method = "Param"), aapc)
})

test_that("fitted lines handle irregular character times on the APC log scale", {
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
  built <- ggplot2::ggplot_build(plot_joinpoint(
    fit,
    show_points = FALSE
  ))$data[[1L]]
  slopes <- diff(built$y) / diff(built$x)
  expect_equal(slopes, rep(slopes[1L], length(slopes)), tolerance = 1e-12)
  expect_equal(built$y, log10(curve$.plot_outcome), tolerance = 1e-12)
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

overlay_plot_fit <- function() {
  data <- expand.grid(year = 2000:2011, sex = 1:2, cancer = c("Lung", "Colon"))
  data$risk <- with(
    data,
    exp(
      -2 +
        0.02 * (year - 2000) +
        0.06 * pmax(year - 2005, 0) +
        0.4 * (sex == 2) +
        0.1 * (cancer == "Colon") +
        0.001 * sin(year)
    )
  )
  joinpoint(
    data,
    "year",
    "risk",
    group_vars = c("sex", "cancer"),
    min_jp = 1,
    max_jp = 1,
    n_eq = 2
  )
}

test_that("AAPC overlay labels use the selected grouping variable", {
  skip_if_not_installed("ggplot2")
  plot <- plot_joinpoint(
    overlay_plot_fit(),
    overlay_by = "cancer",
    show_aapc = TRUE,
    ci_method = "Param"
  )
  aapc_layer <- Filter(
    function(layer) ".label" %in% names(layer$data),
    plot$layers
  )[[1L]]
  labels <- aapc_layer$data
  expect_equal(labels$.label[labels$.aapc_caption], rep("AAPC:", 2L))
  expect_match(
    paste(labels$.label[!labels$.aapc_caption], collapse = " "),
    "Lung"
  )
  expect_match(
    paste(labels$.label[!labels$.aapc_caption], collapse = " "),
    "Colon"
  )
})

test_that("overlay_by shares panels and uses segment shades within each group", {
  skip_if_not_installed("ggplot2")
  fit <- overlay_plot_fit()
  original <- fit
  local_mocked_bindings(fit_candidate = function(...) stop("must not refit"))
  plot <- plot_joinpoint(
    fit,
    overlay_by = "sex",
    show_apc = TRUE,
    show_aapc = TRUE,
    ci_method = "Param"
  )
  built <- ggplot2::ggplot_build(plot)
  expect_equal(nrow(built$layout$layout), 2L)
  expect_named(plot$facet$params$facets, "cancer")
  expect_null(plot$scales$get_scales("colour")$name)
  expect_equal(plot$scales$get_scales("colour")$guide, "legend")
  expect_equal(built$plot$scales$get_scales("colour")$get_labels(), c("1", "2"))
  colour_scale <- built$plot$scales$get_scales("colour")
  colours <- colour_scale$map(colour_scale$get_breaks())
  for (i in seq_len(4L)) {
    expect_setequal(unique(plot$layers[[i]]$data$sex), 1:2)
  }
  expect_setequal(unique(built$data[[1L]]$colour), colours)
  expect_setequal(unique(built$data[[5L]]$colour), c("#999999", colours))
  expect_equal(
    unique(built$data[[5L]]$colour[built$data[[5L]]$label == "AAPC:"]),
    "#999999"
  )
  expect_equal(
    unique(built$data[[5L]]$colour[grepl("^男性", built$data[[5L]]$label)]),
    colours[1L]
  )
  expect_equal(
    unique(built$data[[5L]]$colour[grepl("^女性", built$data[[5L]]$label)]),
    colours[2L]
  )
  expect_setequal(unique(built$data[[3L]]$colour), colours)
  expect_setequal(
    unique(plot$layers[[3L]]$data$.plot_colour),
    colour_scale$get_breaks()
  )
  expect_equal(plot$layers[[3L]]$show.legend, FALSE)
  curves <- plot$layers[[2L]]$data
  labels <- plot$layers[[4L]]$data
  expect_setequal(unique(curves$.plot_colour), unique(labels$.plot_colour))
  for (sex in 1:2) {
    keys <- unique(curves$.plot_colour[curves$sex == sex])
    rgb <- grDevices::col2rgb(colour_scale$map(keys))
    expect_equal(length(keys), 2L)
    expect_gt(sum(rgb[, 1L]), sum(rgb[, 2L]))
  }
  expect_length(unique(built$data[[2L]]$group), 8L)
  for (panel in split(built$data[[2L]], built$data[[2L]]$PANEL)) {
    expect_length(unique(panel$colour), 4L)
    expect_length(unique(panel$group), 4L)
  }
  aapc <- plot$layers[[5L]]$data
  expect_equal(aapc$.label[aapc$.aapc_caption], rep("AAPC:", 2L))
  expect_match(aapc$.label[!aapc$.aapc_caption], "^(男性|女性) ")
  expect_equal(nrow(aapc), 6L)
  custom_aapc_colour <- ggplot2::ggplot_build(plot_joinpoint(
    fit,
    overlay_by = "sex",
    show_aapc = TRUE,
    aapc_text_color = "#6A3D9A",
    ci_method = "Param"
  ))$data[[4L]]$colour
  expect_setequal(unique(custom_aapc_colour), c("#999999", colours))
  expect_identical(fit, original)
  lines_only <- plot_joinpoint(fit, overlay_by = "sex", show_points = FALSE)
  lines_built <- ggplot2::ggplot_build(lines_only)
  expect_equal(
    lines_built$plot$scales$get_scales("colour")$get_labels(),
    c("1", "2")
  )
  expect_equal(
    length(lines_built$plot$scales$get_scales("colour")$get_breaks()),
    2L
  )
  expect_s3_class(ggplot2::ggplotGrob(lines_only), "gtable")
})

test_that("overlay facets respect other variables, filtering and facet FALSE", {
  skip_if_not_installed("ggplot2")
  fit <- overlay_plot_fit()
  default <- plot_joinpoint(fit)
  expect_named(default$facet$params$facets, c("sex", "cancer"))
  one_column <- plot_joinpoint(fit, ncol = 1L)
  expect_equal(max(ggplot2::ggplot_build(one_column)$layout$layout$COL), 1L)
  one_row <- plot_joinpoint(fit, nrow = 1L)
  expect_equal(max(ggplot2::ggplot_build(one_row)$layout$layout$ROW), 1L)
  swapped <- plot_joinpoint(fit, overlay_by = "cancer")
  expect_named(swapped$facet$params$facets, "sex")
  selected <- plot_joinpoint(fit, overlay_by = "sex", group_id = c(1L, 2L))
  expect_setequal(unique(selected$data$group_id), c(1L, 2L))
  selected_curve <- selected$layers[[2L]]$data
  expect_equal(
    selected_curve$.plot_outcome,
    default$layers[[2L]]$data$.plot_outcome[
      default$layers[[2L]]$data$group_id %in% c(1L, 2L)
    ]
  )
  combined <- plot_joinpoint(
    fit,
    overlay_by = "sex",
    facet = FALSE,
    show_aapc = TRUE,
    ci_method = "Param"
  )
  expect_equal(nrow(ggplot2::ggplot_build(combined)$layout$layout), 1L)
  expect_length(combined$layers[[4L]]$data$.label, 5L)

  single_variable <- joinpoint(
    jp_data,
    "year",
    "rate",
    group_vars = "sex",
    max_jp = 0,
    n_eq = 2
  )
  one_panel <- plot_joinpoint(single_variable, overlay_by = "sex")
  expect_equal(nrow(ggplot2::ggplot_build(one_panel)$layout$layout), 1L)
})

test_that("overlay_by requires one existing grouping variable", {
  skip_if_not_installed("ggplot2")
  fit <- overlay_plot_fit()
  expect_snapshot(error = TRUE, plot_joinpoint(fit, overlay_by = "year"))
  expect_snapshot(
    error = TRUE,
    plot_joinpoint(fit, overlay_by = c("sex", "cancer"))
  )
  expect_snapshot(error = TRUE, plot_joinpoint(fit, overlay_by = NA_character_))
  expect_snapshot(error = TRUE, plot_joinpoint(fit, overlay_by = 1))
  expect_snapshot(error = TRUE, plot_joinpoint(fit, overlay_by = character()))
  expect_snapshot(error = TRUE, plot_joinpoint(fit, font_family = ""))
  expect_snapshot(
    error = TRUE,
    plot_joinpoint(fit, font_family = NA_character_)
  )
  expect_snapshot(error = TRUE, plot_joinpoint(fit, apc_text_size = 0))
  expect_snapshot(error = TRUE, plot_joinpoint(fit, aapc_range = "last1"))
  expect_snapshot(error = TRUE, plot_joinpoint(fit, ncol = 2L, nrow = 2L))
  expect_snapshot(error = TRUE, plot_joinpoint(fit, ncol = 0L))
  expect_snapshot(
    error = TRUE,
    plot_joinpoint(fit, aapc_text_color = "not-a-colour")
  )
})
