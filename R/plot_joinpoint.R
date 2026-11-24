#' Plot fitted joinpoint trends
#'
#' @description
#' `plot_joinpoint()` draws the observed outcome, segment-specific fitted
#' lines, and selected joinpoints from a [joinpoint()] result. By default,
#' fitted segments use distinct colours. Optional labels show APC and AAPC
#' values with significance stars.
#'
#' @param x A [joinpoint()] result.
#' @param group_id Optional integer vector selecting `model_info$group_id`.
#'   `NULL` selects every group.
#' @param facet Logical; facet by grouping variables other than `overlay_by`.
#'   The default is `TRUE`.
#' @param scales Scale setting passed to [ggplot2::facet_wrap()].
#' @param ncol Optional positive integer number of facet columns. Specify at
#'   most one of `ncol` and `nrow`; ignored when no facets are drawn.
#' @param nrow Optional positive integer number of facet rows. Specify at most
#'   one of `ncol` and `nrow`; ignored when no facets are drawn.
#' @param show_points Show the original observations.
#' @param show_joinpoints Show vertical dashed lines at selected joinpoints.
#'   Lines are light grey by default and use the corresponding group colour
#'   when `overlay_by` is set.
#' @param show_apc Label every successfully fitted segment with its APC value
#'   and significance stars. Labels are rotated to run parallel to their fitted
#'   lines and use the matching line colour.
#' @param show_aapc Label every successfully fitted group and requested range
#'   with AAPC. AAPC
#'   values append one, two, or three stars for two-sided p-values below 0.05,
#'   0.01, or 0.001, respectively. With `overlay_by`, values for its levels
#'   are shown in their legend colours below a grey `AAPC:` caption per panel.
#' @param aapc_layout Layout for AAPC labels with `overlay_by`: `"column"`
#'   stacks the coloured values below the caption, while `"row"` places them
#'   after the caption on one line.
#' @param aapc_text_size Size of AAPC labels, including the grey `AAPC:`
#'   caption. Defaults to `text_size`.
#' @param aapc_text_color One valid R colour used for AAPC labels when
#'   `overlay_by = NULL`. Defaults to `"black"`; overlay labels retain their
#'   group colours.
#' @param aapc_range AAPC ranges to display. `NULL` uses each group's complete
#'   fitted range, and a two-integer vector defines one common range. Supply
#'   one or more `"lastN"` shortcuts (for example, `"last10"`) to use each
#'   group's final `N` calendar years, or a list of `NULL`, two-year vectors,
#'   and `"lastN"` shortcuts to show multiple ranges. List names are shown in
#'   the labels; unnamed intervals use `AAPC start-end:`. A `"lastN"` range
#'   uses all available years when a group has fewer than `N` years.
#' @param ci_method Confidence-interval and p-value method passed to
#'   [calc_apc()] and the AAPC calculation. The default is `"EQ"`.
#' @param level Confidence level used for APC and AAPC inference.
#' @param n_resamples Number of empirical-quantile resamples. `NULL` reuses the
#'   `n_eq` value stored in `x$control`.
#' @param seed Optional residual-resampling seed.
#' @param estimate_digits Number of decimal places used for APC and AAPC.
#' @param p_digits Retained for backward compatibility and currently unused.
#' @param point_size Size of the observed-data points.
#' @param line_width Width of the fitted segment lines.
#' @param text_size Default size of APC and AAPC annotation text.
#' @param apc_text_size Size of APC labels. Defaults to `text_size`.
#' @param overlay_by Optional single grouping-variable name, for example
#'   `"sex"`. Its levels share a panel and use distinct colour families with
#'   a legend; segments use progressively darker shades within each family.
#'   The remaining grouping variables define facets. `NULL` keeps the default
#'   segment colours and facets by every grouping variable. With `facet = FALSE`,
#'   all selected groups share one panel.
#' @param font_family Font family shared by axis titles, tick labels, facet
#'   labels, legends, and APC/AAPC annotations. The default is the portable
#'   system sans-serif family. For Chinese output, supply an installed family
#'   such as `"Heiti SC"`, `"Microsoft YaHei"`, or `"Noto Sans CJK SC"`.
#'
#' @return A `ggplot` object that can be modified with additional ggplot2
#'   layers, scales, and themes.
#'
#' @details
#' Joinpoint regression is linear in `log(outcome)`. The plot therefore uses a
#' logarithmic y-axis so each fitted segment is a straight line whose slope
#' corresponds to its APC, while axis labels remain in the original outcome
#' units. Each segment contains 101 equally spaced prediction times together
#' with all observed times. Every observed time is included exactly, and
#' adjacent lines share the prediction at their joinpoint.
#' APC and AAPC labels use [calc_apc()] and the cached inference stored by
#' [joinpoint()]; plotting does not refit model selection. The default AAPC is
#' calculated separately over each group's complete fitted range.
#' `aapc_range` can instead request one or more named ranges. For example,
#' `list(overall = NULL, last10 = "last10", early = c(2000L, 2010L))`
#' displays the complete-fit, last-ten-calendar-year, and early-period AAPCs.
#' Consecutive APC labels alternate below and above their fitted segments and
#' use an offset that keeps the text away from the line. Their angle follows
#' the segment slope on the plotted log10 outcome scale. Segment colours are
#' intentionally not shown in a legend because their order follows time.
#' Joinpoint reference lines are light grey by default. With `overlay_by`,
#' colour families identify that variable's levels; joinpoint lines, points,
#' segment lines, and APC labels use matching shades from light to dark in time
#' order. AAPC values for levels of `overlay_by` use their legend colours below
#' a grey `AAPC:` caption in each panel. If no grouping variables remain for
#' faceting, the plot has a single panel. Set `aapc_layout = "row"` to place
#' the caption and values on one line; the default `"column"` uses separated
#' rows for legibility.
#'
#' APC and AAPC confidence intervals describe trend estimates and are not
#' prediction intervals for the fitted line. This function therefore does
#' not draw a confidence ribbon.
#'
#' @examples
#' years <- 2000:2015
#' example_data <- data.frame(
#'   year = rep(years, 2),
#'   sex = rep(c("Male", "Female"), each = length(years))
#' )
#' example_data$risk <- with(
#'   example_data,
#'   exp(-3 + 0.02 * (year - 2000) + 0.2 * (sex == "Female"))
#' )
#' fit <- joinpoint(
#'   example_data,
#'   time = "year",
#'   outcome = "risk",
#'   group_vars = "sex",
#'   max_jp = 0,
#'   n_eq = 19
#' )
#' plot_joinpoint(fit, show_apc = TRUE, show_aapc = TRUE)
#' plot_joinpoint(
#'   fit,
#'   show_aapc = TRUE,
#'   aapc_range = list(overall = NULL, last10 = "last10")
#' )
#' plot_joinpoint(fit, overlay_by = "sex", show_apc = TRUE)
#'
#' @importFrom rlang .data
#' @export
plot_joinpoint <- function(
  x,
  group_id = NULL,
  facet = TRUE,
  scales = c("free_y", "fixed", "free_x", "free"),
  ncol = NULL,
  nrow = NULL,
  show_points = TRUE,
  show_joinpoints = TRUE,
  show_apc = FALSE,
  show_aapc = FALSE,
  aapc_range = NULL,
  ci_method = c("EQ", "Param"),
  level = 0.95,
  n_resamples = NULL,
  seed = NULL,
  estimate_digits = 2L,
  p_digits = 3L,
  point_size = 1.2,
  line_width = 0.6,
  text_size = 3.0,
  apc_text_size = text_size,
  overlay_by = NULL,
  font_family = "sans",
  aapc_layout = c("column", "row"),
  aapc_text_size = text_size,
  aapc_text_color = "black"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package `ggplot2` is required to use `plot_joinpoint()`.",
      call. = FALSE
    )
  }
  options <- validate_plot_joinpoint_inputs(
    x = x,
    group_id = group_id,
    facet = facet,
    scales = scales,
    ncol = ncol,
    nrow = nrow,
    show_points = show_points,
    show_joinpoints = show_joinpoints,
    show_apc = show_apc,
    show_aapc = show_aapc,
    aapc_layout = aapc_layout,
    aapc_range = aapc_range,
    ci_method = ci_method,
    level = level,
    n_resamples = n_resamples,
    seed = seed,
    estimate_digits = estimate_digits,
    p_digits = p_digits,
    point_size = point_size,
    line_width = line_width,
    text_size = text_size,
    apc_text_size = apc_text_size,
    overlay_by = overlay_by,
    font_family = font_family,
    aapc_text_size = aapc_text_size,
    aapc_text_color = aapc_text_color
  )
  models <- x$models[match(options$group_id, x$model_info$group_id)]
  group_vars <- setdiff(x$control$group_vars, options$overlay_by)
  facet_vars <- if (options$facet) group_vars else character()
  observed <- joinpoint_plot_observations(models)
  fitted_segments <- joinpoint_plot_segments(models)
  range_columns <- c(
    "group_id",
    x$control$group_vars,
    ".plot_time",
    ".plot_outcome"
  )
  plot_range_data <- rbind(
    observed[, range_columns, drop = FALSE],
    fitted_segments[, range_columns, drop = FALSE]
  )
  segment_levels <- sort(unique(fitted_segments$segment))
  fitted_segments$.segment <- factor(
    fitted_segments$segment,
    levels = segment_levels,
    labels = paste("Segment", segment_levels)
  )
  joinpoints <- x$joinpoints[
    x$joinpoints$group_id %in% options$group_id,
    ,
    drop = FALSE
  ]

  plot <- ggplot2::ggplot(data = observed)
  if (options$show_points) {
    plot <- plot +
      ggplot2::geom_point(
        data = observed,
        mapping = ggplot2::aes(
          x = .data[[".plot_time"]],
          y = .data[[".plot_outcome"]]
        ),
        size = options$point_size,
        colour = "#333333",
        alpha = 0.75
      )
  }
  plot <- plot +
    ggplot2::geom_line(
      data = fitted_segments,
      mapping = ggplot2::aes(
        x = .data[[".plot_time"]],
        y = .data[[".plot_outcome"]],
        group = interaction(.data[["group_id"]], .data[["segment"]]),
        colour = .data[[".segment"]]
      ),
      linewidth = options$line_width
    )
  if (options$show_joinpoints && nrow(joinpoints) > 0L) {
    plot <- plot +
      ggplot2::geom_vline(
        data = joinpoints,
        mapping = ggplot2::aes(xintercept = .data[["joinpoint"]]),
        colour = "#D0D0D0",
        linetype = "dashed",
        linewidth = 0.4
      )
  }
  if (options$show_apc) {
    apc_labels <- joinpoint_plot_apc_labels(x, models, options)
    if (nrow(apc_labels) > 0L) {
      apc_labels <- joinpoint_plot_apc_label_angles(
        apc_labels,
        plot_range_data,
        facet_vars,
        options$scales
      )
      apc_labels$.segment <- factor(
        apc_labels$segment,
        levels = segment_levels,
        labels = paste("Segment", segment_levels)
      )
      plot <- plot +
        ggplot2::geom_text(
          data = apc_labels,
          mapping = ggplot2::aes(
            x = .data[[".plot_time"]],
            y = .data[[".plot_outcome"]],
            label = .data[[".label"]],
            colour = .data[[".segment"]],
            vjust = .data[[".vjust"]],
            size = .data[[".apc_text_size"]],
            angle = .data[[".line_angle"]]
          ),
          family = options$font_family,
          lineheight = 0.95,
          check_overlap = TRUE,
          show.legend = FALSE
        )
    }
  }
  if (options$show_aapc) {
    aapc_labels <- joinpoint_plot_aapc_labels(models, options, facet_vars)
    if (!is.null(options$overlay_by)) {
      aapc_labels <- joinpoint_plot_overlay_aapc_labels(
        aapc_labels,
        options$overlay_by,
        facet_vars,
        options$aapc_layout
      )
    }
    if (nrow(aapc_labels) > 0L) {
      plot <- plot +
        ggplot2::geom_text(
          data = aapc_labels,
          mapping = ggplot2::aes(
            x = .data[[".plot_time"]],
            y = .data[[".plot_outcome"]],
            label = .data[[".label"]],
            vjust = .data[[".vjust"]],
            hjust = .data[[".hjust"]],
            size = .data[[".aapc_text_size"]]
          ),
          family = options$font_family,
          lineheight = 0.95,
          fontface = "bold",
          colour = options$aapc_text_color,
          show.legend = FALSE
        )
    }
  }
  if (!is.null(options$overlay_by)) {
    plot <- joinpoint_plot_overlay(
      plot,
      options$overlay_by,
      facet_vars
    )
  }
  if (length(facet_vars) > 0L) {
    plot <- plot +
      ggplot2::facet_wrap(
        stats::reformulate(facet_vars),
        scales = options$scales,
        ncol = options$ncol,
        nrow = options$nrow
      )
  }
  plot +
    ggplot2::labs(
      x = x$control$time,
      y = x$control$outcome
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.04, 0.05))
    ) +
    ggplot2::scale_y_log10(
      expand = ggplot2::expansion(mult = c(0.12, 0.22))
    ) +
    (if (is.null(options$overlay_by)) {
      ggplot2::scale_colour_discrete(guide = "none")
    }) +
    ggplot2::scale_size_identity(guide = "none") +
    joinpoint_plot_theme(options$font_family, options$text_size)
}

#' Internal helper: Validate joinpoint plotting arguments.
#'
#' @rdname crtrend-internals
#' @keywords internal
validate_plot_joinpoint_inputs <- function(
  x,
  group_id,
  facet,
  scales,
  ncol,
  nrow,
  show_points,
  show_joinpoints,
  show_apc,
  show_aapc,
  aapc_layout,
  aapc_range,
  ci_method,
  level,
  n_resamples,
  seed,
  estimate_digits,
  p_digits,
  point_size,
  line_width,
  text_size,
  apc_text_size = text_size,
  overlay_by = NULL,
  font_family = "sans",
  aapc_text_size = text_size,
  aapc_text_color = "black"
) {
  if (!inherits(x, "crtrend_joinpoint")) {
    stop("`x` must be a `crtrend_joinpoint` object.", call. = FALSE)
  }
  if (
    !is.null(overlay_by) &&
      (!is.character(overlay_by) ||
        length(overlay_by) != 1L ||
        is.na(overlay_by) ||
        !overlay_by %in% x$control$group_vars)
  ) {
    stop(
      "`overlay_by` must be NULL or the name of one grouping variable in `x`.",
      call. = FALSE
    )
  }
  available_group_id <- x$model_info$group_id
  if (is.null(group_id)) {
    group_id <- available_group_id
  }
  if (
    !is.numeric(group_id) ||
      length(group_id) == 0L ||
      anyNA(group_id) ||
      is.complex(group_id) ||
      any(!is.finite(group_id)) ||
      any(abs(group_id) > .Machine$integer.max) ||
      any(group_id != as.integer(group_id)) ||
      any(!group_id %in% available_group_id)
  ) {
    stop("`group_id` must select available integer group IDs.", call. = FALSE)
  }
  group_id <- unique(as.integer(group_id))
  scales <- match.arg(scales, c("free_y", "fixed", "free_x", "free"))
  if (!is.null(ncol) && !is.null(nrow)) {
    stop("Specify at most one of `ncol` and `nrow`.", call. = FALSE)
  }
  for (argument in c("ncol", "nrow")) {
    value <- get(argument)
    if (
      !is.null(value) &&
        (!is.numeric(value) ||
          length(value) != 1L ||
          is.na(value) ||
          !is_bounded_integer(value) ||
          value < 1L ||
          value != as.integer(value))
    ) {
      stop(
        "`",
        argument,
        "` must be NULL or one positive integer.",
        call. = FALSE
      )
    }
  }
  aapc_layout <- match.arg(aapc_layout, c("column", "row"))
  ci_method <- match.arg(
    ci_method[1L],
    c("EQ", "Param", "empirical_quantile", "parametric")
  )
  ci_method <- normalize_ci_method(ci_method)
  if (is.null(n_resamples)) {
    n_resamples <- x$control$n_eq
  }
  validate_calc_inputs(x, level, ci_method, n_resamples, seed)
  for (argument in c(
    "facet",
    "show_points",
    "show_joinpoints",
    "show_apc",
    "show_aapc"
  )) {
    value <- get(argument)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  aapc_ranges <- normalize_plot_aapc_ranges(aapc_range)
  for (argument in c("estimate_digits", "p_digits")) {
    value <- get(argument)
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !is_bounded_integer(value) ||
        value < 0L ||
        value != as.integer(value)
    ) {
      stop("`", argument, "` must be one non-negative integer.", call. = FALSE)
    }
  }
  for (argument in c(
    "point_size",
    "line_width",
    "text_size",
    "apc_text_size",
    "aapc_text_size"
  )) {
    value <- get(argument)
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !is.finite(value) ||
        value <= 0
    ) {
      stop("`", argument, "` must be one positive number.", call. = FALSE)
    }
  }
  if (
    !is.character(font_family) ||
      length(font_family) != 1L ||
      is.na(font_family) ||
      !nzchar(font_family)
  ) {
    stop("`font_family` must be one non-empty character string.", call. = FALSE)
  }
  valid_aapc_text_color <-
    is.character(aapc_text_color) &&
    length(aapc_text_color) == 1L &&
    !is.na(aapc_text_color) &&
    nzchar(aapc_text_color) &&
    !inherits(
      try(grDevices::col2rgb(aapc_text_color), silent = TRUE),
      "try-error"
    )
  if (!valid_aapc_text_color) {
    stop("`aapc_text_color` must be one valid R colour.", call. = FALSE)
  }
  list(
    group_id = group_id,
    overlay_by = overlay_by,
    facet = facet,
    scales = scales,
    ncol = if (is.null(ncol)) NULL else as.integer(ncol),
    nrow = if (is.null(nrow)) NULL else as.integer(nrow),
    show_points = show_points,
    show_joinpoints = show_joinpoints,
    show_apc = show_apc,
    show_aapc = show_aapc,
    aapc_layout = aapc_layout,
    aapc_ranges = aapc_ranges,
    ci_method = ci_method,
    level = level,
    n_resamples = as.integer(n_resamples),
    seed = seed,
    estimate_digits = as.integer(estimate_digits),
    p_digits = as.integer(p_digits),
    point_size = point_size,
    line_width = line_width,
    text_size = text_size,
    apc_text_size = apc_text_size,
    font_family = font_family,
    aapc_text_size = aapc_text_size,
    aapc_text_color = aapc_text_color
  )
}

joinpoint_plot_theme <- function(font_family, text_size) {
  ggplot2::theme_minimal(base_size = 11, base_family = font_family) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(12, 16, 12, 12),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 10)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 10)
      ),
      axis.text = ggplot2::element_text(colour = "#4D4D4D"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "#E5E5E5",
        linewidth = 0.35
      ),
      panel.spacing = grid::unit(14, "pt"),
      aspect.ratio = 0.58,
      strip.text = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(b = 7)
      ),
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.margin = ggplot2::margin(0, 0, 6, 0),
      legend.spacing.x = grid::unit(6, "pt"),
      plot.caption = ggplot2::element_text(
        colour = "#999999",
        size = max(1, text_size * 2.845276 - 1)
      )
    )
}

joinpoint_plot_overlay <- function(
  plot,
  overlay_by,
  facet_vars
) {
  observed <- plot$data
  observed$.plot_group <- as.factor(observed[[overlay_by]])
  group_levels <- levels(droplevels(observed$.plot_group))
  observed$.plot_group <- factor(observed$.plot_group, levels = group_levels)
  group_keys <- paste0("group_", seq_along(group_levels))
  base_colours <- grDevices::hcl(
    h = seq(15, 375, length.out = length(group_levels) + 1L)[seq_along(
      group_levels
    )],
    c = 65,
    l = 55
  )
  names(base_colours) <- group_keys
  segments <- sort(unique(unlist(lapply(plot$layers, function(layer) {
    layer$data$segment
  }))))
  palette <- c(aapc_caption = "#999999", base_colours)
  shades <- if (length(segments) <= 1L) {
    0
  } else {
    seq(0.3, -0.25, length.out = length(segments))
  }
  for (i in seq_along(segments)) {
    colours <- joinpoint_plot_shade(base_colours, shades[i])
    names(colours) <- paste0(group_keys, "_segment_", segments[i])
    palette <- c(palette, colours)
  }
  plot$data <- observed
  for (layer in plot$layers) {
    rows <- match(layer$data$group_id, observed$group_id)
    layer$data$.plot_group <- observed$.plot_group[rows]
    layer$data$.plot_colour <- group_keys[as.integer(layer$data$.plot_group)]
    if (".aapc_caption" %in% names(layer$data)) {
      layer$data$.plot_colour[layer$data$.aapc_caption] <- "aapc_caption"
    }
    if ("segment" %in% names(layer$data) && nrow(layer$data) > 0L) {
      layer$data$.plot_colour <- paste0(
        layer$data$.plot_colour,
        "_segment_",
        layer$data$segment
      )
    }
    layer$mapping$colour <- ggplot2::aes(
      colour = .data[[".plot_colour"]]
    )$colour
    layer$aes_params$colour <- NULL
    layer$show.legend <- inherits(layer$geom, "GeomLine")
  }
  plot +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = palette,
      limits = names(palette),
      breaks = group_keys,
      labels = group_levels,
      guide = "legend"
    )
}

joinpoint_plot_shade <- function(colours, amount) {
  rgb <- grDevices::col2rgb(colours) / 255
  rgb <- if (amount >= 0) rgb + (1 - rgb) * amount else rgb * (1 + amount)
  grDevices::rgb(rgb[1L, ], rgb[2L, ], rgb[3L, ])
}

#' Internal helper: Build observed data for a joinpoint plot.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_observations <- function(models) {
  rows <- lapply(models, function(model) {
    result <- model$model
    result$group_id <- model$group_id
    result$.plot_time <- coerce_numeric_column(
      result[[model$control$time]],
      model$control$time
    )
    result$.plot_outcome <- coerce_numeric_column(
      result[[model$control$outcome]],
      model$control$outcome
    )
    result
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Internal helper: Build fitted prediction lines by segment.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_segments <- function(models) {
  rows <- lapply(models, function(model) {
    if (model$status != "ok" || nrow(model$segments) == 0L) {
      return(NULL)
    }
    n <- nrow(model$segments)
    observed_time <- coerce_numeric_column(
      model$model[[model$control$time]],
      model$control$time
    )
    grids <- lapply(seq_len(n), function(index) {
      start <- model$segments$start[index]
      end <- model$segments$end[index]
      sort(unique(c(
        seq(start, end, length.out = 101L),
        observed_time[observed_time >= start & observed_time <= end]
      )))
    })
    time <- unlist(grids, use.names = FALSE)
    newdata <- model$model[rep(1L, length(time)), , drop = FALSE]
    newdata[[model$control$time]] <- time
    group <- model$group[rep(1L, length(time)), , drop = FALSE]
    rownames(group) <- NULL
    data.frame(
      group_id = rep(model$group_id, length(time)),
      group,
      segment = rep(model$segments$segment, times = lengths(grids)),
      .plot_time = time,
      .plot_outcome = as.numeric(stats::predict(model, newdata = newdata)),
      check.names = FALSE
    )
  })
  joinpoint_plot_bind_rows(rows)
}

#' Internal helper: Build APC labels at fitted segment midpoints.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_apc_labels <- function(x, models, options) {
  result <- calc_apc(
    x,
    level = options$level,
    ci_method = options$ci_method,
    n_resamples = options$n_resamples,
    seed = options$seed
  )
  result <- result[
    result$group_id %in% options$group_id & result$status == "ok",
    ,
    drop = FALSE
  ]
  if (nrow(result) == 0L) {
    return(result)
  }
  result$.plot_time <- (result$start + result$end) / 2
  result$.plot_outcome <- NA_real_
  result$.vjust <- ifelse(result$segment %% 2L == 1L, 1.5, -0.6)
  result$.apc_text_size <- options$apc_text_size
  for (group_id in unique(result$group_id)) {
    selected <- result$group_id == group_id
    model <- models[[match(
      group_id,
      vapply(models, `[[`, integer(1L), "group_id")
    )]]
    newdata <- model$model[rep(1L, sum(selected)), , drop = FALSE]
    newdata[[model$control$time]] <- result$.plot_time[selected]
    result$.plot_outcome[selected] <- as.numeric(
      stats::predict(model, newdata = newdata)
    )
  }
  result$.label <- paste0(
    format_joinpoint_estimate(result$apc, options$estimate_digits),
    "%",
    format_joinpoint_significance(result$p_value)
  )
  result
}

#' Internal helper: Calculate display angles for APC labels.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_apc_label_angles <- function(
  labels,
  plot_data,
  facet_vars,
  scales,
  panel_aspect = 0.58
) {
  plot_data <- plot_data[
    is.finite(plot_data$.plot_time) &
      is.finite(plot_data$.plot_outcome) &
      plot_data$.plot_outcome > 0,
    ,
    drop = FALSE
  ]
  panels <- if (length(facet_vars) == 0L) {
    list(plot_data)
  } else {
    split_group_rows(plot_data, facet_vars)
  }
  free_x <- scales %in% c("free_x", "free")
  free_y <- scales %in% c("free_y", "free")
  labels$.line_angle <- 0
  for (panel in panels) {
    selected <- labels$group_id %in% panel$group_id
    x_data <- if (free_x) panel else plot_data
    y_data <- if (free_y) panel else plot_data
    x_span <- diff(range(x_data$.plot_time)) * 1.09
    y_span <- diff(range(log10(y_data$.plot_outcome))) * 1.34
    if (x_span > 0 && y_span > 0) {
      labels$.line_angle[selected] <- atan(
        log1p(labels$apc[selected] / 100) /
          log(10) *
          panel_aspect *
          x_span /
          y_span
      ) *
        180 /
        pi
    }
  }
  labels
}

#' Internal helper: Normalize requested AAPC plot ranges.
#'
#' @rdname crtrend-internals
#' @keywords internal
normalize_plot_aapc_ranges <- function(aapc_range) {
  invalid <- function() {
    stop(
      paste0(
        "`aapc_range` must be NULL, two integer years, one or more ",
        "'lastN' shortcuts, or a list of these."
      ),
      call. = FALSE
    )
  }
  make_entry <- function(value, name = NULL) {
    if (is.null(value)) {
      return(list(name = name, interval = NULL, last_n = NULL, full = TRUE))
    }
    if (is.numeric(value) && length(value) == 2L) {
      validate_interval(value[1L], value[2L])
      return(list(
        name = name,
        interval = as.integer(value),
        last_n = NULL,
        full = FALSE
      ))
    }
    if (is.character(value) && length(value) == 1L && !is.na(value)) {
      match <- regexec("^last([1-9][0-9]*)$", value)
      parts <- regmatches(value, match)[[1L]]
      if (length(parts) == 2L) {
        last_n <- suppressWarnings(as.integer(parts[2L]))
        if (!is.na(last_n) && last_n >= 2L) {
          if (is.null(name)) {
            name <- value
          }
          return(list(
            name = name,
            interval = NULL,
            last_n = last_n,
            full = FALSE
          ))
        }
      }
    }
    invalid()
  }
  if (is.null(aapc_range)) {
    return(list(make_entry(NULL)))
  }
  if (is.numeric(aapc_range) && length(aapc_range) == 2L) {
    return(list(make_entry(aapc_range)))
  }
  if (is.character(aapc_range) && length(aapc_range) > 0L) {
    return(lapply(aapc_range, make_entry))
  }
  if (!is.list(aapc_range) || length(aapc_range) == 0L) {
    invalid()
  }
  range_names <- names(aapc_range)
  if (is.null(range_names)) {
    range_names <- rep(NA_character_, length(aapc_range))
  }
  lapply(seq_along(aapc_range), function(index) {
    name <- range_names[index]
    if (is.na(name) || !nzchar(name)) {
      name <- NULL
    }
    make_entry(aapc_range[[index]], name)
  })
}

#' Internal helper: Resolve one AAPC plot range for a fitted group.
#'
#' @rdname crtrend-internals
#' @keywords internal
resolve_plot_aapc_range <- function(range, model) {
  if (range$full) {
    return(c(start = model$start, end = model$end))
  }
  if (!is.null(range$interval)) {
    return(stats::setNames(range$interval, c("start", "end")))
  }
  c(
    start = max(model$start, model$end - range$last_n + 1L),
    end = model$end
  )
}

#' Internal helper: Format an AAPC annotation caption.
#'
#' @rdname crtrend-internals
#' @keywords internal
format_plot_aapc_caption <- function(range, interval) {
  if (range$full) {
    return("AAPC:")
  }
  if (!is.null(range$name)) {
    return(paste0("AAPC (", range$name, "):"))
  }
  paste0("AAPC ", interval[1L], "-", interval[2L], ":")
}

#' Internal helper: Stack AAPC annotations in each panel.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_stack_aapc_labels <- function(labels, facet_vars) {
  panels <- split_group_rows(labels, facet_vars)
  result <- lapply(panels, function(panel) {
    panel <- panel[order(panel$.aapc_range_id, panel$group_id), , drop = FALSE]
    panel$.vjust <- 1.2 + 2.2 * (seq_len(nrow(panel)) - 1L)
    panel$.hjust <- -0.02
    panel
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result
}

#' Internal helper: Build whole-range or custom-range AAPC labels.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_aapc_labels <- function(models, options, facet_vars) {
  rows <- lapply(models, function(model) {
    do.call(
      rbind,
      lapply(seq_along(options$aapc_ranges), function(index) {
        range <- options$aapc_ranges[[index]]
        interval <- resolve_plot_aapc_range(range, model)
        result <- calc_aapc_group(
          fit = model,
          start = as.integer(interval[1L]),
          end = as.integer(interval[2L]),
          level = options$level,
          ci_method = options$ci_method,
          n_resamples = options$n_resamples,
          seed = options$seed
        )
        result$.aapc_range_id <- index
        result$.aapc_caption <- format_plot_aapc_caption(range, interval)
        result
      })
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- result[result$status == "ok", , drop = FALSE]
  if (nrow(result) == 0L) {
    return(result)
  }
  result$.plot_time <- -Inf
  result$.plot_outcome <- Inf
  result$.aapc_text_size <- options$aapc_text_size
  result$.aapc_value <- paste0(
    format_joinpoint_estimate(result$aapc, options$estimate_digits),
    "%",
    format_joinpoint_significance(result$p_value)
  )
  result$.label <- paste(result$.aapc_caption, result$.aapc_value)
  joinpoint_plot_stack_aapc_labels(result, facet_vars)
}

#' Internal helper: Build coloured AAPC labels for an overlay variable.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_overlay_aapc_labels <- function(
  labels,
  overlay_by,
  facet_vars,
  aapc_layout
) {
  if (nrow(labels) == 0L) {
    return(labels)
  }
  panels <- split_group_rows(labels, facet_vars)
  result <- lapply(panels, function(panel) {
    ranges <- split(panel, panel$.aapc_range_id, drop = TRUE)
    next_vjust <- 1.2
    blocks <- vector("list", length(ranges))
    for (index in seq_along(ranges)) {
      range_panel <- ranges[[index]]
      range_panel <- range_panel[order(range_panel$group_id), , drop = FALSE]
      values <- vapply(
        seq_len(nrow(range_panel)),
        function(index) {
          paste(
            format_joinpoint_group_label(
              range_panel[index, , drop = FALSE],
              overlay_by
            ),
            range_panel$.aapc_value[index]
          )
        },
        character(1L)
      )
      caption <- range_panel[1L, , drop = FALSE]
      caption$group_id <- NA_integer_
      caption$.label <- caption$.aapc_caption
      caption$.aapc_caption <- TRUE
      range_panel$.label <- values
      range_panel$.aapc_caption <- FALSE
      result <- rbind(caption, range_panel)
      if (aapc_layout == "row") {
        widths <- nchar(result$.label, type = "width") + 1L
        positions <- c(0L, cumsum(utils::head(widths, -1L))) / sum(widths)
        result$.plot_time <- min(range_panel$start) +
          positions * diff(range(range_panel$start, range_panel$end))
        result$.vjust <- next_vjust
        result$.hjust <- 0
        next_vjust <- next_vjust + 2.2
      } else {
        result$.vjust <- next_vjust + 2.2 * (seq_len(nrow(result)) - 1L)
        result$.hjust <- -0.02
        next_vjust <- max(result$.vjust) + 2.2
      }
      blocks[[index]] <- result
    }
    result <- do.call(rbind, blocks)
    rownames(result) <- NULL
    result
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result
}

#' Internal helper: Format an overlay level for an AAPC annotation.
#'
#' @rdname crtrend-internals
#' @keywords internal
format_joinpoint_group_label <- function(data, overlay_by) {
  value <- as.character(data[[overlay_by]][1L])
  if (tolower(overlay_by) != "sex") {
    return(value)
  }
  labels <- c(
    "0" = "\u5408\u8ba1",
    "1" = "\u7537\u6027",
    "2" = "\u5973\u6027",
    "all" = "\u5408\u8ba1",
    "total" = "\u5408\u8ba1",
    "male" = "\u7537\u6027",
    "female" = "\u5973\u6027",
    "\u7537" = "\u7537\u6027",
    "\u5973" = "\u5973\u6027"
  )
  translated <- unname(labels[tolower(value)])
  if (is.na(translated)) value else translated
}

#' Internal helper: Bind optional plot-data rows.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_bind_rows <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      group_id = integer(),
      segment = integer(),
      .plot_time = numeric(),
      .plot_outcome = numeric()
    ))
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Internal helper: Format APC and AAPC point estimates.
#'
#' @rdname crtrend-internals
#' @keywords internal
format_joinpoint_estimate <- function(x, digits) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format = "f", digits = digits)
  )
}

#' Internal helper: Format AAPC significance stars.
#'
#' @rdname crtrend-internals
#' @keywords internal
format_joinpoint_significance <- function(x) {
  vapply(
    x,
    function(value) {
      if (is.na(value) || value >= 0.05) {
        return("")
      }
      if (value < 0.001) {
        return("***")
      }
      if (value < 0.01) {
        return("**")
      }
      "*"
    },
    character(1L)
  )
}

#' Internal helper: Format p-values for plot annotations.
#'
#' @rdname crtrend-internals
#' @keywords internal
format_joinpoint_p_value <- function(x, digits) {
  threshold <- 10^-digits
  vapply(
    x,
    function(value) {
      if (is.na(value)) {
        return("p = NA")
      }
      if (value < threshold) {
        return(paste0(
          "p < ",
          formatC(threshold, format = "f", digits = digits)
        ))
      }
      paste0("p = ", formatC(value, format = "f", digits = digits))
    },
    character(1L)
  )
}
