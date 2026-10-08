#' Plot fitted joinpoint trends
#'
#' @description
#' `plot_joinpoint()` draws the observed outcome, segment-specific fitted
#' curves, and selected joinpoints from a [joinpoint()] result. Fitted segments
#' use distinct colours. Optional labels show each segment's APC and p-value
#' and each group's overall or custom-range AAPC and p-value.
#'
#' @param x A [joinpoint()] result.
#' @param group_id Optional integer vector selecting `model_info$group_id`.
#'   `NULL` selects every group.
#' @param facet Logical; facet by all grouping variables. The default is
#'   `TRUE`.
#' @param scales Scale setting passed to [ggplot2::facet_wrap()].
#' @param show_points Show the original observations.
#' @param show_joinpoints Show vertical dashed lines at selected joinpoints.
#' @param show_apc Label every successfully fitted segment with APC and its
#'   two-sided p-value.
#' @param show_aapc Label every successfully fitted group with AAPC and its
#'   two-sided p-value.
#' @param aapc_range `NULL` uses each group's complete fitted range. Otherwise,
#'   supply two integer years defining one common AAPC range.
#' @param ci_method Confidence-interval and p-value method passed to
#'   [calc_apc()] and the AAPC calculation. The default is `"EQ"`.
#' @param level Confidence level used for APC and AAPC inference.
#' @param n_resamples Number of empirical-quantile resamples. `NULL` reuses the
#'   `n_eq` value stored in `x$control`.
#' @param seed Optional residual-resampling seed.
#' @param estimate_digits Number of decimal places used for APC and AAPC.
#' @param p_digits Number of decimal places used for p-values; smaller values
#'   are displayed with a less-than sign.
#' @param point_size Size of the observed-data points.
#' @param line_width Width of the fitted segment lines.
#' @param text_size Size of APC and AAPC annotation text.
#'
#' @return A `ggplot` object that can be modified with additional ggplot2
#'   layers, scales, and themes.
#'
#' @details
#' Each fitted segment uses 101 equally spaced prediction times together with
#' all observed times in that segment. Predictions are on the original outcome
#' scale, so log-linear segments are drawn as exponential curves rather than
#' straight lines joining their endpoints. Every observed time is included
#' exactly, and adjacent curves share the prediction at their joinpoint.
#' APC and AAPC labels use [calc_apc()] and the cached inference stored by
#' [joinpoint()]; plotting does not refit model selection. The default AAPC is
#' calculated separately over each group's complete fitted range.
#' Consecutive APC labels alternate below and above their fitted segments and
#' use an offset that keeps the text away from the line. Segment colours are
#' intentionally not shown in a legend because their order follows time.
#'
#' APC and AAPC confidence intervals describe trend estimates and are not
#' prediction intervals for the fitted curve. This function therefore does
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
#'
#' @importFrom rlang .data
#' @export
plot_joinpoint <- function(
  x,
  group_id = NULL,
  facet = TRUE,
  scales = c("free_y", "fixed", "free_x", "free"),
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
  text_size = 3.0
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
    show_points = show_points,
    show_joinpoints = show_joinpoints,
    show_apc = show_apc,
    show_aapc = show_aapc,
    aapc_range = aapc_range,
    ci_method = ci_method,
    level = level,
    n_resamples = n_resamples,
    seed = seed,
    estimate_digits = estimate_digits,
    p_digits = p_digits,
    point_size = point_size,
    line_width = line_width,
    text_size = text_size
  )
  models <- x$models[match(options$group_id, x$model_info$group_id)]
  observed <- joinpoint_plot_observations(models)
  fitted_segments <- joinpoint_plot_segments(models)
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
        colour = "#D55E00",
        linetype = "dashed",
        linewidth = 0.5
      )
  }
  if (options$show_apc) {
    apc_labels <- joinpoint_plot_apc_labels(x, models, options)
    if (nrow(apc_labels) > 0L) {
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
            vjust = .data[[".vjust"]]
          ),
          size = options$text_size,
          check_overlap = TRUE,
          show.legend = FALSE
        )
    }
  }
  if (options$show_aapc) {
    aapc_labels <- joinpoint_plot_aapc_labels(models, options)
    if (nrow(aapc_labels) > 0L) {
      plot <- plot +
        ggplot2::geom_text(
          data = aapc_labels,
          mapping = ggplot2::aes(
            x = .data[[".plot_time"]],
            y = .data[[".plot_outcome"]],
            label = .data[[".label"]]
          ),
          size = options$text_size,
          hjust = -0.05,
          vjust = 1.2,
          fontface = "bold",
          colour = "#009E73"
        )
    }
  }
  group_vars <- x$control$group_vars
  if (options$facet && length(group_vars) > 0L) {
    plot <- plot +
      ggplot2::facet_wrap(
        stats::reformulate(group_vars),
        scales = options$scales
      )
  }
  plot +
    ggplot2::labs(
      x = x$control$time,
      y = x$control$outcome
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.18))
    ) +
    ggplot2::scale_colour_discrete(guide = "none") +
    ggplot2::theme_minimal()
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
  show_points,
  show_joinpoints,
  show_apc,
  show_aapc,
  aapc_range,
  ci_method,
  level,
  n_resamples,
  seed,
  estimate_digits,
  p_digits,
  point_size,
  line_width,
  text_size
) {
  if (!inherits(x, "crtrend_joinpoint")) {
    stop("`x` must be a `crtrend_joinpoint` object.", call. = FALSE)
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
  if (!is.null(aapc_range)) {
    if (!is.numeric(aapc_range) || length(aapc_range) != 2L) {
      stop("`aapc_range` must be NULL or two integer years.", call. = FALSE)
    }
    validate_interval(aapc_range[1L], aapc_range[2L])
    aapc_range <- as.integer(aapc_range)
  }
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
  for (argument in c("point_size", "line_width", "text_size")) {
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
  list(
    group_id = group_id,
    facet = facet,
    scales = scales,
    show_points = show_points,
    show_joinpoints = show_joinpoints,
    show_apc = show_apc,
    show_aapc = show_aapc,
    aapc_range = aapc_range,
    ci_method = ci_method,
    level = level,
    n_resamples = as.integer(n_resamples),
    seed = seed,
    estimate_digits = as.integer(estimate_digits),
    p_digits = as.integer(p_digits),
    point_size = point_size,
    line_width = line_width,
    text_size = text_size
  )
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

#' Internal helper: Build dense response-scale prediction curves by segment.
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
    "APC = ",
    format_joinpoint_estimate(result$apc, options$estimate_digits),
    "%\n",
    format_joinpoint_p_value(result$p_value, options$p_digits)
  )
  result
}

#' Internal helper: Build whole-range or custom-range AAPC labels.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_plot_aapc_labels <- function(models, options) {
  rows <- lapply(models, function(model) {
    interval <- options$aapc_range
    if (is.null(interval)) {
      interval <- c(model$start, model$end)
    }
    calc_aapc_group(
      fit = model,
      start = as.integer(interval[1L]),
      end = as.integer(interval[2L]),
      level = options$level,
      ci_method = options$ci_method,
      n_resamples = options$n_resamples,
      seed = options$seed
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
  result$.label <- paste0(
    "AAPC = ",
    format_joinpoint_estimate(result$aapc, options$estimate_digits),
    "%\n",
    format_joinpoint_p_value(result$p_value, options$p_digits)
  )
  result
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
