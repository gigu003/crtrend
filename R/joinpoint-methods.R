#' Internal helper: Construct standard model components from grouped fits.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_model_components <- function(fits, n, group_vars, time) {
  response_fitted <- response_residuals <- rep(NA_real_, n)
  link_fitted <- link_residuals <- rep(NA_real_, n)
  coefficient_rows <- vector("list", length(fits))
  rank <- 0L
  df_residual <- 0L

  for (index in seq_along(fits)) {
    fit <- fits[[index]]
    if (fit$status != "ok") {
      next
    }
    rows <- fit$data$trend_row_id
    fitted_link <- as.numeric(stats::fitted(fit$model))
    link_fitted[rows] <- fitted_link
    link_residuals[rows] <- as.numeric(stats::residuals(fit$model))
    response_fitted[rows] <- exp(fitted_link)
    response_residuals[rows] <- fit$data$trend_outcome - exp(fitted_link)
    rank <- rank + fit$model$rank
    df_residual <- df_residual + stats::df.residual(fit$model)

    estimates <- stats::coef(fit$model)
    raw_terms <- names(estimates)
    hinge_number <- suppressWarnings(
      as.integer(sub("^joinpoint_", "", raw_terms))
    )
    is_hinge <- grepl("^joinpoint_[0-9]+$", raw_terms)
    joinpoint <- rep(NA_real_, length(estimates))
    joinpoint[is_hinge] <- fit$joinpoints[hinge_number[is_hinge]]
    term <- raw_terms
    term[raw_terms == "trend_time"] <- time
    term[is_hinge] <- paste0("slope_change_at_", joinpoint[is_hinge])
    component <- rep("slope_change", length(estimates))
    component[raw_terms == "(Intercept)"] <- "intercept"
    component[raw_terms == "trend_time"] <- "baseline_slope"
    group_key <- fit$group_key[rep(1L, length(estimates)), , drop = FALSE]
    rownames(group_key) <- NULL
    coefficient_rows[[index]] <- data.frame(
      group_id = rep(index, length(estimates)),
      group_key,
      term = term,
      estimate = unname(estimates),
      component = component,
      joinpoint = joinpoint,
      check.names = FALSE
    )
  }

  coefficient_rows <- Filter(Negate(is.null), coefficient_rows)
  if (length(coefficient_rows) > 0L) {
    coefficients <- do.call(rbind, coefficient_rows)
    rownames(coefficients) <- NULL
  } else {
    coefficients <- fits[[1L]]$group_key[FALSE, , drop = FALSE]
    coefficients <- data.frame(
      group_id = integer(),
      coefficients,
      check.names = FALSE
    )
    coefficients$term <- character()
    coefficients$estimate <- numeric()
    coefficients$component <- character()
    coefficients$joinpoint <- numeric()
  }

  list(
    coefficients = coefficients,
    response_fitted = response_fitted,
    response_residuals = response_residuals,
    link_fitted = link_fitted,
    link_residuals = link_residuals,
    rank = rank,
    df_residual = df_residual
  )
}

#' Internal helper: Return analysis data from a computational or public fit.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_analysis_data <- function(fit) {
  fit$.analysis %||% fit$data
}

#' Internal helper: Return the underlying least-squares model.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_engine <- function(fit) {
  fit$.engine %||% fit$model
}

#' Internal helper: Return the group key from a computational or public fit.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_group <- function(fit) {
  group <- fit$group %||% fit$group_key
  if (!is.null(fit$group_id)) {
    group <- data.frame(
      group_id = fit$group_id,
      group,
      check.names = FALSE
    )
  }
  group
}

#' Internal helper: Return the group weighting setting.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_is_weighted <- function(fit) {
  isTRUE(fit$control$weighted %||% fit$weighted)
}

#' Internal helper: Return the minimum segment observation setting.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_min_observations <- function(fit) {
  fit$control$min_obs %||% fit$min_observations
}

#' Internal helper: Convert computational fits to standard group model objects.
#'
#' @rdname crtrend-internals
#' @keywords internal
build_joinpoint_models <- function(
  fits,
  model_frame,
  call,
  formula,
  time,
  outcome,
  group_vars
) {
  lapply(
    seq_along(fits),
    function(index) {
      build_joinpoint_model(
        fit = fits[[index]],
        group_id = index,
        model_frame = model_frame,
        call = call,
        formula = formula,
        time = time,
        outcome = outcome,
        group_vars = group_vars
      )
    }
  )
}

#' Internal helper: Construct one standard group model object.
#'
#' @rdname crtrend-internals
#' @keywords internal
build_joinpoint_model <- function(
  fit,
  group_id,
  model_frame,
  call,
  formula,
  time,
  outcome,
  group_vars
) {
  rows <- fit$data$trend_row_id
  group_model_frame <- model_frame[rows, , drop = FALSE]
  rownames(group_model_frame) <- row.names(model_frame)[rows]
  n <- length(rows)
  fitted_values <- residuals <- linear_predictors <- link_residuals <-
    rep(NA_real_, n)
  coefficients <- numeric()
  rank <- 0L
  df_residual <- 0L

  if (fit$status == "ok") {
    linear_predictors <- as.numeric(stats::fitted(fit$model))
    link_residuals <- as.numeric(stats::residuals(fit$model))
    fitted_values <- exp(linear_predictors)
    residuals <- fit$data$trend_outcome - fitted_values
    coefficients <- rename_joinpoint_coefficients(
      stats::coef(fit$model),
      fit$joinpoints,
      time
    )
    rank <- fit$model$rank
    df_residual <- stats::df.residual(fit$model)
  }
  observation_names <- row.names(group_model_frame)
  names(fitted_values) <- observation_names
  names(residuals) <- observation_names
  names(linear_predictors) <- observation_names
  names(link_residuals) <- observation_names

  segments <- group_segment_table(fit)
  selection <- list(
    bic = fit$bic,
    bic3 = fit$bic3,
    modified_bic = fit$modified_bic,
    weighted_bic = fit$weighted_bic,
    partial_r_squared = fit$partial_r_squared,
    dds_delta = fit$dds_delta,
    dds_method = fit$dds_method,
    dds_sigma_squared = fit$dds_sigma_squared,
    dds_variance_df = fit$dds_variance_df,
    dds_bic_delta = fit$dds_bic_delta,
    dds_bic3_delta = fit$dds_bic3_delta,
    dds_delta_max = fit$dds_delta_max,
    dds_delta_min = fit$dds_delta_min,
    dds_decision_rule = fit$dds_decision_rule,
    score = fit$selection_score,
    p_value = fit$selection_p_value,
    p_value_lower = fit$selection_p_value_lower,
    p_value_upper = fit$selection_p_value_upper,
    confidence_level = fit$selection_p_value_confidence_level,
    n_permutations = fit$selection_n_permutations,
    stopped_early = fit$selection_stopped_early,
    effective_max_jp = fit$effective_max_joinpoints,
    tests = fit$selection_tests,
    candidates = fit$candidates %||% empty_model_candidates()
  )
  structure(
    list(
      call = call,
      formula = formula,
      terms = stats::terms(formula),
      coefficients = coefficients,
      residuals = residuals,
      fitted.values = fitted_values,
      linear.predictors = linear_predictors,
      link.residuals = link_residuals,
      rank = rank,
      df.residual = df_residual,
      model = group_model_frame,
      na.action = NULL,
      weights = fit$model$weights %||% NULL,
      covariance = fit$model$covariance %||% NULL,
      effects = fit$model$effects %||% NULL,
      assign = fit$model$assign %||% NULL,
      qr = fit$model$qr %||% NULL,
      contrasts = fit$model$contrasts %||% NULL,
      xlevels = fit$model$xlevels %||% NULL,
      group_id = group_id,
      group = fit$group_key,
      status = fit$status,
      n_joinpoints = fit$n_joinpoints,
      joinpoints = fit$joinpoints,
      segments = segments,
      selection = selection,
      inference = fit$inference,
      resamples = list(
        method = "EQ",
        n = fit$eq_n_resamples,
        draws = fit$eq_draws
      ),
      control = list(
        time = time,
        outcome = outcome,
        group_vars = group_vars,
        min_obs = fit$min_observations %||% NA_integer_,
        weighted = fit$weighted %||% FALSE,
        wbic_weight = attr(fit$data, "crtrend_wbic_weight", exact = TRUE) %||%
          "qr",
        error = fit$error_structure %||% "Ind",
        rho = fit$rho %||% 0,
        cov_scale = fit$covariance_scale %||% "outcome"
      ),
      start = fit$start,
      end = fit$end,
      .engine = fit$model,
      .analysis = fit$data
    ),
    class = "crtrend_joinpoint_model"
  )
}

#' Internal helper: Rename computational coefficient terms for public output.
#'
#' @rdname crtrend-internals
#' @keywords internal
rename_joinpoint_coefficients <- function(coefficients, joinpoints, time) {
  terms <- names(coefficients)
  terms[terms == "trend_time"] <- time
  for (index in seq_along(joinpoints)) {
    terms[terms == paste0("joinpoint_", index)] <- paste0(
      "slope_change_at_",
      joinpoints[index]
    )
  }
  names(coefficients) <- terms
  coefficients
}

#' Internal helper: Construct segment metadata for one fitted group.
#'
#' @rdname crtrend-internals
#' @keywords internal
group_segment_table <- function(fit) {
  if (fit$status != "ok") {
    return(data.frame(
      segment = integer(),
      start = numeric(),
      end = numeric(),
      slope = numeric()
    ))
  }
  slopes <- fit$inference$slope_estimates %||% nci_grid_slope_estimates(fit)
  data.frame(
    segment = seq_len(fit$n_joinpoints + 1L),
    start = c(fit$start, fit$joinpoints),
    end = c(fit$joinpoints, fit$end),
    slope = slopes
  )
}

#' Internal helper: Combine numeric joinpoint locations across groups.
#'
#' @rdname crtrend-internals
#' @keywords internal
combine_joinpoint_table <- function(models, group_vars) {
  rows <- lapply(models, function(model) {
    n <- length(model$joinpoints)
    if (n == 0L) {
      return(NULL)
    }
    group <- model$group[rep(1L, n), , drop = FALSE]
    rownames(group) <- NULL
    data.frame(
      group_id = rep(model$group_id, n),
      group,
      index = seq_len(n),
      joinpoint = model$joinpoints,
      check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) > 0L) {
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    return(result)
  }
  result <- models[[1L]]$group[FALSE, group_vars, drop = FALSE]
  result <- data.frame(group_id = integer(), result, check.names = FALSE)
  result$index <- integer()
  result$joinpoint <- numeric()
  result
}

#' Internal helper: Combine segment metadata across groups.
#'
#' @rdname crtrend-internals
#' @keywords internal
combine_segment_table <- function(models, group_vars) {
  rows <- lapply(models, function(model) {
    n <- nrow(model$segments)
    if (n == 0L) {
      return(NULL)
    }
    group <- model$group[rep(1L, n), , drop = FALSE]
    rownames(group) <- NULL
    data.frame(
      group_id = rep(model$group_id, n),
      group,
      model$segments,
      check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) > 0L) {
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    return(result)
  }
  result <- models[[1L]]$group[FALSE, group_vars, drop = FALSE]
  result <- data.frame(group_id = integer(), result, check.names = FALSE)
  result$segment <- integer()
  result$start <- numeric()
  result$end <- numeric()
  result$slope <- numeric()
  result
}

#' Internal helper: Combine permutation-selection audits across groups.
#'
#' @rdname crtrend-internals
#' @keywords internal
combine_selection_test_table <- function(models, group_vars) {
  rows <- lapply(models, function(model) {
    tests <- model$selection$tests
    if (is.null(tests) || nrow(tests) == 0L) {
      return(NULL)
    }
    group <- model$group[rep(1L, nrow(tests)), , drop = FALSE]
    rownames(group) <- NULL
    data.frame(
      group_id = rep(model$group_id, nrow(tests)),
      group,
      tests,
      check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) > 0L) {
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    return(result)
  }
  group <- models[[1L]]$group[FALSE, group_vars, drop = FALSE]
  data.frame(
    group_id = integer(),
    group,
    empty_selection_tests(),
    check.names = FALSE
  )
}

#' @rdname joinpoint
#'
#' @param object,x,formula A fitted `crtrend_joinpoint` or
#'   `crtrend_joinpoint_model` object.
#' @param type Scale of the returned values: `"response"` for the original
#'   outcome scale or `"link"` for the log scale.
#' @param newdata Optional data frame containing the time variable and all
#'   grouping variables. Groups not present among the successful fits return
#'   `NA` predictions.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return `coef()` returns a tidy coefficient data frame. `fitted()`,
#'   `residuals()`, and `predict()` return numeric vectors. `formula()` returns
#'   the model formula, `nobs()` returns the input row count, `model.frame()`
#'   returns the original model variables, and `summary()` returns a concise
#'   grouped-model summary.
#' @importFrom stats nobs
#' @export
coef.crtrend_joinpoint <- function(object, ...) {
  object$coefficients
}

#' @rdname joinpoint
#' @export
fitted.crtrend_joinpoint <- function(
  object,
  type = c("response", "link"),
  ...
) {
  type <- match.arg(type)
  if (type == "response") object$fitted.values else object$linear.predictors
}

#' @rdname joinpoint
#' @export
residuals.crtrend_joinpoint <- function(
  object,
  type = c("response", "link"),
  ...
) {
  type <- match.arg(type)
  if (type == "response") object$residuals else object$link.residuals
}

#' @rdname joinpoint
#' @export
formula.crtrend_joinpoint <- function(x, ...) {
  x$formula
}

#' @rdname joinpoint
#' @export
nobs.crtrend_joinpoint <- function(object, ...) {
  nrow(object$model)
}

#' @rdname joinpoint
#' @export
model.frame.crtrend_joinpoint <- function(formula, ...) {
  formula$model
}

#' @rdname joinpoint
#' @export
predict.crtrend_joinpoint <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  ...
) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    return(stats::fitted(object, type = type))
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  required <- c(object$control$group_vars, object$control$time)
  missing_columns <- setdiff(required, names(newdata))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing columns in `newdata`: ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  time <- coerce_numeric_column(
    newdata[[object$control$time]],
    object$control$time
  )
  prediction <- rep(NA_real_, nrow(newdata))
  names(prediction) <- row.names(newdata)

  for (model in object$models) {
    if (model$status != "ok") {
      next
    }
    selected <- rep(TRUE, nrow(newdata))
    for (group_var in object$control$group_vars) {
      selected <- selected &
        !is.na(newdata[[group_var]]) &
        newdata[[group_var]] == model$group[[group_var]][1L]
    }
    if (!any(selected)) {
      next
    }
    prediction[selected] <- stats::predict(
      model,
      newdata = newdata[selected, , drop = FALSE],
      type = type
    )
  }
  prediction
}

#' @rdname joinpoint
#' @export
summary.crtrend_joinpoint <- function(object, ...) {
  structure(
    list(
      call = object$call,
      formula = object$formula,
      method = object$method,
      nobs = nobs(object),
      n_groups = nrow(object$model_info),
      df.residual = object$df.residual,
      model_info = object$model_info,
      coefficients = object$coefficients
    ),
    class = "summary.crtrend_joinpoint"
  )
}

#' @rdname joinpoint
#' @export
print.crtrend_joinpoint <- function(x, ...) {
  cat(x$method, "\n")
  cat("Formula: ")
  print(x$formula)
  cat(
    "Groups:",
    nrow(x$model_info),
    "  Observations:",
    nobs(x),
    "  Selection:",
    x$control$selection_method,
    "\n\n"
  )
  print(joinpoint_model_info_display(x), row.names = FALSE)
  invisible(x)
}

#' @rdname joinpoint
#' @export
coef.crtrend_joinpoint_model <- function(object, ...) {
  object$coefficients
}

#' @rdname joinpoint
#' @export
fitted.crtrend_joinpoint_model <- function(
  object,
  type = c("response", "link"),
  ...
) {
  type <- match.arg(type)
  if (type == "response") object$fitted.values else object$linear.predictors
}

#' @rdname joinpoint
#' @export
residuals.crtrend_joinpoint_model <- function(
  object,
  type = c("response", "link"),
  ...
) {
  type <- match.arg(type)
  if (type == "response") object$residuals else object$link.residuals
}

#' @rdname joinpoint
#' @export
formula.crtrend_joinpoint_model <- function(x, ...) {
  x$formula
}

#' @rdname joinpoint
#' @export
nobs.crtrend_joinpoint_model <- function(object, ...) {
  nrow(object$model)
}

#' @rdname joinpoint
#' @export
model.frame.crtrend_joinpoint_model <- function(formula, ...) {
  formula$model
}

#' @rdname joinpoint
#' @export
predict.crtrend_joinpoint_model <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  ...
) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    return(stats::fitted(object, type = type))
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  if (object$status != "ok") {
    prediction <- rep(NA_real_, nrow(newdata))
    names(prediction) <- row.names(newdata)
    return(prediction)
  }
  time_name <- object$control$time
  if (!time_name %in% names(newdata)) {
    stop("Missing time column in `newdata`: ", time_name, ".", call. = FALSE)
  }
  time <- coerce_numeric_column(newdata[[time_name]], time_name)
  model_data <- data.frame(trend_time = time)
  for (index in seq_along(object$joinpoints)) {
    model_data[[paste0("joinpoint_", index)]] <- pmax(
      time - object$joinpoints[index],
      0
    )
  }
  design <- joinpoint_design(object$joinpoints, model_data)
  link_prediction <- as.numeric(
    design %*% stats::coef(object$.engine)
  )
  prediction <- if (type == "link") link_prediction else exp(link_prediction)
  names(prediction) <- row.names(newdata)
  prediction
}

#' @rdname joinpoint
#' @export
summary.crtrend_joinpoint_model <- function(object, ...) {
  structure(
    list(
      call = object$call,
      formula = object$formula,
      group_id = object$group_id,
      group = object$group,
      status = object$status,
      nobs = nobs(object),
      df.residual = object$df.residual,
      coefficients = object$coefficients,
      joinpoints = object$joinpoints,
      segments = object$segments,
      selection = object$selection
    ),
    class = "summary.crtrend_joinpoint_model"
  )
}

#' @rdname joinpoint
#' @export
print.crtrend_joinpoint_model <- function(x, ...) {
  cat("Joinpoint model for group", x$group_id, "\n")
  if (ncol(x$group) > 0L) {
    print(x$group, row.names = FALSE)
  }
  cat("Status:", x$status, "  Observations:", nobs(x), "\n")
  if (x$status == "ok") {
    cat("Joinpoints:", format_joinpoints(x$joinpoints), "\n")
    print(x$coefficients)
  }
  invisible(x)
}

#' @rdname joinpoint
#' @export
print.summary.crtrend_joinpoint_model <- function(x, ...) {
  cat("Joinpoint model summary for group", x$group_id, "\n")
  if (ncol(x$group) > 0L) {
    print(x$group, row.names = FALSE)
  }
  cat(
    "Status:",
    x$status,
    "  Observations:",
    x$nobs,
    "  Residual df:",
    x$df.residual,
    "\n"
  )
  if (x$status == "ok") {
    cat("Joinpoints:", format_joinpoints(x$joinpoints), "\n\n")
    print(x$coefficients)
  }
  invisible(x)
}

#' @rdname joinpoint
#' @export
print.summary.crtrend_joinpoint <- function(x, ...) {
  cat(x$method, "\n")
  cat("Call:\n")
  print(x$call)
  cat("\nFormula: ")
  print(x$formula)
  cat(
    "Groups:",
    x$n_groups,
    "  Observations:",
    x$nobs,
    "  Residual df:",
    x$df.residual,
    "\n\nModel selection:\n"
  )
  print(joinpoint_model_info_display(x), row.names = FALSE)
  cat("\nCoefficients:\n")
  print(x$coefficients, row.names = FALSE)
  invisible(x)
}

#' Internal helper: Add formatted joinpoints for model printing only.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_model_info_display <- function(object) {
  result <- object$model_info
  result$joinpoints <- vapply(
    object$models,
    function(model) format_joinpoints(model$joinpoints),
    character(1L)
  )
  position <- match("n_joinpoints", names(result))
  result[c(
    names(result)[seq_len(position)],
    "joinpoints",
    names(result)[seq.int(position + 1L, ncol(result) - 1L)]
  )]
}
