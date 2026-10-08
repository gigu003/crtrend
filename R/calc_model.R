#' Extract model data, coefficients, and candidate-selection summaries
#'
#' @description
#' `calc_model()` assembles the fitted data, model statistics, standard and
#' general coefficient parameterizations, parameter covariance and correlation
#' matrices, and candidate-model scores from a [joinpoint()] result.
#'
#' @param x A [joinpoint()] result, including grouped models.
#' @param level Confidence level used to flag significant yearly APC values.
#'   The default is `0.95`.
#' @param ci_method Method for the yearly APC p-values, `"EQ"` (default) or
#'   `"Param"`. The aliases `"empirical_quantile"` and `"parametric"` are
#'   accepted. Coefficient SEs, t statistics, and p-values are always parametric.
#' @param n_resamples Number of EQ resamples. `NULL` reuses `x$control$n_eq`.
#' @param seed Optional seed passed to [calc_apc()] for EQ resampling.
#'
#' @return A named list with these components:
#' \describe{
#'   \item{data}{One row per input observation, in original input order.
#'     Columns include `row_id`, `group_id`, grouping variables, `status`,
#'     numeric `time`, `observed`, original-scale `se` (missing if not supplied),
#'     `fitted`, `residual`, `linear_predictor`, `link_residual`, `joinpoint`
#'     (the index of a joinpoint at this observation), `segment`, `apc`,
#'     `apc_p_value`, `apc_significant`, and `apc_method`. Segment and APC fields
#'     are missing at joinpoint observations, matching the NCI Data tab.}
#'   \item{statistics}{One row per group with fit and inference status,
#'     selection method, error structure, `n_observations`, `n_joinpoints`,
#'     `n_parameters`, `model_df`, `fixed_grid_df`, `inference_df`,
#'     `effective_n`, `sse`, `mse`, `inference_residual_variance`, and
#'     `candidates_available`. SSE is the whitened log-scale residual sum of
#'     squares. MSE is SSE divided by `model_df`.}
#'   \item{coefficients}{A long table with `parameterization` (`"standard"`
#'     or `"general"`), `term`, `estimate`, `se`, `statistic`,
#'     `statistic_type`, `df`, and `p_value`, plus group and inference status.
#'     Standard terms are the first intercept, first slope, and changes in
#'     slope. General terms alternate each segment's intercept and slope.}
#'   \item{vcov, correlation}{Lists named by `group_id`, containing named
#'     general-parameterization matrices in the same interleaved order as the
#'     coefficient table. These are parameter matrices, not observation-error
#'     matrices. Entries are `NULL` when inference is unavailable. Correlations
#'     involving a zero-variance coefficient are `NA`.}
#'   \item{candidates}{One row per successfully fitted candidate count, with
#'     grouping variables, selection method, `n_joinpoints`, formatted
#'     `joinpoints`, observation and parameter counts, `model_df`, `sse`,
#'     `mse`, `bic`, `bic3`, `modified_bic`, `weight`, `weighted_bic`, and
#'     `selected`, and `weight_method` (`"qr"` or `"joinpoint_6.1"`).
#'     Weight is missing for zero joinpoints. The selected flag
#'     follows the actual selector, including DDS and permutation selection;
#'     the displayed BIC-family scores do not replace their decision rules.}
#'   \item{selection}{One row per group with the selected model, selection
#'     method, score, and DDS decision diagnostics when applicable, including
#'     both BIC-selected Delta values, the variance estimate and degrees of
#'     freedom, and the rule used. For permutation selection, Monte Carlo
#'     p-values, interval, resample count, and stopping status are included.}
#'   \item{selection_tests}{The sequential permutation comparisons and their
#'     status, with group identifiers. It has zero rows for non-permutation
#'     methods.}
#'   \item{joinpoints}{The existing numeric joinpoint table from `x`.
#'     No joinpoint confidence intervals are added.}
#' }
#'
#' @details
#' Model selection is not rerun. Candidate summaries and the full parameter
#' covariance are cached when fitting new models. For older serialized fits,
#' missing coefficient covariance is reconstructed at the existing joinpoints;
#' missing candidate history produces a warning and an empty candidate table
#' for that group. Refit the original data to obtain that history.
#' Set `wbic_weight` in [joinpoint()] before fitting to choose the candidate
#' Weight/WBIC numerical method; extraction never changes the cached scores.
#'
#' For `k` joinpoints and `n` observations, `n_parameters = 2*k + 2` and
#' `model_df = n - n_parameters`. The fixed-grid continuous regression has
#' `k + 2` coefficients and its own `fixed_grid_df`. NCI-style coefficient
#' inference removes observations at joinpoints and fits unconstrained segment
#' lines to estimate covariance, giving a third, separate `inference_df`.
#' Coefficient point estimates still come from the continuous fitted model.
#' Intercepts refer to time zero on the log-outcome scale.
#'
#' Failed groups remain in `data` and `statistics`, with missing fitted
#' statistics, no coefficient/candidate rows, and `NULL` matrix entries.
#' If fitting succeeds but inference fails, coefficient estimates are retained
#' with missing inference fields and status `"inference_failed"`.
#' Missing or infinite t statistics follow the conventions of [calc_apc()].
#'
#' The coefficient calculations share the package's covariance-whitening path
#' for correlated errors. Agreement with NCI 6.1.0 has been checked for the
#' bundled independent-error `jp_data` cohort-zero selected model, not for all
#' error structures or all candidate scores. Extraction does not change the
#' existing fitting, EQ, or parametric algorithms. APC/AAPC interval tables
#' remain available through [calc_apc()] and [calc_aapc()].
#'
#' @references
#' National Cancer Institute. Estimated Regression Coefficients (Beta).
#' \url{https://surveillance.cancer.gov/help/joinpoint/statistical-notes/statistics-related-to-the-k-joinpoint-model/estimated-regression-coefficients-beta}.
#'
#' @examples
#' fit <- joinpoint(
#'   subset(jp_data, sex == 0), "year", "rate", se = "se",
#'   max_jp = 1, n_eq = 19
#' )
#' result <- calc_model(fit, ci_method = "Param")
#' result$statistics
#' result$coefficients
#' result$candidates
#' result$vcov[["1"]]
#' @export
calc_model <- function(
  x,
  level = 0.95,
  ci_method = c("EQ", "Param"),
  n_resamples = NULL,
  seed = NULL
) {
  if (!inherits(x, "crtrend_joinpoint")) {
    stop("`x` must be a `crtrend_joinpoint` object.", call. = FALSE)
  }
  ci_method <- normalize_ci_method(match.arg(
    ci_method[1L],
    c("EQ", "Param", "empirical_quantile", "parametric")
  ))
  n_resamples <- n_resamples %||% x$control$n_eq %||% 5001L
  validate_calc_inputs(x, level, ci_method, n_resamples, seed)
  apc <- calc_apc(x, level, ci_method, n_resamples, seed)
  missing_history <- vapply(
    x$models,
    function(model) {
      model$status == "ok" && is.null(model$selection$candidates)
    },
    logical(1)
  )
  if (any(missing_history)) {
    warning(
      "Candidate history is unavailable for some older fits; refit the original data to obtain it.",
      call. = FALSE
    )
  }
  rows <- lapply(x$models, function(model) {
    extract_model_group(
      model,
      apc,
      x$control$selection_method,
      level,
      ci_method
    )
  })
  bind <- function(name) {
    result <- do.call(rbind, lapply(rows, `[[`, name))
    rownames(result) <- NULL
    result
  }
  data <- bind("data")
  data <- data[order(data$row_id), , drop = FALSE]
  rownames(data) <- NULL
  ids <- as.character(vapply(x$models, `[[`, integer(1), "group_id"))
  list(
    data = data,
    statistics = bind("statistics"),
    coefficients = bind("coefficients"),
    vcov = stats::setNames(lapply(rows, `[[`, "vcov"), ids),
    correlation = stats::setNames(lapply(rows, `[[`, "correlation"), ids),
    candidates = bind("candidates"),
    selection = bind("selection"),
    selection_tests = bind("selection_tests"),
    joinpoints = x$joinpoints
  )
}

empty_model_candidates <- function() {
  data.frame(
    n_joinpoints = integer(),
    joinpoints = character(),
    n_observations = integer(),
    n_parameters = integer(),
    model_df = integer(),
    sse = numeric(),
    mse = numeric(),
    bic = numeric(),
    bic3 = numeric(),
    modified_bic = numeric(),
    weight = numeric(),
    weighted_bic = numeric(),
    dds_delta = numeric(),
    dds_score = numeric(),
    selected = logical()
  )
}

joinpoint_candidate_summary <- function(candidates, n, selected_k) {
  rows <- lapply(candidates, function(candidate) {
    k <- candidate$n_joinpoints
    p <- 2L * k + 2L
    data.frame(
      n_joinpoints = k,
      joinpoints = format_joinpoints(candidate$joinpoints),
      n_observations = n,
      n_parameters = p,
      model_df = n - p,
      sse = candidate$sse,
      mse = candidate$sse / (n - p),
      bic = candidate$bic,
      bic3 = candidate$bic3,
      modified_bic = candidate$modified_bic,
      weight = if (k == 0L) NA_real_ else candidate$partial_r_squared,
      weighted_bic = candidate$weighted_bic,
      dds_delta = candidate$dds_delta %||% NA_real_,
      dds_score = candidate$dds_score %||% NA_real_,
      selected = k == selected_k
    )
  })
  if (!length(rows)) {
    return(empty_model_candidates())
  }
  do.call(rbind, rows)
}

model_coefficient_results <- function(model, inference) {
  empty <- data.frame(
    parameterization = character(),
    term = character(),
    estimate = numeric(),
    se = numeric(),
    statistic = numeric(),
    statistic_type = character(),
    df = numeric(),
    p_value = numeric(),
    status = character()
  )
  if (model$status != "ok") {
    return(list(coefficients = empty, vcov = NULL, correlation = NULL))
  }
  k <- length(model$joinpoints)
  b <- as.numeric(stats::coef(model))
  changes <- if (k) b[-c(1L, 2L)] else numeric()
  slopes <- c(b[2L], b[2L] + cumsum(changes))
  intercepts <- c(b[1L], b[1L] - cumsum(changes * model$joinpoints))
  terms <- as.vector(rbind(
    paste("Intercept", seq_len(k + 1L)),
    paste("Slope", seq_len(k + 1L))
  ))
  general <- as.vector(rbind(intercepts, slopes))
  covariance <- inference$coefficient_covariance
  if (!is.null(covariance)) {
    dimnames(covariance) <- list(terms, terms)
  }
  transform <- matrix(0, k + 2L, 2L * (k + 1L))
  transform[1L, 1L] <- transform[2L, 2L] <- 1
  for (j in seq_len(k)) {
    transform[j + 2L, 2L * j] <- -1
    transform[j + 2L, 2L * (j + 1L)] <- 1
  }
  standard_covariance <- if (is.null(covariance)) {
    NULL
  } else {
    transform %*% covariance %*% t(transform)
  }
  table <- function(estimates, covariance, terms, parameterization) {
    n <- length(estimates)
    available <- !is.null(covariance)
    se <- if (available) sqrt(pmax(0, diag(covariance))) else rep(NA_real_, n)
    statistic <- ifelse(se == 0 & estimates == 0, 0, estimates / se)
    df <- if (available) inference$df else NA_real_
    data.frame(
      parameterization,
      term = terms,
      estimate = estimates,
      se,
      statistic,
      statistic_type = if (available) "t" else NA_character_,
      df,
      p_value = 2 * stats::pt(-abs(statistic), df),
      status = if (available) "ok" else "inference_failed",
      row.names = NULL
    )
  }
  correlation <- covariance
  if (!is.null(covariance)) {
    scales <- sqrt(pmax(0, diag(covariance)))
    denominator <- outer(scales, scales)
    correlation <- covariance / denominator
    correlation[denominator == 0] <- NA_real_
  }
  list(
    coefficients = rbind(
      table(
        b,
        standard_covariance,
        c(
          "Intercept 1",
          "Slope 1",
          if (k) paste0("Slope ", 2:(k + 1L), " - Slope ", seq_len(k))
        ),
        "standard"
      ),
      table(general, covariance, terms, "general")
    ),
    vcov = covariance,
    correlation = correlation
  )
}

extract_model_group <- function(
  model,
  apc,
  selection_method,
  level,
  ci_method
) {
  group <- joinpoint_group(model)
  with_group <- function(table) {
    result <- cbind(group[rep(1L, nrow(table)), , drop = FALSE], table)
    rownames(result) <- NULL
    result
  }
  ok <- model$status == "ok"
  analysis <- joinpoint_analysis_data(model)
  inference <- if (ok) cached_joinpoint_inference(model) else NULL
  if (ok && is.null(inference$coefficient_covariance)) {
    inference <- nci_joinpoint_inference(model)
  }
  coefficients <- model_coefficient_results(model, inference)
  candidates <- model$selection$candidates %||% empty_model_candidates()
  n <- nrow(model$model)
  k <- if (ok) model$n_joinpoints else NA_integer_
  p <- 2L * k + 2L
  sse <- if (ok) {
    model_sse(joinpoint_engine(model), analysis, joinpoint_is_weighted(model))
  } else {
    NA_real_
  }
  statistics <- with_group(data.frame(
    status = model$status,
    inference_status = if (!ok) {
      model$status
    } else if (is.null(inference)) {
      "inference_failed"
    } else {
      "ok"
    },
    selection_method,
    error = model$control$error,
    n_observations = n,
    n_joinpoints = k,
    n_parameters = p,
    model_df = n - p,
    fixed_grid_df = if (ok) model$df.residual else NA_integer_,
    inference_df = inference$df %||% NA_integer_,
    effective_n = inference$effective_n %||% NA_integer_,
    sse,
    mse = sse / (n - p),
    inference_residual_variance = inference$residual_variance %||% NA_real_,
    candidates_available = ok && !is.null(model$selection$candidates)
  ))
  time <- analysis$trend_time
  joinpoint <- match(time, model$joinpoints)
  segment <- if (ok) {
    findInterval(time, model$joinpoints) + 1L
  } else {
    rep(NA_integer_, n)
  }
  segment[!is.na(joinpoint)] <- NA_integer_
  group_apc <- apc[apc$group_id == model$group_id, , drop = FALSE]
  index <- match(segment, group_apc$segment)
  index[is.na(segment)] <- NA_integer_
  data <- with_group(data.frame(
    row_id = analysis$trend_row_id,
    status = model$status,
    time,
    observed = analysis$trend_outcome,
    se = analysis$trend_se %||% rep(NA_real_, n),
    fitted = as.numeric(stats::fitted(model)),
    residual = as.numeric(stats::residuals(model)),
    linear_predictor = as.numeric(model$linear.predictors),
    link_residual = as.numeric(model$link.residuals),
    joinpoint,
    segment,
    apc = group_apc$apc[index],
    apc_p_value = group_apc$p_value[index],
    apc_significant = group_apc$p_value[index] < (1 - level),
    apc_method = ci_method
  ))
  candidates$selection_method <- rep(selection_method, nrow(candidates))
  candidates$weight_method <- rep(
    model$control$wbic_weight %||% "qr",
    nrow(candidates)
  )
  selection <- with_group(data.frame(
    selection_method = selection_method,
    selected_n_joinpoints = k,
    selection_score = model$selection$score %||% NA_real_,
    dds_method = model$selection$dds_method %||% NA_character_,
    dds_decision_rule = model$selection$dds_decision_rule %||% NA_character_,
    dds_bic_delta = model$selection$dds_bic_delta %||% NA_real_,
    dds_bic3_delta = model$selection$dds_bic3_delta %||% NA_real_,
    dds_delta_min = model$selection$dds_delta_min %||% NA_real_,
    dds_delta_max = model$selection$dds_delta_max %||% NA_real_,
    dds_sigma_squared = model$selection$dds_sigma_squared %||% NA_real_,
    dds_variance_df = model$selection$dds_variance_df %||% NA_integer_,
    selected_joinpoints = format_joinpoints(model$joinpoints),
    p_value = model$selection$p_value %||% NA_real_,
    p_value_lower = model$selection$p_value_lower %||% NA_real_,
    p_value_upper = model$selection$p_value_upper %||% NA_real_,
    n_permutations = model$selection$n_permutations %||% NA_integer_,
    stopped_early = model$selection$stopped_early %||% NA
  ))
  selection_tests <- model$selection$tests %||% data.frame()
  if (nrow(selection_tests) > 0L) {
    selection_tests <- with_group(selection_tests)
  } else {
    selection_tests <- data.frame(
      group_id = integer(),
      selection_tests,
      check.names = FALSE
    )
    if (nrow(group) > 0L) {
      for (column in names(group)) {
        selection_tests[[column]] <- group[[column]][FALSE]
      }
      selection_tests <- selection_tests[c(
        "group_id",
        names(group),
        setdiff(names(selection_tests), c("group_id", names(group)))
      )]
    }
  }
  list(
    data = data,
    statistics = statistics,
    coefficients = with_group(coefficients$coefficients),
    vcov = coefficients$vcov,
    correlation = coefficients$correlation,
    candidates = with_group(candidates),
    selection = selection,
    selection_tests = selection_tests
  )
}
