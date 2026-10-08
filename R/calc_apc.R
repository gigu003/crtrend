#' Calculate annual percent change for each fitted segment
#'
#' @param x A [joinpoint()] result.
#' @param level Confidence level. The default is `0.95`.
#' @param ci_method Confidence interval method. The default, `"EQ"`, is the
#'   Empirical Quantile method targeting NCI Joinpoint 6.1.0.
#'   `"Param"` uses the parametric method. The full values
#'   `"empirical_quantile"` and `"parametric"` are accepted as aliases.
#' @param n_resamples Number of EQ estimates for `ci_method = "EQ"`, including
#'   the original estimate and `n_resamples - 1` residual-resampled refits.
#'   `NULL` (default) uses the fitted object's `n_eq`, or `5001` for older fits.
#' @param seed Optional seed for the residual resampling.
#'   Must be NULL or a finite integer in `[-2147483647, 2147483647]`.
#'   Resample counts must be integers in `[2, 2147483646]`.
#'
#' @return A data frame with one row per fitted segment. It contains the
#'   segment boundaries, APC, confidence interval, and the two-sided p-value
#'   for a zero log-scale slope. The `start` and `end` columns are integer
#'   years. For `ci_method = "Param"`, `statistic` is the signed t statistic
#'   (log-scale slope divided by its standard error), and `statistic_type`
#'   is `"t"`. Both columns are `NA` for EQ or unavailable inference.
#'
#' @details
#' An EQ distribution with all log-slope draws exactly zero has an undefined
#' sign test: its p-value is `NA`, while its interval is retained. This is a
#' degeneracy safeguard, not an externally benchmarked NCI convention.
#' Parametric extraction does not change R's random-number state. Seeded EQ
#' uses its own generator; unseeded, uncached EQ consumes R's current stream.
#'
#' With the default empirical-quantile method, the selected number of
#' joinpoints is held fixed while each residual resample is refitted by the
#' integer-year grid search. Joinpoint 6.1.0 includes the original fitted
#' estimate as the first value of the EQ distribution. The confidence limits
#' are empirical quantiles of these estimates, and the p-value is the two-sided
#' empirical sign probability on the log-slope scale. When the requested number of
#' resamples matches `joinpoint()`'s `n_eq` value and no `seed` is supplied,
#' these draws are reused from the fitted object. For a different number of
#' resamples, `joinpoint()`'s `eq_seed` is used unless `seed` is supplied.
#' With a non-NULL seed, EQ uses the Joinpoint 6.1.0-compatible uniform
#' generator, robust residual endpoints, and one-based Method 2 interpolation
#' convention. The implementation matches saved 6.1.0 exports at their
#' reported precision for the bundled independent-error, zero- and
#' one-joinpoint audit cases and the two-joinpoint DDS reference case. This does
#' not establish equivalence for other
#' Joinpoint versions or for correlated-error models.
#' For a zero standard error, the parametric statistic is zero when the
#' slope is zero and signed infinity otherwise, consistent with the existing
#' p-value convention of one and zero, respectively.
#'
#' @references
#' Kim, H.-J., Luo, J., Chen, H.-S., Green, D., Buckman, D., Byrne, J., and
#' Feuer, E. J. (2017). Improved confidence interval for average annual
#' percent change in trend analysis. *Statistics in Medicine*, 36(19),
#' 3059-3074.
#'
#' Kim, H.-J., Chen, H.-S., Byrne, J., Wheeler, B., and Feuer, E. J. (2022).
#' Twenty years since Joinpoint 1.0: Two major enhancements, their
#' justification, and impact. *Statistics in Medicine*, 41(16), 3102-3130.
#'
#' National Cancer Institute. Empirical Quantile Confidence Interval.
#' Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/setting-parameters/method-and-parameters-tab/apc-aapc-tau-confidence-intervals/empirical-quantile}.
#' @export
calc_apc <- function(
  x,
  level = 0.95,
  ci_method = c("EQ", "Param"),
  n_resamples = NULL,
  seed = NULL
) {
  ci_method <- match.arg(
    ci_method[1L],
    c("EQ", "Param", "empirical_quantile", "parametric")
  )
  ci_method <- normalize_ci_method(ci_method)
  n_resamples <- resolve_calc_resamples(x, n_resamples)
  validate_calc_inputs(x, level, ci_method, n_resamples, seed)
  resampling_seed <- resolve_empirical_seed(x, n_resamples, seed)
  result <- do.call(
    rbind,
    lapply(
      x$models,
      calc_apc_group,
      level = level,
      ci_method = ci_method,
      n_resamples = n_resamples,
      seed = resampling_seed
    )
  )
  rownames(result) <- NULL
  result
}

#' Internal helper: Calculate APC results for one fitted group.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Luo, J., Chen, H.-S., Green, D., Buckman, D., Byrne, J., and
#' Feuer, E. J. (2017). Improved confidence interval for average annual
#' percent change in trend analysis. *Statistics in Medicine*, 36(19),
#' 3059-3074.
#'
#' National Cancer Institute. Empirical Quantile Confidence Interval.
#' Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/setting-parameters/method-and-parameters-tab/apc-aapc-tau-confidence-intervals/empirical-quantile}.
calc_apc_group <- function(fit, level, ci_method, n_resamples, seed) {
  if (fit$status != "ok") {
    return(empty_apc_result(fit))
  }
  inference <- cached_joinpoint_inference(fit)
  if (is.null(inference)) {
    fit$status <- "inference_failed"
    return(empty_apc_result(fit))
  }
  critical_value <- stats::qt((1 + level) / 2, df = inference$df)
  start <- c(fit$start, fit$joinpoints)
  end <- c(fit$joinpoints, fit$end)
  slope <- inference$slope_estimates
  slope_se <- inference$slope_se
  p_value <- vapply(
    seq_along(slope),
    function(i) {
      if (slope_se[i] > 0) {
        2 * stats::pt(-abs(slope[i] / slope_se[i]), df = inference$df)
      } else if (slope[i] == 0) {
        1
      } else {
        0
      }
    },
    numeric(1)
  )
  group <- joinpoint_group(fit)
  result <- group[rep(1L, length(start)), , drop = FALSE]
  result$status <- "ok"
  result$n_joinpoints <- fit$n_joinpoints
  result$joinpoints <- format_joinpoints(fit$joinpoints)
  result$segment <- seq_along(start)
  result$period <- paste0(
    format_period_value(start),
    "-",
    format_period_value(end)
  )
  result$start <- as.integer(round(start))
  result$end <- as.integer(round(end))
  result$apc <- 100 * (exp(slope) - 1)
  if (ci_method == "EQ") {
    slope_draws <- empirical_slope_draws(
      fit,
      n_resamples = n_resamples,
      seed = seed
    )
    slope_interval <- empirical_interval_matrix(slope_draws, level)
    result$apc_lower <- 100 * (exp(slope_interval[1L, ]) - 1)
    result$apc_upper <- 100 * (exp(slope_interval[2L, ]) - 1)
    result$p_value <- empirical_two_sided_p(slope_draws)
    result$statistic <- NA_real_
    result$statistic_type <- NA_character_
  } else {
    result$apc_lower <- 100 * (exp(slope - critical_value * slope_se) - 1)
    result$apc_upper <- 100 * (exp(slope + critical_value * slope_se) - 1)
    result$p_value <- p_value
    result$statistic <- ifelse(slope_se == 0 & slope == 0, 0, slope / slope_se)
    result$statistic_type <- "t"
  }
  result
}

#' Internal helper: Create an empty APC result for a failed group.
#'
#' @rdname crtrend-internals
#' @keywords internal
empty_apc_result <- function(fit) {
  result <- joinpoint_group(fit)
  result$status <- fit$status
  result$n_joinpoints <- fit$n_joinpoints
  result$joinpoints <- format_joinpoints(fit$joinpoints)
  result$segment <- NA_integer_
  result$period <- paste0(
    format_period_value(fit$start),
    "-",
    format_period_value(fit$end)
  )
  result$start <- as.integer(round(fit$start))
  result$end <- as.integer(round(fit$end))
  result$apc <- NA_real_
  result$apc_lower <- NA_real_
  result$apc_upper <- NA_real_
  result$p_value <- NA_real_
  result$statistic <- NA_real_
  result$statistic_type <- NA_character_
  result
}
