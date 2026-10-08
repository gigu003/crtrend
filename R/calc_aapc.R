#' Calculate average annual percent change over a specified interval
#'
#' @param x A [joinpoint()] result.
#' @param start Integer start year of the requested interval.
#' @param end Integer end year of the requested interval.
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
#' @return A data frame with one row per group containing the requested
#'   interval, duration-weighted AAPC, its standard error, confidence interval,
#'   and the two-sided p-value for a zero weighted log-slope. The weighting is
#'   performed on the log-slope scale before transformation to percent change,
#'   following Clegg et al. (2009). The `start` and `end` columns are integer
#'   years. For `ci_method = "Param"`, `statistic` is the signed weighted
#'   log-slope divided by its standard error. `statistic_type` is `"t"`
#'   when the interval is within one fitted segment, and `"z"` when it spans
#'   multiple segments. Both columns are `NA` for EQ or unavailable inference.
#'
#' @details
#' An EQ distribution with all weighted log-slope draws exactly zero returns
#' an `NA` p-value, retaining its interval, as a degeneracy safeguard rather
#' than an externally benchmarked NCI convention. Parametric extraction does
#' not change R's random-number state. Seeded EQ uses its own generator;
#' unseeded, uncached EQ consumes R's current stream.
#'
#' For an interval spanning multiple fitted segments, the segment log-slopes
#' are weighted by their overlapping durations. The standard error is the
#' quadratic form of those weights and the complete segment-slope covariance
#' matrix, as in Clegg et al. (2009). The confidence interval is formed on the
#' log-slope scale and then exponentiated. It uses a standard-normal critical
#' value for multiple segments, and a residual-degrees-of-freedom t critical
#' value when the requested interval lies within one fitted segment, matching
#' Joinpoint's parametric AAPC convention. For the default empirical-quantile
#' method,
#' residual-resampled models are refitted with the selected number of
#' joinpoints. Joinpoint 6.1.0 includes the original fitted estimate in the EQ
#' distribution; the two-sided p-value is its empirical sign probability on
#' the log-AAPC scale. When the requested number of resamples matches
#' `joinpoint()`'s `n_eq` value and no `seed` is supplied, these draws are
#' reused from the fitted object. For a different number of resamples,
#' `joinpoint()`'s `eq_seed` is used unless `seed` is supplied.
#' EQ uses the same version-specific residual simulation as [calc_apc()].
#' Saved independent-error, zero- and one-joinpoint 6.1.0 audit exports and a
#' two-joinpoint DDS case are reproduced at their reported precision; other
#' versions and correlated
#' error models have not been externally benchmarked.
#' For a zero log-scale standard error, the parametric statistic is zero
#' when the weighted log-slope is zero and signed infinity otherwise,
#' consistent with the existing p-value convention of one and zero,
#' respectively. This statistic is not `aapc / aapc_se`, because those
#' returned columns are on the percent-change scale.
#'
#' @references
#' Clegg, L. X., Hankey, B. F., Tiwari, R., Feuer, E. J., and Edwards, B. K.
#' (2009). Estimating average annual percent change in trend analysis.
#' *Statistics in Medicine*, 28(29), 3670-3682.
#'
#' Kim, H.-J., Luo, J., Chen, H.-S., Green, D., Buckman, D., Byrne, J., and
#' Feuer, E. J. (2017). Improved confidence interval for average annual
#' percent change in trend analysis. *Statistics in Medicine*, 36(19),
#' 3059-3074.
#'
#' National Cancer Institute. Empirical Quantile Confidence Interval.
#' Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/setting-parameters/method-and-parameters-tab/apc-aapc-tau-confidence-intervals/empirical-quantile}.
#' @export
calc_aapc <- function(
  x,
  start,
  end,
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
  validate_interval(start, end)
  start <- as.integer(start)
  end <- as.integer(end)
  resampling_seed <- resolve_empirical_seed(x, n_resamples, seed)
  result <- do.call(
    rbind,
    lapply(
      x$models,
      calc_aapc_group,
      start = start,
      end = end,
      level = level,
      ci_method = ci_method,
      n_resamples = n_resamples,
      seed = resampling_seed
    )
  )
  rownames(result) <- NULL
  result
}

#' Internal helper: Calculate AAPC results for one fitted group.
#'
#' @title Internal joinpoint implementation helpers
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Clegg, L. X., Hankey, B. F., Tiwari, R., Feuer, E. J., and Edwards, B. K.
#' (2009). Estimating average annual percent change in trend analysis.
#' *Statistics in Medicine*, 28(29), 3670-3682.
#'
#' Kim, H.-J., Luo, J., Chen, H.-S., Green, D., Buckman, D., Byrne, J., and
#' Feuer, E. J. (2017). Improved confidence interval for average annual
#' percent change in trend analysis. *Statistics in Medicine*, 36(19),
#' 3059-3074.
calc_aapc_group <- function(
  fit,
  start,
  end,
  level,
  ci_method,
  n_resamples,
  seed
) {
  result <- joinpoint_group(fit)
  result$status <- fit$status
  result$n_joinpoints <- fit$n_joinpoints
  result$joinpoints <- format_joinpoints(fit$joinpoints)
  result$start <- as.integer(start)
  result$end <- as.integer(end)
  result$aapc <- NA_real_
  result$aapc_se <- NA_real_
  result$aapc_lower <- NA_real_
  result$aapc_upper <- NA_real_
  result$p_value <- NA_real_
  result$statistic <- NA_real_
  result$statistic_type <- NA_character_
  if (fit$status != "ok") {
    return(result)
  }
  if (start < fit$start || end > fit$end) {
    result$status <- "interval_out_of_range"
    return(result)
  }
  inference <- cached_joinpoint_inference(fit)
  if (is.null(inference)) {
    result$status <- "inference_failed"
    return(result)
  }
  weights <- aapc_weights(
    joinpoints = fit$joinpoints,
    fit_start = fit$start,
    fit_end = fit$end,
    start = start,
    end = end
  )
  aapc_log <- sum(weights * inference$slope_estimates)
  aapc_variance <- drop(
    crossprod(weights, inference$slope_covariance %*% weights)
  )
  aapc_se <- sqrt(max(0, aapc_variance))
  positive_weights <- weights > 0
  critical_value <- if (sum(positive_weights) == 1L) {
    stats::qt((1 + level) / 2, df = inference$df)
  } else {
    stats::qnorm((1 + level) / 2)
  }
  result$status <- "ok"
  result$aapc <- 100 * (exp(aapc_log) - 1)
  result$aapc_se <- 100 * exp(aapc_log) * aapc_se
  if (ci_method == "EQ") {
    aapc_draws <- empirical_aapc_draws(
      fit,
      start = start,
      end = end,
      n_resamples = n_resamples,
      seed = seed
    )
    aapc_interval <- empirical_interval_matrix(
      matrix(aapc_draws, ncol = 1L),
      level
    )[, 1L]
    result$aapc_lower <- 100 * (exp(aapc_interval[1L]) - 1)
    result$aapc_upper <- 100 * (exp(aapc_interval[2L]) - 1)
    result$p_value <- empirical_two_sided_p(aapc_draws)
  } else {
    result$aapc_lower <- 100 *
      (exp(aapc_log - critical_value * aapc_se) - 1)
    result$aapc_upper <- 100 *
      (exp(aapc_log + critical_value * aapc_se) - 1)
    result$p_value <- if (aapc_se > 0 && sum(positive_weights) == 1L) {
      2 * stats::pt(-abs(aapc_log / aapc_se), df = inference$df)
    } else if (aapc_se > 0) {
      2 * stats::pnorm(-abs(aapc_log / aapc_se))
    } else if (aapc_log == 0) {
      1
    } else {
      0
    }
    result$statistic <- if (aapc_se == 0 && aapc_log == 0) {
      0
    } else {
      aapc_log / aapc_se
    }
    result$statistic_type <- if (sum(positive_weights) == 1L) "t" else "z"
  }
  result
}
