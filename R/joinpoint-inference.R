#' Internal helper: Normalize the APC/AAPC confidence-interval method.
#'
#' @rdname crtrend-internals
#' @keywords internal
normalize_ci_method <- function(ci_method) {
  switch(
    ci_method,
    EQ = "EQ",
    Param = "Param",
    empirical_quantile = "EQ",
    parametric = "Param"
  )
}

#' Internal helper: Reuse coefficient inference cached by `joinpoint()`.
#'
#' @rdname crtrend-internals
#' @keywords internal
cached_joinpoint_inference <- function(fit) {
  fit$inference %||% nci_joinpoint_inference(fit)
}

#' Internal helper: Calculate empirical confidence limits for resampled estimates.
#'
#' @rdname crtrend-internals
#' @keywords internal
empirical_interval_matrix <- function(draws, level) {
  if (is.null(draws)) {
    return(matrix(NA_real_, nrow = 2L, ncol = 0L))
  }
  quantiles <- c((1 - level) / 2, (1 + level) / 2)
  vapply(
    seq_len(ncol(draws)),
    function(index) {
      values <- draws[, index]
      if (any(!is.finite(values))) {
        return(rep(NA_real_, 2L))
      }
      stats::quantile(
        values,
        probs = quantiles,
        names = FALSE,
        type = 7L
      )
    },
    numeric(2L)
  )
}

#' Internal helper: Calculate the NCI empirical two-sided sign p-value.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' National Cancer Institute. Empirical Quantile Confidence Interval.
#' Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/setting-parameters/method-and-parameters-tab/apc-aapc-tau-confidence-intervals/empirical-quantile}.
empirical_two_sided_p <- function(draws) {
  if (is.null(draws)) {
    return(NA_real_)
  }
  if (is.null(dim(draws))) {
    draws <- matrix(draws, ncol = 1L)
  }
  vapply(
    seq_len(ncol(draws)),
    function(index) {
      values <- draws[, index]
      if (length(values) == 0L || any(!is.finite(values))) {
        return(NA_real_)
      }
      if (all(values == 0)) {
        return(NA_real_)
      }
      2 * min(sum(values < 0), sum(values >= 0)) / length(values)
    },
    numeric(1L)
  )
}

#' Internal helper: Generate resampled segment slope estimates.
#'
#' @rdname crtrend-internals
#' @keywords internal
empirical_slope_draws <- function(fit, n_resamples, seed = NULL) {
  draws <- cached_empirical_model_draws(fit, n_resamples, seed)
  if (is.null(draws)) {
    return(NULL)
  }
  draws$slopes
}

#' Internal helper: Generate resampled AAPC log-slope estimates.
#'
#' @rdname crtrend-internals
#' @keywords internal
empirical_aapc_draws <- function(fit, start, end, n_resamples, seed = NULL) {
  draws <- cached_empirical_model_draws(fit, n_resamples, seed)
  if (is.null(draws)) {
    return(numeric())
  }
  result <- rep(NA_real_, n_resamples)
  for (iteration in seq_len(n_resamples)) {
    joinpoints <- draws$joinpoints[iteration, , drop = TRUE]
    weights <- aapc_weights(
      joinpoints = joinpoints,
      fit_start = fit$start,
      fit_end = fit$end,
      start = start,
      end = end
    )
    result[iteration] <- sum(weights * draws$slopes[iteration, ])
  }
  result
}

#' Internal helper: Resolve the seed for uncached empirical resamples.
#'
#' @rdname crtrend-internals
#' @keywords internal
resolve_empirical_seed <- function(x, n_resamples, seed) {
  if (!is.null(seed)) {
    return(seed)
  }
  if (
    !is.null(x$control) &&
      !identical(as.integer(n_resamples), as.integer(x$control$n_eq))
  ) {
    return(x$control$eq_seed %||% NULL)
  }
  NULL
}

#' Internal helper: Reuse empirical draws cached by `joinpoint()` when possible.
#'
#' @rdname crtrend-internals
#' @keywords internal
cached_empirical_model_draws <- function(fit, n_resamples, seed = NULL) {
  resamples <- fit$resamples %||%
    list(
      n = fit$eq_n_resamples,
      draws = fit$eq_draws
    )
  if (
    is.null(seed) &&
      identical(as.integer(n_resamples), resamples$n) &&
      !is.null(resamples$draws)
  ) {
    return(resamples$draws)
  }
  nci_empirical_model_draws(fit, n_resamples, seed = seed)
}

#' Internal helper: Refit the selected joinpoint count for each empirical resample.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Luo, J., Chen, H.-S., Green, D., Buckman, D., Byrne, J., and
#' Feuer, E. J. (2017). Improved confidence interval for average annual
#' percent change in trend analysis. *Statistics in Medicine*, 36(19),
#' 3059-3074.
#'
#' Kim, H.-J., Chen, H.-S., Byrne, J., Wheeler, B., and Feuer, E. J. (2022).
#' Twenty years since Joinpoint 1.0: Two major enhancements, their
#' justification, and impact. *Statistics in Medicine*, 41(16), 3102-3130.
nci_empirical_model_draws <- function(fit, n_resamples, seed = NULL) {
  # Joinpoint 6.1.0 counts the original estimate as the first EQ value.
  # The remaining values come from residual-resampled model refits.
  n_boot <- n_resamples - 1L
  responses <- nci_resampled_responses(fit, n_boot, seed = seed)
  if (is.null(responses)) {
    return(NULL)
  }
  data <- joinpoint_analysis_data(fit)
  transformed_response <- joinpoint_whiten(data, responses)
  grid <- if (fit$n_joinpoints == 0L) {
    list(numeric())
  } else {
    integer_joinpoint_grid(
      data$trend_time,
      fit$n_joinpoints,
      min_points = joinpoint_min_observations(fit)
    )
  }
  if (is.null(grid)) {
    return(NULL)
  }
  grid_solver <- compile_joinpoint_grid_solver(grid, data, fit$weighted)
  if (is.null(grid_solver)) {
    return(NULL)
  }
  best_sse <- rep(Inf, n_boot)
  best_index <- rep(NA_integer_, n_boot)
  batch_size <- max(
    1L,
    min(
      128L,
      as.integer(floor(2e6 / ncol(grid_solver$basis)))
    )
  )
  first <- 1L
  while (first <= n_boot) {
    last <- min(n_boot, first + batch_size - 1L)
    columns <- first:last
    projection <- grid_projection_minima(
      grid_solver,
      transformed_response[, columns, drop = FALSE]
    )
    best_sse[columns] <- projection$sse
    valid <- !is.na(projection$model_index)
    best_index[columns[valid]] <- grid_solver$model_indices[
      projection$model_index[valid]
    ]
    first <- last + 1L
  }
  if (all(is.na(best_index))) {
    return(NULL)
  }
  n_segments <- fit$n_joinpoints + 1L
  slopes <- matrix(NA_real_, nrow = n_resamples, ncol = n_segments)
  selected_joinpoints <- matrix(
    NA_real_,
    nrow = n_resamples,
    ncol = fit$n_joinpoints
  )
  slopes[1L, ] <- cached_joinpoint_inference(fit)$slope_estimates
  if (fit$n_joinpoints > 0L) {
    selected_joinpoints[1L, ] <- fit$joinpoints
  }
  for (index in unique(best_index[!is.na(best_index)])) {
    selected <- which(best_index == index)
    candidate_joinpoints <- grid[[index]]
    design <- joinpoint_design(candidate_joinpoints, data)
    solver <- qr(joinpoint_whiten(data, design))
    candidate_slopes <- nci_grid_resampled_slopes(
      solver = solver,
      transformed_response = transformed_response[, selected, drop = FALSE],
      n_joinpoints = fit$n_joinpoints
    )
    slopes[selected + 1L, ] <- candidate_slopes
    if (fit$n_joinpoints > 0L) {
      selected_joinpoints[selected + 1L, ] <- matrix(
        candidate_joinpoints,
        nrow = length(selected),
        ncol = fit$n_joinpoints,
        byrow = TRUE
      )
    }
  }
  list(
    slopes = slopes,
    joinpoints = selected_joinpoints
  )
}

#' Internal helper: Extract segment slopes from a resampled grid-search fit.
#'
#' @rdname crtrend-internals
#' @keywords internal
nci_grid_resampled_slopes <- function(
  solver,
  transformed_response,
  n_joinpoints
) {
  coefficients <- tryCatch(
    as.matrix(qr.coef(solver, transformed_response)),
    error = function(error) NULL
  )
  if (
    is.null(coefficients) ||
      nrow(coefficients) < n_joinpoints + 2L ||
      any(!is.finite(coefficients))
  ) {
    return(matrix(
      NA_real_,
      nrow = ncol(transformed_response),
      ncol = n_joinpoints + 1L
    ))
  }
  baseline <- coefficients[2L, ]
  if (n_joinpoints == 0L) {
    return(matrix(baseline, ncol = 1L))
  }
  delta <- coefficients[seq.int(3L, n_joinpoints + 2L), , drop = FALSE]
  cumulative <- matrix(
    0,
    nrow = ncol(delta),
    ncol = n_joinpoints
  )
  for (index in seq_len(ncol(delta))) {
    cumulative[index, ] <- cumsum(delta[, index])
  }
  cbind(
    matrix(baseline, ncol = 1L),
    matrix(baseline, nrow = length(baseline), ncol = n_joinpoints) +
      cumulative
  )
}

#' Internal helper: Generate residual-resampled log-scale responses.
#'
#' @rdname crtrend-internals
#' @keywords internal
nci_resampled_responses <- function(fit, n_resamples, seed = NULL) {
  data <- joinpoint_analysis_data(fit)
  engine <- joinpoint_engine(fit)
  if (is.null(data) || is.null(engine)) {
    return(NULL)
  }
  residuals <- stats::residuals(engine)
  fitted <- stats::fitted(engine)
  residuals <- as.numeric(joinpoint_whiten(data, residuals))
  if (any(!is.finite(c(residuals, fitted)))) {
    return(NULL)
  }
  sampled_residuals <- matrix(
    nci_empirical_residuals(
      residuals,
      nrow(data) * n_resamples,
      seed = seed
    ),
    nrow = nrow(data),
    ncol = n_resamples
  )
  sampled_residuals <- joinpoint_unwhiten(data, sampled_residuals)
  matrix(fitted, nrow = nrow(data), ncol = n_resamples) +
    sampled_residuals
}

#' Internal helper: Generate NCI EmpQ Method 2 residual quantiles.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Luo, J., Chen, H.-S., Green, D., Buckman, D., Byrne, J., and
#' Feuer, E. J. (2017). Improved confidence interval for average annual
#' percent change in trend analysis. *Statistics in Medicine*, 36(19),
#' 3059-3074.
#'
#' Kim, H.-J., Chen, H.-S., Byrne, J., Wheeler, B., and Feuer, E. J. (2022).
#' Twenty years since Joinpoint 1.0: Two major enhancements, their
#' justification, and impact. *Statistics in Medicine*, 41(16), 3102-3130.
nci_empirical_residuals <- function(residuals, size, seed = NULL) {
  residuals <- sort(as.numeric(residuals))
  n <- length(residuals)
  if (n == 0L || any(!is.finite(residuals))) {
    return(rep(NA_real_, size))
  }
  endpoints <- nci_empirical_endpoints(residuals)
  knots <- c(endpoints, residuals)
  knots <- sort(knots)
  if (!is.null(seed)) {
    uniform <- nci_mt_uniforms(size, seed)
  } else {
    uniform <- stats::runif(size)
  }
  interval <- pmin(floor(uniform * (n + 1L)), n)
  # The 6.1.0 SimulatingEmpiricalCDF Method 2 branch subtracts the
  # one-based interval index. Preserve this version-specific convention.
  fraction <- uniform * (n + 1L) - (interval + 1L)
  knots[interval + 1L] +
    (knots[interval + 2L] - knots[interval + 1L]) * fraction
}

#' Internal helper: Reproduce Joinpoint 6.1.0's robust EmpQ endpoints.
#'
#' The 6.1.0 executable uses the EmpQ2 endpoint rule described by Kim et al.
#' (2022), with final limits based on the residual median and 5.189 times the
#' unscaled MAD. Observed extremes at or beyond those limits receive the
#' `2 * IQR / n` extension. The preliminary log-log fence in the executable
#' is overwritten before interpolation and does not affect the final limits.
#'
#' @rdname crtrend-internals
#' @keywords internal
nci_empirical_endpoints <- function(residuals) {
  residuals <- sort(as.numeric(residuals))
  n <- length(residuals)
  spread <- stats::IQR(residuals, type = 7L)
  center <- stats::median(residuals)
  mad <- stats::median(abs(residuals - center))
  lower_candidate <- center - 5.189 * mad
  upper_candidate <- center + 5.189 * mad
  lower <- if (lower_candidate >= residuals[1L]) {
    residuals[1L] - 2 * spread / n
  } else {
    lower_candidate
  }
  upper <- if (residuals[n] >= upper_candidate) {
    residuals[n] + 2 * spread / n
  } else {
    upper_candidate
  }
  c(lower, upper)
}

#' Internal helper: Generate Joinpoint 6.1.0's MT19937-compatible uniforms.
#'
#' Joinpoint 6.1.0 uses a 624-word generator with a legacy seed expansion and
#' the standard MT19937 twist.  Keeping it local avoids changing R's global RNG
#' state and makes the EQ stream reproducible across R sessions.
#'
#' @rdname crtrend-internals
#' @keywords internal
nci_mt_uniforms <- function(size, seed) {
  nci_mt_stream(seed)(size)
}

nci_mt_stream <- function(seed) {
  state <- numeric(624L)
  current <- nci_uint32(seed)
  for (index in seq_len(624L)) {
    product <- (current * 0x10dcd + 1) %% 2^32
    state[index] <- nci_u32_or(
      nci_u32_and(current, 4294901760),
      floor(product / 2^16)
    )
    current <- (product * 0x10dcd + 1) %% 2^32
  }
  buffer <- numeric()
  # Each block reads only old words or words already updated by prior blocks.
  blocks <- list(1:227, 228:454, 455:623, 624L)
  function(size) {
    output <- numeric(size)
    written <- 0L
    while (written < size) {
      if (length(buffer) > 0L) {
        count <- min(length(buffer), size - written)
        output[written + seq_len(count)] <- buffer[seq_len(count)]
        buffer <<- buffer[-seq_len(count)]
        written <- written + count
        next
      }
      next_state <- state
      for (position in blocks) {
        next_position <- position %% 624L + 1L
        value <- nci_u32_or(
          nci_u32_and(next_state[position], 2147483648),
          nci_u32_and(next_state[next_position], 2147483647)
        )
        updated <- nci_u32_xor(
          next_state[(position + 396L) %% 624L + 1L],
          floor(value / 2)
        )
        odd <- value %% 2 != 0
        if (any(odd)) {
          updated[odd] <- nci_u32_xor(updated[odd], 2567483615)
        }
        next_state[position] <- updated
      }
      state <<- next_state
      value <- state
      value <- nci_u32_xor(value, floor(value / 2^11))
      value <- nci_u32_xor(
        value,
        (nci_u32_and(value, 4282013869) * 2^7) %% 2^32
      )
      value <- nci_u32_xor(
        value,
        (nci_u32_and(value, 4294958988) * 2^15) %% 2^32
      )
      value <- nci_u32_xor(value, floor(value / 2^18))
      buffer <<- value / 2^32
    }
    output
  }
}

nci_permutation_shuffler <- function(residuals, uniform_stream) {
  n <- length(residuals)
  function(size) {
    draws <- matrix(uniform_stream((n - 1L) * size), nrow = n - 1L)
    output <- matrix(NA_real_, n, size)
    current <- residuals
    for (draw in seq_len(size)) {
      for (i in seq_len(n - 1L)) {
        j <- i + floor(draws[i, draw] * (n - i + 1L))
        temporary <- current[i]
        current[i] <- current[j]
        current[j] <- temporary
      }
      output[, draw] <- current
    }
    residuals <<- current
    output
  }
}

#' @rdname crtrend-internals
#' @keywords internal
nci_uint32 <- function(value) {
  as.double(value) %% 2^32
}

#' @rdname crtrend-internals
#' @keywords internal
nci_u32_and <- function(left, right) {
  low <- bitwAnd(
    as.integer(left %% 2^16),
    as.integer(right %% 2^16)
  )
  high <- bitwAnd(
    as.integer(floor(left / 2^16)),
    as.integer(floor(right / 2^16))
  )
  as.double(low) + 2^16 * as.double(high)
}

#' @rdname crtrend-internals
#' @keywords internal
nci_u32_or <- function(left, right) {
  low <- bitwOr(
    as.integer(left %% 2^16),
    as.integer(right %% 2^16)
  )
  high <- bitwOr(
    as.integer(floor(left / 2^16)),
    as.integer(floor(right / 2^16))
  )
  as.double(low) + 2^16 * as.double(high)
}

#' @rdname crtrend-internals
#' @keywords internal
nci_u32_xor <- function(left, right) {
  low <- bitwXor(
    as.integer(left %% 2^16),
    as.integer(right %% 2^16)
  )
  high <- bitwXor(
    as.integer(floor(left / 2^16)),
    as.integer(floor(right / 2^16))
  )
  as.double(low) + 2^16 * as.double(high)
}

#' Internal helper: Calculate duration weights for an AAPC interval.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Clegg, L. X., Hankey, B. F., Tiwari, R., Feuer, E. J., and Edwards, B. K.
#' (2009). Estimating average annual percent change in trend analysis.
#' *Statistics in Medicine*, 28(29), 3670-3682.
aapc_weights <- function(joinpoints, fit_start, fit_end, start, end) {
  segment_start <- c(fit_start, joinpoints)
  segment_end <- c(joinpoints, fit_end)
  overlap <- pmax(
    0,
    pmin(end, segment_end) - pmax(start, segment_start)
  )
  overlap / sum(overlap)
}

#' Internal helper: Validate APC/AAPC calculation arguments.
#'
#' @rdname crtrend-internals
#' @keywords internal
validate_calc_inputs <- function(x, level, ci_method, n_resamples, seed) {
  validate_rng_seed(seed)
  if (!is_bounded_integer(n_resamples, 2)) {
    stop(
      "`n_resamples` must be a finite integer in [2, 2147483646].",
      call. = FALSE
    )
  }
  if (!inherits(x, "crtrend_joinpoint")) {
    stop("`x` must be a `crtrend_joinpoint` object.", call. = FALSE)
  }
  if (
    !is.numeric(level) ||
      is.complex(level) ||
      length(level) != 1L ||
      is.na(level) ||
      level <= 0 ||
      level >= 1
  ) {
    stop("`level` must be one number between 0 and 1.", call. = FALSE)
  }
  if (!ci_method %in% c("EQ", "Param")) {
    stop("`ci_method` must be `EQ` or `Param`.", call. = FALSE)
  }
  if (
    !is.numeric(n_resamples) ||
      length(n_resamples) != 1L ||
      is.na(n_resamples) ||
      n_resamples < 2L ||
      n_resamples != as.integer(n_resamples)
  ) {
    stop("`n_resamples` must be one integer no smaller than 2.", call. = FALSE)
  }
  if (
    !is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
  ) {
    stop("`seed` must be NULL or one non-missing number.", call. = FALSE)
  }
}

#' Internal helper: Validate an AAPC interval.
#'
#' @rdname crtrend-internals
#' @keywords internal
validate_interval <- function(start, end) {
  if (
    !is_bounded_integer(start, -.Machine$integer.max, .Machine$integer.max) ||
      !is_bounded_integer(end, -.Machine$integer.max, .Machine$integer.max)
  ) {
    stop(
      "`start` and `end` must be finite integer years in [-2147483647, 2147483647].",
      call. = FALSE
    )
  }
  if (
    !is.numeric(start) ||
      length(start) != 1L ||
      is.na(start) ||
      !is.finite(start) ||
      !is.numeric(end) ||
      length(end) != 1L ||
      is.na(end) ||
      !is.finite(end) ||
      abs(start - round(start)) > sqrt(.Machine$double.eps) ||
      abs(end - round(end)) > sqrt(.Machine$double.eps) ||
      start >= end
  ) {
    stop(
      "`start` and `end` must be finite integer years with `start < end`.",
      call. = FALSE
    )
  }
}

#' Internal helper: Format an integer-year period value.
#'
#' @rdname crtrend-internals
#' @keywords internal
format_period_value <- function(x) {
  if (all(abs(x - round(x)) < sqrt(.Machine$double.eps))) {
    return(format(round(x), trim = TRUE, scientific = FALSE, nsmall = 0))
  }
  format(round(x, 2), trim = TRUE, scientific = FALSE, nsmall = 2)
}
