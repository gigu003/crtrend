#' Fit grouped joinpoint regression models
#'
#' @description
#' `joinpoint()` fits log-linear segmented regression models to a time-series
#' outcome. Data are split by `group_vars`, and the number of joinpoints is
#' selected by `selection_method` from `min_jp` through `max_jp`. APC and
#' AAPC are calculated separately with `calc_apc()` and `calc_aapc()`.
#'
#' @param data A data frame containing the time-series data.
#' @param time Name of the numeric or numeric-like time variable.
#' @param outcome Name of the positive numeric outcome variable. The model is
#'   fitted to `log(outcome)`.
#' @param group_vars Optional character vector of grouping-variable names. If
#'   `NULL`, one model is fitted to all observations.
#'   Output column names (for example `status`, `group_id`, or `segment`) and
#'   internal prefixes `trend_`, `selection_`, and `mc_` are reserved; rename
#'   conflicting grouping columns before fitting.
#' @param max_jp Maximum number of candidate joinpoints. The default is four.
#' @param min_jp Minimum number of joinpoints to consider. The default is zero.
#' @param min_obs Minimum number of observations allowed in each
#'   end segment and between adjacent Joinpoints, excluding the Joinpoint
#'   observation itself. The default is two, matching the NCI Grid Search
#'   setting.
#' @param seed Integer seed, default `7160`. Seeded fixed permutations use
#'   the native-compatible MT stream, reset for each group and continued
#'   across its model comparisons. Other random fitting uses [base::set.seed()].
#'   With `NULL`, permutations consume R's current random-number stream.
#'   Seeded fitting restores the caller's RNG state, including on error.
#'   Seeds must be finite integers in `[-2147483647, 2147483647]`.
#'   Count controls must fit in `[0, 2147483646]`, subject to their stated minima.
#' @param se Optional name of a column containing the standard error of
#'   `outcome` on its original scale. If supplied, weighted least squares is
#'   used with weights `(outcome / se)^2`, the delta-method weight for
#'   the log-linear model.
#' @param error Error structure. `"Ind"` assumes independent observations,
#'   `"AR1"` uses a lag-1 correlation structure, `"GLS"` uses a supplied
#'   general covariance matrix, and `"Survey"` uses a supplied design-based
#'   covariance matrix for aggregate complex-survey estimates.
#' @param rho Lag-1 correlation for `error = "AR1"`. It must be strictly
#'   between -1 and 1. The default is zero.
#' @param cov Optional numeric covariance matrix aligned with the rows of
#'   `data`. It is required for `error = "GLS"` or `"Survey"`. Cross-group
#'   covariance entries are ignored because groups are fitted independently.
#' @param cov_scale Scale of `cov`: `"outcome"` applies the delta-method log
#'   transformation, while `"log"` treats it as the covariance of
#'   `log(outcome)` directly.
#' @param selection_method Model selection method. The default, `"WBIC"`, uses
#'   the NCI-style data-dependent weighted BIC. `"Perm"` uses the Monte Carlo
#'   permutation test, `"BIC"` selects the classical BIC, `"BIC3"` uses the
#'   more conservative BIC3 penalty, `"MBIC"` uses NCI's modified BIC, and
#'   `"DDS"` applies NCI's data-dependent choice between BIC and BIC3. The
#'   former lower-case values remain accepted as compatibility aliases.
#' @param n_perm Number of Monte Carlo permutations used when
#'   `selection_method = "Perm"`. The default is `4499`.
#' @param alpha Significance level for the sequential permutation tests. The
#'   default is `0.05`.
#' @param perm_method Monte Carlo stopping method. The default, `"fixed"`,
#'   uses all `n_perm` resamples, matching NCI Grid Search. `"tsprt"` uses
#'   the truncated sequential probability ratio test of Fay, Kim, and Hachey
#'   (2007) as a package extension; NCI currently does not offer early
#'   stopping with Grid Search.
#' @param perm_alpha0 Type-I boundary error used by the tSPRT. The
#'   default is `0.0001`.
#' @param perm_beta0 Type-II boundary error used by the tSPRT. The
#'   default is `0.0001`.
#' @param perm_p1 Alternative p-value used by the tSPRT. If `NULL`, it
#'   is set to `0.8` times the individual sequential-test alpha. If supplied,
#'   it must be smaller than the smallest individual alpha.
#' @param perm_conf Confidence level for the Monte Carlo
#'   p-value interval. The default is `0.99`.
#' @param n_eq Number of empirical-quantile estimates cached in the fitted
#'   object for [calc_apc()] and [calc_aapc()]. For Joinpoint 6.1.0
#'   compatibility, this includes the original estimate and `n_eq - 1`
#'   residual-resampled refits. The default is `5001`.
#' @param eq_seed Integer seed used only for the cached
#'   empirical-quantile residual resamples. This is independent of `seed`,
#'   which controls random model-selection procedures. Use the same value to
#'   reproduce the cached empirical intervals. The default is `10000`,
#'   using a Joinpoint 6.1.0-compatible generator reset for each group.
#'   Use `NULL` to draw from R's current random-number stream without
#'   separately resetting the empirical-quantile seed.
#'   If `eq_seed = NULL` and `seed` is supplied, uncached EQ during fitting
#'   uses the locally seeded fitting stream; otherwise it consumes the caller's
#'   stream. A non-NULL `eq_seed` never changes the caller's RNG state.
#' @param wbic_weight Numerical method for the local partial R-squared used
#'   in Weight and WBIC. `"qr"` (default) preserves the stable QR calculation.
#'   `"joinpoint_6.1"` reproduces the explicit projection and scalar arithmetic
#'   order audited against Joinpoint 6.1.0. This opt-in compatibility method
#'   requires `error = "Ind"` and an input `se` column. It affects candidate
#'   Weight/WBIC scores, including those extracted by [calc_model()], not SSE,
#'   BIC, BIC3, or the fitting/inference algorithms. It can change the selected
#'   model when WBIC scores are close.
#'
#' @return An object of class `crtrend_joinpoint`, a list with these elements:
#'   \describe{
#'     \item{coefficients}{A tidy data frame containing the grouping variables,
#'       stable `group_id`, coefficient term, estimate, component type, and
#'       associated joinpoint.}
#'     \item{residuals, fitted.values}{Response-scale residuals and fitted
#'       values, aligned to the original input-row order. Unfitted groups have
#'       `NA` values.}
#'     \item{linear.predictors, link.residuals}{Fitted values and residuals on
#'       the log scale, aligned to the original input-row order.}
#'     \item{rank, df.residual}{Total rank and residual degrees of freedom
#'       summed over successfully fitted group models.}
#'     \item{formula, terms, model, na.action, method}{Standard model metadata.
#'       `model` preserves the requested grouping, time, outcome, and optional
#'       standard-error columns in their original row order.}
#'     \item{models}{A named list containing one independent
#'       `crtrend_joinpoint_model` object per smallest group. Each element uses
#'       standard model component names and additionally stores its group,
#'       numeric joinpoints, segments, selection results, inference, and
#'       empirical-quantile resamples. `selection$candidates` retains a
#'       lightweight summary of each successfully fitted candidate count;
#'       [calc_model()] extracts these summaries and full coefficient
#'       inference without repeating model selection.}
#'     \item{model_info}{A data frame with `group_id`, grouping variables,
#'       model status, selected Joinpoint count, time range, BIC, BIC3,
#'       Modified BIC, Weighted BIC, DDS diagnostics,
#'       the final permutation selection p-value, its Monte Carlo interval,
#'       the number of resamples actually used, and the largest feasible
#'       candidate count in `effective_max_jp`. For permutation selection,
#'       `mc_p_value` is an explicit alias for the final Monte Carlo p-value;
#'       `selection_p_value` is retained for compatibility.}
#'     \item{selection_tests}{A long-format audit table containing every
#'       sequential permutation comparison, its adjusted significance level,
#'       F-type statistic, Monte Carlo p-value and interval, number of
#'       permutations, and test decision. It has zero rows for non-permutation
#'       selection methods.}
#'     \item{joinpoints}{A numeric long-format table containing one row per
#'       selected joinpoint.}
#'     \item{segments}{A long-format table containing segment number, start,
#'       end, and fitted log slope for every successfully fitted group.}
#'     \item{control}{A named list containing the normalized fitting options.}
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' DDS uses the largest fitted candidate's SSE divided by `n - 2*k - 2`
#' to estimate variance, counting estimated joinpoint locations. MBIC uses
#' the whitened design matrix and the Joinpoint 6.1.0 score normalization:
#' the published expression minus the model-independent `2*log(n)/n` term.
#' This normalization does not change which candidate minimizes MBIC.
#' Seven-decimal agreement has been verified for the supplied independent,
#' SE-weighted 22-observation audit case, not for every error structure.
#'
#' For permutation selection, `selection_tests$test_status` is `"unavailable"`
#' if the observed or any simulated statistic is non-finite. The comparison
#' stops with missing p-value and decision; invalid draws are not discarded.
#' Selection stops at the current lower candidate, which is a fallback, not
#' a successful test accepting that model. `status = "ok"` describes fitting,
#' not availability of the selection test. Near-zero resampling SSEs use
#' explicit projection residuals to avoid subtractive cancellation.
#' The `tsprt` option is an extension for grid search, not a verified replica
#' of Joinpoint 6.1.0 early stopping. Seeded fixed permutations use forward
#' cumulative shuffling of whitened residuals and a native-compatible MT
#' stream. All four selection p-values match the supplied 22-observation,
#' independent SE-weighted 6.1.0 case at seven decimals. The three jp_data
#' cohorts also match native single-group, grouped and relabelled-group runs
#' (36 comparisons), supporting per-group seed resets for this setting.
#' Other error structures and tSPRT remain unbenchmarked against native output.
#' Unseeded fixed permutations and tSPRT retain R's `sample()` stream.
#'
#' The input must contain one observation for each time point within each
#' group. Missing values in the time, outcome, or grouping columns are not
#' silently removed. Groups with all-zero outcomes are returned with
#' `status = "zero_risk"`. Zero outcomes may have zero SE. These groups do
#' not prevent valid groups from being fitted. No pseudocount is added.
#' Groups containing a mixture of zero and positive
#' outcomes are returned with `status = "nonpositive_outcome"` because a
#' log-linear model is not defined for them.
#'
#' Joinpoints are selected from the observed integer-valued time points,
#' matching the annual-data use of the NCI Joinpoint software. When
#' `se` is supplied, it is interpreted on the original outcome scale
#' and transformed to a log-scale standard error using the first-order delta
#' method. The standard errors affect the weighted least-squares fit; `weighted_bic`
#' additionally uses a data-dependent penalty. When
#' `error = "AR1"`, correlations are `rho^abs(time_i - time_j)`; supplied
#' standard errors define the marginal variances, otherwise the AR(1)
#' correlation matrix is used up to an estimated residual scale. For `"GLS"`
#' and `"Survey"`, `cov` is treated as a known covariance matrix. Every grid
#' fit, information criterion, permutation refit, and empirical-quantile refit
#' is calculated after whitening by that covariance matrix. `"Survey"` does
#' not estimate a design covariance from raw strata, cluster, or weight
#' variables: users must first calculate the covariance matrix of the aggregate
#' time-specific estimates with appropriate survey software and supply it in
#' `cov`.
#'
#' When `cov_scale = "outcome"`, the covariance is transformed to the log scale
#' by the first-order delta method, `V_log[i,j] = V[i,j] / (y[i] * y[j])`.
#' Covariance entries between different `group_vars` combinations are ignored,
#' because each group is fitted as an independent model.
#'
#' When
#' `selection_method = "Perm"`, the selection p-value is the final
#' residual-permutation comparison of adjacent candidate models. It is not the
#' APC p-value returned by `calc_apc()`. The `mc_p_value_lower` and
#' `mc_p_value_upper` columns quantify Monte Carlo uncertainty at
#' `perm_conf`; they do not represent uncertainty in the
#' fitted regression coefficients. The individual sequential-test alpha is
#' adjusted as `alpha / (effective_max_jp - k)` for a null model with `k`
#' joinpoints, following the NCI procedure. `effective_max_jp` is the largest
#' candidate count compatible with the data and `min_obs` setting.
#'
#' For every candidate count, grid search first selects the joinpoint locations
#' that minimize SSE. BIC, BIC3, and WBIC are then used only to select among
#' these least-squares models. Every permutation data set is likewise refitted
#' over the full admissible grid for both models in the comparison; cached QR
#' projection bases avoid repeated calls to [stats::lm()].
#'
#' `wbic_weight = "joinpoint_6.1"` has been checked against the seven-decimal
#' Weight and WBIC output for `jp_data` cohort zero, with zero through four
#' joinpoints. It is sensitive to the time origin and is not a general
#' numerical improvement over QR. Singular or invalid projection calculations
#' stop with an error recommending `wbic_weight = "qr"`; no silent fallback
#' is performed. Compatibility for other datasets and error structures is not
#' established. This option does not resolve remaining SSE/BIC last-digit
#' differences. APC/AAPC results are unchanged if the selected model is unchanged.
#'
#' The returned object follows the usual R model-object conventions while
#' retaining grouped Joinpoint metadata. Methods are available for [coef()],
#' [fitted()], [residuals()], [predict()], [formula()], [nobs()],
#' [model.frame()], and [summary()]. Response-scale residuals are defined as
#' observed outcome minus `exp(fitted log outcome)`; use `type = "link"` for
#' quantities on the fitted log scale.
#'
#' `models` is the sole per-group model container. The previous computational
#' `fits` structure is not exposed. Top-level `model_info`, `joinpoints`,
#' `segments`, and `coefficients` are constructed once from these group models;
#' `group_id` provides a stable key between these tables.
#'
#' Point estimates come from the selected continuous grid-search model. For
#' standard errors and tests, the implementation follows the NCI convention
#' of fitting one block-diagonal unconstrained model for all segments with a
#' common residual variance. Observations located exactly at a selected
#' joinpoint are excluded from this calculation. Thus, for `k` joinpoints the
#' residual degrees of freedom are `n - 3 * k - 2`. When empirical-quantile
#' intervals are
#' requested by `calc_apc()` or `calc_aapc()`, residuals are resampled and the
#' selected number of joinpoints is refitted by an integer-year grid search
#' for every resample.
#'
#' @references
#' Fay, M. P., Kim, H.-J., and Hachey, M. (2007). On using truncated
#' sequential probability ratio test boundaries for Monte Carlo implementation
#' of hypothesis tests. *Journal of Computational and Graphical Statistics*,
#' 16(4), 946-967.
#'
#' Hudson, D. J. (1966). Fitting segmented curves whose join points have to be
#' estimated. *Journal of the American Statistical Association*, 61(316),
#' 1097-1129.
#'
#' Feder, P. I. (1975). On asymptotic distribution theory in segmented
#' regression problems: Identified case. *The Annals of Statistics*, 3(1),
#' 49-83.
#'
#' Kim, H.-J., Fay, M. P., Feuer, E. J., and Midthune, D. N. (2000).
#' Permutation tests for Joinpoint regression with applications to cancer
#' rates. *Statistics in Medicine*, 19(3), 335-351.
#'
#' Kim, H.-J., Chen, H.-S., Midthune, D., Wheeler, B., Buckman, D. W., Green,
#' D., Byrne, J., Luo, J., and Feuer, E. J. (2022). Data-driven choice of a
#' model selection method in Joinpoint regression. *Journal of Applied
#' Statistics*, 1-22.
#'
#' Zhang, N. R., and Siegmund, D. O. (2007). A modified Bayes information
#' criterion with applications to the analysis of comparative genomic
#' hybridization data. *Biometrics*, 63(1), 22-32.
#'
#' Lerman, P. M. (1980). Fitting segmented regression models by grid search.
#' *Applied Statistics*, 29(1), 77-84.
#'
#' Liu, B., Kim, H.-J., Feuer, E. J., and Graubard, B. I. (2022). Joinpoint
#' regression methods of aggregate outcomes for complex survey data.
#' *Journal of Survey Statistics and Methodology*, smac014.
#' \doi{10.1093/jssam/smac014}.
#'
#' @examples
#' years <- 1990:2023
#' example_data <- expand.grid(
#'   year = years,
#'   sex = c(1, 2),
#'   cancer = c("A", "B")
#' )
#' example_data$risk <- with(
#'   example_data,
#'   exp(-3 + 0.02 * (year - 1990) + 0.01 * (sex == 2))
#' )
#'
#' fit <- joinpoint(
#'   example_data,
#'   time = "year",
#'   outcome = "risk",
#'   group_vars = c("sex", "cancer"),
#'   max_jp = 1,
#'   seed = 1
#' )
#' calc_apc(fit)
#' calc_aapc(fit, start = 2000, end = 2020)
#'
#' @export
joinpoint <- function(
  data,
  time,
  outcome,
  group_vars = NULL,
  max_jp = 4L,
  min_jp = 0L,
  min_obs = 2L,
  seed = 7160L,
  se = NULL,
  error = c("Ind", "AR1", "GLS", "Survey"),
  rho = 0,
  cov = NULL,
  cov_scale = c("outcome", "log"),
  selection_method = c("WBIC", "Perm", "BIC", "BIC3", "MBIC", "DDS"),
  n_perm = 4499L,
  alpha = 0.05,
  perm_method = c("fixed", "tsprt"),
  perm_alpha0 = 0.0001,
  perm_beta0 = 0.0001,
  perm_p1 = NULL,
  perm_conf = 0.99,
  n_eq = 5001L,
  eq_seed = 10000L,
  wbic_weight = c("qr", "joinpoint_6.1")
) {
  options <- normalize_joinpoint_options(
    data = data,
    time = time,
    outcome = outcome,
    group_vars = group_vars,
    max_jp = max_jp,
    min_jp = min_jp,
    min_obs = min_obs,
    seed = seed,
    se = se,
    error = error,
    rho = rho,
    cov = cov,
    cov_scale = cov_scale,
    selection_method = selection_method,
    n_perm = n_perm,
    alpha = alpha,
    perm_method = perm_method,
    perm_alpha0 = perm_alpha0,
    perm_beta0 = perm_beta0,
    perm_p1 = perm_p1,
    perm_conf = perm_conf,
    n_eq = n_eq,
    eq_seed = eq_seed,
    wbic_weight = match.arg(wbic_weight)
  )
  analysis <- prepare_joinpoint_data(
    data = data,
    time = time,
    outcome = outcome,
    group_vars = options$group_vars,
    std_error = options$std_error
  )
  fits <- with_joinpoint_seed(
    seed,
    fit_joinpoint_groups(
      analysis,
      options,
      cache_empirical = is.null(options$eq_seed)
    )
  )
  if (!is.null(options$eq_seed)) {
    fits <- cache_joinpoint_empirical_resamples(
      fits,
      options$n_empirical_quantile_resamples,
      seed = options$eq_seed
    )
  }
  assemble_joinpoint_result(
    fits = fits,
    group_vars = options$group_vars,
    time = time,
    outcome = outcome,
    model_frame = as.data.frame(data, optional = TRUE)[
      c(options$group_vars, time, outcome, options$std_error)
    ],
    call = match.call(),
    options = options
  )
}

#' Internal helper: Normalize public arguments into the internal NCI analysis options.
#'
#' @rdname crtrend-internals
#' @keywords internal
normalize_joinpoint_options <- function(
  data,
  time,
  outcome,
  group_vars,
  max_jp,
  min_jp,
  min_obs,
  seed,
  se,
  error,
  rho,
  cov,
  cov_scale,
  selection_method,
  n_perm,
  alpha,
  perm_method,
  perm_alpha0,
  perm_beta0,
  perm_p1,
  perm_conf,
  n_eq,
  eq_seed,
  wbic_weight = "qr"
) {
  max_joinpoints <- max_jp
  min_joinpoints <- min_jp
  min_observations <- min_obs
  std_error <- se
  error_structure <- match.arg(error[1L], c("Ind", "AR1", "GLS", "Survey"))
  covariance_scale <- match.arg(cov_scale[1L], c("outcome", "log"))
  wbic_weight <- match.arg(wbic_weight, c("qr", "joinpoint_6.1"))
  if (
    wbic_weight == "joinpoint_6.1" &&
      (error_structure != "Ind" || is.null(se))
  ) {
    stop(
      '`wbic_weight = "joinpoint_6.1"` requires `error = "Ind"` and `se`.',
      call. = FALSE
    )
  }
  n_permutations <- n_perm
  permutation_method <- perm_method
  permutation_alpha0 <- perm_alpha0
  permutation_beta0 <- perm_beta0
  permutation_p1 <- perm_p1
  permutation_confidence_level <- perm_conf
  n_empirical_quantile_resamples <- n_eq
  validate_joinpoint_inputs(
    data = data,
    time = time,
    outcome = outcome,
    group_vars = group_vars,
    max_joinpoints = max_joinpoints,
    min_joinpoints = min_joinpoints,
    min_observations = min_observations,
    seed = seed,
    std_error = std_error,
    error_structure = error_structure,
    rho = rho,
    covariance = cov,
    covariance_scale = covariance_scale,
    selection_method = selection_method,
    n_permutations = n_permutations,
    alpha = alpha,
    eq_seed = eq_seed,
    permutation_method = permutation_method,
    permutation_alpha0 = permutation_alpha0,
    permutation_beta0 = permutation_beta0,
    permutation_p1 = permutation_p1,
    permutation_confidence_level = permutation_confidence_level,
    n_empirical_quantile_resamples = n_empirical_quantile_resamples
  )
  selection_method <- match.arg(
    selection_method[1L],
    c(
      "WBIC",
      "Perm",
      "BIC",
      "BIC3",
      "MBIC",
      "DDS",
      "weighted_bic",
      "permutation",
      "bic",
      "bic3",
      "modified_bic",
      "dds"
    )
  )
  permutation_method <- match.arg(
    permutation_method[1L],
    c("fixed", "tsprt")
  )
  selection_method_label <- normalize_selection_method_label(selection_method)
  selection_method <- normalize_selection_method(selection_method)
  result <- list(
    max_joinpoints = as.integer(max_joinpoints),
    min_joinpoints = as.integer(min_joinpoints),
    min_observations = as.integer(min_observations),
    std_error = std_error,
    weighted = !is.null(std_error) || error_structure != "Ind",
    error_structure = error_structure,
    rho = rho,
    covariance = cov,
    covariance_scale = covariance_scale,
    group_vars = group_vars %||% character(),
    selection_method = selection_method,
    selection_method_label = selection_method_label,
    wbic_weight = wbic_weight,
    n_permutations = as.integer(n_permutations),
    alpha = alpha,
    permutation_method = permutation_method,
    permutation_alpha0 = permutation_alpha0,
    permutation_beta0 = permutation_beta0,
    permutation_p1 = permutation_p1,
    permutation_confidence_level = permutation_confidence_level,
    permutation_seed = seed,
    eq_seed = eq_seed,
    n_empirical_quantile_resamples = as.integer(
      n_empirical_quantile_resamples
    ),
    permutation_boundaries = make_permutation_boundaries(
      selection_method = selection_method,
      max_joinpoints = max_joinpoints,
      n_permutations = n_permutations,
      alpha = alpha,
      permutation_method = permutation_method,
      permutation_alpha0 = permutation_alpha0,
      permutation_beta0 = permutation_beta0,
      permutation_p1 = permutation_p1
    )
  )
}

#' Internal helper: Return the left-hand value unless it is `NULL`.
#'
#' @rdname crtrend-internals
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Internal helper: Construct cached tSPRT boundaries for the sequential permutation tests.
#'
#' @rdname crtrend-internals
#' @keywords internal
make_permutation_boundaries <- function(
  selection_method,
  max_joinpoints,
  n_permutations,
  alpha,
  permutation_method,
  permutation_alpha0,
  permutation_beta0,
  permutation_p1
) {
  if (
    selection_method != "permutation" ||
      permutation_method != "tsprt" ||
      max_joinpoints < 1L
  ) {
    return(NULL)
  }
  null_joinpoints <- seq_len(max_joinpoints) - 1L
  stats::setNames(
    lapply(
      null_joinpoints,
      function(null_joinpoints) {
        make_tsprt_boundary(
          n_permutations = n_permutations,
          alpha = alpha / (max_joinpoints - null_joinpoints),
          alpha0 = permutation_alpha0,
          beta0 = permutation_beta0,
          p1 = permutation_p1
        )
      }
    ),
    as.character(null_joinpoints)
  )
}

#' Internal helper: Validate and construct the internal grouped trend data.
#'
#' @rdname crtrend-internals
#' @keywords internal
prepare_joinpoint_data <- function(data, time, outcome, group_vars, std_error) {
  analysis <- data.frame(
    data[group_vars],
    trend_time = coerce_numeric_column(data[[time]], time),
    trend_outcome = coerce_numeric_column(data[[outcome]], outcome),
    check.names = FALSE
  )
  analysis$trend_row_id <- seq_len(nrow(analysis))
  if (!is.null(std_error)) {
    analysis$trend_se <- coerce_numeric_column(data[[std_error]], std_error)
    if (anyNA(analysis$trend_se)) {
      stop(
        "Missing values are not allowed in time, outcome, standard error, or grouping columns.",
        call. = FALSE
      )
    }
    if (
      any(
        analysis$trend_se < 0 |
          (analysis$trend_se == 0 & analysis$trend_outcome != 0)
      )
    ) {
      stop("`se` must contain positive values.", call. = FALSE)
    }
    analysis$trend_weight <-
      (analysis$trend_outcome / analysis$trend_se)^2
    analysis$trend_weight[analysis$trend_se == 0] <- NA_real_
  }
  required <- c(group_vars, "trend_time", "trend_outcome")
  if (!is.null(std_error)) {
    required <- c(required, "trend_se")
  }
  if (anyNA(analysis[required])) {
    stop(
      "Missing values are not allowed in time, outcome, standard error, or grouping columns.",
      call. = FALSE
    )
  }
  if (
    any(
      abs(analysis$trend_time - round(analysis$trend_time)) >
        sqrt(.Machine$double.eps)
    )
  ) {
    stop("`time` must contain integer-valued years.", call. = FALSE)
  }
  duplicate_key <- joinpoint_duplicate_key(analysis, group_vars)
  if (anyDuplicated(duplicate_key)) {
    stop(
      "Each group must have at most one observation per time point.",
      call. = FALSE
    )
  }
  analysis
}

#' Internal helper: Fit the joinpoint model independently within each group.
#'
#' @rdname crtrend-internals
#' @keywords internal
fit_joinpoint_groups <- function(data, options, cache_empirical = TRUE) {
  attr(data, "crtrend_wbic_weight") <- options$wbic_weight %||% "qr"
  groups <- split_group_rows(data, options$group_vars)
  grid_cache <- new.env(parent = emptyenv())
  lapply(
    groups,
    fit_joinpoint_group,
    group_vars = options$group_vars,
    max_joinpoints = options$max_joinpoints,
    min_joinpoints = options$min_joinpoints,
    min_observations = options$min_observations,
    weighted = options$weighted,
    error_structure = options$error_structure,
    rho = options$rho,
    covariance = options$covariance,
    covariance_scale = options$covariance_scale,
    selection_method = options$selection_method,
    n_permutations = options$n_permutations,
    alpha = options$alpha,
    permutation_method = options$permutation_method,
    permutation_alpha0 = options$permutation_alpha0,
    permutation_beta0 = options$permutation_beta0,
    permutation_p1 = options$permutation_p1,
    permutation_boundaries = options$permutation_boundaries,
    permutation_confidence_level = options$permutation_confidence_level,
    permutation_seed = options$permutation_seed,
    n_empirical_quantile_resamples = options$n_empirical_quantile_resamples,
    cache_empirical = cache_empirical,
    grid_cache = grid_cache
  )
}

#' Internal helper: Cache empirical-quantile resamples after model fitting.
#'
#' @rdname crtrend-internals
#' @keywords internal
cache_joinpoint_empirical_resamples <- function(
  fits,
  n_resamples,
  seed = NULL
) {
  lapply(fits, function(fit) {
    if (fit$status != "ok") {
      return(fit)
    }
    fit$eq_draws <- tryCatch(
      nci_empirical_model_draws(
        fit,
        n_resamples,
        seed = seed
      ),
      error = function(error) NULL
    )
    fit$eq_n_resamples <- n_resamples
    fit
  })
}

#' Internal helper: Assemble the public `crtrend_joinpoint` result.
#'
#' @rdname crtrend-internals
#' @keywords internal
assemble_joinpoint_result <- function(
  fits,
  group_vars,
  time,
  outcome,
  model_frame,
  call,
  options
) {
  model_info <- do.call(rbind, lapply(fits, fit_info_row))
  rownames(model_info) <- NULL
  model_info$joinpoints <- NULL
  model_info <- data.frame(
    group_id = seq_len(nrow(model_info)),
    model_info,
    check.names = FALSE
  )
  formula <- stats::reformulate(time, response = paste0("log(", outcome, ")"))
  environment(formula) <- baseenv()
  components <- joinpoint_model_components(
    fits = fits,
    n = nrow(model_frame),
    group_vars = group_vars,
    time = time
  )
  observation_names <- row.names(model_frame)
  names(components$response_residuals) <- observation_names
  names(components$response_fitted) <- observation_names
  names(components$link_residuals) <- observation_names
  names(components$link_fitted) <- observation_names
  models <- build_joinpoint_models(
    fits = fits,
    model_frame = model_frame,
    call = call,
    formula = formula,
    time = time,
    outcome = outcome,
    group_vars = group_vars
  )
  names(models) <- names(fits)
  joinpoints <- combine_joinpoint_table(models, group_vars)
  segments <- combine_segment_table(models, group_vars)
  selection_tests <- combine_selection_test_table(models, group_vars)
  structure(
    list(
      coefficients = components$coefficients,
      residuals = components$response_residuals,
      fitted.values = components$response_fitted,
      linear.predictors = components$link_fitted,
      link.residuals = components$link_residuals,
      rank = components$rank,
      df.residual = components$df_residual,
      formula = formula,
      terms = stats::terms(formula),
      model = model_frame,
      na.action = NULL,
      method = "NCI-style grid-search joinpoint regression",
      models = models,
      model_info = model_info,
      selection_tests = selection_tests,
      joinpoints = joinpoints,
      segments = segments,
      control = list(
        group_vars = group_vars,
        time = time,
        outcome = outcome,
        se = options$std_error,
        error = options$error_structure,
        rho = options$rho,
        cov_scale = options$covariance_scale,
        min_jp = options$min_joinpoints,
        max_jp = options$max_joinpoints,
        min_obs = options$min_observations,
        selection_method = options$selection_method_label,
        wbic_weight = options$wbic_weight,
        n_perm = options$n_permutations,
        alpha = options$alpha,
        perm_method = options$permutation_method,
        perm_alpha0 = options$permutation_alpha0,
        perm_beta0 = options$permutation_beta0,
        perm_p1 = options$permutation_p1,
        perm_conf = options$permutation_confidence_level,
        n_eq = options$n_empirical_quantile_resamples,
        eq_seed = options$eq_seed
      ),
      call = call
    ),
    class = "crtrend_joinpoint"
  )
}

#' Internal helper: Validate the public `joinpoint()` arguments.
#'
#' @rdname crtrend-internals
#' @keywords internal
validate_joinpoint_inputs <- function(
  data,
  time,
  outcome,
  group_vars,
  max_joinpoints,
  min_joinpoints,
  min_observations,
  seed,
  std_error,
  error_structure,
  rho,
  covariance,
  covariance_scale,
  selection_method,
  n_permutations,
  alpha,
  eq_seed,
  permutation_method,
  permutation_alpha0,
  permutation_beta0,
  permutation_p1,
  permutation_confidence_level,
  n_empirical_quantile_resamples
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  check_single_name(time, "time")
  check_single_name(outcome, "outcome")
  integer_options <- list(
    max_jp = max_joinpoints,
    min_jp = min_joinpoints,
    min_obs = min_observations,
    n_perm = n_permutations,
    n_eq = n_empirical_quantile_resamples
  )
  lower <- c(max_jp = 0, min_jp = 0, min_obs = 2, n_perm = 1, n_eq = 2)
  for (name in names(integer_options)) {
    if (!is_bounded_integer(integer_options[[name]], lower[[name]])) {
      stop(
        "`",
        name,
        "` must be a finite integer in [",
        lower[[name]],
        ", 2147483646].",
        call. = FALSE
      )
    }
  }
  validate_rng_seed(seed)
  validate_rng_seed(eq_seed, "eq_seed")
  real_options <- list(
    rho = rho,
    alpha = alpha,
    perm_alpha0 = permutation_alpha0,
    perm_beta0 = permutation_beta0,
    perm_conf = permutation_confidence_level
  )
  if (!is.null(permutation_p1)) {
    real_options$perm_p1 <- permutation_p1
  }
  for (name in names(real_options)) {
    value <- real_options[[name]]
    if (
      !is.numeric(value) ||
        is.complex(value) ||
        length(value) != 1L ||
        !is.finite(value)
    ) {
      stop("`", name, "` must be one finite real number.", call. = FALSE)
    }
  }
  if (!is.null(std_error)) {
    check_single_name(std_error, "se")
  }
  if (
    !is.numeric(rho) ||
      length(rho) != 1L ||
      is.na(rho) ||
      !is.finite(rho) ||
      rho <= -1 ||
      rho >= 1
  ) {
    stop("`rho` must be one number strictly between -1 and 1.", call. = FALSE)
  }
  if (error_structure %in% c("GLS", "Survey") && is.null(covariance)) {
    stop(
      "`cov` is required for `error = \"GLS\"` or `\"Survey\"`.",
      call. = FALSE
    )
  }
  if (!is.null(covariance)) {
    if (
      !is.matrix(covariance) ||
        !is.numeric(covariance) ||
        any(dim(covariance) != nrow(data)) ||
        any(!is.finite(covariance))
    ) {
      stop(
        "`cov` must be a finite numeric square matrix with one row per input row.",
        call. = FALSE
      )
    }
    if (!isTRUE(all.equal(covariance, t(covariance), tolerance = 1e-10))) {
      stop("`cov` must be symmetric.", call. = FALSE)
    }
    if (!is.null(std_error)) {
      stop("Supply either `se` or `cov`, not both.", call. = FALSE)
    }
  }
  if (error_structure == "Ind" && !is.null(covariance)) {
    stop(
      "Use `error = \"GLS\"` or `\"Survey\"` when supplying `cov`.",
      call. = FALSE
    )
  }
  if (error_structure == "AR1" && !is.null(covariance)) {
    stop("`cov` cannot be combined with `error = \"AR1\"`.", call. = FALSE)
  }
  if (!is.character(selection_method) || length(selection_method) < 1L) {
    stop(
      "`selection_method` must contain at least one character value.",
      call. = FALSE
    )
  }
  if (!is.character(permutation_method) || length(permutation_method) < 1L) {
    stop(
      "`perm_method` must contain at least one character value.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(n_permutations) ||
      length(n_permutations) != 1L ||
      is.na(n_permutations) ||
      n_permutations < 1 ||
      n_permutations != as.integer(n_permutations)
  ) {
    stop("`n_perm` must be one positive integer.", call. = FALSE)
  }
  if (
    !is.numeric(n_empirical_quantile_resamples) ||
      length(n_empirical_quantile_resamples) != 1L ||
      is.na(n_empirical_quantile_resamples) ||
      n_empirical_quantile_resamples < 2L ||
      n_empirical_quantile_resamples !=
        as.integer(n_empirical_quantile_resamples)
  ) {
    stop("`n_eq` must be one integer no smaller than 2.", call. = FALSE)
  }
  if (
    !is.numeric(min_observations) ||
      length(min_observations) != 1L ||
      is.na(min_observations) ||
      min_observations < 2 ||
      min_observations != as.integer(min_observations)
  ) {
    stop(
      "`min_obs` must be one integer no smaller than 2.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop("`alpha` must be one number between 0 and 1.", call. = FALSE)
  }
  for (argument in c("permutation_alpha0", "permutation_beta0")) {
    value <- get(argument)
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        is.na(value) ||
        value <= 0 ||
        value >= 1
    ) {
      stop(
        "`",
        if (argument == "permutation_alpha0") "perm_alpha0" else "perm_beta0",
        "` must be one number between 0 and 1.",
        call. = FALSE
      )
    }
  }
  if (
    !is.null(permutation_p1) &&
      (!is.numeric(permutation_p1) ||
        length(permutation_p1) != 1L ||
        is.na(permutation_p1) ||
        permutation_p1 <= 0 ||
        permutation_p1 >= alpha / max(1, max_joinpoints))
  ) {
    stop(
      "`perm_p1` must be NULL or smaller than `alpha / max_jp`.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(permutation_confidence_level) ||
      length(permutation_confidence_level) != 1L ||
      is.na(permutation_confidence_level) ||
      permutation_confidence_level <= 0 ||
      permutation_confidence_level >= 1
  ) {
    stop(
      "`perm_conf` must be one number between 0 and 1.",
      call. = FALSE
    )
  }
  if (is.null(group_vars)) {
    group_vars <- character()
  }
  if (!is.character(group_vars) || anyNA(group_vars)) {
    stop("`group_vars` must be a character vector or NULL.", call. = FALSE)
  }
  validate_joinpoint_group_names(group_vars)
  if (anyDuplicated(c(time, outcome, group_vars, std_error))) {
    stop("`time`, `outcome`, and `group_vars` must be unique.", call. = FALSE)
  }
  missing_columns <- setdiff(
    c(time, outcome, group_vars, std_error),
    names(data)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "Missing columns: ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (
    !is.numeric(max_joinpoints) ||
      length(max_joinpoints) != 1L ||
      is.na(max_joinpoints) ||
      max_joinpoints < 0 ||
      max_joinpoints != as.integer(max_joinpoints)
  ) {
    stop("`max_jp` must be one non-negative integer.", call. = FALSE)
  }
  if (
    !is.numeric(min_joinpoints) ||
      length(min_joinpoints) != 1L ||
      is.na(min_joinpoints) ||
      min_joinpoints < 0 ||
      min_joinpoints != as.integer(min_joinpoints) ||
      min_joinpoints > max_joinpoints
  ) {
    stop(
      "`min_jp` must be a non-negative integer no greater than `max_jp`.",
      call. = FALSE
    )
  }
  if (
    !is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
  ) {
    stop("`seed` must be NULL or one non-missing number.", call. = FALSE)
  }
  if (
    !is.null(eq_seed) &&
      (!is.numeric(eq_seed) || length(eq_seed) != 1L || is.na(eq_seed))
  ) {
    stop("`eq_seed` must be NULL or one non-missing number.", call. = FALSE)
  }
}

#' Internal helper: Validate a single non-empty column name.
#'
#' @rdname crtrend-internals
#' @keywords internal
check_single_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", arg, "` must be one non-empty column name.", call. = FALSE)
  }
}

#' Internal helper: Coerce a time-series column to finite numeric values.
#'
#' @rdname crtrend-internals
#' @keywords internal
coerce_numeric_column <- function(x, name) {
  if (is.numeric(x)) {
    value <- as.numeric(x)
  } else {
    value <- suppressWarnings(as.numeric(as.character(x)))
  }
  if (anyNA(value) && !anyNA(x)) {
    stop("`", name, "` must be numeric or numeric-like.", call. = FALSE)
  }
  if (any(!is.finite(value))) {
    stop("`", name, "` cannot contain infinite values.", call. = FALSE)
  }
  value
}

#' Internal helper: Construct keys used to detect duplicate group-time observations.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_duplicate_key <- function(data, group_vars) {
  key_columns <- c(group_vars, "trend_time")
  if (length(key_columns) == 1L) {
    return(as.character(data[[key_columns]]))
  }
  do.call(paste, c(unname(data[key_columns]), sep = "\r"))
}

#' Internal helper: Split the analysis data into groups.
#'
#' @rdname crtrend-internals
#' @keywords internal
split_group_rows <- function(data, group_vars) {
  if (length(group_vars) == 0L) {
    return(stats::setNames(list(data), "all"))
  }
  key <- do.call(
    interaction,
    c(
      unname(data[group_vars]),
      list(drop = TRUE, lex.order = TRUE, sep = "\r")
    )
  )
  split(data, key, drop = TRUE)
}

#' Internal helper: Attach a log-scale covariance model to one fitted group.
#'
#' @rdname crtrend-internals
#' @keywords internal
prepare_group_error_model <- function(
  data,
  error_structure,
  rho,
  covariance,
  covariance_scale
) {
  n <- nrow(data)
  if (error_structure %in% c("GLS", "Survey")) {
    row_id <- data$trend_row_id
    covariance <- covariance[row_id, row_id, drop = FALSE]
    if (covariance_scale == "outcome") {
      inverse_outcome <- 1 / data$trend_outcome
      covariance <- covariance * tcrossprod(inverse_outcome)
    }
  } else if (error_structure == "AR1") {
    lag <- abs(outer(data$trend_time, data$trend_time, `-`))
    covariance <- rho^lag
    if ("trend_se" %in% names(data)) {
      log_se <- data$trend_se / data$trend_outcome
      covariance <- covariance * tcrossprod(log_se)
    }
  } else if ("trend_se" %in% names(data)) {
    log_se <- data$trend_se / data$trend_outcome
    covariance <- diag(log_se^2, nrow = n)
  } else {
    covariance <- diag(n)
  }
  covariance <- (covariance + t(covariance)) / 2
  factorization <- tryCatch(chol(covariance), error = function(error) NULL)
  if (is.null(factorization)) {
    stop(
      "The within-group log-scale covariance matrix must be positive definite.",
      call. = FALSE
    )
  }
  ids <- as.character(data$trend_row_id)
  dimnames(covariance) <- list(ids, ids)
  attr(data, "crtrend_covariance") <- covariance
  attr(data, "crtrend_error") <- error_structure
  attr(data, "crtrend_known_covariance") <- error_structure %in%
    c("GLS", "Survey")
  data
}

#' Internal helper: Extract the covariance matrix for the current data rows.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_covariance <- function(data) {
  covariance <- attr(data, "crtrend_covariance", exact = TRUE)
  if (is.null(covariance)) {
    if ("trend_weight" %in% names(data)) {
      return(diag(1 / data$trend_weight, nrow = nrow(data)))
    }
    return(diag(nrow(data)))
  }
  ids <- as.character(data$trend_row_id)
  covariance[ids, ids, drop = FALSE]
}

#' Internal helper: Construct a GLS whitening matrix.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_whitener <- function(data) {
  factorization <- chol(joinpoint_covariance(data))
  forwardsolve(t(factorization), diag(nrow(data)))
}

#' Internal helper: Transform values to independent unit-variance coordinates.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_whiten <- function(data, value) {
  joinpoint_whitener(data) %*% value
}

#' Internal helper: Transform whitened values back to observation coordinates.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_unwhiten <- function(data, value) {
  solve(joinpoint_whitener(data), value)
}

#' Internal helper: Fit and select a model for one group.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Hudson, D. J. (1966). Fitting segmented curves whose join points have to be
#' estimated. *Journal of the American Statistical Association*, 61(316),
#' 1097-1129.
#'
#' Feder, P. I. (1975). On asymptotic distribution theory in segmented
#' regression problems: Identified case. *The Annals of Statistics*, 3(1),
#' 49-83.
#'
#' Lerman, P. M. (1980). Fitting segmented regression models by grid search.
#' *Applied Statistics*, 29(1), 77-84.
fit_joinpoint_group <- function(
  group,
  group_vars,
  max_joinpoints,
  min_joinpoints,
  min_observations,
  weighted,
  error_structure,
  rho,
  covariance,
  covariance_scale,
  selection_method,
  n_permutations,
  alpha,
  permutation_method,
  permutation_alpha0,
  permutation_beta0,
  permutation_p1,
  permutation_boundaries,
  permutation_confidence_level,
  n_empirical_quantile_resamples,
  cache_empirical,
  grid_cache,
  permutation_seed = NULL
) {
  group <- group[order(group$trend_time), , drop = FALSE]
  group_key <- group[1L, group_vars, drop = FALSE]
  year_min <- min(group$trend_time)
  year_max <- max(group$trend_time)

  if (all(group$trend_outcome == 0)) {
    return(empty_joinpoint_group(
      group_key,
      group_vars,
      status = "zero_risk",
      start = year_min,
      end = year_max,
      data = group
    ))
  }
  if (any(group$trend_outcome <= 0)) {
    return(empty_joinpoint_group(
      group_key,
      group_vars,
      status = "nonpositive_outcome",
      start = year_min,
      end = year_max,
      data = group
    ))
  }
  if (nrow(group) < 4L) {
    return(empty_joinpoint_group(
      group_key,
      group_vars,
      status = "too_few_observations",
      start = year_min,
      end = year_max,
      data = group
    ))
  }

  group <- prepare_group_error_model(
    group,
    error_structure = error_structure,
    rho = rho,
    covariance = covariance,
    covariance_scale = covariance_scale
  )
  group$trend_log <- log(group$trend_outcome)
  candidates <- lapply(
    min_joinpoints:max_joinpoints,
    fit_candidate,
    data = group,
    weighted = weighted,
    min_observations = min_observations,
    selection_method = selection_method,
    grid_cache = grid_cache
  )
  candidates <- Filter(Negate(is.null), candidates)
  if (length(candidates) == 0L) {
    return(empty_joinpoint_group(
      group_key,
      group_vars,
      status = "model_failed",
      start = year_min,
      end = year_max,
      data = group
    ))
  }

  selection <- select_joinpoint_candidate(
    candidates,
    data = group,
    weighted = weighted,
    selection_method = selection_method,
    min_joinpoints = min_joinpoints,
    n_permutations = n_permutations,
    alpha = alpha,
    permutation_method = permutation_method,
    permutation_alpha0 = permutation_alpha0,
    permutation_beta0 = permutation_beta0,
    permutation_p1 = permutation_p1,
    max_joinpoints = max_joinpoints,
    permutation_boundaries = permutation_boundaries,
    permutation_confidence_level = permutation_confidence_level,
    permutation_seed = permutation_seed
  )
  best <- selection$best
  joinpoints <- best$joinpoints
  boundaries <- c(year_min, joinpoints, year_max)
  if (any(!is.finite(boundaries)) || any(diff(boundaries) <= 0)) {
    return(empty_joinpoint_group(
      group_key,
      group_vars,
      status = "invalid_joinpoints",
      start = year_min,
      end = year_max,
      data = group
    ))
  }

  result <- list(
    group_key = group_key,
    status = "ok",
    n_joinpoints = length(joinpoints),
    joinpoints = joinpoints,
    start = year_min,
    end = year_max,
    bic = best$bic,
    bic3 = best$bic3,
    modified_bic = best$modified_bic,
    weighted_bic = best$weighted_bic,
    partial_r_squared = best$partial_r_squared,
    dds_delta = best$dds_delta,
    dds_method = selection$dds_method %||% NA_character_,
    dds_sigma_squared = selection$dds_sigma_squared %||% NA_real_,
    dds_variance_df = selection$dds_variance_df %||% NA_integer_,
    dds_bic_delta = selection$dds_bic_delta %||% NA_real_,
    dds_bic3_delta = selection$dds_bic3_delta %||% NA_real_,
    dds_delta_max = selection$dds_delta_max %||% NA_real_,
    dds_delta_min = selection$dds_delta_min %||% NA_real_,
    dds_decision_rule = selection$dds_decision_rule %||% NA_character_,
    selection_score = best[[selection_score_name(selection_method)]],
    selection_p_value = selection$p_value,
    selection_p_value_lower = selection$p_value_lower,
    selection_p_value_upper = selection$p_value_upper,
    selection_p_value_confidence_level = selection$confidence_level,
    selection_n_permutations = selection$n_permutations,
    selection_stopped_early = selection$stopped_early,
    selection_tests = selection$tests,
    candidates = joinpoint_candidate_summary(
      selection$candidates %||% candidates,
      nrow(group),
      best$n_joinpoints
    ),
    effective_max_joinpoints = selection$effective_max_joinpoints,
    mc_p_value = selection$p_value,
    mc_p_value_lower = selection$p_value_lower,
    mc_p_value_upper = selection$p_value_upper,
    mc_p_value_confidence_level = selection$confidence_level,
    mc_n_permutations = selection$n_permutations,
    mc_stopped_early = selection$stopped_early,
    min_observations = min_observations,
    weighted = weighted,
    error_structure = error_structure,
    rho = rho,
    covariance_scale = covariance_scale,
    data = group,
    model = best$model
  )
  result$eq_draws <- if (cache_empirical) {
    tryCatch(
      nci_empirical_model_draws(result, n_empirical_quantile_resamples),
      error = function(error) NULL
    )
  } else {
    NULL
  }
  result$inference <- nci_joinpoint_inference(result)
  result$eq_n_resamples <- n_empirical_quantile_resamples
  result
}

#' Internal helper: Fit the best candidate model for a fixed number of joinpoints.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Lerman, P. M. (1980). Fitting segmented regression models by grid search.
#' *Applied Statistics*, 29(1), 77-84.
#'
#' Hudson, D. J. (1966). Fitting segmented curves whose join points have to be
#' estimated. *Journal of the American Statistical Association*, 61(316),
#' 1097-1129.
fit_candidate <- function(
  n_joinpoints,
  data,
  weighted,
  min_observations,
  selection_method,
  grid_cache = NULL
) {
  if (n_joinpoints == 0L) {
    candidate_joinpoints <- list(numeric())
  } else {
    candidate_joinpoints <- cached_integer_joinpoint_grid(
      data$trend_time,
      n_joinpoints,
      min_points = min_observations,
      cache = grid_cache
    )
    if (is.null(candidate_joinpoints)) {
      return(NULL)
    }
  }
  if (n_joinpoints > 0L) {
    scores <- lapply(
      candidate_joinpoints,
      fixed_joinpoint_scores,
      data = data,
      weighted = weighted,
      calculate_partial = FALSE
    )
    criterion <- vapply(
      scores,
      `[[`,
      numeric(1),
      "sse"
    )
    best_index <- if (all(!is.finite(criterion))) {
      1L
    } else {
      which.min(criterion)
    }
    model <- fit_fixed_joinpoint_model(
      candidate_joinpoints[[best_index]],
      data,
      weighted = weighted
    )
    return(make_candidate(
      model,
      n_joinpoints,
      joinpoints = candidate_joinpoints[[best_index]],
      data = data,
      weighted = weighted,
      grid_joinpoints = if (selection_method == "permutation") {
        candidate_joinpoints[is.finite(criterion)]
      } else {
        NULL
      }
    ))
  }

  model <- fit_lm_model(data, weighted = weighted)
  make_candidate(
    model,
    n_joinpoints,
    joinpoints = if (n_joinpoints == 0L) numeric() else NULL,
    data = data,
    weighted = weighted,
    grid_joinpoints = if (selection_method == "permutation") {
      list(numeric())
    } else {
      NULL
    }
  )
}

#' Internal helper: Map a public selection label to an internal score column.
#'
#' @rdname crtrend-internals
#' @keywords internal
selection_score_name <- function(selection_method) {
  switch(
    selection_method,
    permutation = "sse",
    bic = "bic",
    bic3 = "bic3",
    weighted_bic = "weighted_bic",
    modified_bic = "modified_bic",
    dds = "dds_score"
  )
}

#' Internal helper: Normalize the model-selection method label.
#'
#' @rdname crtrend-internals
#' @keywords internal
normalize_selection_method <- function(selection_method) {
  switch(
    selection_method,
    WBIC = "weighted_bic",
    Perm = "permutation",
    BIC = "bic",
    BIC3 = "bic3",
    MBIC = "modified_bic",
    DDS = "dds",
    weighted_bic = "weighted_bic",
    permutation = "permutation",
    bic = "bic",
    bic3 = "bic3",
    modified_bic = "modified_bic",
    dds = "dds"
  )
}

#' Internal helper: Return the short public model-selection label.
#'
#' @rdname crtrend-internals
#' @keywords internal
normalize_selection_method_label <- function(selection_method) {
  switch(
    normalize_selection_method(selection_method),
    weighted_bic = "WBIC",
    permutation = "Perm",
    bic = "BIC",
    bic3 = "BIC3",
    modified_bic = "MBIC",
    dds = "DDS"
  )
}

#' Internal helper: Construct an empty permutation-selection audit table.
#'
#' @rdname crtrend-internals
#' @keywords internal
empty_selection_tests <- function() {
  data.frame(
    test_status = character(),
    step = integer(),
    k_null = integer(),
    k_alt = integer(),
    alpha = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    p_value_lower = numeric(),
    p_value_upper = numeric(),
    confidence_level = numeric(),
    n_permutations = integer(),
    stopped_early = logical(),
    reject = logical()
  )
}

#' Internal helper: Record one sequential permutation comparison.
#'
#' @rdname crtrend-internals
#' @keywords internal
selection_test_row <- function(step, k_null, k_alt, alpha, result) {
  data.frame(
    test_status = if (is.finite(result$p_value)) "ok" else "unavailable",
    step = as.integer(step),
    k_null = as.integer(k_null),
    k_alt = as.integer(k_alt),
    alpha = alpha,
    statistic = result$statistic %||% NA_real_,
    p_value = result$p_value,
    p_value_lower = result$p_value_lower,
    p_value_upper = result$p_value_upper,
    confidence_level = result$confidence_level,
    n_permutations = result$n_permutations,
    stopped_early = result$stopped_early,
    reject = result$reject
  )
}

#' Internal helper: Select among candidate joinpoint counts.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Fay, M. P., Feuer, E. J., and Midthune, D. N. (2000).
#' Permutation tests for Joinpoint regression with applications to cancer
#' rates. *Statistics in Medicine*, 19(3), 335-351.
#'
#' Kim, H.-J., Chen, H.-S., Midthune, D., Wheeler, B., Buckman, D. W., Green,
#' D., Byrne, J., Luo, J., and Feuer, E. J. (2022). Data-driven choice of a
#' model selection method in Joinpoint regression. *Journal of Applied
#' Statistics*, 1-22.
select_joinpoint_candidate <- function(
  candidates,
  data,
  weighted,
  selection_method,
  min_joinpoints,
  n_permutations,
  alpha,
  permutation_method,
  permutation_alpha0,
  permutation_beta0,
  permutation_p1,
  max_joinpoints,
  permutation_boundaries,
  permutation_confidence_level,
  permutation_seed = NULL
) {
  if (selection_method == "dds") {
    return(select_dds_candidate(candidates, data, weighted))
  }
  if (selection_method != "permutation") {
    score_name <- selection_score_name(selection_method)
    return(list(
      best = candidates[[which.min(vapply(
        candidates,
        `[[`,
        numeric(1),
        score_name
      ))]],
      p_value = NA_real_,
      p_value_lower = NA_real_,
      p_value_upper = NA_real_,
      confidence_level = NA_real_,
      n_permutations = NA_integer_,
      stopped_early = NA,
      reject = NA,
      tests = empty_selection_tests(),
      effective_max_joinpoints = max(vapply(
        candidates,
        `[[`,
        integer(1),
        "n_joinpoints"
      ))
    ))
  }

  by_n <- split(
    candidates,
    vapply(candidates, `[[`, integer(1), "n_joinpoints")
  )
  available <- sort(unique(vapply(
    candidates,
    `[[`,
    integer(1),
    "n_joinpoints"
  )))
  if (length(available) == 0L) {
    return(list(
      best = candidates[[1L]],
      p_value = NA_real_,
      p_value_lower = NA_real_,
      p_value_upper = NA_real_,
      confidence_level = NA_real_,
      n_permutations = NA_integer_,
      stopped_early = NA,
      reject = NA,
      tests = empty_selection_tests(),
      effective_max_joinpoints = NA_integer_
    ))
  }
  effective_max_joinpoints <- max(available)
  low_index <- match(min_joinpoints, available)
  if (is.na(low_index)) {
    low_index <- 1L
  }
  high_index <- length(available)
  if (low_index >= high_index) {
    return(list(
      best = by_n[[as.character(available[low_index])]][[1L]],
      p_value = NA_real_,
      p_value_lower = NA_real_,
      p_value_upper = NA_real_,
      confidence_level = NA_real_,
      n_permutations = NA_integer_,
      stopped_early = NA,
      reject = NA,
      tests = empty_selection_tests(),
      effective_max_joinpoints = effective_max_joinpoints
    ))
  }

  selection <- list(
    p_value = NA_real_,
    p_value_lower = NA_real_,
    p_value_upper = NA_real_,
    confidence_level = permutation_confidence_level,
    n_permutations = NA_integer_,
    stopped_early = FALSE,
    reject = FALSE
  )
  tests <- list()
  uniform_stream <- if (
    permutation_method == "fixed" && !is.null(permutation_seed)
  ) {
    nci_mt_stream(permutation_seed)
  } else {
    NULL
  }
  permutation_boundary_for <- function(null_joinpoints) {
    if (
      is.null(permutation_boundaries) ||
        effective_max_joinpoints != max_joinpoints
    ) {
      return(NULL)
    }
    permutation_boundaries[[as.character(null_joinpoints)]]
  }
  while (high_index - low_index > 1L) {
    k_null <- available[low_index]
    k_alt <- available[high_index]
    test_alpha <- alpha / (effective_max_joinpoints - k_null)
    permutation_boundary <- permutation_boundary_for(k_null)
    selection <- monte_carlo_p_value(
      reduced = by_n[[as.character(k_null)]][[1L]],
      complex = by_n[[as.character(k_alt)]][[1L]],
      data = data,
      weighted = weighted,
      n_permutations = n_permutations,
      alpha = test_alpha,
      method = permutation_method,
      alpha0 = permutation_alpha0,
      beta0 = permutation_beta0,
      p1 = permutation_p1,
      boundary = permutation_boundary,
      confidence_level = permutation_confidence_level,
      uniform_stream = uniform_stream
    )
    tests[[length(tests) + 1L]] <- selection_test_row(
      length(tests) + 1L,
      k_null,
      k_alt,
      test_alpha,
      selection
    )
    if (is.na(selection$p_value)) {
      break
    }
    if (isTRUE(selection$reject)) {
      low_index <- low_index + 1L
    } else {
      high_index <- high_index - 1L
    }
  }
  if (high_index - low_index == 1L) {
    k_null <- available[low_index]
    k_alt <- available[high_index]
    test_alpha <- alpha / (effective_max_joinpoints - k_null)
    permutation_boundary <- permutation_boundary_for(k_null)
    selection <- monte_carlo_p_value(
      reduced = by_n[[as.character(k_null)]][[1L]],
      complex = by_n[[as.character(k_alt)]][[1L]],
      data = data,
      weighted = weighted,
      n_permutations = n_permutations,
      alpha = test_alpha,
      method = permutation_method,
      alpha0 = permutation_alpha0,
      beta0 = permutation_beta0,
      p1 = permutation_p1,
      boundary = permutation_boundary,
      confidence_level = permutation_confidence_level,
      uniform_stream = uniform_stream
    )
    tests[[length(tests) + 1L]] <- selection_test_row(
      length(tests) + 1L,
      k_null,
      k_alt,
      test_alpha,
      selection
    )
  }
  selected_index <- if (isTRUE(selection$reject)) {
    high_index
  } else {
    low_index
  }
  list(
    best = by_n[[as.character(available[selected_index])]][[1L]],
    p_value = selection$p_value,
    p_value_lower = selection$p_value_lower,
    p_value_upper = selection$p_value_upper,
    confidence_level = selection$confidence_level,
    n_permutations = selection$n_permutations,
    stopped_early = selection$stopped_early,
    reject = selection$reject,
    tests = if (length(tests) == 0L) {
      empty_selection_tests()
    } else {
      do.call(rbind, tests)
    },
    effective_max_joinpoints = effective_max_joinpoints
  )
}

#' Internal helper: Apply NCI data-dependent BIC/BIC3 selection.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' National Cancer Institute. Data Dependent Selection (DDS). Joinpoint Help
#' System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/setting-parameters/method-and-parameters-tab/model-selection-method/data-dependent-selection}
select_dds_candidate <- function(candidates, data, weighted) {
  bic_index <- which.min(vapply(candidates, `[[`, numeric(1), "bic"))
  bic3_index <- which.min(vapply(candidates, `[[`, numeric(1), "bic3"))
  max_candidate <- candidates[[which.max(vapply(
    candidates,
    `[[`,
    integer(1),
    "n_joinpoints"
  ))]]
  # Include estimated joinpoint locations in the variance degrees of freedom.
  sigma_squared <- max_candidate$sse /
    (nrow(data) - 2L * max_candidate$n_joinpoints - 2L)
  deltas <- vapply(
    candidates,
    dds_minimum_delta,
    numeric(1),
    data = data,
    weighted = weighted,
    sigma_squared = sigma_squared
  )
  for (index in seq_along(candidates)) {
    candidates[[index]]$dds_delta <- deltas[index]
    candidates[[index]]$dds_score <- NA_real_
  }
  if (bic_index == bic3_index) {
    selected_index <- bic_index
    selected_method <- "BIC=BIC3"
  } else {
    selected_delta <- deltas[c(bic_index, bic3_index)]
    delta_max <- max(selected_delta)
    delta_min <- min(selected_delta)
    if (is.na(delta_max - delta_min)) {
      stop(
        "DDS cannot resolve BIC/BIC3 disagreement: local Delta statistics are undefined.",
        call. = FALSE
      )
    }
    use_bic <- delta_max <= 10 || delta_max - delta_min > 200
    selected_index <- if (use_bic) bic_index else bic3_index
    selected_method <- if (use_bic) "BIC" else "BIC3"
    decision_rule <- if (delta_max <= 10) {
      "max_delta_le_10"
    } else if (delta_max - delta_min > 200) {
      "delta_range_gt_200"
    } else {
      "default_BIC3"
    }
  }
  if (bic_index == bic3_index) {
    selected_delta <- deltas[c(bic_index, bic3_index)]
    delta_max <- max(selected_delta)
    delta_min <- min(selected_delta)
    decision_rule <- "BIC_and_BIC3_agree"
  }
  candidates[[selected_index]]$dds_score <- if (selected_method == "BIC3") {
    candidates[[selected_index]]$bic3
  } else {
    candidates[[selected_index]]$bic
  }
  list(
    best = candidates[[selected_index]],
    p_value = NA_real_,
    p_value_lower = NA_real_,
    p_value_upper = NA_real_,
    confidence_level = NA_real_,
    n_permutations = NA_integer_,
    stopped_early = NA,
    reject = NA,
    tests = empty_selection_tests(),
    effective_max_joinpoints = max_candidate$n_joinpoints,
    dds_method = selected_method,
    candidates = candidates,
    dds_sigma_squared = sigma_squared,
    dds_variance_df = nrow(data) - 2L * max_candidate$n_joinpoints - 2L,
    dds_bic_delta = deltas[bic_index],
    dds_bic3_delta = deltas[bic3_index],
    dds_delta_max = delta_max,
    dds_delta_min = delta_min,
    dds_decision_rule = decision_rule
  )
}

#' Internal helper: Calculate the DDS minimum local slope-change statistic.
#'
#' @rdname crtrend-internals
#' @keywords internal
dds_minimum_delta <- function(candidate, data, weighted, sigma_squared) {
  if (candidate$n_joinpoints == 0L) {
    return(0)
  }
  if (!is.finite(sigma_squared) || sigma_squared <= 0) {
    return(Inf)
  }
  joinpoints <- candidate$joinpoints
  coefficients <- stats::coef(candidate$model)
  changes <- coefficients[seq.int(3L, candidate$n_joinpoints + 2L)]
  values <- vapply(
    seq_along(joinpoints),
    function(index) {
      previous <- if (index == 1L) {
        min(data$trend_time) - 1
      } else {
        joinpoints[index - 1L]
      }
      following <- if (index == length(joinpoints)) {
        max(data$trend_time)
      } else {
        joinpoints[index + 1L]
      }
      keep <- data$trend_time > previous & data$trend_time <= following
      local_data <- data[keep, , drop = FALSE]
      x <- local_data$trend_time
      z <- pmax(x - joinpoints[index], 0)
      design <- joinpoint_whiten(local_data, cbind(1, x))
      z <- joinpoint_whiten(local_data, z)
      solver <- qr(design)
      z_residual <- z - qr.fitted(solver, z)
      changes[index]^2 * sum(z_residual^2) / sigma_squared
    },
    numeric(1)
  )
  min(values)
}

#' Internal helper: Calculate the Monte Carlo permutation p-value and interval.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Fay, M. P., Feuer, E. J., and Midthune, D. N. (2000).
#' Permutation tests for Joinpoint regression with applications to cancer
#' rates. *Statistics in Medicine*, 19(3), 335-351.
#'
#' Fay, M. P., Kim, H.-J., and Hachey, M. (2007). On using truncated
#' sequential probability ratio test boundaries for Monte Carlo implementation
#' of hypothesis tests. *Journal of Computational and Graphical Statistics*,
#' 16(4), 946-967.
monte_carlo_p_value <- function(
  reduced,
  complex,
  data,
  weighted,
  n_permutations,
  alpha,
  method,
  alpha0,
  beta0,
  p1,
  boundary = NULL,
  confidence_level = 0.99,
  uniform_stream = NULL
) {
  observed_reduced <- model_sse(reduced$model, data, weighted)
  observed_complex <- model_sse(complex$model, data, weighted)
  observed_statistic <- permutation_f_statistic(
    observed_reduced,
    observed_complex,
    reduced$n_joinpoints,
    complex$n_joinpoints,
    nrow(data)
  )
  unavailable <- function(n_used = 0L) {
    list(
      p_value = NA_real_,
      statistic = observed_statistic,
      p_value_lower = NA_real_,
      p_value_upper = NA_real_,
      confidence_level = confidence_level,
      n_permutations = as.integer(n_used),
      stopped_early = n_used < n_permutations,
      reject = NA
    )
  }
  if (!is.finite(observed_statistic)) {
    return(unavailable())
  }
  fitted <- stats::fitted(reduced$model)
  residuals <- stats::residuals(reduced$model)
  standardized_residuals <- as.numeric(joinpoint_whiten(data, residuals))
  if (method == "fixed") {
    shuffle <- if (!is.null(uniform_stream)) {
      nci_permutation_shuffler(standardized_residuals, uniform_stream)
    } else {
      NULL
    }
    batch_size <- permutation_batch_size(reduced, complex)
    null_statistics <- numeric(n_permutations)
    first <- 1L
    while (first <= n_permutations) {
      last <- min(n_permutations, first + batch_size - 1L)
      n_batch <- last - first + 1L
      permuted_residuals <- if (is.null(shuffle)) {
        replicate(
          n_batch,
          sample(standardized_residuals, replace = FALSE)
        )
      } else {
        shuffle(n_batch)
      }
      permuted_log <- fitted + joinpoint_unwhiten(data, permuted_residuals)
      reduced_sse <- candidate_sse_batch(
        reduced,
        permuted_log,
        data,
        weighted
      )
      complex_sse <- candidate_sse_batch(
        complex,
        permuted_log,
        data,
        weighted
      )
      statistic <- permutation_f_statistic(
        reduced_sse,
        complex_sse,
        reduced$n_joinpoints,
        complex$n_joinpoints,
        nrow(data)
      )
      if (any(!is.finite(statistic))) {
        return(unavailable(last))
      }
      null_statistics[first:last] <- statistic
      first <- last + 1L
    }
    if (length(null_statistics) == 0L) {
      return(list(
        p_value = NA_real_,
        statistic = observed_statistic,
        p_value_lower = NA_real_,
        p_value_upper = NA_real_,
        confidence_level = confidence_level,
        n_permutations = 0L,
        stopped_early = FALSE,
        reject = NA
      ))
    }
    successes <- sum(null_statistics >= observed_statistic)
    n_used <- length(null_statistics)
    p_value <- (1 + successes) / (n_used + 1)
    p_value_interval <- stats::binom.test(
      successes,
      n_used,
      conf.level = confidence_level
    )$conf.int
    return(list(
      p_value = p_value,
      statistic = observed_statistic,
      p_value_lower = p_value_interval[1L],
      p_value_upper = p_value_interval[2L],
      confidence_level = confidence_level,
      n_permutations = n_used,
      stopped_early = FALSE,
      reject = p_value <= alpha
    ))
  }

  if (is.null(boundary)) {
    boundary <- make_tsprt_boundary(
      n_permutations = n_permutations,
      alpha = alpha,
      alpha0 = alpha0,
      beta0 = beta0,
      p1 = p1
    )
  }
  successes <- 0L
  for (iteration in seq_len(n_permutations)) {
    permuted_residuals <- sample(standardized_residuals, replace = FALSE)
    permuted_log <- fitted + joinpoint_unwhiten(data, permuted_residuals)
    reduced_sse <- candidate_sse(reduced, permuted_log, data, weighted)
    complex_sse <- candidate_sse(complex, permuted_log, data, weighted)
    null_statistic <- permutation_f_statistic(
      reduced_sse,
      complex_sse,
      reduced$n_joinpoints,
      complex$n_joinpoints,
      nrow(data)
    )
    if (!is.finite(null_statistic)) {
      return(unavailable(iteration))
    }
    if (is.finite(null_statistic) && null_statistic >= observed_statistic) {
      successes <- successes + 1L
    }
    stopping <- tsprt_stopping_decision(
      successes,
      iteration,
      boundary
    )
    if (!is.null(stopping)) {
      p_value <- tsprt_valid_p_value(
        successes,
        iteration,
        boundary
      )
      p_value_interval <- tsprt_p_value_interval(
        successes,
        iteration,
        boundary,
        confidence_level
      )
      return(list(
        p_value = p_value,
        statistic = observed_statistic,
        p_value_lower = p_value_interval[1L],
        p_value_upper = p_value_interval[2L],
        confidence_level = confidence_level,
        n_permutations = iteration,
        stopped_early = iteration < n_permutations,
        reject = p_value <= alpha
      ))
    }
  }
  p_value <- tsprt_valid_p_value(successes, n_permutations, boundary)
  p_value_interval <- tsprt_p_value_interval(
    successes,
    n_permutations,
    boundary,
    confidence_level
  )
  list(
    p_value = p_value,
    statistic = observed_statistic,
    p_value_lower = p_value_interval[1L],
    p_value_upper = p_value_interval[2L],
    confidence_level = confidence_level,
    n_permutations = n_permutations,
    stopped_early = FALSE,
    reject = p_value <= alpha
  )
}

#' Internal helper: Calculate the Joinpoint permutation F-type statistic.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Fay, M. P., Feuer, E. J., and Midthune, D. N. (2000).
#' Permutation tests for Joinpoint regression with applications to cancer
#' rates. *Statistics in Medicine*, 19(3), 335-351.
permutation_f_statistic <- function(
  reduced_sse,
  complex_sse,
  reduced_joinpoints,
  complex_joinpoints,
  n_observations
) {
  numerator_df <- 2 * (complex_joinpoints - reduced_joinpoints)
  denominator_df <- n_observations - 2 * complex_joinpoints - 2
  result <- rep(NA_real_, max(length(reduced_sse), length(complex_sse)))
  valid <- is.finite(reduced_sse) &
    is.finite(complex_sse) &
    complex_sse > 0 &
    numerator_df > 0 &
    denominator_df > 0
  result[valid] <- pmax(
    0,
    ((reduced_sse[valid] - complex_sse[valid]) / numerator_df) /
      (complex_sse[valid] / denominator_df)
  )
  result
}

#' Internal helper: Choose a memory-bounded permutation matrix batch size.
#'
#' @rdname crtrend-internals
#' @keywords internal
permutation_batch_size <- function(reduced, complex) {
  basis_columns <- vapply(
    list(reduced, complex),
    function(candidate) {
      if (is.null(candidate$grid_solver)) {
        return(ncol(candidate$design))
      }
      ncol(candidate$grid_solver$basis)
    },
    integer(1)
  )
  max(
    1L,
    min(128L, as.integer(floor(2e6 / max(1L, basis_columns))))
  )
}

#' Internal helper: Build the truncated sequential probability ratio test boundary.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Fay, M. P., Kim, H.-J., and Hachey, M. (2007). On using truncated
#' sequential probability ratio test boundaries for Monte Carlo implementation
#' of hypothesis tests. *Journal of Computational and Graphical Statistics*,
#' 16(4), 946-967.
make_tsprt_boundary <- function(
  n_permutations,
  alpha,
  alpha0,
  beta0,
  p1 = NULL
) {
  if (is.null(p1)) {
    p1 <- 0.8 * alpha
  }
  p0 <- tsprt_p0_given_p1(p1, alpha)
  likelihood_ratio <- p1 * (1 - p0) / (p0 * (1 - p1))
  a <- (1 - beta0) / alpha0
  b <- beta0 / (1 - alpha0)
  c0 <- log((1 - p0) / (1 - p1)) / log(likelihood_ratio)
  c1 <- log(b) / log(likelihood_ratio)
  c2 <- log(a) / log(likelihood_ratio)

  lower <- floor(c2 + seq_len(n_permutations) * c0)
  upper <- ceiling(c1 + seq_len(n_permutations) * c0)
  continuation_log_mass <- 0
  boundary_s <- integer()
  boundary_n <- integer()
  boundary_log_kstar <- numeric()
  boundary_log_path_count <- numeric()

  for (n in seq_len(n_permutations)) {
    successes <- 0:n
    log_mass <- rep(-Inf, n + 1L)
    if (n == 1L) {
      log_mass <- rep.int(-log(2), 2L)
    } else {
      log_mass[1L] <- continuation_log_mass[1L] - log(2)
      log_mass[n + 1L] <- continuation_log_mass[n] - log(2)
      log_mass[2:n] <- logspace_add(
        continuation_log_mass[1:(n - 1L)] - log(2),
        continuation_log_mass[2:n] - log(2)
      )
    }
    is_boundary <- if (n == n_permutations) {
      is.finite(log_mass)
    } else {
      successes %in% c(lower[n], upper[n])
    }
    is_boundary <- is_boundary & is.finite(log_mass)
    if (any(is_boundary)) {
      boundary_s <- c(boundary_s, successes[is_boundary])
      boundary_n <- c(boundary_n, rep.int(n, sum(is_boundary)))
      boundary_log_kstar <- c(
        boundary_log_kstar,
        log_mass[is_boundary] +
          n * log(2) +
          lbeta(successes[is_boundary] + 1, n - successes[is_boundary] + 1)
      )
      boundary_log_path_count <- c(
        boundary_log_path_count,
        log_mass[is_boundary] + n * log(2)
      )
    }
    continuation_log_mass <- log_mass
    if (n < n_permutations) {
      continuation_log_mass[is_boundary] <- -Inf
    }
  }

  kstar <- exp(boundary_log_kstar - max(boundary_log_kstar))
  kstar <- kstar / sum(kstar)
  list(
    s = boundary_s,
    n = boundary_n,
    ratio = boundary_s / boundary_n,
    kstar = kstar,
    log_path_count = boundary_log_path_count,
    p0 = p0,
    p1 = p1,
    c0 = c0,
    lower = lower,
    upper = upper
  )
}

#' Internal helper: Add two values on the log scale.
#'
#' @rdname crtrend-internals
#' @keywords internal
logspace_add <- function(x, y) {
  maximum <- pmax(x, y)
  result <- maximum
  finite <- is.finite(maximum)
  result[finite] <- maximum[finite] +
    log(
      exp(x[finite] - maximum[finite]) +
        exp(y[finite] - maximum[finite])
    )
  result[!finite] <- -Inf
  result
}

#' Internal helper: Solve the tSPRT null probability boundary.
#'
#' @rdname crtrend-internals
#' @keywords internal
tsprt_p0_given_p1 <- function(p1, alpha) {
  lower <- alpha + sqrt(.Machine$double.eps)
  upper <- 1 - sqrt(.Machine$double.eps)
  objective <- function(p0) {
    ratio <- p1 * (1 - p0) / (p0 * (1 - p1))
    log((1 - p0) / (1 - p1)) / log(ratio) - alpha
  }
  root <- tryCatch(
    stats::uniroot(objective, lower = lower, upper = upper)$root,
    error = function(error) NA_real_
  )
  if (!is.finite(root)) {
    root <- (alpha + 1) / 2
  }
  root
}

#' Internal helper: Determine whether a tSPRT path has stopped.
#'
#' @rdname crtrend-internals
#' @keywords internal
tsprt_stopping_decision <- function(successes, n, boundary) {
  if (n >= max(boundary$n)) {
    return("truncated")
  }
  lower <- boundary$lower[n]
  upper <- boundary$upper[n]
  if (successes <= lower) {
    return("reject")
  }
  if (successes >= upper) {
    return("accept")
  }
  NULL
}

#' Internal helper: Calculate the valid tSPRT p-value at a stopping boundary.
#'
#' @rdname crtrend-internals
#' @keywords internal
tsprt_valid_p_value <- function(successes, n, boundary) {
  index <- which(boundary$s == successes & boundary$n == n)
  if (length(index) != 1L) {
    return(NA_real_)
  }
  value <- sum(boundary$kstar[boundary$ratio <= boundary$ratio[index]])
  max(0, min(1, value))
}

#' Internal helper: Sum probabilities stably on the log scale.
#'
#' @rdname crtrend-internals
#' @keywords internal
logspace_sum <- function(x) {
  finite <- is.finite(x)
  if (!any(finite)) {
    return(-Inf)
  }
  maximum <- max(x[finite])
  maximum + log(sum(exp(x[finite] - maximum)))
}

#' Internal helper: Calculate a binomial log probability.
#'
#' @rdname crtrend-internals
#' @keywords internal
binomial_log_probability <- function(successes, trials, probability) {
  if (probability == 0) {
    return(ifelse(successes == 0, 0, -Inf))
  }
  if (probability == 1) {
    return(ifelse(successes == trials, 0, -Inf))
  }
  successes * log(probability) + (trials - successes) * log1p(-probability)
}

#' Internal helper: Calculate the tSPRT boundary cumulative distribution function.
#'
#' @rdname crtrend-internals
#' @keywords internal
tsprt_boundary_cdf <- function(probability, boundary, ratio_cutoff) {
  keep <- boundary$ratio <= ratio_cutoff
  if (!any(keep)) {
    return(0)
  }
  log_probability <- boundary$log_path_count[keep] +
    binomial_log_probability(
      boundary$s[keep],
      boundary$n[keep],
      probability
    )
  exp(logspace_sum(log_probability))
}

#' Internal helper: Calculate the Monte Carlo p-value interval for tSPRT.
#'
#' @rdname crtrend-internals
#' @keywords internal
tsprt_p_value_interval <- function(
  successes,
  n,
  boundary,
  confidence_level
) {
  index <- which(boundary$s == successes & boundary$n == n)
  if (length(index) != 1L) {
    return(c(NA_real_, NA_real_))
  }
  ratio_cutoff <- boundary$ratio[index]
  cdf <- function(probability) {
    tsprt_boundary_cdf(probability, boundary, ratio_cutoff)
  }
  lower_target <- 1 - (1 - confidence_level) / 2
  upper_target <- (1 - confidence_level) / 2
  solve_target <- function(target) {
    at_zero <- cdf(0)
    at_one <- cdf(1)
    if (at_zero < target) {
      return(0)
    }
    if (at_one > target) {
      return(1)
    }
    tryCatch(
      stats::uniroot(
        function(probability) cdf(probability) - target,
        interval = c(0, 1),
        tol = 1e-10
      )$root,
      error = function(error) NA_real_
    )
  }
  c(solve_target(lower_target), solve_target(upper_target))
}

#' Internal helper: Calculate a candidate model SSE for a resampled response.
#'
#' @rdname crtrend-internals
#' @keywords internal
candidate_sse <- function(candidate, response, data, weighted) {
  grid_solver <- candidate$grid_solver
  if (!is.null(grid_solver)) {
    return(candidate_sse_batch(
      candidate,
      matrix(response, ncol = 1L),
      data,
      weighted
    )[1L])
  }
  solver <- candidate$solver
  if (!is.null(solver)) {
    transformed_response <- joinpoint_whiten(data, response)
    residuals <- transformed_response -
      qr.fitted(
        solver$qr,
        transformed_response
      )
    return(sum(residuals^2))
  }
  permuted_data <- data
  permuted_data$trend_log <- response
  model <- if (candidate$n_joinpoints == 0L) {
    fit_lm_model(permuted_data, weighted = weighted)
  } else {
    fit_fixed_joinpoint_model(
      candidate$joinpoints,
      permuted_data,
      weighted = weighted
    )
  }
  model_sse(model, permuted_data, weighted)
}

#' Internal helper: Calculate minimum-grid SSEs for response columns.
#'
#' @rdname crtrend-internals
#' @keywords internal
candidate_sse_batch <- function(candidate, response, data, weighted) {
  response <- as.matrix(response)
  grid_solver <- candidate$grid_solver
  if (is.null(grid_solver)) {
    return(vapply(
      seq_len(ncol(response)),
      function(index) {
        candidate_sse(candidate, response[, index], data, weighted)
      },
      numeric(1)
    ))
  }
  transformed_response <- joinpoint_whiten(data, response)
  projected <- crossprod(grid_solver$basis, transformed_response)
  projection_squares <- colSums(matrix(
    projected^2,
    nrow = grid_solver$n_coefficients
  ))
  projection_squares <- matrix(
    projection_squares,
    nrow = grid_solver$n_models,
    ncol = ncol(response)
  )
  residual_sum_squares <- rep(
    colSums(transformed_response^2),
    each = grid_solver$n_models
  ) -
    projection_squares
  residual_sum_squares <- matrix(
    residual_sum_squares,
    nrow = grid_solver$n_models,
    ncol = ncol(response)
  )
  threshold <- 64 * .Machine$double.eps * colSums(transformed_response^2)
  unstable <- which(
    !is.finite(residual_sum_squares) |
      residual_sum_squares <= rep(threshold, each = grid_solver$n_models),
    arr.ind = TRUE
  )
  for (i in seq_len(nrow(unstable))) {
    model_index <- unstable[i, 1L]
    response_index <- unstable[i, 2L]
    columns <- (model_index - 1L) *
      grid_solver$n_coefficients +
      seq_len(grid_solver$n_coefficients)
    basis <- grid_solver$basis[, columns, drop = FALSE]
    y <- transformed_response[, response_index]
    residual <- y - as.numeric(basis %*% crossprod(basis, y))
    residual_sum_squares[model_index, response_index] <- sum(residual^2)
  }
  if (grid_solver$n_models == 1L) {
    return(as.numeric(residual_sum_squares))
  }
  apply(residual_sum_squares, 2L, min)
}

#' Internal helper: Calculate the weighted or unweighted model SSE.
#'
#' @rdname crtrend-internals
#' @keywords internal
model_sse <- function(model, data, weighted) {
  residuals <- stats::residuals(model)
  sum(joinpoint_whiten(data, residuals)^2)
}

#' Internal helper: Generate the admissible observed integer-year joinpoint grid.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Lerman, P. M. (1980). Fitting segmented regression models by grid search.
#' *Applied Statistics*, 29(1), 77-84.
integer_joinpoint_grid <- function(time, n_joinpoints, min_points = 2L) {
  time <- sort(unique(time))
  n_time <- length(time)
  if (any(abs(time - round(time)) > sqrt(.Machine$double.eps))) {
    stop("`time` must contain integer-valued years.", call. = FALSE)
  }
  if (n_joinpoints == 0L) {
    return(list(numeric()))
  }
  minimum_time_points <- (as.double(n_joinpoints) + 1) *
    min_points +
    n_joinpoints
  if (n_time < minimum_time_points) {
    return(NULL)
  }
  positions <- seq.int(min_points + 1L, n_time - min_points)
  if (length(positions) < n_joinpoints) {
    return(NULL)
  }
  combinations <- utils::combn(positions, n_joinpoints, simplify = FALSE)
  combinations <- Filter(
    function(position) {
      counts <- c(
        position[1L] - 1L,
        diff(position) - 1L,
        n_time - position[length(position)]
      )
      all(counts >= min_points)
    },
    combinations
  )
  lapply(combinations, function(position) time[position])
}

#' Internal helper: Reuse integer-year grids within one grouped fit.
#'
#' @rdname crtrend-internals
#' @keywords internal
cached_integer_joinpoint_grid <- function(
  time,
  n_joinpoints,
  min_points = 2L,
  cache = NULL
) {
  if (is.null(cache)) {
    return(integer_joinpoint_grid(time, n_joinpoints, min_points))
  }
  normalized_time <- sort(unique(time))
  key <- paste(
    min_points,
    n_joinpoints,
    paste(normalized_time, collapse = ","),
    sep = "|"
  )
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  result <- integer_joinpoint_grid(
    normalized_time,
    n_joinpoints,
    min_points
  )
  assign(key, result, envir = cache)
  result
}

#' Internal helper: Fit a zero-joinpoint log-linear model.
#'
#' @rdname crtrend-internals
#' @keywords internal
fit_lm_model <- function(data, weighted) {
  fit_gls_design(
    joinpoint_design(numeric(), data),
    data$trend_log,
    data
  )
}

#' Internal helper: Fit a linear model after covariance whitening.
#'
#' @rdname crtrend-internals
#' @keywords internal
fit_gls_design <- function(design, response, data) {
  transformed_design <- joinpoint_whiten(data, design)
  transformed_response <- as.numeric(joinpoint_whiten(data, response))
  fit <- tryCatch(
    stats::lm.fit(transformed_design, transformed_response),
    error = function(error) NULL
  )
  if (is.null(fit) || fit$rank < ncol(design)) {
    return(NULL)
  }
  coefficients <- fit$coefficients
  names(coefficients) <- colnames(design)
  fitted_values <- as.numeric(design %*% coefficients)
  fit$coefficients <- coefficients
  fit$fitted.values <- fitted_values
  fit$residuals <- as.numeric(response) - fitted_values
  fit$effects <- as.numeric(fit$effects)
  fit$call <- match.call()
  fit$terms <- stats::terms(trend_log ~ trend_time)
  fit$model <- data.frame(
    trend_log = as.numeric(response),
    trend_time = data$trend_time
  )
  covariance <- joinpoint_covariance(data)
  off_diagonal <- covariance
  diag(off_diagonal) <- 0
  fit$weights <- if (all(abs(off_diagonal) <= sqrt(.Machine$double.eps))) {
    unname(1 / diag(covariance))
  } else {
    NULL
  }
  fit$covariance <- covariance
  fit$df.residual <- nrow(design) - fit$rank
  class(fit) <- "lm"
  fit
}

#' Internal helper: Fit a continuous segmented log-linear model.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Hudson, D. J. (1966). Fitting segmented curves whose join points have to be
#' estimated. *Journal of the American Statistical Association*, 61(316),
#' 1097-1129.
#'
#' Feder, P. I. (1975). On asymptotic distribution theory in segmented
#' regression problems: Identified case. *The Annals of Statistics*, 3(1),
#' 49-83.
fit_fixed_joinpoint_model <- function(joinpoints, data, weighted) {
  model <- fit_gls_design(
    joinpoint_design(joinpoints, data),
    data$trend_log,
    data
  )
  if (is.null(model)) {
    return(NULL)
  }
  attr(model, "crtrend_joinpoints") <- joinpoints
  model
}

#' Internal helper: Calculate model fit and BIC-family scores.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Yu, B., and Feuer, E. J. (2009). Selecting the number of
#' change-points in segmented line regression. *Statistica Sinica*, 19(2),
#' 597-609. This reference supports the classical BIC comparison used here.
#'
#' Kim, J., and Kim, H.-J. (2016). Consistent model selection in segmented
#' line regression. *Journal of Statistical Planning and Inference*, 170,
#' 106-116. This reference supports the BIC3 penalty used here.
#'
#' Kim, H.-J., Chen, H.-S., Midthune, D., Wheeler, B., Buckman, D. W., Green,
#' D., Byrne, J., Luo, J., and Feuer, E. J. (2022). Data-driven choice of a
#' model selection method in Joinpoint regression. *Journal of Applied
#' Statistics*, 1-22. This reference supports the data-dependent choice of
#' the BIC-family method.
#'
#' National Cancer Institute. Bayesian Information Criterion (BIC) Details.
#' Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/tech-help/frequently-asked-questions/bayesian-information-criterion-bic-details}
#'
#' National Cancer Institute. BIC3. Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/setting-parameters/method-and-parameters-tab/model-selection-method/bic-3}
#'
#' National Cancer Institute. Weighted BIC (WBIC). Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/setting-parameters/method-and-parameters-tab/model-selection-method/weighted-bic-wbic}
fixed_joinpoint_scores <- function(
  joinpoints,
  data,
  weighted,
  calculate_partial = TRUE
) {
  model_data <- data
  design <- joinpoint_design(joinpoints, model_data)
  transformed_design <- joinpoint_whiten(model_data, design)
  transformed_response <- as.numeric(joinpoint_whiten(
    model_data,
    model_data$trend_log
  ))
  fit <- tryCatch(
    stats::lm.fit(transformed_design, transformed_response),
    error = function(error) NULL
  )
  if (
    is.null(fit) ||
      fit$rank < ncol(design) ||
      any(!is.finite(fit$coefficients))
  ) {
    return(invalid_joinpoint_scores())
  }
  residual_sum_squares <- sum(fit$residuals^2)
  if (!is.finite(residual_sum_squares)) {
    return(invalid_joinpoint_scores())
  }
  hinge_coefficients <- fit$coefficients[-seq_len(2L)]
  if (
    length(hinge_coefficients) > 0L &&
      all(abs(hinge_coefficients) <= sqrt(.Machine$double.eps))
  ) {
    return(invalid_joinpoint_scores(residual_sum_squares))
  }
  n_observations <- nrow(design)
  n_joinpoints <- length(joinpoints)
  log_n <- log(n_observations)
  fit_term <- log(residual_sum_squares / n_observations)
  information <- crossprod(transformed_design)
  log_determinant <- determinant(information, logarithm = TRUE)
  gamma_argument <- (n_observations - n_joinpoints - 3) / 2
  modified_bic <- if (
    log_determinant$sign > 0 &&
      is.finite(log_determinant$modulus) &&
      gamma_argument > 0 &&
      residual_sum_squares > 0
  ) {
    fit_term +
      (2 * n_joinpoints) * log_n / n_observations +
      as.numeric(log_determinant$modulus) / n_observations -
      2 * lgamma(gamma_argument) / n_observations -
      (n_joinpoints + 3) * log(residual_sum_squares) / n_observations
  } else {
    Inf
  }
  partial_r_squared <- if (calculate_partial) {
    max_adjacent_partial_r_squared(joinpoints, data, weighted)
  } else {
    0
  }
  bic <- fit_term + (2 * n_joinpoints + 2) * log_n / n_observations
  bic3 <- fit_term + (3 * n_joinpoints + 2) * log_n / n_observations
  weighted_bic <- if (
    identical(
      attr(data, "crtrend_wbic_weight", exact = TRUE),
      "joinpoint_6.1"
    )
  ) {
    if (n_joinpoints == 0L) {
      bic
    } else {
      bic * (1 - partial_r_squared) + bic3 * partial_r_squared
    }
  } else {
    fit_term +
      (2 * n_joinpoints + 2 + n_joinpoints * partial_r_squared) *
        log_n /
        n_observations
  }
  list(
    sse = residual_sum_squares,
    bic = bic,
    bic3 = bic3,
    modified_bic = modified_bic,
    weighted_bic = weighted_bic,
    partial_r_squared = partial_r_squared
  )
}

#' Internal helper: Return sentinel scores for an invalid candidate.
#'
#' @rdname crtrend-internals
#' @keywords internal
invalid_joinpoint_scores <- function(sse = Inf) {
  list(
    sse = sse,
    bic = Inf,
    bic3 = Inf,
    modified_bic = Inf,
    weighted_bic = Inf,
    partial_r_squared = NA_real_
  )
}

#' Internal helper: Construct the continuous segmented-regression design matrix.
#'
#' @rdname crtrend-internals
#' @keywords internal
joinpoint_design <- function(joinpoints, data) {
  hinge_names <- if (length(joinpoints) > 0L) {
    paste0("joinpoint_", seq_along(joinpoints))
  } else {
    character()
  }
  model_data <- data
  for (i in seq_along(joinpoints)) {
    model_data[[hinge_names[i]]] <- pmax(
      model_data$trend_time - joinpoints[i],
      0
    )
  }
  design <- as.matrix(model_data[c("trend_time", hinge_names)])
  cbind(`(Intercept)` = 1, design)
}

#' Internal helper: Calculate the maximum adjacent partial R-squared.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Chen, H.-S., Midthune, D., Wheeler, B., Buckman, D. W., Green,
#' D., Byrne, J., Luo, J., and Feuer, E. J. (2022). Data-driven choice of a
#' model selection method in Joinpoint regression. *Journal of Applied
#' Statistics*, 1-22.
max_adjacent_partial_r_squared <- function(joinpoints, data, weighted) {
  if (length(joinpoints) == 0L) {
    return(0)
  }
  values <- vapply(
    seq_along(joinpoints),
    partial_joinpoint_r_squared,
    numeric(1),
    joinpoints = joinpoints,
    data = data,
    weighted = weighted
  )
  if (all(!is.finite(values))) {
    return(0)
  }
  max(0, min(1, max(values, na.rm = TRUE)))
}

#' Internal helper: Calculate the partial R-squared for one joinpoint.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Chen, H.-S., Midthune, D., Wheeler, B., Buckman, D. W., Green,
#' D., Byrne, J., Luo, J., and Feuer, E. J. (2022). Data-driven choice of a
#' model selection method in Joinpoint regression. *Journal of Applied
#' Statistics*, 1-22.
partial_joinpoint_r_squared <- function(index, joinpoints, data, weighted) {
  previous <- if (index == 1L) {
    min(data$trend_time) - 1
  } else {
    joinpoints[index - 1L]
  }
  following <- if (index == length(joinpoints)) {
    max(data$trend_time)
  } else {
    joinpoints[index + 1L]
  }
  keep <- data$trend_time > previous & data$trend_time <= following
  if (sum(keep) < 3L) {
    return(0)
  }
  local_data <- data[keep, , drop = FALSE]
  if (
    !weighted &&
      is.null(attr(local_data, "crtrend_covariance", exact = TRUE))
  ) {
    local_data$trend_weight <- NULL
  }
  x <- local_data$trend_time
  y <- local_data$trend_log
  z <- pmax(x - joinpoints[index], 0)
  if (
    identical(
      attr(data, "crtrend_wbic_weight", exact = TRUE),
      "joinpoint_6.1"
    )
  ) {
    return(joinpoint_61_partial_weight(local_data, z))
  }
  base_design <- joinpoint_whiten(local_data, cbind(1, x))
  y <- as.numeric(joinpoint_whiten(local_data, y))
  z <- as.numeric(joinpoint_whiten(local_data, z))
  solver <- tryCatch(qr(base_design), error = function(error) NULL)
  if (is.null(solver) || solver$rank < 2L) {
    return(0)
  }
  y_residual <- y - qr.fitted(solver, y)
  z_residual <- z - qr.fitted(solver, z)
  yy <- sum(y_residual^2)
  zz <- sum(z_residual^2)
  yz <- sum(y_residual * z_residual)
  denominator <- yy * zz
  if (!is.finite(denominator) || denominator <= 0) {
    return(0)
  }
  value <- yz^2 / denominator
  max(0, min(1, value))
}

#' Internal helper: Create internal metadata for a fitted candidate.
#'
#' @rdname crtrend-internals
#' @keywords internal
make_candidate <- function(
  model,
  n_joinpoints,
  joinpoints = NULL,
  data,
  weighted,
  grid_joinpoints = NULL
) {
  if (is.null(model)) {
    return(NULL)
  }
  residual_df <- tryCatch(stats::df.residual(model), error = function(error) {
    NA_real_
  })
  coefficients <- tryCatch(stats::coef(model), error = function(error) NULL)
  if (
    !is.finite(residual_df) ||
      residual_df <= 0 ||
      is.null(coefficients) ||
      any(!is.finite(coefficients))
  ) {
    return(NULL)
  }
  if (is.null(joinpoints)) {
    joinpoints <- attr(model, "crtrend_joinpoints")
  }
  if (is.null(joinpoints)) {
    joinpoints <- if (n_joinpoints == 0L) {
      numeric()
    } else {
      extract_model_joinpoints(model)
    }
  }
  if (is.null(joinpoints)) {
    return(NULL)
  }
  scores <- fixed_joinpoint_scores(joinpoints, data, weighted)
  design <- joinpoint_design(joinpoints, data)
  solver <- tryCatch(
    {
      transformed_design <- joinpoint_whiten(data, design)
      fit <- qr(transformed_design)
      if (fit$rank < ncol(transformed_design)) {
        NULL
      } else {
        list(qr = fit)
      }
    },
    error = function(error) NULL
  )
  if (is.null(solver)) {
    return(NULL)
  }
  grid_solver <- compile_joinpoint_grid_solver(
    grid_joinpoints,
    data,
    weighted
  )
  list(
    n_joinpoints = n_joinpoints,
    model = model,
    joinpoints = joinpoints,
    sse = scores$sse,
    bic = scores$bic,
    bic3 = scores$bic3,
    modified_bic = scores$modified_bic,
    weighted_bic = scores$weighted_bic,
    partial_r_squared = scores$partial_r_squared,
    dds_delta = NA_real_,
    dds_score = NA_real_,
    dds_sigma_squared = NA_real_,
    dds_variance_df = NA_integer_,
    dds_bic_delta = NA_real_,
    dds_bic3_delta = NA_real_,
    dds_delta_max = NA_real_,
    dds_delta_min = NA_real_,
    dds_decision_rule = NA_character_,
    design = design,
    solver = solver,
    grid_solver = grid_solver
  )
}

#' Internal helper: Compile all fixed-grid least-squares projections.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' Kim, H.-J., Fay, M. P., Feuer, E. J., and Midthune, D. N. (2000).
#' Permutation tests for Joinpoint regression with applications to cancer
#' rates. *Statistics in Medicine*, 19(3), 335-351.
compile_joinpoint_grid_solver <- function(grid_joinpoints, data, weighted) {
  if (is.null(grid_joinpoints) || length(grid_joinpoints) == 0L) {
    return(NULL)
  }
  n_coefficients <- length(grid_joinpoints[[1L]]) + 2L
  bases <- lapply(grid_joinpoints, function(joinpoints) {
    transformed_design <- joinpoint_whiten(
      data,
      joinpoint_design(joinpoints, data)
    )
    factorization <- qr(transformed_design)
    if (factorization$rank < ncol(transformed_design)) {
      return(NULL)
    }
    qr.Q(factorization, complete = FALSE)[,
      seq_len(n_coefficients),
      drop = FALSE
    ]
  })
  bases <- Filter(Negate(is.null), bases)
  if (length(bases) == 0L) {
    return(NULL)
  }
  list(
    basis = do.call(cbind, bases),
    n_coefficients = n_coefficients,
    n_models = length(bases)
  )
}

#' Internal helper: Extract joinpoint estimates from a model object.
#'
#' @rdname crtrend-internals
#' @keywords internal
extract_model_joinpoints <- function(model) {
  psi <- tryCatch(model$psi, error = function(error) NULL)
  if (is.null(psi) || !"Est." %in% colnames(psi)) {
    return(NULL)
  }
  sort(as.numeric(psi[, "Est."]))
}

#' Internal helper: Calculate NCI-style segment-wise standard errors and degrees of freedom.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' National Cancer Institute. Estimated Regression Coefficients (Beta).
#' Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/statistical-notes/statistics-related-to-the-k-joinpoint-model/estimated-regression-coefficients-beta}.
#'
#' Feder, P. I. (1975). On asymptotic distribution theory in segmented
#' regression problems: Identified case. *The Annals of Statistics*, 3(1),
#' 49-83.
#'
#' Lerman, P. M. (1980). Fitting segmented regression models by grid search.
#' *Applied Statistics*, 29(1), 77-84.
nci_joinpoint_inference <- function(fit) {
  engine <- joinpoint_engine(fit)
  data <- joinpoint_analysis_data(fit)
  if (is.null(engine) || is.null(data)) {
    return(NULL)
  }
  joinpoints <- sort(as.numeric(fit$joinpoints))
  if (length(joinpoints) > 0L) {
    data <- data[!(data$trend_time %in% joinpoints), , drop = FALSE]
  }
  segment <- findInterval(data$trend_time, joinpoints) + 1L
  n_segments <- length(joinpoints) + 1L
  slope_estimates <- nci_grid_slope_estimates(fit)
  if (length(slope_estimates) != n_segments) {
    return(NULL)
  }
  segment_counts <- tabulate(segment, nbins = n_segments)
  if (any(segment_counts < 2L)) {
    return(NULL)
  }
  design <- matrix(0, nrow = nrow(data), ncol = 2L * n_segments)
  for (index in seq_len(n_segments)) {
    rows <- segment == index
    columns <- 2L * index - c(1L, 0L)
    design[rows, columns[1L]] <- 1
    design[rows, columns[2L]] <- data$trend_time[rows]
  }
  transformed_design <- joinpoint_whiten(data, design)
  transformed_response <- joinpoint_whiten(data, data$trend_log)
  solver <- tryCatch(qr(transformed_design), error = function(error) NULL)
  if (is.null(solver) || solver$rank < ncol(transformed_design)) {
    return(NULL)
  }
  unconstrained_coefficients <- qr.coef(solver, transformed_response)
  unconstrained_residuals <- transformed_response -
    transformed_design %*% unconstrained_coefficients
  degrees_of_freedom <- nrow(data) - ncol(design)
  if (degrees_of_freedom <= 0L) {
    return(NULL)
  }
  residual_variance <- sum(unconstrained_residuals^2) / degrees_of_freedom
  covariance_scale <- if (
    isTRUE(attr(
      data,
      "crtrend_known_covariance",
      exact = TRUE
    ))
  ) {
    1
  } else {
    residual_variance
  }
  coefficient_covariance <- tryCatch(
    covariance_scale * chol2inv(chol(crossprod(transformed_design))),
    error = function(error) NULL
  )
  slope_columns <- 2L * seq_len(n_segments)
  if (
    is.null(coefficient_covariance) ||
      any(!is.finite(coefficient_covariance))
  ) {
    return(NULL)
  }
  slope_covariance <- coefficient_covariance[
    slope_columns,
    slope_columns,
    drop = FALSE
  ]
  slope_se <- sqrt(pmax(0, diag(slope_covariance)))
  list(
    joinpoints = joinpoints,
    slope_estimates = slope_estimates,
    coefficient_covariance = coefficient_covariance,
    slope_covariance = slope_covariance,
    slope_se = slope_se,
    df = degrees_of_freedom,
    effective_n = nrow(data),
    residual_variance = residual_variance
  )
}

#' Internal helper: Extract continuous grid-search segment slope estimates.
#'
#' @rdname crtrend-internals
#' @keywords internal
#' @references
#' National Cancer Institute. Estimated Regression Coefficients (Beta).
#' Joinpoint Help System.
#' \url{https://surveillance.cancer.gov/help/joinpoint/statistical-notes/statistics-related-to-the-k-joinpoint-model/estimated-regression-coefficients-beta}.
nci_grid_slope_estimates <- function(fit) {
  coefficients <- tryCatch(
    stats::coef(joinpoint_engine(fit)),
    error = function(error) NULL
  )
  if (is.null(coefficients) || !is.finite(coefficients[["trend_time"]])) {
    return(numeric())
  }
  joinpoints <- sort(as.numeric(fit$joinpoints))
  baseline <- coefficients[["trend_time"]]
  if (length(joinpoints) == 0L) {
    return(baseline)
  }
  delta <- coefficients[paste0("joinpoint_", seq_along(joinpoints))]
  if (length(delta) != length(joinpoints) || any(!is.finite(delta))) {
    return(numeric())
  }
  c(baseline, baseline + cumsum(delta))
}

#' Internal helper: Create the result for a group that cannot be fitted.
#'
#' @rdname crtrend-internals
#' @keywords internal
empty_joinpoint_group <- function(
  group_key,
  group_vars,
  status,
  start,
  end,
  data
) {
  list(
    group_key = group_key,
    status = status,
    n_joinpoints = NA_integer_,
    joinpoints = numeric(),
    start = start,
    end = end,
    bic = NA_real_,
    bic3 = NA_real_,
    modified_bic = NA_real_,
    weighted_bic = NA_real_,
    partial_r_squared = NA_real_,
    dds_delta = NA_real_,
    dds_method = NA_character_,
    dds_sigma_squared = NA_real_,
    dds_variance_df = NA_integer_,
    dds_bic_delta = NA_real_,
    dds_bic3_delta = NA_real_,
    dds_delta_max = NA_real_,
    dds_delta_min = NA_real_,
    dds_decision_rule = NA_character_,
    selection_score = NA_real_,
    selection_p_value = NA_real_,
    selection_p_value_lower = NA_real_,
    selection_p_value_upper = NA_real_,
    selection_p_value_confidence_level = NA_real_,
    selection_n_permutations = NA_integer_,
    selection_stopped_early = NA,
    selection_tests = empty_selection_tests(),
    effective_max_joinpoints = NA_integer_,
    mc_p_value = NA_real_,
    mc_p_value_lower = NA_real_,
    mc_p_value_upper = NA_real_,
    mc_p_value_confidence_level = NA_real_,
    mc_n_permutations = NA_integer_,
    mc_stopped_early = NA,
    eq_draws = NULL,
    eq_n_resamples = NA_integer_,
    inference = NULL,
    data = data,
    model = NULL
  )
}

#' Internal helper: Convert one internal fit object to a model-information row.
#'
#' @rdname crtrend-internals
#' @keywords internal
fit_info_row <- function(fit) {
  row <- fit$group_key
  row$status <- fit$status
  row$error <- fit$error_structure %||% NA_character_
  row$n_joinpoints <- fit$n_joinpoints
  row$joinpoints <- if (fit$status == "ok") {
    format_joinpoints(fit$joinpoints)
  } else {
    NA_character_
  }
  row$start <- fit$start
  row$end <- fit$end
  row$bic <- fit$bic
  row$bic3 <- fit$bic3
  row$modified_bic <- fit$modified_bic
  row$weighted_bic <- fit$weighted_bic
  row$partial_r_squared <- fit$partial_r_squared
  row$dds_delta <- fit$dds_delta
  row$dds_method <- fit$dds_method
  row$dds_sigma_squared <- fit$dds_sigma_squared
  row$dds_variance_df <- fit$dds_variance_df
  row$dds_bic_delta <- fit$dds_bic_delta
  row$dds_bic3_delta <- fit$dds_bic3_delta
  row$dds_delta_max <- fit$dds_delta_max
  row$dds_delta_min <- fit$dds_delta_min
  row$dds_decision_rule <- fit$dds_decision_rule
  row$selection_score <- fit$selection_score
  row$selection_p_value <- fit$selection_p_value
  row$selection_p_value_lower <- fit$selection_p_value_lower
  row$selection_p_value_upper <- fit$selection_p_value_upper
  row$selection_p_value_confidence_level <-
    fit$selection_p_value_confidence_level
  row$selection_n_permutations <- fit$selection_n_permutations
  row$selection_stopped_early <- fit$selection_stopped_early
  row$effective_max_jp <- fit$effective_max_joinpoints
  row$mc_p_value <- fit$mc_p_value
  row$mc_p_value_lower <- fit$mc_p_value_lower
  row$mc_p_value_upper <- fit$mc_p_value_upper
  row$mc_p_value_confidence_level <- fit$mc_p_value_confidence_level
  row$mc_n_permutations <- fit$mc_n_permutations
  row$mc_stopped_early <- fit$mc_stopped_early
  row
}

#' Internal helper: Format selected joinpoints for model output.
#'
#' @rdname crtrend-internals
#' @keywords internal
format_joinpoints <- function(joinpoints) {
  if (length(joinpoints) == 0L) {
    return(NA_character_)
  }
  paste(format_period_value(joinpoints), collapse = ";")
}
