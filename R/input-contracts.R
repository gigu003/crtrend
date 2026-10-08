is_bounded_integer <- function(x, lower = 0, upper = .Machine$integer.max - 1) {
  is.numeric(x) &&
    !is.complex(x) &&
    length(x) == 1L &&
    is.finite(x) &&
    x >= lower &&
    x <= upper &&
    x == trunc(x)
}

validate_rng_seed <- function(seed, name = "seed") {
  if (
    !is.null(seed) &&
      !is_bounded_integer(seed, -.Machine$integer.max, .Machine$integer.max)
  ) {
    stop(
      "`",
      name,
      "` must be NULL or a finite integer in [-2147483647, 2147483647].",
      call. = FALSE
    )
  }
}

with_joinpoint_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (existed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(
    {
      if (existed) {
        assign(".Random.seed", previous, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(seed)
  force(code)
}

validate_joinpoint_group_names <- function(group_vars) {
  reserved <- c(
    "group_id",
    "status",
    "test_status",
    "n_joinpoints",
    "joinpoints",
    "segment",
    "period",
    "start",
    "end",
    "apc",
    "apc_lower",
    "apc_upper",
    "aapc",
    "aapc_lower",
    "aapc_upper",
    "aapc_se",
    "p_value",
    "statistic",
    "statistic_type",
    "index",
    "joinpoint",
    "slope",
    "bic",
    "bic3",
    "modified_bic",
    "weighted_bic",
    "partial_r_squared",
    "selection_score",
    "error",
    "effective_max_jp",
    "dds_delta",
    "dds_method",
    "row_id",
    "time",
    "observed",
    "fitted",
    "residual",
    "se",
    "weight",
    "weight_method",
    "n_observations",
    "n_parameters",
    "model_df",
    "sse",
    "mse",
    "selected",
    "parameterization",
    "term",
    "estimate",
    "df",
    "inference_df",
    "fixed_grid_df",
    "selection_method",
    "rho",
    "significant",
    "inference_status",
    "effective_n",
    "inference_residual_variance",
    "candidates_available",
    "linear_predictor",
    "link_residual",
    "apc_p_value",
    "apc_significant",
    "apc_method"
  )
  bad <- group_vars[
    group_vars %in% reserved | grepl("^(trend_|selection_|mc_)", group_vars)
  ]
  if (length(bad)) {
    stop(
      "Reserved grouping column name(s): ",
      paste(bad, collapse = ", "),
      ". Rename these columns before fitting.",
      call. = FALSE
    )
  }
}

validate_standard_age_groups <- function(
  data,
  groups,
  age,
  standard_population,
  expected_ages
) {
  observed_ages <- as.character(data[[age]])
  expected_ages <- expected_ages %||% unique(observed_ages)
  expected_ages <- as.character(expected_ages)
  if (
    !length(expected_ages) ||
      anyNA(expected_ages) ||
      anyDuplicated(expected_ages)
  ) {
    stop(
      "`expected_ages` must contain unique, non-missing age labels.",
      call. = FALSE
    )
  }
  reference <- NULL
  for (i in groups) {
    if (!setequal(observed_ages[i], expected_ages)) {
      stop(
        "Every group must contain exactly the expected age strata.",
        call. = FALSE
      )
    }
    weights <- data[[standard_population]][i][match(
      expected_ages,
      observed_ages[i]
    )]
    if (sum(weights) <= 0) {
      stop(
        "Each group must have a positive standard-population total.",
        call. = FALSE
      )
    }
    weights <- weights / sum(weights)
    if (
      !is.null(reference) &&
        !isTRUE(all.equal(weights, reference, tolerance = 1e-12))
    ) {
      stop(
        "Standard-population weights must be identical across groups.",
        call. = FALSE
      )
    }
    reference <- weights
  }
}

resolve_calc_resamples <- function(x, n_resamples) {
  if (!inherits(x, "crtrend_joinpoint")) {
    stop("`x` must be a `crtrend_joinpoint` object.", call. = FALSE)
  }
  n_resamples %||% x$control$n_eq %||% 5001L
}
