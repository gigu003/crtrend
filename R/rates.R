#' Calculate crude rates from counts and populations
#'
#' `calc_cr()` aggregates counts and populations within each requested
#' group and calculates a rate on the requested scale. The returned `se` is the
#' Poisson standard error on the same rate scale, so the result can be passed
#' directly to [joinpoint()] with `outcome = "rate"` and `se = "se"`.
#' `calc_crude_rate()` is retained as a compatibility alias.
#'
#' @param data A data frame containing counts, populations, and optional
#'   grouping variables.
#' @param count Name of the count column.
#' @param population Name of the population-at-risk column.
#' @param group_vars Optional character vector of grouping-variable names. If
#'   `NULL`, all columns other than `count` and `population` are used.
#' @param rate_per Positive rate multiplier. The default, `100000`, returns
#'   rates per 100,000 people.
#'
#' @return A data frame with the grouping columns, aggregated `count` and
#'   `population`, `rate`, and its Poisson `se`.
#'
#' @details
#' For a group with total count `C` and total population `P`, the returned
#' rate is `C / P * rate_per` and its standard error is
#' `sqrt(C) / P * rate_per`. Counts may be zero, but populations must be
#' strictly positive. Missing values and negative counts are not allowed.
#'
#' @examples
#' data(jp_crude)
#' crude <- calc_cr(
#'   jp_crude,
#'   count = "count",
#'   population = "population",
#'   group_vars = c("sex", "year")
#' )
#'
#' fit <- joinpoint(
#'   crude,
#'   time = "year",
#'   outcome = "rate",
#'   group_vars = "sex",
#'   se = "se",
#'   max_jp = 1
#' )
#'
#' @seealso [calc_asr()], [calc_prop()], [calc_percent()], [joinpoint()]
#' @export
calc_cr <- function(
  data,
  count = "count",
  population = "population",
  group_vars = NULL,
  rate_per = 100000
) {
  validated <- validate_rate_data(
    data = data,
    measure_vars = c(count, population),
    group_vars = group_vars,
    rate_per = rate_per,
    positive_vars = population
  )
  data <- validated$data
  groups <- rate_group_rows(data, validated$group_vars)
  counts <- vapply(groups, function(i) sum(data[[count]][i]), numeric(1))
  populations <- vapply(
    groups,
    function(i) sum(data[[population]][i]),
    numeric(1)
  )
  result <- rate_result_frame(data, validated$group_vars, groups)
  result$count <- counts
  result$population <- populations
  result$rate <- counts / populations * validated$rate_per
  result$se <- sqrt(counts) / populations * validated$rate_per
  rownames(result) <- NULL
  result
}

#' @rdname calc_cr
#' @export
calc_crude_rate <- function(
  data,
  count = "count",
  population = "population",
  group_vars = NULL,
  rate_per = 100000
) {
  calc_cr(
    data = data,
    count = count,
    population = population,
    group_vars = group_vars,
    rate_per = rate_per
  )
}

#' Calculate directly age-standardized rates from age-specific counts
#' and populations
#'
#' `calc_asr()` performs direct age standardization within each
#' requested group. It uses the supplied standard-population counts as weights
#' and returns a Poisson/delta-method standard error on the same rate scale.
#' The result can be passed directly to [joinpoint()] with
#' `outcome = "rate"` and `se = "se"`.
#' `calc_age_adjusted_rate()` is retained as a compatibility alias.
#'
#' @param data A data frame containing one row per group and age stratum.
#' @param count Name of the age-specific count column.
#' @param population Name of the age-specific population-at-risk column.
#' @param age Name of the age-stratum column.
#' @param expected_ages Optional complete vector of age labels. By default,
#'   uses the union of observed ages. Supply this explicitly to detect ages
#'   missing from every group. Restricted-age rates require explicitly
#'   subsetting all groups to the same age range.
#' @param std Name of the standard-population count column.
#' @param group_vars Optional character vector of grouping-variable names. If
#'   `NULL`, all columns other than `count`, `population`, `age`, and
#'   the column selected by `std` are used.
#' @param rate_per Positive rate multiplier. The default, `100000`, returns
#'   rates per 100,000 people.
#'
#' @return A data frame with the grouping columns, aggregated `count` and
#'   `population`, the directly age-standardized `rate`, and its `se`.
#'
#' @details
#' Within each group, the age-specific rates are combined as
#' `sum(w_i * count_i / population_i) * rate_per`, where
#' `w_i = std_i / sum(std_i)`. The returned
#' standard error is
#' `sqrt(sum(w_i^2 * count_i / population_i^2)) * rate_per`, the usual
#' independent-Poisson approximation for direct standardization. Each
#' group-age combination must occur at most once, and every group must have a
#' positive total standard population.
#' All groups must contain the same expected age strata and the same
#' normalized standard-population weights. Missing strata are not silently
#' renormalized. An age absent from every group cannot be detected without
#' supplying `expected_ages`.
#'
#' @examples
#' data(jp_age_adjusted)
#' age_adjusted <- calc_asr(
#'   jp_age_adjusted,
#'   count = "count",
#'   population = "population",
#'   age = "age_group",
#'   std = "standard_population",
#'   group_vars = c("sex", "year")
#' )
#'
#' fit <- joinpoint(
#'   age_adjusted,
#'   time = "year",
#'   outcome = "rate",
#'   group_vars = "sex",
#'   se = "se",
#'   max_jp = 1
#' )
#'
#' @seealso [calc_cr()], [calc_prop()], [calc_percent()], [joinpoint()]
#' @export
calc_asr <- function(
  data,
  count = "count",
  population = "population",
  age = "age",
  std = "std",
  group_vars = NULL,
  rate_per = 100000,
  expected_ages = NULL
) {
  if (
    !is.character(age) ||
      length(age) != 1L ||
      is.na(age) ||
      !nzchar(age)
  ) {
    stop("`age` must be one column name.", call. = FALSE)
  }
  if (
    !is.character(std) ||
      length(std) != 1L ||
      is.na(std) ||
      !nzchar(std)
  ) {
    stop("`std` must be one column name.", call. = FALSE)
  }
  if (!is.data.frame(data) || !age %in% names(data)) {
    stop("Missing age column: `", age, "`.", call. = FALSE)
  }
  measure_vars <- c(count, population, std)
  if (is.null(group_vars)) {
    group_vars <- setdiff(names(data), c(measure_vars, age))
  }
  validated <- validate_rate_data(
    data = data,
    measure_vars = measure_vars,
    group_vars = group_vars,
    rate_per = rate_per,
    positive_vars = population
  )
  data <- validated$data
  if (anyNA(data[[age]])) {
    stop("Missing values are not allowed in the age column.", call. = FALSE)
  }
  if (anyDuplicated(rate_group_key(data, c(validated$group_vars, age)))) {
    stop(
      "Each group must have at most one observation per age stratum.",
      call. = FALSE
    )
  }
  groups <- rate_group_rows(data, validated$group_vars)
  validate_standard_age_groups(
    data,
    groups,
    age,
    std,
    expected_ages
  )
  counts <- vapply(groups, function(i) sum(data[[count]][i]), numeric(1))
  populations <- vapply(
    groups,
    function(i) sum(data[[population]][i]),
    numeric(1)
  )
  rates <- vapply(
    groups,
    function(i) {
      weights <- data[[std]][i]
      weight_total <- sum(weights)
      if (weight_total <= 0) {
        stop(
          "Each group must have a positive standard-population total.",
          call. = FALSE
        )
      }
      weights <- weights / weight_total
      sum(weights * data[[count]][i] / data[[population]][i]) *
        validated$rate_per
    },
    numeric(1)
  )
  standard_errors <- vapply(
    groups,
    function(i) {
      weights <- data[[std]][i]
      weight_total <- sum(weights)
      if (weight_total <= 0) {
        stop(
          "Each group must have a positive standard-population total.",
          call. = FALSE
        )
      }
      weights <- weights / weight_total
      sqrt(sum(weights^2 * data[[count]][i] / data[[population]][i]^2)) *
        validated$rate_per
    },
    numeric(1)
  )
  result <- rate_result_frame(data, validated$group_vars, groups)
  result$count <- counts
  result$population <- populations
  result$rate <- rates
  result$se <- standard_errors
  rownames(result) <- NULL
  result
}

#' @rdname calc_asr
#' @export
calc_age_adjusted_rate <- function(
  data,
  count = "count",
  population = "population",
  age = "age",
  std = "std",
  group_vars = NULL,
  rate_per = 100000,
  expected_ages = NULL
) {
  calc_asr(
    data = data,
    count = count,
    population = population,
    age = age,
    std = std,
    group_vars = group_vars,
    rate_per = rate_per,
    expected_ages = expected_ages
  )
}

#' Calculate proportions or percentages from counts and populations
#'
#' `calc_prop()` and `calc_percent()` aggregate counts and populations within
#' each requested group and calculate a binomial proportion or percentage.
#' The returned `se` is on the same scale as the returned outcome, so the
#' result can be passed directly to [joinpoint()].
#'
#' @param data A data frame containing counts, populations, and optional
#'   grouping variables.
#' @param count Name of the count column.
#' @param population Name of the population or denominator column.
#' @param group_vars Optional character vector of grouping-variable names. If
#'   `NULL`, all columns other than `count` and `population` are used.
#'
#' @return `calc_prop()` returns the grouping columns, aggregated `count` and
#'   `population`, `prop`, and its binomial `se`. `calc_percent()` returns the
#'   same columns with `percent` replacing `prop`; its `se` is on the percentage
#'   scale.
#'
#' @details
#' For an aggregated count `C` and population `N`, `calc_prop()` returns
#' `p = C / N` and `sqrt(p * (1 - p) / N)`. `calc_percent()` multiplies both
#' quantities by 100. The normal approximation used by Joinpoint may be
#' unreliable when `C < 10` or `N - C < 10`; these functions issue a warning
#' when either condition occurs. For complex-survey estimates, calculate a
#' design-based standard error or covariance matrix externally instead.
#'
#' Boundary values zero and one are valid proportions and have a binomial
#' standard error of zero. They are returned without adjustment, but they are
#' not directly compatible with every log-linear [joinpoint()] analysis.
#'
#' @examples
#' d <- data.frame(
#'   year = 2020:2023,
#'   count = c(40, 44, 51, 57),
#'   population = c(100, 102, 105, 108)
#' )
#' proportions <- calc_prop(d)
#' percentages <- calc_percent(d)
#'
#' joinpoint(
#'   proportions,
#'   time = "year",
#'   outcome = "prop",
#'   se = "se",
#'   max_jp = 0
#' )
#'
#' @seealso [calc_cr()], [calc_asr()], [joinpoint()]
#' @export
calc_prop <- function(
  data,
  count = "count",
  population = "population",
  group_vars = NULL
) {
  calc_binomial_outcome(
    data = data,
    count = count,
    population = population,
    group_vars = group_vars,
    scale = "prop"
  )
}

#' @rdname calc_prop
#' @export
calc_percent <- function(
  data,
  count = "count",
  population = "population",
  group_vars = NULL
) {
  calc_binomial_outcome(
    data = data,
    count = count,
    population = population,
    group_vars = group_vars,
    scale = "percent"
  )
}

#' Internal helper: Calculate a grouped binomial outcome.
#'
#' @rdname crtrend-internals
#' @keywords internal
calc_binomial_outcome <- function(
  data,
  count,
  population,
  group_vars,
  scale
) {
  validated <- validate_rate_data(
    data = data,
    measure_vars = c(count, population),
    group_vars = group_vars,
    rate_per = 1,
    positive_vars = population
  )
  data <- validated$data
  if (any(data[[count]] > data[[population]])) {
    stop("`count` must not exceed `population`.", call. = FALSE)
  }
  groups <- rate_group_rows(data, validated$group_vars)
  counts <- vapply(groups, function(i) sum(data[[count]][i]), numeric(1))
  populations <- vapply(
    groups,
    function(i) sum(data[[population]][i]),
    numeric(1)
  )
  prop <- counts / populations
  standard_errors <- sqrt(prop * (1 - prop) / populations)
  multiplier <- if (scale == "percent") 100 else 1
  outcome_name <- if (scale == "percent") "percent" else "prop"
  result <- rate_result_frame(data, validated$group_vars, groups)
  result$count <- counts
  result$population <- populations
  result[[outcome_name]] <- multiplier * prop
  result$se <- multiplier * standard_errors
  rownames(result) <- NULL
  small <- counts < 10 | populations - counts < 10
  if (any(small)) {
    warning(
      "Normal-approximation standard errors may be unreliable for ",
      sum(small),
      " group(s) with `count < 10` or `population - count < 10`.",
      call. = FALSE
    )
  }
  result
}

#' Internal helper: Validate common rate-calculation arguments.
#'
#' @rdname crtrend-internals
#' @keywords internal
validate_rate_data <- function(
  data,
  measure_vars,
  group_vars,
  rate_per,
  positive_vars
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!nrow(data)) {
    stop("`data` must contain at least one row.", call. = FALSE)
  }
  if (anyDuplicated(names(data))) {
    stop("`data` must have unique column names.", call. = FALSE)
  }
  measure_vars <- validate_rate_names(measure_vars, "rate columns")
  if (!all(measure_vars %in% names(data))) {
    missing <- measure_vars[!measure_vars %in% names(data)]
    stop(
      "Missing rate column(s): ",
      paste0("`", missing, "`", collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(group_vars)) {
    group_vars <- setdiff(names(data), measure_vars)
  } else {
    group_vars <- validate_rate_names(group_vars, "`group_vars`")
    if (!all(group_vars %in% names(data))) {
      missing <- group_vars[!group_vars %in% names(data)]
      stop(
        "Missing grouping column(s): ",
        paste0("`", missing, "`", collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (
    any(
      group_vars %in% c("count", "population", "rate", "prop", "percent", "se")
    )
  ) {
    stop(
      paste0(
        "Grouping columns cannot be named `count`, `population`, `rate`, ",
        "`prop`, `percent`, or `se`."
      ),
      call. = FALSE
    )
  }
  if (
    length(rate_per) != 1L ||
      !is.numeric(rate_per) ||
      !is.finite(rate_per) ||
      rate_per <= 0
  ) {
    stop("`rate_per` must be one positive finite number.", call. = FALSE)
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  for (column in measure_vars) {
    data[[column]] <- validate_rate_numeric(data[[column]], column)
  }
  if (anyNA(data[group_vars])) {
    stop("Missing values are not allowed in grouping columns.", call. = FALSE)
  }
  for (column in positive_vars) {
    if (any(data[[column]] <= 0)) {
      stop("`", column, "` must contain positive values.", call. = FALSE)
    }
  }
  list(
    data = data,
    group_vars = group_vars,
    rate_per = as.numeric(rate_per)
  )
}

#' Internal helper: Validate a rate-calculation column name vector.
#'
#' @rdname crtrend-internals
#' @keywords internal
validate_rate_names <- function(x, label) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    stop(label, " must be unique, non-missing column names.", call. = FALSE)
  }
  x
}

#' Internal helper: Convert and validate a numeric rate-calculation column.
#'
#' @rdname crtrend-internals
#' @keywords internal
validate_rate_numeric <- function(x, column, allow_negative = FALSE) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (!is.numeric(x) && !is.character(x)) {
    stop("`", column, "` must be numeric or numeric-like.", call. = FALSE)
  }
  converted <- suppressWarnings(as.numeric(x))
  if (anyNA(converted) || any(!is.finite(converted))) {
    stop("`", column, "` must contain finite numeric values.", call. = FALSE)
  }
  if (!allow_negative && any(converted < 0)) {
    stop("`", column, "` must not contain negative values.", call. = FALSE)
  }
  converted
}

#' Internal helper: Construct stable group keys for rate calculations.
#'
#' @rdname crtrend-internals
#' @keywords internal
rate_group_key <- function(data, group_vars) {
  if (!length(group_vars)) {
    return(rep(".all", nrow(data)))
  }
  columns <- lapply(
    data[group_vars],
    function(x) factor(as.character(x), levels = unique(as.character(x)))
  )
  as.character(do.call(
    interaction,
    c(columns, list(drop = TRUE, lex.order = FALSE, sep = "\r"))
  ))
}

#' Internal helper: Split rate-calculation rows by group.
#'
#' @rdname crtrend-internals
#' @keywords internal
rate_group_rows <- function(data, group_vars) {
  keys <- rate_group_key(data, group_vars)
  group_id <- match(keys, unique(keys))
  split(seq_len(nrow(data)), group_id)
}

#' Internal helper: Create the common rate-calculation output frame.
#'
#' @rdname crtrend-internals
#' @keywords internal
rate_result_frame <- function(data, group_vars, groups) {
  first_rows <- vapply(groups, function(i) i[[1L]], integer(1))
  if (length(group_vars)) {
    data[first_rows, group_vars, drop = FALSE]
  } else {
    data.frame(row.names = seq_along(groups))
  }
}
