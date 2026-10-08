test_that("Joinpoint sample datasets preserve their expected schemas", {
  expect_named(
    jp_crude,
    c("sex", "year", "count", "population")
  )
  expect_named(
    jp_age_adjusted,
    c(
      "sex",
      "year",
      "age_group",
      "count",
      "population",
      "standard_population"
    )
  )
  expect_named(
    jp_data,
    c(
      "sex",
      "year_index",
      "year",
      "rate",
      "se",
      "lower_ci",
      "upper_ci",
      "count",
      "population"
    )
  )
  expect_named(
    jp_jump,
    c(
      "year_index",
      "year",
      "rate",
      "se",
      "lower_ci",
      "upper_ci",
      "count",
      "population"
    )
  )

  expect_equal(nrow(jp_crude), 66L)
  expect_equal(nrow(jp_age_adjusted), 1311L)
  expect_equal(nrow(jp_data), 66L)
  expect_equal(nrow(jp_jump), 31L)
})

test_that("Joinpoint rate samples can be fitted directly", {
  seer_fit <- joinpoint(
    jp_data,
    time = "year",
    outcome = "rate",
    group_vars = "sex",
    se = "se",
    max_jp = 0,
    n_eq = 19L
  )
  jump_fit <- joinpoint(
    jp_jump,
    time = "year",
    outcome = "rate",
    se = "se",
    max_jp = 0,
    n_eq = 19L
  )

  expect_s3_class(seer_fit, "crtrend_joinpoint")
  expect_s3_class(jump_fit, "crtrend_joinpoint")
  expect_equal(nrow(seer_fit$model_info), 3L)
  expect_equal(nrow(jump_fit$model_info), 1L)
  expect_equal(seer_fit$control$se, "se")
  expect_equal(jump_fit$control$se, "se")
})

test_that("crude rates use aggregated counts and populations", {
  crude <- calc_crude_rate(jp_crude)

  expect_named(crude, c("sex", "year", "count", "population", "rate", "se"))
  expect_equal(nrow(crude), nrow(jp_crude))
  expect_equal(crude$count, as.numeric(jp_crude$count))
  expect_equal(crude$population, as.numeric(jp_crude$population))
  expect_equal(
    crude$rate,
    crude$count / crude$population * 100000,
    tolerance = 1e-12
  )
  expect_equal(
    crude$se,
    sqrt(crude$count) / crude$population * 100000,
    tolerance = 1e-12
  )

  fit <- joinpoint(
    crude,
    time = "year",
    outcome = "rate",
    group_vars = "sex",
    se = "se",
    max_jp = 0,
    n_eq = 19L
  )
  expect_s3_class(fit, "crtrend_joinpoint")
  expect_equal(nrow(fit$model_info), 3L)
})

test_that("age-adjusted rates use direct standardization", {
  age_adjusted <- calc_age_adjusted_rate(
    jp_age_adjusted,
    age = "age_group"
  )
  first <- subset(
    jp_age_adjusted,
    sex == jp_age_adjusted$sex[1L] & year == jp_age_adjusted$year[1L]
  )
  weights <- first$standard_population / sum(first$standard_population)
  expected_rate <- sum(weights * first$count / first$population) * 100000
  expected_se <- sqrt(
    sum(weights^2 * first$count / first$population^2)
  ) *
    100000

  expect_named(
    age_adjusted,
    c("sex", "year", "count", "population", "rate", "se")
  )
  expect_equal(nrow(age_adjusted), 3L * 23L)
  expect_equal(age_adjusted$rate[1L], expected_rate, tolerance = 1e-12)
  expect_equal(age_adjusted$se[1L], expected_se, tolerance = 1e-12)
  expect_equal(age_adjusted$count[1L], sum(first$count))
  expect_equal(age_adjusted$population[1L], sum(first$population))

  fit <- joinpoint(
    age_adjusted,
    time = "year",
    outcome = "rate",
    group_vars = "sex",
    se = "se",
    max_jp = 0,
    n_eq = 19L
  )
  expect_s3_class(fit, "crtrend_joinpoint")
  expect_equal(nrow(fit$model_info), 3L)
})

test_that("age-adjusted rates and Joinpoint fits agree with V5.3", {
  reference <- new.env(parent = emptyenv())
  load(
    test_path("fixtures", "joinpoint-v5.3-fitted.rda"),
    envir = reference
  )
  expected <- reference$v5_3_age_adjusted
  age_adjusted <- calc_age_adjusted_rate(
    jp_age_adjusted,
    age = "age_group"
  )
  expected_rates <- expected[expected$model == 0L, ]
  expected_rates <- expected_rates[
    order(expected_rates$sex, expected_rates$year),
  ]
  observed_rates <- age_adjusted[order(age_adjusted$sex, age_adjusted$year), ]
  row.names(expected_rates) <- NULL
  row.names(observed_rates) <- NULL

  expect_equal(
    observed_rates[, c("sex", "year", "rate", "se")],
    expected_rates[, c("sex", "year", "rate", "se")],
    tolerance = 1e-8
  )

  for (n_joinpoints in 0:4) {
    fit <- joinpoint(
      age_adjusted,
      time = "year",
      outcome = "rate",
      group_vars = "sex",
      se = "se",
      min_jp = n_joinpoints,
      max_jp = n_joinpoints,
      n_eq = 19L
    )
    observed <- data.frame(
      model = n_joinpoints,
      sex = age_adjusted$sex,
      year = age_adjusted$year,
      fitted = as.numeric(fitted(fit))
    )
    expected_model <- expected[expected$model == n_joinpoints, ]
    expected_model <- expected_model[
      order(expected_model$sex, expected_model$year),
    ]
    observed <- observed[order(observed$sex, observed$year), ]
    row.names(expected_model) <- NULL
    row.names(observed) <- NULL
    expect_equal(
      observed,
      expected_model[, c("model", "sex", "year", "fitted")],
      tolerance = 1e-8
    )
  }

  selected_fit <- joinpoint(
    age_adjusted,
    time = "year",
    outcome = "rate",
    group_vars = "sex",
    se = "se",
    max_jp = 4,
    n_eq = 19L
  )
  expect_equal(
    selected_fit$joinpoints[, c("sex", "joinpoint")],
    data.frame(
      sex = c(0L, 1L, 2L),
      joinpoint = c(2012, 2012, 2012)
    )
  )
})

test_that("parametric APC and AAPC agree with the Joinpoint 6.1 baseline", {
  fit <- joinpoint(
    jp_data,
    time = "year",
    outcome = "rate",
    group_vars = "sex",
    se = "se",
    min_jp = 1,
    max_jp = 1,
    n_eq = 19L
  )
  apc <- calc_apc(fit, ci_method = "Param")
  aapc <- calc_aapc(
    fit,
    start = 2000,
    end = 2021,
    ci_method = "Param"
  )

  expect_equal(
    apc[apc$sex == 0, c("apc", "apc_lower", "apc_upper")],
    data.frame(
      apc = c(-2.788413, -1.443868),
      apc_lower = c(-3.176572, -2.067780),
      apc_upper = c(-2.398698, -0.815981)
    ),
    tolerance = 1e-6
  )
  expect_lt(
    max(
      abs(
        unname(apc[apc$sex == 0, "p_value"]) -
          c(3.395006e-11, 1.560300e-4)
      )
    ),
    1e-10
  )
  expect_equal(
    aapc[aapc$sex == 0, c("aapc", "aapc_lower", "aapc_upper")],
    data.frame(
      aapc = -2.214440,
      aapc_lower = -2.536801,
      aapc_upper = -1.891013
    ),
    tolerance = 1e-6
  )
  expect_lt(
    abs(unname(aapc[aapc$sex == 0, "p_value"]) - 2.585288e-40),
    1e-44
  )
})

test_that("Joinpoint 6.1 EQ audit inputs retain their diagnostic designs", {
  audit_dir <- system.file(
    "extdata",
    "joinpoint-6.1.0-eq-audit",
    package = "crtrend"
  )
  filenames <- c(
    "01-equal-log-se-high-leverage.txt",
    "02-heteroscedastic-whitening.txt",
    "03-endpoint-outlier.txt"
  )
  inputs <- lapply(file.path(audit_dir, filenames), read.delim)

  expect_length(inputs, 3L)
  expect_identical(vapply(inputs, nrow, integer(1L)), c(12L, 14L, 10L))
  expect_identical(
    lapply(inputs, names),
    rep(list(c("Year", "Rate", "Standard.Error")), 3L)
  )
  expect_lt(
    max(abs(inputs[[1L]]$Standard.Error / inputs[[1L]]$Rate - 0.02)),
    1e-10
  )
  expect_gt(
    diff(range(inputs[[2L]]$Standard.Error / inputs[[2L]]$Rate)),
    0.05
  )
  expect_lt(
    max(abs(inputs[[3L]]$Standard.Error / inputs[[3L]]$Rate - 0.018)),
    1e-10
  )

  fits <- lapply(inputs, function(data) {
    joinpoint(
      data,
      time = "Year",
      outcome = "Rate",
      se = "Standard.Error",
      min_jp = 1L,
      max_jp = 1L,
      n_eq = 19L
    )
  })
  expect_identical(
    vapply(fits, \(fit) fit$joinpoints$joinpoint, numeric(1L)),
    c(2, 5, 2)
  )
})

test_that("max_jp changes the candidate range without changing the selected trend", {
  age_adjusted <- calc_age_adjusted_rate(
    jp_age_adjusted,
    age = "age_group"
  )
  max_values <- c(0L, 1L, 4L, 8L)
  fits <- lapply(max_values, function(max_jp) {
    joinpoint(
      age_adjusted,
      time = "year",
      outcome = "rate",
      group_vars = "sex",
      se = "se",
      max_jp = max_jp,
      n_eq = 19L
    )
  })

  expect_equal(
    lapply(fits, function(x) x$model_info$n_joinpoints),
    list(c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 1L, 1L), c(1L, 1L, 1L))
  )
  expect_equal(
    lapply(fits, function(x) unique(x$model_info$effective_max_jp)),
    list(0L, 1L, 4L, 7L)
  )
  expect_equal(
    lapply(fits[-1L], function(x) x$joinpoints$joinpoint),
    list(c(2012, 2012, 2012), c(2012, 2012, 2012), c(2012, 2012, 2012))
  )
})

test_that("rate calculations reject invalid or duplicated age inputs", {
  expect_error(
    calc_crude_rate(
      data.frame(count = 1, population = 0),
      group_vars = character()
    ),
    "positive"
  )
  expect_error(
    calc_age_adjusted_rate(
      data.frame(
        sex = 1,
        age = c(1, 1),
        count = c(1, 2),
        population = c(10, 10),
        standard_population = c(1, 1)
      ),
      group_vars = "sex"
    ),
    "at most one observation per age stratum"
  )
})
