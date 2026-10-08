test_that("seeded fitting restores RNG and unseeded work consumes it", {
  d <- subset(jp_data, sex == 0)
  set.seed(91)
  before <- .Random.seed
  for (method in c("fixed", "tsprt")) {
    fit <- joinpoint(
      d,
      "year",
      "rate",
      max_jp = 1,
      n_eq = 3,
      selection_method = "Perm",
      perm_method = method,
      n_perm = 19
    )
    expect_identical(.Random.seed, before)
  }
  joinpoint(d, "year", "rate", max_jp = 0, n_eq = 3, eq_seed = NULL)
  expect_identical(.Random.seed, before)
  joinpoint(d, "year", "rate", max_jp = 0, n_eq = 3, seed = NULL)
  expect_identical(.Random.seed, before)
  joinpoint(
    d,
    "year",
    "rate",
    max_jp = 0,
    n_eq = 3,
    seed = NULL,
    eq_seed = NULL
  )
  expect_identical(identical(.Random.seed, before), FALSE)
  set.seed(91)
  joinpoint(
    d,
    "year",
    "rate",
    max_jp = 1,
    n_eq = 3,
    seed = NULL,
    selection_method = "Perm",
    n_perm = 19
  )
  expect_identical(identical(.Random.seed, before), FALSE)
})

test_that("local seed restores absent RNG and errors", {
  set.seed(19)
  old <- .Random.seed
  on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  expect_snapshot(
    error = TRUE,
    with_joinpoint_seed(41, stop("deliberate failure"))
  )
  expect_identical(.Random.seed, old)
  rm(list = ".Random.seed", envir = .GlobalEnv)
  d <- subset(jp_data, sex == 0)
  joinpoint(d, "year", "rate", max_jp = 0, n_eq = 3)
  expect_identical(
    exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    FALSE
  )
  expect_snapshot(
    error = TRUE,
    with_joinpoint_seed(41, stop("absent RNG failure"))
  )
  expect_identical(
    exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    FALSE
  )
})

test_that("numeric controls fail before unsafe coercion", {
  d <- subset(jp_data, sex == 0)
  expect_snapshot({
    for (name in c(
      "max_jp",
      "min_jp",
      "min_obs",
      "n_perm",
      "n_eq",
      "seed",
      "eq_seed"
    )) {
      for (value in list(Inf, NaN, 2^32, 2.5, numeric(), c(1, 2), 1i)) {
        args <- list(
          data = d,
          time = "year",
          outcome = "rate",
          max_jp = 0,
          n_eq = 2
        )
        args[[name]] <- value
        cat(
          name,
          ": ",
          tryCatch(
            {
              do.call(joinpoint, args)
              "ACCEPTED"
            },
            error = function(e) conditionMessage(e)
          ),
          "\n",
          sep = ""
        )
      }
    }
  })
  expect_identical(is_bounded_integer(.Machine$integer.max - 1), TRUE)
  expect_identical(is_bounded_integer(.Machine$integer.max), FALSE)
  expect_null(validate_rng_seed(-.Machine$integer.max))
  expect_null(validate_rng_seed(.Machine$integer.max))
  expect_null(integer_joinpoint_grid(2000:2021, 2L, .Machine$integer.max - 1L))
  fit <- joinpoint(d, "year", "rate", max_jp = 0, n_eq = 3)
  expect_snapshot(error = TRUE, calc_apc(fit, n_resamples = Inf))
  expect_snapshot(error = TRUE, calc_aapc(fit, seed = 1.5))
  expect_snapshot(error = TRUE, calc_aapc(fit, start = -2^32, end = 2021))
})

test_that("fitting errors restore caller RNG", {
  d <- subset(jp_data, sex == 0)
  set.seed(55)
  before <- .Random.seed
  local_mocked_bindings(fit_joinpoint_groups = function(...) {
    runif(3)
    stop("injected fitting error", call. = FALSE)
  })
  expect_snapshot(error = TRUE, joinpoint(d, "year", "rate", n_eq = 2))
  expect_identical(.Random.seed, before)
})

test_that("real controls reject non-finite and non-scalar input", {
  d <- subset(jp_data, sex == 0)
  expect_snapshot({
    for (name in c(
      "rho",
      "alpha",
      "perm_alpha0",
      "perm_beta0",
      "perm_conf",
      "perm_p1"
    )) {
      args <- list(data = d, time = "year", outcome = "rate", n_eq = 2)
      args[[name]] <- Inf
      cat(
        name,
        ": ",
        tryCatch(do.call(joinpoint, args), error = function(e) {
          conditionMessage(e)
        }),
        "\n",
        sep = ""
      )
    }
  })
})

test_that("seeded AAPC plot labels preserve RNG", {
  skip_if_not_installed("ggplot2")
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    max_jp = 0,
    n_eq = 3
  )
  set.seed(19)
  old <- .Random.seed
  plot_joinpoint(fit, show_aapc = TRUE, seed = 123, n_resamples = 5)
  expect_identical(.Random.seed, old)
})
