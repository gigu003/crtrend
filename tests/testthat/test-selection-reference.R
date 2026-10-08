test_that("DDS and weighted MBIC match the saved 6.1 exports", {
  ref <- read.csv(test_path("fixtures", "joinpoint-6.1-selection.csv"))
  d <- subset(jp_data, sex == 0)
  dds <- joinpoint(
    d,
    "year",
    "rate",
    se = "se",
    max_jp = 4,
    selection_method = "DDS",
    n_eq = 2
  )
  mbic <- joinpoint(
    d,
    "year",
    "rate",
    se = "se",
    max_jp = 4,
    selection_method = "MBIC",
    n_eq = 2
  )
  expect_equal(dds$model_info$n_joinpoints, 1L)
  expect_equal(round(dds$model_info$dds_delta, 7), ref$dds_delta[2])
  expect_equal(mbic$model_info$n_joinpoints, 0L)
  candidates <- calc_model(mbic, ci_method = "Param")$candidates
  expect_equal(candidates$n_joinpoints, ref$k)
  expect_equal(round(candidates$modified_bic, 7), ref$mbic)
})

test_that("MBIC normalization is a constant shift independent of k", {
  d <- data.frame(trend_time = 1:20, trend_log = sin(1:20) + (1:20) / 10)
  for (knots in list(numeric(), 10, c(6, 14))) {
    score <- fixed_joinpoint_scores(knots, d, FALSE)
    k <- length(knots)
    x <- joinpoint_design(knots, d)
    published <- score$bic +
      as.numeric(determinant(crossprod(x), logarithm = TRUE)$modulus) / 20 -
      2 * lgamma((20 - k - 3) / 2) / 20 -
      (k + 3) * log(score$sse) / 20
    expect_equal(score$modified_bic, published - 2 * log(20) / 20)
  }
})

test_that("fixed seeded permutations reproduce the saved 6.1 selection path", {
  ref <- read.csv(test_path("fixtures", "joinpoint-6.1-permutation.csv"))
  fit <- joinpoint(
    subset(jp_data, sex == 0),
    "year",
    "rate",
    se = "se",
    max_jp = 4,
    selection_method = "Perm",
    perm_method = "fixed",
    seed = 7160,
    n_perm = 4499,
    n_eq = 2
  )
  actual <- fit$selection_tests
  expect_equal(actual$k_null, ref$k_null)
  expect_equal(actual$k_alt, ref$k_alt)
  expect_equal(round(actual$p_value, 7), ref$p_value)
  expect_equal(round(actual$alpha, 7), ref$alpha)
  expect_equal(actual$n_permutations, rep(4499L, 4))
  expect_equal(actual$reject, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(fit$model_info$n_joinpoints, 1L)
})

test_that("native MT draws and cumulative shuffles are batch invariant", {
  set.seed(41)
  before <- .Random.seed
  stream <- nci_mt_stream(7160)
  actual <- c(stream(1), stream(623), stream(0), stream(627))
  expect_identical(actual, nci_mt_uniforms(1251, 7160))
  one <- nci_permutation_shuffler(1:22, nci_mt_stream(7160))
  batched <- nci_permutation_shuffler(1:22, nci_mt_stream(7160))
  expected <- one(65)
  expect_identical(cbind(batched(1), batched(31), batched(33)), expected)
  expect_equal(apply(expected, 2, sort), matrix(rep(1:22, 65), 22))
  expect_identical(.Random.seed, before)
})
