test_that("short rate names and compatibility aliases agree", {
  expect_equal(names(formals(calc_asr))[5L], "std")
  expect_equal("standard_population" %in% names(formals(calc_asr)), FALSE)
  expect_equal(calc_cr(jp_crude), calc_crude_rate(jp_crude))
  expect_equal(
    calc_asr(jp_age_adjusted, age = "age_group", std = "standard_population"),
    calc_age_adjusted_rate(
      jp_age_adjusted,
      age = "age_group",
      std = "standard_population"
    )
  )
})

test_that("proportions and percentages use grouped binomial standard errors", {
  data <- data.frame(
    year = c(2000, 2000, 2001, 2001),
    stratum = c("a", "b", "a", "b"),
    count = c(20, 30, 24, 36),
    population = c(50, 50, 60, 60)
  )
  prop <- calc_prop(data, group_vars = "year")
  percent <- calc_percent(data, group_vars = "year")
  expected_prop <- c(0.5, 0.5)
  expected_se <- sqrt(expected_prop * (1 - expected_prop) / c(100, 120))

  expect_named(prop, c("year", "count", "population", "prop", "se"))
  expect_named(
    percent,
    c("year", "count", "population", "percent", "se")
  )
  expect_equal(prop$count, c(50, 60))
  expect_equal(prop$population, c(100, 120))
  expect_equal(prop$prop, expected_prop)
  expect_equal(prop$se, expected_se)
  expect_equal(percent$percent, 100 * expected_prop)
  expect_equal(percent$se, 100 * expected_se)
})

test_that("proportion inputs and normal approximation are checked", {
  expect_snapshot(
    error = TRUE,
    calc_prop(data.frame(count = 2, population = 1), group_vars = character())
  )
  expect_snapshot(
    calc_prop(data.frame(count = 1, population = 20), group_vars = character())
  )
})
