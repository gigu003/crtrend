test_that("seeded fixed selection resets per group independent of group order", {
  reference <- read.csv(test_path(
    "fixtures",
    "joinpoint-6.1-permutation-groups.csv"
  ))
  d <- jp_data[, c("sex", "year", "rate", "se")]
  run <- function(x, group = "sex") {
    joinpoint(
      x,
      "year",
      "rate",
      se = "se",
      group_vars = group,
      selection_method = "Perm",
      perm_method = "fixed",
      seed = 7160,
      n_perm = 4499,
      max_jp = 4,
      n_eq = 2
    )
  }
  grouped <- run(d)
  reverse <- d[nrow(d):1, ]
  reverse$sex <- 2 - reverse$sex
  reordered <- run(reverse)
  fields <- c(
    "step",
    "k_null",
    "k_alt",
    "alpha",
    "statistic",
    "p_value",
    "n_permutations",
    "reject"
  )
  for (sex in 0:2) {
    single <- run(d[d$sex == sex, ], NULL)
    a <- grouped$selection_tests[grouped$selection_tests$sex == sex, fields]
    b <- reordered$selection_tests[
      reordered$selection_tests$sex == 2 - sex,
      fields
    ]
    rownames(a) <- rownames(b) <- NULL
    expect_identical(a, single$selection_tests[, fields])
    expect_identical(a, b)
    for (run_label in c("Original", "relabelled", paste0("single-", sex))) {
      ref <- reference[reference$sex == sex & reference$run == run_label, ]
      expect_equal(nrow(ref), 4L)
      expect_equal(round(a$p_value, 7), ref$p_value)
      expect_equal(round(a$alpha, 7), ref$alpha)
      expect_equal(a$k_null, ref$k_null)
      expect_equal(a$k_alt, ref$k_alt)
      expect_equal(ifelse(a$reject, a$k_alt, a$k_null), ref$selected)
      expect_equal(a$n_permutations + 1L, ref$n_total)
      expect_equal(rep(single$model_info$n_joinpoints, nrow(ref)), ref$final_k)
    }
    expect_equal(
      grouped$model_info$n_joinpoints[grouped$model_info$sex == sex],
      single$model_info$n_joinpoints
    )
  }
})
