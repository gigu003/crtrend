# R-side scheduling checks and inputs for native validation; not NCI results.
devtools::load_all(quiet = TRUE)
destination <- "inst/extdata/joinpoint-6.1.0-permutation-groups"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
d <- jp_data[, c("sex", "year", "rate", "se")]
d <- d[order(d$sex, d$year), ]
d$cohort <- d$sex
write_input <- function(x, name) {
  write.table(
    x[, c("cohort", "year", "rate", "se")],
    file.path(destination, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
write_input(d, "groups-original.txt")
reversed <- d
reversed$cohort <- 2 - reversed$sex
reversed <- reversed[order(reversed$cohort, reversed$year), ]
write_input(reversed, "groups-relabelled.txt")
for (sex in 0:2) {
  write_input(d[d$sex == sex, ], paste0("single-", sex, ".txt"))
}
fit <- function(x, grouped = TRUE) {
  joinpoint(
    x,
    "year",
    "rate",
    se = "se",
    group_vars = if (grouped) "cohort" else NULL,
    selection_method = "Perm",
    perm_method = "fixed",
    seed = 7160,
    n_perm = 4499,
    min_jp = 0,
    max_jp = 4,
    min_obs = 2,
    n_eq = 2
  )
}
original <- fit(d)
relabelled <- fit(reversed)
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
  single <- fit(d[d$sex == sex, ], FALSE)
  a <- original$selection_tests[original$selection_tests$cohort == sex, fields]
  b <- relabelled$selection_tests[
    relabelled$selection_tests$cohort == 2 - sex,
    fields
  ]
  c <- single$selection_tests[, fields]
  rownames(a) <- rownames(b) <- rownames(c) <- NULL
  stopifnot(identical(a, b), identical(a, c))
  write.csv(
    single$selection_tests,
    file.path(destination, paste0("r-single-", sex, ".csv")),
    row.names = FALSE
  )
}
write.csv(
  original$selection_tests,
  file.path(destination, "r-groups-original.csv"),
  row.names = FALSE
)
write.csv(
  relabelled$selection_tests,
  file.path(destination, "r-groups-relabelled.csv"),
  row.names = FALSE
)
write.csv(
  original$model_info,
  file.path(destination, "r-model-info.csv"),
  row.names = FALSE
)
message(
  "R-only audit passed: all three cohorts agree across grouped, single and relabelled runs."
)
