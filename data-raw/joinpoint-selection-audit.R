# Local diagnostics only: these outputs are NOT NCI reference values.
devtools::load_all(quiet = TRUE)
destination <- "inst/extdata/joinpoint-6.1.0-selection-audit"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
d <- subset(jp_data, sex == 0, c("year", "rate", "se"))
write.table(
  d,
  file.path(destination, "input.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
fit <- joinpoint(d, "year", "rate", se = "se", max_jp = 4, n_eq = 2)
data <- crtrend:::joinpoint_analysis_data(fit$models[[1]])
candidates <- lapply(
  0:4,
  crtrend:::fit_candidate,
  data = data,
  weighted = TRUE,
  min_observations = 2L,
  selection_method = "bic"
)
n <- nrow(data)
largest <- candidates[[5]]
fixed_variance <- largest$sse / stats::df.residual(largest$model)
model_variance <- largest$sse / (n - 2 * largest$n_joinpoints - 2)
rows <- lapply(candidates, function(candidate) {
  k <- candidate$n_joinpoints
  x <- crtrend:::joinpoint_design(candidate$joinpoints, data)
  raw_mbic <- candidate$bic +
    as.numeric(determinant(crossprod(x), logarithm = TRUE)$modulus) / n -
    2 * lgamma((n - k - 3) / 2) / n -
    (k + 3) * log(candidate$sse) / n
  raw_mbic <- raw_mbic - 2 * log(n) / n
  data.frame(
    k,
    sse = candidate$sse,
    bic = candidate$bic,
    bic3 = candidate$bic3,
    mbic_current_whitened = candidate$modified_bic,
    mbic_raw_design_alternative = raw_mbic,
    dds_sigma2_legacy_fixed_grid = fixed_variance,
    dds_sigma2_current = model_variance,
    dds_delta_legacy_fixed_grid = crtrend:::dds_minimum_delta(
      candidate,
      data,
      TRUE,
      fixed_variance
    ),
    dds_delta_current = crtrend:::dds_minimum_delta(
      candidate,
      data,
      TRUE,
      model_variance
    )
  )
})
write.csv(
  do.call(rbind, rows),
  file.path(destination, "diagnostic-alternatives.csv"),
  row.names = FALSE
)
permutation <- joinpoint(
  d,
  "year",
  "rate",
  se = "se",
  max_jp = 4,
  selection_method = "Perm",
  perm_method = "fixed",
  n_perm = 4499,
  seed = 7160,
  n_eq = 2
)
write.csv(
  permutation$selection_tests,
  file.path(destination, "r-permutation-tests.csv"),
  row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(destination, "R-session.txt"))
