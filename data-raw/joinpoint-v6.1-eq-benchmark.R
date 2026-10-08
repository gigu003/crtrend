# Compare crtrend empirical-quantile inference with the audited Joinpoint 6.1
# exports. Run with Rscript from the package root; no Windows runtime is needed.

devtools::load_all(quiet = TRUE)

audit_dir <- file.path(
  "inst",
  "extdata",
  "joinpoint-6.1.0-eq-audit"
)
reference <- utils::read.csv(
  file.path(audit_dir, "nci-6.1-eq-reference.csv"),
  na.strings = "NA"
)
input_files <- list.files(
  audit_dir,
  pattern = "^[0-9].*\\.txt$",
  full.names = TRUE
)
input_files <- input_files[
  sub("\\.txt$", "", basename(input_files)) %in% unique(reference$case)
]

benchmark_rows <- lapply(input_files, function(input_file) {
  data <- utils::read.delim(input_file, check.names = FALSE)
  case <- sub("\\.txt$", "", basename(input_file))
  fit <- joinpoint(
    data,
    time = "Year",
    outcome = "Rate",
    se = "Standard Error",
    min_jp = 0,
    max_jp = 1,
    selection_method = "WBIC",
    n_eq = 5001,
    eq_seed = 10000
  )
  apc <- calc_apc(fit, ci_method = "EQ")
  aapc <- calc_aapc(
    fit,
    start = min(data$Year),
    end = max(data$Year),
    ci_method = "EQ"
  )
  selected_model <- fit$models[[1L]]
  tau_interval <- if (selected_model$n_joinpoints == 0L) {
    c(NA_real_, NA_real_)
  } else {
    tau <- selected_model$resamples$draws$joinpoints[, 1L]
    stats::quantile(
      tau,
      probs = c(0.025, 0.975),
      names = FALSE,
      type = 7L
    )
  }
  rbind(
    data.frame(
      case = case,
      statistic = "APC",
      segment = apc$segment,
      start = apc$start,
      end = apc$end,
      estimate = apc$apc,
      lower = apc$apc_lower,
      upper = apc$apc_upper,
      p_value = apc$p_value,
      tau_lower = tau_interval[1L],
      tau_upper = tau_interval[2L]
    ),
    data.frame(
      case = case,
      statistic = "AAPC",
      segment = NA_integer_,
      start = aapc$start,
      end = aapc$end,
      estimate = aapc$aapc,
      lower = aapc$aapc_lower,
      upper = aapc$aapc_upper,
      p_value = aapc$p_value,
      tau_lower = tau_interval[1L],
      tau_upper = tau_interval[2L]
    )
  )
})

benchmark <- do.call(rbind, benchmark_rows)
comparison <- merge(
  reference,
  benchmark,
  by = c("case", "statistic", "segment", "start", "end"),
  suffixes = c("_nci", "_crtrend"),
  all = TRUE
)
for (name in c(
  "estimate",
  "lower",
  "upper",
  "p_value",
  "tau_lower",
  "tau_upper"
)) {
  comparison[[paste0(name, "_difference")]] <-
    comparison[[paste0(name, "_crtrend")]] -
    comparison[[paste0(name, "_nci")]]
}
comparison <- comparison[
  order(
    comparison$case,
    comparison$statistic,
    comparison$segment
  ),
]
rownames(comparison) <- NULL
print(comparison, digits = 7)
for (name in c(
  "estimate",
  "lower",
  "upper",
  "p_value",
  "tau_lower",
  "tau_upper"
)) {
  digits <- if (name == "p_value") 6L else 4L
  comparison[[paste0(name, "_matches")]] <-
    (is.na(comparison[[paste0(name, "_nci")]]) &
      is.na(comparison[[paste0(name, "_crtrend")]])) |
    round(comparison[[paste0(name, "_crtrend")]], digits) ==
      comparison[[paste0(name, "_nci")]]
}
utils::write.csv(
  comparison,
  file.path(audit_dir, "crtrend-6.1-eq-comparison.csv"),
  row.names = FALSE
)
matches <- comparison[, endsWith(names(comparison), "_matches")]
stopifnot(nrow(comparison) == nrow(reference), all(as.matrix(matches)))
message(
  "All ",
  nrow(comparison),
  " APC/AAPC rows match the saved export precision."
)
