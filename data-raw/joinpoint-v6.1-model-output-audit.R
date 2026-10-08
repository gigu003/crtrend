# Run from the package root. Screenshot reference: jp_data, sex = 0,
# Joinpoint 6.1.0, WBIC, 0-4 joinpoints, independent errors, Parametric.
# This audit reconstructs currently unexposed outputs; it does not change APIs.
devtools::load_all(quiet = TRUE)
fit <- joinpoint(
  subset(jp_data, sex == 0),
  "year",
  "rate",
  se = "se",
  max_jp = 4,
  n_eq = 5001,
  eq_seed = 10000
)
model <- fit$models[[1L]]
analysis <- model$.analysis
out_dir <- "inst/extdata/joinpoint-6.1.0-model-audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
comparisons <- list()
compare <- function(section, field, actual, reference, digits) {
  stopifnot(length(actual) == length(reference))
  comparisons[[length(comparisons) + 1L]] <<- data.frame(
    section,
    field,
    row = seq_along(actual),
    actual,
    reference,
    digits,
    difference = actual - reference,
    rounded_match = round(actual, digits) == reference,
    # One unit in the final screenshot decimal, not bitwise equality.
    within_display_unit = abs(actual - reference) <= 10^(-digits)
  )
}
observed <- c(
  55.71,
  54.86,
  54.06,
  52.48,
  51.11,
  49.54,
  48.12,
  47.38,
  46.09,
  44.10,
  42.35,
  41.27,
  40.04,
  39.33,
  39.44,
  38.94,
  38.49,
  37.62,
  37.59,
  37.44,
  33.31,
  36.67
)
se <- c(
  .2113,
  .2082,
  .2050,
  .2004,
  .1963,
  .1917,
  .1876,
  .1844,
  .1802,
  .1746,
  .1694,
  .1659,
  .1617,
  .1586,
  .1572,
  .1550,
  .1528,
  .1497,
  .1484,
  .1470,
  .1381,
  .1444
)
modeled <- c(
  56.84,
  55.26,
  53.71,
  52.22,
  50.76,
  49.35,
  47.97,
  46.63,
  45.33,
  44.07,
  42.84,
  41.64,
  40.48,
  39.90,
  39.32,
  38.75,
  38.19,
  37.64,
  37.10,
  36.56,
  36.04,
  35.52
)
compare("data", "observed", model$model$rate, observed, 2)
compare("data", "se", model$model$se, se, 4)
compare("data", "fitted", as.numeric(fitted(model)), modeled, 2)
compare("joinpoints", "estimate", model$joinpoints, 2012, 0)

candidates <- lapply(0:4, function(k) {
  crtrend:::fit_candidate(k, analysis, TRUE, 2L, "weighted_bic")
})
candidate_table <- do.call(
  rbind,
  lapply(candidates, function(x) {
    k <- x$n_joinpoints
    data.frame(
      k,
      n = 22L,
      parameters = 2L * k + 2L,
      df = 22L - (2L * k + 2L),
      sse = x$sse,
      bic = x$bic,
      bic3 = x$bic3,
      weight = if (k == 0L) NA_real_ else x$partial_r_squared,
      wbic = x$weighted_bic,
      joinpoints = paste(x$joinpoints, collapse = ";")
    )
  })
)
reference_scores <- data.frame(
  sse = c(1220.6072632, 587.7385961, 523.1575449, 503.7452244, 488.6752202),
  bic = c(4.2970652, 3.8472475, 4.0118518, 4.2550436, 4.5056749),
  bic3 = c(4.2970652, 3.9877495, 4.2928556, 4.6765494, 5.0676827),
  weight = c(NA, .5184867, .6844505, .9146089, .9116191),
  wbic = c(4.2970652, 3.9200959, 4.2041850, 4.6405566, 5.0180119)
)
for (field in names(reference_scores)) {
  keep <- !is.na(reference_scores[[field]])
  compare(
    "candidates",
    field,
    candidate_table[[field]][keep],
    reference_scores[[field]][keep],
    7
  )
}
selected <- candidate_table[which.min(candidate_table$wbic), ]
compare(
  "statistics",
  "n_parameters_df",
  unlist(selected[c("n", "parameters", "df")]),
  c(22, 4, 18),
  0
)
compare(
  "statistics",
  "sse_mse",
  c(selected$sse, selected$sse / selected$df),
  c(587.7386, 32.6521),
  4
)

# Reproduce the package's NCI inference covariance before its slope-only subset.
# Continuous-grid coefficient estimates are NOT replaced by unconstrained ones.
d <- analysis[!analysis$trend_time %in% model$joinpoints, , drop = FALSE]
segment <- findInterval(d$trend_time, model$joinpoints) + 1L
design <- matrix(0, nrow(d), 4L)
for (j in 1:2) {
  design[segment == j, 2L * j - 1L] <- 1
  design[segment == j, 2L * j] <- d$trend_time[segment == j]
}
xw <- crtrend:::joinpoint_whiten(d, design)
yw <- crtrend:::joinpoint_whiten(d, d$trend_log)
solver <- qr(xw)
df <- nrow(d) - ncol(design)
variance <- sum((yw - xw %*% qr.coef(solver, yw))^2) / df
covariance <- variance * chol2inv(chol(crossprod(xw)))
b <- as.numeric(coef(model))
general <- c(b[1], b[2], b[1] - b[3] * 2012, b[2] + b[3])
standard_transform <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, -1, 0, 1))
coef_table <- function(est, cov, terms) {
  se <- sqrt(diag(cov))
  data.frame(
    term = terms,
    estimate = est,
    se,
    statistic = est / se,
    p_value = 2 * pt(-abs(est / se), df)
  )
}
general_table <- coef_table(
  general,
  covariance,
  c("Intercept 1", "Slope 1", "Intercept 2", "Slope 2")
)
standard_table <- coef_table(
  b,
  standard_transform %*% covariance %*% t(standard_transform),
  c("Intercept 1", "Slope 1", "Slope 2 - Slope 1")
)
reference_general <- data.frame(
  estimate = c(60.600790, -.028280, 32.963267, -.014544),
  se = c(3.802859, .001896, 6.071286, .003010),
  statistic = c(15.935587, -14.913115, 5.429371, -4.831790),
  p_value = c(0, 0, .000045, .000156)
)
reference_standard <- data.frame(
  estimate = c(60.600790, -.028280, .013736),
  se = c(3.802859, .001896, .003558),
  statistic = c(15.935587, -14.913115, 3.861130),
  p_value = c(0, 0, .001253)
)
for (field in names(reference_general)) {
  compare(
    "general_coefficients",
    field,
    general_table[[field]],
    reference_general[[field]],
    6
  )
  compare(
    "standard_coefficients",
    field,
    standard_table[[field]],
    reference_standard[[field]],
    6
  )
}
reference_covariance <- matrix(
  c(
    14.461736,
    -.007211,
    0,
    0,
    -.007211,
    .000004,
    0,
    0,
    0,
    0,
    36.860516,
    -.018275,
    0,
    0,
    -.018275,
    .000009
  ),
  4,
  4
)
reference_correlation <- diag(4)
reference_correlation[1, 2] <- reference_correlation[2, 1] <- -.999999
reference_correlation[3, 4] <- reference_correlation[4, 3] <- -.999999
compare(
  "covariance",
  "value",
  as.vector(covariance),
  as.vector(reference_covariance),
  6
)
compare(
  "correlation",
  "value",
  as.vector(cov2cor(covariance)),
  as.vector(reference_correlation),
  6
)

apc <- calc_apc(fit, ci_method = "Param")
aapc <- calc_aapc(fit, start = 2000, end = 2021, ci_method = "Param")
compare("trends", "apc", apc$apc, c(-2.7884, -1.4439), 4)
compare("trends", "apc_lower", apc$apc_lower, c(-3.1766, -2.0678), 4)
compare("trends", "apc_upper", apc$apc_upper, c(-2.3987, -.8160), 4)
compare("trends", "apc_statistic", apc$statistic, c(-14.9131, -4.8318), 4)
compare("trends", "apc_p_value", apc$p_value[2], .000156, 6)
compare("trends", "aapc", aapc$aapc, -2.2144, 4)
compare("trends", "aapc_lower", aapc$aapc_lower, -2.5368, 4)
compare("trends", "aapc_upper", aapc$aapc_upper, -1.8910, 4)
compare("trends", "aapc_statistic", aapc$statistic, -13.2917, 4)
stopifnot(apc$p_value[1] < 1e-6, aapc$p_value < 1e-6)

# Tau profile is a diagnostic, not a validated general-purpose CI implementation.
grid <- crtrend:::cached_integer_joinpoint_grid(
  analysis$trend_time,
  1L,
  min_points = 2L
)
tau_profile <- data.frame(
  tau = unlist(grid),
  sse = vapply(
    grid,
    function(tau) {
      crtrend:::fixed_joinpoint_scores(tau, analysis, TRUE)$sse
    },
    numeric(1)
  )
)
critical_f <- qf(.95, 1, selected$df)
threshold <- selected$sse * (1 + critical_f / selected$df)
threshold_web_literal <- selected$sse * (1 + 1 / selected$df) * critical_f
tau_profile$inside_profile_ci <- tau_profile$sse <= threshold
tau_profile$inside_literal_web_ci <- tau_profile$sse <= threshold_web_literal
# Investigate outward rounding only; matching this one case is not proof of
# Joinpoint's endpoint algorithm.
profile_difference <- function(tau) {
  crtrend:::fixed_joinpoint_scores(tau, analysis, TRUE)$sse - threshold
}
continuous_bounds <- c(
  uniroot(profile_difference, c(2009, 2010), tol = 1e-9)$root,
  uniroot(profile_difference, c(2014, 2015), tol = 1e-9)$root
)

plot_segments <- crtrend:::joinpoint_plot_segments(fit$models)
years <- model$model$year
plotted <- approx(
  plot_segments$.plot_time,
  plot_segments$.plot_outcome,
  xout = years,
  ties = mean
)$y
endpoint_times <- sort(unique(unlist(
  model$segments[c("start", "end")],
  use.names = FALSE
)))
endpoint_data <- model$model[rep(1L, length(endpoint_times)), , drop = FALSE]
endpoint_data[[model$control$time]] <- endpoint_times
legacy_chord <- approx(
  endpoint_times,
  as.numeric(predict(model, endpoint_data)),
  xout = years
)$y
plot_audit <- data.frame(
  year = years,
  fitted = as.numeric(fitted(model)),
  plotted_curve = plotted,
  difference = plotted - as.numeric(fitted(model)),
  legacy_endpoint_chord = legacy_chord,
  legacy_difference = legacy_chord - as.numeric(fitted(model))
)
stopifnot(max(abs(plot_audit$difference)) < 1e-10)
data_table <- data.frame(
  year = years,
  observed = model$model$rate,
  se = model$model$se,
  fitted = as.numeric(fitted(model)),
  joinpoint = ifelse(years == 2012, 1L, NA_integer_),
  apc = ifelse(
    years == 2012,
    NA_real_,
    ifelse(years < 2012, apc$apc[1], apc$apc[2])
  )
)
keep <- years != 2012
compare(
  "data",
  "yearly_apc",
  data_table$apc[keep],
  c(rep(-2.7884, 12), rep(-1.4439, 9)),
  4
)
stopifnot(all(apc$p_value < .05), is.na(data_table$apc[years == 2012]))

# A time-origin shift is a numerical stability diagnostic, not an NCI oracle.
centered_analysis <- analysis
centered_analysis$trend_time <- centered_analysis$trend_time - 2000
origin_audit <- do.call(
  rbind,
  lapply(0:4, function(k) {
    centered <- crtrend:::fit_candidate(
      k,
      centered_analysis,
      TRUE,
      2L,
      "weighted_bic"
    )
    original <- candidates[[k + 1L]]
    data.frame(
      k,
      sse_difference = centered$sse - original$sse,
      weight_difference = centered$partial_r_squared -
        original$partial_r_squared
    )
  })
)
result <- do.call(rbind, comparisons)
tables <- list(
  comparison = result,
  candidates = candidate_table,
  data = data_table,
  standard_coefficients = standard_table,
  general_coefficients = general_table,
  covariance = covariance,
  correlation = cov2cor(covariance),
  tau_profile = tau_profile,
  plot_audit = plot_audit,
  time_origin_audit = origin_audit
)
for (name in names(tables)) {
  write.csv(
    tables[[name]],
    file.path(out_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
print(aggregate(
  cbind(rounded_match, within_display_unit) ~ section,
  result,
  sum
))
print(result[!result$rounded_match, ])
print(candidate_table, digits = 12)
print(general_table, digits = 12)
print(standard_table, digits = 12)
print(list(
  model_df = selected$df,
  fixed_grid_df = model$df.residual,
  inference_df = df,
  inference_residual_variance = variance,
  tau_profile_threshold = threshold,
  tau_profile_ci = range(tau_profile$tau[tau_profile$inside_profile_ci]),
  literal_web_threshold = threshold_web_literal,
  literal_web_ci = range(tau_profile$tau[tau_profile$inside_literal_web_ci]),
  continuous_profile_bounds = continuous_bounds,
  exploratory_outward_rounded_bounds = c(
    floor(continuous_bounds[1]),
    ceiling(continuous_bounds[2])
  ),
  max_plot_difference = max(abs(plot_audit$difference))
))
# Candidate mismatches are findings, not waived equivalence tests. Keep their
# exact discrepancies in comparison.csv, and assert only established matches.
stopifnot(all(result$rounded_match[result$section != "candidates"]))
stopifnot(selected$k == 1L, selected$joinpoints == "2012")
