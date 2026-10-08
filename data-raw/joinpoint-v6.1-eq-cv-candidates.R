# Static candidate audit for the Joinpoint 6.1.0 capcvErrors vector.

devtools::load_all(quiet = TRUE)
audit_dir <- file.path("inst", "extdata", "joinpoint-6.1.0-eq-audit")
reference <- utils::read.csv(
  file.path(audit_dir, "nci-6.1-eq-reference.csv"),
  na.strings = "NA"
)

draw_summary <- function(
  fit,
  residuals,
  scale,
  correction,
  uniform
) {
  data <- crtrend:::joinpoint_analysis_data(fit)
  design <- crtrend:::joinpoint_design(numeric(), data)
  whitened_design <- crtrend:::joinpoint_whiten(data, design)
  leverage <- stats::hat(whitened_design)
  residuals <- switch(
    correction,
    none = residuals,
    hc2 = residuals / sqrt(1 - leverage),
    hc3 = residuals / (1 - leverage),
    press = residuals / (1 - leverage)
  )
  endpoints <- crtrend:::nci_empirical_endpoints(residuals)
  knots <- sort(c(endpoints, sort(residuals)))
  n <- length(residuals)
  interval <- pmin(floor(uniform * (n + 1L)), n)
  fraction <- uniform * (n + 1L) - interval
  sampled <- matrix(
    knots[interval + 1L] +
      (knots[interval + 2L] - knots[interval + 1L]) * fraction,
    nrow = nrow(data),
    ncol = 5001L
  )
  if (scale == "whitened") {
    sampled <- crtrend:::joinpoint_unwhiten(data, sampled)
  }
  responses <- matrix(
    stats::fitted(crtrend:::joinpoint_engine(fit)),
    nrow = nrow(data),
    ncol = 5001L
  ) +
    sampled
  coefficients <- qr.coef(
    qr(whitened_design),
    crtrend:::joinpoint_whiten(data, responses)
  )
  draws <- 100 * expm1(coefficients[2L, ])
  limits <- stats::quantile(
    draws,
    c(0.025, 0.975),
    names = FALSE,
    type = 7L
  )
  c(
    lower = limits[1L],
    upper = limits[2L],
    p_value = 2 * min(sum(draws < 0), sum(draws >= 0)) / length(draws)
  )
}

cases <- sprintf("%02d", 7:14)
rows <- lapply(cases, function(number) {
  input <- list.files(
    audit_dir,
    pattern = paste0("^", number, ".*\\.txt$"),
    full.names = TRUE
  )
  data <- utils::read.delim(input, check.names = FALSE)
  fit <- joinpoint(
    data,
    time = "Year",
    outcome = "Rate",
    se = "Standard Error",
    min_jp = 0,
    max_jp = 0,
    n_eq = 2,
    eq_seed = NULL
  )$models[[1L]]
  analysis <- crtrend:::joinpoint_analysis_data(fit)
  residuals <- as.numeric(stats::residuals(crtrend:::joinpoint_engine(fit)))
  whitened <- as.numeric(crtrend:::joinpoint_whiten(analysis, residuals))
  uniform <- crtrend:::nci_mt_uniforms(nrow(data) * 5001L, 10000L)
  candidates <- expand.grid(
    scale = c("raw", "whitened"),
    correction = c("none", "hc2", "hc3"),
    stringsAsFactors = FALSE
  )
  values <- lapply(seq_len(nrow(candidates)), function(index) {
    scale <- candidates$scale[index]
    source <- if (scale == "raw") residuals else whitened
    summary <- draw_summary(
      fit,
      source,
      scale,
      candidates$correction[index],
      uniform
    )
    target <- reference[
      reference$case == sub("\\.txt$", "", basename(input)) &
        reference$statistic == "APC",
    ]
    data.frame(
      case = sub("\\.txt$", "", basename(input)),
      scale = scale,
      correction = candidates$correction[index],
      lower = summary[["lower"]],
      upper = summary[["upper"]],
      p_value = summary[["p_value"]],
      absolute_error = sum(abs(
        summary - c(target$lower, target$upper, target$p_value)
      ))
    )
  })
  do.call(rbind, values)
})
result <- do.call(rbind, rows)
result <- result[order(result$case, result$absolute_error), ]
print(do.call(rbind, lapply(split(result, result$case), head, 3L)), digits = 7)
