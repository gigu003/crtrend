# Compare independently decoded EQ endpoint/interpolation branches with exports.
# Run from the package root. This script does not invoke Joinpoint or CrossOver.
devtools::load_all(quiet = TRUE)
audit_dir <- file.path("inst", "extdata", "joinpoint-6.1.0-eq-audit")
reference <- read.csv(file.path(audit_dir, "nci-6.1-eq-reference.csv"))
uniform_cache <- new.env(parent = emptyenv())
rows <- list()
for (number in 7:14) {
  input <- list.files(
    audit_dir,
    paste0("^", sprintf("%02d", number), ".*txt$"),
    full.names = TRUE
  )
  data <- read.delim(input, check.names = FALSE)
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
  engine <- crtrend:::joinpoint_engine(fit)
  residuals <- sort(as.numeric(crtrend:::joinpoint_whiten(
    analysis,
    residuals(engine)
  )))
  n <- length(residuals)
  key <- as.character(n)
  if (!exists(key, uniform_cache, inherits = FALSE)) {
    assign(key, crtrend:::nci_mt_uniforms(n * 5001L, 10000L), uniform_cache)
  }
  uniform <- get(key, uniform_cache, inherits = FALSE)
  target <- subset(
    reference,
    case == sub(".txt$", "", basename(input)) & statistic == "APC"
  )
  design <- crtrend:::joinpoint_whiten(
    analysis,
    crtrend:::joinpoint_design(numeric(), analysis)
  )
  solver <- qr(design)
  for (upper_rule in c("previous", "decoded")) {
    endpoints <- crtrend:::nci_empirical_endpoints(residuals)
    if (upper_rule == "previous") {
      candidate <- median(residuals) +
        5.189 * median(abs(residuals - median(residuals)))
      endpoints[2L] <- if (max(residuals) <= candidate) {
        max(residuals) + 2 * IQR(residuals) / n
      } else {
        candidate
      }
    }
    knots <- sort(c(endpoints[1L], residuals, endpoints[2L]))
    for (offset in c(0, -1)) {
      interval <- floor(uniform * (n + 1L))
      fraction <- uniform * (n + 1L) - interval + offset
      sampled <- matrix(
        knots[interval + 1L] + diff(knots)[interval + 1L] * fraction,
        nrow = n
      )
      # Linear fitting in whitened coordinates avoids repeated refits.
      coefficients <- qr.coef(solver, sampled)
      draws <- 100 * expm1(coef(engine)[2L] + coefficients[2L, ])
      limits <- quantile(draws, c(.025, .975), names = FALSE)
      p <- 2 * min(sum(draws < 0), sum(draws >= 0)) / length(draws)
      rows[[length(rows) + 1L]] <- data.frame(
        case = number,
        upper_rule,
        offset,
        lower = limits[1L],
        upper = limits[2L],
        p_value = p,
        error = sum(abs(
          c(limits, p) - c(target$lower, target$upper, target$p_value)
        ))
      )
    }
  }
}
result <- do.call(rbind, rows)
print(result, digits = 8, row.names = FALSE)
