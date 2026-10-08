# Diagnose residual-distribution details in zero-joinpoint EQ audit cases.
# Run from the package root after devtools::load_all().

audit_dir <- file.path(
  "inst",
  "extdata",
  "joinpoint-6.1.0-eq-audit"
)
reference <- utils::read.csv(
  file.path(audit_dir, "nci-6.1-eq-reference.csv"),
  na.strings = "NA"
)

empq2_residuals <- function(residuals, size, quartile_type = 7L) {
  residuals <- sort(as.numeric(residuals))
  n <- length(residuals)
  quartiles <- if (quartile_type == 10L) {
    residuals[c(
      floor((n + 1) * 0.25),
      ceiling((n + 1) * 0.75)
    )]
  } else {
    stats::quantile(
      residuals,
      probs = c(0.25, 0.75),
      names = FALSE,
      type = quartile_type
    )
  }
  spread <- quartiles[2L] - quartiles[1L]
  spread <- max(spread, sqrt(.Machine$double.eps))
  truncation <- 1.5 * spread
  lower <- if (residuals[1L] <= quartiles[1L] - truncation) {
    residuals[1L] - 2 * spread / n
  } else {
    quartiles[1L] - truncation
  }
  upper <- if (residuals[n] >= quartiles[2L] + truncation) {
    residuals[n] + 2 * spread / n
  } else {
    quartiles[2L] + truncation
  }
  knots <- c(lower, residuals, upper)
  uniform <- stats::runif(size)
  interval <- pmin(floor(uniform * (n + 1L)), n)
  fraction <- uniform * (n + 1L) - interval
  knots[interval + 1L] +
    (knots[interval + 2L] - knots[interval + 1L]) * fraction
}

zero_joinpoint_draws <- function(
  fit,
  n_resamples,
  seed,
  center = FALSE,
  quartile_type = 7L,
  residual_scale = c("whitened", "raw"),
  leverage_adjustment = c("none", "HC2", "HC3")
) {
  residual_scale <- match.arg(residual_scale)
  leverage_adjustment <- match.arg(leverage_adjustment)
  data <- crtrend:::joinpoint_analysis_data(fit)
  engine <- crtrend:::joinpoint_engine(fit)
  residuals <- as.numeric(stats::residuals(engine))
  if (residual_scale == "whitened") {
    residuals <- as.numeric(crtrend:::joinpoint_whiten(data, residuals))
  }
  if (leverage_adjustment != "none") {
    design <- crtrend:::joinpoint_design(numeric(), data)
    whitened_design <- crtrend:::joinpoint_whiten(data, design)
    leverage <- stats::hat(whitened_design)
    exponent <- if (leverage_adjustment == "HC2") 0.5 else 1
    residuals <- residuals / (1 - leverage)^exponent
  }
  if (center) {
    residuals <- residuals - mean(residuals)
  }
  sampled <- matrix(
    if (quartile_type == 7L) {
      crtrend:::nci_empirical_residuals(
        residuals,
        nrow(data) * n_resamples,
        seed = seed
      )
    } else {
      set.seed(seed)
      empq2_residuals(
        residuals,
        nrow(data) * n_resamples,
        quartile_type = quartile_type
      )
    },
    nrow = nrow(data),
    ncol = n_resamples
  )
  sampled <- if (residual_scale == "whitened") {
    crtrend:::joinpoint_unwhiten(data, sampled)
  } else {
    sampled
  }
  responses <- matrix(
    stats::fitted(engine),
    nrow = nrow(data),
    ncol = n_resamples
  ) +
    sampled
  design <- crtrend:::joinpoint_design(numeric(), data)
  solver <- qr(crtrend:::joinpoint_whiten(data, design))
  coefficients <- qr.coef(
    solver,
    crtrend:::joinpoint_whiten(data, responses)
  )
  100 * expm1(coefficients[2L, ])
}

summarize_draws <- function(draws, final_type = 7L) {
  limits <- stats::quantile(
    draws,
    probs = c(0.025, 0.975),
    names = FALSE,
    type = final_type
  )
  c(
    lower = limits[1L],
    upper = limits[2L],
    p_value = 2 * min(sum(draws < 0), sum(draws >= 0)) / length(draws)
  )
}

input_files <- list.files(
  audit_dir,
  pattern = "^(07|08|09|10|11|12|13|14).*\\.txt$",
  full.names = TRUE
)
diagnostic_rows <- lapply(input_files, function(input_file) {
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
  )$models[[1L]]
  stopifnot(fit$n_joinpoints == 0L)
  nci <- reference[
    reference$case == case & reference$statistic == "APC",
  ]
  do.call(
    rbind,
    lapply(c("none", "HC2", "HC3"), function(leverage_adjustment) {
      do.call(
        rbind,
        lapply(c("whitened", "raw"), function(residual_scale) {
          do.call(
            rbind,
            lapply(c(FALSE, TRUE), function(center) {
              do.call(
                rbind,
                lapply(seq_len(10L), function(quartile_type) {
                  values <- summarize_draws(zero_joinpoint_draws(
                    fit,
                    n_resamples = 5001,
                    seed = 10000,
                    center = center,
                    quartile_type = quartile_type,
                    residual_scale = residual_scale,
                    leverage_adjustment = leverage_adjustment
                  ))
                  data.frame(
                    case = case,
                    leverage_adjustment = leverage_adjustment,
                    residual_scale = residual_scale,
                    center = center,
                    quartile_type = quartile_type,
                    lower = values["lower"],
                    upper = values["upper"],
                    p_value = values["p_value"],
                    absolute_error = sum(abs(
                      values -
                        c(
                          nci$lower,
                          nci$upper,
                          nci$p_value
                        )
                    ))
                  )
                })
              )
            })
          )
        })
      )
    })
  )
})

diagnostic <- do.call(rbind, diagnostic_rows)
diagnostic <- diagnostic[
  order(diagnostic$case, diagnostic$absolute_error),
]
best <- do.call(rbind, lapply(split(diagnostic, diagnostic$case), head, 4L))
rownames(best) <- NULL
print(best, digits = 7)
