# Build diagnostic inputs for distinguishing candidate Joinpoint 6.1.0
# empirical-quantile residual corrections.

output_dir <- file.path(
  "inst",
  "extdata",
  "joinpoint-6.1.0-eq-audit"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
old_options <- options(digits = 12)

make_audit_data <- function(
  year,
  joinpoint,
  intercept,
  first_slope,
  second_slope,
  log_se,
  standardized_error
) {
  mean_log_rate <- intercept +
    first_slope * year +
    (second_slope - first_slope) * pmax(year - joinpoint, 0)
  rate <- exp(mean_log_rate + log_se * standardized_error)
  data.frame(
    Year = year,
    Rate = rate,
    `Standard Error` = rate * log_se,
    check.names = FALSE
  )
}

heteroscedastic_log_se <- c(
  0.008,
  0.012,
  0.02,
  0.035,
  0.06,
  0.018,
  0.01,
  0.045,
  0.012,
  0.07,
  0.02,
  0.01,
  0.05,
  0.015
)
heteroscedastic_data <- make_audit_data(
  year = 0:13,
  joinpoint = 5,
  intercept = log(60),
  first_slope = -0.025,
  second_slope = 0.015,
  log_se = heteroscedastic_log_se,
  standardized_error = c(
    1.5,
    -1,
    0.8,
    -0.6,
    0.4,
    -0.5,
    0.7,
    -0.8,
    0.6,
    -0.4,
    0.9,
    -0.7,
    0.5,
    -0.3
  )
)

zero_joinpoint_data <- make_audit_data(
  year = 0:17,
  joinpoint = 8,
  intercept = log(50),
  first_slope = -0.0015,
  second_slope = -0.0015,
  log_se = rep(0.02, 18),
  standardized_error = c(
    -0.8,
    0.6,
    -0.4,
    0.7,
    -0.5,
    0.2,
    -0.3,
    0.4,
    -0.2,
    0.1,
    -0.1,
    0.3,
    -0.4,
    0.5,
    -0.6,
    0.4,
    -0.2,
    0.3
  )
)
zero_joinpoint_log_se <- c(
  0.008,
  0.014,
  0.025,
  0.05,
  0.012,
  0.035,
  0.01,
  0.06,
  0.018,
  0.04,
  0.009,
  0.03,
  0.055,
  0.015,
  0.022,
  0.045,
  0.011,
  0.028
)

audit_data <- list(
  `01-equal-log-se-high-leverage.txt` = make_audit_data(
    year = 0:11,
    joinpoint = 2,
    intercept = log(50),
    first_slope = -0.04,
    second_slope = 0.025,
    log_se = rep(0.02, 12),
    standardized_error = c(
      2.8,
      -1.4,
      0.5,
      -0.6,
      0.2,
      0.4,
      -0.3,
      0.1,
      0.2,
      -0.2,
      0.1,
      -0.1
    )
  ),
  `02-heteroscedastic-whitening.txt` = heteroscedastic_data,
  `03-endpoint-outlier.txt` = make_audit_data(
    year = 0:9,
    joinpoint = 2,
    intercept = log(45),
    first_slope = -0.08,
    second_slope = 0.03,
    log_se = rep(0.018, 10),
    standardized_error = c(
      -4.5,
      1.2,
      -0.2,
      0.3,
      -0.4,
      0.2,
      0.1,
      -0.2,
      0.3,
      -0.1
    )
  ),
  `04-paired-equal-log-se.txt` = transform(
    heteroscedastic_data,
    `Standard Error` = Rate * 0.02
  ),
  `05-paired-reversed-log-se.txt` = transform(
    heteroscedastic_data,
    `Standard Error` = Rate * rev(heteroscedastic_log_se)
  ),
  `06-paired-scaled-log-se.txt` = transform(
    heteroscedastic_data,
    `Standard Error` = `Standard Error` * 10
  ),
  `07-zero-jp-equal-log-se.txt` = zero_joinpoint_data,
  `08-zero-jp-endpoint-leverage.txt` = make_audit_data(
    year = 0:17,
    joinpoint = 8,
    intercept = log(50),
    first_slope = -0.0015,
    second_slope = -0.0015,
    log_se = rep(0.02, 18),
    standardized_error = c(
      2.0,
      -0.9,
      -0.4,
      0.5,
      -0.3,
      0.2,
      -0.2,
      0.3,
      -0.1,
      0.1,
      -0.1,
      0.2,
      -0.3,
      0.4,
      -0.5,
      0.3,
      -0.2,
      0.1
    )
  ),
  `09-zero-jp-heteroscedastic.txt` = transform(
    zero_joinpoint_data,
    `Standard Error` = Rate * zero_joinpoint_log_se
  ),
  `10-zero-jp-scaled-log-se.txt` = transform(
    zero_joinpoint_data,
    `Standard Error` = Rate * zero_joinpoint_log_se * 10
  ),
  `11-zero-jp-17-observations.txt` = make_audit_data(
    year = 0:16,
    joinpoint = 8,
    intercept = log(50),
    first_slope = -0.0015,
    second_slope = -0.0015,
    log_se = rep(0.02, 17),
    standardized_error = c(
      -0.8,
      0.6,
      -0.4,
      0.7,
      -0.5,
      0.2,
      -0.3,
      0.4,
      -0.2,
      0.1,
      -0.1,
      0.3,
      -0.4,
      0.5,
      -0.6,
      0.4,
      -0.2
    )
  ),
  `12-zero-jp-19-observations.txt` = make_audit_data(
    year = 0:18,
    joinpoint = 8,
    intercept = log(50),
    first_slope = -0.0015,
    second_slope = -0.0015,
    log_se = rep(0.02, 19),
    standardized_error = c(
      -0.8,
      0.6,
      -0.4,
      0.7,
      -0.5,
      0.2,
      -0.3,
      0.4,
      -0.2,
      0.1,
      -0.1,
      0.3,
      -0.4,
      0.5,
      -0.6,
      0.4,
      -0.2,
      0.3,
      0.1
    )
  ),
  `13-zero-jp-near-endpoint.txt` = make_audit_data(
    year = 0:17,
    joinpoint = 8,
    intercept = log(50),
    first_slope = -0.0015,
    second_slope = -0.0015,
    log_se = rep(0.02, 18),
    standardized_error = c(
      1.0,
      -0.9,
      -0.4,
      0.5,
      -0.3,
      0.2,
      -0.2,
      0.3,
      -0.1,
      0.1,
      -0.1,
      0.2,
      -0.3,
      0.4,
      -0.5,
      0.3,
      -0.2,
      0.1
    )
  ),
  `14-zero-jp-outside-endpoint.txt` = make_audit_data(
    year = 0:17,
    joinpoint = 8,
    intercept = log(50),
    first_slope = -0.0015,
    second_slope = -0.0015,
    log_se = rep(0.02, 18),
    standardized_error = c(
      1.5,
      -0.9,
      -0.4,
      0.5,
      -0.3,
      0.2,
      -0.2,
      0.3,
      -0.1,
      0.1,
      -0.1,
      0.2,
      -0.3,
      0.4,
      -0.5,
      0.3,
      -0.2,
      0.1
    )
  )
)

invisible(lapply(names(audit_data), function(filename) {
  write.table(
    audit_data[[filename]],
    file = file.path(output_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}))
options(old_options)
