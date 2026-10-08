devtools::load_all(quiet = TRUE)

audit_dir <- file.path("inst", "extdata", "joinpoint-6.1.0-dds-disagreement")
native <- read.csv(file.path(audit_dir, "native-dds-results.csv"))
records <- list()

for (case in native$case) {
  input <- read.csv(file.path(audit_dir, paste0("dds-selects-", case, ".csv")))
  fit <- joinpoint(
    input,
    "year",
    "rate",
    max_jp = 3,
    selection_method = "DDS",
    n_eq = 2
  )
  model <- fit$models[[1L]]
  analysis_data <- joinpoint_analysis_data(model)
  selection <- model$selection
  native_case <- native[native$case == case, ]

  for (criterion in c("bic", "bic3")) {
    candidate_index <- which.min(selection$candidates[[criterion]])
    n_joinpoints <- selection$candidates$n_joinpoints[candidate_index]
    candidate <- fit_candidate(
      n_joinpoints,
      analysis_data,
      weighted = FALSE,
      min_observations = 2L,
      selection_method = "DDS"
    )
    knots <- candidate$joinpoints
    coefficients <- stats::coef(candidate$model)
    changes <- coefficients[seq.int(3L, n_joinpoints + 2L)]

    for (index in seq_len(n_joinpoints)) {
      previous <- if (index == 1L) {
        min(analysis_data$trend_time) - 1
      } else {
        knots[index - 1L]
      }
      following <- if (index == n_joinpoints) {
        max(analysis_data$trend_time)
      } else {
        knots[index + 1L]
      }
      keep <- analysis_data$trend_time > previous &
        analysis_data$trend_time <= following
      time <- analysis_data$trend_time[keep]
      z <- pmax(time - knots[index], 0)
      design <- cbind(1, time)
      decomposition <- qr(design)
      q <- qr.Q(decomposition)
      qr_residual <- z - qr.fitted(decomposition, z)
      q_residual <- z - q %*% crossprod(q, z)
      xz <- crossprod(design, z)
      normal_sse <- sum(z^2) -
        as.numeric(crossprod(xz, solve(crossprod(design), xz)))
      centered_design <- cbind(1, time - mean(time))
      centered_residual <- z - qr.fitted(qr(centered_design), z)
      lm_residual <- stats::residuals(stats::lm(z ~ time))
      local_sse <- c(
        qr = sum(qr_residual^2),
        explicit_q = sum(q_residual^2),
        normal_equations = normal_sse,
        centered_qr = sum(centered_residual^2),
        lm = sum(lm_residual^2)
      )
      ref_delta <- native_case[[paste0("delta_", criterion)]]
      records[[length(records) + 1L]] <- data.frame(
        case = case,
        criterion = criterion,
        n_joinpoints = n_joinpoints,
        local_joinpoint = index,
        joinpoint = knots[index],
        slope_change = unname(changes[index]),
        sigma_squared = selection$dds_sigma_squared,
        local_sse_qr = unname(local_sse["qr"]),
        local_sse_explicit_q = unname(local_sse["explicit_q"]),
        local_sse_normal_equations = unname(local_sse["normal_equations"]),
        local_sse_centered_qr = unname(local_sse["centered_qr"]),
        local_sse_lm = unname(local_sse["lm"]),
        delta_crtrend = changes[index]^2 *
          local_sse["qr"] /
          selection$dds_sigma_squared,
        delta_joinpoint_export = ref_delta,
        stringsAsFactors = FALSE
      )
    }
  }
}

diagnostic <- do.call(rbind, records)
rownames(diagnostic) <- NULL
utils::write.csv(
  diagnostic,
  file.path(audit_dir, "delta-projection-diagnostic.csv"),
  row.names = FALSE,
  quote = FALSE
)
