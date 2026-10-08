# Independently diagnose and validate the opt-in production WBIC calculation.
# Run from the package root: Rscript data-raw/joinpoint-v6.1-weight-audit.R
devtools::load_all(quiet = TRUE)

# Match the scalar double-precision multiply/add order in the native program.
# R's BLAS and sum() are deliberately not used in this diagnostic kernel.
audit_multiply <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  stopifnot(ncol(a) == nrow(b))
  value <- matrix(0, nrow(a), ncol(b))
  for (i in seq_len(nrow(a))) {
    for (j in seq_len(ncol(b))) {
      accumulator <- 0
      for (h in seq_len(ncol(a))) {
        accumulator <- accumulator + a[i, h] * b[h, j]
      }
      value[i, j] <- accumulator
    }
  }
  value
}

# Nonsingular 2 x 2 specialization of native 0x50b550 (via 0x50be70).
# This is NOT a general replacement for solve(), nor a stable recommended solver.
audit_inverse <- function(a) {
  stopifnot(identical(dim(a), c(2L, 2L)), a[1, 1] > 0)
  first <- 1 / a[1, 1]
  lower <- a[1, 2] * first
  pivot <- a[2, 2] - lower * a[1, 2]
  stopifnot(pivot > 0)
  second <- 1 / pivot
  off_diagonal <- 0 - second * lower
  matrix(
    c(first - off_diagonal * lower, off_diagonal, off_diagonal, second),
    nrow = 2L
  )
}

audit_local_weight <- function(data, knots, index, origin = 0) {
  lower <- if (index == 1L) min(data$trend_time) - 1 else knots[index - 1L]
  upper <- if (index == length(knots)) {
    max(data$trend_time)
  } else {
    knots[index + 1L]
  }
  local <- data[
    data$trend_time > lower & data$trend_time <= upper,
    ,
    drop = FALSE
  ]
  time <- local$trend_time - origin
  hinge <- pmax(time - (knots[index] - origin), 0)
  # Native input processing: SE squared, divided by rate squared, reciprocal.
  log_variance <- local$trend_se *
    local$trend_se /
    (local$trend_outcome * local$trend_outcome)
  multiplier <- sqrt(1 / log_variance)
  x <- cbind(1, time) * multiplier
  y <- matrix(local$trend_log * multiplier)
  z <- matrix(hinge * multiplier)
  normal <- audit_multiply(t(x), x)
  inverse <- audit_inverse(normal)
  projection <- audit_multiply(audit_multiply(x, inverse), t(x))
  residual_projection <- diag(nrow(x)) - projection
  yy <- audit_multiply(audit_multiply(t(y), residual_projection), y)[1, 1]
  zz <- audit_multiply(audit_multiply(t(z), residual_projection), z)[1, 1]
  yz <- audit_multiply(audit_multiply(t(y), residual_projection), z)[1, 1]
  solver <- qr(x)
  ry <- y - qr.fitted(solver, y)
  rz <- z - qr.fitted(solver, z)
  qr_native_weights <- sum(ry * rz)^2 / (sum(ry^2) * sum(rz^2))
  data.frame(
    k = length(knots),
    index = index,
    knot = knots[index],
    first_year = min(local$trend_time),
    last_year = max(local$trend_time),
    n = nrow(local),
    origin = origin,
    package_qr = crtrend:::partial_joinpoint_r_squared(
      index,
      knots,
      data,
      TRUE
    ),
    native_order = yz * yz / (yy * zz),
    qr_native_weights = qr_native_weights,
    condition_normal = kappa(normal, exact = TRUE),
    projection_idempotence_error = max(abs(
      audit_multiply(projection, projection) - projection
    ))
  )
}

# Independent kernel checks on a well-conditioned matrix.
test_matrix <- matrix(c(3, 1, 1, 2), 2)
stopifnot(
  max(abs(audit_inverse(test_matrix) - solve(test_matrix))) < 1e-15,
  identical(
    audit_multiply(test_matrix, test_matrix),
    test_matrix %*% test_matrix
  )
)
fit <- joinpoint(
  subset(jp_data, sex == 0),
  "year",
  "rate",
  se = "se",
  max_jp = 4,
  n_eq = 2
)
analysis <- fit$models[[1]]$.analysis
candidates <- lapply(0:4, function(k) {
  crtrend:::fit_candidate(k, analysis, TRUE, 2L, "weighted_bic")
})
local <- do.call(
  rbind,
  lapply(candidates[-1], function(candidate) {
    do.call(
      rbind,
      lapply(seq_along(candidate$joinpoints), function(index) {
        rbind(
          audit_local_weight(analysis, candidate$joinpoints, index),
          audit_local_weight(
            analysis,
            candidate$joinpoints,
            index,
            origin = 2000
          )
        )
      })
    )
  })
)

# Transcribed from the user's Joinpoint 6.1.0 Model tab screenshot (7 decimals).
reference_weight <- c(NA, .5184867, .6844505, .9146089, .9116191)
reference_wbic <- c(4.2970652, 3.9200959, 4.2041850, 4.6405566, 5.0180119)
comparison <- do.call(
  rbind,
  lapply(seq_along(candidates), function(position) {
    candidate <- candidates[[position]]
    k <- position - 1L
    rows <- local[local$k == k & local$origin == 0, , drop = FALSE]
    weight <- if (k == 0L) NA_real_ else max(rows$native_order)
    wbic <- if (k == 0L) {
      candidate$bic
    } else {
      candidate$bic * (1 - weight) + candidate$bic3 * weight
    }
    data.frame(
      k = k,
      joinpoints = paste(candidate$joinpoints, collapse = ";"),
      maximizing_local_index = if (k == 0L) {
        NA_integer_
      } else {
        which.max(rows$native_order)
      },
      screenshot_weight = reference_weight[position],
      package_weight = if (k == 0L) NA_real_ else candidate$partial_r_squared,
      native_order_weight = weight,
      weight_display_match = if (k == 0L) {
        NA
      } else {
        round(weight, 7) == reference_weight[position]
      },
      screenshot_wbic = reference_wbic[position],
      package_wbic = candidate$weighted_bic,
      native_weight_wbic = wbic,
      wbic_display_match = round(wbic, 7) == reference_wbic[position]
    )
  })
)

compatible_fit <- joinpoint(
  subset(jp_data, sex == 0),
  "year",
  "rate",
  se = "se",
  max_jp = 4,
  n_eq = 2,
  wbic_weight = "joinpoint_6.1"
)
production <- calc_model(compatible_fit, ci_method = "Param")$candidates
comparison$production_weight <- production$weight
comparison$production_wbic <- production$weighted_bic

# Validate the reproduction without relaxing the screenshot's rounding criterion.
# With centered years, explicit projection and QR agree to 1e-9 on all locals.
stopifnot(
  identical(comparison$production_weight, comparison$native_order_weight),
  identical(comparison$production_wbic, comparison$native_weight_wbic),
  all(comparison$weight_display_match[-1]),
  all(comparison$wbic_display_match),
  max(abs(local$qr_native_weights - local$package_qr)) < 1e-10,
  max(abs(with(subset(local, origin == 2000), native_order - package_qr))) <
    1e-9,
  which.min(comparison$package_wbic) == which.min(comparison$native_weight_wbic)
)
destination <- "inst/extdata/joinpoint-6.1.0-model-audit"
write.csv(
  local,
  file.path(destination, "weight_local_audit.csv"),
  row.names = FALSE
)
write.csv(
  comparison,
  file.path(destination, "weight_reproduction.csv"),
  row.names = FALSE
)
print(comparison, digits = 13, row.names = FALSE)
cat(
  "Maximum QR change from native weights:",
  max(abs(local$qr_native_weights - local$package_qr)),
  "\n"
)
cat(
  "Maximum centered projection versus QR difference:",
  max(abs(with(subset(local, origin == 2000), native_order - package_qr))),
  "\n"
)
cat(
  "Normal-matrix condition number range (original years):",
  range(subset(local, origin == 0)$condition_normal),
  "\n"
)
