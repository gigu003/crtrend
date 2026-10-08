# Opt-in numerical compatibility kernels; not shared by fitting or inference.
# Evidence: inst/extdata/joinpoint-6.1.0-model-audit/WEIGHT-DIAGNOSIS.md.
joinpoint_61_weight_error <- function() {
  stop(
    'Joinpoint 6.1 Weight projection is singular or invalid; use `wbic_weight = "qr"`.',
    call. = FALSE
  )
}

# Keep each multiplication and addition separate, in increasing index order.
joinpoint_61_multiply <- function(a, b) {
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

# Nonsingular 2 x 2 specialization, native 0x50b550 via 0x50be70.
joinpoint_61_inverse <- function(a) {
  if (any(!is.finite(a)) || a[1, 1] <= 0) {
    joinpoint_61_weight_error()
  }
  first <- 1 / a[1, 1]
  lower <- a[1, 2] * first
  pivot <- a[2, 2] - lower * a[1, 2]
  if (!is.finite(pivot) || pivot <= 0) {
    joinpoint_61_weight_error()
  }
  second <- 1 / pivot
  off_diagonal <- 0 - second * lower
  matrix(
    c(first - off_diagonal * lower, off_diagonal, off_diagonal, second),
    nrow = 2L
  )
}

joinpoint_61_partial_weight <- function(data, hinge) {
  variance <- (data$trend_se * data$trend_se) /
    (data$trend_outcome * data$trend_outcome)
  multiplier <- sqrt(1 / variance)
  if (any(!is.finite(multiplier)) || any(multiplier <= 0)) {
    joinpoint_61_weight_error()
  }
  x <- cbind(1, data$trend_time) * multiplier
  y <- matrix(data$trend_log * multiplier)
  z <- matrix(hinge * multiplier)
  inverse <- joinpoint_61_inverse(joinpoint_61_multiply(t(x), x))
  projection <- joinpoint_61_multiply(
    joinpoint_61_multiply(x, inverse),
    t(x)
  )
  residual_projection <- diag(nrow(x)) - projection
  yy <- joinpoint_61_multiply(
    joinpoint_61_multiply(t(y), residual_projection),
    y
  )[1L, 1L]
  zz <- joinpoint_61_multiply(
    joinpoint_61_multiply(t(z), residual_projection),
    z
  )[1L, 1L]
  yz <- joinpoint_61_multiply(
    joinpoint_61_multiply(t(y), residual_projection),
    z
  )[1L, 1L]
  denominator <- yy * zz
  value <- yz * yz / denominator
  if (!is.finite(value) || !is.finite(denominator) || yy <= 0 || zz <= 0) {
    joinpoint_61_weight_error()
  }
  max(0, min(1, value))
}
