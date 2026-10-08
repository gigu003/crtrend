# Diagnostic only; no production RNG changes.
# Static evidence: Permute 0x471980, shuffle 0x481e40,
# integer draw 0x4679d0, seed expansion 0x4c7b40.
devtools::load_all(quiet = TRUE)
d <- subset(jp_data, sex == 0, c("year", "rate", "se"))
fit <- joinpoint(d, "year", "rate", se = "se", max_jp = 4, n_eq = 2)
d <- crtrend:::joinpoint_analysis_data(fit$models[[1]])
candidates <- lapply(
  0:4,
  crtrend:::fit_candidate,
  data = d,
  weighted = TRUE,
  min_observations = 2L,
  selection_method = "permutation"
)
n <- nrow(d)
b <- 4499L
uniforms <- crtrend:::nci_mt_uniforms((n - 1L) * b * 4L, 7160)
pairs <- list(c(0, 4), c(0, 3), c(0, 2), c(1, 2))
results <- list()
for (reset_stream in c(FALSE, TRUE)) {
  for (step in seq_along(pairs)) {
    pair <- pairs[[step]]
    reduced <- candidates[[pair[1] + 1L]]
    complex <- candidates[[pair[2] + 1L]]
    residual <- as.numeric(crtrend:::joinpoint_whiten(
      d,
      residuals(reduced$model)
    ))
    observed <- crtrend:::permutation_f_statistic(
      reduced$sse,
      complex$sse,
      pair[1],
      pair[2],
      n
    )
    offset <- if (reset_stream) 0L else (step - 1L) * b * (n - 1L)
    shuffled <- matrix(NA_real_, n, b)
    for (draw in seq_len(b)) {
      for (i in seq_len(n - 1L)) {
        offset <- offset + 1L
        j <- i + floor(uniforms[offset] * (n - i + 1L))
        temporary <- residual[i]
        residual[i] <- residual[j]
        residual[j] <- temporary
      }
      shuffled[, draw] <- residual
    }
    count <- 0L
    batch_size <- crtrend:::permutation_batch_size(reduced, complex)
    for (first in seq.int(1L, b, by = batch_size)) {
      index <- first:min(b, first + batch_size - 1L)
      response <- fitted(reduced$model) +
        crtrend:::joinpoint_unwhiten(d, shuffled[, index, drop = FALSE])
      sse0 <- crtrend:::candidate_sse_batch(reduced, response, d, TRUE)
      sse1 <- crtrend:::candidate_sse_batch(complex, response, d, TRUE)
      statistic <- crtrend:::permutation_f_statistic(
        sse0,
        sse1,
        pair[1],
        pair[2],
        n
      )
      stopifnot(all(is.finite(statistic)))
      count <- count + sum(statistic >= observed)
    }
    row <- data.frame(
      reset_stream,
      step,
      k_null = pair[1],
      k_alt = pair[2],
      exceedances = count,
      p_value = (count + 1) / (b + 1)
    )
    print(row)
    results[[length(results) + 1L]] <- row
  }
}
write.csv(
  do.call(rbind, results),
  "inst/extdata/joinpoint-6.1.0-selection-audit/native-permutation-diagnostic.csv",
  row.names = FALSE
)
