devtools::load_all(quiet = TRUE)

output_dir <- file.path(
  "inst",
  "extdata",
  "joinpoint-6.1.0-dds-disagreement"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

time <- 0:30
make_case <- function(first_change, second_change) {
  data.frame(
    year = time,
    rate = exp(
      4 - .02 * time + first_change * pmax(time - 10, 0) +
        second_change * pmax(time - 22, 0) + .01 * sin(2 * time)
    )
  )
}

cases <- list(
  bic = make_case(.02, .003),
  bic3 = make_case(.005, .00325)
)
for (case in names(cases)) {
  utils::write.csv(
    cases[[case]],
    file.path(output_dir, paste0("dds-selects-", case, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
}

results <- lapply(cases, function(data) {
  fit <- joinpoint(
    data,
    time = "year",
    outcome = "rate",
    max_jp = 3,
    selection_method = "DDS",
    n_eq = 2
  )
  output <- calc_model(fit, ci_method = "Param")
  list(fit = fit, output = output)
})
summary <- do.call(rbind, lapply(names(results), function(case) {
  fit <- results[[case]]$fit
  output <- results[[case]]$output
  data.frame(
    case = case,
    bic_minimizer_k = output$candidates$n_joinpoints[
      which.min(output$candidates$bic)
    ],
    bic3_minimizer_k = output$candidates$n_joinpoints[
      which.min(output$candidates$bic3)
    ],
    dds_branch = output$selection$dds_method,
    selected_k = output$selection$selected_n_joinpoints,
    bic_delta = output$selection$dds_bic_delta,
    bic3_delta = output$selection$dds_bic3_delta,
    delta_min = output$selection$dds_delta_min,
    delta_max = output$selection$dds_delta_max,
    sigma_squared = output$selection$dds_sigma_squared,
    variance_df = output$selection$dds_variance_df
  )
}))
utils::write.csv(
  summary,
  file.path(output_dir, "crtrend-reference.csv"),
  row.names = FALSE,
  quote = FALSE
)
