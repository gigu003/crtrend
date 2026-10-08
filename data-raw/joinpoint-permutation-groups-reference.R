# Import actual NCI exports, never generate official values from R fits.
root <- "/Users/qc0824/Documents/IMS/Joinpoint/Selection-audit-results/Perm-groups"
destination <- "inst/extdata/joinpoint-6.1.0-permutation-groups"
runs <- c(
  Original = "groups-original",
  relabelled = "groups-relabelled",
  `single-0` = "single-0",
  `single-1` = "single-1",
  `single-2` = "single-2"
)
rows <- lapply(names(runs), function(run) {
  prefix <- file.path(root, run, runs[[run]])
  read_export <- function(kind) {
    read.delim(paste0(prefix, ".Export.", kind, ".txt"), check.names = FALSE)
  }
  summary <- readLines(
    paste0(prefix, ".Export.Run.Summary.txt"),
    skipNul = TRUE,
    warn = FALSE
  )
  stopifnot(
    any(grepl("Seed for Randomly Permuting Data\\s+7160", summary)),
    any(grepl("Number of Permutations for Each Test\\s+4499", summary))
  )
  x <- read_export("Model.Selection")
  model <- read_export("Final.Selected.Model")
  grouped <- "cohort" %in% names(x)
  sex <- if (grouped) {
    x$cohort
  } else {
    rep(as.integer(sub("single-", "", run)), nrow(x))
  }
  if (run == "relabelled") {
    sex <- 2L - sex
  }
  final <- if (grouped) {
    model$Model[match(x$cohort, model$cohort)]
  } else {
    rep(model$Model, nrow(x))
  }
  ref <- data.frame(
    run,
    sex,
    step = x[["Test #"]] + 1L,
    k_null = x[["Null Hypothesis"]],
    k_alt = x[["Alt Hypothesis"]],
    selected = x[["Selected Hypothesis"]],
    numerator_df = x[["Numerator DF"]],
    denominator_df = x[["Denominator DF"]],
    n_total = x[["# of Permutations"]],
    p_value = x[["P-Value"]],
    alpha = x[["Significance level"]],
    final_k = final
  )
  baseline <- read.csv(file.path(destination, "r-groups-original.csv"))
  idx <- match(paste(ref$sex, ref$step), paste(baseline$cohort, baseline$step))
  stopifnot(!anyNA(idx), nrow(ref) == if (grouped) 12L else 4L)
  r <- baseline[idx, ]
  stopifnot(
    all(ref$k_null == r$k_null),
    all(ref$k_alt == r$k_alt),
    all(ref$p_value == round(r$p_value, 7)),
    all(ref$alpha == round(r$alpha, 7)),
    all(ref$selected == ifelse(r$reject, r$k_alt, r$k_null)),
    all(ref$n_total == r$n_permutations + 1L),
    all(ref$final_k == 1L),
    all(ref$numerator_df == 2 * (r$k_alt - r$k_null)),
    all(ref$denominator_df == 22 - 2 * r$k_alt - 2)
  )
  # Confirm cohort mapping using observed data, not filenames alone.
  actual <- read_export("Data")
  input <- read.delim(file.path(destination, paste0(runs[[run]], ".txt")))
  key <- function(x) paste(if ("cohort" %in% names(x)) x$cohort else 0, x$year)
  if (!grouped) {
    input$cohort <- NULL
  }
  j <- match(key(actual), key(input))
  stopifnot(
    nrow(actual) == nrow(input),
    !anyNA(j),
    !anyDuplicated(key(actual)),
    all(actual[["Observed rate"]] == round(input$rate[j], 2)),
    all(actual[["Standard Error"]] == round(input$se[j], 4))
  )
  ref
})
reference <- do.call(rbind, rows)
stopifnot(nrow(reference) == 36L)
write.csv(
  reference,
  "tests/testthat/fixtures/joinpoint-6.1-permutation-groups.csv",
  row.names = FALSE
)
message(
  "All 36 native comparisons match R; all 9 cohort-level final models match."
)
