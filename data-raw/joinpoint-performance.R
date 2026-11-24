# Run from the package root. Compare snapshots in separate R processes:
# Rscript data-raw/joinpoint-performance.R /tmp/before-R /tmp/before
# Rscript data-raw/joinpoint-performance.R R /tmp/after
# The input remains local and is not included in the benchmark output.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Supply an R source directory and an output prefix.")
}
source_dir <- args[1L]
output_prefix <- args[2L]
implementation <- new.env(parent = globalenv())
for (file in list.files(source_dir, pattern = "[.]R$", full.names = TRUE)) {
  sys.source(file, envir = implementation)
}
data <- readRDS("inst/extdata/risk_inci.rds")
invisible(gc())
Rprof(paste0(output_prefix, ".prof"), interval = 0.02)
timing <- system.time(
  fit <- implementation$joinpoint(
    data,
    time = "year",
    outcome = "risk",
    group_vars = c("sex", "cancer"),
    se = "se"
  )
)
Rprof(NULL)
print(timing)
print(head(summaryRprof(paste0(output_prefix, ".prof"))$by.total, 20L))
saveRDS(
  list(
    timing = timing,
    session = utils::sessionInfo(),
    model_info = fit$model_info,
    joinpoints = fit$joinpoints,
    fitted = fit$fitted.values,
    candidates = lapply(fit$models, function(model) model$selection$candidates),
    resamples = lapply(fit$models, function(model) model$resamples)
  ),
  paste0(output_prefix, ".rds")
)
