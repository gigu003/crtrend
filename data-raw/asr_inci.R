# Build the example data shipped with crtrend.

source_file <- file.path(
  "..",
  "..",
  "..",
  "1_Projects",
  "1990-2023中国患癌风险趋势",
  "output",
  "asr_inci.rds"
)

if (!file.exists(source_file)) {
  stop("Cannot find the source file: ", source_file, call. = FALSE)
}

asr_inci <- readRDS(source_file)

usethis::use_data(
  asr_inci,
  overwrite = TRUE,
  compress = "xz",
  version = 3
)
