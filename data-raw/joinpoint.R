# Build the Joinpoint 6.1.0 sample data shipped with crtrend.

source_dir <- file.path("inst", "extdata", "joinpoint-6.1.0")

if (!dir.exists(source_dir)) {
  stop(
    "Cannot find the Joinpoint sample-data directory: ",
    source_dir,
    call. = FALSE
  )
}

crude <- read.delim(
  file.path(source_dir, "Sample.Crude.Rate.Calculation.Data.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
jp_crude <- data.frame(
  sex = as.integer(crude[["Sex"]]),
  year = as.integer(crude[[2L]]),
  count = as.integer(crude[["Count"]]),
  population = as.integer(crude[["Population"]])
)

age_adjusted <- read.delim(
  file.path(source_dir, "Sample.Age.Adjusted.Rate.Calculation.Data.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
jp_age_adjusted <- data.frame(
  sex = as.integer(age_adjusted[["Sex"]]),
  year = as.integer(age_adjusted[[2L]]),
  age_group = as.integer(age_adjusted[["Age recode"]]),
  count = as.integer(age_adjusted[["Count"]]),
  population = as.integer(age_adjusted[["Population"]]),
  standard_population = as.integer(age_adjusted[["Std Pop"]])
)

seer <- read.table(
  file.path(source_dir, "Sample.Data.txt"),
  sep = ";",
  header = FALSE,
  col.names = c(
    "sex",
    "year_index",
    "rate",
    "se",
    "lower_ci",
    "upper_ci",
    "count",
    "population"
  ),
  stringsAsFactors = FALSE
)
jp_data <- data.frame(
  sex = as.integer(seer$sex),
  year_index = as.integer(seer$year_index),
  year = as.integer(seer$year_index) + 2000L,
  rate = as.numeric(seer$rate),
  se = as.numeric(seer$se),
  lower_ci = as.numeric(seer$lower_ci),
  upper_ci = as.numeric(seer$upper_ci),
  count = as.integer(seer$count),
  population = as.integer(seer$population)
)

jump <- read.table(
  file.path(source_dir, "Sample.Jump.CR.Standard.Data.txt"),
  sep = ";",
  header = TRUE,
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
jp_jump <- data.frame(
  year_index = as.integer(jump[[1L]]),
  year = as.integer(jump[[1L]]) + 1992L,
  rate = as.numeric(jump[[2L]]),
  se = as.numeric(jump[[3L]]),
  lower_ci = as.numeric(jump[[4L]]),
  upper_ci = as.numeric(jump[[5L]]),
  count = as.integer(jump[[6L]]),
  population = as.integer(jump[[7L]])
)

save(jp_crude, file = file.path("data", "jp_crude.rda"), version = 3)
save(
  jp_age_adjusted,
  file = file.path("data", "jp_age_adjusted.rda"),
  version = 3
)
save(jp_data, file = file.path("data", "jp_data.rda"), version = 3)
save(jp_jump, file = file.path("data", "jp_jump.rda"), version = 3)
