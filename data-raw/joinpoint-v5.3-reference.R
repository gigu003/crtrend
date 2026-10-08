# Extract candidate fitted values from the official Joinpoint V5.3.0.0 sample
# output. The .jpo file is a binary format and this narrow reader is only for
# creating a stable regression-test fixture from the published sample output.

load(file.path("data", "jp_data.rda"))

reference_file <- file.path(
  "inst",
  "extdata",
  "joinpoint-5.3.0",
  "Sample.Session.jpo"
)

if (!file.exists(reference_file)) {
  stop(
    "Cannot find the Joinpoint V5.3 reference file: ",
    reference_file,
    call. = FALSE
  )
}

reference_raw <- readBin(
  reference_file,
  what = "raw",
  n = file.info(reference_file)$size
)

read_double <- function(offset) {
  readBin(
    reference_raw[(offset + 1L):(offset + 8L)],
    what = "double",
    n = 1L,
    size = 8L,
    endian = "little"
  )
}

find_double <- function(target, tolerance = 1e-8) {
  offsets <- vapply(
    seq_len(length(reference_raw) - 7L),
    function(offset) {
      value <- read_double(offset - 1L)
      isTRUE(is.finite(value) && abs(value - target) < tolerance)
    },
    logical(1)
  )
  which(offsets) - 1L
}

extract_group <- function(sex) {
  first_row <- which(jp_data$sex == sex)[1L]
  offsets <- find_double(jp_data$rate[first_row])
  starts <- offsets - 8L
  starts <- starts[vapply(
    starts,
    function(start) {
      start >= 0L &&
        abs(read_double(start) - 2000) < 1e-8 &&
        abs(read_double(start + 21L * 32L) - 2021) < 1e-8
    },
    logical(1)
  )]

  if (length(starts) != 5L) {
    stop(
      "Expected five candidate models for sex ",
      sex,
      "; found ",
      length(starts),
      call. = FALSE
    )
  }

  do.call(
    rbind,
    lapply(seq_along(starts), function(model) {
      start <- starts[[model]]
      data.frame(
        model = model - 1L,
        sex = sex,
        year = vapply(
          0:21,
          function(index) read_double(start + index * 32L),
          numeric(1)
        ),
        rate = vapply(
          0:21,
          function(index) read_double(start + index * 32L + 8L),
          numeric(1)
        ),
        fitted = vapply(
          0:21,
          function(index) read_double(start + index * 32L + 16L),
          numeric(1)
        ),
        se = vapply(
          0:21,
          function(index) read_double(start + index * 32L + 24L),
          numeric(1)
        )
      )
    })
  )
}

v5_3_fitted <- do.call(rbind, lapply(c(0L, 1L, 2L), extract_group))

load(file.path("data", "jp_age_adjusted.rda"))

age_groups <- split(
  seq_len(nrow(jp_age_adjusted)),
  interaction(
    jp_age_adjusted$sex,
    jp_age_adjusted$year,
    drop = TRUE,
    lex.order = FALSE
  )
)
age_rates <- do.call(
  rbind,
  lapply(age_groups, function(index) {
    data <- jp_age_adjusted[index, , drop = FALSE]
    weights <- data$standard_population / sum(data$standard_population)
    data.frame(
      sex = data$sex[1L],
      year = data$year[1L],
      rate = sum(weights * data$count / data$population) * 100000,
      se = sqrt(
        sum(weights^2 * data$count / data$population^2)
      ) *
        100000
    )
  })
)

age_reference_file <- file.path(
  "inst",
  "extdata",
  "joinpoint-5.3.0",
  "Sample.Age.Adjusted.Rate.Calculation.Output.jpo"
)
age_raw <- readBin(
  age_reference_file,
  what = "raw",
  n = file.info(age_reference_file)$size
)

read_age_double <- function(offset) {
  readBin(
    age_raw[(offset + 1L):(offset + 8L)],
    what = "double",
    n = 1L,
    size = 8L,
    endian = "little"
  )
}

find_age_double <- function(target, tolerance = 1e-8) {
  offsets <- vapply(
    seq_len(length(age_raw) - 7L),
    function(offset) {
      value <- read_age_double(offset - 1L)
      isTRUE(is.finite(value) && abs(value - target) < tolerance)
    },
    logical(1)
  )
  which(offsets) - 1L
}

extract_age_group <- function(sex) {
  first_row <- which(age_rates$sex == sex & age_rates$year == 2000)[1L]
  offsets <- find_age_double(age_rates$rate[first_row])
  starts <- offsets - 8L
  starts <- starts[vapply(
    starts,
    function(start) {
      start >= 0L &&
        abs(read_age_double(start) - 2000) < 1e-8 &&
        abs(read_age_double(start + 22L * 32L) - 2022) < 1e-8
    },
    logical(1)
  )]

  if (length(starts) != 5L) {
    stop(
      "Expected five candidate age-adjusted models for sex ",
      sex,
      "; found ",
      length(starts),
      call. = FALSE
    )
  }

  do.call(
    rbind,
    lapply(seq_along(starts), function(model) {
      start <- starts[[model]]
      data.frame(
        model = model - 1L,
        sex = sex,
        year = vapply(
          0:22,
          function(index) read_age_double(start + index * 32L),
          numeric(1)
        ),
        rate = vapply(
          0:22,
          function(index) read_age_double(start + index * 32L + 8L),
          numeric(1)
        ),
        fitted = vapply(
          0:22,
          function(index) read_age_double(start + index * 32L + 16L),
          numeric(1)
        ),
        se = vapply(
          0:22,
          function(index) read_age_double(start + index * 32L + 24L),
          numeric(1)
        )
      )
    })
  )
}

v5_3_age_adjusted <- do.call(
  rbind,
  lapply(c(0L, 1L, 2L), extract_age_group)
)

jump_reference_file <- file.path(
  "inst",
  "extdata",
  "joinpoint-5.3.0",
  "Sample.Jump.CR.Standard.Output.jpo"
)
jump_raw <- readBin(
  jump_reference_file,
  what = "raw",
  n = file.info(jump_reference_file)$size
)

read_jump_double <- function(offset) {
  readBin(
    jump_raw[(offset + 1L):(offset + 8L)],
    what = "double",
    n = 1L,
    size = 8L,
    endian = "little"
  )
}

find_jump_double <- function(target, tolerance = 1e-8) {
  offsets <- vapply(
    seq_len(length(jump_raw) - 7L),
    function(offset) {
      value <- read_jump_double(offset - 1L)
      isTRUE(is.finite(value) && abs(value - target) < tolerance)
    },
    logical(1)
  )
  which(offsets) - 1L
}

jump_starts <- find_jump_double(2.705446699) - 8L
jump_starts <- jump_starts[vapply(
  jump_starts,
  function(start) {
    start >= 0L &&
      abs(read_jump_double(start) - 1992) < 1e-8 &&
      abs(read_jump_double(start + 30L * 32L) - 2022) < 1e-8
  },
  logical(1)
)]

if (length(jump_starts) != 6L) {
  stop(
    "Expected six Standard candidate models in the V5.3 Jump output; found ",
    length(jump_starts),
    call. = FALSE
  )
}

v5_3_jump_standard <- do.call(
  rbind,
  lapply(seq_along(jump_starts), function(model) {
    start <- jump_starts[[model]]
    data.frame(
      model = model - 1L,
      year = vapply(
        0:30,
        function(index) read_jump_double(start + index * 32L),
        numeric(1)
      ),
      rate = vapply(
        0:30,
        function(index) read_jump_double(start + index * 32L + 8L),
        numeric(1)
      ),
      fitted = vapply(
        0:30,
        function(index) read_jump_double(start + index * 32L + 16L),
        numeric(1)
      ),
      se = vapply(
        0:30,
        function(index) read_jump_double(start + index * 32L + 24L),
        numeric(1)
      )
    )
  })
)

dir.create(
  file.path("tests", "testthat", "fixtures"),
  recursive = TRUE,
  showWarnings = FALSE
)
save(
  v5_3_fitted,
  v5_3_age_adjusted,
  v5_3_jump_standard,
  file = file.path(
    "tests",
    "testthat",
    "fixtures",
    "joinpoint-v5.3-fitted.rda"
  ),
  version = 3
)
