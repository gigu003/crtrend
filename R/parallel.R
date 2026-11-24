joinpoint_parallel_flag <- function(parallel) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  parallel
}

joinpoint_validate_cluster <- function(cluster) {
  if (is.null(cluster)) {
    return(NULL)
  }
  if (!inherits(cluster, "cluster")) {
    stop("`cluster` must be NULL or a parallel cluster object.", call. = FALSE)
  }
  cluster
}

joinpoint_worker_count <- function(workers, n_groups, use_parallel) {
  if (!use_parallel || n_groups <= 1L) {
    return(1L)
  }
  if (!is.null(workers) && !is_bounded_integer(workers, 1L)) {
    stop("`workers` must be NULL or a positive integer.", call. = FALSE)
  }
  if (is.null(workers)) {
    detected <- parallel::detectCores(logical = TRUE)
    if (is.na(detected)) {
      detected <- 2L
    }
    return(as.integer(min(n_groups, max(1L, detected - 1L))))
  }
  as.integer(min(workers, n_groups))
}

crtrend_initialize_cluster <- function(cluster) {
  namespace_path <- getNamespaceInfo(asNamespace("crtrend"), "path")
  source_files <- list.files(
    file.path(namespace_path, "R"),
    pattern = "[.][Rr]$",
    full.names = TRUE
  )
  source_signature <- if (length(source_files) > 0L) {
    paste(
      normalizePath(namespace_path, winslash = "/", mustWork = FALSE),
      max(file.info(source_files)$mtime),
      sep = "@"
    )
  } else {
    paste0("crtrend@", as.character(utils::packageVersion("crtrend")))
  }
  parallel::clusterCall(
    cluster,
    crtrend_worker_initialize,
    source_files = source_files,
    signature = source_signature
  )
  invisible(NULL)
}

crtrend_worker_initialize <- function(source_files, signature) {
  current <- getOption("crtrend.worker_state")
  if (is.list(current) && identical(current$signature, signature)) {
    return(invisible(NULL))
  }
  use_sources <- length(source_files) > 0L && all(file.exists(source_files))
  if (use_sources) {
    source_environment <- new.env(parent = globalenv())
    for (source_file in source_files) {
      sys.source(source_file, envir = source_environment)
    }
    mode <- "source"
  } else {
    loadNamespace("crtrend")
    source_environment <- NULL
    mode <- "namespace"
  }
  options(
    crtrend.worker_state = list(
      signature = signature,
      mode = mode,
      environment = source_environment
    )
  )
  invisible(NULL)
}
environment(crtrend_worker_initialize) <- baseenv()

crtrend_worker_invoke <- function(name, ...) {
  state <- getOption("crtrend.worker_state")
  target <- if (
    is.list(state) &&
      identical(state$mode, "source") &&
      is.environment(state$environment)
  ) {
    state$environment
  } else {
    asNamespace("crtrend")
  }
  get(name, envir = target, inherits = FALSE)(...)
}
environment(crtrend_worker_invoke) <- baseenv()
