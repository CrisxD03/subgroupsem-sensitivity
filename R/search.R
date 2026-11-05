# Subgroup search with the package's own engine (subgroupsem::subgroupsem, which
# calls pysubgroup), using one of our measures from R/measures.R.

# Must be called once per R process before any search.
init_python <- function() {
  reticulate::py_require(c("pysubgroup==0.7.8", "setuptools<70"), python_version = "3.11")
  ready <- subgroupsem::subgroupsem_ready()
  if (!isTRUE(ready)) stop("subgroupsem is not ready: ", ready)
  invisible(TRUE)
}

# Canonical description of a subgroup from its members: a predictor is part of the
# description if all members share one value of it. Used to label the log rows.
describe_members <- function(sg, dat, predictors) {
  members <- dat[!is.na(sg) & sg == 1, predictors, drop = FALSE]
  parts <- character(0)
  for (p in sort(predictors)) {
    vals <- unique(as.character(members[[p]]))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 1) parts <- c(parts, paste0(p, "=", vals))
  }
  paste(parts, collapse = " & ")
}

# settings: list(algorithm, depth, k, min_size)
run_search <- function(data, measure, predictors, settings) {
  stopifnot(all(predictors %in% names(data)),
            all(c("algorithm", "depth", "k", "min_size") %in% names(settings)))

  log_rows <- list()
  f_fit <- function(sg, dat) {
    res <- measure(sg, dat)
    log_rows[[length(log_rows) + 1L]] <<- data.frame(
      description = describe_members(sg, dat, predictors),
      n_sg = sum(sg == 1, na.rm = TRUE), value = res$value, status = res$status,
      stringsAsFactors = FALSE)
    # The search sees exactly what the package's wrappers would return.
    if (res$status %in% c("ok", "all_in_subgroup")) res$value else -1
  }

  t0 <- Sys.time()
  task <- subgroupsem::subgroupsem(
    f_fit = f_fit, dat = data, columns = predictors,
    algorithm = settings$algorithm, search_depth = as.integer(settings$depth),
    max_n_subgroups = as.integer(settings$k),
    min_subgroup_size = as.integer(settings$min_size))
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  top <- task@summary_statistics
  top$description <- vapply(top$subgroup, canonical_description, character(1))
  log <- if (length(log_rows)) do.call(rbind, log_rows) else
    data.frame(description = character(0), n_sg = integer(0), value = numeric(0),
               status = character(0))
  list(top = top, log = log, seconds = seconds)
}
