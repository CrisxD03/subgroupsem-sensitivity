# SEM tree on the Introductory Example, with the same model and settings as our
# replication of Figure 2 in the report: single-factor RAM model,
# Var(eta) = 1, free loadings, intercepts and residual variances, covariates A and
# B as unordered factors, semtree_control() defaults.

intro_openmx_model <- function(data) {
  y <- paste0("y", 1:5)
  OpenMx::mxModel(
    "cfa", type = "RAM", manifestVars = y, latentVars = "eta",
    OpenMx::mxPath("eta", y, free = TRUE, values = 0.8),
    OpenMx::mxPath("eta", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(y, arrows = 2, free = TRUE, values = 1),
    OpenMx::mxPath("one", y, free = TRUE, values = 0),
    OpenMx::mxData(data[y], type = "raw"))
}

# Returns the leaf (node id) of every row of `data`.
run_semtree <- function(data, control = semtree::semtree_control(verbose = FALSE)) {
  y <- paste0("y", 1:5)
  stopifnot(all(c(y, "A", "B") %in% names(data)), is.factor(data$A), is.factor(data$B))
  tree_data <- data[c(y, "A", "B")]
  model <- OpenMx::mxRun(intro_openmx_model(tree_data), silent = TRUE)

  t0 <- Sys.time()
  tree <- suppressMessages(semtree::semtree(model, data = tree_data, control = control))
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  leaf <- as.vector(stats::predict(tree, tree_data))
  if (anyNA(leaf)) stop("semtree returned no leaf for some rows")
  list(leaf = leaf, n_leaves = length(unique(leaf)), seconds = seconds)
}
