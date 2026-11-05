# Sensitivity study: runs SubgroupSEM (three IMs) and the SEM tree on simulated
# datasets for every condition s x N, R replications each.
#
# Run from the project root:
#   Rscript scripts/20_run_grid.R            full study
#   Rscript scripts/20_run_grid.R timing     one dataset per N at s = 1, for timing
# One .rds file per dataset is written to results/sensitivity/runs/, so the run can
# be stopped and resumed; finished datasets are skipped. A dataset that errors is
# written to results/sensitivity/errors/ and reported by scripts/30_summarise.R.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

settings <- list(
  s_levels = c(0, 0.25, 0.5, 1),
  N_levels = c(500, 2500, 5000, 10000),
  R = 20,
  base_seed = 2022,
  min_size_fraction = 0.02,      # N_min = 0.02 N, see report Sect. 3.3
  search = list(algorithm = "SimpleDFS", depth = 2, k = 10),
  top_n = 3,                     # a planted group is recovered if it is in the top 3
  tree_jaccard = 0.8,            # tree: best leaf must match a planted group this well
  workers = 6
)

# Seed depends on N and r only: for a given replication the four effect sizes use
# the same draws of A, B, eta and residuals and differ only in the parameters.
design_grid <- function(settings) {
  g <- expand.grid(r = seq_len(settings$R), s = settings$s_levels, N = settings$N_levels)
  g$seed <- settings$base_seed + 1000L * match(g$N, settings$N_levels) + g$r
  g$id <- sprintf("N%05d_s%.2f_r%02d", g$N, g$s, g$r)
  g[order(-g$N, g$r, g$s), ]   # large N first: better load balance
}

# run_dataset() and measures_for() live in R/study.R so the parallel workers,
# which source R/, have them.
run_and_save <- function(task, settings, out_dir) {
  file <- file.path(out_dir, "runs", paste0(task$id, ".rds"))
  if (file.exists(file)) return("skipped")
  result <- tryCatch(run_dataset(task, settings), error = function(e) e)
  if (inherits(result, "error")) {
    saveRDS(list(task = task, error = conditionMessage(result)),
            file.path(out_dir, "errors", paste0(task$id, ".rds")))
    return("error")
  }
  saveRDS(result, file)
  "done"
}

main <- function(args) {
  timing <- length(args) > 0 && args[[1]] == "timing"
  out_dir <- file.path("results", if (timing) "timing" else "sensitivity")
  for (sub in c("runs", "errors")) dir.create(file.path(out_dir, sub), recursive = TRUE, showWarnings = FALSE)

  grid <- design_grid(settings)
  if (timing) grid <- grid[grid$s == 1 & grid$r == 1, ]
  jsonlite::write_json(settings, file.path(out_dir, "settings.json"), auto_unbox = TRUE, pretty = TRUE)
  write.csv(grid, file.path(out_dir, "design.csv"), row.names = FALSE)
  tasks <- split(grid, seq_len(nrow(grid)))

  if (timing) {
    init_python()
    status <- vapply(tasks, run_and_save, character(1), settings = settings, out_dir = out_dir)
  } else {
    cl <- parallel::makePSOCKcluster(settings$workers)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, {
      for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
      OpenMx::mxOption(NULL, "Number of Threads", 1)
      init_python()
    })
    status <- unlist(parallel::parLapplyLB(cl, tasks, run_and_save, settings = settings, out_dir = out_dir, chunk.size = 1))
  }
  print(table(status))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
}

main(commandArgs(trailingOnly = TRUE))
