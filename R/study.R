# One dataset of the sensitivity study: simulate, search with every IM, fit the
# SEM tree, score. Called by scripts/20_run_grid.R (also inside parallel workers).

measures_for <- function(data) {
  list(LR   = make_im_lr(model_intro_lr(), data),
       Wald = make_im_wald(model_intro_labelled(), wald_constraints_intro()),
       UD   = make_im_ud(model_intro_labelled(), loading_pairs_intro()))
}

run_dataset <- function(task, settings) {
  data <- simulate_intro(task$N, task$s, task$seed)
  search_data <- data[c(paste0("y", 1:5), "A", "B")]
  search_settings <- c(settings$search,
                       list(min_size = ceiling(settings$min_size_fraction * task$N)))

  searches <- lapply(measures_for(search_data), function(measure) {
    res <- run_search(search_data, measure, c("A", "B"), search_settings)
    list(top = res$top, log = res$log, seconds = res$seconds,
         recovered = recovery_sgsem(res$top$description, settings$top_n))
  })

  tree <- run_semtree(search_data)
  cells <- unique(data.frame(A = data$A, B = data$B, leaf = tree$leaf))
  tree_rec <- recovery_tree(cells, settings$tree_jaccard)

  list(task = task, cell_counts = table(A = data$A, B = data$B), searches = searches,
       tree = list(n_leaves = tree$n_leaves, leaf_by_cell = cells, seconds = tree$seconds,
                   recovered = tree_rec$recovered, best_jaccard = tree_rec$best_jaccard))
}
