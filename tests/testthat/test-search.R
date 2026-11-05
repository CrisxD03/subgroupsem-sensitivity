# Needs Python with pysubgroup 0.7.8 (slow: two searches).
skip_if_not(python_available(), "pysubgroup not available")

d <- simulate_intro(2000, 1, seed = 31)
settings <- list(algorithm = "SimpleDFS", depth = 2, k = 10, min_size = 30)

test_that("our LR inside the package search gives the same top 10 as subsem_lrt", {
  ours <- run_search(d[c(paste0("y", 1:5), "A", "B")], make_im_lr(model_intro_lr(), d),
                     c("A", "B"), settings)
  pkg <- suppressWarnings(subgroupsem::subsem_lrt(
    model = model_intro_lr(), data = d[c(paste0("y", 1:5), "A", "B")], predictors = c("A", "B"),
    lavaan_options = list(warn = FALSE),
    subsem_options = list(algorithm = "SimpleDFS", search_depth = 2L,
                          min_subgroup_size = 30L, max_n_subgroups = 10L)))
  pkg_top <- pkg@summary_statistics
  expect_identical(ours$top$description,
                   vapply(pkg_top$subgroup, canonical_description, character(1), USE.NAMES = FALSE))
  expect_equal(ours$top$quality, pkg_top$quality, tolerance = 1e-6)
})

test_that("the log labels subgroups the same way pysubgroup does", {
  res <- run_search(d[c(paste0("y", 1:5), "A", "B")], make_im_lr(model_intro_lr(), d),
                    c("A", "B"), settings)
  logged <- res$log[res$log$status == "ok", ]
  idx <- match(res$top$description, logged$description)
  expect_false(anyNA(idx))
  expect_equal(logged$value[idx], res$top$quality, tolerance = 1e-8)
})
