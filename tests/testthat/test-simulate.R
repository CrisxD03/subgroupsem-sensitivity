test_that("simulation is reproducible from its seed and depends on it", {
  a <- simulate_intro(500, 1, seed = 11)
  expect_identical(a, simulate_intro(500, 1, seed = 11))
  expect_false(identical(a$y1, simulate_intro(500, 1, seed = 12)$y1))
})

test_that("columns have the types of subsemdata01", {
  d <- simulate_intro(500, 1, seed = 1)
  expect_named(d, c(paste0("y", 1:5), "A", "B", "planted"))
  expect_true(is.factor(d$A) && !is.ordered(d$A))
  expect_identical(levels(d$A), as.character(1:5))
  expect_identical(d$planted, planted_group(d$A, d$B))
})

test_that("sample moments per group match the moments implied by Table 1", {
  # Implied: mean = nu, covariance = lambda lambda' + diag(theta), since Var(eta) = 1.
  d <- simulate_intro(250000, 1, seed = 3)
  p <- table1_params()
  y <- paste0("y", 1:5)
  for (g in names(p)) {
    sub <- as.matrix(d[d$planted == g, y])
    implied <- tcrossprod(p[[g]]$lambda) + diag(p[[g]]$theta)
    expect_equal(unname(colMeans(sub)), p[[g]]$nu, tolerance = 0.05)
    expect_equal(unname(cov(sub)), implied, tolerance = 0.05)
  }
})

test_that("at s = 0 all groups have the baseline distribution", {
  d <- simulate_intro(250000, 0, seed = 4)
  implied <- tcrossprod(table1_params()$baseline$lambda) + diag(5)
  sub <- as.matrix(d[d$planted == "G3", paste0("y", 1:5)])
  expect_equal(unname(cov(sub)), implied, tolerance = 0.08)
})

test_that("A and B are roughly uniform", {
  d <- simulate_intro(50000, 1, seed = 5)
  expect_true(all(abs(prop.table(table(d$A)) - 0.2) < 0.01))
  expect_true(all(abs(prop.table(table(d$B)) - 0.2) < 0.01))
})
