test_that("s = 1 reproduces Table 1 and s = 0 removes every planted deviation", {
  p <- table1_params()
  expect_equal(scale_params(p, 1), p, tolerance = 1e-12)
  p0 <- scale_params(p, 0)
  for (g in names(p0)) expect_equal(p0[[g]], p$baseline, tolerance = 1e-12)
})

test_that("s = 0.5 halves the deviation", {
  p <- scale_params(table1_params(), 0.5)
  expect_equal(p$G12$lambda[5], 0.6 + 0.5 * (1.4 - 0.6), tolerance = 1e-12)
  expect_equal(p$G3$nu[5], 0.4 + 0.5 * (-0.4 - 0.4), tolerance = 1e-12)
  expect_equal(p$G3$theta, rep(1.1, 5), tolerance = 1e-12)
})

test_that("negative or vector s is rejected", {
  expect_error(scale_params(table1_params(), -0.1))
  expect_error(scale_params(table1_params(), c(0.5, 1)))
})

test_that("planted groups follow the paper (G1: A = 3, G2: B = 1, G3: A = 5 and B = 5)", {
  expect_equal(planted_group(A = c(3, 1, 3, 5, 5, 2), B = c(2, 1, 1, 5, 4, 3)),
               c("G12", "G12", "G12", "G3", "baseline", "baseline"))
  m <- planted_membership(A = c(3, 1, 5), B = c(1, 1, 5))
  expect_equal(m$G1, c(TRUE, FALSE, FALSE))
  expect_equal(m$G2, c(TRUE, TRUE, FALSE))
  expect_equal(m$G3, c(FALSE, FALSE, TRUE))
})
