test_that("pysubgroup descriptions are parsed into canonical form", {
  expect_identical(canonical_description("A=='5' AND B=='5'"), "A=5 & B=5")
  expect_identical(canonical_description("B=='1' AND A=='3'"), "A=3 & B=1")
  expect_identical(canonical_description("A=='3'"), "A=3")
  expect_identical(canonical_description("grade==7 AND school=='Grant-White' AND sex==2"),
                   "grade=7 & school=Grant-White & sex=2")
  expect_error(canonical_description("A > 3"))
})

test_that("recovery counts planted descriptions in the top 3 only", {
  top <- c("A=3", "A=5 & B=5", "A=3 & B=3", "B=1")
  expect_identical(recovery_sgsem(top), c(G1 = TRUE, G2 = FALSE, G3 = TRUE))
  expect_identical(recovery_sgsem(top, top_n = 4), c(G1 = TRUE, G2 = TRUE, G3 = TRUE))
})

grid <- expand.grid(A = 1:5, B = 1:5)

test_that("Jaccard index on known sets", {
  g2 <- grid$B == 1
  expect_equal(jaccard(g2, g2), 1)
  expect_equal(jaccard(g2, grid$B == 2), 0)
  # G2 without the cell it shares with G1: 4 of 5 cells
  expect_equal(jaccard(g2 & grid$A != 3, g2), 0.8)
})

test_that("tree recovery: one leaf per cell finds only G3", {
  leaf <- paste(grid$A, grid$B)
  r <- recovery_tree(data.frame(A = grid$A, B = grid$B, leaf = leaf))
  expect_identical(r$recovered, c(G1 = FALSE, G2 = FALSE, G3 = TRUE))
  expect_equal(unname(r$best_jaccard), c(0.2, 0.2, 1))
})

test_that("tree recovery: leaves G1, G2 minus overlap, G3, rest find all three", {
  leaf <- ifelse(grid$A == 3, "G1", ifelse(grid$B == 1, "G2", ifelse(grid$A == 5 & grid$B == 5, "G3", "rest")))
  r <- recovery_tree(data.frame(A = grid$A, B = grid$B, leaf = leaf))
  expect_identical(r$recovered, c(G1 = TRUE, G2 = TRUE, G3 = TRUE))
})

test_that("tree recovery counts cells, so cell sizes cannot push a best-possible leaf below 0.8", {
  # The input has one row per cell; a cell listed twice is an error, not extra weight.
  cells <- data.frame(A = grid$A, B = grid$B, leaf = ifelse(grid$B == 1 & grid$A != 3, "G2", "rest"))
  expect_equal(unname(recovery_tree(cells)$best_jaccard[["G2"]]), 0.8)
  expect_error(recovery_tree(rbind(cells, cells[1, ])))
})
