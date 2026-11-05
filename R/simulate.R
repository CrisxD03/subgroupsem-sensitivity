# Data-generating mechanism of the Introductory Example.
#
# What the paper states (pp. 6-8): N pupils, drawn uniformly from five schools (A)
# and five grades (B), a single-factor model with five indicators, parameters per
# group as in Table 1.
# What the paper does NOT state, and we therefore assume:
#   * A and B are drawn independently and uniformly (the bundled subsemdata01
#     has unbalanced cells, 175-232 per cell, consistent with random draws);
#   * eta ~ N(0, 1) in every group (Table 1 gives no factor variance; with
#     lambda_1 = 1 fixed, the estimated variance in subsemdata01 is close to 1);
#   * residuals are normal and independent with variances theta.
# scripts/10_check_inputs.R checks these assumptions against subsemdata01.

simulate_intro <- function(N, s, seed, params = table1_params()) {
  stopifnot(length(N) == 1, N == round(N), N >= 25,
            length(seed) == 1, is.finite(seed))
  p <- scale_params(params, s)

  withr::with_seed(seed, {
    A <- sample.int(5, N, replace = TRUE)
    B <- sample.int(5, N, replace = TRUE)
    eta <- stats::rnorm(N)
    eps <- matrix(stats::rnorm(N * 5), N, 5)
  })

  group <- planted_group(A, B)
  y <- matrix(NA_real_, N, 5, dimnames = list(NULL, paste0("y", 1:5)))
  for (g in names(p)) {
    idx <- group == g
    if (!any(idx)) next
    pg <- p[[g]]
    y[idx, ] <- sweep(outer(eta[idx], pg$lambda), 2, pg$nu, "+") +
      sweep(eps[idx, , drop = FALSE], 2, sqrt(pg$theta), "*")
  }
  stopifnot(!anyNA(y))

  # Same column types as subsemdata01: A and B are unordered factors "1".."5".
  data.frame(y, A = factor(A, levels = 1:5), B = factor(B, levels = 1:5),
             planted = group, stringsAsFactors = FALSE)
}
