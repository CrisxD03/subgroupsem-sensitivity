d <- simulate_intro(2000, 1, seed = 21)
y <- paste0("y", 1:5)
sg_A3 <- d$A == "3"

single_logl <- function(model, data)
  lavaan::lavInspect(lavaan::sem(model, data = data, meanstructure = TRUE, se = "none"),
                     "fit")[["logl"]]

test_that("LR equals baseline -2logL minus the sum of two separate group fits", {
  lr <- make_im_lr(model_intro_lr(), d)(sg_A3, d)
  expect_identical(lr$status, "ok")
  by_hand <- -2 * single_logl(model_intro_lr(), d) -
    -2 * (single_logl(model_intro_lr(), d[sg_A3, ]) + single_logl(model_intro_lr(), d[!sg_A3, ]))
  expect_equal(lr$value, by_hand, tolerance = 1e-6)
})

test_that("LR does not depend on the identification (lambda_1 = 1 vs Var(eta) = 1)", {
  marker <- "eta =~ y1 + y2 + y3 + y4 + y5"
  expect_equal(make_im_lr(model_intro_lr(), d)(sg_A3, d)$value,
               make_im_lr(marker, d)(sg_A3, d)$value, tolerance = 1e-6)
})

test_that("Wald with one constraint equals diff^2 / var(diff) from vcov", {
  w <- make_im_wald(model_intro_labelled(), "l2_0 == l2_1")(sg_A3, d)
  fit <- fit_two_group(model_intro_labelled(), d, sg_A3, se = "standard")
  est <- lavaan::coef(fit); V <- lavaan::vcov(fit)
  by_hand <- (est[["l2_0"]] - est[["l2_1"]])^2 /
    (V["l2_0", "l2_0"] + V["l2_1", "l2_1"] - 2 * V["l2_0", "l2_1"])
  expect_equal(w$value, by_hand, tolerance = 1e-8)
})

test_that("user-defined IM with the subgroup count equals the paper's exp(gw1) route", {
  ud <- make_im_ud(model_intro_labelled(), loading_pairs_intro())(sg_A3, d)
  expect_identical(ud$status, "ok")
  # Paper's route (p. 34-35): estimate the group weight and take exp(gw1).
  paper_model <- paste0(model_intro_labelled(), "\ngroup % c(gw0, gw1)*w")
  fit <- fit_two_group(paper_model, d, sg_A3, se = "standard",
                       extra_options = list(group.w.free = TRUE))
  est <- lavaan::coef(fit)
  expect_equal(exp(est[["gw1"]]), sum(sg_A3), tolerance = 1e-4)
  paper_value <- sqrt(exp(est[["gw1"]])) * sum(abs(est[paste0("l", 1:5, "_0")] - est[paste0("l", 1:5, "_1")]))
  expect_equal(ud$value, unname(paper_value), tolerance = 1e-4)
})

test_that("a subgroup too small to fit returns a marked failure, not a number", {
  tiny <- seq_len(nrow(d)) <= 4
  res <- make_im_lr(model_intro_lr(), d)(tiny, d)
  expect_true(is.na(res$value))
  expect_false(identical(res$status, "ok"))
})

test_that("a subgroup holding everyone returns 0, as in the package", {
  res <- make_im_lr(model_intro_lr(), d)(rep(TRUE, nrow(d)), d)
  expect_equal(res, list(value = 0, status = "all_in_subgroup"))
})
