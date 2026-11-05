# Interestingness measures (IMs) for the sensitivity study.
#
# Each make_im_*() returns a function(sg, dat) that fits the two-group SEM for one
# candidate subgroup and returns list(value, status). It mirrors the package's
# own wrappers in subgroupsem/R/subsem_userfun.R (commit fd24703):
#   subsem_lrt : fit with se = "none", check post.check,      LR = |-2 logL_2g - (-2 logL_base)|
#   subsem_wald: fit with default se,  check post.check,      W  = lavTestWald(fit, constraints)$stat
#   subsem     : fit with default se,  check post.check and converged, value from user formula
# The one difference is on purpose: the package turns every failure into -1, which
# looks the same as "not interesting". We return the failure with a status so it
# can be counted. R/search.R still hands -1 to the search, so the search itself
# behaves exactly like the package.

# Single-factor model, Var(eta) = 1, all loadings free (as in our replication).
model_intro_lr <- function() {
  "eta =~ NA*y1 + y2 + y3 + y4 + y5
   eta ~~ 1*eta
   eta ~ 0*1"
}

# Same model with labelled loadings and intercepts per group:
# suffix _0 = complement (group 1), _1 = subgroup (group 2), see fit_two_group().
# c(NA, NA)*y1 frees the first loading, as in the paper's Example 1 syntax (p. 34).
model_intro_labelled <- function() {
  loadings <- paste0("c(l", 1:5, "_0, l", 1:5, "_1)*y", 1:5, collapse = " + ")
  intercepts <- paste0("y", 1:5, " ~ c(n", 1:5, "_0, n", 1:5, "_1)*1", collapse = "\n")
  paste0("eta =~ ", loadings, " + c(NA, NA)*y1\n", intercepts,
         "\neta ~~ 1*eta\neta ~ 0*1")
}

# Wald constraints: equal loadings and intercepts in subgroup and complement
# (measurement invariance, Meredith 1993).
wald_constraints_intro <- function() {
  paste(c(paste0("l", 1:5, "_0 == l", 1:5, "_1"),
          paste0("n", 1:5, "_0 == n", 1:5, "_1")), collapse = "\n")
}

loading_pairs_intro <- function() {
  data.frame(complement = paste0("l", 1:5, "_0"), subgroup = paste0("l", 1:5, "_1"))
}

# ---------------------------------------------------------------------------

measure_ok <- function(value) list(value = value, status = "ok")
measure_failed <- function(status) list(value = NA_real_, status = status)

# Fit the two-group SEM exactly as the package does: subgroup indicator as group
# variable, rows sorted so that the complement (0) is group 1.
fit_two_group <- function(model, dat, sg, se, extra_options = list()) {
  dat$subgroup <- as.numeric(sg)
  dat <- dat[order(dat$subgroup), ]
  args <- c(list(model = model, data = dat, group = "subgroup", se = se, warn = FALSE),
            extra_options)
  do.call(lavaan::sem, args)
}

# Status of a fitted model, in the order the package checks it.
fit_status <- function(fit, require_converged) {
  if (!isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check"))))
    return("post_check_failed")
  if (require_converged && !isTRUE(lavaan::lavInspect(fit, "converged")))
    return("not_converged")
  "ok"
}

# Run one measure computation. lavaan errors become a marked failure carrying the
# error message, never a number.
evaluate_measure <- function(sg, compute) {
  if (all(sg == 1, na.rm = TRUE)) return(list(value = 0, status = "all_in_subgroup"))
  tryCatch(suppressWarnings(compute()),
           error = function(e) measure_failed(paste0("error: ", conditionMessage(e))))
}

make_im_lr <- function(model, data) {
  base <- suppressWarnings(lavaan::sem(model, data = data, se = "none",
                                       meanstructure = TRUE, warn = FALSE))
  if (!isTRUE(lavaan::lavInspect(base, "converged")))
    stop("Baseline model did not converge; the LR is undefined for this dataset.")
  base_m2ll <- -2 * lavaan::lavInspect(base, "fit")[["logl"]]

  function(sg, dat) {
    evaluate_measure(sg, function() {
      fit <- fit_two_group(model, dat, sg, se = "none")
      status <- fit_status(fit, require_converged = FALSE)
      if (status != "ok") return(measure_failed(status))
      measure_ok(abs(-2 * lavaan::lavInspect(fit, "fit")[["logl"]] - base_m2ll))
    })
  }
}

make_im_wald <- function(model, constraints) {
  function(sg, dat) {
    evaluate_measure(sg, function() {
      fit <- fit_two_group(model, dat, sg, se = "standard")
      status <- fit_status(fit, require_converged = FALSE)
      if (status != "ok") return(measure_failed(status))
      measure_ok(lavaan::lavTestWald(fit, constraints = constraints)$stat)
    })
  }
}

# User-defined IM: sqrt(n_sg) * sum_j |lambda_0j - lambda_1j|.
# The paper obtains n_sg as exp(gw1) from an estimated group weight
# (group.w.free = TRUE, p. 35). The ML estimate of that weight is the observed
# group count, so we use the count directly; tests/testthat/test-measures.R
# checks that both routes agree.
make_im_ud <- function(model, loading_pairs) {
  function(sg, dat) {
    evaluate_measure(sg, function() {
      fit <- fit_two_group(model, dat, sg, se = "standard")
      status <- fit_status(fit, require_converged = TRUE)
      if (status != "ok") return(measure_failed(status))
      est <- lavaan::coef(fit)
      stopifnot(all(c(loading_pairs$complement, loading_pairs$subgroup) %in% names(est)))
      diff <- abs(est[loading_pairs$complement] - est[loading_pairs$subgroup])
      measure_ok(sqrt(sum(sg == 1, na.rm = TRUE)) * sum(diff))
    })
  }
}
