# Generating parameters of the Introductory Example.
# Source: Kiefer et al. (2022), Table 1 (p. 8). G1 (A = 3) and G2 (B = 1) have
# identical parameters, so the cell where they overlap (A = 3, B = 1) is
# unambiguous and both are stored as "G12".

table1_params <- function() {
  list(
    baseline = list(lambda = c(1.0, 0.9, 0.8, 0.7, 0.6),
                    nu     = c(0.0, 0.1, 0.2, 0.3, 0.4),
                    theta  = rep(1.0, 5)),
    G12      = list(lambda = c(1.0, 1.1, 1.2, 1.3, 1.4),
                    nu     = c(0.0, 0.1, 0.2, 0.3, 0.4),
                    theta  = rep(1.0, 5)),
    G3       = list(lambda = c(1.0, 0.7, 0.6, 0.5, 0.4),
                    nu     = c(0.0, -0.1, -0.2, -0.3, -0.4),
                    theta  = rep(1.2, 5))
  )
}

# Shrink every planted deviation by factor s:
#   theta_g(s) = theta_base + s * (theta_g - theta_base)
# s = 1 is the paper's setting, s = 0 makes all groups equal to the baseline.
scale_params <- function(params, s) {
  stopifnot(is.numeric(s), length(s) == 1, s >= 0, "baseline" %in% names(params))
  base <- params$baseline
  lapply(params, function(group) {
    stopifnot(identical(names(group), names(base)))
    Map(function(b, g) b + s * (g - b), base, group)
  })
}

# Which planted group a pupil belongs to, from school A and grade B (1..5).
planted_group <- function(A, B) {
  A <- as.integer(as.character(A))
  B <- as.integer(as.character(B))
  stopifnot(!anyNA(A), !anyNA(B))
  ifelse(A == 5 & B == 5, "G3", ifelse(A == 3 | B == 1, "G12", "baseline"))
}

# Planted groups as subgroup descriptions, in the canonical form used by
# R/score.R ("attr=value", conjunctions joined by " & ", sorted by attribute).
planted_descriptions <- function() {
  c(G1 = "A=3", G2 = "B=1", G3 = "A=5 & B=5")
}

# Membership of each planted group as logical vectors.
planted_membership <- function(A, B) {
  A <- as.integer(as.character(A))
  B <- as.integer(as.character(B))
  list(G1 = A == 3, G2 = B == 1, G3 = A == 5 & B == 5)
}
