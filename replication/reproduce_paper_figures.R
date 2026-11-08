# =============================================================================
# Reproduce the paper's figures and tables for the datasets that are available
# (Kiefer, Lemmerich, Langenberg & Mayer, "Subgroup Discovery in SEM")
#
#   Figure 2  Simulated example: Original / Mixture model / SEM tree / SubgroupSEM
#   Table 1   (check) group parameters estimated from the authors' simulated data
#   Table 5   Example 1 (Holzinger-Swineford): top-5 subgroups, IM1 and IM2
#   Table 6   Example 1: loadings and intercepts per subgroup vs baseline
#
# NOT reproducible: Example 2 (Figure 5, Tables 7-9) needs restricted NEPS data.
# Tables 2-4 and Figures 1, 3, 4 are conceptual (no data).
#
# HOW TO RUN: restart R, run top to bottom. Takes ~5 minutes
# (SEM tree ~1 min, mixture model ~30 s, the rest is fast).
# Output: results/paper_figures/  (PDF + PNG figures, CSV + LaTeX tables)
# =============================================================================
#Setup 
for (p in c("lavaan", "remotes", "ggplot2", "patchwork", "OpenMx", "semtree"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("subgroupsem", quietly = TRUE))
  remotes::install_github("chkiefer/subgroupsem@fd24703", upgrade = "never")   # pinned commit
t_start <- Sys.time()
timing <- list()

library(reticulate)
if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  # A Python 3.11 env with pysubgroup==0.7.8, setuptools<70, numpy<2, pandas<2.2
  use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
} else {
  if (packageVersion("reticulate") < "1.41") {
    install.packages("reticulate"); stop("reticulate updated: restart R and rerun.")
  }
  py_require(c("pysubgroup==0.7.8", "setuptools<70"), python_version = "3.11")
}
suppressPackageStartupMessages({
  library(lavaan); library(subgroupsem); library(ggplot2); library(patchwork)
  library(OpenMx); library(semtree)
})
stopifnot(isTRUE(subgroupsem_ready()))
out <- "results/paper_figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(2022)




## PART A: Example 1  Holzinger & Swineford (Tables 5 and 6)

data("HolzingerSwineford1939", package = "lavaan")
hs <- HolzingerSwineford1939
covariates <- c("school", "grade", "sex")
dfs <- list(algorithm = "SimpleDFS", search_depth = 3L, min_subgroup_size = 30L,
            max_n_subgroups = 10L)

model_im1 <- "
visual  =~ NA*x1 + x2 + x3
textual =~ NA*x4 + x5 + x6
speed   =~ NA*x7 + x8 + x9
visual ~~ 1*visual
textual ~~ 1*textual
speed ~~ 1*speed
visual + textual + speed ~ 0*1
"
model_im2 <- "
visual  =~ c(la01,la11)*x1 + c(NA,NA)*x1 + c(la02,la12)*x2 + c(la03,la13)*x3
textual =~ NA*x4 + x5 + x6
speed   =~ NA*x7 + x8 + x9
visual ~~ 1*visual
textual ~~ 1*textual
speed ~~ 1*speed
visual + textual + speed ~ 0*1
group % c(gw0,gw1)*w
"
qf_im2 <- "
N1 := exp(gw1)
diff1 := abs(la01 - la11)
diff2 := abs(la02 - la12)
diff3 := abs(la03 - la13)
subsem_qf := sqrt(N1) * (diff1 + diff2 + diff3)
"

t0 <- Sys.time()
fit_im1 <- suppressWarnings(subsem_lrt(model = model_im1, data = hs, predictors = covariates,
                                       lavaan_options = list(warn = FALSE), subsem_options = dfs))
timing$example1_IM1 <- Sys.time() - t0; t0 <- Sys.time()
fit_im2 <- suppressWarnings(subsem(model = model_im2, data = hs, qf = qf_im2, predictors = covariates,
                                   lavaan_options = list(warn = FALSE, group.w.free = TRUE),   # as on p. 35
                                   subsem_options = dfs))
timing$example1_IM2 <- Sys.time() - t0

# Readable descriptions. lavaan codes sex 1 = male, 2 = female; the paper labels
# code 2 as "Male" in its Table 5. We use the correct labels.
pretty <- function(s) {
  s <- gsub("school=='Grant-White'", "Grant-White", s)
  s <- gsub("school=='Pasteur'", "Pasteur", s)
  s <- gsub("grade==([0-9])", "Grade \\1", s)
  s <- gsub("sex==2", "Female", s); s <- gsub("sex==1", "Male", s)
  gsub(" AND ", " & ", s)
}

#Table 5

paper_t5 <- data.frame(
  IM1_paper_subgroup = c("Pasteur", "Grant-White", "Grant-White & Grade 7",
                         "Grant-White & Female & Grade 7", "Grade 7"),   # paper: "Male"
  IM1_paper_N  = c(156, 145, 79, 43, 157),
  IM1_paper    = c(111.095, 111.095, 107.820, 99.360, 97.032),
  IM2_paper_subgroup = c("Pasteur & Female", "Pasteur", "Grant-White", "Female", "Male"),  # paper: swapped
  IM2_paper_N  = c(82, 156, 145, 155, 146),
  IM2_paper    = c(7.025, 6.890, 6.643, 5.533, 5.370))

s1 <- fit_im1@summary_statistics; s2 <- fit_im2@summary_statistics
t5 <- data.frame(
  Rank = 1:5,
  `IM1 subgroup (ours)` = pretty(s1$subgroup[1:5]), `N` = s1$size_sg[1:5],
  `LR ours` = sprintf("%.3f", s1$quality[1:5]), `LR paper` = sprintf("%.3f", paper_t5$IM1_paper),
  `IM2 subgroup (ours)` = pretty(s2$subgroup[1:5]), `N ` = s2$size_sg[1:5],
  `Delta LR` = sprintf("%+.3f", round(s1$quality[1:5], 3) - paper_t5$IM1_paper),
  `IM2 ours` = sprintf("%.3f", s2$quality[1:5]), `IM2 paper` = sprintf("%.3f", paper_t5$IM2_paper),
  `Delta IM2` = sprintf("%+.3f", round(s2$quality[1:5], 3) - paper_t5$IM2_paper),
  check.names = FALSE)
write.csv(t5, file.path(out, "table5_example1.csv"), row.names = FALSE)

print(t5)

# Table 6: CFA parameters within each IM1 subgroup 
# Paper's Table 6 order. 
sg_def <- list(
  "Subgroup 1" = hs$school == "Grant-White",
  "Subgroup 2" = hs$school == "Pasteur",
  "Subgroup 3" = hs$school == "Grant-White" & hs$grade == 7,
  "Subgroup 4" = hs$school == "Grant-White" & hs$grade == 7 & hs$sex == 2,
  "Subgroup 5" = !is.na(hs$grade) & hs$grade == 7)
sg_def <- lapply(sg_def, function(x) !is.na(x) & x)

cfa_pars <- function(dat) {
  f <- cfa(model_im1, data = dat, meanstructure = TRUE)
  pe <- parameterEstimates(f)
  c(setNames(pe$est[pe$op == "~1" & pe$lhs %in% paste0("x", 1:9)], paste0("nu", 1:9)),
    setNames(pe$est[pe$op == "=~"], paste0("lambda", 1:9)))
}
ours_t6 <- cbind(Baseline = cfa_pars(hs), sapply(sg_def, function(idx) cfa_pars(hs[idx, ])))

paper_t6 <- rbind(
  nu1 = c(4.936, 4.930, 4.941, 4.772, 4.810, 4.746), nu2 = c(6.088, 6.200, 5.984, 6.054, 6.070, 5.932),
  nu3 = c(2.250, 1.996, 2.487, 1.862, 1.805, 2.097), nu4 = c(3.061, 3.317, 2.823, 3.148, 3.372, 2.828),
  nu5 = c(4.341, 4.712, 3.995, 4.446, 4.605, 4.123), nu6 = c(2.186, 2.469, 1.922, 2.329, 2.399, 2.020),
  nu7 = c(4.186, 3.921, 4.432, 3.542, 3.534, 3.823), nu8 = c(5.527, 5.488, 5.563, 5.223, 5.167, 5.238),
  nu9 = c(5.374, 5.327, 5.418, 5.057, 5.121, 5.166),
  lambda1 = c(0.900, 0.777, 1.047, 0.887, 0.954, 0.909), lambda2 = c(0.498, 0.572, 0.412, 0.576, 0.725, 0.462),
  lambda3 = c(0.656, 0.719, 0.597, 0.636, 0.683, 0.662), lambda4 = c(0.990, 0.971, 0.946, 0.964, 1.137, 0.943),
  lambda5 = c(1.102, 0.961, 1.119, 1.075, 1.189, 1.078), lambda6 = c(0.917, 0.935, 0.827, 1.010, 1.075, 0.900),
  lambda7 = c(0.619, 0.679, 0.591, 0.362, 0.462, 0.448), lambda8 = c(0.731, 0.833, 0.665, 0.592, 0.617, 0.513),
  lambda9 = c(0.670, 0.719, 0.545, 0.992, 1.074, 0.855))
colnames(paper_t6) <- colnames(ours_t6)

diff_t6 <- ours_t6 - paper_t6
cat("\nTable 6: largest absolute difference ours vs paper:", round(max(abs(diff_t6)), 3), "\n")
print(round(diff_t6, 3))

t6 <- data.frame(Parameter = sub("lambda([0-9])", "$\\\\lambda_\\1$", sub("^nu([0-9])", "$\\\\nu_\\1$", rownames(ours_t6))),
                 apply(ours_t6, 2, function(x) sprintf("%.3f", x)), check.names = FALSE)
write.csv(data.frame(parameter = rownames(ours_t6), round(ours_t6, 3), check.names = FALSE),
          file.path(out, "table6_example1_ours.csv"), row.names = FALSE)
write.csv(data.frame(parameter = rownames(diff_t6), round(diff_t6, 3), check.names = FALSE),
          file.path(out, "table6_example1_diff_vs_paper.csv"), row.names = FALSE)



# Missing grade and inadmissible solutions 
# subgroupsem sets a case to NA in the two-group model when a selector variable of the
# description is missing; the baseline keeps all 301 pupils. Compare with keeping the
# pupil (id 351, Grant-White) in the complement, and check lavaan's post.check, which
# the package uses to score inadmissible fits as -1.
base_hs <- sem(model_im1, data = hs, se = "none", meanstructure = TRUE)
m2ll_base <- -2 * lavInspect(base_hs, "fit")[["logl"]]
lr_of <- function(sg) {
  dd <- hs; dd$subgroup <- as.numeric(sg); dd <- dd[order(dd$subgroup), ]
  f <- sem(model_im1, data = dd, group = "subgroup", se = "none", warn = FALSE)
  c(LR = abs(-2 * lavInspect(f, "fit")[["logl"]] - m2ll_base),
    admissible = suppressWarnings(lavInspect(f, "post.check")), n_used = sum(!is.na(sg)))
}
na_g <- is.na(hs$grade); g7 <- !na_g & hs$grade == 7
cands <- list("Pasteur" = list(hs$school == "Pasteur", FALSE),
              "Grant-White & Female" = list(hs$sex == 2 & hs$school == "Grant-White", FALSE),
              "Grant-White & Grade 7" = list(hs$school == "Grant-White" & g7, TRUE),
              "Pasteur & Grade 7" = list(hs$school == "Pasteur" & g7, TRUE),
              "Grant-White & Female & Grade 7" = list(hs$school == "Grant-White" & g7 & hs$sex == 2, TRUE),
              "Grade 7" = list(g7, TRUE))
missing_grade <- do.call(rbind, lapply(names(cands), function(k) {
  sg <- cands[[k]][[1]]; a <- lr_of(ifelse(cands[[k]][[2]] & na_g, NA, sg)); b <- lr_of(sg)
  data.frame(subgroup = k, N = sum(sg), LR_as_package = round(a[["LR"]], 3), n_used = a[["n_used"]],
             LR_same_cases = round(b[["LR"]], 3), delta = round(b[["LR"]] - a[["LR"]], 3),
             admissible = as.logical(a[["admissible"]]))
}))
write.csv(missing_grade, file.path(out, "example1_missing_grade.csv"), row.names = FALSE)
print(missing_grade, row.names = FALSE)
# Full ranking with the same 301 cases in both models (pupil kept in the complement)
sel <- list(sex = c(1, 2), school = c("Pasteur", "Grant-White"), grade = c(7, 8))
same_rank <- do.call(rbind, unlist(lapply(1:3, function(k)
  lapply(combn(names(sel), k, simplify = FALSE), function(vs) {
    grid <- expand.grid(sel[vs], stringsAsFactors = FALSE)
    do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      sg <- Reduce(`&`, lapply(vs, function(v) !is.na(hs[[v]]) & hs[[v]] == grid[i, v]))
      if (sum(sg) < 30) return(NULL)
      r <- lr_of(sg)
      data.frame(subgroup = paste(vs, grid[i, ], sep = "=", collapse = " & "), N = sum(sg),
                 LR_same_cases = round(r[["LR"]], 3), admissible = as.logical(r[["admissible"]]))
    }))
  })), recursive = FALSE))
same_rank <- same_rank[order(-same_rank$LR_same_cases), ]
write.csv(same_rank, file.path(out, "example1_same_cases_ranking.csv"), row.names = FALSE)


# PART B: Simulated example (Figure 2, Table 1 check)

data(subsemdata01)
sim <- subsemdata01
y <- paste0("y", 1:5)
sim$planted <- ifelse(sim$A == "5" & sim$B == "5", "G3",
                      ifelse(sim$A == "3" | sim$B == "1", "G1/G2", "Baseline"))

# Table 1 check: estimate each group's parameters

m1 <- "eta =~ NA*y1 + y2 + y3 + y4 + y5
eta ~~ 1*eta"
t1_fit <- lapply(c("Baseline", "G1/G2", "G3"), function(g) {
  pe <- parameterEstimates(sem(m1, data = sim[sim$planted == g, ], meanstructure = TRUE))
  keep <- (pe$op == "=~") | (pe$op == "~1" & pe$lhs %in% y) | (pe$op == "~~" & pe$lhs %in% y & pe$lhs == pe$rhs)
  pe <- pe[keep, ]; pe[order(match(pe$op, c("=~", "~1", "~~"))), c("est", "se")]
})
specified <- cbind(Baseline = c(1, .9, .8, .7, .6, 0, .1, .2, .3, .4, rep(1, 5)),
                   `G1/G2`  = c(1, 1.1, 1.2, 1.3, 1.4, 0, .1, .2, .3, .4, rep(1, 5)),
                   G3       = c(1, .7, .6, .5, .4, 0, -.1, -.2, -.3, -.4, rep(1.2, 5)))
pnames <- c(paste0("lambda", 1:5), paste0("nu", 1:5), paste0("theta", 1:5))
t1_out <- data.frame(parameter = pnames, check.names = FALSE)
for (k in 1:3) {
  g <- colnames(specified)[k]; e <- t1_fit[[k]]
  t1_out[[paste(g, "spec")]] <- specified[, k]
  t1_out[[paste(g, "est")]]  <- round(e$est, 2)
  t1_out[[paste(g, "z")]]    <- round((e$est - specified[, k]) / e$se, 1)
}
t1_out$n_Baseline <- sum(sim$planted == "Baseline"); t1_out$n_G1G2 <- sum(sim$planted == "G1/G2")
t1_out$n_G3 <- sum(sim$planted == "G3")
write.csv(t1_out, file.path(out, "table1_check_simulated.csv"), row.names = FALSE)
print(t1_out)

# Figure 2, panel 4: SubgroupSEM 
model_sim <- "
eta =~ NA*y1 + y2 + y3 + y4 + y5
eta ~~ 1*eta
eta ~ 0*1
"
t0 <- Sys.time()
fit_sim <- suppressWarnings(subsem_lrt(model = model_sim, data = sim, predictors = c("A", "B"),
                                       lavaan_options = list(warn = FALSE),
                                       subsem_options = list(algorithm = "SimpleDFS", search_depth = 2L, min_subgroup_size = 30L,
                                                             max_n_subgroups = 10L)))
timing$intro_SubgroupSEM <- Sys.time() - t0
sg_top5 <- fit_sim@summary_statistics[1:5, ]
print(sg_top5)
write.csv(fit_sim@summary_statistics, file.path(out, "figure2_subgroupsem_top10.csv"), row.names = FALSE)

# Figure 2, panel 2: two-class factor mixture model (OpenMx)
mix_class <- function(k, lam) mxModel(paste0("class", k),
                                      mxMatrix("Full", 5, 1, TRUE, lam, name = "L", labels = paste0("l", k, 1:5)),
                                      mxMatrix("Diag", 5, 5, TRUE, .8, name = "Th", lbound = 1e-3, labels = paste0("t", k, 1:5)),
                                      mxMatrix("Full", 1, 5, TRUE, colMeans(sim[y]), name = "M", labels = paste0("n", k, 1:5)),
                                      mxAlgebra(L %*% t(L) + Th, name = "S"),
                                      mxExpectationNormal("S", "M", dimnames = y), mxFitFunctionML(vector = TRUE))
build_mix <- function(lam1, lam2) mxModel("mix", mix_class(1, lam1), mix_class(2, lam2),
                                          mxData(sim[y], "raw"),
                                          mxMatrix("Full", 2, 1, c(FALSE, TRUE), c(1, 1), lbound = 1e-3, name = "w"),
                                          mxExpectationMixture(c("class1", "class2"), weights = "w", scale = "sum"),
                                          mxFitFunctionML())
class_prob <- function(f) {
  w <- as.vector(mxEval(w, f)); w <- w / sum(w)
  dens <- function(k) {
    S <- mxEvalByName(paste0("class", k, ".S"), f, compute = TRUE)
    M <- as.vector(mxEvalByName(paste0("class", k, ".M"), f, compute = TRUE))
    X <- sweep(as.matrix(sim[y]), 2, M)
    exp(-0.5 * rowSums((X %*% solve(S)) * X)) / sqrt((2 * pi)^5 * det(S))
  }
  p1 <- w[1] * dens(1) / (w[1] * dens(1) + w[2] * dens(2))
  # class 1 = the class resembling G1/G2 (larger loading on y5); relabel if needed
  if (mxEval(class1.L, f)[5] < mxEval(class2.L, f)[5]) { p1 <- 1 - p1; w <- rev(w) }
  list(p1 = p1, w = w)
}
# Robustness to local maxima (paper p. 11): 10 seeds, random starting loadings.
t0 <- Sys.time()
mix_runs <- do.call(rbind, lapply(2022:2031, function(seed) {
  set.seed(seed)
  f <- mxTryHard(build_mix(runif(5, .5, 1.5), runif(5, .5, 1.5)), extraTries = 10, silent = TRUE)
  cp <- class_prob(f); m <- tapply(cp$p1, sim$planted, mean)
  data.frame(seed = seed, minus2LL = f$output$minimum, weight_class1 = cp$w[1],
             p_Baseline = m[["Baseline"]], p_G1G2 = m[["G1/G2"]], p_G3 = m[["G3"]])
}))
timing$mixture_10_seeds <- Sys.time() - t0
write.csv(mix_runs, file.path(out, "figure2_mixture_seeds.csv"), row.names = FALSE)
print(mix_runs, row.names = FALSE)
# Fit used for the figure: start at the generating values, seed 2022
set.seed(2022)
fit_mix <- mxTryHard(build_mix(c(1, 1.1, 1.2, 1.3, 1.4), c(1, .9, .8, .7, .6)), extraTries = 10, silent = TRUE)
cp <- class_prob(fit_mix); w <- cp$w; sim$p_class1 <- cp$p1
cat("\nMixture: class weights", round(w, 3), "\n")
print(round(tapply(sim$p_class1, sim$planted, mean), 3))
print(table(assigned_class1 = sim$p_class1 > .5, sim$planted))

# Figure 2, panel 3: SEM tree (semtree) 
cfa_mx <- mxModel("cfa", type = "RAM", manifestVars = y, latentVars = "eta",
                  mxPath("eta", y, free = TRUE, values = .8),
                  mxPath("eta", arrows = 2, free = FALSE, values = 1),
                  mxPath(y, arrows = 2, free = TRUE, values = 1),
                  mxPath("one", y, free = TRUE, values = 0),
                  mxData(sim[y], "raw"))
cfa_mx <- mxRun(cfa_mx, silent = TRUE)
t0 <- Sys.time()
tree <- semtree(cfa_mx, data = sim[c(y, "A", "B")], control = semtree_control(verbose = FALSE))
print(tree)
sim$leaf <- predict(tree, sim[c(y, "A", "B")])
timing$semtree <- Sys.time() - t0
capture.output(print(tree), file = file.path(out, "figure2_semtree.txt"))
leaf_table <- unique(data.frame(A = sim$A, B = sim$B, leaf = sim$leaf))
write.csv(leaf_table[order(leaf_table$A, leaf_table$B), ], file.path(out, "figure2_semtree_leaf_by_cell.csv"), row.names = FALSE)

# Figure 2: draw the four panels 
# Covariate space: B (grade) on the x-axis, A (school) on the y-axis, A = 1 at the top

ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e0"
cat_pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
cell <- function(a, b) data.frame(xmin = b - .5, xmax = b + .5, ymin = a - .5, ymax = a + .5)
theme_map <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = ink2),
        axis.title = element_text(colour = ink2), legend.position = "none",
        plot.title = element_text(face = "bold", hjust = .5, colour = ink),
        panel.border = element_rect(fill = NA, colour = ink2, linewidth = .4))
base_map <- function(title) ggplot() +
  scale_x_continuous("B (grade)", breaks = 1:5, limits = c(.5, 5.5), expand = c(0, 0)) +
  scale_y_reverse("A (school)", breaks = 1:5, limits = c(5.5, .5), expand = c(0, 0)) +
  coord_fixed() + ggtitle(title) + theme_map +
  geom_hline(yintercept = 1.5:4.5, colour = grid_col, linewidth = .3) +
  geom_vline(xintercept = 1.5:4.5, colour = grid_col, linewidth = .3)

# Panel 1: original (planted) groups
orig_rects <- rbind(cbind(cell(3, 1:5), g = "G1"), cbind(cell(1:5, 1), g = "G2"), cbind(cell(5, 5), g = "G3"))
lab_orig <- data.frame(x = c(3, 1, 5), y = c(3, 4.2, 5), lab = c("G[1]", "G[2]", "G[3]"))
p_orig <- base_map("Original") +
  geom_rect(data = orig_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#8a8984", alpha = .25, colour = NA) +
  geom_label(data = lab_orig, aes(x, y, label = lab), parse = TRUE, size = 3, label.size = .2)

# Panel 2: mixture model. Cell tint = average probability of belonging to class 1
# (the class resembling G1/G2); the number is that average in %. The paper plots
# every person instead; per-cell averages show the same pattern more legibly.
mixdat <- aggregate(p ~ A + B, data.frame(A = as.integer(as.character(sim$A)),
                                          B = as.integer(as.character(sim$B)), p = sim$p_class1), mean)
p_mix <- base_map("Mixture model (2 classes)") +
  geom_rect(data = cbind(mixdat, cell(mixdat$A, mixdat$B)),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = p), colour = "white") +
  scale_fill_gradient(low = "#eef3fb", high = "#2a78d6", limits = c(0, 1)) +
  geom_text(data = mixdat, aes(B, A, label = sprintf("%.0f%%", 100 * p)), size = 2.7, colour = ink)

# Panel 3: SEM tree leaves
leaf_cells <- unique(data.frame(A = as.integer(as.character(sim$A)),
                                B = as.integer(as.character(sim$B)), leaf = sim$leaf))
leaf_ids <- sort(unique(leaf_cells$leaf))
leaf_cells$R <- factor(match(leaf_cells$leaf, leaf_ids))
# label each leaf in the cell of that leaf closest to the leaf's centre
leaf_lab <- do.call(rbind, lapply(split(leaf_cells, leaf_cells$R), function(z) {
  ctr <- colMeans(z[c("A", "B")]); z[which.min((z$A - ctr[1])^2 + (z$B - ctr[2])^2), ] }))
p_tree <- base_map(sprintf("SEM tree (%d leaves)", length(leaf_ids))) +
  geom_rect(data = cbind(leaf_cells, cell(leaf_cells$A, leaf_cells$B)),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = R), alpha = .55, colour = "white") +
  scale_fill_manual(values = cat_pal) +
  geom_label(data = leaf_lab, aes(B, A, label = paste0("hat(R)[", R, "]")), parse = TRUE,
             size = 3, label.size = .2)

# Panel 4: top-5 SubgroupSEM subgroups, intensity = LR
parse_sg <- function(desc) {
  a <- if (grepl("A==", desc)) as.integer(sub(".*A=='([0-9])'.*", "\\1", desc)) else 1:5
  b <- if (grepl("B==", desc)) as.integer(sub(".*B=='([0-9])'.*", "\\1", desc)) else 1:5
  data.frame(xmin = min(b) - .5, xmax = max(b) + .5, ymin = min(a) - .5, ymax = max(a) + .5)
}
sg_rects <- do.call(rbind, lapply(seq_len(5), function(i)
  cbind(parse_sg(sg_top5$subgroup[i]), S = i, lr = sg_top5$quality[i])))
sg_rects$lx <- ifelse(sg_rects$xmax - sg_rects$xmin > 1, 4.3 - .2 * (sg_rects$S > 3), (sg_rects$xmin + sg_rects$xmax) / 2)
sg_rects$ly <- ifelse(sg_rects$ymax - sg_rects$ymin > 1, 4.2, (sg_rects$ymin + sg_rects$ymax) / 2)
sg_rects$lx[sg_rects$S == 4] <- 3   # A=3 & B=3 sits inside the A=3 row
p_sg <- base_map("SubgroupSEM (top 5)") +
  geom_rect(data = sg_rects, aes(xmin = xmin + .04, xmax = xmax - .04, ymin = ymin + .04, ymax = ymax - .04,
                                 alpha = lr), fill = "#4a3aa7", colour = "#4a3aa7", linewidth = .3) +
  scale_alpha_continuous(range = c(.12, .45)) +
  geom_label(data = sg_rects, aes(lx, ly, label = paste0("hat(S)[", S, "]")), parse = TRUE,
             size = 3, label.size = .2)

fig2 <- (p_orig | p_mix) / (p_tree | p_sg)
ggsave(file.path(out, "figure2_simulated.pdf"), fig2, width = 6.4, height = 6.6)
ggsave(file.path(out, "figure2_simulated.png"), fig2, width = 6.4, height = 6.6, dpi = 200)

timing$total <- Sys.time() - t_start
writeLines(sprintf("%-22s %8.2f min", names(timing), sapply(timing, as.numeric, units = "mins")),
           file.path(out, "runtime.txt"))
writeLines(c(capture.output(sessionInfo()), "", capture.output(py_config()),
             paste("pysubgroup", py_eval("__import__('importlib.metadata').metadata.version('pysubgroup')"))),
           file.path(out, "sessionInfo.txt"))
cat("\nDone. Everything is in", out, "\n")