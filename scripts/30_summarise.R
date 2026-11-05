# Summarises the sensitivity study and evaluates H1-H4.
# Run from the project root after scripts/20_run_grid.R:
#   Rscript scripts/30_summarise.R
# Output: results/sensitivity/summary/*.csv, hypotheses.txt, figures (PDF + PNG)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
in_dir <- file.path("results", "sensitivity")
out_dir <- file.path(in_dir, "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
settings <- jsonlite::read_json(file.path(in_dir, "settings.json"), simplifyVector = TRUE)
design <- read.csv(file.path(in_dir, "design.csv"))
robust <- 0.8   # report Sect. 3.3: recovery is robust if the full recovery rate >= 0.8

## Completeness: never summarise a partial or failed run silently ------------------
errors <- list.files(file.path(in_dir, "errors"), full.names = TRUE)
if (length(errors)) {
  for (e in errors) { x <- readRDS(e); message("ERROR in ", x$task$id, ": ", x$error) }
  stop(length(errors), " dataset(s) failed; fix or rerun them before summarising.")
}
files <- file.path(in_dir, "runs", paste0(design$id, ".rds"))
if (!all(file.exists(files)))
  stop(sum(!file.exists(files)), " of ", length(files), " datasets are missing; the grid has not finished.")

## One row per dataset and method ----------------------------------------------------
not_failed <- c("ok", "all_in_subgroup")
rows <- lapply(files, function(f) {
  x <- readRDS(f)
  search_rows <- lapply(names(x$searches), function(im) {
    s <- x$searches[[im]]
    failed <- s$log[!s$log$status %in% not_failed, ]
    data.frame(method = paste0("SubgroupSEM (", im, ")"),
               G1 = s$recovered[["G1"]], G2 = s$recovered[["G2"]], G3 = s$recovered[["G3"]],
               top1_value = s$top$quality[1], top1 = s$top$description[1],
               n_failed = nrow(failed),
               planted_failed = paste(intersect(planted_descriptions(), failed$description), collapse = "; "),
               n_leaves = NA_integer_, seconds = s$seconds)
  })
  # Recomputed from the saved leaf of every cell (cell-based Jaccard, see R/score.R).
  tree_rec <- recovery_tree(x$tree$leaf_by_cell, settings$tree_jaccard)$recovered
  tree_row <- data.frame(method = "SEM tree",
                         G1 = tree_rec[["G1"]], G2 = tree_rec[["G2"]], G3 = tree_rec[["G3"]],
                         top1_value = NA_real_, top1 = NA_character_, n_failed = NA_integer_,
                         planted_failed = NA_character_, n_leaves = x$tree$n_leaves, seconds = x$tree$seconds)
  task <- x$task[c("id", "N", "s", "r", "seed")]
  rownames(task) <- NULL
  cbind(task, rbind(do.call(rbind, search_rows), tree_row))
})
runs <- do.call(rbind, rows)
runs$n_recovered <- runs$G1 + runs$G2 + runs$G3
runs$full <- runs$n_recovered == 3
write.csv(runs, file.path(out_dir, "runs.csv"), row.names = FALSE)

## Aggregate per condition and method -----------------------------------------------
agg <- do.call(rbind, lapply(split(runs, list(runs$method, runs$N, runs$s), drop = TRUE), function(d)
  data.frame(method = d$method[1], N = d$N[1], s = d$s[1], R = nrow(d),
             full_k = sum(d$full), full_rate = mean(d$full),
             mean_recovered = mean(d$n_recovered), sd_recovered = sd(d$n_recovered),
             rate_G1 = mean(d$G1), rate_G2 = mean(d$G2), rate_G3 = mean(d$G3),
             top1_mean = mean(d$top1_value), top1_sd = sd(d$top1_value),
             failed_fits_mean = mean(d$n_failed),
             datasets_planted_failed = sum(!is.na(d$planted_failed) & d$planted_failed != ""),
             leaves_mean = mean(d$n_leaves))))
agg <- agg[order(agg$method, agg$N, agg$s), ]
stopifnot(all(agg$R == settings$R))
write.csv(agg, file.path(out_dir, "by_condition.csv"), row.names = FALSE)

## Hypotheses --------------------------------------------------------------------------
sink(file.path(out_dir, "hypotheses.txt"), split = TRUE)
lr <- agg[agg$method == "SubgroupSEM (LR)", ]
show <- function(d) print(d[c("method", "N", "s", "full_k", "R", "full_rate", "mean_recovered",
                              "rate_G1", "rate_G2", "rate_G3", "failed_fits_mean")], row.names = FALSE)

cat("H1 (effect size): LR at N = 5000, robust (full rate >=", robust, ") at s = 0.5 and 0.25?\n")
h1 <- lr[lr$N == 5000 & lr$s > 0, ]; show(h1)
cat("-> H1", if (all(h1$full_rate[h1$s %in% c(0.25, 0.5)] >= robust)) "SUPPORTED" else "NOT SUPPORTED", "\n\n")

cat("H2 (sample size): LR at s = 1, robust at N = 2500 and 500?\n")
h2 <- lr[lr$s == 1, ]; show(h2)
cat("-> H2", if (all(h2$full_rate[h2$N %in% c(500, 2500)] >= robust)) "SUPPORTED" else "NOT SUPPORTED", "\n\n")

cat("H3 (IM): per condition (s > 0) and planted group, is the group robustly recovered (rate >=",
    robust, ") by each IM? H3 holds if the three IMs agree everywhere.\n")
ims <- agg[grepl("SubgroupSEM", agg$method) & agg$s > 0, ]
long <- do.call(rbind, lapply(c("G1", "G2", "G3"), function(g)
  data.frame(N = ims$N, s = ims$s, group = g, method = ims$method,
             robust = ims[[paste0("rate_", g)]] >= robust)))
wide <- reshape(long, idvar = c("N", "s", "group"), timevar = "method", direction = "wide")
names(wide) <- sub("robust.SubgroupSEM ", "", names(wide))
wide$agree <- apply(wide[grep("^\\(", names(wide))], 1, function(x) length(unique(x)) == 1)
print(wide[order(wide$group, wide$N, wide$s), ], row.names = FALSE)
cat("-> H3", if (all(wide$agree)) "SUPPORTED" else
  paste0("NOT SUPPORTED (", sum(!wide$agree), " of ", nrow(wide), " condition-group pairs disagree)"), "\n\n")

cat("H4 (superiority): mean number of recovered groups, LR vs SEM tree, s > 0\n")
tree <- agg[agg$method == "SEM tree", ]
h4 <- merge(lr[lr$s > 0, c("N", "s", "mean_recovered")], tree[tree$s > 0, c("N", "s", "mean_recovered")],
            by = c("N", "s"), suffixes = c("_LR", "_tree"))
h4$LR_at_least_tree <- h4$mean_recovered_LR >= h4$mean_recovered_tree
print(h4, row.names = FALSE)
cat("-> H4", if (all(h4$LR_at_least_tree)) "SUPPORTED" else "NOT SUPPORTED", "\n\n")

cat("Control s = 0 (nothing planted): IM value of the top-ranked subgroup, mean (sd)\n")
print(agg[agg$s == 0 & grepl("SubgroupSEM", agg$method), c("method", "N", "top1_mean", "top1_sd", "full_k")], row.names = FALSE)
cat("\nFailed fits per search (mean) and datasets in which a planted group itself failed to fit\n")
print(agg[grepl("SubgroupSEM", agg$method) & agg$failed_fits_mean > 0,
          c("method", "N", "s", "failed_fits_mean", "datasets_planted_failed")], row.names = FALSE)
sink()

## Figures ------------------------------------------------------------------------------
agg$method <- factor(agg$method, levels = c("SubgroupSEM (LR)", "SubgroupSEM (Wald)", "SubgroupSEM (UD)", "SEM tree"))
fig_full <- ggplot(agg[agg$s > 0, ], aes(factor(N), factor(s), fill = full_rate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%d/%d", full_k, R)), size = 3) +
  scale_fill_gradient("Full recovery rate", low = "#eef3fb", high = "#2a78d6", limits = c(0, 1)) +
  facet_wrap(~ method, nrow = 1) + labs(x = "N", y = "Effect size s") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(file.path(out_dir, "fig_full_recovery.pdf"), fig_full, width = 8, height = 3.2)
ggsave(file.path(out_dir, "fig_full_recovery.png"), fig_full, width = 8, height = 3.2, dpi = 200)

groups <- do.call(rbind, lapply(c("G1", "G2", "G3"), function(g)
  data.frame(agg[agg$s > 0, c("method", "N", "s")], group = g, rate = agg[agg$s > 0, paste0("rate_", g)])))
fig_groups <- ggplot(groups, aes(factor(N), factor(s), fill = rate)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rate)), size = 2.6) +
  scale_fill_gradient("Recovery rate", low = "#eef3fb", high = "#2a78d6", limits = c(0, 1)) +
  facet_grid(group ~ method) + labs(x = "N", y = "Effect size s") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(file.path(out_dir, "fig_group_recovery.pdf"), fig_groups, width = 8, height = 6)
ggsave(file.path(out_dir, "fig_group_recovery.png"), fig_groups, width = 8, height = 6, dpi = 200)
