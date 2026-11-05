# Scoring of search and tree results against the planted groups.

# pysubgroup prints descriptions like "A=='5' AND B=='5'" or "grade==7".
# Canonical form: "A=5 & B=5", selectors sorted by attribute.
canonical_description <- function(desc) {
  stopifnot(length(desc) == 1, is.character(desc))
  if (desc %in% c("", "Dataset")) return("")
  sel <- strsplit(desc, " AND ", fixed = TRUE)[[1]]
  m <- regmatches(sel, regexec("^\\s*([^=]+?)==\\s*'?([^']*)'?\\s*$", sel))
  ok <- vapply(m, length, integer(1)) == 3
  if (!all(ok)) stop("Cannot parse subgroup description: ", desc)
  attr <- vapply(m, `[`, character(1), 2)
  val <- vapply(m, `[`, character(1), 3)
  o <- order(attr)
  paste(paste0(attr[o], "=", val[o]), collapse = " & ")
}

# Which planted groups are among the top_n subgroups of a search result.
recovery_sgsem <- function(top_descriptions, top_n = 3,
                           planted = planted_descriptions()) {
  found <- head(top_descriptions, top_n)
  vapply(planted, function(p) p %in% found, logical(1))
}

jaccard <- function(a, b) {
  stopifnot(is.logical(a), is.logical(b), length(a) == length(b), !anyNA(a), !anyNA(b))
  union <- sum(a | b)
  if (union == 0) return(NA_real_)
  sum(a & b) / union
}

# A planted group counts as recovered by the tree if its best-matching leaf has a
# Jaccard index of at least `threshold` with it, counted in cells of the A x B grid.
# The tree only splits on A and B, so every cell lies in one leaf. Counting cells
# (not pupils) makes "G2 without the cell it shares with G1" exactly 0.8, the best
# a tree can do, whatever the random cell sizes.
# leaf_by_cell: one row per occupied cell, columns A, B, leaf.
recovery_tree <- function(leaf_by_cell, threshold = 0.8) {
  stopifnot(all(c("A", "B", "leaf") %in% names(leaf_by_cell)), !anyNA(leaf_by_cell$leaf),
            !anyDuplicated(leaf_by_cell[c("A", "B")]))
  leaf <- leaf_by_cell$leaf
  members <- planted_membership(leaf_by_cell$A, leaf_by_cell$B)
  best <- vapply(members, function(g)
    max(vapply(unique(leaf), function(l) jaccard(leaf == l, g), numeric(1))),
    numeric(1))
  list(recovered = best >= threshold, best_jaccard = best)
}
