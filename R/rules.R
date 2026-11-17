#' Root-to-leaf decision rules of a forest
#'
#' Enumerates every root-to-leaf path of selected trees by depth-first
#' traversal. Each path is one rule: the conjunction of the split conditions
#' met on the way down, with the prediction of the leaf. A case goes to the
#' left child when its value is at most the split value, as in `ranger`.
#'
#' Thresholds are kept at full precision in `lower` and `upper` (the interval
#' each feature must fall in, `lower < value <= upper`), so the rules reproduce
#' the trees exactly. The text in `rule` is rounded for display and should not
#' be used to assign cases: rounding a threshold can move it across a data
#' value.
#'
#' Rules of a forest grown to purity are long and each covers few cases. How
#' often two features occur in the same rule is driven mainly by how often
#' each is split on, so co-occurrence counts are not evidence beyond feature
#' importance.
#'
#' @param fit Output of [fit_forest()], or a `ranger` object.
#' @param trees Indices of the trees to traverse.
#' @param digits Significant digits of the thresholds shown in `rule`.
#' @return A data frame with one row per leaf: `tree`, `leaf` (node ID),
#'   `rule` (display text), `prediction`, `n_conditions`, and the list columns
#'   `features`, `lower` and `upper` (one element per distinct feature on the
#'   path).
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs[1:120, , ])
#' fit <- fit_forest(x, sim$slides$subtype[1:120], num_trees = 20,
#'                   tune_trees = 20, num_threads = 1)
#' r <- extract_rules(fit, trees = 1:2)
#' table(rule_membership(r[r$tree == 1, ], x))  # each case in exactly one leaf
#' @export
extract_rules <- function(fit, trees = 1:20, digits = 6) {
  forest <- if (inherits(fit, "ranger")) fit else fit$forest
  out <- vector("list", length(trees))
  for (t in seq_along(trees)) {
    ti <- ranger::treeInfo(forest, trees[t])
    rownames(ti) <- ti$nodeID
    pred_cols <- names(ti)[startsWith(names(ti), "pred.")]
    leaf_pred <- function(row) {
      if (length(pred_cols)) {
        substring(pred_cols[which.max(unlist(row[pred_cols]))], 6)
      } else as.character(row$prediction)
    }
    rows <- list()
    # explicit stack: node, conditions so far (feature, threshold, side)
    stack <- list(list(node = 0, f = character(0), thr = numeric(0),
                       left = logical(0)))
    while (length(stack)) {
      s <- stack[[length(stack)]]
      stack[[length(stack)]] <- NULL
      row <- ti[as.character(s$node), ]
      if (row$terminal) {
        feats <- unique(s$f)
        lower <- sapply(feats, function(v) max(c(-Inf, s$thr[s$f == v & !s$left])))
        upper <- sapply(feats, function(v) min(c(Inf, s$thr[s$f == v & s$left])))
        txt <- if (!length(s$f)) "<root>" else paste(
          paste0(s$f, ifelse(s$left, " <= ", " > "), signif(s$thr, digits)),
          collapse = " & ")
        rows[[length(rows) + 1]] <- data.frame(
          tree = trees[t], leaf = s$node, rule = txt,
          prediction = leaf_pred(row), n_conditions = length(s$f),
          features = I(list(feats)), lower = I(list(unname(lower))),
          upper = I(list(unname(upper))))
        next
      }
      for (side in c(FALSE, TRUE)) {   # right pushed first, so left is visited first
        stack[[length(stack) + 1]] <- list(
          node = if (side) row$leftChild else row$rightChild,
          f = c(s$f, row$splitvarName), thr = c(s$thr, row$splitval),
          left = c(s$left, side))
      }
    }
    out[[t]] <- do.call(rbind, rows)
  }
  do.call(rbind, out)
}

#' Rule satisfied by each case
#'
#' @param rules Rows of [extract_rules()] for a single tree.
#' @param x Feature matrix with named columns.
#' @return For each row of `x`, the number of rules it satisfies. For the
#'   rules of one tree this is one for every case.
#' @export
rule_membership <- function(rules, x) {
  rowSums(rule_matrix(rules, x))
}

#' Case-by-rule membership matrix
#'
#' @inheritParams rule_membership
#' @return A logical matrix, cases by rules.
#' @export
rule_matrix <- function(rules, x) {
  sapply(seq_len(nrow(rules)), function(r) {
    ok <- rep(TRUE, nrow(x))
    f <- rules$features[[r]]
    for (k in seq_along(f)) {
      v <- x[, f[k]]
      ok <- ok & v > rules$lower[[r]][k] & v <= rules$upper[[r]][k]
    }
    ok
  })
}
