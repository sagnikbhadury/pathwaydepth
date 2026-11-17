# Plotting helpers (base graphics). One fixed colour per feature set or class.

.pd_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6", "#e34948")
.pd_ink <- "#0b0b0b"; .pd_muted <- "#52514e"; .pd_grid <- "#e6e5e1"

.pd_frame <- function(xaxis = TRUE, yaxis = TRUE) {
  graphics::box(col = .pd_grid)
  if (xaxis) graphics::axis(1, col = .pd_grid, col.axis = .pd_muted, cex.axis = 0.85)
  if (yaxis) graphics::axis(2, col = .pd_grid, col.axis = .pd_muted, cex.axis = 0.85, las = 1)
}

.pd_colours <- function(names, col = NULL) {
  if (is.null(col)) col <- .pd_cols[(seq_along(names) - 1) %% length(.pd_cols) + 1]
  stats::setNames(col, names)
}

#' Colours used by the plotting functions
#'
#' @param n Number of colours (at most six).
#' @return A character vector of hex colours, in a fixed order.
#' @examples
#' depth_colours(3)
#' @export
depth_colours <- function(n = 6) .pd_cols[seq_len(n)]

#' ROC curves of several classifiers on one test set
#'
#' @param probs A named list (or a matrix with named columns) of predicted
#'   probabilities of the positive class, one element per classifier.
#' @param y Two-level factor of true labels, the second level positive.
#' @param col Optional colours, one per classifier.
#' @param main Title.
#' @return Invisibly, the AUC of each classifier.
#' @examples
#' set.seed(1)
#' y <- factor(rep(c("a", "b"), each = 60))
#' plot_roc(list(strong = c(rnorm(60), rnorm(60, 1.5)),
#'               weak = c(rnorm(60), rnorm(60, 0.5))), y)
#' @export
plot_roc <- function(probs, y, col = NULL, main = "") {
  if (is.matrix(probs)) probs <- as.list(as.data.frame(probs))
  y <- factor(y)
  col <- .pd_colours(names(probs), col)
  op <- graphics::par(pty = "s"); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
                 xlab = "1 - specificity", ylab = "Sensitivity", main = main,
                 font.main = 1, col.lab = .pd_ink)
  .pd_frame(); graphics::abline(0, 1, col = .pd_grid, lty = 2)
  auc <- numeric(0)
  for (k in names(probs)) {
    r <- pROC::roc(y, probs[[k]], levels = levels(y), direction = "<", quiet = TRUE)
    graphics::lines(1 - r$specificities, r$sensitivities, col = col[k], lwd = 2)
    auc[k] <- as.numeric(pROC::auc(r))
  }
  graphics::legend("bottomright", bty = "n", lwd = 2, col = col, cex = 0.8,
                   text.col = .pd_ink, title = "Features (AUC)", title.col = .pd_muted,
                   legend = sprintf("%s (%.2f)", names(probs), auc))
  invisible(auc)
}

#' Estimates with intervals for several classifiers
#'
#' A dot-and-interval display, one row per classifier. Use it for the output
#' of [evaluate()] (one metric at a time) or for the mean and percentiles of
#' the AUC over repeated splits.
#'
#' @param estimate,lower,upper Numeric vectors, named by classifier.
#' @param xlab Axis label.
#' @param col Optional colours, one per classifier.
#' @param xlim Axis limits.
#' @param main Title.
#' @return `NULL`, invisibly.
#' @examples
#' plot_intervals(c(A = 0.80, B = 0.74), c(A = 0.75, B = 0.68),
#'                c(A = 0.85, B = 0.80), xlab = "Test AUC")
#' @export
plot_intervals <- function(estimate, lower, upper, xlab = "", col = NULL,
                           xlim = range(c(lower, upper)) + c(-0.02, 0.02), main = "") {
  k <- names(estimate); n <- length(k)
  col <- .pd_colours(k, col)
  op <- graphics::par(mar = c(4, max(nchar(k)) * 0.48 + 1.5, 2.2, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5), axes = FALSE, xlab = xlab,
                 ylab = "", main = main, font.main = 1, col.lab = .pd_ink)
  .pd_frame(yaxis = FALSE)
  graphics::abline(v = graphics::axTicks(1), col = .pd_grid)
  yy <- rev(seq_len(n))
  graphics::segments(lower, yy, upper, yy, col = col, lwd = 2)
  graphics::points(estimate, yy, pch = 16, col = col, cex = 1.4)
  graphics::mtext(k, side = 2, at = yy, las = 1, line = 0.5, cex = 0.75, col = .pd_ink)
  graphics::text(estimate, yy, sprintf("%.2f", estimate), pos = 3, cex = 0.7, col = .pd_muted)
  invisible(NULL)
}

#' Volcano plot of a feature-wise group comparison
#'
#' @param cmp Output of [compare_groups()].
#' @param fdr False discovery rate threshold.
#' @param strong Absolute effect size above which a significant feature is
#'   drawn as a strong effect.
#' @param label Number of features to label, by decreasing absolute effect.
#' @param main Title.
#' @return `NULL`, invisibly.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs)
#' plot_volcano(compare_groups(x, sim$slides$subtype))
#' @export
plot_volcano <- function(cmp, fdr = 0.05, strong = 0.3, label = 5, main = "") {
  cmp <- cmp[!is.na(cmp$d), ]
  sig <- cmp$fdr < fdr
  cls <- ifelse(sig & abs(cmp$d) > strong, 2, ifelse(sig, 1, 0))
  col <- c("#b9b8b2", .pd_cols[4], .pd_cols[6])[cls + 1]
  yv <- -log10(cmp$fdr)
  graphics::plot(cmp$d, yv, pch = 16, col = col, axes = FALSE, main = main, font.main = 1,
                 ylim = c(0, max(yv) * 1.22),
                 xlab = sprintf("Cohen's d (%s - %s)", sub("mean_", "", names(cmp)[3]),
                                sub("mean_", "", names(cmp)[2])),
                 ylab = expression(-log[10] ~ FDR), col.lab = .pd_ink)
  .pd_frame()
  graphics::abline(h = -log10(fdr), v = 0, col = .pd_grid, lty = c(2, 1))
  top <- utils::head(order(-abs(cmp$d)), label)
  graphics::text(cmp$d[top], yv[top], cmp$feature[top], cex = 0.6, col = .pd_muted,
                 pos = ifelse(cmp$d[top] > 0, 2, 4), xpd = NA)
  graphics::legend("top", horiz = TRUE, bty = "n", pch = 16, cex = 0.7, text.col = .pd_ink,
                   x.intersp = 0.6, text.width = NA,
                   col = c(.pd_cols[6], .pd_cols[4], "#b9b8b2"),
                   legend = c(sprintf("FDR < %g, |d| > %g", fdr, strong),
                              sprintf("FDR < %g", fdr), "not significant"))
  invisible(NULL)
}

#' Effect size of a depth feature against that of the edge weight
#'
#' Depth is unsigned, so a difference in depth between two groups need not
#' have the sign of the difference in edge weight. This plot shows the two
#' effect sizes against each other, one point per edge.
#'
#' @param depth_cmp Output of [compare_groups()] on depth features.
#' @param weight_cmp Output of [compare_groups()] on the edge weights.
#' @param fdr Threshold for marking a depth feature as significant.
#' @param ylab Axis label for the depth effect size.
#' @return Invisibly, the Pearson correlation of the two effect sizes.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs); y <- sim$slides$subtype
#' plot_depth_vs_weight(compare_groups(within_depth(x, B = 0), y),
#'                      compare_groups(x, y))
#' @export
plot_depth_vs_weight <- function(depth_cmp, weight_cmp, fdr = 0.05,
                                 ylab = "Depth: Cohen's d") {
  w <- weight_cmp$d[match(depth_cmp$feature, weight_cmp$feature)]
  ok <- !is.na(w) & !is.na(depth_cmp$d)
  d <- depth_cmp$d[ok]; w <- w[ok]; sig <- depth_cmp$fdr[ok] < fdr
  lim <- max(abs(c(d, w))) * c(-1.1, 1.1)
  graphics::plot(w, d, xlim = lim, ylim = lim, pch = 16, axes = FALSE,
                 col = ifelse(sig, .pd_ink, "#b9b8b2"), ylab = ylab,
                 xlab = "Edge weight: Cohen's d", col.lab = .pd_ink)
  .pd_frame(); graphics::abline(h = 0, v = 0, col = .pd_grid)
  r <- stats::cor(w, d)
  graphics::mtext(sprintf("r = %.2f", r), side = 3, adj = 1, line = 0.2, cex = 0.75,
                  col = .pd_muted)
  graphics::legend("topright", bty = "n", pch = 16, cex = 0.75, text.col = .pd_ink,
                   col = c(.pd_ink, "#b9b8b2"),
                   legend = c(sprintf("depth FDR < %g", fdr), "not significant"))
  invisible(r)
}

#' Share of graphs in which each node is absent, by class
#'
#' @param graphs Array of directed proximity matrices, graphs by nodes by
#'   nodes.
#' @param y Two-level factor, one label per graph.
#' @return Invisibly, a matrix of absence shares (nodes by classes).
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' plot_presence(sim$graphs, sim$slides$subtype)
#' @export
plot_presence <- function(graphs, y) {
  y <- factor(y)
  absent <- !presence_matrix(graphs)
  s <- sapply(levels(y), function(l) colMeans(absent[y == l, , drop = FALSE]))
  s <- s[order(rowSums(s)), , drop = FALSE]
  col <- stats::setNames(.pd_cols[c(5, 6)], levels(y))
  op <- graphics::par(mar = c(4, 6, 1.5, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0.5, nrow(s) + 0.5), axes = FALSE, ylab = "",
                 xlab = "Share of graphs in which the node is absent", col.lab = .pd_ink)
  .pd_frame(yaxis = FALSE); graphics::abline(v = seq(0, 1, 0.25), col = .pd_grid)
  for (i in seq_len(nrow(s))) {
    graphics::segments(s[i, 1], i, s[i, 2], i, col = "#b9b8b2", lwd = 2)
    graphics::points(s[i, ], c(i, i), pch = 16, col = col, cex = 1.3)
  }
  graphics::mtext(rownames(s), side = 2, at = seq_len(nrow(s)), las = 1, line = 0.5,
                  cex = 0.75, col = .pd_ink)
  graphics::legend("bottomright", bty = "n", pch = 16, col = col, legend = levels(y),
                   cex = 0.8, text.col = .pd_ink)
  invisible(s)
}

#' Permutation importance of the leading features
#'
#' @param fit Output of [fit_forest()], or a `ranger` object fitted with
#'   `importance = "permutation"`.
#' @param top Number of features to show.
#' @param col Bar colour.
#' @param main Title.
#' @return Invisibly, the importance of the features shown.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' fit <- fit_forest(edge_matrix(sim$graphs), sim$slides$subtype,
#'                   num_trees = 100, tune_trees = 50, num_threads = 1)
#' plot_importance(fit, top = 10)
#' @export
plot_importance <- function(fit, top = 20, col = depth_colours(1), main = "") {
  forest <- if (inherits(fit, "ranger")) fit else fit$forest
  imp <- utils::head(sort(forest$variable.importance, decreasing = TRUE), top)
  op <- graphics::par(mar = c(4, max(nchar(names(imp))) * 0.5 + 1.5, 2.2, 1))
  on.exit(graphics::par(op))
  n <- length(imp); yy <- rev(seq_len(n))
  graphics::plot(NA, xlim = c(0, max(imp) * 1.05), ylim = c(0.5, n + 0.5), axes = FALSE,
                 xlab = "Permutation importance", ylab = "", main = main, font.main = 1,
                 col.lab = .pd_ink, xaxs = "i")
  .pd_frame(yaxis = FALSE)
  graphics::rect(0, yy - 0.3, imp, yy + 0.3, col = col, border = NA)
  graphics::mtext(names(imp), side = 2, at = yy, las = 1, line = 0.5, cex = 0.7, col = .pd_ink)
  invisible(imp)
}

#' Distribution of the number of conditions per rule
#'
#' @param rules A named list of outputs of [extract_rules()], one per model.
#' @param col Optional colours, one per model.
#' @return Invisibly, the median number of conditions per model.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' fit <- fit_forest(edge_matrix(sim$graphs), sim$slides$subtype,
#'                   num_trees = 30, tune_trees = 30, num_threads = 1)
#' plot_rule_lengths(list(`Edge weights` = extract_rules(fit, 1:10)))
#' @export
plot_rule_lengths <- function(rules, col = NULL) {
  col <- .pd_colours(names(rules), col)
  mx <- max(sapply(rules, function(r) max(r$n_conditions)))
  op <- graphics::par(mfrow = c(1, length(rules)), mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  med <- sapply(names(rules), function(k) {
    n <- rules[[k]]$n_conditions
    tab <- table(factor(n, levels = seq_len(mx)))
    graphics::plot(NA, xlim = c(0.5, mx + 0.5), ylim = c(0, max(tab) * 1.08), axes = FALSE,
                   xlab = "Conditions per rule", ylab = "Number of rules", yaxs = "i",
                   main = sprintf("%s: %s rules", k, format(length(n), big.mark = ",")),
                   font.main = 1, cex.main = 0.95, col.lab = .pd_ink)
    .pd_frame()
    graphics::rect(seq_len(mx) - 0.4, 0, seq_len(mx) + 0.4, tab, col = col[k], border = NA)
    m <- stats::median(n)
    graphics::abline(v = m, lty = 2, col = .pd_ink)
    graphics::text(m, max(tab), sprintf("median %g", m), pos = 4, cex = 0.75, col = .pd_muted)
    m
  })
  invisible(med)
}

#' Co-occurrence of features in decision rules
#'
#' Counts, for every pair of features, the number of rules in which both
#' occur. A feature that is split near the root of a tree is on every rule of
#' that tree, so these counts mainly reflect how often each feature is used;
#' `nodes$splits` is returned so that this can be checked.
#'
#' @param rules Output of [extract_rules()].
#' @return A list with `pairs` (data frame `feature1`, `feature2`, `count`,
#'   by decreasing count) and `nodes` (data frame `feature`, `degree`,
#'   `total`, `rules`: the number of other features it occurs with, the sum of
#'   its pair counts, and the number of rules that use it).
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' fit <- fit_forest(edge_matrix(sim$graphs), sim$slides$subtype,
#'                   num_trees = 30, tune_trees = 30, num_threads = 1)
#' co <- cooccurrence(extract_rules(fit, 1:10))
#' head(co$pairs, 3)
#' @export
cooccurrence <- function(rules) {
  feats <- sort(unique(unlist(rules$features)))
  m <- matrix(0L, length(feats), length(feats), dimnames = list(feats, feats))
  used <- stats::setNames(integer(length(feats)), feats)
  for (f in rules$features) {
    used[f] <- used[f] + 1L
    if (length(f) > 1) m[f, f] <- m[f, f] + 1L
  }
  diag(m) <- 0L
  idx <- which(upper.tri(m) & m > 0, arr.ind = TRUE)
  pairs <- data.frame(feature1 = feats[idx[, 1]], feature2 = feats[idx[, 2]],
                      count = m[idx], stringsAsFactors = FALSE)
  pairs <- pairs[order(-pairs$count), ]; rownames(pairs) <- NULL
  nodes <- data.frame(feature = feats, degree = rowSums(m > 0), total = rowSums(m),
                      rules = as.integer(used), row.names = NULL)
  list(pairs = pairs, nodes = nodes[order(-nodes$total), ])
}

#' Circular network of the most frequent feature pairs
#'
#' @param co Output of [cooccurrence()].
#' @param top Number of pairs to draw.
#' @param col Node colour.
#' @param main Title.
#' @return `NULL`, invisibly.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' fit <- fit_forest(edge_matrix(sim$graphs), sim$slides$subtype,
#'                   num_trees = 30, tune_trees = 30, num_threads = 1)
#' plot_cooccurrence(cooccurrence(extract_rules(fit, 1:10)), top = 25)
#' @export
plot_cooccurrence <- function(co, top = 40, col = depth_colours(1), main = "") {
  p <- utils::head(co$pairs, top)
  nodes <- unique(c(p$feature1, p$feature2))
  deg <- table(factor(c(p$feature1, p$feature2), levels = nodes))
  ang <- seq(0, 2 * pi, length.out = length(nodes) + 1)[-1]
  xy <- cbind(cos(ang), sin(ang)); rownames(xy) <- nodes
  op <- graphics::par(mar = c(1, 1, 2.2, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7), asp = 1, axes = FALSE,
                 xlab = "", ylab = "", main = main, font.main = 1)
  w <- 0.5 + 4 * p$count / max(p$count)
  graphics::segments(xy[p$feature1, 1], xy[p$feature1, 2], xy[p$feature2, 1],
                     xy[p$feature2, 2], lwd = w, col = grDevices::adjustcolor(.pd_muted, 0.35))
  graphics::points(xy, pch = 21, bg = col, col = "white", lwd = 1.5,
                   cex = 1 + 2.2 * as.numeric(deg) / max(deg))
  for (i in seq_along(nodes)) {            # labels run outwards along the radius
    left <- xy[i, 1] < 0
    graphics::text(1.09 * xy[i, 1], 1.09 * xy[i, 2], nodes[i], cex = 0.62, col = .pd_ink,
                   xpd = NA, srt = ang[i] * 180 / pi + if (left) 180 else 0,
                   adj = c(if (left) 1 else 0, 0.5))
  }
  invisible(NULL)
}
