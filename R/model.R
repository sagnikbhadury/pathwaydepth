#' Patient identifier from a TCGA slide barcode
#'
#' @param slide Character vector of slide barcodes such as
#'   `"TCGA-05-4244-01Z-00-DX1"`.
#' @return The first twelve characters, the TCGA participant barcode.
#' @examples
#' patient_id("TCGA-05-4244-01Z-00-DX1")
#' @export
patient_id <- function(slide) substr(slide, 1, 12)

#' Stratified split that keeps every patient on one side
#'
#' Patients, not slides, are assigned to the training or the test set, within
#' each class. All slides of a patient follow the patient.
#'
#' @param group Vector of patient identifiers, one per slide.
#' @param y Factor of class labels, one per slide. A patient must have a single
#'   class.
#' @param p Proportion of patients assigned to the training set.
#' @return A logical vector, `TRUE` for training slides.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' tr <- grouped_split(sim$slides$patient, sim$slides$subtype)
#' any(sim$slides$patient[tr] %in% sim$slides$patient[!tr])
#' @export
grouped_split <- function(group, y, p = 0.7) {
  first <- !duplicated(group)
  g <- group[first]
  cls <- y[first]
  train_groups <- character(0)
  for (k in unique(cls)) {
    gk <- g[cls == k]
    train_groups <- c(train_groups, sample(gk, round(p * length(gk))))
  }
  group %in% train_groups
}

#' Random forest classifier tuned on out-of-bag predictions
#'
#' Fits probability forests with `ranger` over a small grid of `mtry` and
#' `min.node.size`, selects the pair with the largest out-of-bag area under the
#' ROC curve, and refits with `num_trees` trees. Only the training data are
#' used. The classification threshold maximises Youden's index on the
#' out-of-bag probabilities of the final forest, which, unlike in-bag
#' probabilities, are not separated perfectly by a forest grown to purity.
#'
#' @param x Numeric matrix or data frame of training features.
#' @param y Two-level factor. The second level is the positive class.
#' @param num_trees Number of trees in the final forest.
#' @param tune_trees Number of trees in each tuning forest.
#' @param mtry Candidate values of `mtry`. Defaults to 2, the square root and
#'   one third of the number of features.
#' @param min_node_size Candidate minimal node sizes.
#' @param seed Seed passed to `ranger`.
#' @param num_threads Threads passed to `ranger`.
#' @return A list with the final forest (`forest`), the selected parameters
#'   (`tune`), the tuning table (`grid`), the out-of-bag probabilities of the
#'   positive class (`oob_prob`), their AUC (`oob_auc`) and the threshold
#'   (`threshold`).
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs[1:120, , ])
#' fit <- fit_forest(x, sim$slides$subtype[1:120], num_trees = 100,
#'                   tune_trees = 50, num_threads = 1)
#' fit$oob_auc
#' @export
fit_forest <- function(x, y, num_trees = 1500, tune_trees = 500, mtry = NULL,
                       min_node_size = c(1, 5), seed = 1, num_threads = NULL) {
  y <- factor(y)
  stopifnot(nlevels(y) == 2)
  p <- ncol(x)
  if (is.null(mtry)) {
    mtry <- unique(pmax(1, pmin(p, c(2, floor(sqrt(p)), floor(p / 3)))))
  }
  w <- c(1, sum(y == levels(y)[1]) / sum(y == levels(y)[2]))
  names(w) <- levels(y)
  one <- function(m, s, trees) {
    ranger::ranger(x = x, y = y, num.trees = trees, mtry = m,
                   min.node.size = s, probability = TRUE,
                   importance = "permutation", class.weights = w,
                   seed = seed, num.threads = num_threads)
  }
  auc_of <- function(f) {
    as.numeric(pROC::auc(pROC::roc(y, f$predictions[, 2], levels = levels(y),
                                   direction = "<", quiet = TRUE)))
  }
  grid <- expand.grid(mtry = mtry, min_node_size = min_node_size)
  grid$oob_auc <- NA_real_
  for (r in seq_len(nrow(grid))) {
    grid$oob_auc[r] <- auc_of(one(grid$mtry[r], grid$min_node_size[r],
                                  tune_trees))
  }
  best <- grid[which.max(grid$oob_auc), ]
  forest <- one(best$mtry, best$min_node_size, num_trees)
  oob <- forest$predictions[, 2]
  roc <- pROC::roc(y, oob, levels = levels(y), direction = "<", quiet = TRUE)
  thr <- pROC::coords(roc, "best", ret = "threshold",
                      best.method = "youden")$threshold[1]
  list(forest = forest, tune = best, grid = grid, oob_prob = oob,
       oob_auc = as.numeric(pROC::auc(roc)), threshold = thr,
       levels = levels(y))
}

#' Predicted probability of the positive class
#'
#' @param fit Output of [fit_forest()].
#' @param newdata Feature matrix with the training columns.
#' @return Numeric vector of probabilities.
#' @export
predict_forest <- function(fit, newdata) {
  stats::predict(fit$forest, data = newdata)$predictions[, 2]
}

#' Test-set performance with confidence intervals
#'
#' @param prob Predicted probabilities of the positive class.
#' @param y Two-level factor of true labels, the second level positive.
#' @param threshold Classification threshold, fixed before seeing `y`.
#' @param conf Confidence level.
#' @return A data frame with one row per metric (AUC, accuracy, sensitivity,
#'   specificity) and columns `estimate`, `lower`, `upper`. The AUC interval is
#'   DeLong's; the others are exact binomial intervals.
#' @examples
#' set.seed(1)
#' y <- factor(rep(c("a", "b"), each = 50))
#' evaluate(c(rnorm(50), rnorm(50, 1)), y, threshold = 0.5)
#' @export
evaluate <- function(prob, y, threshold, conf = 0.95) {
  y <- factor(y)
  roc <- pROC::roc(y, prob, levels = levels(y), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc, conf.level = conf, method = "delong"))
  pred <- prob >= threshold
  pos <- y == levels(y)[2]
  bin <- function(k, n) {
    b <- stats::binom.test(k, n, conf.level = conf)
    c(k / n, b$conf.int[1], b$conf.int[2])
  }
  m <- rbind(AUC = ci[c(2, 1, 3)],
             Accuracy = bin(sum(pred == pos), length(pos)),
             Sensitivity = bin(sum(pred & pos), sum(pos)),
             Specificity = bin(sum(!pred & !pos), sum(!pos)))
  data.frame(metric = rownames(m), estimate = m[, 1], lower = m[, 2],
             upper = m[, 3], row.names = NULL)
}

#' Paired comparison of two classifiers on one test set
#'
#' DeLong's test for two correlated ROC curves, with the difference in AUC and
#' its confidence interval.
#'
#' @param prob1,prob2 Predicted probabilities from the two classifiers for the
#'   same test cases.
#' @param y Two-level factor of true labels.
#' @param conf Confidence level.
#' @return A data frame with the two AUCs, their difference (`prob1` minus
#'   `prob2`), its confidence limits and the p-value.
#' @examples
#' set.seed(1)
#' y <- factor(rep(c("a", "b"), each = 50))
#' s <- c(rnorm(50), rnorm(50, 1))
#' compare_auc(s, s + rnorm(100), y)
#' @export
compare_auc <- function(prob1, prob2, y, conf = 0.95) {
  y <- factor(y)
  r1 <- pROC::roc(y, prob1, levels = levels(y), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(y, prob2, levels = levels(y), direction = "<", quiet = TRUE)
  t <- pROC::roc.test(r1, r2, paired = TRUE, method = "delong",
                      conf.level = conf)
  data.frame(auc1 = as.numeric(pROC::auc(r1)), auc2 = as.numeric(pROC::auc(r2)),
             difference = as.numeric(pROC::auc(r1) - pROC::auc(r2)),
             lower = t$conf.int[1], upper = t$conf.int[2], p_value = t$p.value)
}

#' Feature-wise comparison of two groups
#'
#' Welch's two-sample t-test for every column, Cohen's d with the pooled
#' standard deviation (second level minus first level), and Benjamini-Hochberg
#' adjustment over the columns.
#'
#' @param x Numeric matrix, cases by features.
#' @param y Two-level factor.
#' @return A data frame with one row per feature: the two group means, `d`,
#'   `p_value` and `fdr`, ordered by decreasing `|d|`.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs[1:200, , ])
#' head(compare_groups(x, sim$slides$subtype[1:200]), 3)
#' @export
compare_groups <- function(x, y) {
  y <- factor(y)
  a <- y == levels(y)[1]
  b <- !a
  res <- lapply(seq_len(ncol(x)), function(j) {
    u <- x[a, j]; v <- x[b, j]
    sp <- sqrt(((length(u) - 1) * stats::var(u) +
                  (length(v) - 1) * stats::var(v)) /
                 (length(u) + length(v) - 2))
    constant <- !is.finite(sp) || sp == 0
    p <- if (constant) NA_real_ else stats::t.test(v, u)$p.value
    data.frame(feature = colnames(x)[j], mean1 = mean(u), mean2 = mean(v),
               d = if (constant) NA_real_ else (mean(v) - mean(u)) / sp,
               p_value = p)
  })
  res <- do.call(rbind, res)
  names(res)[2:3] <- paste0("mean_", levels(y))
  res$fdr <- stats::p.adjust(res$p_value, method = "BH")
  res <- res[order(-abs(res$d)), ]
  rownames(res) <- NULL
  res
}
