# Two checks requested in review.
#  (a) Rule co-occurrence: is a feature's co-occurrence total in root-to-leaf
#      rules explained by how often the forest splits on it?
#  (b) The depth implemented here equals fda.usc::depth.FM on single-point
#      curves (the function the first submission called).
#
# Usage: Rscript 04_rules_and_checks.R <output folder>

suppressMessages(library(pathwaydepth))
out_dir <- commandArgs(trailingOnly = TRUE)[1]
prim <- readRDS(file.path(out_dir, "primary.rds"))

rule_stats <- function(fit, n_trees = 20) {
  feats <- fit$forest$forest$independent.variable.names
  split_n <- stats::setNames(numeric(length(feats)), feats)
  cooc <- matrix(0, length(feats), length(feats), dimnames = list(feats, feats))
  n_rules <- 0; len <- c(); leaf_n <- c()
  for (k in seq_len(n_trees)) {
    ti <- ranger::treeInfo(fit$forest, k)
    tab <- table(ti$splitvarName[!ti$terminal])
    split_n[names(tab)] <- split_n[names(tab)] + tab
    walk <- function(node, path) {
      r <- ti[ti$nodeID == node, ]
      if (r$terminal) {
        n_rules <<- n_rules + 1; len <<- c(len, length(path))
        u <- unique(path)
        if (length(u) > 1) cooc[u, u] <<- cooc[u, u] + 1
        return(invisible())
      }
      walk(r$leftChild, c(path, r$splitvarName))
      walk(r$rightChild, c(path, r$splitvarName))
    }
    walk(0, character(0))
  }
  diag(cooc) <- 0
  data.frame(feature = feats, splits = split_n, cooc_total = rowSums(cooc),
             cooc_degree = rowSums(cooc > 0),
             importance = fit$forest$variable.importance[feats],
             n_rules = n_rules, median_rule_length = stats::median(len))
}

res <- lapply(c(popdepth = "popdepth", within = "within"), function(k) {
  s <- rule_stats(prim$fits[[k]])
  s$model <- k
  s
})
all <- do.call(rbind, res)
write.csv(all, file.path(out_dir, "rule_cooccurrence_vs_splits.csv"), row.names = FALSE)
for (k in names(res)) {
  s <- res[[k]]
  cat(sprintf("%s: rules=%d, median length=%g, median training cases per rule=%.1f\n",
              k, s$n_rules[1], s$median_rule_length[1],
              sum(prim$train) * 0.632 / (s$n_rules[1] / 20)))
  cat(sprintf("  Spearman(co-occurrence total, split count) = %.3f\n",
              stats::cor(s$cooc_total, s$splits, method = "spearman")))
  cat(sprintf("  Spearman(co-occurrence total, permutation importance) = %.3f\n",
              stats::cor(s$cooc_total, s$importance, method = "spearman")))
  cat(sprintf("  Spearman(co-occurrence degree, split count) = %.3f\n",
              stats::cor(s$cooc_degree, s$splits, method = "spearman")))
}

## (b) agreement with fda.usc
if (requireNamespace("fda.usc", quietly = TRUE)) {
  set.seed(3)
  x <- prim$feats$raw$tr[, "JAKSTAT_TNFa"]
  fd <- fda.usc::fdata(matrix(x, ncol = 1), argvals = 0, rangeval = c(0, 0))
  a <- fda.usc::depth.FM(fd, fdataori = fd, trim = 0.25)$dep
  b <- fda.usc::depth.FM(fd, fdataori = fd, trim = 0.15)$dep
  ours <- halfspace_depth(x, x)
  fm1 <- 1 - abs(0.5 - stats::ecdf(x)(x))
  cat("\nmax |fda.usc depth.FM - halfspace_depth| =", max(abs(a - ours)), "\n")
  cat("max |trim 0.25 - trim 0.15| =", max(abs(a - b)), "\n")
  cat("max |halfspace_depth - 2 * (FM1 - 1/2)| =", max(abs(ours - 2 * (fm1 - 0.5))), "\n")
  cat("range of depth:", range(ours), "\n")
}

## (c) coverage of the rules of the first 20 trees: distinct rules, longest
##     rule, features used, and how many of the 50 most important features
##     (permutation importance of the whole forest) occur in the rules
cov <- do.call(rbind, lapply(c("popdepth", "within"), function(k) {
  fit <- prim$fits[[k]]
  r <- extract_rules(fit, trees = 1:20)
  key <- vapply(seq_len(nrow(r)), function(i) paste(
    r$features[[i]], format(r$lower[[i]], digits = 17),
    format(r$upper[[i]], digits = 17), r$prediction[i], collapse = "|"), "")
  imp <- sort(fit$forest$variable.importance, decreasing = TRUE)
  used <- unique(unlist(r$features))
  data.frame(model = k, trees_total = fit$forest$num.trees, trees_used = 20,
             n_rules = nrow(r), n_distinct = length(unique(key)),
             max_conditions = max(r$n_conditions),
             n_features = length(imp), n_features_used = length(used),
             top50_used = sum(names(imp)[1:50] %in% used))
}))
write.csv(cov, file.path(out_dir, "rule_coverage.csv"), row.names = FALSE)
print(cov)
