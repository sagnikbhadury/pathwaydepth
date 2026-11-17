# Regenerates the figures of the first submission (performance, differential
# analysis, co-occurrence networks, importance, decision rules) from the
# patient-level analysis, using the original figure script unchanged.
#
# Usage: Rscript 06_original_figures.R <output folder of 01_primary.R>
#          <path to integrated_figure_generation.R> <folder to write into>
# The figures are written to <folder>/figures.

suppressMessages({library(pathwaydepth); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
out_dir <- normalizePath(args[1]); fig_script <- normalizePath(args[2])
dir.create(args[3], showWarnings = FALSE, recursive = TRUE)
prim <- readRDS(file.path(out_dir, "primary.rds"))

# The original script labels LUSC "Control" and LUAD "Disease".
relabel <- function(y) factor(ifelse(y == "LUSC", "Control", "Disease"),
                              levels = c("Control", "Disease"))
y_train <- relabel(prim$y[prim$train]); y_test <- relabel(prim$y[!prim$train])

# The evaluated forests, refitted with the original labels (same data, seed and
# tuning parameters, hence the same trees).
refit <- function(key) {
  f <- prim$fits[[key]]
  w <- c(Control = 1, Disease = sum(y_train == "Control") / sum(y_train == "Disease"))
  rf <- ranger::ranger(x = prim$feats[[key]]$tr, y = y_train, num.trees = 1500,
                       mtry = f$tune$mtry, min.node.size = f$tune$min_node_size,
                       probability = TRUE, importance = "permutation",
                       class.weights = w, seed = 123)
  p_new <- predict(rf, prim$feats[[key]]$te)$predictions[, "Disease"]
  stopifnot(max(abs(p_new - prim$probs[, key])) < 1e-12)   # identical forest
  list(rf_final = rf, threshold = f$threshold)
}
fit_A <- refit("popdepth"); fit_B <- refit("within")

# Root-to-leaf rules of the first 20 trees of the evaluated forests.
rules_of <- function(fit, model) {
  r <- extract_rules(fit$rf_final, trees = 1:20)
  data.frame(tree_id = r$tree, rule_id = ave(r$tree, r$tree, FUN = seq_along),
             rule = r$rule, prediction = r$prediction,
             n_conditions = r$n_conditions, Model = model, stringsAsFactors = FALSE)
}
rules_A <- rules_of(fit_A, "Population Referenced Depth")
rules_B <- rules_of(fit_B, "Patient Referenced Depth")

# Pairwise co-occurrence of features within rules (as in the first submission).
pairs_of <- function(rules, model) {
  feats <- lapply(rules$rule, function(r)
    unique(unlist(regmatches(r, gregexpr("[A-Za-z0-9]+_[A-Za-z0-9]+", r)))))
  all <- sort(unique(unlist(feats)))
  m <- matrix(0L, length(all), length(all), dimnames = list(all, all))
  for (f in feats) if (length(f) > 1) m[f, f] <- m[f, f] + 1L
  diag(m) <- 0L
  idx <- which(upper.tri(m) & m > 0, arr.ind = TRUE)
  d <- data.frame(Feature1 = all[idx[, 1]], Feature2 = all[idx[, 2]],
                  CoOccurrence = m[idx], Model = model, stringsAsFactors = FALSE)
  list(pairs = d[order(-d$CoOccurrence), ], matrix = m)
}
pA <- pairs_of(rules_A, "Population Referenced Depth")
pB <- pairs_of(rules_B, "Patient Referenced Depth")

setwd(args[3])
source(fig_script)
output_dir <- normalizePath(output_dir)   # the cairo device needs an absolute path
# The heatmap function of the original script reads these from the global
# environment instead of from its arguments.
Depth_A_train <- prim$feats$popdepth$tr; Depth_A_test <- prim$feats$popdepth$te
Depth_B_train <- prim$feats$within$tr;   Depth_B_test <- prim$feats$within$te
results <- generate_all_figures(
  fit_depthA = fit_A, fit_depthB = fit_B,
  X_depthA_train = prim$feats$popdepth$tr, X_depthA_test = prim$feats$popdepth$te,
  X_depthB_train = prim$feats$within$tr, X_depthB_test = prim$feats$within$te,
  y_train = y_train, y_test = y_test,
  feature_pairs_depthA = pA$pairs, feature_pairs_depthB = pB$pairs,
  network_layout = "circular",
  tree_rules_df_depthA = rules_A, tree_rules_df_depthB = rules_B)

## ---- numbers quoted with these figures -------------------------------------
perf <- function(fit, key) {
  m <- compute_performance_metrics(fit, prim$feats[[key]]$te, y_test)   # threshold 0.5
  tr <- predict(fit$rf_final, prim$feats[[key]]$tr)$predictions[, "Disease"]
  data.frame(model = key, auc = m$auc, accuracy = m$accuracy,
             sensitivity = m$sensitivity, specificity = m$specificity, f1 = m$f1,
             train_auc_inbag = as.numeric(pROC::auc(pROC::roc(
               y_train, tr, levels = levels(y_train), direction = "<", quiet = TRUE))))
}
write.csv(rbind(perf(fit_A, "popdepth"), perf(fit_B, "within")),
          "figures/Figure1_metrics_threshold_0.5.csv", row.names = FALSE)

net <- function(p, rules, key) {
  m <- p$matrix
  data.frame(model = key, feature = rownames(m), degree = rowSums(m > 0),
             cooccurrence_total = rowSums(m), row.names = NULL)
}
nets <- rbind(net(pA, rules_A, "popdepth"), net(pB, rules_B, "within"))
write.csv(nets, "figures/cooccurrence_node_statistics.csv", row.names = FALSE)
summ <- function(p, rules, key) data.frame(
  model = key, rules = nrow(rules), median_conditions = median(rules$n_conditions),
  share_le3 = mean(rules$n_conditions <= 3), share_gt6 = mean(rules$n_conditions > 6),
  features_in_network = sum(rowSums(p$matrix) > 0), pairs = nrow(p$pairs),
  median_pair_count = median(p$pairs$CoOccurrence), max_pair_count = max(p$pairs$CoOccurrence),
  top_pair = paste(p$pairs$Feature1[1], p$pairs$Feature2[1], sep = " + "))
write.csv(rbind(summ(pA, rules_A, "popdepth"), summ(pB, rules_B, "within")),
          "figures/cooccurrence_summary.csv", row.names = FALSE)
write.csv(rbind(rules_A, rules_B), "figures/decision_rules_first20trees.csv", row.names = FALSE)
cat("done\n")
