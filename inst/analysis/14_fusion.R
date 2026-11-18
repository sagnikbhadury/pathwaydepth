# Fusion of the three graph representations (edge weights, Patient Referenced
# Depth, Population Referenced Depth).
#
# Two questions: (1) does combining the representations change the separation
# of the subtypes, and (2) do the representations carry complementary
# information (different patients classified correctly, different features
# used, different pathway pairs distinguished)?
#
# Early fusion: one forest on the concatenated features.
# Late fusion: (a) the mean of the three forests' probabilities; (b) stacking,
#   a logistic regression fitted to the out-of-bag probabilities of the
#   training patients and applied to the test probabilities.
#
# Usage: Rscript 14_fusion.R <output folder> <graphs.rds> <R>
#   The output folder must hold primary.rds (from 01_primary.R) and
#   features_*.csv. graphs.rds is the local copy written by 01_primary.R.

suppressMessages(library(pathwaydepth))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; graphs_file <- args[2]; R <- as.integer(args[3])

auc_ci <- function(p, y) {
  roc <- pROC::roc(y, p, levels = levels(y), direction = "<", quiet = TRUE)
  a <- as.numeric(pROC::auc(roc))
  se <- sqrt(pROC::var(roc, method = "delong"))
  c(auc = a, se = se, p_vs_half = stats::pnorm((a - 0.5) / se, lower.tail = FALSE))
}
youden <- function(p, y) {
  roc <- pROC::roc(y, p, levels = levels(y), direction = "<", quiet = TRUE)
  pROC::coords(roc, "best", ret = "threshold", best.method = "youden")$threshold[1]
}
prefix <- function(m, s) { colnames(m) <- paste0(s, colnames(m)); m }
reps <- c(raw = "w:", within = "pat:", popdepth = "pop:")
combos <- list(raw_within = c("raw", "within"), raw_popdepth = c("raw", "popdepth"),
               within_popdepth = c("within", "popdepth"),
               all3 = c("raw", "within", "popdepth"))

## ---- primary split ------------------------------------------------------
pr <- readRDS(file.path(out_dir, "primary.rds"))
y <- pr$y; train <- pr$train; yte <- y[!train]; ytr <- y[train]
feats <- pr$feats; fits <- pr$fits; probs <- pr$probs

early <- lapply(combos, function(k) {
  tr <- do.call(cbind, lapply(k, function(r) prefix(feats[[r]]$tr, reps[[r]])))
  te <- do.call(cbind, lapply(k, function(r) prefix(feats[[r]]$te, reps[[r]])))
  fit <- fit_forest(tr, ytr, seed = 123)
  list(fit = fit, prob = predict_forest(fit, te))
})
k3 <- c("raw", "within", "popdepth")
oob <- sapply(k3, function(k) fits[[k]]$oob_prob)
late_mean <- rowMeans(probs[, k3])
thr_mean <- youden(rowMeans(oob), ytr)
stack <- stats::glm(ytr ~ ., data = data.frame(qlogis(pmin(pmax(oob, 1e-3), 1 - 1e-3))),
                    family = stats::binomial())
late_stack <- stats::predict(stack, newdata = data.frame(qlogis(pmin(pmax(probs[, k3], 1e-3), 1 - 1e-3))),
                             type = "response")
thr_stack <- youden(stats::fitted(stack), ytr)

all_prob <- cbind(probs[, k3], sapply(early, `[[`, "prob"), late_mean = late_mean,
                  late_stack = late_stack)
thr <- c(sapply(k3, function(k) fits[[k]]$threshold),
         sapply(early, function(e) e$fit$threshold), late_mean = thr_mean,
         late_stack = thr_stack)
perf <- do.call(rbind, lapply(colnames(all_prob), function(k) {
  e <- evaluate(all_prob[, k], yte, thr[[k]])
  a <- auc_ci(all_prob[, k], yte)
  vs <- if (k == "raw") c(NA, NA, NA, NA) else
    unlist(compare_auc(all_prob[, k], all_prob[, "raw"], yte)[c("difference", "lower", "upper", "p_value")])
  data.frame(model = k, auc = e$estimate[1], lower = e$lower[1], upper = e$upper[1],
             p_vs_half = a[["p_vs_half"]], accuracy = e$estimate[2],
             sensitivity = e$estimate[3], specificity = e$estimate[4],
             diff_vs_raw = vs[1], diff_lower = vs[2], diff_upper = vs[3], diff_p = vs[4])
}))
write.csv(perf, file.path(out_dir, "fusion_performance_primary.csv"), row.names = FALSE)
print(perf, digits = 3)

# complementarity of the predictions
cor_s <- stats::cor(probs[, k3], method = "spearman")
correct <- sapply(k3, function(k) (probs[, k] >= fits[[k]]$threshold) == (yte == levels(y)[2]))
pattern <- apply(correct, 1, function(z) paste(ifelse(z, "1", "0"), collapse = ""))
pat_tab <- as.data.frame(table(pattern = pattern))   # order raw, within, popdepth
pat_tab$meaning <- c("000" = "none correct", "001" = "Population depth only",
                     "010" = "Patient depth only", "011" = "both depths, not weights",
                     "100" = "edge weights only", "101" = "weights and Population depth",
                     "110" = "weights and Patient depth", "111" = "all three correct")[as.character(pat_tab$pattern)]
mcn <- function(a, b) {
  t <- table(factor(correct[, a], c(FALSE, TRUE)), factor(correct[, b], c(FALSE, TRUE)))
  data.frame(model1 = a, model2 = b, only1 = t["TRUE", "FALSE"], only2 = t["FALSE", "TRUE"],
             p_mcnemar = stats::mcnemar.test(t)$p.value)
}
mc <- rbind(mcn("raw", "within"), mcn("raw", "popdepth"), mcn("within", "popdepth"))
write.csv(data.frame(model = rownames(cor_s), cor_s, row.names = NULL),
          file.path(out_dir, "fusion_prediction_correlation.csv"), row.names = FALSE)
write.csv(pat_tab, file.path(out_dir, "fusion_correct_patterns.csv"), row.names = FALSE)
write.csv(mc, file.path(out_dir, "fusion_mcnemar.csv"), row.names = FALSE)
print(cor_s, digits = 3); print(pat_tab); print(mc, digits = 3)

# what the fused forest uses
imp <- early$all3$fit$forest$variable.importance
rep_of <- sub(":.*", "", names(imp))
pos <- pmax(imp, 0)
share <- data.frame(representation = c("w", "pat", "pop"),
                    importance_share = sapply(c("w", "pat", "pop"), function(r) sum(pos[rep_of == r]) / sum(pos)),
                    in_top15 = sapply(c("w", "pat", "pop"), function(r) sum(rep_of[order(-imp)][1:15] == r)),
                    in_top30 = sapply(c("w", "pat", "pop"), function(r) sum(rep_of[order(-imp)][1:30] == r)))
top <- data.frame(rank = 1:30, feature = names(imp)[order(-imp)][1:30], importance = sort(imp, decreasing = TRUE)[1:30])
write.csv(share, file.path(out_dir, "fusion_importance_share.csv"), row.names = FALSE)
write.csv(top, file.path(out_dir, "fusion_importance_top.csv"), row.names = FALSE)
print(share, digits = 3); print(head(top, 15), digits = 3)

# which pathway pairs each representation distinguishes (training set, FDR 0.05)
sig <- lapply(c(raw = "raw", within = "within", popdepth = "popdepth"), function(k) {
  d <- utils::read.csv(file.path(out_dir, paste0("features_", k, ".csv")))
  d$feature[!is.na(d$fdr) & d$fdr < 0.05]
})
allf <- colnames(feats$raw$tr)
venn <- data.frame(feature = allf,
                   weights = allf %in% sig$raw, patient_depth = allf %in% sig$within,
                   population_depth = allf %in% sig$popdepth)
venn$pattern <- paste0(ifelse(venn$weights, "W", "-"), ifelse(venn$patient_depth, "P", "-"),
                       ifelse(venn$population_depth, "C", "-"))
write.csv(venn, file.path(out_dir, "fusion_feature_overlap.csv"), row.names = FALSE)
print(table(venn$pattern))

## ---- repeated splits ----------------------------------------------------
dat <- readRDS(graphs_file)
keep <- dat$slides$primary
graphs <- dat$graphs[keep, , ]; slides <- dat$slides[keep, ]
yy <- slides$subtype
x_raw <- edge_matrix(graphs)
set.seed(7)
w_sub <- within_depth(x_raw)          # same RNG use as 02_repeated_splits.R
auc <- function(p, yv) as.numeric(pROC::auc(pROC::roc(yv, p, levels = levels(yv),
                                                     direction = "<", quiet = TRUE)))
one_split <- function(r) {
  set.seed(1000 + r)
  tr <- grouped_split(slides$patient, yy, p = 0.7)
  f <- list(raw = list(x_raw[tr, ], x_raw[!tr, ]),
            popdepth = list(population_depth(x_raw[tr, ], x_raw[tr, ]),
                            population_depth(x_raw[!tr, ], x_raw[tr, ])),
            within = list(w_sub[tr, ], w_sub[!tr, ]))
  single <- lapply(f, function(z) fit_forest(z[[1]], yy[tr], num_trees = 1000, tune_trees = 300, seed = r))
  p_single <- sapply(names(f), function(k) predict_forest(single[[k]], f[[k]][[2]]))
  o_single <- sapply(names(f), function(k) single[[k]]$oob_prob)
  out <- sapply(names(f), function(k) auc(p_single[, k], yy[!tr]))
  for (cn in names(combos)) {
    k <- combos[[cn]]
    xtr <- do.call(cbind, lapply(k, function(q) prefix(f[[q]][[1]], reps[[q]])))
    xte <- do.call(cbind, lapply(k, function(q) prefix(f[[q]][[2]], reps[[q]])))
    fit <- fit_forest(xtr, yy[tr], num_trees = 1000, tune_trees = 300, seed = r)
    out[cn] <- auc(predict_forest(fit, xte), yy[!tr])
  }
  out["late_mean"] <- auc(rowMeans(p_single[, k3]), yy[!tr])
  lg <- function(m) data.frame(qlogis(pmin(pmax(m[, k3], 1e-3), 1 - 1e-3)))
  st <- stats::glm(yy[tr] ~ ., data = lg(o_single), family = stats::binomial())
  out["late_stack"] <- auc(stats::predict(st, newdata = lg(p_single), type = "response"), yy[!tr])
  out
}
res <- t(sapply(seq_len(R), function(r) {
  a <- one_split(r)
  cat(sprintf("split %d: %s\n", r, paste(sprintf("%s=%.3f", names(a), a), collapse = " ")))
  a
}))
write.csv(res, file.path(out_dir, "fusion_repeated.csv"), row.names = FALSE)
summ <- data.frame(model = colnames(res), mean = colMeans(res), sd = apply(res, 2, stats::sd),
                   q025 = apply(res, 2, stats::quantile, 0.025),
                   q975 = apply(res, 2, stats::quantile, 0.975),
                   min = apply(res, 2, min))
dif <- do.call(rbind, lapply(setdiff(colnames(res), "raw"), function(k) {
  d <- res[, k] - res[, "raw"]
  data.frame(model = k, mean_diff_vs_raw = mean(d), q025 = stats::quantile(d, 0.025),
             q975 = stats::quantile(d, 0.975), share_above_raw = mean(d > 0))
}))
write.csv(summ, file.path(out_dir, "fusion_repeated_summary.csv"), row.names = FALSE)
write.csv(dif, file.path(out_dir, "fusion_repeated_difference.csv"), row.names = FALSE)
print(summ, digits = 3); print(dif, digits = 3)
