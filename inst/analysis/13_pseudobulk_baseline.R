# Baseline from the predicted expression itself: the pseudo-bulk of the DeepSpot
# predictions (mean over the spots of a slide, script 12), on the same
# patient-level splits as the graph-based models.
#   pb_genes    random forest on the 2,000 predicted genes with the largest
#               variance among the training patients
#   pb_markers  logistic regression on the canonical markers that DeepSpot
#               predicts (TP63, NKX2-1, NAPSA; SOX2 and KRT5 are not predicted)
#
# Usage: Rscript 13_pseudobulk_baseline.R <output folder of 01_primary.R> <pseudo-bulk .csv.gz> <splits>

suppressMessages(library(pathwaydepth))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; R <- as.integer(args[3])
prim <- readRDS(file.path(out_dir, "primary.rds"))
pb <- read.csv(args[2], check.names = FALSE, row.names = 1)
slides <- prim$slides$slide
stopifnot(all(slides %in% rownames(pb)))
E <- as.matrix(pb[slides, setdiff(colnames(pb), "n_spots")])
y <- prim$y; patient <- prim$slides$patient
markers <- intersect(c("TP63", "SOX2", "KRT5", "NKX2-1", "NAPSA"), colnames(E))
cat("slides", nrow(E), "genes", ncol(E), "markers predicted:", markers, "\n")
x_raw <- rbind(prim$feats$raw$tr, prim$feats$raw$te)[order(c(which(prim$train), which(!prim$train))), ]
auc <- function(p, yy) as.numeric(pROC::auc(pROC::roc(yy, p, levels = levels(yy), direction = "<", quiet = TRUE)))

run <- function(tr, seed, trees = 1500) {
  v <- apply(E[tr, ], 2, stats::var); top <- names(sort(v, decreasing = TRUE))[1:2000]
  f <- fit_forest(E[tr, top], y[tr], seed = seed, num_trees = trees, mtry = c(44, 200))
  p_genes <- predict_forest(f, E[!tr, top])
  d <- data.frame(yy = as.integer(y == "LUAD"), E[, markers, drop = FALSE], check.names = FALSE)
  names(d) <- make.names(names(d))
  g <- stats::glm(yy ~ ., data = d[tr, ], family = binomial())
  p_mark <- stats::predict(g, d[!tr, ], type = "response")
  fr <- fit_forest(x_raw[tr, ], y[tr], seed = seed, num_trees = trees)
  list(pb_genes = p_genes, pb_markers = p_mark, raw_rf = predict_forest(fr, x_raw[!tr, ]), glm = g)
}

## primary split
p1 <- run(prim$train, 123)
yt <- y[!prim$train]
perf <- do.call(rbind, lapply(c("pb_genes", "pb_markers", "raw_rf"), function(k) {
  e <- evaluate(p1[[k]], yt, 0.5)
  data.frame(model = k, auc = e$estimate[1], lower = e$lower[1], upper = e$upper[1])
}))
cmp <- rbind(cbind(a = "pb_genes", b = "raw_rf", compare_auc(p1$pb_genes, p1$raw_rf, yt)),
             cbind(a = "pb_markers", b = "raw_rf", compare_auc(p1$pb_markers, p1$raw_rf, yt)))
write.csv(perf, file.path(out_dir, "pseudobulk_performance_primary.csv"), row.names = FALSE)
write.csv(cmp, file.path(out_dir, "pseudobulk_comparisons_primary.csv"), row.names = FALSE)
print(perf, digits = 3); print(cmp, digits = 3); print(summary(p1$glm)$coefficients, digits = 3)

## repeated patient-level splits (same seeds as script 02)
if (R > 0) {
  res <- t(sapply(seq_len(R), function(r) {
    set.seed(1000 + r)
    tr <- grouped_split(patient, y, p = 0.7)
    p <- run(tr, r, trees = 1000)
    a <- c(pb_genes = auc(p$pb_genes, y[!tr]), pb_markers = auc(p$pb_markers, y[!tr]), raw_rf = auc(p$raw_rf, y[!tr]))
    cat(sprintf("split %d: %s\n", r, paste(sprintf("%s=%.3f", names(a), a), collapse = " ")))
    a
  }))
  write.csv(res, file.path(out_dir, "pseudobulk_repeated.csv"), row.names = FALSE)
  s <- data.frame(features = colnames(res), mean = colMeans(res), q025 = apply(res, 2, stats::quantile, 0.025),
                  q975 = apply(res, 2, stats::quantile, 0.975))
  d <- res[, "pb_genes"] - res[, "raw_rf"]
  s2 <- data.frame(mean_diff = mean(d), q025 = stats::quantile(d, 0.025), q975 = stats::quantile(d, 0.975), share_positive = mean(d > 0))
  write.csv(s, file.path(out_dir, "pseudobulk_repeated_summary.csv"), row.names = FALSE)
  write.csv(s2, file.path(out_dir, "pseudobulk_repeated_difference.csv"), row.names = FALSE)
  print(s, digits = 3); print(s2, digits = 3)
}
