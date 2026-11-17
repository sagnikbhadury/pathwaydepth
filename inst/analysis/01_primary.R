# Primary analysis: one slide per patient (900 patients), patient-level
# stratified 70/30 split, five feature sets, test-set intervals and paired
# comparisons.
#
# Usage: Rscript 01_primary.R <output folder> <folder of graph CSV files> <annotation CSV>
# The graphs are not distributed with the package; see ?read_graphs.

suppressMessages(library(pathwaydepth))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_graphs(args[2], args[3])
lung_graphs <- dat$graphs; lung_slides <- dat$slides
saveRDS(dat, file.path(out_dir, "graphs.rds"))   # local copy for the later scripts
keep <- lung_slides$primary
graphs <- lung_graphs[keep, , ]
slides <- lung_slides[keep, ]
y <- slides$subtype                       # levels LUSC, LUAD; LUAD is positive

## ---- description of the graphs -------------------------------------------
present <- presence_matrix(graphs)
presence_tab <- data.frame(
  pathway = colnames(present),
  absent_LUSC = colSums(!present[y == "LUSC", ]),
  absent_LUAD = colSums(!present[y == "LUAD", ]),
  n_LUSC = sum(y == "LUSC"), n_LUAD = sum(y == "LUAD"))
presence_tab$p_fisher <- sapply(seq_len(ncol(present)), function(j) {
  if (length(unique(present[, j])) < 2) return(NA_real_)   # never absent
  stats::fisher.test(table(present[, j], y))$p.value
})
presence_tab$fdr <- p.adjust(presence_tab$p_fisher, "BH")
write.csv(presence_tab, file.path(out_dir, "node_absence_by_subtype.csv"),
          row.names = FALSE)

x_raw <- edge_matrix(graphs)
x_sym <- edge_matrix(graphs, complete = FALSE)   # symmetrised, no completion
imputation <- c(
  slides = nrow(x_raw),
  slides_with_absent_node = sum(rowSums(!present) > 0),
  median_absent_nodes = stats::median(rowSums(!present)),
  edges_total = length(x_raw),
  edges_touching_absent_node = sum(x_sym == 0 & {
    ut <- which(upper.tri(matrix(0, 14, 14)), arr.ind = TRUE)
    a <- !present[, ut[, 1]] | !present[, ut[, 2]]
    a
  }),
  edges_changed_by_completion = sum(abs(x_raw - x_sym) > 1e-12),
  max_abs_change = max(abs(x_raw - x_sym)),
  edges_zero_after_completion = sum(x_raw == 0))
write.csv(data.frame(quantity = names(imputation), value = imputation),
          file.path(out_dir, "imputation_summary.csv"), row.names = FALSE)

## ---- split ----------------------------------------------------------------
set.seed(20261117)
train <- grouped_split(slides$patient, y, p = 0.7)
stopifnot(!any(slides$patient[train] %in% slides$patient[!train]))
write.csv(as.data.frame(table(set = ifelse(train, "train", "test"), subtype = y)),
          file.path(out_dir, "split_summary.csv"), row.names = FALSE)

## ---- features --------------------------------------------------------------
set.seed(1991)
feats <- list(
  raw      = list(tr = x_raw[train, ], te = x_raw[!train, ]),
  popdepth = list(tr = population_depth(x_raw[train, ], x_raw[train, ]),
                  te = population_depth(x_raw[!train, ], x_raw[train, ])),
  within   = { w <- within_depth(x_raw); list(tr = w[train, ], te = w[!train, ]) },
  rank     = { r <- within_rank(x_raw);  list(tr = r[train, ], te = r[!train, ]) },
  presence = list(tr = present[train, ] + 0, te = present[!train, ] + 0))

## ---- models ----------------------------------------------------------------
fits <- lapply(feats, function(f) fit_forest(f$tr, y[train], seed = 123))
probs <- mapply(function(fit, f) predict_forest(fit, f$te), fits, feats)
perf <- do.call(rbind, lapply(names(fits), function(k) {
  e <- evaluate(probs[, k], y[!train], fits[[k]]$threshold)
  cbind(features = k, n_features = ncol(feats[[k]]$tr),
        mtry = fits[[k]]$tune$mtry, min_node = fits[[k]]$tune$min_node_size,
        oob_auc = fits[[k]]$oob_auc, threshold = fits[[k]]$threshold, e)
}))
write.csv(perf, file.path(out_dir, "performance_primary.csv"), row.names = FALSE)

pairs <- list(c("within", "popdepth"), c("raw", "within"), c("raw", "popdepth"),
              c("rank", "within"), c("within", "presence"),
              c("popdepth", "presence"), c("raw", "presence"))
cmp <- do.call(rbind, lapply(pairs, function(p)
  cbind(model1 = p[1], model2 = p[2],
        compare_auc(probs[, p[1]], probs[, p[2]], y[!train]))))
write.csv(cmp, file.path(out_dir, "paired_comparisons_primary.csv"),
          row.names = FALSE)

## ---- feature-level description (training set only) -------------------------
ut <- which(upper.tri(matrix(0, 14, 14)), arr.ind = TRUE)
edge_absent <- !present[, ut[, 1]] | !present[, ut[, 2]]
colnames(edge_absent) <- colnames(x_raw)
describe <- function(key) {
  d <- compare_groups(feats[[key]]$tr, y[train])
  imp <- fits[[key]]$forest$variable.importance
  d$importance <- imp[d$feature]
  d$importance_rank <- rank(-imp)[d$feature]
  raw <- compare_groups(x_raw[train, ], y[train])
  d$raw_mean_LUSC <- raw$mean_LUSC[match(d$feature, raw$feature)]
  d$raw_mean_LUAD <- raw$mean_LUAD[match(d$feature, raw$feature)]
  d$raw_d <- raw$d[match(d$feature, raw$feature)]
  d$absent_LUSC <- colMeans(edge_absent[train & y == "LUSC", ])[d$feature]
  d$absent_LUAD <- colMeans(edge_absent[train & y == "LUAD", ])[d$feature]
  d
}
for (k in c("popdepth", "within", "raw")) {
  write.csv(describe(k), file.path(out_dir, paste0("features_", k, ".csv")),
            row.names = FALSE)
}

saveRDS(list(train = train, slides = slides, y = y, feats = feats, fits = fits,
             probs = probs, perf = perf, cmp = cmp, present = present),
        file.path(out_dir, "primary.rds"))

cat("\nTraining:", table(y[train]), " Test:", table(y[!train]), "\n\n")
print(presence_tab, digits = 3)
print(imputation)
print(perf[, c("features", "oob_auc", "threshold", "metric", "estimate", "lower", "upper")],
      digits = 3)
print(cmp, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
