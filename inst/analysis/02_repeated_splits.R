# Repeated patient-level splits. For each of R stratified 70/30 splits of the
# patients, every feature set is rebuilt from the training set, a forest is
# tuned and fitted on the training set and its AUC is computed on the test set.
#
# Usage: Rscript 02_repeated_splits.R <output folder> <cohort> <R>
#   cohort = "primary" (one slide per patient) or "all" (all 996 slides, split
#   by patient; test AUC at slide level)

suppressMessages(library(pathwaydepth))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; cohort <- args[2]; R <- as.integer(args[3])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- readRDS(file.path(out_dir, "graphs.rds"))   # written by 01_primary.R
lung_graphs <- dat$graphs; lung_slides <- dat$slides
keep <- if (cohort == "primary") lung_slides$primary else rep(TRUE, nrow(lung_slides))
graphs <- lung_graphs[keep, , ]
slides <- lung_slides[keep, ]
y <- slides$subtype
x_raw <- edge_matrix(graphs)
present <- presence_matrix(graphs) + 0
set.seed(7)
w_sub <- within_depth(x_raw)              # does not depend on the split
w_full <- within_depth(x_raw, B = 0)
w_rank <- within_rank(x_raw)

auc <- function(p, y) as.numeric(pROC::auc(pROC::roc(
  y, p, levels = levels(y), direction = "<", quiet = TRUE)))

one_split <- function(r) {
  set.seed(1000 + r)
  tr <- grouped_split(slides$patient, y, p = 0.7)
  # population rank of each edge, learned on the training set, then the
  # within-sample depth of those ranks (edges put on a common scale first)
  ecdf_scale <- function(new) {
    out <- new
    for (j in seq_len(ncol(new))) out[, j] <- stats::ecdf(x_raw[tr, j])(new[, j])
    out
  }
  s_all <- within_depth(ecdf_scale(x_raw), B = 0)
  f <- list(
    raw = list(x_raw[tr, ], x_raw[!tr, ]),
    popdepth = list(population_depth(x_raw[tr, ], x_raw[tr, ]),
                    population_depth(x_raw[!tr, ], x_raw[tr, ])),
    popdepth_nosub = list(population_depth(x_raw[tr, ], x_raw[tr, ], B = 0),
                          population_depth(x_raw[!tr, ], x_raw[tr, ], B = 0)),
    within = list(w_sub[tr, ], w_sub[!tr, ]),
    within_nosub = list(w_full[tr, ], w_full[!tr, ]),
    within_scaled = list(s_all[tr, ], s_all[!tr, ]),
    rank = list(w_rank[tr, ], w_rank[!tr, ]),
    presence = list(present[tr, ], present[!tr, ]))
  sapply(f, function(z) {
    fit <- fit_forest(z[[1]], y[tr], num_trees = 1000, tune_trees = 300,
                      seed = r)
    auc(predict_forest(fit, z[[2]]), y[!tr])
  })
}

res <- t(sapply(seq_len(R), function(r) {
  a <- one_split(r)
  cat(sprintf("split %d: %s\n", r, paste(sprintf("%s=%.3f", names(a), a), collapse = " ")))
  a
}))
write.csv(res, file.path(out_dir, sprintf("repeated_splits_%s.csv", cohort)),
          row.names = FALSE)

summ <- data.frame(features = colnames(res), mean = colMeans(res),
                   sd = apply(res, 2, stats::sd),
                   q025 = apply(res, 2, stats::quantile, 0.025),
                   q975 = apply(res, 2, stats::quantile, 0.975))
dif <- function(a, b) {
  d <- res[, a] - res[, b]
  data.frame(model1 = a, model2 = b, mean_diff = mean(d), sd = stats::sd(d),
             q025 = stats::quantile(d, 0.025), q975 = stats::quantile(d, 0.975),
             share_positive = mean(d > 0))
}
dd <- rbind(dif("within", "popdepth"), dif("raw", "within"), dif("raw", "popdepth"),
            dif("rank", "within"), dif("within", "within_nosub"),
            dif("popdepth", "popdepth_nosub"), dif("within_scaled", "within"),
            dif("within", "presence"))
write.csv(summ, file.path(out_dir, sprintf("repeated_summary_%s.csv", cohort)), row.names = FALSE)
write.csv(dd, file.path(out_dir, sprintf("repeated_differences_%s.csv", cohort)), row.names = FALSE)
print(summ, digits = 3); print(dd, digits = 3)
