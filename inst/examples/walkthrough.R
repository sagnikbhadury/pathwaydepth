# The full analysis on simulated graphs. This script produces every figure of
# the README. The data are random numbers (see ?simulate_graphs); the results
# illustrate what the package returns and say nothing about lung cancer.
#
# Usage: Rscript inst/examples/walkthrough.R <folder for the figures>

library(pathwaydepth)
fig_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(fig_dir)) fig_dir <- "man/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- normalizePath(fig_dir)
png_ <- function(name, w, h) grDevices::png(file.path(fig_dir, name), width = w, height = h,
                                            units = "in", res = 150)

## 1. Data: one slide per patient ----------------------------------------------
sim <- simulate_graphs(660, seed = 2026)
keep <- sim$slides$primary
graphs <- sim$graphs[keep, , ]
y <- sim$slides$subtype[keep]
patient <- sim$slides$patient[keep]

png_("README-presence.png", 6, 4.2)
plot_presence(graphs, y)
dev.off()

## 2. Edge weights, patient-level split, depth features -------------------------
x <- edge_matrix(graphs)                       # slides x 91 edge weights
set.seed(1)
train <- grouped_split(patient, y, p = 0.7)

set.seed(2)
features <- list(
  "Edge weights" = list(train = x[train, ], test = x[!train, ]),
  "Population depth" = list(train = population_depth(x[train, ], x[train, ]),
                            test  = population_depth(x[!train, ], x[train, ])),
  "Within-sample depth" = list(train = within_depth(x[train, ]),
                               test  = within_depth(x[!train, ])),
  "Node presence" = list(train = presence_matrix(graphs[train, , ]) + 0,
                         test  = presence_matrix(graphs[!train, , ]) + 0))

## 3. Random forests, evaluation, paired comparison ----------------------------
fits  <- lapply(features, function(f) fit_forest(f$train, y[train], seed = 1))
probs <- mapply(function(fit, f) predict_forest(fit, f$test), fits, features, SIMPLIFY = FALSE)
evals <- mapply(function(p, fit) evaluate(p, y[!train], fit$threshold), probs, fits,
                SIMPLIFY = FALSE)
print(evals[["Within-sample depth"]])
print(compare_auc(probs[["Within-sample depth"]], probs[["Population depth"]], y[!train]))

get <- function(metric, col) sapply(evals, function(e) e[e$metric == metric, col])
png_("README-performance.png", 11, 4.2)
layout(matrix(1:2, 1), widths = c(1, 1.25))
par(mar = c(4, 4, 2.2, 1))
plot_roc(probs, y[!train], main = "Held-out patients")
plot_intervals(get("AUC", "estimate"), get("AUC", "lower"), get("AUC", "upper"),
               xlab = "Test AUC with 95% interval")
dev.off()

## 4. Repeated patient-level splits --------------------------------------------
w_all <- within_depth(x)                       # does not depend on the split
res <- t(sapply(1:20, function(r) {
  set.seed(100 + r)
  tr <- grouped_split(patient, y, p = 0.7)
  f <- list("Edge weights" = list(x[tr, ], x[!tr, ]),
            "Population depth" = list(population_depth(x[tr, ], x[tr, ], B = 0),
                                      population_depth(x[!tr, ], x[tr, ], B = 0)),
            "Within-sample depth" = list(w_all[tr, ], w_all[!tr, ]))
  sapply(f, function(z) {
    fit <- fit_forest(z[[1]], y[tr], num_trees = 500, tune_trees = 200, seed = r)
    evaluate(predict_forest(fit, z[[2]]), y[!tr], fit$threshold)$estimate[1]
  })
}))
png_("README-repeated.png", 6.5, 3)
plot_intervals(colMeans(res), apply(res, 2, quantile, 0.025), apply(res, 2, quantile, 0.975),
               xlab = "Test AUC over 20 patient-level splits")
dev.off()

## 5. Which edges differ between the classes ------------------------------------
cmp_weight <- compare_groups(features[["Edge weights"]]$train, y[train])
cmp_pop    <- compare_groups(features[["Population depth"]]$train, y[train])
cmp_within <- compare_groups(features[["Within-sample depth"]]$train, y[train])
print(head(cmp_weight, 5))

png_("README-volcano.png", 11, 4.2)
par(mfrow = c(1, 2), mar = c(4, 4.2, 2.2, 1))
plot_volcano(cmp_pop, main = "Population depth")
plot_volcano(cmp_within, main = "Within-sample depth")
dev.off()

png_("README-depth-vs-weight.png", 9.5, 4.4)
par(mfrow = c(1, 2), mar = c(4, 4.2, 2.2, 1))
plot_depth_vs_weight(cmp_pop, cmp_weight, ylab = "Population depth: Cohen's d")
plot_depth_vs_weight(cmp_within, cmp_weight, ylab = "Within-sample depth: Cohen's d")
dev.off()

## 6. Feature importance ---------------------------------------------------------
png_("README-importance.png", 11, 4.6)
layout(matrix(1:2, 1))
plot_importance(fits[["Population depth"]], top = 15, col = depth_colours(2)[2],
                main = "Population depth")
plot_importance(fits[["Within-sample depth"]], top = 15, col = depth_colours(3)[3],
                main = "Within-sample depth")
dev.off()

## 7. Decision rules and their co-occurrence -------------------------------------
rules <- list("Population depth" = extract_rules(fits[["Population depth"]], trees = 1:20),
              "Within-sample depth" = extract_rules(fits[["Within-sample depth"]], trees = 1:20))
r1 <- rules[[2]][rules[[2]]$tree == 1, ]
stopifnot(all(rule_membership(r1, features[["Within-sample depth"]]$train) == 1))
print(head(rules[[2]][order(rules[[2]]$n_conditions), c("rule", "prediction")], 3))

png_("README-rules.png", 11, 3.8)
plot_rule_lengths(rules, col = depth_colours(3)[2:3])
dev.off()

co <- lapply(rules, cooccurrence)
print(head(co[[2]]$nodes, 5))
cat("Spearman correlation of co-occurrence total with number of rules using the feature:",
    round(cor(co[[2]]$nodes$total, co[[2]]$nodes$rules, method = "spearman"), 2), "\n")

png_("README-cooccurrence.png", 11, 5.6)
layout(matrix(1:2, 1))
plot_cooccurrence(co[[1]], top = 30, col = depth_colours(2)[2], main = "Population depth")
plot_cooccurrence(co[[2]], top = 30, col = depth_colours(3)[3], main = "Within-sample depth")
dev.off()
cat("figures written to", fig_dir, "\n")
