# Pathway pairs distinguished by one representation only, and what the
# difference is: a shift in the pair's position within the slide (signed
# within-sample rank), or a difference in spread across patients
# (Fligner-Killeen test on the edge weights). Training patients only.
#
# Usage: Rscript 15_unique_pairs.R <output folder>
#   needs primary.rds and features_*.csv from 01_primary.R

suppressMessages(library(pathwaydepth))
out_dir <- commandArgs(trailingOnly = TRUE)[1]
pr <- readRDS(file.path(out_dir, "primary.rds"))
y <- pr$y[pr$train]
x <- pr$feats$raw$tr
rd <- function(k) utils::read.csv(file.path(out_dir, paste0("features_", k, ".csv")))
w <- rd("raw"); pat <- rd("within"); pop <- rd("popdepth")
rk <- compare_groups(pr$feats$rank$tr, y)
f <- colnames(x)
get <- function(d, col) d[[col]][match(f, d$feature)]
res <- data.frame(
  feature = f,
  weight_d = get(w, "d"), weight_fdr = get(w, "fdr"),
  patient_depth_d = get(pat, "d"), patient_depth_fdr = get(pat, "fdr"),
  population_depth_d = get(pop, "d"), population_depth_fdr = get(pop, "fdr"),
  signed_rank_d = get(rk, "d"), signed_rank_fdr = get(rk, "fdr"),
  median_LUSC = apply(x[y == "LUSC", ], 2, stats::median),
  median_LUAD = apply(x[y == "LUAD", ], 2, stats::median),
  sd_LUSC = apply(x[y == "LUSC", ], 2, stats::sd),
  sd_LUAD = apply(x[y == "LUAD", ], 2, stats::sd),
  fligner_p = sapply(f, function(j) stats::fligner.test(list(x[y == "LUSC", j], x[y == "LUAD", j]))$p.value))
res$fligner_fdr <- stats::p.adjust(res$fligner_p, "BH")
s <- function(v) !is.na(v) & v < 0.05
res$pattern <- paste0(ifelse(s(res$weight_fdr), "W", "-"), ifelse(s(res$patient_depth_fdr), "P", "-"),
                      ifelse(s(res$population_depth_fdr), "C", "-"))
res <- res[order(res$pattern, res$patient_depth_fdr), ]
utils::write.csv(res, file.path(out_dir, "fusion_unique_pairs.csv"), row.names = FALSE)
print(res[res$pattern %in% c("-P-", "--C", "-PC"), ], digits = 3)
cat("pairs with different spread (Fligner FDR < 0.05):", sum(res$fligner_fdr < 0.05), "\n")
