# Baselines from measured bulk RNA-seq of the same TCGA tumours, evaluated on
# the patients of the primary test set that have RNA-seq.
#   (i)  five canonical markers: TP63, SOX2, KRT5 (squamous), NKX2-1, NAPSA
#        (adenocarcinoma)
#   (ii) the 2,000 most variable genes of the training patients
# The graph-based models are re-evaluated on the same subset of test patients.
#
# Usage: Rscript 03_bulk_baselines.R <output folder> <folder with
#        raw_counts_{luad,lusc}.csv>
# The count files are STAR gene counts from the GDC (TCGA-LUAD, TCGA-LUSC),
# genes in rows, RNA-seq aliquot barcodes in columns. They are not distributed
# with the package.

suppressMessages(library(pathwaydepth))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; rna_dir <- args[2]
prim <- readRDS(file.path(out_dir, "primary.rds"))

read_counts <- function(f) {
  d <- read.csv(f, check.names = FALSE, row.names = 1)
  as.matrix(d)
}
luad <- read_counts(file.path(rna_dir, "raw_counts_luad.csv"))
lusc <- read_counts(file.path(rna_dir, "raw_counts_lusc.csv"))
stopifnot(identical(rownames(luad), rownames(lusc)))
counts <- cbind(luad, lusc)
# one primary-tumour aliquot per patient (first in barcode order)
counts <- counts[, substr(colnames(counts), 14, 15) == "01"]
counts <- counts[, order(colnames(counts))]
counts <- counts[, !duplicated(substr(colnames(counts), 1, 12))]
colnames(counts) <- substr(colnames(counts), 1, 12)
logcpm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)

pat <- prim$slides$patient
has <- pat %in% colnames(logcpm)
tr <- prim$train & has; te <- !prim$train & has
y <- prim$y
E <- t(logcpm[, pat[has]])                 # patients x genes, order of pat[has]
rownames(E) <- pat[has]
Etr <- E[pat[tr], ]; Ete <- E[pat[te], ]
cat("patients with RNA-seq: train", sum(tr), " test", sum(te), "\n")
print(table(y[te]))

markers <- c("TP63", "SOX2", "KRT5", "NKX2-1", "NAPSA")
stopifnot(all(markers %in% colnames(E)))
fit_m <- fit_forest(Etr[, markers], y[tr], seed = 123)
p_m <- predict_forest(fit_m, Ete[, markers])
glm_m <- stats::glm(y[tr] ~ ., data = as.data.frame(Etr[, markers]), family = "binomial")
p_g <- stats::predict(glm_m, as.data.frame(Ete[, markers]), type = "response")

v <- apply(Etr, 2, stats::var)
top <- names(sort(v, decreasing = TRUE))[1:2000]
fit_v <- fit_forest(Etr[, top], y[tr], seed = 123, mtry = c(44, 200))
p_v <- predict_forest(fit_v, Ete[, top])

sub <- which(!prim$train) %in% which(te)   # rows of prim$probs that have RNA
rows <- list(
  "Bulk RNA-seq, 5 markers (logistic)" = p_g,
  "Bulk RNA-seq, 5 markers (forest)" = p_m,
  "Bulk RNA-seq, 2000 variable genes (forest)" = p_v,
  "Graph: raw edge weights" = prim$probs[sub, "raw"],
  "Graph: population depth" = prim$probs[sub, "popdepth"],
  "Graph: within-sample depth" = prim$probs[sub, "within"])
tab <- do.call(rbind, lapply(names(rows), function(k) {
  e <- evaluate(rows[[k]], y[te], 0.5)
  data.frame(model = k, n_test = sum(te), auc = e$estimate[1],
             lower = e$lower[1], upper = e$upper[1])
}))
write.csv(tab, file.path(out_dir, "bulk_baselines.csv"), row.names = FALSE)
print(tab, digits = 3)
print(summary(glm_m)$coefficients, digits = 3)
