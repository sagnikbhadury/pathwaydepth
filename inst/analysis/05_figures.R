# Figures of the revised manuscript, from the saved outputs of scripts 01-04.
#
# Usage: Rscript 05_figures.R <output folder> <figure folder>

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; fig_dir <- args[2]
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- normalizePath(fig_dir)   # the cairo device needs an absolute path
prim <- readRDS(file.path(out_dir, "primary.rds"))

# One colour per feature set, fixed across panels (colour-blind checked).
col <- c(raw = "#2a78d6", popdepth = "#eb6834", within = "#1baf7a",
         presence = "#eda100", rank = "#6250d6")
lab <- c(raw = "Edge weights", popdepth = "Population Referenced Depth",
         within = "Patient Referenced Depth", presence = "Pathway presence only",
         rank = "Within-sample rank (signed)")
ink <- "#0b0b0b"; muted <- "#52514e"; grid <- "#e6e5e1"
sub_col <- c(LUSC = "#6250d6", LUAD = "#e34948")

panel <- function(letter) mtext(letter, side = 3, adj = -0.12, line = 0.6, font = 2, cex = 1)
axes <- function() { box(col = grid); axis(1, col = grid, col.axis = muted, cex.axis = 0.8)
  axis(2, col = grid, col.axis = muted, cex.axis = 0.8, las = 1) }

## ---- Figure 2: classification ------------------------------------------------
y <- prim$y[!prim$train]
perf <- read.csv(file.path(out_dir, "performance_primary.csv"))
rep_s <- read.csv(file.path(out_dir, "repeated_summary_primary.csv"))
rep_a <- read.csv(file.path(out_dir, "repeated_summary_all.csv"))
bulk <- read.csv(file.path(out_dir, "bulk_baselines.csv"))

grDevices::cairo_pdf(file.path(fig_dir, "FigureBaselines_Classification.pdf"), width = 10.5, height = 3.9)
layout(matrix(1:3, 1), widths = c(1, 1.15, 1.25))
par(mar = c(4, 4, 2.2, 1), mgp = c(2.4, 0.6, 0), tcl = -0.25, col.lab = ink)
# A: ROC curves on the held-out patients
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - specificity",
     ylab = "Sensitivity", axes = FALSE, asp = 1)
axes(); abline(0, 1, col = grid, lty = 2)
keys <- c("raw", "popdepth", "within", "presence")
for (k in keys) {
  r <- pROC::roc(y, prim$probs[, k], levels = levels(y), direction = "<", quiet = TRUE)
  lines(1 - r$specificities, r$sensitivities, col = col[k], lwd = 2)
}
auc <- perf$estimate[perf$metric == "AUC"][match(keys, perf$features[perf$metric == "AUC"])]
legend("bottomright", bty = "n", cex = 0.72, lwd = 2, col = col[keys], text.col = ink,
       legend = sprintf("%s (%.2f)", lab[keys], auc), title = "Features (test AUC)",
       title.col = muted)
panel("A")
# B: test AUC over 50 repeated patient-level splits
ord <- c("raw", "rank", "within", "popdepth", "presence")
par(mar = c(4, 11.5, 2.2, 1))
plot(NA, xlim = c(0.55, 0.9), ylim = c(0.5, length(ord) + 0.5), yaxt = "n",
     xlab = "Test AUC over 50 patient-level splits", ylab = "", axes = FALSE)
box(col = grid); axis(1, col = grid, col.axis = muted, cex.axis = 0.8)
abline(v = seq(0.6, 0.9, 0.1), col = grid)
for (i in seq_along(ord)) {
  k <- ord[i]; yy <- length(ord) + 1 - i
  s <- rep_s[rep_s$features == k, ]; a <- rep_a[rep_a$features == k, ]
  segments(s$q025, yy + 0.12, s$q975, yy + 0.12, col = col[k], lwd = 2)
  points(s$mean, yy + 0.12, pch = 16, col = col[k], cex = 1.3)
  segments(a$q025, yy - 0.16, a$q975, yy - 0.16, col = col[k], lwd = 2, lty = 3)
  points(a$mean, yy - 0.16, pch = 21, col = col[k], bg = "white", cex = 1.2, lwd = 1.5)
  mtext(lab[k], side = 2, at = yy, las = 1, line = 0.5, cex = 0.62, col = ink)
}
legend("bottomright", bty = "n", cex = 0.68, pch = c(16, 21), lty = c(1, 3), col = muted,
       pt.bg = "white", text.col = ink,
       legend = c("900 patients, one slide each", "996 slides, split by patient"))
panel("B")
# C: the same test patients with measured bulk RNA-seq
b <- bulk[c(3, 1, 4, 6, 5), ]
bl <- c("Bulk RNA-seq, 2,000 genes", "Bulk RNA-seq, 5 markers", "Edge weights",
        "Patient Referenced Depth", "Population Referenced Depth")
bc <- c(muted, muted, col["raw"], col["within"], col["popdepth"])
plot(NA, xlim = c(0.6, 1), ylim = c(0.5, 5.5), yaxt = "n", axes = FALSE,
     xlab = sprintf("Test AUC, %d patients with RNA-seq", b$n_test[1]), ylab = "")
box(col = grid); axis(1, col = grid, col.axis = muted, cex.axis = 0.8)
abline(v = seq(0.6, 1, 0.1), col = grid)
for (i in 1:5) {
  yy <- 6 - i
  segments(b$lower[i], yy, b$upper[i], yy, col = bc[i], lwd = 2)
  points(b$auc[i], yy, pch = 16, col = bc[i], cex = 1.3)
  mtext(bl[i], side = 2, at = yy, las = 1, line = 0.5, cex = 0.62, col = ink)
}
panel("C")
dev.off()

## ---- Figure 3: what the depth features measure -----------------------------
fw <- read.csv(file.path(out_dir, "features_within.csv"))
fp <- read.csv(file.path(out_dir, "features_popdepth.csv"))
ab <- read.csv(file.path(out_dir, "node_absence_by_subtype.csv"))
pretty_edge <- function(s) {
  s <- gsub("_", " - ", s); s <- gsub("JAKSTAT", "JAK-STAT", s)
  s <- gsub("TNFa", "TNF\u03b1", s); s <- gsub("TGFb", "TGF\u03b2", s)
  s <- gsub("NFkB", "NF\u03baB", s); s <- gsub("Trail", "TRAIL", s); s
}
scatter <- function(f, ylab, letter, n_lab) {
  sig <- !is.na(f$fdr) & f$fdr < 0.05
  lim <- c(-0.48, 0.48)
  plot(f$raw_d, f$d, xlim = lim, ylim = lim, pch = 16, cex = 0.9,
       col = ifelse(sig, ink, "#b9b8b2"), axes = FALSE,
       xlab = "Edge weight: Cohen's d (LUAD - LUSC)", ylab = ylab)
  axes(); abline(h = 0, v = 0, col = grid)
  # label the largest differences; alternate above and below so labels do not collide
  top <- order(-abs(f$d))[seq_len(n_lab)]
  top <- top[order(f$raw_d[top])]
  text(f$raw_d[top], f$d[top], pretty_edge(f$feature[top]), cex = 0.6, col = muted,
       pos = rep(c(3, 1), length.out = length(top)), offset = 0.4, xpd = NA)
  r <- stats::cor(f$raw_d, f$d, use = "complete.obs")
  mtext(sprintf("r = %.2f", r), side = 3, adj = 1, line = 0.2, cex = 0.7, col = muted)
  panel(letter)
}
grDevices::cairo_pdf(file.path(fig_dir, "FigureDepthVsWeight.pdf"), width = 10.5, height = 4)
layout(matrix(1:3, 1), widths = c(1, 1, 1.05))
par(mar = c(4, 4.2, 2.2, 1), mgp = c(2.5, 0.6, 0), tcl = -0.25, col.lab = ink)
scatter(fw, "Patient Referenced Depth: Cohen's d", "A", 4)
scatter(fp, "Population Referenced Depth: Cohen's d", "B", 1)
ab$LUSC <- ab$absent_LUSC / ab$n_LUSC; ab$LUAD <- ab$absent_LUAD / ab$n_LUAD
ab <- ab[order(ab$LUSC + ab$LUAD), ]
par(mar = c(4, 5.2, 2.2, 1))
plot(NA, xlim = c(0, 1), ylim = c(0.5, nrow(ab) + 0.5), yaxt = "n", axes = FALSE,
     xlab = "Share of patients with no spot assigned to the pathway", ylab = "")
box(col = grid); axis(1, col = grid, col.axis = muted, cex.axis = 0.8)
abline(v = seq(0, 1, 0.25), col = grid)
for (i in seq_len(nrow(ab))) {
  segments(ab$LUSC[i], i, ab$LUAD[i], i, col = "#b9b8b2", lwd = 2)
  points(ab$LUSC[i], i, pch = 16, col = sub_col["LUSC"], cex = 1.2)
  points(ab$LUAD[i], i, pch = 16, col = sub_col["LUAD"], cex = 1.2)
  star <- if (!is.na(ab$fdr[i]) && ab$fdr[i] < 0.05) " *" else ""
  mtext(paste0(pretty_edge(ab$pathway[i]), star), side = 2, at = i, las = 1, line = 0.5, cex = 0.62, col = ink)
}
legend("bottomright", bty = "n", pch = 16, col = sub_col, legend = names(sub_col),
       cex = 0.75, text.col = ink)
panel("C")
dev.off()
cat("figures written to", fig_dir, "\n")
