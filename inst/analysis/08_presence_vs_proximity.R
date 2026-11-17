# Presence or proximity? A Bayesian multilevel decomposition of the subtype
# differences in edge weight.
#
# An edge weight is zero when one of its two pathways is absent from the slide
# (no spot is assigned to it), and often also when both are present. A
# difference in mean edge weight between subtypes can therefore arise in three
# ways, which this script separates with three multilevel models fitted to the
# training patients of the primary split (brms):
#   A  presence of a pathway on a slide       (Bernoulli, 14 pathways)
#   B  a positive weight given both pathways present   (Bernoulli, 91 pairs)
#   C  the log weight given that it is positive        (Student-t, 91 pairs)
# Each model has a subtype effect that varies by pathway (A) or by pair (B, C)
# around a common mean (partial pooling, which also handles the multiplicity
# of 91 pairs) and a patient intercept for the overall level of a slide.
#
# Usage: Rscript 08_presence_vs_proximity.R <output folder>

suppressMessages({library(pathwaydepth); library(brms)})
out_dir <- commandArgs(trailingOnly = TRUE)[1]
prim <- readRDS(file.path(out_dir, "primary.rds"))
tr <- prim$train
x <- prim$feats$raw$tr; y <- prim$y[tr]; pres <- prim$present[tr, ]
patient <- prim$slides$patient[tr]
ut <- which(upper.tri(matrix(0, 14, 14)), arr.ind = TRUE)
both <- pres[, ut[, 1]] & pres[, ut[, 2]]
tol <- 1e-6                                   # weights below this are numerical zeros of the completion step
pos <- x > tol
luad <- as.integer(y == "LUAD")

long <- data.frame(patient = rep(patient, ncol(x)), luad = rep(luad, ncol(x)),
                   edge = rep(colnames(x), each = nrow(x)), w = as.vector(x),
                   both = as.vector(both), pos = as.integer(as.vector(pos)))
dB <- long[long$both, ]
dC <- dB[dB$pos == 1, ]; dC$lw <- log(dC$w)
dA <- data.frame(patient = rep(patient, 14), luad = rep(luad, 14),
                 pathway = rep(colnames(pres), each = nrow(pres)), present = as.integer(as.vector(pres)))
cat("rows: A", nrow(dA), " B", nrow(dB), " C", nrow(dC), "\n")

fit <- function(formula, data, family, prior, name) {
  t0 <- Sys.time()
  f <- brm(formula, data = data, family = family, prior = prior, chains = 4, cores = 4,
           iter = 2000, warmup = 1000, seed = 2026, backend = "cmdstanr", threads = threading(6),
           control = list(adapt_delta = 0.9), refresh = 0, silent = 2,
           file = file.path(out_dir, paste0("brms_", name)), file_refit = "on_change")
  cat(name, "fitted in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  f
}
# Pathway- and pair-specific intercepts are unpooled (fixed); the subtype effects are partially pooled.
pr <- c(prior(normal(0, 5), class = "b"), prior(normal(0, 1), class = "b", coef = "luad"),
        prior(exponential(2), class = "sd"))
fA <- fit(present ~ 0 + pathway + luad + (0 + luad | pathway) + (1 | patient), dA, bernoulli(), pr, "A_presence")
fB <- fit(pos ~ 0 + edge + luad + (0 + luad | edge) + (1 | patient), dB, bernoulli(), pr, "B_positive")
fC <- fit(lw ~ 0 + edge + luad + (0 + luad | edge) + (1 | patient), dC, student(), pr, "C_logweight")

effects <- function(f, group, part) {
  d <- as.matrix(f, variable = "b_luad")[, 1] +
    as.matrix(f, variable = paste0("r_", group), regex = TRUE)[, grepl(",luad]", colnames(as.matrix(f, variable = paste0("r_", group), regex = TRUE)), fixed = TRUE)]
  nm <- sub(paste0("r_", group, "\\["), "", sub(",luad\\]", "", colnames(d)))
  data.frame(part = part, unit = nm, median = apply(d, 2, stats::median),
             q025 = apply(d, 2, stats::quantile, 0.025), q975 = apply(d, 2, stats::quantile, 0.975),
             p_positive = colMeans(d > 0), row.names = NULL)
}
eff <- rbind(effects(fA, "pathway", "A_presence"), effects(fB, "edge", "B_positive"), effects(fC, "edge", "C_logweight"))
eff$excludes_zero <- eff$q025 > 0 | eff$q975 < 0
raw <- compare_groups(x, y)
eff$raw_d <- raw$d[match(eff$unit, raw$feature)]
eff$raw_fdr <- raw$fdr[match(eff$unit, raw$feature)]
eff$n <- c(rep(nrow(pres), 14), as.vector(table(dB$edge)[eff$unit[eff$part == "B_positive"]]),
           as.vector(table(dC$edge)[eff$unit[eff$part == "C_logweight"]]))
write.csv(eff, file.path(out_dir, "presence_vs_proximity_effects.csv"), row.names = FALSE)

diag <- do.call(rbind, lapply(list(A_presence = fA, B_positive = fB, C_logweight = fC), function(f) {
  s <- posterior::summarise_draws(f, "rhat", "ess_bulk")
  np <- nuts_params(f)
  e <- s[grepl(",luad]", s$variable, fixed = TRUE) | s$variable == "b_luad", ]
  data.frame(max_rhat = max(s$rhat, na.rm = TRUE), min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
             max_rhat_effects = max(e$rhat), min_ess_effects = min(e$ess_bulk),
             divergences = sum(np$Value[np$Parameter == "divergent__"]),
             pop_effect = fixef(f)["luad", "Estimate"], pop_q025 = fixef(f)["luad", "Q2.5"],
             pop_q975 = fixef(f)["luad", "Q97.5"])
}))
diag$part <- rownames(diag)
write.csv(diag, file.path(out_dir, "presence_vs_proximity_diagnostics.csv"), row.names = FALSE)
print(diag, digits = 3)
cat("units whose 95% interval excludes zero:\n"); print(table(eff$part, eff$excludes_zero))
for (p in unique(eff$part)) {
  e <- eff[eff$part == p, ]
  cat("\n", p, "\n"); print(head(e[order(-abs(e$median)), c("unit", "median", "q025", "q975", "p_positive", "raw_d", "n")], 10), digits = 3, row.names = FALSE)
}
eC <- eff[eff$part == "C_logweight", ]
cat("\ncorrelation of the proximity effect (C) with the raw Cohen's d over 91 pairs:", round(stats::cor(eC$median, eC$raw_d), 2), "\n")
