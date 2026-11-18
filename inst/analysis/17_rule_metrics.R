# Rule metrics, pruning and held-out evaluation of the decision rules of the
# two depth forests, following the defaults of inTrees (Deng 2019): rules are
# the root-to-node paths of the first 100 trees, truncated at depth 6; each
# rule gets a length, a frequency (support) and an error (one minus
# confidence) recomputed on all training patients, with the prediction
# reassigned to the majority class of the patients it covers; conditions whose
# removal changes the error by little are pruned (maxDecay 0.05); rules with
# support of at least 5% and training error of at most 0.25 are kept and
# evaluated once on the held-out patients. For every condition
# the tail of the reference distribution that the covered patients occupy is
# reported, since a depth condition alone does not give a direction.
# Co-occurrence of features in the pruned rules is compared with its
# expectation under independence (lift), which removes the dependence of raw
# co-occurrence counts on how often each feature is used.
#
# Usage: Rscript 17_rule_metrics.R <output folder>
# Writes rule_metrics_all.csv, rule_metrics_selected.csv, rule_lift.csv,
#        rule_metrics_summary.csv

suppressMessages(library(pathwaydepth))
out_dir <- commandArgs(trailingOnly = TRUE)[1]
prim <- readRDS(file.path(out_dir, "primary.rds"))
stopifnot(identical(rownames(prim$feats$within$tr), prim$slides$slide[prim$train]),
          identical(rownames(prim$feats$within$te), prim$slides$slide[!prim$train]))
y_tr <- prim$y[prim$train]; y_te <- prim$y[!prim$train]
lev <- levels(y_tr)

N_TREES <- 100      # inTrees default (extractRules, ntree = 100)
MAX_DEPTH <- 6      # inTrees default (extractRules, maxdepth = 6)
MAX_DECAY <- 0.05   # inTrees default (pruneRule, maxDecay = 0.05)
MIN_SUPPORT <- 0.05 # smallest share of training patients covered by a kept rule
MAX_ERR <- 0.25     # largest training error of a kept rule (as in Fuh-Ngwa et al. 2022)
set.seed(17)

# membership of all patients in one rule given as (features, lower, upper)
covers <- function(x, f, lo, up) {
  ok <- rep(TRUE, nrow(x))
  for (k in seq_along(f)) ok <- ok & x[, f[k]] > lo[k] & x[, f[k]] <= up[k]
  ok
}
# inTrees metrics: prediction = majority class of covered patients,
# err = 1 - share of that class, freq = share of patients covered
metric <- function(m, y) {
  n <- sum(m)
  if (n == 0) return(list(n = 0, freq = 0, pred = NA_character_, err = 1))
  tab <- table(factor(y[m], levels = lev))
  list(n = n, freq = n / length(y), pred = names(which.max(tab)),
       err = 1 - max(tab) / n)
}
# inTrees pruneRule (typeDecay = 2): drop a condition when the relative
# increase in error, (err_new - err) / max(err, 1e-3), is at most MAX_DECAY;
# conditions are visited from the leaf end of the path towards the root
prune <- function(f, lo, up, x, y) {
  cur <- metric(covers(x, f, lo, up), y)
  k <- length(f)
  while (k >= 1 && length(f) > 1) {
    keep <- setdiff(seq_along(f), k)
    new <- metric(covers(x, f[keep], lo[keep], up[keep]), y)
    if ((new$err - cur$err) / max(cur$err, 1e-3) <= MAX_DECAY) {
      f <- f[keep]; lo <- lo[keep]; up <- up[keep]; cur <- new
    }
    k <- k - 1
  }
  list(f = f, lo = lo, up = up)
}
# root-to-node paths truncated at MAX_DEPTH (inTrees extractRules); each
# feature's conditions are merged into one interval lower < x <= upper
trunc_rules <- function(forest, trees) {
  out <- list()
  for (t in trees) {
    ti <- ranger::treeInfo(forest, t); rownames(ti) <- ti$nodeID
    stack <- list(list(node = 0, f = character(0), thr = numeric(0), left = logical(0)))
    while (length(stack)) {
      st <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
      row <- ti[as.character(st$node), ]
      if (row$terminal || length(st$f) == MAX_DEPTH) {
        if (!length(st$f)) next
        fe <- unique(st$f)
        lo <- vapply(fe, function(v) max(c(-Inf, st$thr[st$f == v & !st$left])), 0)
        up <- vapply(fe, function(v) min(c(Inf, st$thr[st$f == v & st$left])), 0)
        out[[length(out) + 1]] <- list(tree = t, f = fe, lo = unname(lo), up = unname(up), len_path = length(st$f))
        next
      }
      for (side in c(FALSE, TRUE)) stack[[length(stack) + 1]] <- list(
        node = if (side) row$leftChild else row$rightChild,
        f = c(st$f, row$splitvarName), thr = c(st$thr, row$splitval), left = c(st$left, side))
    }
  }
  out
}
fmt_cond <- function(f, lo, up) paste(vapply(seq_along(f), function(k) {
  if (is.finite(lo[k]) && is.finite(up[k])) sprintf("%.3g < %s <= %.3g", lo[k], f[k], up[k])
  else if (is.finite(up[k])) sprintf("%s <= %.3g", f[k], up[k])
  else sprintf("%s > %.3g", f[k], lo[k])
}, ""), collapse = " & ")

all_rules <- list(); sel_rules <- list(); lift_tab <- list(); summ <- list()
for (k in c("within", "popdepth")) {
  fit <- prim$fits[[k]]
  xtr <- prim$feats[[k]]$tr; xte <- prim$feats[[k]]$te
  # direction: signed position of the edge weight in the same reference
  # (within-slide rank for Patient Referenced Depth; training-cohort
  # percentile of the edge for Population Referenced Depth)
  if (k == "within") {
    dir_tr <- prim$feats$rank$tr
  } else {
    raw <- prim$feats$raw$tr
    dir_tr <- sapply(colnames(raw), function(j) stats::ecdf(raw[, j])(raw[, j]))
  }
  r <- trunc_rules(fit$forest, seq_len(N_TREES))
  res <- vector("list", length(r))
  for (i in seq_along(r)) {
    f <- r[[i]]$f; lo <- r[[i]]$lo; up <- r[[i]]$up
    m0 <- metric(covers(xtr, f, lo, up), y_tr)
    p <- prune(f, lo, up, xtr, y_tr)
    mtr <- covers(xtr, p$f, p$lo, p$up); mte <- covers(xte, p$f, p$lo, p$up)
    a <- metric(mtr, y_tr); b <- metric(mte, y_te)
    # tail occupied by the covered training patients, per condition:
    # share of covered patients whose signed position is above one half
    upper_share <- vapply(p$f, function(v) mean(dir_tr[mtr, v] > 0.5), 0)
    res[[i]] <- data.frame(
      model = k, tree = r[[i]]$tree,
      len_raw = length(f), n_raw = m0$n, err_raw = m0$err,
      len = length(p$f), condition = fmt_cond(p$f, p$lo, p$up),
      features = paste(p$f, collapse = ";"),
      upper_tail_share = paste(sprintf("%.2f", upper_share), collapse = ";"),
      pred = a$pred, n_train = a$n, freq_train = a$freq, err_train = a$err,
      n_test = b$n, freq_test = b$freq,
      acc_test = if (b$n > 0) mean(y_te[mte] == a$pred) else NA_real_,
      n_test_correct = if (b$n > 0) sum(y_te[mte] == a$pred) else 0L)
  }
  d <- do.call(rbind, res)
  all_rules[[k]] <- d
  # distinct pruned rules with adequate support
  s <- d[!duplicated(d$condition) & d$freq_train >= MIN_SUPPORT & d$err_train <= MAX_ERR, ]
  base <- mean(y_te == "LUAD")
  s$base_test <- ifelse(s$pred == "LUAD", base, 1 - base)
  s$p_test <- mapply(function(x, n, p0) if (n > 0) stats::binom.test(x, n, p0, alternative = "greater")$p.value else NA_real_,
                     s$n_test_correct, s$n_test, s$base_test)
  s$fdr_test <- stats::p.adjust(s$p_test, "BH")
  s <- s[order(s$err_train, -s$freq_train), ]
  sel_rules[[k]] <- s
  # lift of feature co-occurrence in the distinct supported rules:
  # observed number of rules containing both features over the number
  # expected if features entered rules independently
  feats <- strsplit(s$features, ";")
  vars <- sort(unique(unlist(feats))); R <- length(feats)
  inc <- t(vapply(feats, function(v) vars %in% v, logical(length(vars))))
  colnames(inc) <- vars
  co <- crossprod(inc); fr <- colSums(inc) / R
  pr <- which(upper.tri(co), arr.ind = TRUE)
  lt <- data.frame(model = k, a = vars[pr[, 1]], b = vars[pr[, 2]],
                   n_both = co[pr], expected = R * fr[pr[, 1]] * fr[pr[, 2]])
  lt$lift <- lt$n_both / lt$expected
  # exact one-sided test of the overlap (hypergeometric, conditional on margins)
  lt$p <- stats::phyper(lt$n_both - 1, colSums(inc)[pr[, 2]], R - colSums(inc)[pr[, 2]],
                        colSums(inc)[pr[, 1]], lower.tail = FALSE)
  lt$fdr <- stats::p.adjust(lt$p, "BH")
  lift_tab[[k]] <- lt[order(-lt$n_both), ]
  # does lift remove the frequency dependence?
  tot <- rowSums(co) - diag(co); use <- diag(co)
  summ[[k]] <- data.frame(
    model = k, n_trees = N_TREES, max_depth = MAX_DEPTH, min_support = MIN_SUPPORT, max_err = MAX_ERR,
    n_rules = nrow(d), median_len_raw = stats::median(d$len_raw),
    median_n_raw = stats::median(d$n_raw), median_len_pruned = stats::median(d$len),
    n_distinct_pruned = sum(!duplicated(d$condition)), n_selected = nrow(s),
    median_freq_selected = stats::median(s$freq_train),
    median_err_selected = stats::median(s$err_train),
    median_acc_test_selected = stats::median(s$acc_test, na.rm = TRUE),
    n_selected_fdr05 = sum(s$fdr_test < 0.05, na.rm = TRUE),
    rho_cooc_use = stats::cor(tot, use, method = "spearman"),
    rho_lift_use = {
      ml <- tapply(c(lt$lift, lt$lift), c(lt$a, lt$b), mean)
      stats::cor(ml[names(use)], use, method = "spearman", use = "complete.obs") },
    n_pairs_lift_fdr05 = sum(lt$fdr < 0.05))
}
write.csv(do.call(rbind, all_rules), file.path(out_dir, "rule_metrics_all.csv"), row.names = FALSE)
write.csv(do.call(rbind, sel_rules), file.path(out_dir, "rule_metrics_selected.csv"), row.names = FALSE)
write.csv(do.call(rbind, lift_tab), file.path(out_dir, "rule_lift.csv"), row.names = FALSE)
write.csv(do.call(rbind, summ), file.path(out_dir, "rule_metrics_summary.csv"), row.names = FALSE)
print(do.call(rbind, summ))
for (k in names(sel_rules)) {
  cat("\n==", k, "top selected rules\n")
  print(utils::head(sel_rules[[k]][, c("condition", "upper_tail_share", "pred", "n_train", "err_train",
                                       "n_test", "acc_test", "fdr_test")], 12), row.names = FALSE)
  cat("\n==", k, "top co-occurring pairs by lift (n_both >= 5)\n")
  lt <- lift_tab[[k]]; lt <- lt[lt$n_both >= 5, ]
  print(utils::head(lt[order(-lt$lift), ], 8), row.names = FALSE)
}
