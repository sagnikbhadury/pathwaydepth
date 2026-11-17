test_that("no patient is on both sides of a grouped split", {
  sim <- simulate_graphs(200, seed = 1)
  set.seed(1)
  tr <- grouped_split(sim$slides$patient, sim$slides$subtype)
  expect_false(any(sim$slides$patient[tr] %in% sim$slides$patient[!tr]))
  expect_equal(mean(tr), 0.7, tolerance = 0.03)
})

test_that("forest fit returns out-of-bag threshold and evaluates", {
  set.seed(2)
  y <- factor(rep(c("a", "b"), each = 60))
  x <- matrix(rnorm(120 * 6), 120, 6, dimnames = list(NULL, paste0("f", 1:6)))
  x[y == "b", 1] <- x[y == "b", 1] + 1.5
  fit <- fit_forest(x, y, num_trees = 100, tune_trees = 50, num_threads = 1)
  expect_gt(fit$oob_auc, 0.7)
  expect_true(fit$threshold > 0 && fit$threshold < 1)
  expect_length(predict_forest(fit, x), 120)
  e <- evaluate(fit$oob_prob, y, fit$threshold)
  expect_equal(e$metric, c("AUC", "Accuracy", "Sensitivity", "Specificity"))
  cmp <- compare_auc(x[, 1], x[, 2], y)
  expect_gt(cmp$difference, 0)
})

test_that("group comparison reports Cohen's d as second level minus first", {
  y <- factor(rep(c("a", "b"), each = 30))
  set.seed(3)
  x <- cbind(up = rnorm(60) + rep(c(0, 1), each = 30), flat = rnorm(60))
  r <- compare_groups(x, y)
  expect_equal(r$feature[1], "up")
  expect_gt(r$d[1], 0)
})
