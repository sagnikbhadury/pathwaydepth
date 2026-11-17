test_that("extracted rules reproduce every tree exactly", {
  sim <- simulate_graphs(200, seed = 1)
  x <- edge_matrix(sim$graphs[1:150, , ])
  y <- sim$slides$subtype[1:150]
  set.seed(1)
  f <- ranger::ranger(x = x, y = y, num.trees = 5, seed = 1, num.threads = 1)
  rules <- extract_rules(f, trees = 1:5)
  per_tree <- stats::predict(f, data = x, predict.all = TRUE)$predictions
  for (k in 1:5) {
    rk <- rules[rules$tree == k, ]
    ti <- ranger::treeInfo(f, k)
    expect_equal(nrow(rk), sum(ti$terminal))            # one rule per leaf
    m <- rule_matrix(rk, x)
    expect_true(all(rowSums(m) == 1))                   # a partition of the cases
    pred <- rk$prediction[apply(m, 1, which)]
    expect_equal(pred, levels(y)[per_tree[, k]])        # same prediction as the tree
  }
})

test_that("rules work for the probability forest of fit_forest", {
  sim <- simulate_graphs(200, seed = 1)
  x <- edge_matrix(sim$graphs[1:100, , ])
  fit <- fit_forest(x, sim$slides$subtype[1:100], num_trees = 10,
                    tune_trees = 10, num_threads = 1)
  r <- extract_rules(fit, trees = 1:3)
  expect_true(all(r$prediction %in% c("LUSC", "LUAD")))
  expect_true(all(rule_membership(r[r$tree == 2, ], x) == 1))
})
