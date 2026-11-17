test_that("halfspace depth is 2 * min(F, 1 - F) and lies in [0, 1]", {
  set.seed(1)
  ref <- rnorm(101)
  x <- c(min(ref) - 1, stats::median(ref), max(ref) + 1)
  d <- halfspace_depth(x, ref)
  expect_equal(d[1], 0)
  expect_equal(d[3], 0)
  expect_gt(d[2], 0.97)
  expect_true(all(halfspace_depth(ref, ref) >= 0 & halfspace_depth(ref, ref) <= 1))
})

test_that("depth is an affine function of the Fraiman-Muniz univariate depth", {
  set.seed(2)
  x <- rexp(60)
  fm1 <- 1 - abs(0.5 - stats::ecdf(x)(x))
  expect_equal(halfspace_depth(x, x), 2 * (fm1 - 0.5))
})

test_that("depth agrees with fda.usc::depth.FM on single-point curves", {
  skip_if_not_installed("fda.usc")
  set.seed(3)
  x <- c(rexp(40), rep(0, 20))
  fd <- fda.usc::fdata(matrix(x, ncol = 1), argvals = 0, rangeval = c(0, 0))
  expect_equal(unname(fda.usc::depth.FM(fd, fdataori = fd)$dep),
               halfspace_depth(x, x))
})

test_that("depth is folded: it does not distinguish the two tails", {
  ref <- 1:99
  expect_equal(halfspace_depth(10, ref), halfspace_depth(89, ref), tolerance = 0.03)
})

test_that("held-out rows do not enter the population reference", {
  set.seed(4)
  train <- matrix(rnorm(200), 50, 4)
  new1 <- matrix(rnorm(40), 10, 4)
  new2 <- rbind(new1[1:5, ], matrix(100, 5, 4))
  d1 <- population_depth(new1, train, B = 0)
  d2 <- population_depth(new2, train, B = 0)
  expect_equal(d1[1:5, ], d2[1:5, ])
})

test_that("within-sample depth uses only the graph itself and is rank based", {
  set.seed(5)
  x <- matrix(rexp(60), 3, 20)
  a <- within_depth(x, B = 0)
  b <- within_depth(rbind(x[1, , drop = FALSE], x[1, , drop = FALSE] * 0 + 7), B = 0)
  expect_equal(a[1, ], b[1, ])
  expect_equal(within_depth(log1p(x), B = 0), a)
  expect_equal(dim(within_rank(x)), dim(x))
})
