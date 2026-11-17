test_that("completed graphs are symmetric with unit diagonal and no negatives", {
  sim <- simulate_graphs(30, seed = 1)
  for (i in 1:30) {
    w <- complete_graph(sim$graphs[i, , ])
    expect_true(isSymmetric(unname(w)))
    expect_equal(unname(diag(w)), rep(1, 14))
    expect_true(all(w >= 0))
  }
})

test_that("an edge between two absent nodes stays at zero", {
  set.seed(2)
  a <- matrix(rexp(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
  a[c(2, 5), ] <- 0; a[, c(2, 5)] <- 0
  expect_equal(names(which(!node_presence(a))), c("b", "e"))
  w <- complete_graph(a)
  expect_equal(w["b", "e"], 0)
})

test_that("edge matrix has one named column per node pair", {
  sim <- simulate_graphs(3, seed = 1)
  x <- edge_matrix(sim$graphs)
  expect_equal(dim(x), c(3L, 91L))
  expect_true("JAKSTAT_TNFa" %in% colnames(x))
})

test_that("simulated graphs follow the documented format", {
  sim <- simulate_graphs(50, seed = 3)
  expect_equal(dim(sim$graphs), c(50L, 14L, 14L))
  expect_equal(rownames(sim$graphs), sim$slides$slide)
  expect_true(any(duplicated(sim$slides$patient)))
  expect_equal(sum(sim$slides$primary), length(unique(sim$slides$patient)))
  expect_equal(nrow(unique(sim$slides[, c("patient", "subtype")])),
               length(unique(sim$slides$patient)))
})

test_that("graphs written as CSV files are read back", {
  sim <- simulate_graphs(4, seed = 4)
  d <- file.path(tempdir(), "graphs_io"); dir.create(d, showWarnings = FALSE)
  for (i in 1:4) write.csv(sim$graphs[i, , ], file.path(d, paste0(sim$slides$slide[i], ".csv")))
  ann <- file.path(d, "names.txt")
  write.csv(data.frame(pat = sim$slides$slide, disease = tolower(sim$slides$subtype)), ann,
            row.names = FALSE)
  back <- read_graphs(d, ann)
  o <- order(sim$slides$slide)                 # read_graphs sorts by slide
  expect_equal(back$graphs, sim$graphs[o, , ])
  expect_equal(back$slides$subtype, sim$slides$subtype[o])
})
