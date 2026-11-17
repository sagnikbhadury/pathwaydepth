#' Scaled univariate halfspace depth
#'
#' For a value `x` and a reference sample with empirical distribution function
#' `F`, the depth is `2 * min(F(x), 1 - F(x))`. It equals one at the median of
#' the reference sample, decreases towards both tails and carries no sign. It
#' is the Fraiman and Muniz univariate depth `1 - |1/2 - F(x)|` rescaled from
#' `[1/2, 1]` to the unit interval, `2 * (1 - |1/2 - F(x)| - 1/2)`, and it is
#' the quantity returned by `fda.usc::depth.FM()` with its defaults
#' (`dfunc = "FM1"`, `par.dfunc = list(scale = TRUE)`) on curves observed at a
#' single argument value. The `trim` argument of `fda.usc::depth.FM()` selects
#' the curves entering a trimmed mean and does not change the depth.
#'
#' @param x Numeric vector of values to score.
#' @param ref Numeric vector, the reference sample.
#' @return A numeric vector of depths in the unit interval, the same length as
#'   `x`.
#' @examples
#' halfspace_depth(c(-3, 0, 3), ref = rnorm(200))
#' @export
halfspace_depth <- function(x, ref) {
  ref <- sort(ref[is.finite(ref)])
  if (!length(ref)) return(rep(NA_real_, length(x)))
  f <- findInterval(x, ref) / length(ref)
  2 * pmin(f, 1 - f)
}

#' Subsampled halfspace depth
#'
#' Draws `B` subsamples of the reference sample without replacement, computes
#' [halfspace_depth()] of every value of `x` against each subsample, and
#' returns the median over subsamples.
#'
#' @param x Numeric vector of values to score.
#' @param ref Numeric vector, the reference sample.
#' @param frac Fraction of the reference sample in each subsample.
#' @param B Number of subsamples. `B = 0` uses the full reference sample once.
#' @return A numeric vector the same length as `x`.
#' @examples
#' set.seed(1)
#' subsample_depth(c(-1, 0, 1), ref = rnorm(100), B = 50)
#' @export
subsample_depth <- function(x, ref, frac = 0.8, B = 200) {
  if (B < 1) return(halfspace_depth(x, ref))
  n <- length(ref)
  m <- max(2L, floor(frac * n))
  d <- matrix(NA_real_, length(x), B)
  for (b in seq_len(B)) {
    d[, b] <- halfspace_depth(x, ref[sample.int(n, m)])
  }
  apply(d, 1, stats::median)
}

#' Population-referenced depth features
#'
#' For each edge (column), scores every row of `new` against the distribution
#' of that edge over the rows of `reference`. To score a held-out set, pass the
#' training edge matrix as `reference`; the held-out rows then play no part in
#' the reference distribution.
#'
#' @param new Numeric matrix, graphs by edges, to be scored.
#' @param reference Numeric matrix with the same columns as `new`.
#' @inheritParams subsample_depth
#' @return A numeric matrix with the dimensions and dimnames of `new`.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs[1:40, , ])
#' d <- population_depth(x[31:40, ], reference = x[1:30, ], B = 20)
#' range(d)
#' @export
population_depth <- function(new, reference = new, frac = 0.8, B = 200) {
  stopifnot(ncol(new) == ncol(reference))
  out <- new
  for (j in seq_len(ncol(new))) {
    out[, j] <- subsample_depth(new[, j], reference[, j], frac, B)
  }
  out
}

#' Within-sample depth features
#'
#' For each graph (row), scores every edge weight against the distribution of
#' all edge weights of the same graph. No other graph is used, so a held-out
#' graph is scored exactly as a training graph is.
#'
#' The feature is a rank within the graph. It compares weights of different
#' node pairs on one scale, which is meaningful when all edges are measured in
#' the same unit, as they are for a proximity statistic computed with one
#' radius.
#'
#' @param x Numeric matrix, graphs by edges.
#' @inheritParams subsample_depth
#' @return A numeric matrix with the dimensions and dimnames of `x`.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs[1:5, , ])
#' d <- within_depth(x, B = 20)
#' range(d)
#' @export
within_depth <- function(x, frac = 0.8, B = 200) {
  out <- x
  for (i in seq_len(nrow(x))) {
    out[i, ] <- subsample_depth(x[i, ], x[i, ], frac, B)
  }
  out
}

#' Within-sample rank features (unfolded)
#'
#' The empirical distribution function of each edge weight among the edge
#' weights of the same graph. Unlike [within_depth()], the value keeps the
#' direction: it is near one for the strongest edges of the graph and near zero
#' for the weakest.
#'
#' @param x Numeric matrix, graphs by edges.
#' @return A numeric matrix with the dimensions and dimnames of `x`.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' r <- within_rank(edge_matrix(sim$graphs[1:5, , ]))
#' range(r)
#' @export
within_rank <- function(x) {
  out <- x
  for (i in seq_len(nrow(x))) {
    out[i, ] <- findInterval(x[i, ], sort(x[i, ])) / ncol(x)
  }
  out
}
