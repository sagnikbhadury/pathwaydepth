#' Complete and symmetrise a directed spatial proximity matrix
#'
#' A node that is absent from a slide has an all-zero row and an all-zero
#' column in the directed proximity matrix. Rows and columns that are entirely
#' zero are set to missing and the matrix is completed with Soft-Impute. The
#' diagonal is then set to one, negative values are set to zero, and the matrix
#' is symmetrised: where both directions are positive the edge weight is their
#' geometric mean, otherwise it is the larger of the two.
#'
#' Soft-Impute has no information on a row or a column with no observed entry,
#' so the completed value for an edge between two absent nodes is zero. An edge
#' between an absent node and a present node can receive a small positive
#' value only when the absent node has an observed row or column.
#' [node_presence()] reports which nodes were observed.
#'
#' @param adj A square numeric matrix of non-negative directed weights.
#' @param lambda Soft-Impute nuclear norm penalty.
#' @param rank_max Maximum rank of the Soft-Impute fit. Defaults to one less
#'   than the number of nodes.
#' @return A symmetric matrix with unit diagonal and the dimnames of `adj`.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' w <- complete_graph(sim$graphs[1, , ])
#' isSymmetric(w)
#' @export
complete_graph <- function(adj, lambda = 0.1, rank_max = nrow(adj) - 1) {
  adj <- as.matrix(adj)
  stopifnot(is.numeric(adj), nrow(adj) == ncol(adj))
  zero_rows <- apply(adj == 0, 1, all)
  zero_cols <- apply(adj == 0, 2, all)
  filled <- adj
  if (any(zero_rows) || any(zero_cols)) {
    with_na <- adj
    with_na[, zero_cols] <- NA
    with_na[zero_rows, ] <- NA
    if (all(is.na(with_na))) {
      completed <- adj
    } else {
      fit <- softImpute::softImpute(with_na, rank.max = rank_max,
                                    lambda = lambda, type = "svd")
      completed <- softImpute::complete(with_na, fit)
    }
    filled[adj == 0] <- completed[adj == 0]
  }
  diag(filled) <- 1
  filled[filled < 0] <- 0
  both <- (filled > 0) & (t(filled) > 0)
  out <- ifelse(both, sqrt(filled * t(filled)), pmax(filled, t(filled)))
  dimnames(out) <- dimnames(adj)
  out
}

#' Which nodes are observed in a directed proximity matrix
#'
#' @param adj A square numeric matrix.
#' @return A named logical vector, `TRUE` when the node has at least one
#'   non-zero entry in its row or its column.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' node_presence(sim$graphs[1, , ])
#' @export
node_presence <- function(adj) {
  adj <- as.matrix(adj)
  present <- !(apply(adj == 0, 1, all) & apply(adj == 0, 2, all))
  names(present) <- colnames(adj)
  present
}

#' Names of the upper-triangular edges of a graph
#'
#' @param nodes Character vector of node names.
#' @param sep Separator placed between the two node names.
#' @return A character vector of length `choose(length(nodes), 2)`, in the
#'   column-major order used by [upper.tri()].
#' @examples
#' edge_names(c("A", "B", "C"))
#' @export
edge_names <- function(nodes, sep = "_") {
  m <- matrix(FALSE, length(nodes), length(nodes))
  idx <- which(upper.tri(m), arr.ind = TRUE)
  paste0(nodes[idx[, 1]], sep, nodes[idx[, 2]])
}

#' Edge weight matrix for a set of graphs
#'
#' Applies [complete_graph()] to every slice of an array of directed proximity
#' matrices and returns the upper-triangular edge weights, one row per graph.
#'
#' @param graphs A numeric array with dimensions graphs by nodes by nodes.
#' @param complete If `FALSE`, the matrices are only symmetrised by the rule of
#'   [complete_graph()] and no completion is attempted.
#' @param ... Passed to [complete_graph()].
#' @return A numeric matrix with one row per graph and one column per edge.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' x <- edge_matrix(sim$graphs[1:5, , ])
#' dim(x)
#' @export
edge_matrix <- function(graphs, complete = TRUE, ...) {
  stopifnot(length(dim(graphs)) == 3, dim(graphs)[2] == dim(graphs)[3])
  n <- dim(graphs)[1]
  nodes <- dimnames(graphs)[[2]]
  if (is.null(nodes)) nodes <- paste0("V", seq_len(dim(graphs)[2]))
  nodes <- gsub("[^A-Za-z0-9]", "", nodes)
  ut <- upper.tri(matrix(0, length(nodes), length(nodes)))
  out <- matrix(NA_real_, n, sum(ut),
                dimnames = list(dimnames(graphs)[[1]], edge_names(nodes)))
  for (i in seq_len(n)) {
    a <- graphs[i, , ]
    if (complete) {
      w <- complete_graph(a, ...)
    } else {
      diag(a) <- 1
      both <- (a > 0) & (t(a) > 0)
      w <- ifelse(both, sqrt(a * t(a)), pmax(a, t(a)))
    }
    out[i, ] <- w[ut]
  }
  out
}

#' Node presence matrix for a set of graphs
#'
#' @param graphs A numeric array with dimensions graphs by nodes by nodes.
#' @return A logical matrix with one row per graph and one column per node.
#' @examples
#' sim <- simulate_graphs(200, seed = 1)
#' colMeans(presence_matrix(sim$graphs))
#' @export
presence_matrix <- function(graphs) {
  out <- t(apply(graphs, 1, node_presence))
  dimnames(out) <- list(dimnames(graphs)[[1]], dimnames(graphs)[[2]])
  out
}
