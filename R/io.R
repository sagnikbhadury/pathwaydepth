#' Read directed proximity matrices from a folder of CSV files
#'
#' Reads one square matrix per slide. Each file is named `<slide>.csv`, has the
#' node names in its first row and first column, and holds the directed
#' proximity of every ordered pair of nodes. The annotation file has a column
#' `pat` with the slide identifiers and a column `disease` with the class
#' label.
#'
#' The package does not include the lung cancer graphs analysed in the
#' accompanying article. They are derived from the DeepSpot predicted spatial
#' transcriptomics of TCGA slides, whose licence (CC BY-NC-ND 4.0) does not
#' allow derivatives to be distributed. Users who have obtained the DeepSpot
#' data can build the matrices and read them with this function.
#'
#' @param graph_dir Folder that contains the CSV files.
#' @param annotation Path of the annotation CSV file.
#' @param levels Class labels, in the order of the factor levels; the second
#'   level is the positive class.
#' @return A list with `graphs`, a numeric array of dimension slides by nodes
#'   by nodes, and `slides`, a data frame with the columns `slide`, `patient`
#'   (first twelve characters of the slide identifier), `subtype` and
#'   `primary` (`TRUE` for the first slide of each patient).
#' @export
read_graphs <- function(graph_dir, annotation, levels = c("LUSC", "LUAD")) {
  info <- utils::read.csv(annotation, stringsAsFactors = FALSE)
  info <- info[order(info$pat), ]
  files <- file.path(graph_dir, paste0(info$pat, ".csv"))
  stopifnot(all(file.exists(files)))
  first <- as.matrix(utils::read.csv(files[1], check.names = FALSE, row.names = 1))
  nodes <- colnames(first)
  graphs <- array(NA_real_, c(length(files), length(nodes), length(nodes)),
                  dimnames = list(info$pat, nodes, nodes))
  for (i in seq_along(files)) {
    a <- as.matrix(utils::read.csv(files[i], check.names = FALSE, row.names = 1))
    stopifnot(identical(colnames(a), nodes), identical(rownames(a), nodes))
    graphs[i, , ] <- a
  }
  slides <- data.frame(slide = info$pat, patient = patient_id(info$pat),
                       subtype = factor(toupper(info$disease), levels = levels),
                       stringsAsFactors = FALSE)
  slides$primary <- !duplicated(slides$patient)
  list(graphs = graphs, slides = slides)
}

#' Simulated proximity graphs
#'
#' Generates artificial directed proximity matrices for two classes, for use
#' in examples and tests. Nodes are absent from a slide at random (all-zero
#' row and column), weights are drawn from gamma distributions, and a few
#' edges differ in mean between the classes. Some patients contribute two
#' slides. The values are random numbers and are not derived from any data.
#'
#' @param n Number of slides.
#' @param nodes Node names.
#' @param effect Multiplicative shift of the mean weight of the edges that
#'   differ between classes.
#' @param seed Optional seed.
#' @return A list with `graphs` and `slides`, in the format of
#'   [read_graphs()].
#' @examples
#' sim <- simulate_graphs(20, seed = 1)
#' dim(sim$graphs)
#' head(sim$slides)
#' @export
simulate_graphs <- function(n = 200,
                            nodes = c("Androgen", "EGFR", "Estrogen", "Hypoxia",
                                      "JAK-STAT", "MAPK", "NFkB", "PI3K", "TGFb",
                                      "TNFa", "Trail", "VEGF", "WNT", "p53"),
                            effect = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- length(nodes)
  n_pat <- ceiling(0.9 * n)
  pat <- c(seq_len(n_pat), seq_len(n - n_pat))          # the first patients have two slides
  subtype <- factor(c("LUSC", "LUAD")[pat %% 2 + 1], levels = c("LUSC", "LUAD"))  # one class per patient
  patient <- sprintf("SIM-%02d-%04d", as.integer(subtype), pat)
  slide <- paste0(patient, "-01Z-00-DX", stats::ave(pat, pat, FUN = seq_along))
  absent_prob <- seq(0, 0.8, length.out = p)
  base <- matrix(stats::rgamma(p * p, shape = 2, rate = 0.2), p, p)
  differ <- matrix(FALSE, p, p); differ[1:4, 1:4] <- TRUE
  graphs <- array(0, c(n, p, p), dimnames = list(slide, nodes, nodes))
  for (i in seq_len(n)) {
    mu <- base * ifelse(differ & subtype[i] == "LUAD", effect, 1)
    a <- matrix(stats::rgamma(p * p, shape = 2, rate = 2 / mu), p, p)
    a[stats::runif(p * p) < 0.2] <- 0
    gone <- stats::runif(p) < absent_prob
    a[gone, ] <- 0; a[, gone] <- 0
    graphs[i, , ] <- a
  }
  slides <- data.frame(slide = slide, patient = patient, subtype = subtype,
                       stringsAsFactors = FALSE)
  slides$primary <- !duplicated(slides$patient)
  rownames(slides) <- NULL
  list(graphs = graphs, slides = slides)
}
