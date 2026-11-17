# Small tables for the overview figure (Figure 1), from the saved primary analysis.
# Usage: Rscript 09_overview_inputs.R <output folder of 01_primary.R>
# Writes <output folder>/pipeline_inputs/.
suppressMessages(library(pathwaydepth))
out_dir <- commandArgs(trailingOnly = TRUE)[1]
dir.create(file.path(out_dir, "pipeline_inputs"), showWarnings = FALSE)
prim <- readRDS(file.path(out_dir, "primary.rds"))
dat <- readRDS(file.path(out_dir, "graphs.rds"))
x <- prim$feats$raw$tr; y <- prim$y[prim$train]
# example patient: all 14 pathways present, first such training patient
present <- prim$present[prim$train, ]
i <- which(rowSums(present) == 14)[1]
slide <- prim$slides$slide[prim$train][i]
w <- complete_graph(dat$graphs[slide, , ])
write.csv(w, file.path(out_dir, "pipeline_inputs", "example_graph.csv"))
write.csv(data.frame(edge = colnames(x), weight = x[i, ],
                     pop = prim$feats$popdepth$tr[i, ], within = prim$feats$within$tr[i, ]),
          file.path(out_dir, "pipeline_inputs", "example_edges.csv"), row.names = FALSE)
write.csv(data.frame(JAKSTAT_TNFa = x[, "JAKSTAT_TNFa"], subtype = y),
          file.path(out_dir, "pipeline_inputs", "edge_across_patients.csv"), row.names = FALSE)
set.seed(1); s <- c(sample(which(y == "LUSC"), 30), sample(which(y == "LUAD"), 30))
write.csv(prim$feats$within$tr[s, ], file.path(out_dir, "pipeline_inputs", "within_matrix.csv"), row.names = FALSE)
write.csv(prim$feats$popdepth$tr[s, ], file.path(out_dir, "pipeline_inputs", "pop_matrix.csv"), row.names = FALSE)
cat("example:", slide, as.character(y[i]), "\n")
