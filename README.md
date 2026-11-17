# pathwaydepth

Depth features for small weighted graphs, and the analysis of spatial pathway
graphs of lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC)
reported in

> Bhadury S, Tsang AP, Ray D, Lagisetty K, Rao A. Depth features of spatial
> pathway graphs from histology-predicted transcriptomics in lung
> adenocarcinoma and squamous cell carcinoma. Revised manuscript, 2026.

## What the package does

Each slide is summarised by a graph on 14 PROGENy pathways. The 91 edge
weights are turned into edge-level depth features under two reference frames:

| Function | Reference distribution for edge *j* of slide *i* |
|---|---|
| `population_depth()` | edge *j* across the training slides |
| `within_depth()` | the 91 edges of slide *i* |

The depth is `2 * min(F(x), 1 - F(x))`, the Fraiman and Muniz univariate depth
rescaled to the unit interval. It is one at the median of the reference
distribution and falls towards both tails, so it is unsigned: a high and a low
edge weight can have the same depth. `within_rank()` gives the signed version
of the within-sample feature.

Also included: graph completion and symmetrisation (`complete_graph()`),
patient-level splitting (`grouped_split()`), a random forest tuned on
out-of-bag predictions (`fit_forest()`), test-set intervals (`evaluate()`),
DeLong's paired comparison (`compare_auc()`) and feature-wise group
comparisons (`compare_groups()`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/pathwaydepth")
```

## Example

```r
library(pathwaydepth)
sim <- simulate_graphs(300, seed = 1)          # artificial graphs, for illustration

keep   <- sim$slides$primary                   # one slide per patient
x      <- edge_matrix(sim$graphs[keep, , ])    # edge weights, one row per slide
y      <- sim$slides$subtype[keep]
train  <- grouped_split(sim$slides$patient[keep], y)

d_tr <- within_depth(x[train, ])
d_te <- within_depth(x[!train, ])
fit  <- fit_forest(d_tr, y[train])
evaluate(predict_forest(fit, d_te), y[!train], fit$threshold)
```

## Data

This repository contains code only. The spatial pathway graphs analysed in the
article (one 14 x 14 directed G-cross proximity matrix for each of 996
diagnostic slides of TCGA-LUAD and TCGA-LUSC) are **not distributed**. They
are derived from the DeepSpot predicted spatial transcriptomics of TCGA slides
(Nonchev et al., 2025, doi:10.1101/2025.02.09.25321567), which is released
under CC BY-NC-ND 4.0, a licence that does not allow derivatives to be
distributed. The DeepSpot data can be requested from its authors at
<https://huggingface.co/datasets/nonchev/TCGA_digital_spatial_transcriptomics>.

`read_graphs()` reads such matrices from a folder with one CSV file per slide
and an annotation file with the columns `pat` (slide barcode) and `disease`
(`luad` or `lusc`). `simulate_graphs()` produces artificial graphs in the same
format.

## Reproducing the manuscript

The scripts are in `inst/analysis`. They need the graph files (see Data) and
write to an output folder of your choice.

```sh
Rscript inst/analysis/01_primary.R          results <folder of graph CSV files> <annotation CSV>
Rscript inst/analysis/02_repeated_splits.R  results primary 50
Rscript inst/analysis/02_repeated_splits.R  results all 50
Rscript inst/analysis/03_bulk_baselines.R   results <folder with GDC count files>
Rscript inst/analysis/04_rules_and_checks.R results
Rscript inst/analysis/05_figures.R          results figures
Rscript inst/analysis/06_original_figures.R results <figure script> figures_original
```

Script 03 needs the STAR gene counts of TCGA-LUAD and TCGA-LUSC from the
Genomic Data Commons; they are not redistributed here.

## Licence

MIT.
