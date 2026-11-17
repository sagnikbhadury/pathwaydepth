"""Workflow figure of the manuscript, assembled from the analysis itself:
the slide-to-graph panel, an example graph and its two depth references
(tables written by export_pipeline_inputs.R), and panels of the result figures.

Usage (run in the folder that holds pipeline_inputs/ written by 09_overview_inputs.R,
TCGA_ST_graph_fig.pdf and the result figures of scripts 05 and 06):
    python 10_overview_figure.py <output pdf>
"""
import sys
import numpy as np
import pandas as pd
import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
ORANGE, AQUA = "#eb6834", "#1baf7a"

def crop(pdf, box, dpi=240):
    """Render a region of a figure PDF; box = (x0, y0, x1, y1) as fractions of the page."""
    page = fitz.open(pdf)[0]
    r = page.rect
    clip = fitz.Rect(r.x0 + box[0] * r.width, r.y0 + box[1] * r.height,
                     r.x0 + box[2] * r.width, r.y0 + box[3] * r.height)
    pix = page.get_pixmap(dpi=dpi, clip=clip)
    return np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[..., :3]

fig = plt.figure(figsize=(11.5, 9.4))

def image(rect, img, anchor="N"):
    ax = fig.add_axes(rect); ax.imshow(img); ax.set_axis_off(); ax.set_anchor(anchor)
    return ax

def step(x, y, n, text):
    fig.text(x, y, str(n), ha="center", va="center", fontsize=8.5, fontweight="bold", color="white",
             bbox=dict(boxstyle="round,pad=0.25", fc=INK, ec="none"))
    fig.text(x + 0.016, y, text, ha="left", va="center", fontsize=8, color=INK, linespacing=1.15)

def title(y, text):
    fig.text(0.02, y, text, ha="left", va="center", fontsize=10.5, fontweight="bold", color=INK)

def clean(ax):
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=6.5, length=2)

def arrow(y0, y1):
    fig.add_artist(FancyArrowPatch((0.5, y0), (0.5, y1), transform=fig.transFigure,
                                   arrowstyle="-|>", mutation_scale=18, color=MUTED, lw=1.6))

g = pd.read_csv("pipeline_inputs/example_graph.csv", index_col=0)
e = pd.read_csv("pipeline_inputs/example_edges.csv")
pop = pd.read_csv("pipeline_inputs/edge_across_patients.csv")
wm = pd.read_csv("pipeline_inputs/within_matrix.csv").to_numpy()
ex = e[e.edge == "JAKSTAT_TNFa"].iloc[0]

# ---- A: slide -> predicted expression -> graph ------------------------------
title(0.975, "A   From slide to graph")
image([0.03, 0.725, 0.94, 0.225], crop("TCGA_ST_graph_fig.pdf", (0, 0, 1, 1)))
step(0.075, 0.705, 1, "H&E slide (TCGA)")
step(0.345, 0.705, 2, "DeepSpot-predicted expression; each spot labelled\nby its highest-scoring PROGENy pathway")
step(0.735, 0.705, 3, "G-cross proximity between\npathway-labelled spots")
arrow(0.682, 0.655)

# ---- B: graph -> edge weights -> two depths -> feature matrix ----------------
title(0.635, "B   From graph to depth features")
y0, h = 0.405, 0.165
ax = fig.add_axes([0.075, y0, 0.135, h])
m = g.to_numpy().copy(); np.fill_diagonal(m, np.nan)
im = ax.imshow(m, cmap="Blues", vmin=0)
ax.set_xticks(range(14)); ax.set_yticks(range(14))
ax.set_xticklabels(g.columns, rotation=90); ax.set_yticklabels(g.index)
clean(ax); ax.tick_params(labelsize=5.2)
cax = fig.add_axes([0.214, y0, 0.006, h]); cb = fig.colorbar(im, cax=cax)
cb.ax.tick_params(labelsize=6, colors=MUTED); cb.outline.set_edgecolor(GRID)
step(0.045, 0.595, 4, "Edge weights of one slide\n(14 pathways, 91 pairs)")

ax = fig.add_axes([0.305, y0, 0.175, h])
ax.hist(pop.JAKSTAT_TNFa, bins=30, color=ORANGE, alpha=0.9)
ax.axvline(ex.weight, color=INK, lw=1.4)
ax.text(ex.weight + 1, ax.get_ylim()[1] * 0.96, f"this slide\ndepth {ex['pop']:.2f}", fontsize=7, va="top", color=INK)
ax.set_xlabel("JAK-STAT ↔ TNFα weight in the training patients", fontsize=7, color=INK)
ax.set_ylabel("Patients", fontsize=7, color=INK); clean(ax)
step(0.285, 0.595, "5a", "Population Referenced Depth:\nthe same edge across patients")

ax = fig.add_axes([0.545, y0, 0.175, h])
ax.hist(e.weight, bins=30, color=AQUA, alpha=0.95)
ax.axvline(ex.weight, color=INK, lw=1.4)
ax.text(ex.weight + 1, ax.get_ylim()[1] * 0.96, f"JAK-STAT ↔ TNFα\ndepth {ex['within']:.2f}", fontsize=7,
        va="top", color=INK)
ax.set_xlabel("Weights of the 91 edges of this slide", fontsize=7, color=INK)
ax.set_ylabel("Edges", fontsize=7, color=INK); clean(ax)
step(0.525, 0.595, "5b", "Patient Referenced Depth:\nall edges within the slide")

ax = fig.add_axes([0.80, y0, 0.15, h])
im = ax.imshow(wm, aspect="auto", cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
ax.axhline(29.5, color="white", lw=1.5)
ax.set_yticks([15, 45]); ax.set_yticklabels(["LUSC", "LUAD"], rotation=90, va="center")
ax.set_xticks([]); ax.set_xlabel("91 pathway pairs", fontsize=7, color=INK); clean(ax)
cax = fig.add_axes([0.955, y0, 0.006, h]); cb = fig.colorbar(im, cax=cax)
cb.ax.tick_params(labelsize=6, colors=MUTED); cb.outline.set_edgecolor(GRID)
step(0.775, 0.595, 6, "Depth feature matrix\n(60 training patients shown)")
arrow(0.362, 0.335)

# ---- C: classification and interpretation ------------------------------------
title(0.315, "C   Classification, evaluation and interpretation")
panels = [
    ([0.02, 0.02, 0.185, 0.245], "Figure1_Performance.pdf", (0.69, 0.0, 1.0, 1.0), 7,
     "Random forest, patient-level\nsplit: ROC on held-out patients"),
    ([0.215, 0.02, 0.19, 0.245], "FigureBaselines_Classification.pdf", (0.294, 0.0, 0.632, 1.0), 8,
     "Repeated splits; edge-weight\nand presence baselines"),
    ([0.415, 0.02, 0.19, 0.245], "Figure2_DifferentialExpression.pdf", (0.5, 0.0, 1.0, 0.465), 9,
     "Differential analysis\nbetween subtypes"),
    ([0.61, 0.02, 0.19, 0.245], "Figure4_Networks_Circular.pdf", (0.5, 0.0, 0.875, 1.0), 10,
     "Decision rules and their\nco-occurrence"),
    ([0.805, 0.02, 0.185, 0.245], "Figure5_Importance.pdf", (0.5, 0.0, 1.0, 0.405), 11,
     "Permutation importance"),
]
for rect, pdf, box, n, lab in panels:
    image(rect, crop(pdf, box))
    step(rect[0] + 0.012, 0.283, n, lab)

fig.savefig(sys.argv[1])
print("written", sys.argv[1])
