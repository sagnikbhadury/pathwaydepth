"""Figure of the Bayesian two-part model (script 08 of the package).

Usage: python 11_twopart_figure.py <revision_results folder> <output pdf>
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, GRID, GREY = "#0b0b0b", "#52514e", "#e6e5e1", "#b9b8b2"
LUSC, LUAD = "#6250d6", "#e34948"            # the subtype colours of Figure 6C

res, out = sys.argv[1], sys.argv[2]
e = pd.read_csv(os.path.join(res, "presence_vs_proximity_effects.csv"))
TABLE2 = {"JAKSTAT_TNFa", "JAKSTAT_TGFb", "Hypoxia_JAKSTAT", "JAKSTAT_p53",
          "TNFa_Trail", "EGFR_TNFa", "EGFR_PI3K", "TGFb_Trail"}

def pretty(u):
    for a, b in [("JAKSTAT", "JAK-STAT"), ("TNFa", "TNFα"), ("TGFb", "TGFβ"),
                 ("NFkB", "NFκB"), ("Trail", "TRAIL")]:
        u = u.replace(a, b)
    return u.replace("_", " ↔ ")

def forest(ax, d, xlabel, title, ticks, mark_table2=False):
    d = d.sort_values("median").reset_index(drop=True)
    y = np.arange(len(d))
    for i, r in d.iterrows():
        col = GREY if not r.excludes_zero else (LUAD if r["median"] > 0 else LUSC)
        ax.plot(np.exp([r.q025, r.q975]), [i, i], color=col, lw=2, solid_capstyle="round")
        ax.plot(np.exp(r["median"]), i, "o", color=col, ms=6, mec="white", mew=0.8, zorder=3)
    ax.axvline(1, color=MUTED, lw=1, ls="--")
    ax.set_xscale("log")
    lab = [pretty(u) + (" *" if mark_table2 and u in TABLE2 else "") for u in d.unit]
    ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=7.5, color=INK)
    ax.set_ylim(-0.7, len(d) - 0.3)
    ax.set_xlabel(xlabel, fontsize=8, color=INK)
    ax.set_title(title, fontsize=9, color=INK, loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelsize=7.5, colors=MUTED, length=2); ax.tick_params(axis="y", length=0)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks(ticks); ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.grid(axis="x", color=GRID, lw=0.8, which="major"); ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color(GRID)

A = e[e.part == "A_presence"]
B = e[(e.part == "B_positive") & e.excludes_zero]
C = e[(e.part == "C_logweight") & e.excludes_zero]

fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1, 1, 1], "wspace": 0.75})
forest(axes[0], A, "Odds ratio of presence, LUAD vs LUSC", "A   Pathway present on the slide", [0.25, 0.5, 1, 2, 4])
forest(axes[1], B, "Odds ratio of a positive weight, LUAD vs LUSC", "B   Positive weight, both present", [0.25, 0.5, 1, 2, 4])
forest(axes[2], C, "Ratio of the positive weight, LUAD vs LUSC", "C   Size of the weight (proximity)", [0.7, 0.8, 1, 1.25, 1.5], mark_table2=True)
h = [plt.Line2D([], [], color=c, marker="o", lw=2, mec="white", label=l) for c, l in
     [(LUAD, "higher in LUAD"), (LUSC, "higher in LUSC"), (GREY, "95% interval includes 1")]]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.01))
fig.subplots_adjust(left=0.12, right=0.99, top=0.92, bottom=0.2)
fig.savefig(out)
print("written", out, "| panels:", len(A), len(B), len(C))
