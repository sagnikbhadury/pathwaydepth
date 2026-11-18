"""Figure of the fusion analysis (scripts 14 and 15 of the package).

Usage: python 16_fusion_figure.py <revision_results folder> <output pdf>
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, GRID, GREY = "#0b0b0b", "#52514e", "#e6e5e1", "#b9b8b2"
COL = {"raw": "#2a78d6", "popdepth": "#eb6834", "within": "#1baf7a"}   # as in Figure 3
FUSED = "#3d3c39"
LAB = {"raw": "Edge weights", "within": "Patient Referenced Depth", "popdepth": "Population Referenced Depth",
       "all3": "Early fusion, all three", "raw_within": "Early fusion, weights + Patient",
       "raw_popdepth": "Early fusion, weights + Population", "within_popdepth": "Early fusion, Patient + Population",
       "late_mean": "Late fusion, mean probability", "late_stack": "Late fusion, stacking"}

res, out = sys.argv[1], sys.argv[2]
rs = pd.read_csv(os.path.join(res, "fusion_repeated_summary.csv")).set_index("model")
pp = pd.read_csv(os.path.join(res, "fusion_performance_primary.csv")).set_index("model")
pat = pd.read_csv(os.path.join(res, "fusion_correct_patterns.csv"), dtype={"pattern": str})
share = pd.read_csv(os.path.join(res, "fusion_importance_share.csv")).set_index("representation")
up = pd.read_csv(os.path.join(res, "fusion_unique_pairs.csv"))


def pretty(u):
    for a, b in [("JAKSTAT", "JAK-STAT"), ("TNFa", "TNFα"), ("TGFb", "TGFβ"), ("NFkB", "NFκB"), ("Trail", "TRAIL")]:
        u = u.replace(a, b)
    return u.replace("_", " ↔ ")


def style(ax, title):
    ax.set_title(title, fontsize=9, color=INK, loc="left", fontweight="bold")
    ax.tick_params(labelsize=7.5, colors=MUTED, length=2)
    for s in ax.spines.values():
        s.set_color(GRID)


fig, axes = plt.subplots(2, 2, figsize=(11, 7.6), gridspec_kw={"hspace": 0.55, "wspace": 0.75})

# A: test AUC over repeated splits, single and fused
ax = axes[0, 0]
order = ["raw", "within", "popdepth", "raw_within", "raw_popdepth", "within_popdepth", "all3", "late_mean", "late_stack"]
for i, k in enumerate(order):
    y = len(order) - 1 - i
    c = COL.get(k, FUSED)
    ax.plot([rs.loc[k, "q025"], rs.loc[k, "q975"]], [y, y], color=c, lw=2, solid_capstyle="round")
    ax.plot(rs.loc[k, "mean"], y, "o", color=c, ms=6, mec="white", mew=0.8, zorder=3)
    ax.text(0.905, y, f"{rs.loc[k, 'mean']:.3f}", va="center", fontsize=7.5, color=INK)
ax.axvline(0.5, color=MUTED, lw=1, ls="--")
ax.text(0.505, -0.85, "chance", fontsize=7, color=MUTED)
ax.set_ylim(-1.1, len(order) - 0.4)
ax.axhline(5.5, color=GRID, lw=0.8)
ax.set_yticks(range(len(order))); ax.set_yticklabels([LAB[k] for k in order][::-1], fontsize=7.5, color=INK)
ax.set_xlim(0.45, 0.94); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
n_rep = len(pd.read_csv(os.path.join(res, "fusion_repeated.csv")))
ax.set_xlabel(f"Test AUC over {n_rep} repeated patient-level splits (mean, 2.5th-97.5th pct.)",
              fontsize=8, color=INK)
ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True); ax.tick_params(axis="y", length=0)
style(ax, "A   Single representations and their fusion")

# B: which representation classifies each test patient correctly
ax = axes[0, 1]
names = {"111": "all three", "110": "weights + Patient", "101": "weights + Population",
         "100": "weights only", "011": "both depths only", "010": "Patient only", "001": "Population only",
         "000": "none"}
pat = pat.set_index("pattern").reindex(list(names)).fillna(0)
cols = [FUSED, GREY, GREY, COL["raw"], "#8a8985", COL["within"], COL["popdepth"], "#dddcd7"]
yy = np.arange(len(names))[::-1]
ax.barh(yy, pat.Freq.values, color=cols, height=0.62)
for y, v in zip(yy, pat.Freq.values):
    ax.text(v + 2, y, f"{int(v)}", va="center", fontsize=7.5, color=INK)
ax.set_yticks(yy); ax.set_yticklabels(list(names.values()), fontsize=7.5, color=INK)
ax.set_xlabel(f"Test patients of the primary split (n = {int(pat.Freq.sum())})", fontsize=8, color=INK)
ax.set_xlim(0, pat.Freq.max() * 1.15); ax.tick_params(axis="y", length=0)
ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
style(ax, "B   Correctly classified by")

# C: importance share in the fused forest
ax = axes[1, 0]
reps = [("w", "raw", "Edge weights"), ("pat", "within", "Patient Referenced Depth"),
        ("pop", "popdepth", "Population Referenced Depth")]
yy = np.arange(3)[::-1]
for y, (r, k, l) in zip(yy, reps):
    v = share.loc[r, "importance_share"]
    ax.barh(y, 100 * v, color=COL[k], height=0.55)
    ax.text(100 * v + 1, y, f"{100 * v:.0f}%   ({int(share.loc[r, 'in_top30'])} of the top 30 features)",
            va="center", fontsize=7.5, color=INK)
ax.set_yticks(yy); ax.set_yticklabels([l for _, _, l in reps], fontsize=7.5, color=INK)
ax.set_xlim(0, 100); ax.set_xlabel("Share of the permutation importance of the early-fusion forest (%)",
                                   fontsize=8, color=INK)
ax.tick_params(axis="y", length=0); ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
style(ax, "C   What the fused forest uses")

# D: pairs distinguished by one depth representation only
ax = axes[1, 1]
d = up[up.pattern.isin(["-P-", "--C"])].copy().sort_values(["pattern", "patient_depth_fdr"], ascending=[False, True])
cols_d = [("weight_d", "weight_fdr", COL["raw"], "Edge weight", "o"),
          ("signed_rank_d", "signed_rank_fdr", FUSED, "Signed within-slide rank", "D"),
          ("patient_depth_d", "patient_depth_fdr", COL["within"], "Patient Referenced Depth", "s"),
          ("population_depth_d", "population_depth_fdr", COL["popdepth"], "Population Referenced Depth", "^")]
yy = np.arange(len(d))[::-1]
off = np.linspace(0.27, -0.27, len(cols_d))
for (c, f, colr, lab, m), o in zip(cols_d, off):
    sig = d[f].values < 0.05
    ax.scatter(d[c].values[sig], yy[sig] + o, color=colr, marker=m, s=30, zorder=3, label=lab)
    ax.scatter(d[c].values[~sig], yy[~sig] + o, facecolor="white", edgecolor=colr, marker=m, s=30, zorder=3)
ax.axvline(0, color=MUTED, lw=1, ls="--")
ax.set_yticks(yy); ax.set_yticklabels([pretty(u) for u in d.feature], fontsize=7.5, color=INK)
ax.set_xlabel("Cohen's d, LUAD minus LUSC (filled: FDR < 0.05)", fontsize=8, color=INK)
ax.tick_params(axis="y", length=0); ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
style(ax, "D   Pairs distinguished by a depth only")

fig.subplots_adjust(left=0.2, right=0.8, top=0.95, bottom=0.08)
fig.savefig(out)
print("written", out)
