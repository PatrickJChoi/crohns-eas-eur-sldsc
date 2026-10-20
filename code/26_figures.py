"""Manuscript figures, drawn directly from the result files.

Figure 1: tau* per lineage with 95% CIs for the three CD datasets (registered primary model).
Figure 2: East Asian minus European difference D per lineage: registered marginal and simultaneous (Bonferroni, 8)
          95% CIs, and the post hoc aligned-block intervals (common SNPs) for comparison.
Colours: reference categorical slots 1-3 (validated all-pairs for colour-vision deficiency); each dataset also has its
own marker shape so identity is not carried by colour alone.
Writes results/manuscript/fig1_tau_star.{png,pdf} and fig2_difference.{png,pdf}.
Run from the repo root:  .venv/Scripts/python code/26_figures.py
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

RES = pathlib.Path("results/main_analysis")
OUT = pathlib.Path("results/manuscript")
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SERIES = [("EAS", "East Asian", "#2a78d6", "o"),
          ("EUR", "European (de Lange)", "#eb6834", "s"),
          ("FinnGen", "European (FinnGen)", "#1baf7a", "D")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK, "axes.facecolor": SURF, "figure.facecolor": SURF,
                     "savefig.facecolor": SURF})


def style(ax):
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.axvline(0, color=INK2, linewidth=1)


def fig1(t):
    lineages = t["lineage"].tolist()[::-1]
    tt = t.set_index("lineage").loc[lineages]
    y = np.arange(len(lineages))
    fig, ax = plt.subplots(figsize=(6.3, 4.2))
    for i, (key, label, col, mk) in enumerate(SERIES):
        est = tt[f"tau_star_{key if key != 'EUR' else 'EUR'}"] if key != "FinnGen" else tt["tau_star_FinnGen"]
        se = tt[f"se_{key}"] if key != "FinnGen" else tt["se_FinnGen"]
        yy = y + (1 - i) * 0.22
        ax.errorbar(est, yy, xerr=1.96 * se, fmt=mk, color=col, ms=6, mfc=col, mec=SURF, mew=1.2,
                    elinewidth=2, capsize=0, label=label)
    ax.set_yticks(y, lineages)
    ax.set_xlabel("τ* (standardized per-SNP heritability effect of lineage-specific genes)")
    style(ax)
    # legend in one row above the plot so it never covers data; the title lives in the manuscript caption
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="upper center", ncol=3, fontsize=8,
               labelcolor=INK, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig1_tau_star.{ext}", dpi=300)
    plt.close(fig)


def fig2(t, aligned):
    lineages = t["lineage"].tolist()[::-1]
    tt = t.set_index("lineage").loc[lineages]
    y = np.arange(len(lineages))
    zs = stats.norm.ppf(1 - 0.05 / (2 * len(t)))
    fig, ax = plt.subplots(figsize=(6.3, 4.2))
    D, se = tt["D_EAS_minus_EUR"].values, tt["se_D"].values
    ax.hlines(y + 0.12, D - zs * se, D + zs * se, color="#2a78d6", linewidth=1, alpha=0.6,
              label="registered, simultaneous 95% CI")
    ax.errorbar(D, y + 0.12, xerr=1.96 * se, fmt="o", color="#2a78d6", ms=6, mec=SURF, mew=1.2, elinewidth=2.5,
                capsize=0, label="registered, 95% CI")
    a = aligned["EAS_vs_EUR_deLange"]["direct_aligned"]
    order = [t["lineage"].tolist().index(L) for L in lineages]
    Da, sa = np.array(a["D"])[order], np.array(a["se_D"])[order]
    ax.errorbar(Da, y - 0.12, xerr=1.96 * sa, fmt="o", color=INK2, ms=5, mfc=SURF, mew=1.4, elinewidth=1.5,
                capsize=0, label="post hoc aligned-block, 95% CI")
    ax.set_yticks(y, lineages)
    ax.set_xlabel("D = τ* East Asian − τ* European")
    style(ax)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="upper center", ncol=2, fontsize=8,
               labelcolor=INK, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig2_difference.{ext}", dpi=300)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(RES / "primary_table.tsv", sep="\t")
    fig1(t)
    fig2(t, json.load(open(RES / "aligned_block_check.json")))
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
