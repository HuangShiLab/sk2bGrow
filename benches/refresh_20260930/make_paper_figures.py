#!/usr/bin/env python3
"""Generate publication figures from archived refresh and mixed-strain TSVs."""
from __future__ import annotations
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REF = Path(__file__).resolve().parent
MIX = ROOT / "benches/mixedstrain_20261001"
FIG = ROOT / "docs/paper/figures"
FIG.mkdir(parents=True, exist_ok=True)
mpl.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "figure.dpi": 180, "savefig.dpi": 300, "axes.spines.top": False,
    "axes.spines.right": False,
})
ARM_LABEL = {"A": "anchors + V-fit", "E": "FracMinHash + V-fit", "B": "anchors + sorted",
             "C_relaxed": "Pilea gates-off", "C_default": "Pilea defaults"}
ARM_COLOR = {"A": "#1f77b4", "E": "#d62728", "B": "#7f7f7f",
             "C_relaxed": "#9467bd", "C_default": "#8c564b"}

def save(fig, stem):
    for suffix in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{suffix}", bbox_inches="tight")
    plt.close(fig)

def add_ci(ax, x, lo, hi, color):
    y = np.vstack([lo, hi])
    ax.fill_between(x, y[0], y[1], color=color, alpha=0.16, lw=0)

def pipeline_figure():
    fig, ax = plt.subplots(figsize=(7.2, 2.35))
    ax.axis("off")
    boxes = [
        (0.015, "reference genomes\n16 Type-II/IIS enzymes"),
        (0.235, "TGT / anchor DB\ndeterministic loci"),
        (0.455, "reads\nWGS or route-B"),
        (0.655, "integer counts\n+ EM diagnostics"),
        (0.855, "window rates\nV-fit + QC"),
    ]
    for x, label in boxes:
        ax.add_patch(plt.Rectangle((x, 0.36), 0.17, 0.30, facecolor="#e8f0fb",
                                   edgecolor="#2f5597", lw=1))
        ax.text(x + 0.085, 0.51, label, ha="center", va="center")
    for x in [0.185, 0.405, 0.625, 0.825]:
        ax.annotate("", xy=(x + 0.03, 0.51), xytext=(x, 0.51),
                    arrowprops=dict(arrowstyle="->", lw=1.1, color="#444444"))
    ax.text(0.5, 0.15, "Per-stratum windows; shared anchors retained for diagnostics but excluded from integer ZTP PTR fitting",
            ha="center", fontsize=7, style="italic")
    save(fig, "fig1_pipeline")

def zheng_figure():
    qc = pd.read_csv(REF / "results/zheng_all_seeds_qc.tsv", sep="\t")
    boot = pd.read_csv(REF / "results/zheng_bootstrap_media.tsv", sep="\t")
    delta = pd.read_csv(REF / "results/zheng_arm_delta_bootstrap.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 2.9), constrained_layout=True)
    ax = axes[0]
    show = ["A", "E", "B", "C_relaxed"]
    for arm in show:
        d = qc[(qc["arm"] == arm) & (qc["view"] == "all finite estimates")]
        agg = d.groupby("coverage", as_index=False)["pearson_r"].agg(["mean", "min", "max"])
        ax.plot(agg.coverage, agg["mean"], "-o", color=ARM_COLOR[arm], label=ARM_LABEL[arm], ms=3)
        ax.fill_between(agg.coverage, agg["min"], agg["max"], color=ARM_COLOR[arm], alpha=0.10, lw=0)
    ax.set(xscale="log", xlabel="Sequencing coverage", ylabel="Pearson r vs measured λ",
           title="A  Zheng estimator accuracy", xticks=[0.5, 1, 2, 5, 10])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.legend(frameon=False, loc="lower right")
    ax = axes[1]
    a = qc[(qc["arm"] == "A") & (qc["view"] == "default QC passed")]
    agg = a.groupby("coverage", as_index=False)["n_pass_qc"].mean()
    ax.bar(agg.coverage.astype(str), agg.n_pass_qc, color="#1f77b4")
    ax.set(ylim=(0, 16), ylabel="Mean default-QC passes /16", xlabel="Sequencing coverage",
           title="B  Shipped-QC deployability")
    ax = axes[2]
    d = delta[delta.comparison == "A_minus_E"]
    x = np.arange(len(d))
    ax.errorbar(x, d.pearson_r_delta,
                yerr=[d.pearson_r_delta-d.pearson_r_ci_low, d.pearson_r_ci_high-d.pearson_r_delta],
                fmt="o", color="#1f77b4", capsize=3, ms=4)
    ax.axhline(0, color="#444444", lw=0.7)
    ax.set(xticks=x, xticklabels=[f"{v:g}×" for v in d.coverage],
           xlabel="Sequencing coverage", ylabel="Δ Pearson r (A − E)",
           title="C  Paired anchor advantage")
    save(fig, "fig2_zheng_performance")

def a4_figure():
    slopes = pd.read_csv(REF / "results/a4_slope_bootstrap.tsv", sep="\t")
    alls = pd.read_csv(REF / "a4/results/genome_slopes_by_depth.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
    ax = axes[0]
    colors = {"pois": "#1f77b4", "nb": "#d62728"}
    labels = {"pois": "Poisson", "nb": "overdispersed NB"}
    for arm in ["pois", "nb"]:
        d = slopes[slopes.arm == arm].sort_values("depth")
        e = alls[(alls.arm == arm) & (alls.variant == "true_all_windows")].sort_values("depth")
        ax.errorbar(d.depth, d.slope, yerr=[d.slope-d.slope_ci_low, d.slope_ci_high-d.slope],
                    fmt="-o", color=colors[arm], label=f"recovered: {labels[arm]}", capsize=2, ms=3)
        ax.plot(e.depth, e.genome_slope_est_on_truth, "--", color=colors[arm], alpha=0.7,
                label=f"exact rates: {labels[arm]}")
    ax.axhline(1, color="#444444", lw=0.7)
    ax.set(xscale="log", xticks=[0.5, 1, 2, 5], xlabel="Sequencing coverage",
           ylabel="Recovered / true log2PTR slope", title="A  Low-depth compression")
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax = axes[1]
    for arm in ["pois", "nb"]:
        d = slopes[slopes.arm == arm].sort_values("depth")
        ax.errorbar(d.depth, d.b0_floor,
                    yerr=[d.b0_floor-d.b0_ci_low, d.b0_ci_high-d.b0_floor],
                    fmt="-o", color=colors[arm], label=labels[arm], capsize=2, ms=3)
    ax.axhline(0, color="#444444", lw=0.7)
    ax.set(xscale="log", xticks=[0.5, 1, 2, 5], xlabel="Sequencing coverage",
           ylabel="Mean estimate at true b=0", title="B  Stationary floor")
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.legend(frameon=False)
    save(fig, "fig3_a4_compression")

def r3_figure():
    stats = pd.read_csv(REF / "r3/r3_stats.tsv", sep="\t")
    sens = pd.read_csv(REF / "results/r3_threshold_sensitivity.tsv", sep="\t")
    stats["family"] = np.where(stats.ref.str.contains("scr"), "scrambled",
                      np.where(stats.ref.str.contains("ord"), "ordered", "complete"))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
    ax = axes[0]
    colors = {1: "#4c78a8", 5: "#f58518", 10: "#54a24b"}
    for depth in [1, 5, 10]:
        d = (stats[(stats.family == "scrambled") & (stats.depth == depth)]
             .groupby("log2ptr", as_index=False)["flag"].mean())
        ax.plot(d.log2ptr, d.flag, "-o", label=f"{depth}×", color=colors[depth], ms=4)
    ax.set(xticks=[0, 0.5, 1, 1.5], xlabel="True log2PTR", ylabel="WCG detection fraction",
           ylim=(-0.03, 1.03), title="A  Scrambled references")
    ax.legend(frameon=False, title="depth")
    ax = axes[1]
    ax.plot(sens.z_short_threshold, sens.stationary_false_fires, "-o", label="stationary false fires",
            color="#d62728")
    ax.plot(sens.z_short_threshold, sens.strong_5_10x_fires, "-o", label="strong 5–10× detections",
            color="#1f77b4")
    ax.axvline(2.5, color="#666666", ls=":", lw=0.8)
    ax.text(2.52, 36, "post-hoc\ncandidate", fontsize=6, color="#666666")
    ax.set(xlabel="Slope-branch z-short threshold", xticks=sens.z_short_threshold,
           ylabel="Runs", title="B  Threshold sensitivity")
    ax.legend(frameon=False)
    save(fig, "fig4_r3_qc")

def mixed_figure():
    path = MIX / "results" / "sim_summary.tsv"
    if not path.exists():
        path = MIX / "sim_results.tsv"
    if not path.exists():
        return
    d = pd.read_csv(path, sep="\t")
    if d.empty or "arm" not in d:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
    palette = {4: "#4c78a8", 8: "#f58518", 16: "#54a24b"}
    for ns, grp in d.groupby("n_strains"):
        agg = grp.groupby("coverage", as_index=False).agg(recall=("recall", "mean"),
                                                          rmse=("rmse", "mean"))
        ax = axes[0]
        ax.plot(agg.coverage, agg.recall, "-o", label=f"{ns} strains", color=palette.get(ns), ms=4)
        ax = axes[1]
        ax.plot(agg.coverage, agg.rmse, "-o", label=f"{ns} strains", color=palette.get(ns), ms=4)
    axes[0].set(xlabel="Coverage per genome", ylabel="Recall", ylim=(-0.03, 1.03),
                title="A  Reported-strain recall", xscale="log")
    axes[1].set(xlabel="Coverage per genome", ylabel="RMSE (log2PTR)",
                title="B  Estimate error", xscale="log")
    for ax in axes:
        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        ax.set_xticks([1, 2, 4])
        ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
        ax.legend(frameon=False)
    save(fig, "fig5_mixed_strain")

pipeline_figure()
zheng_figure()
a4_figure()
r3_figure()
mixed_figure()
print(f"wrote figures to {FIG}")
