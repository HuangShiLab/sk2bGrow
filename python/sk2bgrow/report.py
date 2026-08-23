"""Output assembly and QC figures.

The main table keeps Pilea's column names — ``coverage``, ``dispersion``,
``fraction``, ``containment``, ``PTR``, ``log2(PTR)`` — so an existing benchmark
script can swap tools without rewriting its parser, which is exactly what the P0
and P1 validation stages need. Three columns are appended:
``enzyme_consistency``, ``n_anchors``, ``ori_confidence``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import io as sk_io

__all__ = ["OUTPUT_COLUMNS", "assemble", "apply_qc", "write_output", "plot_profile", "plot_forest", "write_qc_figures"]

#: Column order of ``output.tsv``. The first block is Pilea-compatible.
OUTPUT_COLUMNS = [
    "sample",
    "genome",
    "taxonomy",
    "coverage",
    "dispersion",
    "fraction",
    "containment",
    "PTR",
    "log2(PTR)",
    # --- sk2bGrow additions ------------------------------------------------
    "enzyme_consistency",
    "n_anchors",
    "ori_confidence",
    "se",
    "ci_low",
    "ci_high",
    "n_enzymes",
    "n_enzymes_attempted",
    "enzyme_fit_rate",
    "enzymes_used",
    "enzyme_i2",
    "fusion_model",
    "ori",
    "method",
    "n_windows",
    "pass_qc",
    "qc_reason",
    "excluded",
    "note",
]


def assemble(fused: pd.DataFrame, manifest: sk_io.DbManifest, stats: dict[str, dict] | None = None) -> pd.DataFrame:
    """Join the fused estimates with manifest and per-sample counting stats.

    ``containment`` comes from the EM step in the Rust layer (fraction of a
    genome's unique anchors that were detected at all), which is the same
    quantity Pilea's ANI filter thresholds on.
    """
    df = fused.copy()
    df["taxonomy"] = [manifest.genome(int(g)).taxonomy or "-" for g in df["genome_id"]]
    df["PTR"] = df["ptr"]
    df["log2(PTR)"] = df["log2_ptr"]

    containment = np.full(len(df), np.nan)
    if stats:
        for i, (sample, gid) in enumerate(zip(df["sample"], df["genome_id"])):
            entry = stats.get(str(sample))
            if not entry:
                continue
            for g in entry.get("em", {}).get("genomes", []):
                if int(g.get("genome_id", -1)) == int(gid):
                    containment[i] = float(g.get("containment", np.nan))
                    break
    df["containment"] = containment
    return df


def apply_qc(
    df: pd.DataFrame,
    min_coverage: float = 1.0,
    min_fraction: float = 0.75,
    max_dispersion: float = 5.0,
    min_containment: float = 0.5,
    consistency_alpha: float = 0.05,
    min_ori_confidence: float = 0.0,
    min_enzyme_fit_rate: float = 0.8,
) -> pd.DataFrame:
    """Flag rows that fail quality gates.

    Nothing is deleted: a failing row keeps its estimate and gains a reason. A
    tool that silently drops rows makes "this genome was not growing" and "this
    genome was filtered out" indistinguishable downstream.

    ``min_enzyme_fit_rate`` defaults to 0.8: on a 16-enzyme panel that tolerates
    three sparse enzymes dropping out (PpiI, PsrI and BplI are legitimately thin
    on some genomes) while flagging a sample where a quarter of the panel saw no
    gradient at all.

    The coverage floor defaults to 1x rather than Pilea's 5x. That is the point
    of the whole design — the report's simulations put the deterministic union's
    usable boundary near 1x, not 0.2x — but it is a *default*, and a study that
    wants Pilea-comparable strictness should pass ``min_coverage=5``.
    """
    out = df.copy()
    reasons: list[str] = []
    for _, r in out.iterrows():
        why = []
        cov = r.get("coverage", np.nan)
        if np.isfinite(cov) and cov < min_coverage:
            why.append(f"coverage {cov:.2f} < {min_coverage}")
        frac = r.get("fraction", np.nan)
        if np.isfinite(frac) and frac < min_fraction:
            why.append(f"fraction {frac:.2f} < {min_fraction}")
        disp = r.get("dispersion", np.nan)
        if np.isfinite(disp) and disp > max_dispersion:
            why.append(f"dispersion {disp:.2f} > {max_dispersion}")
        cont = r.get("containment", np.nan)
        if np.isfinite(cont) and cont < min_containment:
            why.append(f"containment {cont:.2f} < {min_containment}")
        cons = r.get("enzyme_consistency", np.nan)
        if np.isfinite(cons) and cons < consistency_alpha:
            why.append(f"enzymes disagree (p={cons:.2g})")
        rate = r.get("enzyme_fit_rate", np.nan)
        if np.isfinite(rate) and rate < min_enzyme_fit_rate:
            why.append(
                f"only {int(r.get('n_enzymes', 0))}/{int(r.get('n_enzymes_attempted', 0))} enzymes produced a fit"
            )
        oc = r.get("ori_confidence", np.nan)
        if min_ori_confidence > 0 and np.isfinite(oc) and oc < min_ori_confidence:
            why.append(f"ori_confidence {oc:.2f} < {min_ori_confidence}")
        if not np.isfinite(r.get("log2_ptr", np.nan)):
            why.append("no PTR estimate")
        reasons.append("; ".join(why) if why else "")
    out["qc_reason"] = reasons
    out["pass_qc"] = [r == "" for r in reasons]
    return out


def write_output(df: pd.DataFrame, path: str | Path) -> Path:
    """Write ``output.tsv`` in the canonical column order."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    for c in OUTPUT_COLUMNS:
        if c not in out.columns:
            out[c] = np.nan
    out = out[OUTPUT_COLUMNS].sort_values(["sample", "genome"])
    # 10 significant digits, not 6: the ``ori`` column holds a genome coordinate,
    # and %.6g would round 3,923,883 to 3,923,880. Statistics do not need the
    # extra digits, but coordinates do.
    out.to_csv(path, sep="\t", index=False, float_format="%.10g", na_rep="NA")
    return path


# --------------------------------------------------------------------------
# figures (matplotlib is optional; missing it degrades to "no figures")
# --------------------------------------------------------------------------

def _plt():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except Exception:  # pragma: no cover - depends on the environment
        return None


def plot_profile(windows: pd.DataFrame, fits: pd.DataFrame, sample: str, genome: str, path: str | Path) -> Path | None:
    """Coverage profile: log2 window rate against genome coordinate, per enzyme.

    This is the figure that shows whether a PTR number means anything — a real
    replication gradient is visible as a V, and 16 enzymes tracing the same V is
    the visual form of the consistency test.
    """
    plt = _plt()
    if plt is None:
        return None
    sel = windows[(windows["sample"] == sample) & (windows["genome"] == genome)]
    if sel.empty:
        return None
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for enzyme, grp in sel.groupby("enzyme", sort=True):
        ax.scatter(grp["global_mid"] / 1e6, grp["log2_rate"], s=8, alpha=0.55, label=enzyme)
    frow = fits[(fits["sample"] == sample) & (fits["genome"] == genome)]
    if not frow.empty and np.isfinite(frow["ori"].iloc[0]):
        ori = float(frow["ori"].iloc[0])
        ax.axvline(ori / 1e6, color="crimson", lw=1.2, ls="--", label="ori")
    ax.set_xlabel("genome coordinate (Mb)")
    ax.set_ylabel("log2 window rate")
    ax.set_title(f"{sample} / {genome}")
    ax.legend(fontsize=6, ncol=4, frameon=False)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_forest(per_enzyme: pd.DataFrame, fused_row: pd.Series, path: str | Path) -> Path | None:
    """Forest plot of the 16 per-enzyme estimates against the fused value."""
    plt = _plt()
    if plt is None:
        return None
    sel = per_enzyme[
        (per_enzyme["sample"] == fused_row["sample"]) & (per_enzyme["genome"] == fused_row["genome"])
    ].sort_values("enzyme")
    if sel.empty:
        return None
    fig, ax = plt.subplots(figsize=(6, 0.35 * len(sel) + 2))
    y = np.arange(len(sel))
    ax.errorbar(sel["log2_ptr"], y, xerr=sel["se"], fmt="o", ms=4, capsize=2, lw=1)
    ax.axvline(fused_row["log2_ptr"], color="crimson", lw=1.4, label="fused")
    if np.isfinite(fused_row.get("ci_low", np.nan)):
        ax.axvspan(fused_row["ci_low"], fused_row["ci_high"], color="crimson", alpha=0.12)
    ax.set_yticks(y)
    ax.set_yticklabels(sel["enzyme"])
    ax.set_xlabel("log2(PTR)")
    p = fused_row.get("enzyme_consistency", np.nan)
    ax.set_title(f"{fused_row['sample']} / {fused_row['genome']}  (Q p={p:.3g})" if np.isfinite(p) else f"{fused_row['sample']} / {fused_row['genome']}")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def write_qc_figures(windows: pd.DataFrame, per_enzyme: pd.DataFrame, fused: pd.DataFrame, outdir: str | Path) -> list[Path]:
    """Write one profile and one forest plot per (sample, genome)."""
    outdir = Path(outdir)
    written: list[Path] = []
    for _, row in fused.iterrows():
        tag = f"{row['sample']}__{row['genome']}".replace("/", "_")
        p = plot_profile(windows, fused, row["sample"], row["genome"], outdir / f"{tag}.profile.png")
        if p:
            written.append(p)
        f = plot_forest(per_enzyme, row, outdir / f"{tag}.forest.png")
        if f:
            written.append(f)
    return written
