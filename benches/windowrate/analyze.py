#!/usr/bin/env python3
"""C7/A4 analysis: recovered window rates vs exactly known truth.

Reads every sample under ``work/out/`` (GC control: ``work/out_gc/``),
reconstructs the exact expected count of every anchor from the simulation
parameters (see ``simulate_tent.py``), aggregates it to per-window truth using
the window boundaries the pipeline itself wrote, and answers:

  (a) is there a floor/ceiling in recovered rates?
  (b) regression of log2(recovered) on log2(true) per depth — the
      *window-rate slope*, which is NOT the genome-level slope of A4;
  (c) does the BIC-chosen component count / model branch track the bias?
  (d) genome-level slope of estimated log2PTR on true b, three ways:
      recovered rates (the pipeline as shipped), true rates at the same
      surviving windows (isolates rate bias from window dropout), and true
      rates at all windows (the no-window-layer baseline). All three use the
      identical fit + fusion machinery, so differences are attributable.

Two distinct slopes are reported and labelled throughout, per RESEARCH_PLAN
Part 2.2: the per-window ``log2(recovered) ~ log2(true)`` slope, and the
genome-level ``estimated log2PTR ~ true b`` slope.

Outputs TSVs under ``results/`` and prints a compact summary.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from simulate_tent import ORI, READ_LEN, anchor_expected, start_weights  # noqa: E402

from sk2bgrow import fit as sk_fit  # noqa: E402
from sk2bgrow import fusion as sk_fusion  # noqa: E402
from sk2bgrow import io as sk_io  # noqa: E402
from sk2bgrow.ztp import estimate_window_rate  # noqa: E402

#: Tag lengths from crates/sk2bgrow-core/src/enzyme.rs PANEL (digest position is
#: the 0-based forward-strand start of the tag; the tag spans tag_len bases).
TAG_LEN = {
    "BcgI": 32, "AlfI": 32, "AloI": 27, "BaeI": 28, "BplI": 27, "BsaXI": 27,
    "BslFI": 25, "Bsp24I": 27, "CjeI": 28, "CjePI": 27, "CspCI": 33, "FalI": 27,
    "HaeIV": 27, "Hin4I": 27, "PpiI": 27, "PsrI": 27,
}

SAMPLE_RE = re.compile(r"(pois|nb)_d([\d.]+)_b([\d.]+)_r(\d+)$")

#: True-rate bins for the near-zero stratification, in expected counts/anchor.
RATE_BINS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, np.inf]
RATE_LABELS = ["<0.25", "0.25-0.5", "0.5-1", "1-2", "2-4", ">=4"]


def parse_sample(name: str, gc_on: bool) -> dict | None:
    m = SAMPLE_RE.match(name)
    if not m:
        return None
    arm, depth, b, rep = m.group(1), float(m.group(2)), float(m.group(3)), int(m.group(4))
    return {"sample": name, "arm": arm + ("_gc" if gc_on else ""), "depth": depth,
            "true_log2_ptr": b, "rep": rep}


def window_truth(counts: pd.DataFrame, lam: np.ndarray, windows: pd.DataFrame, check: bool = True) -> pd.DataFrame:
    """Mean true expected count over each window's usable anchors.

    Window membership is recovered from the ``start``/``end`` (min/max anchor
    position) the pipeline wrote for each (enzyme, window) — blocks are
    consecutive in position order within an enzyme, so the ranges partition the
    usable anchors exactly. The member count is checked against the pipeline's
    own ``n_anchors``.
    """
    df = counts[counts["usable"]].copy()
    df["true_lambda"] = lam[counts["usable"].to_numpy()]
    out = []
    for enzyme, wgrp in windows.groupby("enzyme", sort=True):
        agrp = df[df["enzyme"] == enzyme].sort_values("position")
        pos = agrp["position"].to_numpy()
        wg = wgrp.sort_values("start")
        starts = wg["start"].to_numpy()
        ends = wg["end"].to_numpy()
        idx = np.searchsorted(ends, pos, side="left")
        if not (np.all(idx < len(wg)) and np.all(starts[idx] <= pos)):
            raise RuntimeError(f"window ranges do not partition usable anchors for {enzyme}")
        truth = pd.Series(agrp["true_lambda"].to_numpy()).groupby(idx).mean().to_numpy()
        n = pd.Series(pos).groupby(idx).size().to_numpy()
        if len(truth) != len(wg):
            raise RuntimeError(f"window count mismatch for {enzyme}")
        wg = wg.copy()
        wg["true_rate"] = truth
        if not np.array_equal(n, wg["n_anchors"].to_numpy()):
            raise RuntimeError(f"member count mismatch for {enzyme}")
        # Refit each window from its member counts with the same public function
        # the pipeline used (``estimate_window_rate``, model "auto"). The TSV
        # does not carry n_components/BIC, so this recovers them; the refit rate
        # must equal the pipeline's, which doubles as a consistency check.
        anchor_counts = agrp["count"].to_numpy()
        refit = [estimate_window_rate(anchor_counts[idx == j], model="auto") for j in range(len(wg))]
        wg["refit_rate"] = [r.rate for r in refit]
        wg["n_components"] = [r.n_components for r in refit]
        wg["bic"] = [r.bic for r in refit]
        wg["model"] = [r.model for r in refit]
        same = np.isfinite(wg["rate"]) & np.isfinite(wg["refit_rate"])
        if check and not np.allclose(wg["rate"][same], wg["refit_rate"][same], rtol=1e-9, atol=1e-12):
            raise RuntimeError(f"refit disagrees with pipeline rates for {enzyme}")
        out.append(wg)
    return pd.concat(out, ignore_index=True)


def estimate_genome_level(windows: pd.DataFrame, manifest, variant: str) -> pd.DataFrame:
    """Run the real fit + fusion on a window table; one row per sample."""
    pe = sk_fit.fit_windows(windows, manifest, method="auto")
    fused = sk_fusion.fuse_table(pe)
    fused = fused[["sample", "log2_ptr", "se", "n_enzymes", "enzyme_fit_rate", "enzyme_consistency"]]
    fused["variant"] = variant
    return fused


def ols(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, int]:
    """Slope, intercept, Pearson r, n for finite pairs."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3 or np.var(x) < 1e-12:
        return np.nan, np.nan, np.nan, int(x.size)
    slope, intercept = np.polyfit(x, y, 1)
    r = float(np.corrcoef(x, y)[0, 1])
    return float(slope), float(intercept), r, int(x.size)


def main() -> None:
    root = Path(__file__).resolve().parent
    work = root / "work"
    results = root / "results"
    results.mkdir(exist_ok=True)

    manifest = sk_io.read_manifest(work / "db")
    ginfo = next(iter(manifest.genomes.values()))
    L = ginfo.genome_len

    per_window_rows = []
    for gc_on, outroot in ((False, work / "out"), (True, work / "out_gc")):
        if not outroot.is_dir():
            continue
        for outdir in sorted(outroot.iterdir()):
            meta = parse_sample(outdir.name, gc_on)
            if meta is None or not (outdir / "windows.rates.tsv").exists():
                continue
            sample = meta["sample"]
            counts = sk_io.read_counts(outdir / f"{sample}.counts.tsv")
            windows = pd.read_csv(outdir / "windows.rates.tsv", sep="\t")

            eff_path = work / "eff" / f"{sample}.npy"
            eff = np.load(eff_path) if eff_path.exists() else None
            w = start_weights(L, ORI, meta["true_log2_ptr"], eff)
            n_reads = int(round(meta["depth"] * L / READ_LEN))
            lam = anchor_expected(
                w, n_reads,
                counts["position"].to_numpy(),
                counts["enzyme"].map(TAG_LEN).to_numpy(),
            )

            win = window_truth(counts, lam, windows, check=not gc_on)
            for k, v in meta.items():
                win[k] = v if k != "sample" else (v + ("__gc" if gc_on else ""))
            per_window_rows.append(win)
            print(f"  {win['sample'].iloc[0]}: {len(win)} windows", file=sys.stderr)

    perwin = pd.concat(per_window_rows, ignore_index=True)
    perwin["log2_true"] = np.where(perwin["true_rate"] > 0, np.log2(perwin["true_rate"]), np.nan)
    perwin["log2_rec"] = np.where(perwin["rate"] > 0, np.log2(perwin["rate"]), np.nan)
    perwin["log2_ratio"] = perwin["log2_rec"] - perwin["log2_true"]
    perwin.to_csv(results / "per_window.tsv", sep="\t", index=False, float_format="%.10g")

    # ------------------------------------------------------------------
    # genome-level estimates: recovered vs the two true-rate counterfactuals
    # ------------------------------------------------------------------
    fit_cols = ["sample", "genome_id", "genome", "enzyme", "window", "global_mid",
                "rate", "se", "log2_rate", "log2_se", "detected_fraction", "dispersion",
                "n_anchors", "true_rate"]
    est_frames = []
    for variant in ("recovered", "true_same_windows", "true_all_windows"):
        frames = []
        for sample, grp in perwin.groupby("sample", sort=True):
            wf = grp[fit_cols].copy()
            if variant != "recovered":
                wf["rate"] = wf["true_rate"]
                wf["log2_rate"] = np.where(wf["true_rate"] > 0, np.log2(wf["true_rate"]), np.nan)
                if variant == "true_same_windows":
                    # Rate bias only: windows the ZTP layer dropped stay dropped.
                    wf.loc[~np.isfinite(grp["log2_rec"]), "log2_rate"] = np.nan
                # Missing error bars get the enzyme median, so weights stay sane.
                med = wf.groupby("enzyme")["log2_se"].transform("median")
                wf["log2_se"] = wf["log2_se"].fillna(med).fillna(1.0)
            frames.append(wf)
        est_frames.append(estimate_genome_level(pd.concat(frames, ignore_index=True), manifest, variant))
    est = pd.concat(est_frames, ignore_index=True)
    meta_df = perwin[["sample", "arm", "depth", "true_log2_ptr", "rep"]].drop_duplicates("sample")
    est = est.merge(meta_df, on="sample", how="left")
    est.to_csv(results / "genome_estimates.tsv", sep="\t", index=False, float_format="%.10g")

    # ------------------------------------------------------------------
    # genome-level slope: estimated log2PTR ~ true b, per arm x depth x variant
    # ------------------------------------------------------------------
    rows = []
    for (arm, depth, variant), grp in est.groupby(["arm", "depth", "variant"]):
        slope, intercept, r, n = ols(grp["true_log2_ptr"].to_numpy(), grp["log2_ptr"].to_numpy())
        rows.append({"arm": arm, "depth": depth, "variant": variant,
                     "genome_slope_est_on_truth": slope, "intercept": intercept,
                     "pearson_r": r, "n_samples": n})
    genome_slopes = pd.DataFrame(rows)
    genome_slopes.to_csv(results / "genome_slopes_by_depth.tsv", sep="\t", index=False, float_format="%.10g")

    # ------------------------------------------------------------------
    # window-rate slope: log2(recovered) ~ log2(true), per arm x depth
    # ------------------------------------------------------------------
    rows = []
    pw = perwin[np.isfinite(perwin["log2_rec"]) & np.isfinite(perwin["log2_true"])]
    for (arm, depth), grp in pw.groupby(["arm", "depth"]):
        # Per-(replicate, enzyme) slopes, b = 0 cells excluded (no x-variance).
        slopes = []
        for (sample, _enzyme), sg in grp[grp["true_log2_ptr"] > 0].groupby(["sample", "enzyme"]):
            s, _, _, n = ols(sg["log2_true"].to_numpy(), sg["log2_rec"].to_numpy())
            if np.isfinite(s) and n >= 5:
                slopes.append(s)
        # Pooled slope on enzyme-and-sample demeaned values (removes the
        # per-enzyme intercept so only the gradient is regressed).
        g2 = grp[grp["true_log2_ptr"] > 0].copy()
        keys = ["sample", "enzyme"]
        g2["x"] = g2["log2_true"] - g2.groupby(keys)["log2_true"].transform("mean")
        g2["y"] = g2["log2_rec"] - g2.groupby(keys)["log2_rec"].transform("mean")
        ps, _, pr, pn = ols(g2["x"].to_numpy(), g2["y"].to_numpy())
        rows.append({"arm": arm, "depth": depth,
                     "window_rate_slope_pooled": ps, "pooled_r": pr, "n_windows": pn,
                     "window_rate_slope_mean": float(np.mean(slopes)) if slopes else np.nan,
                     "window_rate_slope_median": float(np.median(slopes)) if slopes else np.nan,
                     "n_fits": len(slopes)})
    wr_slopes = pd.DataFrame(rows)
    wr_slopes.to_csv(results / "window_rate_slopes.tsv", sep="\t", index=False, float_format="%.10g")

    # ------------------------------------------------------------------
    # (a)+(b) bias by true-rate regime, including the windows the layer drops
    # ------------------------------------------------------------------
    perwin["rate_regime"] = pd.cut(perwin["true_rate"], RATE_BINS, labels=RATE_LABELS, right=False)
    rows = []
    for (arm, depth, regime), grp in perwin.groupby(["arm", "depth", "rate_regime"], observed=True):
        finite = np.isfinite(grp["log2_ratio"])
        rows.append({"arm": arm, "depth": depth, "true_rate_regime": str(regime),
                     "n_windows": int(len(grp)),
                     "frac_dropped": float(1.0 - finite.mean()),
                     "median_log2_ratio": float(grp["log2_ratio"].median()) if finite.any() else np.nan,
                     "median_true_rate": float(grp["true_rate"].median()),
                     "median_recovered_rate": float(grp["rate"].median())})
    bias = pd.DataFrame(rows)
    bias.to_csv(results / "bias_by_regime.tsv", sep="\t", index=False, float_format="%.10g")

    # ------------------------------------------------------------------
    # (c) model selection vs bias; floor/ceiling
    # ------------------------------------------------------------------
    rows = []
    for (arm, depth), grp in perwin.groupby(["arm", "depth"]):
        finite = grp[np.isfinite(grp["log2_ratio"])]
        rows.append({"arm": arm, "depth": depth,
                     "frac_ztnb": float((grp["model"] == "ztnb").mean()),
                     "frac_empty": float((grp["model"] == "empty").mean()),
                     "mean_n_components": float(grp["n_components"].mean()),
                     "corr_components_absbias": (
                         float(np.corrcoef(finite["n_components"], np.abs(finite["log2_ratio"]))[0, 1])
                         if len(finite) > 10 else np.nan),
                     "min_true_rate": float(grp["true_rate"].min()),
                     "min_recovered_rate": float(grp["rate"].min()),
                     "p05_recovered_rate": float(np.nanpercentile(grp["rate"], 5)),
                     })
    modelsel = pd.DataFrame(rows)
    modelsel.to_csv(results / "model_selection.tsv", sep="\t", index=False, float_format="%.10g")

    # ------------------------------------------------------------------
    # summary to stdout
    # ------------------------------------------------------------------
    pd.set_option("display.width", 200)
    print("\n=== window-rate slope: log2(recovered) ~ log2(true), per depth ===")
    print(wr_slopes.to_string(index=False))
    print("\n=== genome-level slope: estimated log2PTR ~ true b, per depth x variant ===")
    print(genome_slopes.pivot_table(index=["arm", "depth"], columns="variant",
                                    values="genome_slope_est_on_truth").to_string())
    print("\n=== bias by true-rate regime (median log2 recovered/true; frac_dropped) ===")
    b0 = bias[bias["arm"].isin(["pois", "nb"])]
    print(b0.pivot_table(index=["arm", "depth", "true_rate_regime"],
                         values=["median_log2_ratio", "frac_dropped"]).to_string())
    print("\n=== model selection ===")
    print(modelsel.to_string(index=False))


if __name__ == "__main__":
    main()
