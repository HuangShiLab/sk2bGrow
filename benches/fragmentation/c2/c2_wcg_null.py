#!/usr/bin/env python3
"""C2 Task 1: WCG null-margin validation on real data (R3_QC.md caveat).

Every C1 / C1b cell is a NULL for the WCG gate (correct coordinate), so the
fire rate with the tuned thresholds should be ~0 and the measured margins
(max J vs the 0.03 fire line, max z_J vs 3.5, max A_loc - b_hat on growing
cells vs 0.5) are the numbers the R3 report said must be re-measured at
scale on real reads.

Scored offline from the windows.rates.tsv / output.tsv the C1/C1b profile
runs already wrote -- no estimator re-runs. Statistic: qc_within_contig.py,
imported verbatim (same directory).

Inputs (all under $BASE):
  bench/C1/counts/<cell>/{windows.rates.tsv,output.tsv}   315 cells, K-12
  bench/C1/run_medium.tsv        run -> medium (RUN_OUT = stationary)
  bench/C1/growth_rates.tsv      medium -> measured lambda, pred_log2ptr
  bench/C1b/counts/<cell>/{windows.rates.tsv,output.tsv}  100 cells, 4 spp
  bench/C1b/run_species.tsv      run -> species, ref, glen
  bench/C1b/run_condition.tsv    run -> concentration, mu_per_h, excluded

Output:
  bench/C2/c2_wcg_null.tsv  one row per cell
  stdout: distributions per depth band, fire rates, margins, firing cells
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
C2 = os.path.join(BASE, "bench/C2")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_within_contig import wcg_stat  # noqa: E402

K12_LEN = 4_641_652
DEPTH_BANDS = [(0.25, 1.0, "0.25-1x"), (2.0, 5.0, "2-5x"), (10.0, 20.0, "10-20x")]
WINDOWS_COLS = ["enzyme", "contig_id", "start", "end", "log2_rate", "log2_se"]


def band(depth):
    for lo, hi, name in DEPTH_BANDS:
        if lo <= depth <= hi:
            return name
    return "other"


def read_tsv(path, **kw):
    return pd.read_csv(path, sep="\t", **kw)


def c1_meta():
    """run -> (medium, condition, pred_log2ptr, lambda)."""
    rm = read_tsv(f"{BASE}/bench/C1/run_medium.tsv", header=None,
                  names=["run", "medium"])
    gr = read_tsv(f"{BASE}/bench/C1/growth_rates.tsv")
    gr = gr.set_index("medium")
    out = {}
    for run, medium in zip(rm["run"], rm["medium"]):
        cond = "stationary" if medium == "RUN_OUT" else "growing"
        pred = lam = np.nan
        if medium in gr.index:
            pred = gr.loc[medium, "pred_log2ptr"]
            lam = gr.loc[medium, "growth_rate"]
        out[run] = dict(medium=medium, condition=cond, truth_pred=pred,
                        growth_rate=lam, species="Escherichia_coli_K12",
                        glen=K12_LEN)
    return out


def c1b_meta():
    """run -> (species, concentration, mu_per_h, excluded, glen)."""
    rs = read_tsv(f"{BASE}/bench/C1b/run_species.tsv")
    rc = read_tsv(f"{BASE}/bench/C1b/run_condition.tsv")
    cond = rc.set_index("run")
    out = {}
    for _, r in rs.iterrows():
        run = r["run"]
        c = cond.loc[run] if run in cond.index else {}
        out[run] = dict(medium=str(c.get("conc", "")),
                        condition="growing",
                        truth_pred=np.nan,
                        growth_rate=c.get("mu_per_h", np.nan),
                        species=r["ref"], glen=int(r["glen"]),
                        excluded=str(c.get("excluded_negative", "")).lower() == "true")
    return out


def score_cell(counts_dir, cell, meta, dataset):
    wpath = os.path.join(counts_dir, cell, "windows.rates.tsv")
    opath = os.path.join(counts_dir, cell, "output.tsv")
    if not (os.path.isfile(wpath) and os.path.isfile(opath)):
        return None
    run, dtag = cell.split(".", 1)  # SRR.... + "0.25x" (run has no dot)
    depth = float(dtag.rstrip("x"))
    m = meta.get(run)
    if m is None:
        return None
    out = read_tsv(opath).iloc[0]
    b_hat = float(out["log2(PTR)"]) if pd.notna(out["log2(PTR)"]) else np.nan
    se_b = float(out["se"]) if pd.notna(out["se"]) else np.nan
    windows = read_tsv(wpath, usecols=WINDOWS_COLS)
    st = wcg_stat(windows, m["glen"], b_hat, se_b)
    row = dict(cell=cell, dataset=dataset, run=run, depth=depth,
               depth_band=band(depth), species=m["species"],
               condition=m["condition"], medium=m["medium"],
               growth_rate=m["growth_rate"], truth_pred_log2ptr=m["truth_pred"],
               excluded_negative=m.get("excluded", False),
               log2ptr=b_hat, se=se_b,
               coverage=float(out["coverage"]), pass_qc=bool(out["pass_qc"]),
               qc_reason=str(out["qc_reason"]) if pd.notna(out["qc_reason"]) else "",
               k_contigs=st["k_contigs"], A_loc=st["A_loc"], se_A=st["se_A"],
               slope_usable=st["slope_usable"], z_loc=st["z_loc"],
               z_short=st["z_short"],
               A_minus_b=(st["A_loc"] - b_hat
                          if np.isfinite(st["A_loc"]) and np.isfinite(b_hat)
                          else np.nan),
               J=st["J"], se_J=st["se_J"], n_boundary=st["n_boundary"],
               z_J=st["z_J"], flag_slope=st["flag_slope"],
               flag_jump=st["flag_jump"], flag=st["flag"])
    return row


def summarise(df):
    pd.set_option("display.width", 200)
    print(f"\n== cells scored: {len(df)} "
          f"(C1 {(df.dataset == 'C1').sum()}, C1b {(df.dataset == 'C1b').sum()})")
    cols = ["J", "z_J", "A_loc", "z_loc", "z_short", "A_minus_b"]
    for (ds, bd, cond), g in df.groupby(["dataset", "depth_band", "condition"]):
        print(f"\n-- {ds} | {bd} | {cond} | n={len(g)}")
        q = g[cols].describe(percentiles=[0.5, 0.9, 0.99]).T
        print(q[["count", "mean", "50%", "90%", "99%", "max"]].to_string(
            float_format=lambda x: f"{x:.4f}"))
        fs, fj = int(g["flag_slope"].sum()), int(g["flag_jump"].sum())
        print(f"   fires: slope {fs}/{len(g)}  jump {fj}/{len(g)}  "
              f"either {int(g['flag'].sum())}/{len(g)}")

    print("\n== margins against the fire lines (null side, all cells)")
    print(f"   max J            = {df['J'].max():.4f}   (fire at J > 0.03)")
    print(f"   max z_J          = {df['z_J'].max():.2f}   (fire at z_J > 3.5)")
    grow = df[df["condition"] == "growing"]
    print(f"   max A_loc-b_hat  = {grow['A_minus_b'].max():.4f} on growing cells "
          f"(fire at > 0.5)")
    stat = df[df["condition"] == "stationary"]
    print(f"   stationary cells: n={len(stat)}, max A_loc = {stat['A_loc'].max():.4f}, "
          f"max z_loc = {stat['z_loc'].max():.2f}, max J = {stat['J'].max():.4f}")

    print("\n== slope branch at 0.25-1x (the A4 window-SE question): "
          "is 'slope branch only at >= 2x' supported?")
    c1 = df[df["dataset"] == "C1"]
    for bd in ["0.25-1x", "2-5x", "10-20x"]:
        g = c1[(c1["depth_band"] == bd) & (c1["condition"] == "growing")]
        u = g[g["slope_usable"]]
        print(f"   C1 growing {bd}: n={len(g)}, slope-usable {len(u)}, "
              f"A_loc median {g['A_loc'].median():.3f} max {g['A_loc'].max():.3f}, "
              f"A_minus_b median {g['A_minus_b'].median():.3f} "
              f"max {g['A_minus_b'].max():.3f}, "
              f"z_short>2: {(g['z_short'] > 2).sum()}, "
              f"flag_slope fires: {g['flag_slope'].sum()}")

    fires = df[df["flag"]]
    print(f"\n== firing cells: {len(fires)}")
    for _, r in fires.iterrows():
        print(f"   {r['cell']}  {r['dataset']} {r['condition']} "
              f"A_loc={r['A_loc']:.3f} b={r['log2ptr']:.3f} "
              f"z_loc={r['z_loc']:.2f} z_short={r['z_short']:.2f} "
              f"J={r['J']:.4f} z_J={r['z_J']:.2f} | qc: {r['qc_reason']}")


def main():
    rows = []
    c1m, c1bm = c1_meta(), c1b_meta()
    for dataset, cdir, meta in [("C1", f"{BASE}/bench/C1/counts", c1m),
                                ("C1b", f"{BASE}/bench/C1b/counts", c1bm)]:
        cells = sorted(os.listdir(cdir))
        for i, cell in enumerate(cells):
            row = score_cell(cdir, cell, meta, dataset)
            if row is not None:
                rows.append(row)
            if (i + 1) % 50 == 0:
                print(f"  [{dataset}] {i + 1}/{len(cells)}", flush=True)
    df = pd.DataFrame(rows)
    out = os.path.join(C2, "c2_wcg_null.tsv")
    df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {out}: {len(df)} rows")
    summarise(df)


if __name__ == "__main__":
    main()
