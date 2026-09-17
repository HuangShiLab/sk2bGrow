#!/usr/bin/env python3
"""C2 scorer: accuracy + WCG per condition per depth (Tasks 2.1, 2.2, 2.4).

Reads counts/<condition>/<cell>/{output.tsv,windows.rates.tsv} for every
condition in refs.tsv, joins the Zheng ground truth (growth_rates.tsv via
runs.tsv), and reports per condition per depth:

  n, recall (fraction with a finite estimate; QC-pass fraction beside it),
  Pearson r and OLS slope of estimate on measured lambda (table2-style),
  RMSE and bias vs pred_log2ptr (a non-independent target -- the source
  paper's own marker-frequency analysis of the same reads; label kept),
  RUN_OUT stationary estimate,
  WCG fire rate (both branches), median J / A_loc.

Writes c2_cells.tsv (one row per cell x condition, incl. all WCG fields) and
c2_summary.tsv; prints the summary. Slope is labelled: OLS of estimate on
measured lambda.
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
C2 = os.path.join(BASE, "bench/C2")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_within_contig import wcg_stat  # noqa: E402

WINDOWS_COLS = ["enzyme", "contig_id", "start", "end", "log2_rate", "log2_se"]
K12_LEN = 4_641_652


def ref_len(fasta):
    n = 0
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def main():
    refs = pd.read_csv(f"{C2}/refs.tsv", sep="\t")
    runs = pd.read_csv(f"{C2}/runs.tsv", sep="\t")
    gr = pd.read_csv(f"{BASE}/bench/C1/growth_rates.tsv", sep="\t")
    gr = gr.set_index("medium")
    glens = {c: ref_len(f) for c, f in zip(refs["condition"], refs["fasta"])}

    rows = []
    for cond in refs["condition"]:
        cdir = f"{C2}/counts/{cond}"
        if not os.path.isdir(cdir):
            continue
        for cell in sorted(os.listdir(cdir)):
            opath = os.path.join(cdir, cell, "output.tsv")
            wpath = os.path.join(cdir, cell, "windows.rates.tsv")
            if not os.path.isfile(opath):
                continue
            run, dtag = cell.split(".", 1)  # SRR.... + "0.25x"
            depth = float(dtag.rstrip("x"))
            med = runs.loc[runs["run"] == run, "medium"]
            medium = med.iloc[0] if len(med) else ""
            lam = gr.loc[medium, "growth_rate"] if medium in gr.index else np.nan
            pred = (gr.loc[medium, "pred_log2ptr"]
                    if medium in gr.index else np.nan)
            out = pd.read_csv(opath, sep="\t").iloc[0]
            b_hat = (float(out["log2(PTR)"])
                     if pd.notna(out["log2(PTR)"]) else np.nan)
            se_b = float(out["se"]) if pd.notna(out["se"]) else np.nan
            st = (wcg_stat(pd.read_csv(wpath, sep="\t", usecols=WINDOWS_COLS),
                           glens[cond], b_hat, se_b)
                  if os.path.isfile(wpath) else {})
            rows.append(dict(
                condition=cond, cell=cell, run=run, depth=depth, medium=medium,
                growth_condition=("stationary" if medium == "RUN_OUT"
                                  else "growing"),
                lam=lam, pred_log2ptr=pred, log2ptr=b_hat, se=se_b,
                pass_qc=bool(out["pass_qc"]),
                coverage=float(out["coverage"]),
                k_contigs=st.get("k_contigs", np.nan),
                A_loc=st.get("A_loc", np.nan), z_loc=st.get("z_loc", np.nan),
                z_short=st.get("z_short", np.nan),
                J=st.get("J", np.nan), z_J=st.get("z_J", np.nan),
                flag_slope=bool(st.get("flag_slope", False)),
                flag_jump=bool(st.get("flag_jump", False)),
                wcg_flag=bool(st.get("flag", False))))
    df = pd.DataFrame(rows)
    df.to_csv(f"{C2}/c2_cells.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"wrote {C2}/c2_cells.tsv: {len(df)} rows")

    summ = []
    for (cond, depth), g in df.groupby(["condition", "depth"]):
        grow = g[g["growth_condition"] == "growing"].dropna(subset=["log2ptr"])
        stat = g[g["growth_condition"] == "stationary"]
        r = slope = rmse = bias = np.nan
        ok = grow.dropna(subset=["lam"])
        if len(ok) >= 3 and ok["log2ptr"].std() > 0:
            r = float(np.corrcoef(ok["log2ptr"], ok["lam"])[0, 1])
            slope = float(np.polyfit(ok["lam"], ok["log2ptr"], 1)[0])
        okp = grow.dropna(subset=["pred_log2ptr"])
        if len(okp) >= 1:
            err = okp["log2ptr"] - okp["pred_log2ptr"]
            rmse = float(np.sqrt((err ** 2).mean()))
            bias = float(err.mean())
        summ.append(dict(
            condition=cond, depth=depth, n=len(g),
            n_growing_est=len(grow),
            recall_qc=float(grow["pass_qc"].mean()) if len(grow) else np.nan,
            r=r, slope_est_on_lambda=slope, rmse_vs_pred=rmse, bias_vs_pred=bias,
            runout_est=float(stat["log2ptr"].median())
            if len(stat) else np.nan,
            wcg_fire=float(g["wcg_flag"].mean()),
            wcg_fire_growing=float(grow["wcg_flag"].mean()) if len(grow) else np.nan,
            J_med=float(g["J"].median()), A_loc_med=float(g["A_loc"].median())))
    s = pd.DataFrame(summ).sort_values(["condition", "depth"])
    s.to_csv(f"{C2}/c2_summary.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"wrote {C2}/c2_summary.tsv")
    pd.set_option("display.width", 250)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
