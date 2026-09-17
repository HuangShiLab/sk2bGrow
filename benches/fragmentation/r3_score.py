#!/usr/bin/env python3
"""Score the R3 grid: WCG statistic on every (reference, log2PTR, depth) cell.

    WORK=/path/to/work python3 r3_score.py

Writes $WORK/r3_stats.tsv and prints the evaluation table: one row per
reference condition x simulated log2PTR x depth, with the fused V-fit
estimate, both WCG branches, and whether the proposed gate fires.
References suffixed R are the rotated-origin replicate family (ori at a
generic interior position rather than at a fragment cut point).
"""
import glob
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_within_contig import wcg_stat

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("WORK") or os.path.join(HERE, "..", "work", "r3")

REF_LABEL = {
    "complete": "complete (N=1)",
    "ord10": "ordered N=10",
    "ord100": "ordered N=100",
    "scr2": "scrambled N=2",
    "scr5": "scrambled N=5",
    "scr10": "scrambled N=10",
    "scr20": "scrambled N=20",
    "scr50": "scrambled N=50",
    "scr100": "scrambled N=100",
    "completeR": "complete, rot",
    "ord100R": "ordered N=100, rot",
    "scr2R": "scrambled N=2, rot",
    "scr20R": "scrambled N=20, rot",
    "scr100R": "scrambled N=100, rot",
    "completeRX": "complete, ori-int",
    "ord100RX": "ordered N=100, ori-int",
    "scr2RX": "scrambled N=2, ori-int",
    "scr20RX": "scrambled N=20, ori-int",
    "scr100RX": "scrambled N=100, ori-int",
}
REF_ORDER = list(REF_LABEL)


def ref_length(ref):
    total = 0
    subdir = "refsR" if ref.endswith(("R", "RX")) else "refs"
    with open(os.path.join(WORK, subdir, ref.removesuffix("X") + ".fna")) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def main():
    rows = []
    for d in sorted(glob.glob(os.path.join(WORK, "out", "*/"))
                    + glob.glob(os.path.join(WORK, "outX", "*/"))):
        name = os.path.basename(d.rstrip("/"))
        m = re.match(r"^(completeR?X?|ord\d+R?X?|scr\d+R?X?)_lp([\d.]+)_(\d+)x$", name)
        if not m:
            continue
        ref, lp, cov = m.group(1), float(m.group(2)), float(m.group(3))
        wf, of = os.path.join(d, "windows.rates.tsv"), os.path.join(d, "output.tsv")
        if not (os.path.exists(wf) and os.path.exists(of)):
            continue
        odf = pd.read_csv(of, sep="\t", na_values=["NA", "n/a"])
        if odf.empty:
            continue
        b = float(odf.iloc[0]["log2(PTR)"])
        se_b = float(odf.iloc[0].get("se", np.nan))
        passed = bool(odf.iloc[0].get("pass_qc", True))
        win = pd.read_csv(wf, sep="\t")
        st = wcg_stat(win, ref_length(ref), b, se_b)
        rows.append(dict(ref=ref, log2ptr=lp, depth=cov, vf_est=b, vf_se=se_b,
                         pass_qc=passed, **st))
    res = pd.DataFrame(rows)
    if res.empty:
        sys.exit("no runs found under $WORK/out")
    res.to_csv(os.path.join(WORK, "r3_stats.tsv"), sep="\t", index=False,
               float_format="%.4f")

    hdr = (f"{'reference':>20}{'truth':>6}{'depth':>6}{'V-fit':>7}{'A_loc':>7}"
           f"{'z_loc':>7}{'z_sht':>7}{'J':>7}{'z_J':>6}{'QC':>6}{'WCG':>6}")
    print(hdr)
    print("-" * len(hdr))
    for ref in REF_ORDER:
        sub = res[res["ref"] == ref]
        if sub.empty:
            continue
        for lp in sorted(sub["log2ptr"].unique()):
            for cov in sorted(sub["depth"].unique()):
                r = sub[(sub["log2ptr"] == lp) & (sub["depth"] == cov)].iloc[0]
                f = lambda v, w=7: f"{v:>{w}.3f}" if np.isfinite(v) else f"{'--':>{w}}"
                print(f"{REF_LABEL[ref]:>20}{lp:>6.1f}{cov:>5g}x"
                      f"{f(r['vf_est'])}{f(r['A_loc'])}{f(r['z_loc'])}"
                      f"{f(r['z_short'])}{f(r['J'])}{f(r['z_J'], 6)}"
                      f"{'pass' if r['pass_qc'] else 'FAIL':>6}"
                      f"{'FIRE' if r['flag'] else '':>6}")
        print()

    grow = res[res["log2ptr"] > 0]
    stat = res[res["log2ptr"] == 0]
    print("DETECTION SUMMARY (growing cells: fire rate; stationary: false-fire rate)")
    for ref in REF_ORDER:
        g, s = grow[grow["ref"] == ref], stat[stat["ref"] == ref]
        if g.empty and s.empty:
            continue
        print(f"  {REF_LABEL[ref]:>20}  growing fires {int(g['flag'].sum())}/{len(g)}"
              f"   stationary fires {int(s['flag'].sum())}/{len(s)}")


if __name__ == "__main__":
    main()
