#!/usr/bin/env python3
"""Score the fragmentation benchmark: complete vs 100 contigs vs scaffolded.

    WORK=/path/to/work python3 analyze.py
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("WORK") or sys.exit("set WORK")
GT = os.path.join(HERE, "..", "zheng2020", "growth_rates.tsv")

COND = {
    "complete": "complete chromosome",
    "frag": "100 contigs",
    "scafSelf": "scaffolded vs itself",
    "scafRel": "scaffolded vs a relative",
    "pileaFrag": "Pilea on 100 contigs",
}

rows = []
for d in sorted(glob.glob(os.path.join(WORK, "out", "*/"))):
    name = os.path.basename(d.rstrip("/"))
    m = re.match(r"^(complete|frag|scafSelf|scafRel|pileaFrag)_(.+?)_([\d.]+)x$", name)
    if not m:
        continue
    cond, medium, cov = m.group(1), m.group(2), float(m.group(3))
    # Pilea writes the same filename and the same log2(PTR) column, but has no
    # pass_qc of its own -- the gates-off run reports everything it fits.
    f = os.path.join(d, "output.tsv")
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, sep="\t", na_values=["NA", "n/a"])
    if df.empty:
        continue
    val = df.iloc[0]["log2(PTR)"]
    qc = bool(df.iloc[0].get("pass_qc", True))
    rows.append(dict(cond=cond, medium=medium, cov=cov, log2ptr=val, passed=qc))

res = pd.DataFrame(rows)
if res.empty:
    sys.exit("no results under $WORK/out")
gt = pd.read_csv(GT, sep="\t").set_index("medium")
res["growth_rate"] = res["medium"].map(gt["growth_rate"])
res["pred_log2ptr"] = res["medium"].map(gt["pred_log2ptr"])
res.to_csv(os.path.join(WORK, "results_raw.tsv"), sep="\t", index=False,
           float_format="%.4f")

grow = res[res["medium"] != "RUN_OUT"]
print("=" * 92)
print("REFERENCE FRAGMENTATION  (Zheng E. coli, 16 media, same reads throughout)")
print("=" * 92)
print(f"{'coverage':>9}  {'reference':26}{'n':>4}{'Pearson r':>11}{'RMSE':>9}"
      f"{'bias':>9}{'slope':>8}{'QC pass':>9}")
for cov in sorted(grow["cov"].unique()):
    for cond, lab in COND.items():
        s = grow[(grow["cov"] == cov) & (grow["cond"] == cond)].dropna(subset=["log2ptr"])
        if s.empty:
            continue
        if len(s) < 3 or s["log2ptr"].nunique() < 2:
            print(f"{cov:>8}x  {lab:26}{len(s):>4}{'--':>11}{'--':>9}{'--':>9}{'--':>8}")
            continue
        ok = s.dropna(subset=["pred_log2ptr"])
        e = ok["log2ptr"] - ok["pred_log2ptr"]
        print(f"{cov:>8}x  {lab:26}{len(s):>4}"
              f"{stats.pearsonr(s['growth_rate'], s['log2ptr'])[0]:>11.4f}"
              f"{np.sqrt((e ** 2).mean()):>9.3f}{e.mean():>9.3f}"
              f"{np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0]:>8.3f}"
              f"{100 * s['passed'].mean():>8.0f}%")
    print()

print("PER-CONDITION at 10x (predicted log2PTR in the first column)")
top = grow["cov"].max()
piv = grow[grow["cov"] == top].pivot_table(index="medium", columns="cond", values="log2ptr")
piv.insert(0, "predicted", gt["pred_log2ptr"])
print(piv.sort_values("predicted", ascending=False)
      .to_string(float_format=lambda v: f"{v:.3f}"))
