#!/usr/bin/env python3
"""Score current sk2bGrow mixed-strain cells without absent external tools."""
from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent

def estimate(path):
    f = path / "output.tsv"
    if not f.exists(): return {}
    df = pd.read_csv(f, sep="\t", na_values=["NA", "n/a"])
    return {r["genome"]: r["log2(PTR)"] for _, r in df.iterrows()
            if pd.notna(r.get("log2(PTR)"))}

rows=[]
for truth_file in sorted((ROOT/"res").glob("*.truth")):
    tag=truth_file.stem
    m=re.fullmatch(r"s(\d+)_c([\d.]+)_r(\d+)",tag)
    if not m: continue
    truth=pd.read_csv(truth_file,sep="\t").set_index("genome")["true_log2ptr"].to_dict()
    est=estimate(ROOT/"res"/f"A_{tag}")
    got={g:v for g,v in est.items() if g in truth}
    err=np.array([got[g]-truth[g] for g in got],dtype=float) if got else np.array([])
    rows.append(dict(n_strains=int(m.group(1)),coverage=float(m.group(2)),rep=int(m.group(3)),
        n_true=len(truth),n_reported=len(got),recall=len(got)/len(truth),
        rmse=float(np.sqrt(np.mean(err**2))) if len(err) else np.nan,
        bias=float(err.mean()) if len(err) else np.nan,
        l2=float(np.sqrt(np.sum(err**2))) if len(err) else np.nan,
        spurious=sum(g not in truth for g in est)))
df=pd.DataFrame(rows).sort_values(["n_strains","coverage","rep"])
df.to_csv(ROOT/"results/sim_cells.tsv",sep="\t",index=False,float_format="%.6g")
summary=df.groupby(["n_strains","coverage"],as_index=False).agg(recall=("recall","mean"),rmse=("rmse","mean"),bias=("bias","mean"),n=("rep","count"))
summary["arm"]="sk2bGrow"
summary.to_csv(ROOT/"results/sim_summary.tsv",sep="\t",index=False,float_format="%.6g")
print(summary.to_string(index=False))
