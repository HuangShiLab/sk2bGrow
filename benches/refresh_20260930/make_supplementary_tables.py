#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
REF=Path(__file__).resolve().parent
OUT=REF.parents[1]/"docs/paper/supplementary_tables.md"
items=[
 ("S1. Zheng all-seed all-finite and default-QC metrics","results/zheng_all_seeds_qc.tsv"),
 ("S2. Zheng medium bootstrap intervals","results/zheng_bootstrap_media.tsv"),
 ("S3. Zheng paired A-minus-E bootstrap deltas","results/zheng_arm_delta_bootstrap.tsv"),
 ("S4. Zheng QC pass counts by seed","results/zheng_qc_pass_counts.tsv"),
 ("S5. A4 recovered-slope and b=0 bootstrap intervals","results/a4_slope_bootstrap.tsv"),
 ("S6. R3 threshold sensitivity","results/r3_threshold_sensitivity.tsv"),
 ("S7. Mixed-strain cell-level results","../mixedstrain_20261001/results/sim_cells.tsv"),
 ("S8. Mixed-strain summary","../mixedstrain_20261001/results/sim_summary.tsv"),
]
lines=["# Supplementary tables",""]
for title, rel in items:
    path=(REF/rel).resolve()
    lines += [f"## {title}", "", f"Source: `{Path(rel).as_posix()}`.", ""]
    d=pd.read_csv(path,sep="\t")
    lines.append(d.to_markdown(index=False))
    lines.append("")
OUT.write_text("\n".join(lines))
print(OUT)
