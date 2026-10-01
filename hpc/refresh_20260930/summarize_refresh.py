#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def read_tsv(path):
    return pd.read_csv(path, sep='\t') if path.exists() else pd.DataFrame()

ap=argparse.ArgumentParser()
ap.add_argument('--zheng-root',required=True); ap.add_argument('--a4',required=True); ap.add_argument('--r3',required=True); ap.add_argument('--out',required=True)
a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
frames=[]
for p in sorted(Path(a.zheng_root).glob('s*/results_qc.tsv')):
    df=read_tsv(p); df['seed']=p.parent.name; frames.append(df)
zheng=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
if not zheng.empty:
    zheng.to_csv(out/'zheng_all_seeds_qc.tsv',sep='\t',index=False,float_format='%.8g')
    # Aggregate all-finite metrics across seeds by arm/coverage; n is media x finite seeds.
    grow=zheng[(zheng['view']=='all finite estimates') & np.isfinite(zheng['pearson_r'])]
    grow.groupby(['arm','coverage']).size().to_csv(out/'zheng_cell_counts.tsv',sep='\t')
a4p=Path(a.a4)
for n in ['window_rate_slopes.tsv','genome_slopes_by_depth.tsv','bias_by_regime.tsv','decomposition.tsv','model_selection.tsv']:
    df=read_tsv(a4p/n)
    if not df.empty: df.to_csv(out/f'a4_{n}',sep='\t',index=False,float_format='%.8g')
r3p=Path(a.r3)/'r3_stats.tsv'; r3=read_tsv(r3p)
if not r3.empty:
    r3.to_csv(out/'r3_stats.tsv',sep='\t',index=False,float_format='%.8g')
    r3.to_csv(out/'r3_summary.tsv',sep='\t',index=False,float_format='%.8g')
print('zheng rows',len(zheng),'a4 files',len(list(a4p.glob('*.tsv'))),'r3 rows',len(r3))
