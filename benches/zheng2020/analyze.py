"""Zheng et al. 2020 E. coli benchmark: sk2bGrow vs Pilea."""
import re, sys, glob
import numpy as np, pandas as pd
from scipy import stats

ARM = {'A':'sk2bGrow', 'B':'sk2bGrow (Pilea-parity est.)',
       'C_default':'Pilea (defaults)', 'C_relaxed':'Pilea (gates off)'}

gt = pd.read_csv('growth_rates.tsv', sep='\t').set_index('medium')

rows = []
for d in sorted(glob.glob('out/*/')):
    name = d.split('/')[1]
    m = re.match(r'^(A|B|C_default|C_relaxed)_(.+?)_([\d.]+)x$', name)
    if not m: continue
    arm, medium, cov = m.group(1), m.group(2), float(m.group(3))
    try:
        df = pd.read_csv(d + 'output.tsv', sep='\t', na_values=['NA','n/a'])
    except Exception:
        rows.append(dict(arm=arm, medium=medium, cov=cov, log2ptr=np.nan)); continue
    if df.empty:
        rows.append(dict(arm=arm, medium=medium, cov=cov, log2ptr=np.nan)); continue
    r = df.iloc[0]
    rows.append(dict(arm=arm, medium=medium, cov=cov,
                     log2ptr=r.get('log2(PTR)', np.nan),
                     est_cov=r.get('coverage', np.nan),
                     passed=bool(r.get('pass_qc', True))))
res = pd.DataFrame(rows)
res['growth_rate'] = res['medium'].map(gt['growth_rate'])
res['pred_log2ptr'] = res['medium'].map(gt['pred_log2ptr'])
res.to_csv('results_raw.tsv', sep='\t', index=False, float_format='%.4f')

grow = res[res['medium'] != 'RUN_OUT'].copy()
print("=" * 96)
print("BENCHMARK: log2(PTR) vs measured growth rate  (Zheng et al. 2020, E. coli K-12 MG1655)")
print("=" * 96)
print(f"{'coverage':>9}  {'arm':30}{'n':>4}{'Pearson r':>11}{'Spearman':>10}"
      f"{'RMSE vs pred':>14}{'slope':>8}")
for cov in sorted(grow['cov'].unique()):
    for arm in ['A','B','C_default','C_relaxed']:
        s = grow[(grow['cov']==cov) & (grow['arm']==arm) & np.isfinite(grow['log2ptr'])]
        if len(s) < 3:
            print(f"{cov:>8}x  {ARM[arm]:30}{len(s):>4}{'--':>11}{'--':>10}{'--':>14}{'--':>8}")
            continue
        r,_ = stats.pearsonr(s['growth_rate'], s['log2ptr'])
        rho,_ = stats.spearmanr(s['growth_rate'], s['log2ptr'])
        ok = s[np.isfinite(s['pred_log2ptr'])]
        rmse = np.sqrt(np.mean((ok['log2ptr']-ok['pred_log2ptr'])**2)) if len(ok) else np.nan
        slope = np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0]
        print(f"{cov:>8}x  {ARM[arm]:30}{len(s):>4}{r:>11.4f}{rho:>10.4f}{rmse:>14.4f}{slope:>8.3f}")
    print()

print("NEGATIVE CONTROL — RUN_OUT (stationary phase, expect log2(PTR) ~ 0)")
ctl = res[res['medium']=='RUN_OUT']
for arm in ['A','B','C_default','C_relaxed']:
    s = ctl[ctl['arm']==arm].sort_values('cov')
    if s.empty: continue
    vals = "  ".join(f"{c:g}x={v:.3f}" if np.isfinite(v) else f"{c:g}x=--"
                     for c,v in zip(s['cov'], s['log2ptr']))
    print(f"  {ARM[arm]:30} {vals}")

print("\nPER-CONDITION at the highest coverage")
top = grow['cov'].max()
piv = grow[grow['cov']==top].pivot_table(index='medium', columns='arm', values='log2ptr')
piv.insert(0,'growth_rate', gt['growth_rate']); piv.insert(1,'Zheng pred', gt['pred_log2ptr'])
print(piv.sort_values('growth_rate', ascending=False).to_string(float_format=lambda v:f"{v:.3f}"))
