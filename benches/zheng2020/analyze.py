"""Zheng et al. 2020 E. coli benchmark: sk2bGrow vs Pilea.

A finite estimate and a deployable estimate are different outcomes. This script
therefore reports both the all-finite estimator view and the default-QC view.
"""
import os
import re
import glob

import numpy as np
import pandas as pd
from scipy import stats

# Ordered as the 2x2 of sketch x estimator, then Pilea at its shipped gates.
ARM = {'A':'anchors + V-fit (sk2bGrow)', 'E':'FracMinHash + V-fit',
       'B':'anchors + rank regression', 'C_relaxed':'Pilea (gates off)',
       'C_default':'Pilea (defaults)'}

gt = pd.read_csv('growth_rates.tsv', sep='\t').set_index('medium')

rows = []
for d in sorted(glob.glob('out/*/')):
    name = d.split('/')[1]
    m = re.match(r'^(A|B|E|C_default|C_relaxed)_(.+?)_s(\d+)_([\d.]+)x$', name)
    if not m: continue
    arm, medium, seed_id, cov = m.group(1), m.group(2), int(m.group(3)), float(m.group(4))
    try:
        df = pd.read_csv(d + 'output.tsv', sep='\t', na_values=['NA','n/a'])
    except Exception:
        rows.append(dict(arm=arm, medium=medium, cov=cov, seed=seed_id, log2ptr=np.nan)); continue
    if df.empty:
        rows.append(dict(arm=arm, medium=medium, cov=cov, seed=seed_id, log2ptr=np.nan)); continue
    r = df.iloc[0]
    rows.append(dict(arm=arm, medium=medium, cov=cov, seed=seed_id,
                     log2ptr=r.get('log2(PTR)', np.nan),
                     est_cov=r.get('coverage', np.nan),
                     passed=bool(r.get('pass_qc', True))))
res = pd.DataFrame(rows)
if rows:
    res['commit'] = os.environ.get('SK2BGROW_COMMIT', 'unrecorded')
    res['pilea_version'] = os.environ.get('PILEA_VERSION', 'unrecorded')
else:
    # Allow the QC view to be regenerated from an archived raw table without
    # pretending that the raw profile runs happened again.
    res = pd.read_csv('results_raw.tsv', sep='\t')
    if 'commit' not in res:
        res['commit'] = 'unrecorded'
    if 'pilea_version' not in res:
        res['pilea_version'] = 'unrecorded'
res['growth_rate'] = res['medium'].map(gt['growth_rate'])
res['pred_log2ptr'] = res['medium'].map(gt['pred_log2ptr'])
res.to_csv('results_raw.tsv', sep='\t', index=False, float_format='%.4f')

grow = res[res['medium'] != 'RUN_OUT'].copy()
qc_rows = []


def metrics(s):
    """Correlation/rmse/slope for a selection of finite per-condition rows."""
    if len(s) < 3 or s['log2ptr'].nunique(dropna=True) < 2:
        return len(s), np.nan, np.nan, np.nan, np.nan
    r, _ = stats.pearsonr(s['growth_rate'], s['log2ptr'])
    rho, _ = stats.spearmanr(s['growth_rate'], s['log2ptr'])
    ok = s[np.isfinite(s['pred_log2ptr'])]
    rmse = np.sqrt(np.mean((ok['log2ptr'] - ok['pred_log2ptr'])**2)) if len(ok) else np.nan
    try:
        slope = np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0]
    except (np.linalg.LinAlgError, ValueError):
        slope = np.nan
    return len(s), r, rho, rmse, slope


def emit(view, label):
    print("=" * 112)
    print(f"{label} (Pearson/Spearman vs measured lambda; RMSE vs non-independent Zheng prediction)")
    print("=" * 112)
    print(f"{'coverage':>9}  {'arm':30}{'n':>4}{'n_pass':>8}{'Pearson r':>11}"
          f"{'Spearman':>10}{'RMSE vs pred':>14}{'slope':>8}")
    for cov in sorted(view['cov'].unique()):
        for arm in ARM:
            finite = view[(view['cov'] == cov) & (view['arm'] == arm) & np.isfinite(view['log2ptr'])]
            passed = finite[finite['passed'].fillna(False).astype(bool)]
            n, r, rho, rmse, slope = metrics(finite)
            n_pass = int(passed.shape[0])
            if np.isnan(r):
                print(f"{cov:>8}x  {ARM[arm]:30}{n:>4}{n_pass:>8}{'--':>11}{'--':>10}"
                      f"{'--':>14}{'--':>8}")
            else:
                print(f"{cov:>8}x  {ARM[arm]:30}{n:>4}{n_pass:>8}{r:>11.4f}{rho:>10.4f}"
                      f"{rmse:>14.4f}{slope:>8.3f}")
            qc_rows.append(dict(arm=arm, coverage=cov, view=label, n_finite=n,
                                n_pass_qc=n_pass, pearson_r=r, spearman=rho,
                                rmse_vs_pred=rmse, slope_vs_growth=slope))
        print()


emit(grow, "all finite estimates")
emit(grow[grow['passed'].fillna(False).astype(bool)], "default QC passed")
pd.DataFrame(qc_rows).to_csv('results_qc.tsv', sep='\t', index=False, float_format='%.6g')

print("NEGATIVE CONTROL — RUN_OUT (stationary phase, expect log2(PTR) ~ 0)")
ctl = res[res['medium']=='RUN_OUT']
for arm in ARM:
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
