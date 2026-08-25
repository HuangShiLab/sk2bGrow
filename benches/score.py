#!/usr/bin/env python3
"""Score the multi-strain grid: accuracy and computational cost.

Accuracy uses Pilea's own metric (L2 distance between estimated and true PTR
vectors) but reports **recall alongside it**, because L2 over only the genomes a
tool chose to report rewards a tool for staying silent. A tool that reports 2 of
16 strains accurately is not better than one that reports 16 with small errors.
"""
import re, glob
from pathlib import Path
import numpy as np, pandas as pd

ARMS = {'A': 'sk2bGrow', 'C': 'Pilea (defaults)', 'D': 'Pilea (gates off)'}


def read_out(d):
    f = Path(d) / 'output.tsv'
    if not f.exists():
        return {}
    try:
        df = pd.read_csv(f, sep='\t', na_values=['NA', 'n/a'])
    except Exception:
        return {}
    if df.empty or 'log2(PTR)' not in df.columns:
        return {}
    return {r['genome']: r['log2(PTR)'] for _, r in df.iterrows()
            if pd.notna(r.get('log2(PTR)'))}


def peak_rss_and_time(f):
    try:
        t = Path(f).read_text()
    except Exception:
        return np.nan, np.nan
    rss = re.search(r'(\d+)\s+maximum resident set size', t)
    real = re.search(r'([\d.]+)\s+real', t)
    return (int(rss.group(1)) / 1e6 if rss else np.nan,
            float(real.group(1)) if real else np.nan)


rows = []
for truth_f in sorted(glob.glob('res/*.truth')):
    tag = Path(truth_f).stem
    m = re.match(r's(\d+)_c([\d.]+)_r(\d+)', tag)
    if not m:
        continue
    ns, cov, rep = int(m.group(1)), float(m.group(2)), int(m.group(3))
    truth = pd.read_csv(truth_f, sep='\t').set_index('genome')['true_log2ptr'].to_dict()
    for arm, label in ARMS.items():
        est = read_out(f'res/{arm}_{tag}')
        rss, secs = peak_rss_and_time(f'res/{arm}_{tag}.time')
        present = set(truth)
        got = {g: v for g, v in est.items() if g in present}
        spurious = len([g for g in est if g not in present])   # false positives
        if got:
            e = np.array([got[g] - truth[g] for g in got])
            l2 = float(np.sqrt((e ** 2).sum()))
            rmse = float(np.sqrt((e ** 2).mean()))
            bias = float(e.mean())
        else:
            l2 = rmse = bias = np.nan
        rows.append(dict(arm=label, n_strains=ns, coverage=cov, rep=rep,
                         n_true=len(truth), n_reported=len(got), spurious=spurious,
                         recall=len(got) / len(truth), l2=l2, rmse=rmse, bias=bias,
                         seconds=secs, peak_rss_mb=rss))

df = pd.DataFrame(rows)
df.to_csv('sim_results.tsv', sep='\t', index=False, float_format='%.4f')

print('=' * 92)
print('MULTI-STRAIN SIMULATION  (Pilea Fig-3 design, laptop scale)')
print('=' * 92)
print(f"{'strains':>8}{'cov':>6}  {'arm':20}{'recall':>8}{'RMSE':>8}{'L2':>8}{'bias':>8}{'sec':>8}{'RSS MB':>9}")
for ns in sorted(df.n_strains.unique()):
    for cov in sorted(df.coverage.unique()):
        for label in ARMS.values():
            s = df[(df.n_strains == ns) & (df.coverage == cov) & (df.arm == label)]
            if s.empty:
                continue
            print(f"{ns:>8}{cov:>6g}  {label:20}"
                  f"{s.recall.mean():>8.2f}{s.rmse.mean():>8.3f}{s.l2.mean():>8.3f}"
                  f"{s.bias.mean():>8.3f}{s.seconds.mean():>8.1f}{s.peak_rss_mb.mean():>9.0f}")
        print()

print('\nAGGREGATE  (mean over the whole grid)')
agg = df.groupby('arm').agg(recall=('recall', 'mean'), rmse=('rmse', 'mean'),
                            bias=('bias', 'mean'), spurious=('spurious', 'mean'),
                            sec=('seconds', 'mean'), rss=('peak_rss_mb', 'mean'))
print(agg.to_string(float_format=lambda v: f'{v:.3f}'))
