#!/usr/bin/env python3
"""Audit what the per-enzyme and fused QC flags actually reject.

Run on the Zheng grid it reproduces the numbers quoted in
``tests/python/test_fusion_qc.py``. Run on a C5 table it says how much of
recall = 1.00 each mechanism accounts for.

    python3 benches/qc_audit.py PER_ENZYME.tsv [--fused FUSED.tsv]

``PER_ENZYME.tsv`` is the output of ``sk2bgrow.fit.fit_windows`` (or the
committed ``sk2bGrow-paper/data/per_enzyme_zheng.tsv``, whose columns are
abbreviated -- both spellings are accepted). ``--fused`` is the output of
``sk2bgrow.fusion.fuse_table``, one row per (sample, genome).
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

# The committed Zheng table abbreviates several columns.
ALIASES = {
    "log2ptr": "log2_ptr",
    "nw": "n_windows_used",
    "na": "n_anchors",
    "medium": "sample",
}


def load(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={k: v for k, v in ALIASES.items() if k in df.columns and v not in df.columns})
    if "ok" in df.columns:
        df["ok"] = df["ok"].astype(str).str.strip().str.lower().isin(("true", "1", "yes"))
    for c in ("log2_ptr", "r2", "se", "depth", "mean_rate"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def rule(title: str) -> None:
    print()
    print("=" * 74)
    print(title)
    print("=" * 74)


def audit_per_enzyme(df: pd.DataFrame) -> None:
    keys = [c for c in ("depth",) if c in df.columns] or [None]
    cell = [c for c in ("sample", "genome") if c in df.columns]

    rule("1. How many enzymes survive per (sample, genome)")
    print("   k = 1 is the case where Cochran's Q cannot be computed, so")
    print("   `consistent` is True for a test that never ran.")
    print(f"\n{'depth':>7} {'cells':>7} {'mean k':>8} {'k=0':>6} {'k=1':>6} {'fit_rate':>10}")
    for key, g in (df.groupby(keys[0]) if keys[0] else [(None, df)]):
        if not cell:
            continue
        surv = g.groupby(cell)["ok"].sum()
        att = g.groupby(cell).size()
        lab = f"{key:>7}" if key is not None else "    all"
        print(f"{lab} {len(surv):>7} {surv.mean():>8.1f} {(surv == 0).sum():>6} "
              f"{(surv == 1).sum():>6} {(surv / att).mean():>10.3f}")

    rule("2. What the per-enzyme `ok` flag rejects")
    if "ok" not in df.columns:
        print("   no `ok` column")
        return
    n, nbad = len(df), int((~df["ok"]).sum())
    print(f"   fits {n}, ok=False {nbad} ({nbad / n * 100:.1f}%)")
    if nbad and "log2_ptr" in df.columns:
        bad = df[~df["ok"]]
        nonfinite = int((~np.isfinite(bad["log2_ptr"])).sum())
        negative = int((bad["log2_ptr"] < 0).sum())
        print(f"     non-finite estimate : {nonfinite}")
        print(f"     negative estimate   : {negative}")
        print("   fit.py sets `ok = log2_ptr >= 0`. If `negative` is 0 that")
        print("   condition rejected nothing -- every rejection came from an")
        print("   earlier return path, so the flag is a sign check, not a")
        print("   quality check.")

    if "r2" in df.columns:
        rule("3. Accepted fits that are worse than a horizontal line (r2 < 0)")
        acc = df[df["ok"]]
        print(f"   overall: {int((acc['r2'] < 0).sum())}/{len(acc)} "
              f"({(acc['r2'] < 0).mean() * 100:.1f}%)")
        if keys[0]:
            print(f"\n{'depth':>7} {'r2<0':>8} {'of':>6} {'share':>8} {'mean r2':>9}")
            for key, g in acc.groupby(keys[0]):
                print(f"{key:>7} {int((g['r2'] < 0).sum()):>8} {len(g):>6} "
                      f"{(g['r2'] < 0).mean() * 100:>7.1f}% {g['r2'].mean():>9.4f}")
        print("\n   Nothing gates on r2: `fuse` never reads it, and before the")
        print("   min_r2 / n_enzymes_negative_r2 columns were added it did not")
        print("   survive fusion at all.")


def audit_fused(df: pd.DataFrame) -> None:
    rule("4. Fused table: how many estimates rest on an unrun check")
    n = len(df)
    print(f"   rows (sample x genome): {n}")
    if "consistency_checked" in df.columns:
        unchecked = int((~df["consistency_checked"].astype(bool)).sum())
        print(f"   consistency_checked=False : {unchecked} ({unchecked / n * 100:.1f}%)")
        if "consistent" in df.columns:
            both = df["consistent"].astype(bool) & df["consistency_checked"].astype(bool)
            print(f"   consistent alone          : {int(df['consistent'].astype(bool).sum())}")
            print(f"   checked AND consistent    : {int(both.sum())}   <-- the honest QC pass count")
    else:
        print("   no `consistency_checked` column -- table predates the fix;")
        print("   re-fuse to get it. Rows with n_enzymes == 1 are the affected set:")
        if "n_enzymes" in df.columns:
            k1 = int((df["n_enzymes"] == 1).sum())
            print(f"     n_enzymes == 1: {k1} ({k1 / n * 100:.1f}%)")

    if "coverage" in df.columns:
        rule("5. The coverage column, which any suspicious detector thresholds on")
        cov = pd.to_numeric(df["coverage"], errors="coerce")
        print("   quantiles:", ", ".join(
            f"p{q}={cov.quantile(q / 100):.3f}" for q in (1, 5, 25, 50, 75, 95, 99)))
        print("   NOTE this is the mean over *surviving* enzymes, unweighted, so it")
        print("   is biased upward exactly where coverage is low. Compare against")
        print("   enzyme_fit_rate before trusting a threshold on it.")
        if "enzyme_fit_rate" in df.columns:
            fr = pd.to_numeric(df["enzyme_fit_rate"], errors="coerce")
            print(f"   enzyme_fit_rate: median {fr.median():.3f}, "
                  f"share below 0.5 = {(fr < 0.5).mean() * 100:.1f}%")
            lowfit = fr < 0.5
            if lowfit.any():
                print(f"   rows where most of the panel saw nothing yet coverage reads "
                      f"{cov[lowfit].median():.3f} (median): {int(lowfit.sum())}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("per_enzyme")
    ap.add_argument("--fused")
    a = ap.parse_args()

    df = load(a.per_enzyme)
    print(f"per-enzyme table: {df.shape[0]} rows, "
          f"{df['enzyme'].nunique() if 'enzyme' in df else '?'} enzymes")
    audit_per_enzyme(df)
    if a.fused:
        audit_fused(load(a.fused))
    else:
        print("\n(no --fused given; sections 4-5 skipped)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
