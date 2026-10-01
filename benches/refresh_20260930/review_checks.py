#!/usr/bin/env python3
"""Post-run review checks for the 2026-09-30 refresh.

Creates bootstrap intervals at the medium level (Zheng), replicate-level
intervals (A4), Wilson intervals and threshold sensitivity (R3). These are
diagnostics for review response, not a replacement for the primary analyzer.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

RNG = np.random.default_rng(20260930)


def metric_rows(fin: pd.DataFrame) -> tuple[float, float, float, float]:
    if len(fin) < 3 or fin["log2ptr"].nunique(dropna=True) < 2:
        return np.nan, np.nan, np.nan, np.nan
    pearson = pearsonr(fin["growth_rate"], fin["log2ptr"])[0]
    rho = spearmanr(fin["growth_rate"], fin["log2ptr"])[0]
    rmse = np.sqrt(np.mean((fin["log2ptr"] - fin["pred_log2ptr"]) ** 2))
    slope = np.polyfit(fin["growth_rate"], fin["log2ptr"], 1)[0]
    return float(pearson), float(rho), float(rmse), float(slope)


def zheng_bootstrap(raw: pd.DataFrame, n_boot: int = 10000):
    media = sorted(raw.loc[raw["medium"] != "RUN_OUT", "medium"].unique())
    out = []
    delta_rows = []
    for arm in ["A", "E", "B", "C_relaxed", "C_default"]:
        for cov in [0.5, 1.0, 2.0, 5.0, 10.0]:
            q = raw[(raw["arm"] == arm) & (raw["cov"] == cov) & (raw["medium"] != "RUN_OUT")]
            per_media = {}
            for medium in media:
                # One technical-seed estimate per medium. Bootstrap media only;
                # treating the three seed rows as independent would deflate r.
                gm = q[q["medium"] == medium].groupby("seed")[["log2ptr", "pred_log2ptr", "growth_rate"]].mean().mean()
                per_media[medium] = gm
            full = pd.DataFrame([
                dict(medium=medium, log2ptr=gm["log2ptr"],
                     pred_log2ptr=gm["pred_log2ptr"], growth_rate=gm["growth_rate"])
                for medium, gm in per_media.items()
            ])
            full = full[np.isfinite(full["log2ptr"])]
            point = metric_rows(full)
            boots = np.full((n_boot, 4), np.nan)
            for b in range(n_boot):
                chosen = RNG.choice(media, size=len(media), replace=True)
                parts = []
                for j, medium in enumerate(chosen):
                    gm = per_media[medium]
                    parts.append(dict(medium=f"{medium}_{j}", log2ptr=gm["log2ptr"],
                                      pred_log2ptr=gm["pred_log2ptr"],
                                      growth_rate=gm["growth_rate"]))
                boot = pd.DataFrame.from_records(parts)
                boot = boot[np.isfinite(boot["log2ptr"])]
                boots[b] = metric_rows(boot)
            names = ["pearson_r", "spearman_rho", "rmse_vs_pred", "slope_vs_growth"]
            for j, name in enumerate(names):
                vals = boots[:, j]
                vals = vals[np.isfinite(vals)]
                lo, hi = (np.quantile(vals, [0.025, 0.975]) if len(vals) >= 100 else (np.nan, np.nan))
                out.append(dict(arm=arm, coverage=cov, metric=name,
                                estimate=point[j], ci_low=lo, ci_high=hi,
                                n_boot_finite=len(vals), method="media bootstrap"))
            # Paired A-E deltas on the same bootstrap media.
            if arm == "A":
                for other in ["E"]:
                    qo = raw[(raw["arm"] == other) & (raw["cov"] == cov) & (raw["medium"] != "RUN_OUT")]
                    per_other = {}
                    for medium in media:
                        gm = qo[qo["medium"] == medium].groupby("seed")[["log2ptr", "pred_log2ptr", "growth_rate"]].mean().mean()
                        per_other[medium] = gm
                    a_full = full.copy()
                    o_full = pd.DataFrame([
                        dict(medium=medium, log2ptr=gm["log2ptr"],
                             pred_log2ptr=gm["pred_log2ptr"], growth_rate=gm["growth_rate"])
                        for medium, gm in per_other.items()
                    ])
                    o_full = o_full[np.isfinite(o_full["log2ptr"])]
                    ap = metric_rows(a_full)
                    op = metric_rows(o_full)
                    point_deltas = {name: ap[j] - op[j] for j, name in enumerate(names)}
                    boot_deltas = boots - np.array([metric_rows(o_full)])
                    row = dict(coverage=cov, comparison=f"A_minus_{other}")
                    for j, name in enumerate(names):
                        vals = boot_deltas[:, j]
                        vals = vals[np.isfinite(vals)]
                        lo, hi = (np.quantile(vals, [0.025, 0.975]) if len(vals) >= 100 else (np.nan, np.nan))
                        greater = np.mean(vals > 0)
                        less = np.mean(vals < 0)
                        row[f"{name}_delta"] = point_deltas[name]
                        row[f"{name}_ci_low"] = lo
                        row[f"{name}_ci_high"] = hi
                        row[f"{name}_bootstrap_p"] = float(min(1.0, 2 * min(greater, less)))
                    delta_rows.append(row)
    return pd.DataFrame(out), pd.DataFrame(delta_rows)


def zheng_qc_counts(raw: pd.DataFrame):
    rows = []
    for arm in ["A", "E", "B", "C_relaxed", "C_default"]:
        for cov in [0.5, 1.0, 2.0, 5.0, 10.0]:
            q = raw[(raw["arm"] == arm) & (raw["cov"] == cov) & (raw["medium"] != "RUN_OUT")]
            row = dict(arm=arm, coverage=cov)
            for seed in sorted(q["seed"].unique()):
                qq = q[q["seed"] == seed]
                row[f"pass_{seed}"] = int(qq[qq["passed"] == True]["medium"].nunique())
            # Strict deployment set: passing in every seed.
            wide = q.pivot_table(index="medium", columns="seed", values="passed", aggfunc="first")
            row["pass_all_seeds"] = int(wide.all(axis=1).sum()) if len(wide) else 0
            row["pass_any_seed"] = int(wide.any(axis=1).sum()) if len(wide) else 0
            rows.append(row)
    return pd.DataFrame(rows)


def a4_bootstrap(est: pd.DataFrame, n_boot: int = 10000):
    rows = []
    rng = np.random.default_rng(20260930)
    for arm in ["pois", "nb"]:
        for depth in [0.5, 1.0, 2.0, 5.0]:
            # Include b=0 so the slope cell matches the primary analyzer.
            base = est[(est["arm"] == arm) & (est["depth"] == depth) &
                       (est["variant"] == "recovered")]
            if base.empty:
                continue
            all_fit = base.copy()
            X = np.column_stack([np.ones(len(all_fit)), all_fit["true_log2_ptr"]])
            point = float(np.linalg.lstsq(X, all_fit["log2_ptr"], rcond=None)[0][1])
            boots = np.empty(n_boot)
            for b in range(n_boot):
                parts = []
                for truth, grp in base.groupby("true_log2_ptr"):
                    take = rng.choice(np.arange(len(grp)), size=len(grp), replace=True)
                    parts.append(grp.iloc[take])
                boot = pd.concat(parts, ignore_index=True)
                X = np.column_stack([np.ones(len(boot)), boot["true_log2_ptr"]])
                boots[b] = np.linalg.lstsq(X, boot["log2_ptr"], rcond=None)[0][1]
            zero = est[(est["arm"] == arm) & (est["depth"] == depth) &
                       (est["variant"] == "recovered") & (est["true_log2_ptr"] == 0)]
            zero_point = float(zero["log2_ptr"].mean())
            zero_boot = np.array([
                zero.iloc[rng.choice(np.arange(len(zero)), size=len(zero), replace=True)]["log2_ptr"].mean()
                for _ in range(n_boot)
            ])
            rows.append(dict(arm=arm, depth=depth, n_samples=len(base),
                             slope=point,
                             slope_ci_low=float(np.quantile(boots, 0.025)),
                             slope_ci_high=float(np.quantile(boots, 0.975)),
                             b0_floor=zero_point,
                             b0_ci_low=float(np.quantile(zero_boot, 0.025)),
                             b0_ci_high=float(np.quantile(zero_boot, 0.975))))
    return pd.DataFrame(rows)


def wilson(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def r3_sensitivity(stats: pd.DataFrame):
    rows = []
    for zshort in [2.0, 2.25, 2.5, 2.75, 3.0]:
        slope = ((stats["z_loc"] > 2.0) & (stats["z_short"] > zshort) &
                 ((stats["A_loc"] - stats["b_hat"]) > 0.5) &
                 stats["slope_usable"].astype(bool))
        jump = stats["flag_jump"]
        flag = slope | jump
        strong = (stats["log2ptr"] >= 1) & (stats["depth"].isin([5, 10])) & stats["ref"].str.contains("scr")
        k_strong = int((strong & flag).sum()); n_strong = int(strong.sum())
        k_false = int(((stats["log2ptr"] == 0) & flag).sum()); n_false = int((stats["log2ptr"] == 0).sum())
        slo, shi = wilson(k_strong, n_strong)
        flo, fhi = wilson(k_false, n_false)
        rows.append(dict(z_short_threshold=zshort, slope_flags=int(slope.sum()),
                         jump_flags=int(jump.sum()), total_flags=int(flag.sum()),
                         strong_5_10x_fires=k_strong, strong_5_10x_runs=n_strong,
                         strong_5_10x_fraction=k_strong / n_strong,
                         strong_ci_low=slo, strong_ci_high=shi,
                         stationary_false_fires=k_false, stationary_runs=n_false,
                         stationary_false_fraction=k_false / n_false,
                         stationary_ci_low=flo, stationary_ci_high=fhi))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args()
    raw = pd.concat([
        pd.read_csv(args.root / f"zheng/s{s}/results_raw.tsv", sep="\t").assign(seed=f"s{s}")
        for s in range(3)
    ], ignore_index=True)
    boot, deltas = zheng_bootstrap(raw)
    qc = zheng_qc_counts(raw)
    est = pd.read_csv(args.root / "a4/results/genome_estimates.tsv", sep="\t")
    a4 = a4_bootstrap(est)
    stats = pd.read_csv(args.root / "r3/r3_stats.tsv", sep="\t")
    r3 = r3_sensitivity(stats)
    out = args.root / "results"
    out.mkdir(exist_ok=True)
    boot.to_csv(out / "zheng_bootstrap_media.tsv", sep="\t", index=False, float_format="%.8g")
    deltas.to_csv(out / "zheng_arm_delta_bootstrap.tsv", sep="\t", index=False, float_format="%.8g")
    qc.to_csv(out / "zheng_qc_pass_counts.tsv", sep="\t", index=False)
    a4.to_csv(out / "a4_slope_bootstrap.tsv", sep="\t", index=False, float_format="%.8g")
    r3.to_csv(out / "r3_threshold_sensitivity.tsv", sep="\t", index=False, float_format="%.8g")
    print(f"wrote review diagnostics to {out}")


if __name__ == "__main__":
    main()
