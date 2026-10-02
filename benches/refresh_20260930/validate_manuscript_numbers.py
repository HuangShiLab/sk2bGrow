#!/usr/bin/env python3
"""Validate headline numbers in docs/paper/manuscript.md against result TSVs."""
from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
REF = Path(__file__).resolve().parent
PAPER = ROOT / "docs/paper/manuscript.md"
RESULTS_DOC = ROOT / "docs/PAPER_RESULTS.md"

def zheng_point(raw, arm, cov):
    q = raw[(raw["arm"] == arm) & (raw["cov"] == cov) & (raw["medium"] != "RUN_OUT")]
    if q.empty:
        return {}
    m = q.groupby("medium")[["log2ptr", "pred_log2ptr", "growth_rate"]].mean().reset_index()
    fin = m[np.isfinite(m["log2ptr"])]
    out = dict(n=len(fin))
    if len(fin) >= 3 and fin.log2ptr.nunique() > 1:
        out["r"] = pearsonr(fin.growth_rate, fin.log2ptr)[0]
        out["rho"] = spearmanr(fin.growth_rate, fin.log2ptr)[0]
        out["rmse"] = np.sqrt(np.mean((fin.log2ptr-fin.pred_log2ptr)**2))
        out["slope"] = np.polyfit(fin.growth_rate, fin.log2ptr, 1)[0]
    qc = q.pivot_table(index="medium", columns="seed", values="passed", aggfunc="first")
    out["qc_mean"] = float(qc.mean().mean()) * 16
    out["qc_sum_seeds"] = int(qc.sum().sum())
    out["qc_all_seeds"] = int(qc.all(axis=1).sum())
    return out

def parse_paper_zheng():
    text = RESULTS_DOC.read_text()
    rows = {}
    pattern = re.compile(r"^\| (A|E|B) \| ([0-9.]+)× \| (\d+) \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \|$", re.M)
    for arm, depth, n, qc, r, rho, rmse, slope in pattern.findall(text):
        rows[(arm, float(depth))] = dict(n=int(n), qc_mean=float(qc), r=float(r),
                                         rho=float(rho), rmse=float(rmse), slope=float(slope))
    return rows

def close(actual, paper, tol=0.0015):
    return abs(float(actual)-float(paper)) <= tol

def main():
    raw = pd.concat([
        pd.read_csv(REF / f"zheng/s{s}/results_raw.tsv", sep="\t").assign(seed=f"s{s}")
        for s in range(3)
    ], ignore_index=True)
    paper = parse_paper_zheng()
    lines = ["# Claim-consistency report", "",
             "Generated from archived TSVs against `docs/paper/manuscript.md`.", ""]
    errors = []
    for (arm, cov), expected in paper.items():
        actual = zheng_point(raw, arm, cov)
        checks = [
            ("n finite", actual.get("n"), expected["n"], 0),
            ("mean QC passes", actual.get("qc_mean"), expected["qc_mean"], 0.051),
            ("Pearson r", actual.get("r"), expected["r"], 0.0015),
            ("Spearman rho", actual.get("rho"), expected["rho"], 0.0015),
            ("RMSE", actual.get("rmse"), expected["rmse"], 0.0015),
            ("slope", actual.get("slope"), expected["slope"], 0.0015),
        ]
        for name, got, want, tol in checks:
            ok = got is not None and close(got, want, tol)
            status = "PASS" if ok else "FAIL"
            if not ok:
                errors.append(f"{arm} {cov}x {name}: paper={want} data={got}")
            lines.append(f"- {status} {arm} {cov:g}x {name}: paper={want}, data={got}")
    # QC seed/all-seed statements.
    qc = pd.read_csv(REF / "results/zheng_qc_pass_counts.tsv", sep="\t")
    a2 = qc[(qc.arm == "A") & (qc.coverage == 2.0)].iloc[0]
    a5 = qc[(qc.arm == "A") & (qc.coverage == 5.0)].iloc[0]
    a10 = qc[(qc.arm == "A") & (qc.coverage == 10.0)].iloc[0]
    expected_qc = [
        ("2x seed counts", [a2.pass_s0, a2.pass_s1, a2.pass_s2], [8,8,5]),
        ("5x seed counts", [a5.pass_s0, a5.pass_s1, a5.pass_s2], [14,11,13]),
        ("10x seed counts", [a10.pass_s0, a10.pass_s1, a10.pass_s2], [14,15,13]),
        ("2x all-seed", a2.pass_all_seeds, 2),
        ("5x all-seed", a5.pass_all_seeds, 8),
        ("10x all-seed", a10.pass_all_seeds, 10),
    ]
    for name, got, want in expected_qc:
        got_values = got if isinstance(got, list) else [got]
        want_values = want if isinstance(want, list) else [want]
        ok = list(map(int, got_values)) == list(map(int, want_values))
        status = "PASS" if ok else "FAIL"
        if not ok: errors.append(f"{name}: paper={want}, data={got}")
        lines.append(f"- {status} A {name}: paper={want}, data={got}")
    # A4 recovered slope and CIs.
    a4 = pd.read_csv(REF / "a4/results/genome_slopes_by_depth.tsv", sep="\t")
    ci = pd.read_csv(REF / "results/a4_slope_bootstrap.tsv", sep="\t")
    expected_a4 = {(arm,d):(s,lo,hi) for arm,d,s,lo,hi in zip(ci.arm,ci.depth,ci.slope,ci.slope_ci_low,ci.slope_ci_high)}
    pattern = re.compile(r"^\| (pois|nb) \| ([0-9.]+)× \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \|$", re.M)
    text = (ROOT / "docs/PAPER_RESULTS.md").read_text()
    for arm, depth, paper_slope, paper_lo, paper_hi in pattern.findall(text):
        row = a4[(a4.arm == arm) & (a4.depth == float(depth)) & (a4.variant == "recovered")].iloc[0]
        slo, lo, hi = expected_a4[(arm,float(depth))]
        checks = [("slope", row.genome_slope_est_on_truth, float(paper_slope)),
                  ("CI low", lo, float(paper_lo)), ("CI high", hi, float(paper_hi))]
        for name, got, want in checks:
            ok = close(got, want, 0.0015)
            if not ok: errors.append(f"A4 {arm} {depth}x {name}: paper={want}, data={got}")
            lines.append(f"- {'PASS' if ok else 'FAIL'} A4 {arm} {depth:g}x {name}: paper={want}, data={got}")
    # Mixed-strain headline.
    mixed_path = ROOT / "benches/mixedstrain_20261001/results/sim_cells.tsv"
    mixed = pd.read_csv(mixed_path, sep="\t")
    mixed_checks = [
        ("mixed profiles", len(mixed), 18),
        ("mixed recall", mixed.recall.min(), 1.0),
        ("mixed spurious", mixed.spurious.sum(), 0),
        ("mixed RMSE", mixed.rmse.mean(), 0.140),
        ("mixed bias", mixed.bias.mean(), 0.015),
    ]
    for name, got, want in mixed_checks:
        ok = close(got, want, 0.0015)
        if not ok: errors.append(f"{name}: paper={want}, data={got}")
        lines.append(f"- {'PASS' if ok else 'FAIL'} {name}: paper={want}, data={got}")
    # Real-community headline claims.
    real_text = PAPER.read_text()
    real_path = ROOT / "benches/realcommunity_20261001/remote_summaries"

    c1 = pd.read_csv(real_path / "c1b_summary_pooled.tsv", sep="\t")
    c1_checks = [
        ("C1b evaluable cells at 0.5x", c1.loc[(c1.arm == "A_sk2bgrow") & (c1.depth == 0.5), "recall"].iloc[0], 18, 0),
        ("C1b evaluable cells at 10x", c1.loc[(c1.arm == "A_sk2bgrow") & (c1.depth == 10.0), "recall"].iloc[0], 18, 0),
    ]
    for depth, want in [(0.5, 0.680), (1.0, 0.278), (2.0, 0.561), (5.0, 0.663), (10.0, 0.675)]:
        got = c1.loc[(c1.arm == "A_sk2bgrow") & (c1.depth == depth), "pearson_mu"].iloc[0]
        c1_checks.append((f"C1b sk2bGrow r at {depth:g}x", got, want, 0.0006))
    for depth, want in [(1.0, 0.321), (2.0, 0.868), (5.0, 0.736), (10.0, 0.787)]:
        got = c1.loc[(c1.arm == "C_relaxed") & (c1.depth == depth), "pearson_mu"].iloc[0]
        c1_checks.append((f"C1b Pilea gates-off r at {depth:g}x", got, want, 0.0006))
    for name, got, want, tol in c1_checks:
        token = str(int(want)) if float(want).is_integer() else f"{want:.3f}"
        ok = close(got, want, tol) and token in real_text
        if not ok: errors.append(f"{name}: data={got}, paper_value={want}, text_present={token in real_text}")
        lines.append(f"- {'PASS' if ok else 'FAIL'} {name}: data={got}, paper={want}")

    c4 = pd.read_csv(real_path / "c4_summary.tsv", sep="\t").set_index("arm")
    c4_checks = [
        ("C4 sk2bGrow samples", c4.loc["A_sk2bgrow", "n_samples_ok"], 20),
        ("C4 sk2bGrow any-estimate MAGs", c4.loc["A_sk2bgrow", "n_mags_any"], 51),
        ("C4 Pilea any-estimate MAGs", c4.loc["C_pilea_default", "n_mags_any"], 64),
        ("C4 sk2bGrow protocol MAGs", c4.loc["A_sk2bgrow", "n_mags_protocol"], 2),
        ("C4 Pilea protocol MAGs", c4.loc["C_pilea_default", "n_mags_protocol"], 18),
        ("C4 sk2bGrow protocol median r", c4.loc["A_sk2bgrow", "median_r_protocol"], 0.438),
        ("C4 Pilea protocol median r", c4.loc["C_pilea_default", "median_r_protocol"], 0.520),
    ]
    for name, got, want in c4_checks:
        ok = close(got, want, 0.0006) and (str(int(want)) in real_text if float(want).is_integer() else f"{want:.3f}" in real_text)
        if not ok: errors.append(f"{name}: data={got}, paper={want}")
        lines.append(f"- {'PASS' if ok else 'FAIL'} {name}: data={got}, paper={want}")

    c5 = pd.read_csv(real_path / "c5_sample_summary.tsv", sep="\t")
    c5cost = pd.read_csv(real_path / "c5_cost_per_sample.tsv", sep="\t")
    sk = c5[c5.arm == "sk2bgrow"]
    pilea_default = c5[c5.arm == "pilea_default"]
    pilea_gatesoff = c5[c5.arm == "pilea_gatesoff"]
    c5_checks = [
        ("C5 sk2bGrow expected cells", len(sk) * 522, 4698),
        ("C5 sk2bGrow estimates", sk.n_est.sum(), 4698),
        ("C5 sk2bGrow QC passes", sk.n_qc_pass.sum(), 484),
        ("C5 sk2bGrow suspicious estimates", sk.n_suspicious.sum(), 0),
        ("C5 sk2bGrow mean recall", sk.recall.mean(), 1.0),
        ("C5 Pilea default estimates", pilea_default.n_est.sum(), 333),
        ("C5 Pilea default mean recall", pilea_default.recall.mean(), 0.071),
        ("C5 Pilea gates-off estimates", pilea_gatesoff.n_est.sum(), 1432),
        ("C5 Pilea gates-off mean recall", pilea_gatesoff.recall.mean(), 0.914),
        ("C5 sk2bGrow mean wall seconds", c5cost.sk2bgrow_total_min.mean() * 60, 76866),
        ("C5 Pilea default mean wall seconds", pilea_default.pilea_default_wall_s.mean(), 536),
        ("C5 Pilea gates-off mean wall seconds", pilea_gatesoff.pilea_gatesoff_wall_s.mean(), 2379),
    ]
    for name, got, want in c5_checks:
        formatted = f"{int(want):,}" if float(want).is_integer() else f"{want:.3f}"
        ok = close(got, want, 0.7) and formatted in real_text
        if not ok: errors.append(f"{name}: data={got}, paper={want}, text_token={formatted}")
        lines.append(f"- {'PASS' if ok else 'FAIL'} {name}: data={got}, paper={want}")

    real_expected_tokens = [
        "275778f350b87e10c6366bf90884d965fbab45a6",
        "PRJNA1280254", "PRJNA551656", "PRJNA974210",
    ]
    for token in real_expected_tokens:
        ok = token in real_text
        if not ok: errors.append(f"real-community provenance token absent: {token}")
        lines.append(f"- {'PASS' if ok else 'FAIL'} real-community provenance token: {token}")

    out = REF / "results/claim_consistency_report.md"
    out.write_text("\n".join(lines) + "\n")
    if errors:
        print("CLAIM CONSISTENCY FAILED")
        for e in errors: print(" -", e)
        raise SystemExit(1)
    print(f"PASS: {len(lines)-3} checked claims; report={out}")

if __name__ == "__main__":
    main()
