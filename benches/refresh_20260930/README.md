# 2026-09-30 current-commit refresh (HPC)

This directory archives the completed current-implementation refresh, not the
historical 2026-08-24 benchmark.

## Scope

- **Zheng:** 17 conditions (16 growing media + RUN_OUT) × 3 deterministic
  subsample seeds (`s0` C1 subsample; `s1`/`s2` paired resamples from the C1
  20× pool) × depths 0.5/1/2/5/10× × arms A/B/E/Pilea-default/Pilea-relaxed.
  Final count: **1,275 output.tsv files**.
- **A4/window-rate:** full 120-cell primary grid plus 10 GC controls.
- **R3/fragmentation:** 27 references × 12 read cells = **324 profiles**.

## Run identity

- Final source archive (including the paired-E and analyzer fixes used for all
  completed cells):
  `sha256=8a157b800887df43d5417b901e229fc302554418c119a682499f129ea01c92a7`.
- Results archive:
  `sha256=4b5840dff3681aa0ca3c639e48d26c0cbcf4de48bc20eae6d7f3404da25315d1`.
- HPC root:
  `/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930`.
- Binary: `sk2bgrow 0.1.0`, built from the uploaded working tree.
- Pilea: vendor environment `pilea138` from the prior C1 setup.
- Pilea version: v1.3.8.

## Main result files

- `results/zheng_all_seeds_qc.tsv`: all-finite and default-QC-passed metrics,
  by seed, arm, and coverage.
- `zheng/s*/results_raw.tsv`: raw per-condition estimates.
- `a4/results/`: current-commit window-rate summaries.
- `r3/r3_stats.tsv`: current-commit WCG/R3 grid.
- `results/`: copies of cross-benchmark summary tables.
- `results/zheng_bootstrap_media.tsv`,
  `results/a4_slope_bootstrap.tsv`, and
  `results/r3_threshold_sensitivity.tsv`: post-run review diagnostics.
- `review_checks.py`: deterministic script that regenerates those diagnostics.

The historical numbers in `benches/zheng2020/RESULTS.txt` are not replaced by
these files; use this refresh for current-commit claims and cite the hashes.

Manuscript-facing summaries and wording constraints are in
[`../../docs/PAPER_RESULTS.md`](../../docs/PAPER_RESULTS.md) and
[`../../docs/PAPER_CLAIMS.md`](../../docs/PAPER_CLAIMS.md).
