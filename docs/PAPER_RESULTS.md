# Current-implementation results for the manuscript

These are the 2026-09-30 HPC refresh results. Use them in preference to the
2026-08-24 numbers in the benchmark notes. The profiling code, statistical
code and the final analyzer fixes are contained in source archive
`sk2bgrow-source-20261001.tar.gz`
(`sha256:8a157b800887df43d5417b901e229fc302554418c119a682499f129ea01c92a7`).
The archived result set is
`sk2bgrow-refresh-results-20260930.tar.gz`
(`sha256:4b5840dff3681aa0ca3c639e48d26c0cbcf4de48bc20eae6d7f3404da25315d1`).
Machine-readable copies are in `benches/refresh_20260930/`.
Review-driven uncertainty diagnostics and accepted scope changes are in
[`PEER_REVIEW_RESPONSE.md`](PEER_REVIEW_RESPONSE.md).
The formal draft, figures, supplement and pre-submission audit are under
[`paper/`](paper/).

## Zheng 2020 E. coli benchmark

### Design

We used 16 growing-media conditions plus a stationary RUN_OUT control from
Zheng et al. (2020) for *E. coli* K-12 MG1655. Condition `s0` is the archived
C1 paired subsample; `s1` and `s2` are deterministic paired resamples of the
C1 20× pool. Every cell was run at 0.5×, 1×, 2×, 5× and 10× with five arms:

| arm | sketch | estimator |
|---|---|---|
| A | 16-enzyme anchors | adaptive windows + V-fit |
| E | FracMinHash | adaptive windows + V-fit |
| B | 16-enzyme anchors | 25 kb windows + sorted/RANSAC |
| C relaxed | FracMinHash | Pilea with gates off |
| C default | FracMinHash | Pilea shipped defaults |

Pilea was the pinned vendor `pilea138` environment used by the prior C1 run
(v1.3.8). The sk2bGrow binary reports version 0.1.0; exact provenance is the
source archive hash below.

The primary accuracy target is the independently measured growth rate
`lambda`. RMSE uses the paper's predicted log2 PTR and is therefore a
same-reads comparison, not absolute accuracy. For each arm/depth we averaged
the three seed estimates within a medium and then computed statistics across
the 16 growing media. `n` is the number of finite estimates. QC-passed count
uses sk2bGrow's shipped default QC.

### Estimator performance

| arm | depth | n finite | mean default-QC passes /16 | Pearson r | Spearman rho | RMSE vs predicted | slope |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 0.5× | 16 | 0 | 0.911 | 0.832 | 0.398 | 0.565 |
| A | 1× | 16 | 0 | 0.960 | 0.932 | 0.280 | 0.850 |
| A | 2× | 16 | 7.0 | 0.966 | 0.932 | 0.161 | 0.862 |
| A | 5× | 16 | 12.7 | 0.974 | 0.965 | 0.052 | 0.914 |
| A | 10× | 16 | 14.0 | 0.967 | 0.962 | 0.055 | 0.940 |
| E | 0.5× | 16 | 0.3 | 0.852 | 0.788 | 0.390 | 0.701 |
| E | 1× | 16 | 0 | 0.923 | 0.950 | 0.312 | 0.964 |
| E | 2× | 16 | 5.0 | 0.962 | 0.915 | 0.131 | 0.896 |
| E | 5× | 16 | 15.7 | 0.968 | 0.965 | 0.051 | 0.924 |
| E | 10× | 16 | 16 | 0.960 | 0.968 | 0.069 | 0.892 |
| B | 0.5× | 16 | 0 | 0.130 | 0.194 | 1.185 | 0.035 |
| B | 1× | 16 | 0 | 0.571 | 0.500 | 0.992 | 0.124 |
| B | 2× | 16 | 3.0 | 0.774 | 0.800 | 0.704 | 0.285 |
| B | 5× | 16 | 5.7 | 0.889 | 0.903 | 0.356 | 0.589 |
| B | 10× | 16 | 7.0 | 0.907 | 0.903 | 0.252 | 0.758 |
| C relaxed | 0.5× | 16 | 16 | undefined | undefined | undefined | 0.000 |
| C relaxed | 1× | 16 | 16 | 0.875 | 0.888 | 0.527 | 0.505 |
| C relaxed | 2× | 16 | 16 | 0.873 | 0.885 | 0.275 | 0.671 |
| C relaxed | 5× | 16 | 16 | 0.980 | 0.979 | 0.101 | 0.807 |
| C relaxed | 10× | 16 | 16 | 0.973 | 0.965 | 0.106 | 0.820 |
| C default | 0.5–5× | 0 | 0 | — | — | — | — |
| C default | 10× | 16 | 16 | 0.973 | 0.965 | 0.106 | 0.820 |

The A-arm Pearson correlation ranges over the three subsample seeds were
0.716–0.899 at 0.5×, 0.891–0.964 at 1×, 0.952–0.972 at 2×, 0.964–0.977 at
5× and 0.963–0.969 at 10×. Thus the low-coverage point is materially more
variable than the 2–10× points.

### Interpretation

1. **Do not report low-coverage sk2bGrow as deployable.** Although the
   all-finite A arm reaches r=0.911 at 0.5× and r=0.960 at 1×, shipped default
   QC passes a mean of 0/16 media at both depths and 7/16 at 2×.
2. **The coordinate fit, not the rank fit, supplies low-depth signal.** At 1×,
   A reaches r=0.960 whereas the otherwise anchor-based B arm reaches r=0.571.
3. **Under identical all-finite rules, Pilea gates-off is competitive.** The
   A-versus-E correlation advantage is largest at 0.5× (0.911 vs 0.852) and
   1× (0.960 vs 0.923), and is small at 2–10×. The main low-depth difference
   is recall/reporting: Pilea defaults return no estimate below 10×, whereas
   sk2bGrow returns a finite estimate, but fails its own default QC below 2×.
4. **Correlation must not be quoted without slope and QC.** The A-arm slope is
   compressed at 0.5× (0.565) and 1× (0.850).

The RUN_OUT control supports this caution. Its A-arm log2 PTR values across
three seeds are 0.113/-0.015/0.083 at 0.5×, 0.033/-0.036/0.076 at 1×, and
0.022/0.013/0.038 at 10×. Sorted regression (B) falsely reports 2.12/2.26/1.99
at 0.5× and 1.94/1.93/1.94 at 1×.

### Review diagnostics and uncertainty

After the initial refresh, we added a medium-level bootstrap diagnostic
(10,000 resamples; seed estimates are averaged within a resampled medium).
This interval expresses variability across media, not biological replication.
Seed-to-seed ranges are reported separately.

| arm | depth | Pearson r (95% bootstrap CI) | slope (95% bootstrap CI) |
|---|---:|---:|---:|
| A | 0.5× | 0.911 (0.820–0.966) | 0.565 (0.438–0.686) |
| A | 1× | 0.960 (0.932–0.980) | 0.850 (0.679–0.986) |
| A | 2× | 0.966 (0.926–0.987) | 0.862 (0.696–0.983) |
| A | 5× | 0.974 (0.939–0.991) | 0.914 (0.808–1.005) |
| A | 10× | 0.967 (0.924–0.990) | 0.940 (0.819–1.060) |
| E | 0.5× | 0.852 (0.713–0.939) | 0.701 (0.438–0.961) |
| E | 1× | 0.923 (0.814–0.977) | 0.964 (0.734–1.155) |
| E | 2× | 0.962 (0.911–0.989) | 0.896 (0.783–0.993) |
| E | 5× | 0.968 (0.920–0.993) | 0.924 (0.812–1.020) |
| E | 10× | 0.960 (0.927–0.986) | 0.892 (0.761–1.057) |

The paired A-minus-E bootstrap warns against over-interpretation. At 0.5×, the
Pearson advantage is 0.059 (CI -0.031 to 0.114); at 1× it is 0.037
(CI 0.010 to 0.057). RMSE differences are small and mostly include zero. At
0.5×, A has the *worse* slope (delta -0.136; CI -0.263 to -0.015). Therefore
the supported statement is better recall and rank recovery at low depth, not a
uniform, significant accuracy win.

QC pass counts also vary by seed. For arm A, exact counts are 8/8/5 at 2×,
14/11/13 at 5× and 14/15/13 at 10×. Only 2, 8 and 10 media respectively pass
in all three seeds. The mean pass counts above are not a deployment guarantee.

## A4 / window-rate experiment

This is the same real-anchor grid as the earlier A4 experiment: Poisson counts
and overdispersed NB counts; true tent amplitude 0–2; depths 0.5/1/2/5×; three
read seeds. The table reports the slope of fused estimated log2 PTR on true
amplitude (`n=15` growing samples per cell). Exact rates passed through the
same windows and fit machinery are shown as a pipeline check.

| arm | depth | recovered slope | exact-rate slope |
|---|---:|---:|---:|
| Poisson | 0.5× | 0.774 | 0.999 |
| Poisson | 1× | 0.958 | 0.999 |
| Poisson | 2× | 1.012 | 0.999 |
| Poisson | 5× | 0.995 | 0.999 |
| NB | 0.5× | 0.701 | 0.967 |
| NB | 1× | 0.853 | 0.997 |
| NB | 2× | 0.857 | 0.978 |
| NB | 5× | 0.982 | 0.992 |

After the fixed-origin signed-slope fix, the stationary b=0 floor is small:
mean estimates are 0.099/0.055/0.035/0.026 for Poisson at
0.5/1/2/5× and 0.033/0.050/0.084/0.036 for NB. Thus the earlier
downhill-only floor is no longer the dominant low-depth defect. Remaining
compression is modest in well-specified Poisson data, but persists in the
overdispersed NB arm at ≤1× and is accompanied by upward bias in low-rate
windows (median log2 recovered/true +0.49 in the true-rate <0.25 bin at 0.5×).

Replicate bootstrap CIs preserve the same conclusion: Poisson recovered slopes
are 0.774 (0.720–0.823) at 0.5×, 0.958 (0.886–1.033) at 1×,
1.012 (0.998–1.026) at 2× and 0.995 (0.981–1.008) at 5×. NB recovered slopes
are 0.701 (0.658–0.741), 0.853 (0.807–0.899), 0.857 (0.807–0.907) and
0.982 (0.951–1.014). Thus the NB defect is clearer than random seed noise,
whereas the Poisson point is close to one by 1×.

## Mixed-strain validation

The current implementation was evaluated on 18 simulated communities from 16
bacterial references: 4, 8 or 16 strains; 1, 2 or 4× coverage per selected
genome; two replicates. Explicit coordinate fitting was used because several
reference assemblies contained multiple contigs or plasmid entries. Recall was
1.00 in all 18 cells, with no spurious genomes. Aggregate RMSE was 0.140 log2
units and mean bias was 0.015. RMSE decreased with coverage; for 16 strains it
was 0.197, 0.100 and 0.056 at 1, 2 and 4×.

This is controlled attribution evidence. It does not replace real community
validation.
Machine-readable cell and summary tables are in
`benches/mixedstrain_20261001/results/`.

## R3 / fragmented-coordinate QC

The current grid has 27 references (three origin-layout families) × 12 read
cells = 324 profiles. Coordinate fits were run explicitly with
`--method v_shape --min-coverage 0`. Complete and correctly ordered controls
are quiet:

| reference class | true log2 PTR | fires / runs |
|---|---:|---:|
| complete | 0, 0.5, 1, 1.5 | 0 / 36 |
| correctly ordered multi-contig | 0, 0.5, 1, 1.5 | 0 / 72 |
| scrambled | 0 | 1 / 54 |
| scrambled | 0.5 | 15 / 54 |
| scrambled | 1 | 40 / 54 |
| scrambled | 1.5 | 40 / 54 |

At 5–10×, detection is complete for scrambled references with true
log2 PTR ≥1 (72/72). It remains weak at 1× and incomplete at 0.5×. One
stationary scrambled control fires (`scr20RX`, 1×, slope branch), so WCG is
not yet ready as an unconditional production gate. It should remain a
prototype flag, preferably restricted to sufficiently deep data and evaluated
with real-data null margins.

Threshold sensitivity shows that raising the slope-branch `z_short` cut from
2.0 to 2.5 removes the stationary false fire while retaining all 72/72 strong
5–10× detections; total flags fall from 96 to 93. This is a post-hoc
sensitivity analysis, not independent validation, so the production gate should
remain disabled pending a pre-registered null/positive evaluation.

### EM scope

The Zheng run is single-isolate, so no anchors are shared across genomes. More
generally, the current ZTP PTR model intentionally consumes raw integer anchor
counts and excludes database-shared anchors; it does **not** regress the
fractional EM weights. `<sample>.em.tsv` is an interface for abundance,
containment and downstream assignment tools. Manuscripts must not describe the
current PTR estimate as EM-reassigned shared-anchor coverage.

## Manuscript wording constraints

* Say “16 deterministic enzyme strata, at most ~15 independent channels”;
  do not say “16 independent enzymes”.
* Label Zheng RMSE as a comparison against the paper's marker-frequency
  prediction, not independent absolute accuracy.
* State whether a number is all-finite or default-QC-passed.
* State the three-seed scheme; these are sequencing subsamples, not
  independent biological replicates.
* State the exact archive hashes above. The historical 2026-08-24 numbers must
  not be mixed with this refresh.
* Treat the R3 result as prototype QC evidence, not a deployed production gate.
