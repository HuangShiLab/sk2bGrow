# R3 — a QC statistic that can see a destroyed coordinate

**Question** (RESEARCH_PLAN Part 4, R3): can a statistic separate "no gradient
because the culture is not growing" (legitimate, pass) from "no gradient
because the coordinate is scrambled" (failure, flag)? Cochran's Q cannot: a
destroyed coordinate makes all sixteen enzymes agree there is no gradient, and
100% of fragmented estimates pass QC at 5–10×.

**Answer: yes at 5–10× for strong gradients, but not production-ready at 1×.**
The original one-seed grid fired on 47/48 scrambled cells with log2PTR ≥1 at
5–10× and stayed quiet on every complete, ordered, and stationary control. The
larger 2026-09-30 current grid confirms those controls and gives 72/72
detection for scrambled log2PTR≥1 at 5–10×, but also finds one stationary
false fire and weak 1× sensitivity (see “Current-refresh result” below). WCG is
therefore best described as a promising prototype QC statistic, not a deployed
gate.

Prototype only; nothing is wired into the pipeline. Implementation:
`qc_within_contig.py`; grid runner `r3_run.sh`; scorer `r3_score.py`;
references `r3_refs.py`; reads `r3_simulate.py`.

Version note (2026-09-29): the numeric grid below predates the fail-closed
fragmented-reference rule and fixed-origin signed-slope change. `r3_run.sh` now
passes `--method v_shape --min-coverage 0` explicitly, but the published R3
numbers must be regenerated at the exact final commit before manuscript use.

## Current-refresh result (2026-09-30)

The current grid has 27 references (three origin-layout families) × 12 read
cells = 324 profiles. Complete and correctly ordered controls are quiet:
0/36 and 0/72, including stationary controls. Scrambled results are:

| true log2PTR | fires / runs | 1× | 5× | 10× |
|---:|---:|---:|---:|---:|
| 0 | 1/54 | 1/18 | 0/18 | 0/18 |
| 0.5 | 15/54 | 1/18 | 6/18 | 8/18 |
| 1.0 | 40/54 | 4/18 | 18/18 | 18/18 |
| 1.5 | 40/54 | 4/18 | 18/18 | 18/18 |

Thus current sensitivity is excellent at 5–10× for true log2PTR≥1, but weak at
1× and incomplete at 0.5×. One scrambled stationary control (`scr20RX`, 1×)
fires on the slope branch. WCG therefore remains a prototype flag; do not wire
it into production QC without deeper control data and threshold recalibration.
The pre-fix grid below is retained for mechanism/provenance only.

## The statistic

Scrambling destroys the coordinate *between* contigs only. Two things survive
inside contigs, and the two branches exploit one each. Both are computed from
the `windows.rates.tsv` every profile run already writes, plus `log2(PTR)` and
`se` from `output.tsv`.

**Slope branch — local amplitude vs global amplitude.** Enzymes are
median-centred genome-wide; per contig c a weighted line
`y = a_c + m_c·x` (x contig-local bp, weights 1/log2_se², slope SE rescaled by
the contig's residual χ²/dof so mis-calibrated window SEs do not produce
artificially tight slopes). The noise-corrected RMS slope, rescaled by half
the genome length, is a local amplitude estimate in log2PTR units:

```
A_loc = sqrt(max(0, Σw m_c²/Σw − k/Σw)) · (L/2),   w_c = 1/s_c²
```

`k/Σw` is E[Σw m²/Σw] under "no slopes", so a stationary sample returns ~0
rather than the spread of its own noise. A_loc uses only within-contig
coordinates — contig order and orientation are irrelevant to it. Compare with
the pipeline's fused V-fit amplitude `b_hat`:

```
z_loc   = A_loc / se(A_loc)                       is there a local gradient?
z_short = (A_loc − b_hat) / sqrt(se_A² + se_b²)    did the global fit lose it?
fire    = usable AND z_loc > 2 AND z_short > 2 AND A_loc − b_hat > 0.5
```

`usable` requires ≥ 2 usable contigs with median ≥ 20 windows per contig
(≈1.25 panel passes; on this 16-enzyme panel that holds up to N ≈ 20). Below
that information density a contig holds ~1 window per enzyme and "slope" is
not separable from per-enzyme offsets — measured, not theoretical: without the
domain restriction A_loc reads 2–3 on *correctly ordered* 100-contig controls.

**Jump branch — boundary continuity.** Within each enzyme series ordered by
reference coordinate, consecutive window pairs give a step d = y₂−y₁ with
noise floor f = s₁²+s₂²; e = min(d², 9f) − f estimates the squared true step
(Winsorized at 3σ so one wild low-count window cannot dominate; pairs with
log2_se > 1 carry no step information and are dropped).

```
J   = mean(e | boundary-crossing) − mean(e | within-contig)
z_J = J / se(J),  se empirical from per-pair scatter
fire = ≥ 8 boundary pairs AND z_J > 3.5 AND J > 0.03
```

A correct coordinate is continuous, so J ≈ 0 regardless of growth. A scrambled
one jumps to an unrelated coverage level at every boundary: E[J] ≈
log2PTR²/6 − (within-contig excess). Power *grows* with N — complementing the
slope branch, which fades. The J > 0.03 margin exists because windows never
straddle contigs: boundary-flanking windows sit ~2× further apart than
within-contig pairs, and spatially correlated residual wiggle (GC correction
leftovers) gives a null floor of ~0.01–0.02 at 10×. Observed max J on correct
references: 0.013.

**Why stationary can never fire:** no gradient means no within-contig slopes
(A_loc ≈ 0, noise-corrected) and no boundary jumps (all differences are
sampling noise, subtracted by the floor). Both branches require a gradient to
exist locally before they can accuse the global coordinate of losing it. That
is precisely the asymmetry Cochran's Q lacks.

## Evaluation

Simulation: reads with a V-shaped gradient (Pilea's Methods model, as in
`benches/simulate.py`; ori at position 0, ter at midpoint) on K-12 at
log2PTR ∈ {0, 0.5, 1.0, 1.5} × depth {1, 5, 10}×, fixed seeds. References: the
same chromosome cut into N ∈ {2, 5, 10, 20, 50, 100} lognormal contigs
(fragment.py protocol, seed 0), shuffled and flipped (`scr`); the same cuts in
true order/orientation (`ord`, the correctly-scaffolded control); the complete
chromosome. Three layout families: `cut` (ori at a fragment boundary — the
simulator's position 0 is always a cut), `rot` (genome rotated before both
steps; same property), `ori-int` (reference rotated, reads not: ori interior
to a contig, like the real-data sweep). 228 profile runs total; statistic
scored offline from run outputs.

Fused V-fit estimate (`*` = WCG fires), 5–10× — the regime where the existing
QC passes fragmented estimates:

| reference | layout | 0.5@5× | 0.5@10× | 1.0@5× | 1.0@10× | 1.5@5× | 1.5@10× |
|---|---|---|---|---|---|---|---|
| complete | | 0.48 | 0.51 | 1.04 | 1.00 | 1.46 | 1.48 |
| ordered N=10 | | 0.49 | 0.51 | 1.04 | 1.00 | 1.46 | 1.48 |
| ordered N=100 | all 3 | 0.47–0.49 | 0.48–0.51 | 1.01–1.04 | 0.98–1.00 | 1.42–1.47 | 1.46–1.48 |
| scr N=2 | cut/rot | 0.16–0.17* | 0.13–0.15* | 0.31–0.32* | 0.28–0.29* | 0.35–0.38* | 0.37–0.40* |
| scr N=2 | ori-int | 0.38 | 0.38 | 0.72* | 0.74* | 1.02 | 1.06* |
| scr N=5 | cut | 0.23 | 0.26 | 0.48* | 0.47* | 0.57* | 0.64* |
| scr N=10 | cut | 0.25* | 0.29* | 0.55* | 0.56* | 0.79* | 0.82* |
| scr N=20 | cut/rot | 0.26–0.28 | 0.29–0.30 | 0.62* | 0.57–0.59* | 0.88–0.90* | 0.90–0.91* |
| scr N=20 | ori-int | 0.05* | 0.08* | 0.08* | 0.08* | 0.07* | 0.09* |
| scr N=50 | cut | 0.06 | 0.06 | 0.11* | 0.05* | 0.12* | 0.06* |
| scr N=100 | all 3 | 0.05–0.16 | 0.07–0.18* | 0.10–0.33* | 0.09–0.32* | 0.09–0.40* | 0.11–0.44* |

Detection summary:

| cell class | fires |
|---|---|
| complete, growing, all depths | 0/27 |
| ordered multi-contig, growing, all depths | 0/36 |
| stationary (log2PTR = 0), any reference, all depths | **0/57** |
| scrambled, growing, 5–10×, log2PTR ≥ 1 | **47/48** |
| scrambled, growing, 5–10×, log2PTR = 0.5 | 11/24 |
| scrambled, growing, 1× | 11/36 — but see below |

The one miss at log2PTR ≥ 1, 5–10× is the mildest cell in the grid (N=2,
ori-int, log2PTR = 1.5, 5×: estimate 1.02 vs truth 1.5, J z = 3.34 against the
3.5 threshold).

## Sensitivity limit

The gate fires when scrambling destroys ≳ 0.3–0.5 log2 of fitted amplitude,
and is quiet below that. Read against the established dose-response (plan Part
3.1: bias −0.06 at N=2 → −0.81 at N=100, smooth and monotone), the statistic
degrades smoothly too, in the same direction:

- **N = 2 is layout-dependent, and the gate tracks the damage, not the contig
  count.** With ori at a cut boundary (both cut/rot families here) N = 2 is
  catastrophic (estimates 25–30% of truth) and fires 15/18. With ori interior
  (matching the real-data sweep, where N=2 cost only slope 0.95→0.88) the
  estimate keeps 68–76% of truth and the gate fires only on the
  most-damaged cells (3/9, all with J z ≥ 3.3). This is the desired behaviour:
  the gate measures coordinate destruction, not fragmentation.
- **log2PTR = 0.5 (PTR 1.41)** is the weak-amplitude floor: the whole gradient
  is 0.5 log2 wide, so even total destruction at N = 5–50 costs ≤ 0.25 log2 —
  inside the margin the gate needs to distinguish destruction from the two
  estimators' ordinary disagreement. Detected at N = 2, 10 and N ≥ 100 at 10×;
  missed at N = 5/20/50 (estimates 0.23–0.30 vs truth 0.5, PTR error ~15%).
- **1× is the A4 floor.** Window SEs at near-zero counts are not calibrated
  (the known A4 suspect layer); A_loc inflates ~20% even with the χ² rescale.
  All 1× cells here also fail the existing coverage gate anyway (anchor
  coverage 0.83 < 1.0), so the blind spot is in practice a ≥ 2× phenomenon.
- **Slope branch domain:** N ≲ 20 on this panel (needs ≥ ~20 windows/contig).
  **Jump branch domain:** N ≳ 3, strengthens with N; its null floor is the
  boundary-spacing artefact, ~0.01–0.02 at 10×, worse at 1×.

## Recommendation (proposal, not applied)

Wire WCG into the QC as a **flag, not a correction** — the fix remains
`sk2bgrow scaffold` (or the spread estimator at ≥ 5×), and R5's two-estimator
design is unaffected:

1. `python/sk2bgrow/wcg.py`: port `per_contig_slopes`, `local_amplitude`,
   `boundary_jumps`, `wcg_stat` from `benches/fragmentation/qc_within_contig.py`
   verbatim (they take the windows DataFrame and the fused estimate; ~120
   lines, no new dependencies).
2. `cli.py`: after `fusion.fuse_table`, compute `wcg_stat(windows,
   total_reference_length, log2_ptr, se)` per genome and append columns
   `wcg_A_loc, wcg_z_short, wcg_J, wcg_z_J, wcg_flag` to output.tsv (additive;
   Pilea-compatible column order untouched).
3. `report.py::apply_qc`: if `wcg_flag`, add reason
   `"local gradient A_loc={A:.2f} not captured by coordinate fit ({b:.2f});
   scaffold the reference (sk2bgrow scaffold) or use the spread estimator"`.
   Fail-closed like the other gates.
4. Thresholds as evaluated: slope branch z_loc > 2, z_short > 2,
   A_loc − b_hat > 0.5, median ≥ 20 windows/contig; jump branch ≥ 8 boundary
   pairs, z_J > 3.5, J > 0.03. Consider enabling the slope branch only at
   coverage ≥ 2× until A4 is resolved.
   Review sensitivity on the 2026-09-30 grid suggests `z_short > 2.5` removes
   the sole stationary false fire without losing any strong 5–10× detection.
   This is post hoc; validate it on an independent grid before changing the
   default.
5. Regression test: profile the committed simulation grid (or one synthetic
   V-gradient fixture) against a 10-contig shuffled reference and assert the
   flag fires, plus a stationary control asserting it does not.

Do **not** replace this with a contig-count guard: the ordered N=10/N=100
controls (0/36 fires) are exactly the cases a count guard would condemn.

## Reproduce

```bash
cargo build --release
WORK=benches/work/r3
python3 benches/fragmentation/r3_refs.py benches/genomes/Escherichia_coli_K12.fna -o $WORK/refs
python3 benches/fragmentation/r3_simulate.py benches/genomes/Escherichia_coli_K12.fna -o $WORK/reads
WORK=$WORK benches/fragmentation/r3_run.sh        # needs --python with statsmodels
WORK=$WORK python3 benches/fragmentation/r3_score.py
```

Record `git rev-parse HEAD`, the simulator seed, and the R3 stats TSV with any
result quoted from this grid.

`r3_score.py` writes `$WORK/r3_stats.tsv` (one row per run: V-fit estimate,
both branches, flags) and prints the evaluation table.

## Caveats carried

- One genome (K-12), one seed per cell; the grid is 228 runs but not
  replicated across fragmentation seeds. The N=2 layout dependence is the
  known-large variance term and was characterised directly (three families).
- Thresholds are tuned on this grid. The null-side margins (J = 0.013 max on
  correct references vs the 0.03 fire margin; A_loc − b_hat = 0.18 max on
  growing controls vs 0.5) are comfortable but should be re-measured on the
  real-data sweep (Zheng reads exist per media; `sweep.sh` outputs can be
  scored offline by `qc_within_contig.py` with no re-run of the estimator).
- The gate says nothing about *why* the coordinate is broken beyond "local
  gradient not captured"; a heavily mis-assembled single contig (inversions)
  would also fire — arguably correctly.
