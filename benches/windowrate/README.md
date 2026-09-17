# C7 — the A4 window-rate experiment

**Question.** Open issue A4 (`../RESEARCH_PLAN.md` Part 3.3): at low coverage the fitted
slope of estimated log2PTR on truth is compressed (0.681 at 0.5×, 1.074 at 10× on the
Zheng data). One candidate remained: does the ZTP/ZTNB window-rate layer
(`python/sk2bgrow/ztp.py`) bias a window's rate at near-zero counts — before the V-fit
ever sees it — and thereby compress the coverage gradient?

**Verdict, one paragraph.** The compression is reproduced almost exactly (genome-level
slope 0.60 at 0.5× here vs 0.681 in A4) and is shown to enter **entirely with the
recovered window-rate values**: exact rates pushed through the identical window
boundaries, origin search, V-fit and fusion give slope 0.998–0.999 at *every* depth.
But the specific A4 hypothesis is **not** what happens. On well-specified (Poisson)
counts the surviving windows are *median-unbiased* — slightly downward, not upward —
and the per-window `log2(recovered) ~ log2(true)` slope is **> 1** at low depth, not
< 1. The genome-level compression is manufactured downstream of the rates: by
inverse-variance **self-weighting** (a window whose Poisson draw fluctuated up gets a
bigger rate *and* a smaller SE; corr(weight, error) = +0.42 at 0.5×), by asymmetric
Tukey trimming of a one-sided error tail, by a one-sided (downhill-only) floor of
≈ +0.2 at true log2PTR = 0, and — a genuine bug found in passing — by an all-ones
boundary in the ZTP path that hands such windows to the ZTNB branch, which returns a
garbage rate of ~1e-6. Under overdispersion the low-rate windows *are* biased upward
(+0.14 to +0.32 log2 median at ≤1×), but that arm's slope is even more compressed, so
the upward bias is not the compressor either.

## Design

- Real E. coli K-12 MG1655 (`../genomes/Escherichia_coli_K12.fna`), indexed with the
  real pipeline: `sk2bgrow index ... --enzymes all --write-tgt` (43,735 anchors,
  42,234 usable). Anchor coordinates, tag lengths, window cutting, counting and the
  stats layer are all the production code; **nothing in the pipeline was modified**.
- `simulate_tent.py` draws 150 bp single-end reads with a known tent gradient:
  read start `s` is sampled ∝ `2^(−b·d(mid)/(L/2))`, `d` = circular distance from
  ori = 3,923,883. Because the gradient is linear in log space, the per-base coverage
  tent has amplitude exactly `b` = true log2PTR. The **exact** expected count of every
  anchor is `N · Σ p(s)` over the starts whose read covers the whole tag
  (`anchor_expected`, prefix sum) — window truth is the mean of that over the window's
  usable anchors, using the window boundaries the pipeline itself wrote.
- Grid: depths 0.5 / 1 / 2 / 5× (5× = high-depth control) × true log2PTR
  b ∈ {0, 0.5, 1, 1.5, 2} × 3 seeds × 2 arms. Arms: `pois` (Poisson counts) and
  `nb` (per-1 kb lognormal efficiency, σ = 0.5, unit mean — genuinely overdispersed
  within-window counts; the efficiency vector is saved and enters the truth).
  Seeds are fixed functions of the cell; everything is deterministic.
- GC correction is **off** in the primary arms (no GC bias is simulated, so the loess
  could only add noise). A GC-on control (pois, 0.5×/1×, rep 0) changes nothing
  (genome slope 0.596 vs 0.601 at 0.5×; window-rate slopes within noise).
- Per sample the full real pipeline runs: `sk2bgrow profile` (count) then
  `python -m sk2bgrow.cli profile --no-gc-correct`. `analyze.py` reconstructs truth,
  joins it to `windows.rates.tsv` (membership verified against the pipeline's own
  `n_anchors`), refits every window with the same public `estimate_window_rate` to
  recover `n_components`/BIC (refit rates match the pipeline's to 1e-9 — used as a
  consistency check), and re-estimates genome-level log2PTR through the real
  `fit_windows` + `fuse_table` for three rate variants. `decompose.py` changes one
  fit factor at a time through the same machinery.

Two slopes are reported and they are **not** interchangeable (RESEARCH_PLAN Part 2.2):
the **window-rate slope** (OLS of window `log2(recovered)` on `log2(true)` rate) and
the **genome-level slope** (OLS of fused estimated log2PTR on true b, 15 samples per
cell). A4's 0.681 is a genome-level slope.

## Results

### Window-rate layer, per depth (question b)

Window-rate slope `log2(recovered) ~ log2(true)` per (replicate, enzyme), b > 0
cells (median over 192 fits; "pooled" = single regression of enzyme-and-sample
demeaned values):

| arm | depth | pooled | median per-fit | mean per-fit |
|---|---|---:|---:|---:|
| pois | 0.5× | 2.16 | 1.35 | 2.02 |
| pois | 1× | 1.21 | 1.10 | 1.20 |
| pois | 2× | 1.14 | 1.06 | 1.15 |
| pois | 5× | 1.01 | 1.01 | 1.00 |
| nb | 0.5× | 1.28 | 0.78 | 1.03 |
| nb | 1× | 1.07 | 0.97 | 1.20 |
| nb | 2× | 1.16 | 0.99 | 1.18 |
| nb | 5× | 0.98 | 0.98 | 0.97 |

The confirmation criterion "slope < 1 at low depth" **does not occur**. If anything
the surviving-window gradient is steepened at 0.5× (the lowest windows read slightly
low in the median).

### Bias by true-count regime (questions a, b; `bias_by_regime.tsv`)

Median `log2(recovered/true)` and fraction of windows the layer drops (rate NaN):

| arm | depth | true < 0.25 | 0.25–0.5 | 0.5–1 | 1–2 | ≥ 2 |
|---|---|---|---|---|---|---|
| pois | 0.5× | −0.17 (16% dropped) | −0.04 (5%) | −0.02 (1%) | — | — |
| pois | 1× | — | −0.06 (1%) | −0.03 (0%) | −0.02 | — |
| pois | 5× | — | — | — | −0.10* | 0.00 |
| nb | 0.5× | **+0.32** (15% dropped) | **+0.21** (4%) | **+0.17** (1%) | — | — |
| nb | 1× | — | **+0.19** (2%) | **+0.17** (0%) | +0.13 | — |
| nb | 5× | — | — | — | +0.02 | −0.02 |

(*the 1–2 bin at pois 5× has n = 57, pooled over b; treat as noise.)

- **Poisson arm: no upward bias.** Surviving near-zero windows are median-unbiased
  (slightly *down*). The ZTP MLE is fine on well-specified counts; the damage is
  scatter (std ≈ 0.7–0.85 log2 below true rate 0.5) plus two one-sided failure modes.
- **Overdispersed arm: upward bias confirmed** at low true rates (+0.13 to +0.32
  log2, i.e. +10–25%, decaying with depth; gone by 5×). The truncation model
  (Poisson/NB) is mis-specified for lognormal site heterogeneity and the excess zeros
  it predicts get compensated by an inflated mean. Note the eff vector is 1 kb
  lognormal, not gamma — this is the realistic mis-specification case.
- **Floor, not ceiling:** recovered rates have no upper pathologies, but ~3.6% of
  windows at 0.5× come back at ~1e-8 (min recovered 1.4e-8 vs min true 0.19) — see
  the bug below. 5% of all 0.5× windows sit below 0.09 while the true floor is 0.16.

### The all-ones / ZTNB boundary bug (question c, and a genuine defect)

At 0.5×, ~16% of windows with true rate < 0.25 have all positive counts equal to 1.
For such a window `solve_ztp_lambda(mean = 1.0)` returns 0.0, and `fit_ztp_mixture`
then hard-codes `loglik = −inf` (`ztp.py:195`) — although the true ZTP likelihood at
the boundary is *finite* (P(X=1 | X≥1) → 1 as λ → 0⁺, so loglik → 0). BIC therefore
always prefers the ZTNB branch, whose MLE drives μ → ~1e-6 with a huge α and returns
that as a **finite** rate (log2 ≈ −20) instead of NaN. 99% of "collapsed" windows
(log2 ratio < −3) are `model == "ztnb"`; 81% of them (pois, 0.5×) are exactly
all-ones. The collapse frequency tracks the all-ones frequency, not the BIC component
count: `mean n_components` is 1.00–1.01 everywhere in the pois arm (mixtures are
essentially never selected; `corr(n_components, |bias|)` ≈ 0.06–0.33 is driven by the
few multi-component fits at high depth, not by the bias).

**Reported, not fixed** (per the experiment's ground rules). Two candidate fixes:
evaluate the ZTP loglik at the λ → 0⁺ limit instead of −inf, or return NaN from
`estimate_window_rate(model="auto")` when the positive counts are all ones (the
window carries no rate information either way — the pure-ZTP path already returns
NaN there).

Downstream impact is *limited but real*: Tukey trimming catches most collapsed
windows, and hand-removing them moves the genome slope only 0.601 → 0.637 at 0.5×.
Without trimming they destroy fits outright (the no-trim arm needs them removed by
hand to be interpretable at all).

### Genome-level slope and the decomposition (question d; `decomposition.tsv`)

Genome-level slope, estimated log2PTR ~ true b (15 samples per cell; b0 floor = mean
estimate at true b = 0):

| variant (pois arm) | 0.5× | 1× | 2× | 5× |
|---|---:|---:|---:|---:|
| **recovered (pipeline as shipped)** | **0.601** | 0.874 | 0.953 | 0.977 |
| true rates, all windows | 0.999 | 0.999 | 0.999 | 0.999 |
| true rates, same surviving windows | 0.998 | 0.999 | 0.999 | 0.999 |
| recovered, flat weights (const SE) | 0.692 | 0.921 | 0.994 | 0.987 |
| recovered, origin fixed at truth | 0.640 | 0.865 | 0.953 | 0.978 |
| recovered, flat SE + ori fixed | 0.741 | 0.922 | 0.990 | 0.989 |
| recovered, flat SE + no Tukey* | 0.750 | 0.967 | 1.004 | 0.985 |
| recovered, flat SE + ori fixed + no Tukey* | 0.784 | 0.973 | 1.005 | 0.989 |

(*collapsed bug-artifact windows hand-removed first.) nb arm, recovered: 0.561 /
0.712 / 0.807 / 0.914; its true-rate counterfactual is 0.960–0.979 (the small
shortfall is Jensen's inequality: log2 of a window's *mean* rate is not exactly the
linear tent when per-kb efficiency varies within a window).

The ledger at 0.5× (pois; deficit vs 1.0 is 0.40):

- **Fit machinery: innocent.** True rates through the identical windows, origin
  search, V-fit, trimming and fusion give 0.998–0.999 at every depth in both arms.
- **Window dropout: innocent.** Same dropout pattern, true values: 0.998 — losing
  16% of the lowest windows costs nothing when the survivors are exact. (Dropout is
  not random — it selects low Poisson draws — but it is not *itself* the compressor.)
- **IV self-weighting: −0.09.** The window SE is a decreasing function of the window's
  own fitted rate, so weights correlate with the error: corr(1/SE², error) = **+0.42**
  at 0.5× (+0.50 within the true < 0.5 regime), +0.36 / +0.27 / +0.16 at 1/2/5×.
  Windows that fluctuated *up* get more weight; the noisy terminus is pulled up more
  than the quiet origin, compressing the amplitude. Flat weights recover 0.692.
  (Note: this is the *within-fit* weighting. A4's "IV beats unweighted 0.681 vs
  0.497" compared fused estimates on real data — a different weighting level, not
  contradicted here.)
- **Origin search on noisy pooled rates: −0.04.** Fixing the ori at truth gives 0.640.
  Consistent with A4's "origin error is minor".
- **Tukey trimming asymmetry: −0.05.** With a one-sided error tail, trimming removes
  more from the low side.
- **b = 0 floor: ≈ −0.07 of slope.** The fit rejects uphill solutions, so at flat
  truth the estimates are one-sided noise: mean estimate +0.223 (pois) / +0.305 (nb)
  at 0.5×, shrinking to +0.05 / +0.10 at 5×. Same floor as the real-data negative
  control (0.260 at 0.5×, RESEARCH_PLAN Part 3.1). With b ∈ [0, 2] evenly spaced, an
  intercept of +0.22 costs ≈ 0.32 × 0.22 ≈ 0.07 of OLS slope.
- **Residual ≈ −0.10:** plain Monte-Carlo scatter of the rates (std 0.7–0.85 log2 in
  the sub-0.5 regime) propagating through per-enzyme fits and fusion; even with every
  factor above neutralised the slope is 0.784, not 0.998.

In the nb arm the same structure holds plus the genuine upward rate bias at low
counts (+0.2 log2), which is why its slope is *worse* despite the bias pushing the
"right" way for the A4 hypothesis: the bias lifts the terminus, i.e. it compresses
the fitted amplitude directly.

## What fix the data supports

1. **Fix the all-ones boundary bug** (above). Clear defect, cheap fix, removes the
   ~1e-8 garbage rates that currently depend on Tukey for containment. Does not by
   itself move the slope much (+0.04) but is unambiguous.
2. **Decouple within-fit weights from the window's own rate.** Self-weighting costs
   −0.09 of slope at 0.5× and is the largest *fixable* fit-side term. Options the
   data motivates: weights from the expected information at a shrunken rate (e.g. the
   enzyme median), or two-stage fitting (fit flat, reweight, refit). Needs a
   real-data check before adopting — flat weights also downweight genuinely bad
   windows, and A4's cross-enzyme fusion finding cuts the other way.
3. **Do not "fix" the ZTP model itself for the Poisson case** — it is median-unbiased
   there. The upward bias under overdispersion (+0.2 log2 at ≤1×) is a
   mis-specification effect; a per-anchor efficiency term (or fitting window rates on
   the raw scale rather than via truncation inversion) would address it, but that is
   a design change, not a bug fix.
4. **The b = 0 floor is a property of the downhill-only fit**, visible identically in
   the real-data negative control. Report it as a known bias floor; it is not
   fixable by the rate layer.

For A4's bottom line: the genome-level compression at 0.5× is **fully explained by
the window-rate layer's output** (values + their SEs + failure modes), with ~nothing
left for the fit machinery to answer for — but the mechanism is noise amplification
and one-sided failures, not a systematic truncation bias of well-specified rates.

## After the fixes (2026-08-26)

The fixes the data supported are implemented in `python/sk2bgrow/` (the sections
above are the pre-fix record; `work/results_baseline/` holds the pre-fix TSVs,
`results/` the post-fix ones, `work/results_boundary_only/` an intermediate
state with only the all-ones fix + two-stage weights):

1. **All-ones boundary** (`ztp.py`): the ZTP loglik at λ = 0 is evaluated at its
   finite limit (0.0) instead of −inf, so BIC keeps the ZTP branch and the
   window returns NaN. Two sibling defects found while chasing the residual
   collapses: `ztp_logpmf` computed `log(1−e^−λ)` as `log1p(−exp(−λ))`, which
   loses the mantissa below λ ~ 1e-8 and returns logpmf(1) ≈ +0.69 — phantom
   likelihood that let EM mixture components collapse onto λ ≈ 0 — now computed
   in `expm1` form; and the ZTNB branch has a non-identifiable ridge
   (μ → 0, α → ∞) at low counts because the truncated likelihood never sees the
   zeros — a NB win is now vetoed when the fit's implied detection probability
   misses the observed detected fraction by >10× (`NB_DETECTION_MISMATCH`),
   falling back to the identified ZTP rate.
2. **Two-stage decoupled weights** (`fit.py`, `fit_v_shape(shrink_se=True)`):
   stage 1 fits with the SEs as given; stage 2 re-expresses each window's SE as
   a smooth function of the stage-1 *fitted* log2 rate (log se regressed on
   observed log2 rate, predicted at the fitted values), keeping the
   rate-to-precision trend while breaking the weight-error correlation.

Before → after, genome-level slope of estimated log2PTR on true b (recovered
variant, n = 15 per cell):

| arm | depth | before | after | b0 floor before | b0 floor after |
|---|---|---:|---:|---:|---:|
| pois | 0.5× | 0.601 | **0.757** | +0.223 | +0.246 |
| pois | 1× | 0.874 | **0.908** | +0.097 | +0.174 |
| pois | 2× | 0.953 | **0.988** | +0.079 | +0.091 |
| pois | 5× | 0.977 | **0.985** | +0.051 | +0.050 |
| nb | 0.5× | 0.561 | **0.667** | +0.305 | +0.249 |
| nb | 1× | 0.712 | **0.783** | +0.204 | +0.160 |
| nb | 2× | 0.807 | **0.843** | +0.160 | +0.142 |
| nb | 5× | 0.914 | **0.954** | +0.097 | +0.106 |

Collapsed windows (log2 recovered/true < −3): pois 3.59% → **0.07%** at 0.5×
(1.13% → 0.07% at 1×); nb 3.74% → 0.25% at 0.5×, 4.13% → 0.70% at 1×.

The b=0 floor was re-checked with 11 seeds on identical count tables against the
pre-fix code (git worktree at HEAD): pois 0.5× 0.292 → 0.336, 1× 0.166 → 0.146;
nb 0.5× 0.300 → 0.258, 1× 0.169 → 0.171 — no systematic inflation; the apparent
3-seed uptick at 1× was small-n noise. The floor is the documented one-sided-fit
behaviour and is unchanged in kind.

Remaining, out of scope here: ~100 of 72k windows (mostly nb arm) still collapse
through the **dominant-component mixture rule** — on a bimodal window (e.g. 31
ones plus a handful of 2s/4s) BIC legitimately selects two components, and
Pilea's "highest-weight component" rule reads the all-ones component's λ ≈ 0 as
the window rate. The likelihood is now correct; the rate-selection heuristic is
the issue. A mixture-mean rate would have read 0.49 against a truth of 0.47 on
the inspected case, but changing the rule alters estimator semantics and needs
its own experiment. The decoupled weights + Tukey fence contain these windows;
they no longer move the genome-level slope measurably.

Test state after the fixes: `cargo test` 86/86; `pytest tests/python` 113/113
(7 new regression tests: all-ones boundary, boundary numerics, mixture collapse,
ZTNB ridge fallback, self-weighting resistance, single-stage parity).


## Caveats

One genome; error-free reads (no sequencing error, no GC bias — mismatch tolerance
and GC correction are exercised only as no-ops); the nb arm's overdispersion is
1 kb lognormal, not gamma; n = 15 samples per genome-level cell (3 seeds × 5 truth
values), so genome-level slopes carry Monte-Carlo error of roughly ±0.05; the
decomposition terms are not exactly additive (fit factors interact).

## Reproduce

```bash
cargo build --release
./run.sh          # index (once) + 130 pipeline runs + analyze.py + gc control
../../.venv/bin/python decompose.py   # fit-factor ledger (needs results/per_window.tsv)
```

Inputs, per-sample FASTQs, count tables and stats outputs live under `work/`
(git-ignored). `results/` holds: `per_window.tsv` (72k windows, recovered vs true),
`genome_estimates.tsv`, `genome_slopes_by_depth.tsv`, `window_rate_slopes.tsv`,
`bias_by_regime.tsv`, `model_selection.tsv`, `decomposition.tsv`.
