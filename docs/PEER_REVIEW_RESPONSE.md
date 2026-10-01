# Simulated peer-review response

This is an internal adversarial review of the 2026-09-30 refresh and the
manuscript material in [`PAPER_RESULTS.md`](PAPER_RESULTS.md). It is not an
external journal review. “Resolved” means the concern has been either corrected
or explicitly scoped; it does not imply that the underlying biological or
engineering question is closed.

## Reviewer 1 — correlations are presented without uncertainty

### Critique

The Zheng comparisons use only 16 media. Point correlations at 0.5× and 1× are
not enough to support a superiority claim, especially because the same media
are used for every arm and the three sequencing seeds are not independent.

### Resolution

Added a 10,000-resample medium-level bootstrap. Within each resampled medium,
the three seed estimates are averaged before the metric is computed, so seeds
are not treated as independent biological observations.

Key all-finite A-arm intervals are:

| depth | Pearson r | slope |
|---:|---:|---:|
| 0.5× | 0.911 (0.820–0.966) | 0.565 (0.438–0.686) |
| 1× | 0.960 (0.932–0.980) | 0.850 (0.679–0.986) |

For the paired A-minus-E comparison, the Pearson delta is 0.059
(-0.031 to 0.114) at 0.5× and 0.037 (0.010 to 0.057) at 1×. Thus the paper now
claims better low-depth recall/rank recovery, not a uniform significant accuracy
win. The bootstrap p-values are exploratory and have no multiplicity correction.

Files: `benches/refresh_20260930/results/zheng_bootstrap_media.tsv`,
`benches/refresh_20260930/results/zheng_arm_delta_bootstrap.tsv`.

## Reviewer 2 — “mean QC-passed” can hide instability

### Critique

Averaging pass/fail over three seeds may make a non-deterministic QC gate look
stable. Deployment users need to know whether a condition passes reproducibly.

### Resolution

Added exact per-seed, all-seed and any-seed pass counts. For arm A, the counts
are 8/8/5 at 2×, 14/11/13 at 5× and 14/15/13 at 10×. Only 2, 8 and 10 media
pass in all three seeds. The manuscript wording now distinguishes the mean pass
count from reproducible deployment.

File: `benches/refresh_20260930/results/zheng_qc_pass_counts.tsv`.

## Reviewer 3 — EM appears overclaimed

### Critique

The pipeline computes EM weights, but the ZTP model consumes raw integer counts
and excludes masked/shared anchors. Describing PTR as EM-reassigned shared-
anchor coverage would overstate the implementation.

### Resolution

Scoped the claim. `<sample>.em.tsv` is an interface for abundance, containment
and downstream assignment tools. The production PTR model intentionally uses raw
integer observations and excludes database-shared anchors. The single-isolate
Zheng benchmark has no cross-genome shared anchors.

Updated: `docs/PAPER_RESULTS.md`, `docs/PAPER_CLAIMS.md`,
`docs/design/03-data-formats.md`.

## Reviewer 4 — the old A4 floor may be stale after signed fixed-origin slopes

### Critique

The historical b=0 floor and compression estimates predate the signed
fixed-origin fit and cannot support current mechanistic claims.

### Resolution

Re-ran the full A4 grid. The b=0 mean is now at most 0.10 in every arm/depth:
Poisson 0.099/0.055/0.035/0.026 and NB
0.033/0.050/0.084/0.036 at 0.5/1/2/5×. Bootstrap CIs agree that the former
floor is no longer the dominant defect. Poisson recovered slopes are close to
one by 1×, while NB remains compressed at ≤1× and has upward bias in low-rate
windows.

Files: `benches/refresh_20260930/results/a4_slope_bootstrap.tsv`,
`benches/refresh_20260930/a4/results/`.

## Reviewer 5 — one stationary false fire undermines WCG as a production gate

### Critique

The current R3 grid contains one false fire (`scr20RX`, 1×) even though the
truth is stationary. A method that flags a non-growing culture as fragmented is
not ready for production QC.

### Resolution

Kept the primary result and added threshold sensitivity. With the original
`z_short > 2` cut, there is 1 stationary false fire among 81 stationary runs
and 96 total flags. Raising the slope-branch cut to `z_short > 2.5` removes that
false fire and retains 72/72 strong scrambled detections at 5–10×; total flags
fall to 93.

This is a post-hoc sensitivity result, not independent validation. The docs now
recommend keeping WCG disabled until an independent null/positive evaluation.

File: `benches/refresh_20260930/results/r3_threshold_sensitivity.tsv`.

## Reviewer 6 — three subsamples are not biological replication

### Critique

The s1/s2 resamples derive from the C1 20× pool. They quantify technical
subsampling variability, not biological or laboratory replication, and cannot
support metagenomic deployment claims.

### Resolution

The methods and wording gates state this explicitly. The refresh supports
current-implementation reproducibility only. Claims about real communities,
bioinformatics pipeline variants, independent biological replicates and
cross-strain references remain out of scope.

## Reviewer 7 — historical and current numbers were at risk of being mixed

### Critique

The repository contains 2026-08-24 results, pre-fix mechanism records and a
2026-09-30 refresh. Without a provenance gate, manuscript text could quote
stale numbers.

### Resolution

Added `docs/PAPER_RESULTS.md` and `docs/PAPER_CLAIMS.md`, marked old benchmark
sections as historical/provenance, and recorded source/result archive hashes.
Current claims must cite the refresh; old numbers must not be mixed with it.

## Remaining open issues

1. Independent biological replicates and real communities.
2. Independent validation of WCG thresholds.
3. A metagenomic PTR design that can use EM-assigned fractional counts without
   violating the integer count model.
4. A multi-tool benchmark beyond Pilea.
