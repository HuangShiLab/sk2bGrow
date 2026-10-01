# Manuscript claim corrections

Use [`PAPER_RESULTS.md`](PAPER_RESULTS.md) for the 2026-09-30 current-commit
numbers and archive hashes. The rules below are wording gates.

## 1. Low coverage: signal is not the same as deployability

**Do not write:** “sk2bGrow is usable at 0.5–1× and passes its gate.”

**Write instead:** in the three-subsample Zheng refresh, the all-finite
estimator view has mean per-medium Pearson r=0.911 at 0.5× and r=0.960 at 1×
for anchors + V-fit. The 0.5× seed range is wide (r=0.716–0.899). Default QC
nevertheless passes a mean of 0/16 media at both depths and 7/16 at 2×. Thus the
coordinate estimator has low-coverage signal, but shipped QC does not call
those estimates deployable.

Pass counts vary by seed: A-arm counts are 8/8/5 at 2×, 14/11/13 at 5× and
14/15/13 at 10×. Only 2, 8 and 10 media pass in all three seeds. Do not report
the mean pass count as if deployment were deterministic.

Under the same all-finite rule, Pilea gates-off is competitive: FracMinHash +
V-fit reaches r=0.852 at 0.5× and r=0.923 at 1×. The principal low-depth
difference is recall/reporting, not a uniform accuracy win.

Do not convert these contrasts into a significant-superiority claim. The media
bootstrap A-minus-E Pearson delta is 0.059 (CI -0.031 to 0.114) at 0.5× and
0.037 (CI 0.010 to 0.057) at 1×; the 0.5× slope delta favours FracMinHash +
V-fit.

Report both views:

1. all finite estimates — estimator performance;
2. default-QC-passed estimates — deployed performance.

Never compare the all-finite sk2bGrow column with Pilea defaults as though the
reporting rules were identical.

## 2. Enzyme panel independence

**Do not write:** “16 independent enzymes” or “16 independent measurements.”

**Write instead:** 16 deterministic enzyme strata. Every Bsp24I tag is also a
CjePI tag, so the panel contains at most 15 independent enzyme channels. Cochran
Q p-values should be described as a heterogeneity diagnostic, not an exact test
of fully independent strata.

## 3. Fragmented references and coordinate QC

Use the 2026-09-30 R3 grid. Complete and correctly ordered controls are quiet
(0/36 and 0/72). Scrambled references fire in 15/54 cells at true
log2PTR=0.5 and 40/54 at each of log2PTR=1 and 1.5; detection is complete
(72/72) at 5–10× for true log2PTR≥1.

However, one scrambled stationary cell fires at 1×, and 1× sensitivity is weak.
WCG is therefore prototype evidence, not a deployed production gate.

A post-hoc sensitivity check suggests `z_short > 2.5` removes that false fire
while retaining all strong 5–10× detections. This threshold must be validated
on independent null/positive data before implementation; do not present it as
a pre-specified result.

Keep the architectural statement: on unscaffolded contigs, a position-free
estimator is structurally advantaged. `scaffold` now emits an indexable
single-contig FASTA, but this does not make coordinate QC production-ready.

## 4. EM is not a PTR shortcut

**Do not write:** “shared anchors are EM-reassigned and then used for PTR.”

**Write instead:** the production ZTP PTR model uses raw integer anchor counts
and conservatively excludes database-shared anchors. EM output is retained in
`<sample>.em.tsv` for abundance, containment and downstream assignment tools.
In the single-isolate Zheng benchmark there are no cross-genome shared anchors.

## 5. Low-coverage compression mechanism

The current A4 refresh supports a mechanism claim, not merely “ZTP causes it”:

- exact window rates through the same fit/fusion machinery recover slopes near
  one;
- after the fixed-origin signed-slope fix, the b=0 floor is small (≤0.10 in
  every arm/depth) and is no longer the dominant low-depth defect;
- recovered slopes are 0.774/0.958/1.012/0.995 for Poisson at
  0.5/1/2/5× and 0.701/0.853/0.857/0.982 for NB;
- residual NB compression at ≤1× is accompanied by upward bias in low-rate
  windows and should be described as a mis-specification/noise effect;
- the all-ones ZTP boundary and ZTNB ridge defects remain fixed, but those
  fixes should not be credited with eliminating all NB compression.

## 6. Ground truth

Use `growth_rate` (independent optical density/microscopy measurement) for
accuracy/correlation. `pred_log2ptr` is a within-reads marker-frequency target
and must be labelled non-independent. Never describe RMSE against
`pred_log2ptr` as absolute accuracy.

## 7. Minimum reporting standard

Every numeric claim should include:

- exact source/archive hash and binary/build profile;
- Pilea version;
- genome/reference accession;
- read count and subsampling scheme;
- estimator (`--method`) and windowing;
- QC view (`all finite`, `default QC passed`, or another explicitly defined set);
- n estimates and n QC-passed;
- uncertainty or seed-to-seed variability;
- negative-control result.

Three sequencing subsamples are sufficient for a current-implementation check,
but not for a biological-replicate claim. The Zheng refresh has no independent
biological replication and no community complexity. Low-coverage deployment
claims require both repeated subsamples and independent biological replicates.
