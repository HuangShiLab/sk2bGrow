# Reference fragmentation (A2 / experiment C2)

The V-shape fit regresses replication signal on **genomic coordinate**. A MAG
has no coordinate: contig order and orientation are unknown. This is the one
place where Pilea is architecturally better — sorted-rank regression discards
position by construction and is fragmentation-proof for free.

So the question is not rhetorical. Either `sk2bgrow scaffold` recovers enough
order that the fit degrades gracefully, or the honest framing of the whole method
is *"for complete references"* and every MAG-based dataset becomes a limitation
section.

## Design

Pilea's Fig 3 protocol: the same genome, cut into 100 contigs with lognormal
lengths (μ = 0, σ = 1), shuffled, each independently reverse-complemented with
probability 0.5. Sequence content is untouched — the 100-contig database holds
43,707 anchors against the complete genome's 43,735, the 28 lost being tags that
straddled a cut — so any change in accuracy is attributable to the lost
coordinate and nothing else.

Four reference conditions over identical reads (Zheng *E. coli*, 16 media plus
the stationary control, subsampled to 0.5–10×):

| condition | reference |
|---|---|
| `complete` | the finished chromosome — upper bound |
| `frag` | 100 contigs — the MAG case |
| `scafSelf` | `frag` scaffolded against the same genome — upper bound on scaffolding |
| `scafRel` | `frag` scaffolded against *E. coli* O157:H7 — the realistic MAG case |

plus Pilea gates-off on `frag`, to check rather than assume that rank regression
is indifferent to fragmentation.

`scafSelf` is circular by construction and is not a result; it exists to separate
"scaffolding cannot work" from "scaffolding cannot work *across strains*".

## Scaffolding accuracy, scored against the known layout

| | vs itself | vs O157:H7 |
|---|---|---|
| contigs placed | 100/100 | 99/100 |
| orientation correct | 100/100 | 99/99 |
| median start error | 0 bp | 712,324 bp |
| after removing a rigid rotation | 0 bp | 73,276 bp |
| contig order (Spearman) | **1.0000** | **1.0000** |

The raw 712 kb against O157:H7 reads worse than it is. Almost all of it is a
single rigid rotation — O157:H7's origin sits ~786 kb away in the coordinate
frame — and the estimator searches for the origin rather than assuming it, so a
rotation is invisible to the fit. What matters is relative order, and that is
recovered *exactly*: Spearman 1.0000 against the truth, across a different
strain. The residual 73 kb median is local stretch from O157:H7's ~865 kb of
strain-specific insertions, which is 1.6% of the chromosome.

**Caveat in the reconstruction, not the scaffolding.** `rescaffold.py` writes each
contig into an N canvas at its inferred start, so where two contigs are placed
overlapping — which happens when the coordinates come from a different strain —
the later write wins and that sequence is lost. Against O157:H7 this drops
89,597 bp, 1.93% of the draft (760 of 43,707 anchors); against itself it is
728 bp. The results below are therefore a slight *under*-estimate of what
scaffolding delivers. Emitting contigs at cumulative order-preserving offsets, or
teaching `index` to read a scaffolded TGT directly, would remove it.

## Result

Accuracy against the independently measured growth rate, n = 16 media per cell:

| coverage | | complete | 100 contigs | scaffolded vs O157:H7 | Pilea on 100 contigs |
|---|---|---:|---:|---:|---:|
| 1× | Pearson r | 0.981 | 0.550 | 0.977 | 0.827 |
| | RMSE | 0.157 | 0.890 | 0.152 | 0.545 |
| | slope | 0.779 | 0.103 | 0.755 | 0.628 |
| 10× | Pearson r | 0.968 | 0.859 | 0.967 | 0.960 |
| | RMSE | 0.063 | 0.862 | 0.086 | 0.077 |
| | slope | 0.951 | 0.210 | 0.984 | 0.820 |
| | QC pass | 75% | **100%** | 88% | — |

**Fragmentation does not add noise; it removes the gradient.** At 10× every
fragmented estimate lands at roughly a fifth of its true value — 0.35 where the
prediction is 1.73, 0.10 where it is 0.55 — for a fitted slope of 0.21 against
0.95 on the complete chromosome. Correlation flatters this badly: r = 0.86 at
10× looks survivable, and it is not. Read the slope.

**Scaffolding restores it completely**, and against a *different strain*, which
is the case that matters. Every scaffolded column matches the complete-genome
column to within about 0.04 at every depth, and the RUN_OUT stationary control
stays near zero throughout (0.05–0.30) in all conditions.

**The QC does not catch the failure — it is anti-correlated with it.** 100% of
the fragmented estimates pass at 5–10×, against 75% for the correct ones. That
is not a bug in the thresholds but a blind spot in what they test: the fusion QC
asks whether the enzymes *agree*, and a destroyed coordinate makes all sixteen
agree that there is no gradient. Cochran's Q cannot see a failure that is
identical across strata.

**On unscaffolded contigs, Pilea wins outright.** Its rank regression discards
position by construction, so fragmentation barely touches it: r falls from 0.889
to 0.827 at 1× and from 0.971 to 0.960 at 10×. Against sk2bGrow's unscaffolded
0.550 and 0.859 that is not close. This is a real architectural advantage of the
sorted estimator and should be stated as one. What reverses it is scaffolding,
not the sketch — and scaffolding is available precisely because the anchors are
deterministic and shared between genomes.

## How fragmented is too fragmented?

`sweep.sh`, at 10×, n = 16 media per row:

| contigs | N50 | Pearson r | RMSE | bias | slope | QC pass |
|---:|---:|---:|---:|---:|---:|---:|
| 1 (complete) | 4.64 Mb | 0.968 | 0.063 | +0.031 | 0.951 | 75% |
| 2 | 2.62 Mb | 0.970 | 0.080 | −0.064 | 0.876 | 94% |
| 5 | 939 kb | 0.971 | 0.308 | −0.285 | 0.634 | 100% |
| 10 | 626 kb | 0.948 | 0.371 | −0.349 | 0.609 | 100% |
| 20 | 305 kb | 0.961 | 0.445 | −0.412 | 0.508 | 100% |
| 50 | 156 kb | 0.962 | 0.609 | −0.574 | 0.444 | 94% |
| 100 | 78 kb | 0.859 | 0.862 | −0.808 | 0.210 | 100% |

**There is no safe threshold above two contigs, and correlation will not tell you
where you are.** r never leaves 0.86–0.97 across the whole sweep while the slope
falls from 0.88 to 0.21 and the bias grows to −0.81. At 50 contigs — an N50 of
156 kb, a draft most people would call good — r reads 0.96 and every estimate is
44% of truth.

The degradation is smooth and monotone in bias, so this is not a cliff to stay
above; it is a tax that starts as soon as the reference is not closed. The
practical rule is simply: **scaffold anything that is not a single contig**, and
report the slope beside r.

A contig-count check is therefore the right guard, and it should fire at 2, not
at some tuned larger number.

## Can Pilea's estimator be borrowed instead?

Partly, and the useful half is the *formulation*, not the code.

Under the standard model log2(coverage) is a tent function of position with equal
|slope| on both replichores. Over a uniformly tiled genome the coverage *values*
are therefore uniform on [log2 c_ter, log2 c_ori], and the width of that uniform
**is** log2(PTR). Only the spread of coverage is needed; position never enters.
That is exactly why sorting the windows costs Pilea nothing here.

What should not be borrowed is how Pilea estimates the width. It sorts the
per-window rates and RANSAC-regresses them on rank, taking the fitted rise as
log2(PTR). Sorted values of *any* sample rise, so sampling noise is counted as
growth — on the stationary control it reports 1.15 at 1× for a culture that is
not growing.

`spread_estimator.py` keeps the formulation and puts the noise in the model.
Each window already carries a standard error, so with mu_w ~ U(a, a+W) and
e_w ~ N(0, s_w²), maximising the marginal likelihood over (a, W) per enzyme and
fusing across enzymes gives an order-free estimate that cannot manufacture a
gradient. `score_spread.py` scores it from the `windows.rates.tsv` the runs
already wrote — no re-run needed.

| coverage | estimator, on the 100-contig reference | RMSE | slope |
|---|---|---:|---:|
| 5× | V-fit (coordinate destroyed) | 0.871 | 0.187 |
| | **spread-MLE** | **0.144** | **0.922** |
| | Pilea | 0.116 | 0.888 |
| 10× | V-fit | 0.862 | 0.210 |
| | **spread-MLE** | **0.227** | **0.919** |
| | Pilea | 0.077 | 0.820 |

| stationary control (truth ≈ 0) | 1× | 2× | 5× | 10× |
|---|---:|---:|---:|---:|
| Pilea | 1.153 | 0.621 | 0.255 | 0.203 |
| spread-MLE | **0.000** | **0.000** | 0.383 | 0.458 |

So it recovers most of what fragmentation destroys — RMSE 0.87 → 0.14 at 5×,
slope 0.19 → 0.92 — with no scaffolding at all, and it does not invent growth
where there is none. Three things it does not do, stated plainly:

1. **It does not beat Pilea at depth on a fragmented reference** (RMSE 0.227
   against 0.077 at 10×). Pilea's advantage there is real.
2. **Below 5× it over-shrinks to zero** (bias −0.98 at 1×). The distribution
   alone carries too little information. Coordinates are not redundant — which is
   why the V-fit wins at 1–2× and why scaffolding stays the better answer in the
   band this project is actually about.
3. **It has a positive floor** (~0.4 on the stationary control at 5–10×). Real
   window scatter includes systematic terms — mappability, anchor density,
   residual GC — that a pure noise model can only attribute to growth. An
   overdispersion parameter is the obvious next step and should fix both this and
   (2).

Pooling all sixteen enzymes into one width was tried and is much worse (8.6 on
the stationary control at 0.5×): the per-enzyme fusion earns its place by
down-weighting sparse enzymes.

**The design this points to is two estimators, not one.** With a coordinate —
complete reference, or a scaffolded draft — use the V-fit, which is the only
thing that works at 1–2×. Without one, and at ≥5×, use the spread estimator.
Without one and below 5×, report nothing: that is the case the current QC passes
at 100% while being wrong by a factor of five.

## Running

```bash
WORK=/scratch/frag SUB=../zheng2020/sub GENOME=../ecoli.fna \
  RELDB=/path/to/db_with_relatives ./run.sh
WORK=/scratch/frag python3 analyze.py
```

`RELDB` is any sk2bGrow database containing the scaffolding reference;
`SELF_NAME` and `REL_NAME` select which genome inside it to use.

## Files

| file | what |
|---|---|
| `spread_estimator.py` | the order-free spread estimator (prototype) |
| `score_spread.py` | scores it against the V-fit and Pilea from existing runs |
| `fragment.py` | cut a complete genome into lognormal contigs, shuffled and flipped; writes a `.layout.tsv` truth file |
| `rescaffold.py` | score a scaffold result against that truth, and re-emit the draft as a coordinate-bearing FASTA |
| `run.sh` | the four reference conditions plus the Pilea control |
| `analyze.py` | accuracy per coverage per condition |
