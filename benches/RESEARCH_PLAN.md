# sk2bGrow — research plan for the HPC phase

For whoever runs the large-scale benchmarks next. Everything here was measured on
an M3 Max laptop against reduced datasets in `benches/`. Measured numbers are
given as measured; extrapolations are labelled as extrapolations. This file
supersedes the earlier `HPC.md` handover and absorbs it.

**Read Part 3 before Part 6.** Two claims in this project have already been
overturned by adding a missing control, and one documented mechanism turned out
to be simply wrong. The most useful thing you can do is keep breaking things that
way.

---

## Part 0 — Orientation

| where | what |
|---|---|
| `crates/sk2bgrow-core` | Rust: digest, anchor DB, counting, windows, origin search |
| `python/sk2bgrow` | Python: ZTP/ZTNB window rates, GC correction, V-fit, fusion, QC, report |
| `benches/zheng2020` | *E. coli* isolate benchmark, all arms, ground truth |
| `benches/fragmentation` | reference-fragmentation benchmark + the order-free prototype |
| `../sk2bGrow-paper/data/*.tsv` | every committed result table |
| `../sk2bGrow-paper/figures/` | `make_figures.py`, `make_tables.py` — pure functions of `data/` |

The integration contract is deliberately narrow: **the figure and table scripts
read only `data/*.tsv`, use no network, and regenerate everything**. Dropping
refreshed TSVs into `data/` with the same column names is the whole handoff. PDF
timestamps are pinned, so a regeneration with no change produces no git diff —
use that to check whether a figure really moved.

```bash
cargo test                                   # 86 tests, keep green
cargo build --release
COVS="0.5 1 2 5 10" ./benches/zheng2020/bench.sh
```

---

## Part 1 — The algorithm

### 1.1 What it estimates

Peak-to-trough ratio. During replication a growing cell has more copies of
sequence near the origin than near the terminus, so read coverage carries a
gradient whose amplitude is the growth rate:

```
log2(PTR) = C / tau = lambda * C / ln 2
```

where *C* is the replication period. Over a uniformly tiled chromosome, log2
coverage is a **tent function** of position — linear from ori down to ter with
equal |slope| on both replichores.

Two consequences are used repeatedly below, and it is worth stating them
explicitly because most of the design follows from them:

1. **The tent's amplitude is log2(PTR).** Fitting the tent in coordinate space is
   the direct estimator.
2. **The tent's *values* are uniformly distributed** on [log2 c_ter, log2 c_ori],
   with width log2(PTR). So PTR is also recoverable from the coverage
   *distribution alone*, with no coordinate. This is why Pilea's sorted-rank
   estimator is indifferent to reference fragmentation, and it is the basis of
   the prototype in `benches/fragmentation/spread_estimator.py`.

### 1.2 Why a deterministic sketch

sk2bGrow counts **Type IIB restriction-enzyme tags** (2bRAD anchors) rather than
a random hash subset. A 16-enzyme panel on *E. coli* yields 42,234 usable anchors
after masking multi-copy and shared sites — mean spacing 110 bp, median 59 bp.
Density varies ~20-fold across the panel (CjeI 1,962/Mb, PpiI 74/Mb).

The point is not that anchors are denser than a FracMinHash sketch. It is that
they are **at known coordinates and are the same in every sample and every
genome**, which buys three things a random sketch cannot: a coordinate for the
tent fit, sixteen *independent strata* for a consistency test, and shared tags
between genomes, which is what makes `sk2bgrow scaffold` possible.

### 1.3 The pipeline

**Stage 1 — index (offline, once).** Digest each reference with the panel; store
each tag with its coordinate, strand, local GC (±250 bp), and flags for
multi-copy / shared-across-genomes. Anchors are held in CSR arrays (sorted keys,
offsets, values) — see Part 3.3 A1.

**Stage 2 — count.** Scan reads for panel tags with pigeonhole seed indexing,
tolerating `--max-mismatch` (default 2). Scanning is **grouped by tag length, not
by enzyme**: two enzymes of the same tag length can satisfy their patterns on the
same physical read window, and the panel guarantees this happens because every
Bsp24I tag is byte-identical to a CjePI tag. One physical tag is one observation,
credited to each stratum the locus belongs to. Getting this wrong inflated
low-depth estimates by a large margin — see Part 3.2.

**Stage 3 — windows and rates.** Windows are cut inside each (genome, enzyme)
series in coordinate order, holding a fixed *number of anchors*, never straddling
a contig boundary. Size per enzyme is `clamp(n_anchors / 25, 25, 100)`, targeting
~25 windows. Counts within a window are fitted as a **zero-truncated Poisson
mixture** (EM, component count by BIC), with a **zero-truncated negative
binomial** branch selected by BIC when overdispersed. Truncation matters because
sequence divergence makes some reference anchors genuinely absent.

**Every window rate is returned with a standard error.** This is the single most
load-bearing output of the stage: it is what makes inverse-variance fusion
possible, and what lets the order-free estimator subtract noise instead of
reading it as growth.

**Stage 4 — GC correction.** A loess curve of log2 efficiency against GC fitted
**within each enzyme** and applied **at anchor resolution**, as an offset in log2
space. Each Type IIB site has its own base composition, so a single global curve
fits none of the sixteen. A constant per-enzyme efficiency cannot bias PTR (it is
absorbed by that enzyme's intercept); only GC *slope within* an enzyme matters.

**Stage 5 — origin and V-fit.** The origin is a property of the chromosome, so it
is estimated **once per genome from all enzymes pooled**, after median-centring
each enzyme. Then log2 coverage is regressed on circular distance *d* from the
origin,

```
log2 mu(x) = a - b1 * min(d, k) - b2 * max(0, d - k)
```

a piecewise-linear V on *real coordinates*, with an optional kink (segmented form
offered only with >=30 windows, selected by BIC). An externally supplied origin
(DoriC / Ori-Finder / *dnaA*) is used as a prior when given.

**Stage 6 — fusion and QC.** Each enzyme is an independent stratum. Estimates are
combined by inverse-variance weighting; **Cochran's Q** tests the null that all
enzymes measure one common value; on rejection the estimator escalates to
DerSimonian–Laird random effects. QC gates on coverage, detected fraction,
dispersion, containment, enzyme consistency, enzyme fit rate, and origin
confidence.

---

## Part 2 — Analysis method

### 2.1 The arms

sk2bGrow differs from Pilea in **two** places — what is counted, and how the
gradient is fitted. Comparing the two pipelines end-to-end cannot attribute the
difference, so the two factors are crossed:

| | coordinate V-fit | sorted-rank regression |
|---|---|---|
| **2bRAD anchors** | sk2bGrow (arm A) | arm B |
| **FracMinHash** | arm E | Pilea, gates off (arm C) |

Arm B runs our anchors through Pilea's windowing and estimator. Arm E rewrites
Pilea's own sketch into our count-table format so the **unmodified** estimator
runs on it (`benches/zheng2020/armE_counts.py`). Pilea is always reported in
**both** configurations — shipped defaults, and `-x 0 -z 0 -c 0` — because at its
defaults it returns nothing below 10×, and a benchmark of "estimate vs no
estimate" is not a comparison of accuracy, while reporting only the gates-off arm
is misleading in the other direction.

**Do not build arm E by feeding raw window rates into the V-fit.** That returns
log2PTR 2.65 against a measured 1.73, because what it removes is the surrounding
machinery, not the fitting geometry.

### 2.2 Metrics, and what each one misses

| metric | detects | misses |
|---|---|---|
| Pearson r | whether the ranking is recovered | a compressed or inflated scale — r is invariant to affine transformation |
| Spearman rho | monotone recovery, robust to outliers | magnitude, as above |
| RMSE | magnitude error in log2PTR units | direction |
| bias | systematic over/under-estimation | scatter |
| **slope** | dynamic-range compression | offset |
| recall | silence — a tool that reports nothing scores no error | whether what *was* reported is right |
| L2 (Pilea's metric) | — | **it rewards silence**: fewer terms, smaller norm |
| spurious | false positives | — |

Three rules, each of which this project learned the hard way:

1. **Recall and error must be read together.** A tool reporting 2 of 16 strains
   accurately scores better on L2 than one reporting all 16 with small errors.
   Report recall beside RMSE everywhere, and prefer per-genome RMSE to per-sample
   L2.
2. **Never quote r alone.** The fragmentation experiment is the extreme case:
   across 1 to 100 contigs r never leaves 0.86–0.97 while the slope falls from
   0.95 to 0.21. At 50 contigs r reads 0.96 and every estimate is 44% of truth.
3. **A constant output has no correlation.** At 0.5× Pilea's gates-off arm
   returns PTR = 1.0 for every sample; report r as *absent*, not zero, or it
   reads as "scored badly" instead of "produced no signal".

**Two different slopes appear in our outputs and they are not interchangeable.**
`table2` reports OLS of estimate on measured λ (0.615 at 0.5×, 0.951 at 10×);
the A4 analysis reports OLS of estimate on the λ·C-derived *predicted* log2PTR
(0.681 at 0.5×, 1.074 at 10×). Label which one you mean every time. RMSE and bias
are always against the predicted log2PTR; correlation always against λ.

**RMSE is against a non-independent target.** The predicted log2PTR comes from
the source paper's own marker-frequency analysis of the same reads. Say so
wherever it is quoted.

### 2.3 Cost

`/usr/bin/time -l` (macOS) or `-v` (Linux). On macOS `ru_maxrss` does propagate
from a grandchild, so the Python layer is included. Report wall time **per
depth**, never averaged: Pilea's cost is non-monotonic (5.6 s at 0.5×, **21.1 s
at 2×**, 5.6 s at 10×) because its ZTP-mixture EM is slowest where the mixture is
least identifiable — exactly the depth band this project is about.

---

## Part 3 — What is established

All on the Zheng *E. coli* dataset (16 media, λ = 0.40–1.72 h⁻¹, independent
ground truth) unless stated. n = 16 per cell.

### 3.1 Headline results

**Accuracy, complete reference** (Pearson r / RMSE vs predicted):

| depth | anchors + V-fit | FracMinHash + V-fit | anchors + rank | Pilea, gates off | Pilea, defaults |
|---|---|---|---|---|---|
| 0.5× | **0.913** / 0.304 | 0.724 / 0.382 | 0.164 / 1.188 | degenerate | no estimate |
| 1× | **0.981** / 0.157 | 0.940 / 0.213 | 0.683 / 1.027 | 0.888 / 0.397 | no estimate |
| 2× | 0.982 / 0.128 | 0.984 / 0.093 | 0.756 / 0.667 | 0.947 / 0.259 | no estimate |
| 5× | 0.979 / **0.039** | 0.977 / 0.051 | 0.914 / 0.346 | 0.953 / 0.103 | no estimate |
| 10× | 0.968 / 0.063 | 0.942 / 0.149 | 0.913 / 0.230 | 0.971 / 0.077 | 0.971 / 0.077 |

**The gain is an interaction, not either factor alone.** At 1× the coordinate fit
is worth +0.30 r on anchors but only +0.05 on a FracMinHash sketch, and the
sketch effect *changes sign* with the estimator (+0.04 under the V-fit, −0.21
under rank regression). At 0.5× only the combination works.

**Which Pilea gate binds.** Re-applying each default threshold to gates-off
output (exact, no extra runs): `--min-cove` is responsible for **68 of 68**
suppressed *E. coli* estimates and **167 of 172** in simulation. A 150 bp read
gives 120 31-mers, so its k-mer threshold of 5 is ≈6.6× read coverage.

**Negative control** (stationary culture, truth ≈ 0):

| arm | 0.5× | 1× | 2× | 5× | 10× |
|---|---:|---:|---:|---:|---:|
| sk2bGrow | 0.260 | 0.077 | 0.097 | 0.060 | 0.046 |
| anchors + rank | 2.212 | 1.816 | 1.421 | 0.843 | 0.679 |
| Pilea, gates off | 0.000 | 0.963 | 0.586 | 0.244 | 0.168 |

Sorted regression reports **log2PTR 2.21 for a culture that is not growing**. It
manufactures a gradient out of rank-ordered noise. This is the clearest single
demonstration in the benchmark of why coordinates matter.

**Multi-strain simulation** (16 genomes, 4/8/16 strains × 1/2/4/8×):

| | recall | RMSE | bias |
|---|---:|---:|---:|
| sk2bGrow | **1.000** | **0.134** | **−0.013** |
| Pilea, gates off | 0.997 | 0.265 | +0.168 |
| Pilea, defaults | 0.224 | 0.083 | −0.045 |

Pilea's defaults have the best RMSE *because they report 22% of the genomes*.
This is metric rule 1 in one line.

**Fragmented references** (100 shuffled lognormal contigs; 43,707 of 43,735
anchors survive, so only the coordinate changes):

| depth | | complete | 100 contigs | scaffolded vs O157:H7 | Pilea on contigs |
|---|---|---:|---:|---:|---:|
| 1× | r | 0.981 | 0.550 | 0.977 | 0.827 |
| | slope | 0.779 | 0.103 | 0.755 | 0.628 |
| 10× | RMSE | 0.063 | 0.862 | 0.086 | 0.077 |
| | slope | 0.951 | 0.210 | 0.984 | 0.820 |
| | QC pass | 75% | **100%** | 88% | — |

- Fragmentation **removes the gradient** rather than adding scatter: every
  estimate lands at about a fifth of truth.
- `sk2bgrow scaffold` restores the complete-reference result **exactly, across
  strains** — contig order returns with Spearman 1.0000 against O157:H7,
  orientation 99/99 correct. The 712 kb raw placement error is a rigid rotation,
  invisible to a fit that searches for the origin.
- **On unscaffolded contigs Pilea wins outright** (r 0.827 vs 0.550 at 1×). Say
  so; it is a real architectural advantage of the sorted estimator.
- **The QC is anti-correlated with the failure**: 100% of fragmented estimates
  pass at 5–10× against 75% of correct ones. Cochran's Q asks whether the enzymes
  *agree*, and a destroyed coordinate makes all sixteen agree there is no
  gradient. Q cannot see a failure identical across strata.
- There is **no safe contig count**. Degradation is smooth and monotone from two
  contigs upward; even 5 contigs (N50 939 kb) costs 0.29 in bias.

**Panel size** (averaged over ≤2×):

| enzymes | anchors | r | RMSE | run-out bias | s/sample |
|---|---:|---:|---:|---:|---:|
| 2 | 17,055 | 0.951 | 0.188 | 0.073 | 3.2 |
| 4 | 24,753 | 0.959 | **0.165** | 0.074 | 4.4 |
| **8** | 37,232 | **0.969** | 0.171 | 0.123 | 6.5 |
| 16 | 43,735 | 0.959 | 0.196 | 0.145 | 8.2 |

**Cost**, seconds per sample, 85 *E. coli* cells. Reported per depth because
Pilea's is non-monotonic and a single averaged number hides that:

| depth | k=2 | k=8 | k=16 | Pilea, gates off | Pilea, defaults |
|---|---:|---:|---:|---:|---:|
| 0.5× | 2.2 | 3.8 | 4.3 | 5.6 | 0.9 |
| 1× | 3.1 | 5.9 | 6.8 | **16.8** | 1.0 |
| 2× | 3.4 | 6.8 | 8.3 | **21.0** | 1.0 |
| 5× | 3.6 | 7.5 | 9.8 | 8.6 | 1.1 |
| 10× | 4.0 | 8.4 | 11.9 | 5.6 | 5.5 |

sk2bGrow is monotone in depth and linear in panel size. Pilea's gates-off arm
peaks at 2× — the depth band this project is about — because its ZTP-mixture EM
is slowest where the mixture is least identifiable. Its default arm is cheap
below 10× only because it does no fitting there and returns nothing.
Peak RSS: 162 MB at k=2, 187 at k=8, 195 at k=16; Pilea 152–160 MB.

**Memory.** CSR index: 352 → **160 B/anchor** marginal at `--max-mismatch 2`,
count tables bit-identical, ~1% slower. GTDB species-rep projection 1.65 TB →
**752 GB** at 16 enzymes, ~672 GB at 8.

### 3.2 Claims that were tested and **refuted**

This section exists because it is the most useful part of the handover.

| claim | status |
|---|---|
| "The gain is the estimator; deterministic anchors are if anything behind" | **wrong** — an artefact of a missing 2×2 cell. The factors interact |
| "Pilea is ~7× faster" | **withdrawn** — at 1–2× it is 2.5× *slower* (16.8 s and 21.0 s against our 6.8 and 8.3); it is cheaper only at 0.5× and 10× |
| "Sparse enzymes are more biased, which explains the panel-size result" | **falsified** — −0.319 vs −0.312, no difference |
| "Fixed vs random effects explains the panel-size result" | **falsified** — RMSE 0.581 vs 0.506, wrong direction |
| "CSR will also be faster (cache-friendly, no pointer chase)" | **wrong** — ~1% slower; it buys memory only |
| "CSR will cost 12–16 B/anchor" | **wrong by an order of magnitude** — that costed only the steady-state arrays |
| "Low-coverage compression is errors-in-variables attenuation" | **wrong** — the predictor in that regression is effectively noise-free, and the slope exceeds 1 at depth, which no attenuation story predicts |
| "Origin misplacement explains the compression" | **mostly no** — a perfect origin recovers 0.06 of a 0.32 shortfall |

### 3.3 Open-issue status

| id | issue | status |
|---|---|---|
| A1 | index memory | **done** — 2.2×, 752 GB at GTDB scale. Remaining halvings: drop the build-time pair buffer (~58 B/anchor), mmap the AnchorDb (~38), ship 8 enzymes, shard by genome |
| A2 | fragmented references | **done** — scaffolding rescues it; **the QC blind spot it exposed is now the top open item** |
| A3 | attribution 2×2 | **done** — changed the claim |
| A4 | low-coverage compression | **open, narrowed**. Excluded: outlier trimming (100% of windows kept), IV-fusion shrinkage (estimate and SE are *negatively* correlated, −0.06 to −0.28; IV beats unweighted 0.681 vs 0.497), origin error (minor). **Remaining candidate: what the ZTP/ZTNB layer does to a window's rate at near-zero counts, before the fit sees it** |
| A5 | enzyme containment (Bsp24I ⊂ CjePI) | **measured, immaterial** — max 0.050 log2 change, Q rejection does not improve. Document; do not code |
| A6 | ship 8 enzymes not 16 | analysis done, **needs the panel sweep crossed with fragmentation** before the default changes. Recommendation stands: ship 8 |
| A7 | origin annotation | **tested** — ~10% RMSE at ≤1× (0.304→0.271, 0.157→0.143), nothing at ≥5×. Not the fix for A4 |
| A8 | strain mixtures | **out of scope by design** — two tents superposed are not identifiable from one sample. State as a limitation |

Also worth carrying: at 0.5× the mean per-enzyme V-fit has **r² = −0.009** —
individually worthless — while fusing sixteen of them correlates 0.913 with
measured growth rate. The panel, not any single fit, is what works there.

---

## Part 4 — Research objectives

Each objective states the question, why it matters, and **what result would
falsify the current position**. Answering "no" to any of these is a result.

### R1 — Does the 1–2× advantage survive real metagenomes?

Every real-data result so far is one strain against a complete reference. The
whole claim is usable estimates at 1–2× where Pilea's defaults report nothing.
**Falsified if** in a real community sk2bGrow's advantage at 1–2× disappears, or
its recall advantage is bought by spurious estimates on absent genomes.

### R2 — Where does scaffolding stop working?

`scaffold` was perfect against a different strain of the same species. Real MAGs
are incomplete and contaminated, and the nearest complete relative may be a
different genus. **The failure boundary is the number a MAG user actually needs**
and it cannot be guessed. Falsified if placement degrades already at the genus
level, in which case the marine and RBC datasets become limitations.

### R3 — Can the QC be made to see a destroyed coordinate?

Currently 100% of fragmented estimates pass. A contig-count guard firing at 2 is
the cheap fix, but it is a proxy: a *correctly ordered* multi-contig scaffold is
fine and would be flagged. Is there a statistic that separates "no gradient
because not growing" from "no gradient because the coordinate is scrambled"?
Candidate: within-contig slopes should have consistent |magnitude| under a real
gradient and be uncorrelated under a scrambled one.

### R4 — What compresses the slope at low coverage?

See A4. One candidate left. The decisive experiment is a simulation with **known
λ per window**: compare the recovered window rates to truth at 0.5–1× and see
whether low-count windows are pulled toward the mean before fitting. If they are,
no change to the fit will help.

### R5 — Does the order-free estimator become competitive with an overdispersion term?

`spread_estimator.py` recovers RMSE 0.87 → 0.14 at 5× on contigs and reports
0.000 on the stationary control where Pilea reports 1.153, but over-shrinks below
5× and has a ~0.4 floor at depth. Both failures point at the same missing term,
`s_eff² = s_w² + tau²`. **If that fixes both**, sk2bGrow gains a principled
fragmentation-proof mode and the two-estimator design of Part 3.1 becomes a real
feature rather than a workaround.

### R6 — Does the method scale, and at what panel size?

C6 below. Also the last blocker on shipping 8 enzymes by default.

### R7 — Do the conclusions hold against more tools?

Everything so far is sk2bGrow vs Pilea. Pilea's own paper compares CoPTR, GRiD,
iRep and DEMIC. **A two-tool comparison is not a benchmark**; add them.

---

## Part 5 — Data

Sizes are FASTQ bytes from the ENA file report, queried 2026-08-25.

| BioProject | what | runs | bases | FASTQ | objective |
|---|---|---|---|---|---|
| **PRJNA615952** | *E. coli* K-12 MG1655, 16 defined media + run-out (Zheng et al. 2020) | 45 | 99.4 Gbp | **54.7 GB** | R1 (full depth) |
| **PRJNA1280254** | *B. subtilis*, *K. pneumoniae*, *M. morganii*, *P. putida*, LB at 0.1–2× nutrient | 20 | 48.0 Gbp | **42.5 GB** | R1 cross-species |
| **PRJNA551656** | 20 marine surface-water metagenomes, 4–5 timepoints over 2 days (Long et al.) | 20 | 100.1 Gbp | **59.6 GB** | R1, R2 |
| **PRJNA974210** | rotating biological contactor biofilms + MAGs | 18 | 210 Gbp | **100.5 GB** | R2, application |

**Total for reads: ~260 GB.** References add ~160 GB.

| reference set | where | size | notes |
|---|---|---|---|
| 101 marine MAGs | figshare `10.6084/m9.figshare.9730628` | small | Pilea got estimates for only 18 |
| 525 MAGs (RBC) | under PRJNA974210 | small | remove 3 archaeal |
| 120 NCBI-Pathogen complete genomes | RefSeq | ~500 Mb | simulation set; completeness ≥99.97, contamination ≤1.54, 1–9 contigs |
| 45,529 *Escherichia* assemblies | GTDB R226 | ~60 GB est. | assembly-quality sweep, one at a time |
| GTDB R226 species reps (136,646) | GTDB | ~100 GB est. | R6 only; needs ~750 GB RAM at 16 enzymes, ~670 GB at 8 |

Pilea's global-sludge survey (4,448 SRA sludge samples) is many TB — a separate
project, not part of this benchmark.

**Download procedure.** Write whole files, then subsample. `curl | gunzip | head`
races with SIGPIPE and silently truncates; a partially written FASTQ once made a
sample look too small and it was skipped at 10× without any error. Subsample only
after the download has finished, and verify the byte count against the ENA report
before using a file.

---

## Part 6 — Experiments

Compare against **Pilea v1.3.8 in both configurations** everywhere. Add CoPTR
v1.1.6, GRiD v1.3, iRep v1.1.14, DEMIC for R7.

### C1 — *E. coli* isolate accuracy at full depth (R1)
Repeat the laptop grid without the 600k-read cap: full ~490× per run, subsampled
to 0.25/0.5/1/2/5/10/20×. **Expect** the 1–2× advantage to hold and the 0.5×
point to firm up (n = 16 is thin there). **Failure**: if the advantage was an
artefact of taking the first 600k reads in flowcell order.

### C1b — cross-species (R1)
PRJNA1280254, four species at graded nutrient. **Expect** the *E. coli* result to
generalise. **Failure**: species-specific behaviour, which would mean the panel's
anchor density interacts with the estimate.

### C2 — fragmentation, generalised (R2)
`benches/fragmentation/run.sh` takes any genome. In priority order:
1. **Genuinely incomplete MAGs**, not just cut ones — drop random contigs and
   splice in foreign ones. `fragment.py` preserves every base; a real MAG does
   not.
2. **Progressively distant scaffolding references** — O157:H7, then *Shigella*,
   *Salmonella*, then a different genus. This is R2's boundary.
3. The 16-genome simulation grid under fragmentation.
4. The enzyme-panel sweep under fragmentation (last blocker on A6).

### C3 — multi-strain simulation at full scale (R1)
Pilea's grid: 120 genomes, up to 32 strains, up to 32×, 400 samples. The laptop
ran 16 genomes / 16 strains / 8× / 24 cells. **Expect** recall to stay at 1.00 and
RMSE to degrade gracefully with strain count. **Failure**: recall falling below
Pilea's gates-off arm at high strain counts.

### C4 — marine metagenome (R1, R2)
PRJNA551656 against the 101 MAGs, with Long et al.'s growth rates. Note Pilea
itself got estimates for only 18 of 101 — **report the same denominator**, or the
comparison is recall-vs-accuracy again.

### C5 — rotating biological contactor (R2)
PRJNA974210, 525 MAGs. This is the application section the paper does not have.

### C6 — GTDB-scale scalability (R6)
136,646 species reps against 100 Gbp. Budget **~750 GB at 16 enzymes, ~670 GB at
8, ~290 GB at 2** from the measured 160 B/anchor. Run at **8 enzymes** unless
there is a reason not to. Worth one 2-enzyme run for a like-for-like sketch-size
comparison: a 2-enzyme panel (17,055 anchors on *E. coli*) is almost exactly the
size of Pilea's FracMinHash sketch (18,261 k-mers).

### C7 — the A4 window-rate experiment (R4)
Not a scale experiment; can be done on any machine but has not been. Simulate
with known λ per window, recover rates at 0.5–1×, compare to truth.

---

## Part 7 — Reporting

Per experiment, send back: raw per-sample TSVs with **unchanged column names**,
the `/usr/bin/time` files, and tool versions. Dropping refreshed TSVs into
`sk2bGrow-paper/data/` regenerates every figure and table.

Report per depth, never averaged across depths. Report recall beside every error
metric. Report slope beside every correlation. Say which slope you mean.

**Flag anything that contradicts Part 3**, especially Part 3.2 — the refuted
claims are refuted on laptop-scale data and some of them could come back at
scale. If a result overturns something, that is the result; write it up as one.

---

## Part 8 — Traps already hit

- **`curl | gunzip | head` races with SIGPIPE** and silently truncates.
- **Subsample only after the download finishes.**
- **`pgrep -f script.sh` matches its own command line.** Wait on a log marker.
- **`%.6g` truncates genome coordinates** (3923883 → 3923880). Use `%.10g`.
- **A pandas column named `cov` shadows `DataFrame.cov()`.** Use `depth`.
- **Two enzymes of the same tag length can match the same read window.** Fixed;
  the regression test is
  `a_shared_locus_is_counted_once_per_enzyme_not_once_per_pass`. If you touch
  `count_read`, keep it green.
- **`profile --enzymes` was once parsed and discarded.** Excluded enzymes must be
  left out of the index, not emitted as zeros — zeros tell the statistics layer
  the enzyme was measured and found empty.
- **Pilea writes `output.tsv`, not `profile.tsv`**, and its log2(PTR) is already
  a column.
- **Pilea's cost is non-monotonic in depth.** Never average it.
- **matplotlib lacks a `→` glyph in Helvetica Neue.** Avoid arrows in figure text.
- **`savefig` with `bbox='tight'` widens the canvas to fit one long caption
  line.** Hard-wrap captions.
- **Contributors:** commits in these repositories must show `HuangShiLab` only.
  No `Co-Authored-By` trailers.
