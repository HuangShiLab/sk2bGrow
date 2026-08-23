# The algorithm

Six steps, following report §7.1. Each names the defect it addresses.

## 1. Index (offline)

Digest each reference with the 16-enzyme panel, record for every anchor its
coordinate, enzyme, strand, ±250 bp GC and uniqueness flags, and store the result
as TGT v2 plus a binary anchor table.

Uniqueness is resolved once, into two independent bits:

* `UNIQUE_IN_GENOME` — the tag occurs once in its own genome. Multi-copy anchors
  contaminate coverage estimates (Pilea masks multi-copy k-mers for the same
  reason).
* `UNIQUE_ACROSS_DB` — the tag occurs in one genome of the database. Shared
  anchors are excluded from coverage modelling but *retained* for step 2.

Budget: ~28 000 anchors × ~44 B ≈ 1.2 MB per E. coli-scale genome.

## 2. Count (online)

Scan reads for the 16 recognition motifs; at each hit extract the implied tag and
match it against the anchor table with a ≤2 mismatch budget (pigeonhole seeding,
then Hamming verification). Both tag orientations are probed — see
[`../enzymes.md`](../enzymes.md) for why that is not redundant.

Tags matching anchors in several genomes are recorded on **all** of them, and
`em.rs` then splits the mass in proportion to current abundance, iterating until
the split stops moving. Abundance is re-estimated from unique anchors only on
every iteration; letting shared anchors feed back into the abundance used to
split them makes the fixed point unidentifiable, and a genome with no unique
anchors could inflate without bound.

> **Addresses D5.** A mismatch budget extends usable reference distance below the
> 95 % ANI wall that exact k-mer matching imposes.

## 3. Window

Cut windows inside each `(genome, enzyme)` anchor series with a fixed **anchor
count** (default 100), not a fixed base-pair span.

Pilea's 25 kb window holds ~100 k-mers *on average*, but the Bernoulli sampling
makes occupancy vary and collapse on fragmented references. Equal-anchor windows
equalise statistical power along the genome; window *width* then varies instead,
and that variation is itself a legible signal — a wide window marks a sparse
region. `--windowing bp` restores fixed-bp windows for A/B benchmarking.

**The anchor count is per enzyme, not global.** Panel densities span 20×
(CjeI 1 910/Mb, PpiI 73/Mb). At a flat 100 anchors/window on E. coli, CjeI gets
88 windows while PpiI, BplI and PsrI get 3–4 — below any usable minimum, so the
three sparsest enzymes drop out of a "16-enzyme" design entirely. `auto_window_size`
targets ~25 windows per enzyme within a 25–100 anchor band, which keeps all
sixteen channels alive; sparse enzymes get smaller, noisier windows and therefore
*earn less weight* in the fusion, which is the right way for them to count less.

> **Addresses D2.**

## 4. Estimate a rate per window

Fit the window's positive counts with a zero-truncated Poisson mixture (EM,
component count by BIC) or a zero-truncated negative binomial, whichever BIC
prefers.

Truncation matters: systematic divergence makes some reference anchors genuinely
absent, and counting those as true zeros drags a plain Poisson fit down. The NB
branch accommodates residual overdispersion instead of discarding the sample the
way Pilea's dispersion filter does; at low coverage NB is unidentifiable, so BIC
falls back to ZTP without needing a hand-set coverage threshold.

Every rate carries a standard error, from the delta method:

```
se(rate) = sqrt(Var_model[X | X ≥ 1] / n_eff) / |d mean / d rate|
```

Same construction for every model, so BIC can switch between them without making
the errors incomparable. Those errors are what step 6 weights by.

> **Addresses D1, D6.**

*Implementation note.* Computing the NB log-pmf as
`gammaln(k+r) − gammaln(r)` cancels catastrophically as α → 0 (r = 1/α → ∞:
two terms of order 10¹³ differing by ~10). The optimiser reads the resulting
rounding noise as free likelihood and NB "wins" on perfectly Poisson data. The
implementation uses the rising-factorial form `Σ log(r+i)`, which is exact for
integer counts.

## 5. Correct GC bias

Fit a loess curve of log2 count against local GC **inside each enzyme**, and
apply it as a log2 offset **per anchor**, averaged into the window.

Each Type IIB site has its own GC composition, so the 16 enzymes sample 16
different narrow GC neighbourhoods; one global curve averages them into a shape
that fits none. And since GC is stored per anchor, correcting on a window's mean
GC would discard exactly the within-window variation being corrected for.

The correction is an offset in log2 space, never a rescaling of counts — the
window models need integer counts, and a multiplicative fudge would break the
zero-truncation.

What this does *not* need to fix: a constant per-enzyme efficiency factor cannot
bias PTR at all, because each enzyme is fitted separately and a constant is
absorbed by that fit's intercept. Only GC *slope* within an enzyme matters. The
factors are still reported, for QC.

**Curves are shrunk toward flat.** A loess will trace Poisson scatter happily,
and a sparse enzyme has few anchors to average over — on a genome with no GC bias
at all, testing produced spurious "corrections" of ~0.5 log2. Because each enzyme
gets its own curve, and each curve then distorts its enzyme differently, that
noise becomes *between-enzyme* disagreement and the consistency test fires on an
artefact of the correction. Each curve is therefore scaled by the variance it
explains beyond the null expectation for fitting noise (`edf/n`): a real gradient
keeps ~all its amplitude, a curve tracing scatter shrinks to exactly zero.

> **Addresses D6.**

## 6. Fit the profile

**Origin known** (DoriC / Ori-Finder / `dnaA`): regress log2 window rate directly
on circular distance from the origin.

```
log2 μ(x) = a − b₁·min(d, k) − b₂·max(0, d−k),    d = dist(x, ori)
```

With `b₁ = b₂` this is the plain V that CoPTR-Ref shows is the maximum-likelihood
model. The two-slope form exists for multi-fork replication: above PTR ≈ 2,
overlapping rounds put a genuine kink in the profile, and a single line through a
kinked profile is biased.

The segmented form is **gated, not merely BIC-selected**. With ~20 windows, two
extra parameters buy a BIC improvement by chance often enough that a plain V gets
reported as kinked — and since the two forms yield different log2(PTR), that
again surfaces as enzymes disagreeing. A kink is a multi-fork phenomenon, so it
is offered only where multi-fork is physically possible (log2 PTR ≳ 1) and there
are ≥ 30 windows to resolve it. Within that gate, BIC decides.

**One origin per genome, not per enzyme.** The origin is a property of the
chromosome; the enzymes are independent measurements *of the same gradient*, not
independent genomes. Searching separately per enzyme wastes power — each search
sees a fraction of the windows — and injects between-enzyme variance that has
nothing to do with biology: two enzymes landing on origins 200 kb apart report
different slopes, and the consistency test reads that as disagreement when in
fact they agree and the *search* diverged. `find_shared_ori` pools all enzymes
(each median-centred first, to remove its efficiency offset) and each enzyme's
slope is then fitted at that fixed coordinate. `--per-enzyme-ori` restores
independent searches as a diagnostic.

**Origin unknown**: grid-search it jointly with the slope, then refine. Candidate
origins whose fit slopes *uphill* are rejected outright — that solution is the
same line read from the terminus, and letting it win reports the ter as the ori.
The returned `ori_confidence` is the circular mean resultant length of the
posterior over origin position: 1 means the data pin it down, near 0 means they
are consistent with it being almost anywhere. A slow grower has little gradient
and therefore little information about the origin, and the output says so instead
of reporting a confident coordinate.

**Fragmented reference**: fall back to sorted regression + RANSAC. A MAG with no
contig order has no x-axis; `sk2bgrow scaffold` supplies one, and until it has,
the coordinate fit would be fitting noise.

Standard errors are scaled by the fit's reduced χ². The window errors from step 4
describe counting noise only; anchor efficiency, residual GC structure and
profile misspecification all add scatter. Taking the residuals at face value
keeps the error bars honest.

> **Addresses D3, D8.**

## 7. Fuse across enzymes, and test them against each other

Each enzyme yields `log2(PTR) ± se`. Combine by inverse-variance weighting — the
minimum-variance linear combination — so enzymes with more anchors and cleaner
fits count for more, automatically.

Then test them against each other. Under the null that all 16 measure one common
value, Cochran's Q is χ² on `k−1` degrees of freedom. A significant Q means the
enzymes *disagree*: an enzyme with too few anchors, a methylation-blocked site
class, a mis-assembled region. **This is a real replicate structure, available at
zero extra sequencing cost, and no single sketch can produce it.**

When Q rejects, the fixed-effect standard error is known to be too small — it
assumes the only scatter is sampling noise. The estimator escalates to
DerSimonian-Laird random-effects weights, which add the between-enzyme variance
component. Reporting a tight interval around a value the enzymes visibly disagree
about would be the worst of both worlds.

Caveat, from [`../enzymes.md`](../enzymes.md): HaeIV's recognition site is a
strict subset of Hin4I's, so those two strata are not fully independent and their
agreement is partly structural. With 16 enzymes this is minor; with a two-enzyme
subset of exactly those two it would be meaningless.

> **Addresses D4.**

## 8. Report

`output.tsv` keeps Pilea's column names so existing benchmark scripts work
unchanged, and appends `enzyme_consistency`, `n_anchors`, `ori_confidence` plus
the fusion diagnostics.

QC gates flag rows; they never delete them. A tool that silently drops rows makes
"this genome was not growing" and "this genome was filtered out"
indistinguishable downstream.

One gate exists specifically because the χ² test cannot see it. Enzyme
discordance has **two** modes: an enzyme that fits to a *different slope*
(caught by Q), and an enzyme that produces *no usable fit at all* — a flat
profile, a digestion failure, a mis-assembled region. The second never reaches
the Q statistic, so without a separate gate the surviving enzymes agree with each
other and a sample where a quarter of the panel saw nothing reads as clean.
`enzyme_fit_rate` closes that hole. The coverage floor defaults to 1× rather than
Pilea's 5× — that is the point of the design — but it is a default, and
`--min-coverage 5` restores Pilea-comparable strictness.

## 9. Dynamics

Because the anchor set is fixed by the enzymes, the same loci are re-observed in
every sample. A time series stops being "estimate a noisy curve per sample, then
compare curves" and becomes a repeated-measures table on fixed loci — the object
mixed-effects models are built for. `dynamics.anchor_matrix` exposes the raw
anchor × sample counts; `delta_ptr` and `trend_test` cover the summary view.

## Known limits

Inherited from PTR itself, not from this implementation:

* PTR measures replication activity, not growth rate. Converting needs a
  species-specific C-period, so PTR is not comparable across distant taxa.
* Relic DNA from dead cells flattens the profile.
* Plasmids carry recognition sites but no ori-ter gradient — hence the
  `NON_CHROMOSOMAL` flag and `ContigKind`.
* Multi-fork replication is only partly captured by the segmented model.
* Long-read libraries violate the V-shape assumption through size selection.

And one specific to this design: closely related strains share most anchors, so
strain-level PTR needs strain-specific anchors and accepts the sensitivity loss.
