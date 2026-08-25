# sk2bGrow on HPC — handover

Written for whoever runs the large-scale benchmarks. Everything below was
measured on an M3 Max laptop against the reduced datasets in `benches/`; the
numbers are real, the extrapolations are marked as extrapolations.

**Read Part A before starting Part C.** A1 (index memory) has been addressed but
only halves the requirement, so the GTDB-scale experiment is now expensive rather
than impossible; A2 and A3 decide whether the headline claim survives. A3 is
settled — see below — and it changed the answer.

---

## Part A — open algorithm/engineering issues

Priority order. A1 and A3 are done; A2 is open and is the one that decides the
paper's scope. A4–A8 are refinements.

### A1. Index memory — **done, but only 2.2× and still the binding constraint**

`AnchorIndex` stored `HashMap<u64, Vec<u32>>` for the exact table and one more
per seed slot. A separate heap-allocated `Vec` per key, holding typically one
`u32`, was the whole problem. It is now a CSR layout: a sorted `Vec<u64>` of
keys, a `Vec<u32>` of offsets, a `Vec<u32>` of anchor indices, binary-searched.

Measured marginal RAM per anchor — peak RSS from `/usr/bin/time -l` on the
506,785-anchor database minus the 43,735-anchor one, median of five runs:

| `--max-mismatch` | HashMap | CSR | |
|---|---|---|---|
| 0 | 174 | 111 | 1.57× |
| 1 | 356 | 173 | 2.06× |
| 2 (default) | 352 | 160 | **2.20×** |

Count tables are bit-identical before and after. Wall time is ~1% *worse*
(6.39 s against 6.31 s for one 10× *E. coli* sample, median of five): binary
search is not faster than a hash lookup here, and the earlier guess in this
document that CSR would also be quicker was wrong.

**The 12–16 bytes/anchor this document originally projected was also wrong**, by
an order of magnitude, because it costed only the steady-state arrays. Of the
160 bytes actually spent per anchor at the default: ~38 are the `AnchorDb`
itself, held in RAM alongside the index; ~64 are the CSR arrays across the exact
table and three seed slots; the remaining ~58 are build transient — the
`(key, value)` pair buffer and its sort. Seed slots are already built one at a
time to keep that transient down.

Extrapolated to GTDB R226 species representatives (136,646 genomes, ~0.55 Gbp,
~4.7 × 10⁹ anchors at 8,642 anchors/Mb), at mm = 2:

| | 16 enzymes | 8 enzymes | 2 enzymes |
|---|---|---|---|
| HashMap (before) | 1.65 TB | 1.48 TB | 645 GB |
| CSR (now) | **752 GB** | 672 GB | 293 GB |

That is the difference between "needs a distributed index" and "needs one fat
node", which is the change that matters for C6.

Pilea's sketch for the same 16 genomes is 4.46 MB against our 19.3 MB on disk
(0.076 vs 0.33 MB per Mb of genome), so it starts 4.3× smaller.

**Where the remaining factor of ~2 is, in priority order:**

1. *Drop the transient.* Build each CSR table by counting keys, prefix-summing,
   then scattering, instead of materialising and sorting a pair buffer. Removes
   ~58 B/anchor — the single largest remaining item.
2. *Ship 8 enzymes, not 16* (A6). Independent of the layout and worth 10%; the
   panel sweep already says 4–8 is at least as accurate.
3. *Memory-map the `AnchorDb`* rather than deserialising it. Removes ~38
   B/anchor and makes the resident set page-cache-backed.
4. *Shard by genome.* The index is embarrassingly partitionable; a GTDB-scale
   run does not need one process holding everything.

With (1) and (3), 16 enzymes projects to ~300 GB, which is a single fat node
rather than a distributed problem. Until then, keep reference sets under
~4 × 10⁹ anchors per process (roughly 100,000 average bacterial genomes at 16
enzymes on a 700 GB node). **Everything in Part C fits, including C6 at 8
enzymes.**

### A2. Fragmented references and MAGs are untested — **decides the paper's scope**

The V-shape fit regresses on genomic coordinate. On a MAG the contig order and
orientation are unknown, so the coordinate does not exist. `scaffold.rs` (398
lines) exists to order contigs against a reference but **no benchmark exercises
it**, and no result in the paper uses anything but a complete genome.

This is the one place where Pilea is architecturally better: sorted-rank
regression needs no coordinates and is fragmentation-proof by construction. Its
Fig 3 runs every cell twice, on complete genomes and on the same genomes split
into 100 contigs.

Experiment C2 settles it. If sk2bGrow collapses on contigs and `scaffold` does
not rescue it, the honest framing is *"a method for complete references"*, and
the marine/RBC datasets (both MAG-based) become a limitation section rather than
a result.

### A3. The attribution 2×2 — **done; it changed the claim**

The missing cell is now arm E (`benches/zheng2020/armE.sh`). Pearson r at 1×:

| | rank regression | coordinate V-fit |
|---|---|---|
| FracMinHash | Pilea, 0.889 | arm E, **0.940** |
| 2bRAD anchors | arm B, 0.683 | sk2bGrow, **0.981** |

**The factors interact.** The coordinate fit is worth +0.30 r on anchors but only
+0.05 on the sketch, and the sketch effect changes sign with the estimator
(+0.04 under the V-fit, −0.21 under rank regression). The previous reading off
three cells — "the gain is the estimator, the anchors are if anything behind" —
was an artefact of the hole. At 0.5× only the combination works at all: 0.913
against 0.724 for the same estimator on a FracMinHash sketch, with Pilea's own
arm degenerate.

Construction note for anyone repeating it: do **not** feed Pilea's per-window
rates straight into `fit_v_shape`. That returns log₂PTR 2.65 against a measured
1.73, because it removes the GC correction, the ZTP/ZTNB standard errors, the
outlier handling and the fusion QC along with the sorting. `armE_counts.py`
instead rewrites the sketch into our count-table format, so the entire estimator
runs unchanged and the only difference from arm A is which loci are counted.
Sketch positions are not in the `.pdb` and are recovered by replaying `hash64`
over the reference.

### A4. Low-coverage compression is unfixed

Slope of estimated on measured log₂PTR: 0.62 at 0.5×, 0.78 at 1×, 0.84 at 2×,
0.92 at 5×, 0.95 at 10×. Systematic underestimation that shrinks with depth —
the signature of errors-in-variables attenuation: window log-rates are noisy
regressors, and noise in the *predictor* (here, effectively, in the rate) biases
the slope toward zero. Pilea shows the same effect in the opposite direction
(it overestimates, bias +0.168 in simulation, against ours of −0.013).

Candidate fixes, none tried: a measurement-error-aware fit (weight by the ZTP
standard errors already computed, in a Deming/total-least-squares sense rather
than OLS), or a depth-conditioned shrinkage correction calibrated on simulation.

### A5. Enzyme containment is documented but not used

`enzyme::CONTAINMENTS` records that Bsp24I's sites are *totally* contained in
CjePI's, and that ~half are also CjeI's. `fusion.py` still treats all 16 as
independent strata, so Cochran's *Q* has ~15 degrees of freedom, not 16, and is
mildly anti-conservative. Either drop Bsp24I from the strata or model the
covariance.

### A6. The panel should probably be 8 enzymes, not 16

On the *E. coli* panel, averaged over ≤2×: r peaks at k = 8 (0.967 vs 0.959 at
k = 16) and RMSE is best at k = 8 (0.167 vs 0.196). The four sparsest enzymes
(AloI, PpiI, BplI, PsrI, all < 450 anchors) add no accuracy and **double the
negative-control bias** (mean |log₂PTR| on a replication run-out: 0.073 at k = 2,
0.145 at k = 16). Cost is linear in k: 3.2 s / 7.0 s / 8.6 s per sample at
k = 2 / 8 / 16.

A2 may change this — sparse enzymes may matter more on fragmented references,
where every window is short. Re-run the sweep there before fixing a default.

### A7. The origin is always fitted, never annotated

`index --ori` accepts a DoriC/Ori-Finder/*dnaA* annotation but no database in any
benchmark supplies one. Median |fitted − reference| origin error is 90–230 kb at
0.5× and 9–15 kb at 10×, and origin error propagates into every enzyme's slope.
Supplying an annotation should mostly help at low depth. Untested.

### A8. One species = one strain = one PTR

No model for strain mixtures within a species. Real metagenomes have them. Out
of scope for this paper; state it as a limitation.

---

## Part B — datasets

Sizes are FASTQ bytes from the ENA file report, queried 2026-08-25.

| BioProject | what | runs | bases | FASTQ | used for |
|---|---|---|---|---|---|
| **PRJNA615952** | *E. coli* K-12 MG1655, 16 defined media + run-out controls (Zheng et al. 2020) | 45 | 99.4 Gbp | **54.7 GB** | C1 isolate accuracy |
| **PRJNA1280254** | *B. subtilis*, *K. pneumoniae*, *M. morganii*, *P. putida* in LB at 0.1–2× nutrient | 20 | 48.0 Gbp | **42.5 GB** | C1b cross-species |
| **PRJNA551656** | 20 marine surface-water metagenomes, 4–5 timepoints over 2 days (Long et al.) | 20 | 100.1 Gbp | **59.6 GB** | C4 real metagenome |
| **PRJNA974210** | rotating biological contactor biofilms + MAGs | 18 | 210 Gbp | **100.5 GB** | C5 application |

Reference sets:

| what | where | size | notes |
|---|---|---|---|
| 101 marine MAGs | figshare `10.6084/m9.figshare.9730628` | small | Pilea got estimates for only 18 of them |
| 525 MAGs (RBC) | under PRJNA974210 | small | 3 archaeal to be removed |
| 120 NCBI-Pathogen complete genomes | RefSeq | ~500 Mb | simulation reference set; completeness ≥ 99.97, contamination ≤ 1.54, 1–9 contigs |
| 45,529 *Escherichia* assemblies | GTDB R226 | **~60 GB est.** | assembly-quality sweep; used one at a time |
| GTDB R226 species reps (136,646) | GTDB | **~100 GB est.** | scalability only; needs a ~700 GB node, or ~300 GB after A1's follow-ups |

Marine growth rates come from Long et al.; *E. coli* growth rates are in
`benches/zheng2020/growth_rates.tsv` (already extracted).

**Total for C1–C5: ~260 GB of reads.** C6 adds ~160 GB of references.
Pilea's global-sludge survey (4,448 SRA sludge samples) is many TB — treat it as
a separate project, not part of this benchmark.

---

## Part C — experiments

Each block gives the objective, what to run, the expected result, and what would
count as a failure. Compare against **Pilea v1.3.8** in both configurations
(shipped defaults, and `-x 0 -z 0 -c 0`) — reporting only one is misleading in
one direction or the other. Where Pilea's paper compares more tools, add CoPTR
v1.1.6, GRiD v1.3, iRep v1.1.14, DEMIC.

Measure cost with `/usr/bin/time -l` (macOS) or `/usr/bin/time -v` (Linux); on
macOS it does propagate a grandchild's `maxrss`, so the Python statistics layer
is included. Pilea's paper used `gtime -v` on a 32-core M3 Ultra with `purge`
before each timing run.

### C1 — *E. coli* isolate accuracy, full depth (PRJNA615952)

**Objective.** Reproduce Pilea's Fig 2a at full depth, which we could not do on
a laptop, and confirm the coverage titration extends to it.

**Run.** All 45 runs at native depth (≈490× — subsample to 100× if that is
cheaper; it is well above every gate). Then the titration 0.5/1/2/5/10/20/50×.
Reference GCF_000005845.2. Keep the three run-out samples as negative controls;
Pilea excluded them.

**Expect.** sk2bGrow r ≈ 0.97–0.98 against measured growth rate, matching
Pilea's published 0.976, and **flat from 1× upward** — our laptop titration
already gives 0.981 at 1× and 0.982 at 2×. RMSE against λ·C/ln2 ≈ 0.03–0.06
above 5×. Run-out |log₂PTR| < 0.1.

**Failure.** r below 0.95 at full depth, or a run-out estimate above 0.3, means
something in the full-depth path differs from the subsampled path — check for
duplicate reads first (`fastp --dedup`; the RBC samples needed it).

### C2 — reference fragmentation — **the important one**

**Objective.** Settle A2. Pilea's Fig 3 runs every simulation cell on complete
genomes and on the same genomes split into 100 contigs (lognormal lengths,
μ = 0, σ = 1).

**Run.** Same design, three reference conditions: complete, 100 contigs, 100
contigs + `sk2bgrow scaffold` against a complete relative. Also re-run the
enzyme-panel sweep (A6) under fragmentation.

**Expect.** Honestly: unknown, and we should not pretend otherwise. The V-fit
needs coordinates the contigs do not carry. Best case, `scaffold` recovers
enough order that sk2bGrow degrades gracefully; worst case it collapses to
chance and Pilea's rank regression wins outright on fragmented references.

**Either outcome is publishable** — the second one just changes the title.

### C3 — multi-strain simulation at Pilea's full scale

**Objective.** Pilea's Fig 3a–f at their scale, not our laptop scale.

**Run.** 120 complete NCBI-Pathogen genomes; 4/8/16/32 strains × 4/8/16/32×
coverage × 5 replicates × 5 samples per replicate = **400 samples**. Origin at
position 0, terminus at midpoint, coverage decaying log₂-linearly,
log₂PTR ~ U[0,2]. `benches/simulate.py` already implements exactly this
generative model — raise the grid, do not rewrite it.
Add 1× and 2× coverage rows: that is the regime this paper is about and Pilea's
grid starts at 4×.

**Expect.** Pilea reports L2 *d* = 10.681 over its 400 samples (GRiD 39.070).
From our 24-cell laptop grid, sk2bGrow gives recall 1.000 / RMSE 0.134 /
bias −0.013 against Pilea-gates-off 0.997 / 0.265 / +0.168 and Pilea-defaults
**recall 0.224**. Expect that ordering to hold, with sk2bGrow's advantage
largest at 1–4× and shrinking above 16×.

**Report recall beside L2.** L2 is computed only over the genomes a tool chose
to report, so it rewards silence — Pilea's flattering 0.083 RMSE at defaults
comes from answering 22% of cases, the easy ones.

### C4 — marine metagenome (PRJNA551656)

**Objective.** Pilea's Fig 3g and the scalability claim, at a size a laptop
cannot hold. 100 Gbp of reads against 101 MAGs.

**Run.** All 20 samples against the 101 MAGs. Compare per-MAG correlation
between estimated log₂PTR and Long et al.'s cell-count-derived growth rates.
Keep only MAGs with > 3 paired observations, discard negative growth rates.
Record wall-clock per 10 Gbp and peak RSS (Pilea's Fig 3h,i).

**Expect.** Pilea returned estimates for only **18 of 101** MAGs and had the
highest median correlation. The interesting number for us is **how many MAGs we
return** — if the anchor panel plus the fusion QC returns 40 with comparable
correlation, that is a strong result. If we return 18 with the same correlation,
we have matched, not beaten.

**Caveat.** These are MAGs, so C2 gates the interpretation. Run C2 first.

### C5 — application: rotating biological contactor (PRJNA974210)

**Objective.** The equivalent of Pilea's Fig 4 — a biological claim, not a
benchmark. Abundance-weighted log₂PTR of ammonia- and nitrite-oxidising
bacteria along the flowpath, against NH₄⁺-N and NO₂⁻-N.

**Run.** `fastp --dedup` first (these samples have > 10% duplicates). 522 MAGs
after removing 3 archaeal, plus a GTDB R226 AOB/NOB/COM reference set
(*Nitrospira* n = 378, *Nitrosomonas* n = 266, *Nitrotoga* n = 89,
*Nitrobacter* n = 22, *Nitrosymbiomonas* n = 21, *Nitronauta* n = 3).

**Expect.** The point is a coherent AOB→NOB gradient along the flowpath that
tracks the nitrogen chemistry. Do this last; it is only meaningful once C2 has
established whether MAG references are usable at all.

### C6 — GTDB-scale scalability — **needs a large-memory node**

**Objective.** Pilea's headline: 136,646 GTDB species representatives profiled
against 100 Gbp in under ten minutes on 32 threads.

Budget from the measured 160 B/anchor: **~750 GB at 16 enzymes, ~670 GB at 8,
~290 GB at 2**. Before the CSR rewrite this was 1.65 TB and there was no node to
run it on; it is now a single fat node, and A1's follow-ups (drop the build
transient, mmap the database) would bring 16 enzymes to ~300 GB.

Run it at **8 enzymes** unless there is a reason not to: the panel sweep (Table
7) shows 8 is at least as accurate as 16, and it saves 10% of the index and 40%
of the wall time.

Note that a 2-enzyme panel (17,055 anchors on *E. coli*) is almost exactly the
size of Pilea's FracMinHash sketch (18,261 k-mers) — that is the honest
like-for-like comparison of the two sketches, and worth one run even though 2
enzymes is not the recommended default.

---

## Part D — traps already hit, do not repeat

- **`curl | gunzip | head` races with SIGPIPE** and silently truncates the
  download. Write the whole file, then subsample.
- **Subsample only after the download has finished.** A partially written FASTQ
  made one sample look too small and it was silently skipped at 10×.
- **`pgrep -f script.sh` matches its own command line.** Wait on a marker in the
  log, not on a process listing.
- **`%.6g` truncates genome coordinates** (3923883 → 3923880). Use `%.10g`.
- **A pandas column named `cov` shadows `DataFrame.cov()`.** Bracket access, or
  name it `depth`.
- **Pilea's cost is non-monotonic in depth**: 5.6 s at 0.5×, **21.1 s at 2×**,
  5.6 s at 10× on our cells, because its ZTP-mixture EM (`--max-iter` defaults to
  infinity) is slowest where the mixture is least identifiable. Do not average
  cost across depths and quote one number.
- **Two enzymes of the same tag length can match the same read window.** Fixed in
  82eda91; the regression test is
  `a_shared_locus_is_counted_once_per_enzyme_not_once_per_pass`. If you touch
  `count_read`, keep it green.
- **Contributors:** commits in these repositories must show `HuangShiLab` only.

---

## Part E — what to send back

Per experiment: the raw per-sample TSVs, the `/usr/bin/time` files, and the tool
versions. `sk2bGrow-paper/figures/make_figures.py` and `make_tables.py` read only
`data/*.tsv` and regenerate every figure and table with no network access, so
dropping refreshed TSVs into `data/` is the whole integration step. Keep the
column names as they are.

Flag anything that contradicts Part A — especially A2 and A3, where the current
answer is "we do not know" and the paper's framing depends on it.
