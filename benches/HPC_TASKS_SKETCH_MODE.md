# HPC task list — adding FracMinHash as a second landmark source

> **Ranks below [`HPC_TASKS_PAIRED_2BRAD.md`](HPC_TASKS_PAIRED_2BRAD.md).** That
> document measures σ_eff — roadmap open question R1, the largest open risk — from
> paired WGS+2bRAD data already available, and its arm A *is* the FracMinHash arm
> below. Run it first; these tasks stay queued behind it.

**Additive to `RESEARCH_PLAN.md`. Do not disturb work already running.** R1–R7 and
C1–C7 stand as written; nothing here changes them. These tasks share the same
datasets and the same reporting contract, so they slot into gaps rather than
replacing anything. Where a task can reuse a run that is already scheduled, it
says so.

Prompted by `Syn2b/docs/LANDMARK_COMPARISON.md`, which measured 2bRAD digestion
against FracMinHash for the *synteny* metric and found the structural mathematics
to be landmark-agnostic. The question here is the narrower one for us: **should
`sk2bgrow index` gain a `--mode fracminhash` alongside the enzyme panel, for WGS
input?**

---

## 1. What we already know, before spending any HPC time

### 1.1 The estimator is already sketch-agnostic — we ran it

Arm E (`benches/zheng2020/armE_counts.py`) rewrites Pilea's FracMinHash sketch
into our count-table format and runs the **unmodified** estimator on it. Pearson
r / RMSE against measured growth rate, Zheng *E. coli*, n = 16:

| depth | 2bRAD anchors | FracMinHash (scale 250) |
|---|---|---|
| 0.5× | **0.913** / 0.304 | 0.724 / 0.382 |
| 1× | **0.981** / 0.157 | 0.940 / 0.213 |
| 2× | 0.982 / 0.128 | **0.984** / 0.093 |
| 5× | **0.979** / 0.039 | 0.977 / 0.051 |
| 10× | **0.968** / 0.063 | 0.942 / 0.149 |

So the question is not *whether* a sketch works in our pipeline. It does. GC
correction, ZTP/ZTNB window rates, adaptive windows, origin search, V-fit and
fusion all run on it unchanged. This is an engineering question, not a research
one — **except for what follows.**

### 1.2 That comparison was not fair, and the unfairness is quantified

Arm E used Pilea's default `--scale 250`, which on *E. coli* yields 18,261
landmarks = **3,934 /Mb**. Our 16-enzyme panel yields 43,735 = **9,422 /Mb**.

**The sketch arm was run at 2.4× lower landmark density than the enzyme arm.** At
0.5× read coverage that is exactly the difference between populated and empty
windows, and low coverage is where the gap is largest (0.913 vs 0.724). Density
is a continuous knob on the sketch side — Syn2b confirms the count tracks
`genome_length / scale` to within 2% from scale 250 down to 50 — so the matched
comparison is one parameter away and has never been run.

Scales that match each panel size on *E. coli*:

| panel | landmarks/Mb | matching `--scale` |
|---|---:|---:|
| k = 2 | 3,674 | 268 |
| k = 4 | 5,333 | 184 |
| k = 8 | 8,021 | 123 |
| k = 16 | 9,422 | 104 |

**This is the single most important untested thing in this document**, and it
needs no new Rust: `armE_counts.py` already takes `-s`.

### 1.3 GC dependence — a real scope limit, and it interacts with A6

Syn2b measured landmark density across GC 0.25–0.75 on synthetic sequence: BcgI
alone varies **9.7×**, a four-enzyme panel **2.9×**, FracMinHash **1.04×**.

Our panel is better than that, and we can show why from `data/density_3genomes.tsv`
(*B. subtilis* 43.5%, *E. coli* 50.8%, *P. putida* 61.5% GC). Individual enzymes
swing 2–3.2× over just those 18 GC points (AlfI 237→762 tags/Mb), but the
**panel total varies only 1.16×** — enzymes with opposite GC preferences cancel.

That cancellation is bought with panel size, and it decays as the panel shrinks:

| panel | density ratio over GC 43.5–61.5% | GC slope, % of mean per GC point |
|---|---:|---:|
| k = 2 | 1.51× | 1.99% |
| k = 4 | 1.31× | 1.31% |
| k = 8 | 1.22× | 0.98% |
| k = 16 | **1.16×** | **0.78%** |

**This is a dimension A6 never considered.** A6 recommends shipping 8 enzymes
because on *E. coli* they are as accurate as 16 and 25% cheaper. *E. coli* is one
GC value. If GC robustness is part of what the other eight enzymes buy, the
recommendation needs re-deriving across GC, not re-confirming at 50.8%.

Two warnings about the table above. It is three genomes over a **narrow** range,
and the GC values are literature values, not measured from the FASTA.
Extrapolating linearly to GTDB's real range (25–72%) gives 1.46× at k = 16 and
1.63× at k = 8, but Syn2b's synthetic sweep is **non-monotonic** — their
four-enzyme panel peaks at GC 0.65 and falls at 0.75 — so a linear extrapolation
is probably optimistic. **Measure it; do not quote the extrapolation.**

### 1.4 Mismatch tolerance is our exposure, and it is larger than Syn2b's

Syn2b §3.3: enzyme landmarks must contain a recognition motif, so they are
crammed into a small region of sequence space and a unique landmark can sit **one
substitution** from a multi-copy family — 0.34% of four-enzyme landmarks are at
risk, against **0.00%** for FracMinHash, whose k-mers are drawn from the whole
4^31 space with no shared constraint.

We match reads to anchors with `--max-mismatch 2` by default. **Our exposure is
at Hamming distance 2, not 1**, over a landmark set that is more motif-constrained
than theirs. This has never been measured in sk2bGrow and it is a
misassignment channel, not a sensitivity one: a read from locus A credited to
locus B puts coverage at the wrong genomic coordinate, which is precisely what
the V-fit integrates.

### 1.5 What we can already tell Syn2b

Their §4 lists as unmeasured: *"Behaviour on reads rather than assemblies.
FracMinHash needs an assembly. Whether a read-level sketch reproduces these
results is untested."*

Arm E **is** that experiment for our metric: FracMinHash landmarks defined on the
assembled reference, counted from reads, through a coordinate-dependent
estimator, on real data with an independent ground truth. It works (r 0.940 at
1×, 0.984 at 2×). That does not answer their adjacency question, but it removes
the read-level sketch from their "might not work at all" column.

---

## 2. The design under test

`sk2bgrow index --mode {2brad|fracminhash}`, with `--scale` and `--kmer` for the
sketch mode. Everything downstream consumes `(landmark identity, coordinate,
contig, strand, local GC, flags)` and does not care where landmarks came from —
arm E proves this empirically, not just by inspection.

The honest framing, following Syn2b: **these are not competitors.** 2bRAD is a
wet-lab protocol that yields landmarks from a mixed sample without assembly;
FracMinHash is an in-silico rule that needs an assembled reference. For **route B
(real 2bRAD reads) the enzyme panel is not optional**, because the enzymes are
what was in the tube. The question is only about **route A (WGS reads against a
reference)**, where today we impose a wet-lab constraint on a purely
computational problem for no stated reason.

---

## 3. Tasks

Priority order. F1 and F2 are worth doing before anything else here; F6 and F7
are optional.

### F1 — Density-matched sketch comparison — **do this first** ★★★

**Objective.** Settle whether arm E's low-coverage deficit is a property of
FracMinHash or an artefact of running it at 2.4× lower density.

**Procedure.** No new code. On the Zheng grid, for `--scale` in
{268, 184, 123, 104, 60, 30}:

```bash
pilea_env/bin/python benches/zheng2020/armE_counts.py \
    -r ecoli.fna -s <scale> -o counts_s<scale> sub/*.fq
# then profile each count table exactly as armE.sh does
```

Score with `benches/zheng2020/analyze.py`. Report against the matched enzyme
panel, not against k = 16 only: scale 268 vs k = 2, 123 vs k = 8, 104 vs k = 16.

**Cost.** ~6 × 85 profile runs plus sketch construction. Hours, not days. The
reference sketch is scanned once and reused across samples.

**Expect.** The 0.5× gap to close substantially. If it closes completely, the
claim "deterministic anchors are what keep windows populated at low depth"
(Fig 8 caption, `outline.md` §4, and the A3 write-up) is **wrong** and must be
rewritten — the anchors' remaining advantages would be wet-lab realisability and
scaffolding, not low-coverage accuracy.

**Falsifies the current position if:** at matched density FracMinHash equals or
beats the panel at 0.5–1×. That is a publishable correction, not a failure.

**Also report:** landmarks actually observed per window at each scale, so the
mechanism is visible rather than inferred.

### F2 — GC sweep, both modes, both panel sizes ★★★

**Objective.** Establish the GC range over which the enzyme panel is a usable
instrument, and whether A6's 8-enzyme recommendation survives outside *E. coli*.

**Procedure.** Pick ~20 complete genomes spanning GC 25–75% — GTDB R226 has
plenty; *Buchnera* (~25%), *Campylobacter* (~30%), *Staphylococcus* (~33%),
*Bacillus* (~43%), *E. coli* (~51%), *Klebsiella* (~57%), *Pseudomonas* (~62%),
*Mycobacterium* (~66%), *Streptomyces* (~72%). For each, measure:

1. landmark density, panel k ∈ {2, 4, 8, 16} and FracMinHash at matched scale;
2. **PTR accuracy**, by simulating reads with a known planted gradient
   (`python/sk2bgrow/simulate.py`) at 0.5/1/2/5×, both modes.

Density alone is not the deliverable — a panel can be dense and still fit badly.

**Expect.** Panel density to fall at both GC extremes and the fall to be steeper
at k = 8 than k = 16. FracMinHash flat.

**Decides:** whether the paper claims a GC range, and whether A6 ships 8 or 16.
If the k = 8 panel degrades below GC 35% while k = 16 holds, **A6's
recommendation is wrong outside *E. coli*** and the default stays at 16.

**Note the non-monotonicity** Syn2b found. Sample GC densely enough to see a
peak, not just the endpoints.

### F3 — Misassignment under mismatch tolerance ★★

**Objective.** Measure the §1.4 risk, which is specific to us because we allow
two mismatches.

**Procedure.** For each reference and each mode, count landmark pairs within
Hamming distance 1 and 2 of each other, splitting unique-vs-unique from
unique-vs-multi-copy-family. Then, empirically: simulate reads from known loci
and measure the fraction credited to the wrong coordinate, at
`--max-mismatch` 0, 1, 2, for both modes.

**Expect.** A higher near-collision rate for enzyme landmarks than for
FracMinHash at matched density, and growth from d = 1 to d = 2.

**Decides:** whether `--max-mismatch 2` should remain the default, and whether
the sketch mode can safely use a *higher* tolerance than the enzyme mode. Note
the memory interaction: seed tables dominate the index (A1 — 174 → 356 B/anchor
going from mm = 0 to mm = 1), so if mm = 1 is safe the GTDB budget nearly halves.

### F4 — Fragmentation × sketch mode ★★

**Objective.** A2 showed fragmentation destroys the coordinate and `scaffold`
restores it. Both are landmark-agnostic in principle; check that in practice.

**Procedure.** `benches/fragmentation/run.sh` with FracMinHash landmarks. Two
distinct questions:

1. Does the V-fit collapse identically on contigs? (It must — the mechanism is
   the lost coordinate, not the landmark type.)
2. **Can `scaffold` place contigs using FracMinHash landmarks?** Scaffolding
   needs landmarks *shared between draft and reference*; FracMinHash landmarks
   are shared exactly as enzyme tags are. If it works, scaffolding stops being a
   reason to keep the enzyme panel in route A.

**Expect.** (1) identical collapse. (2) works, possibly better — Syn2b measured
higher landmark retention under divergence for FracMinHash (94.6% vs 89.5% at
0.1% substitution), and cross-strain scaffolding is exactly a divergence problem.

### F5 — Index size and cost at matched density ★★

**Objective.** Feed the C6 budget. A1 measured 160 B/anchor and projected 752 GB
for GTDB species reps at 16 enzymes.

**Procedure.** Build both modes at matched density on a few hundred genomes;
measure peak RSS, index build time, and profile wall time per sample.

**Expect.** Near-identical — the index stores `(u64 hash, u32 anchor)` either
way. If FracMinHash is materially cheaper, it is because a hash threshold needs
no motif scan at index time, which matters at GTDB scale.

**Also:** at matched density the sketch has a knob the panel does not. If GTDB at
16 enzymes needs 752 GB, `--scale 200` gives half the landmarks and half the
index, with accuracy measured rather than assumed. That may be the practical
answer to C6.

### F6 — Read-level sketch, reported back to Syn2b ★

Package the F1 result as an answer to their §4 open question. Costs nothing extra
once F1 is done; it is a write-up.

### F7 — Hybrid landmark sets ★

Union of enzyme tags and sketch k-mers as separate strata in the existing fusion.
Only worth testing if F1 shows the two modes have **different** failure regimes —
if the sketch wins everywhere at matched density there is nothing to combine.
Do not start this before F1 reports.

---

## 4. What not to do

- **Do not remove the enzyme path.** Route B (real 2bRAD reads) requires it: the
  enzymes are what was in the tube. Any sketch mode is additive, for route A only.
- **Do not compare at unmatched density again.** That is the mistake this
  document exists to correct. Every comparison states landmarks/Mb for both sides.
- **Do not quote the GC extrapolation in §1.3** as a result. It is three genomes
  over a narrow range, linearly extrapolated, against a response Syn2b showed to
  be non-monotonic. F2 replaces it with measurement.
- **Do not re-derive A6 from *E. coli* alone.** That is what made it incomplete.
- **Do not let this delay C1–C7.** The running plan answers whether the method
  works at scale; this answers what the landmarks should be. If HPC time is
  scarce, F1 alone (hours) is worth more than F2–F7 combined.

---

## 5. What would change in the paper

Ordered by how much each would cost us to discover late.

| if F1 shows | consequence |
|---|---|
| matched-density FracMinHash equals the panel at 0.5–1× | **The low-coverage claim for deterministic anchors is wrong.** Fig 8's caption, `outline.md` §4 and the A3 interaction result all rest on arm E's deficit. The contribution becomes the coordinate-aware estimator plus wet-lab realisability, and the paper is *more* defensible for being narrower |
| matched-density FracMinHash still loses at 0.5× | the anchors' advantage is real and is about *where* landmarks sit, not how many. Say so, with the density control that proves it |
| **if F2 shows** | |
| k = 8 degrades outside GC 40–60%, k = 16 holds | A6 is wrong outside *E. coli*; ship 16, and state the GC range |
| both hold across 25–75% | the panel is a uniform instrument; a claim we currently cannot make |
| **if F3 shows** | |
| material misassignment at mm = 2 | lower the default, and the GTDB index budget roughly halves as a side effect |
