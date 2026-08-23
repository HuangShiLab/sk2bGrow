# Overview

Two source documents define this project:

* **`2bRAD-TGT_PTR算法评估与设计报告.docx`** — the feasibility and algorithm
  argument: why deterministic 2bRAD anchors are a better sketch for PTR than
  FracMinHash, with measured anchor densities on three genomes and Monte-Carlo
  simulations of both routes.
* **`sk2bGrow_代码架构设计.docx`** — the architecture: language split, module
  boundaries, reuse map onto Syn2b/Syn2bANI and Pilea, CLI shape, milestones.

This directory records how those translate into code, and what changed when they
met an implementation.

## The one-sentence version

A 2bRAD tag is a k-mer sketch whose sampling rule is *a restriction motif* rather
than *a hash*. Everything follows from that substitution:

| property | FracMinHash (Pilea) | Type IIB anchors (here) |
|---|---|---|
| sampling rule | `h(x) < H/s`, Bernoulli 1/250 | motif-conditioned, deterministic |
| spacing | geometric, unbounded tail | motif density; 16-enzyme union caps the worst gap near 1.5 kb |
| position | discarded (sort-and-fit) | first-class — the fit's x-axis |
| strata | one | 16, each an independent channel |
| cross-sample loci | same reference ⇒ same | same, *and* physically guaranteed at the bench |
| blind spots | hidden per realisation | auditable at build time |

## Two routes, one kernel

**Route A — pure computation.** Treat the anchor set as a deterministic sketch and
mine existing shotgun metagenomes. Zero marginal cost: it adds a PTR readout to
data that already exists.

**Route B — wet lab.** Real 2bRAD-M libraries already sequence exactly these loci
for taxonomic profiling. Adding a PTR readout costs nothing extra at the bench.

The two share one algorithmic kernel and differ only in a flag (`--mode 2brad`)
and in which noise dominates: Poisson sampling in route A, per-site efficiency in
route B.

## Reading order

1. [`01-architecture.md`](01-architecture.md) — layers, module map, why the
   Rust/Python split falls where it does.
2. [`02-algorithm.md`](02-algorithm.md) — the six pipeline steps and the
   statistical choices in each.
3. [`03-data-formats.md`](03-data-formats.md) — TGT v2, the anchor database, the
   count table. The interface contracts.
4. [`04-roadmap.md`](04-roadmap.md) — milestones and validation gates.
5. [`../enzymes.md`](../enzymes.md) — the panel, and two of its properties that
   the code has to handle explicitly.

## What implementation changed

Three things the design documents did not anticipate, all found by tests:

* **Palindromic enzymes need two-orientation matching.** AlfI, BplI, FalI and
  HaeIV would otherwise lose ~half their reads silently.
  ([`../enzymes.md`](../enzymes.md))
* **HaeIV's recognition site is a strict subset of Hin4I's**, so those two strata
  are not independent. Handled in counting; noted as a caveat for the χ² test.
* **Shared anchors must keep their counts** through the counting layer, or the EM
  reassignment step has nothing to reassign. Discarding multimappers at count
  time silently zeroed the entire shared fraction.

And four statistical adjustments, each found by running the pipeline end to end
rather than by reasoning about it:

* **Window size must adapt per enzyme.** A flat 100 anchors/window drops the three
  sparsest enzymes (PpiI, BplI, PsrI) from a 16-enzyme design even on E. coli.
* **The origin belongs to the genome, not the enzyme.** Searching it per enzyme
  turns search divergence into apparent enzyme disagreement.
* **Fitted-to-noise corrections must shrink.** Both the GC loess and the
  segmented V-shape will fit scatter given the chance, and because each enzyme is
  corrected separately, that noise becomes cross-enzyme discordance.
* **A failed enzyme fit is evidence.** Cochran's Q only sees enzymes that
  produced a number, so a separate gate is needed for those that produced none.

The pattern is worth naming: in a stratified design, *anything fitted per stratum
adds between-stratum variance when it fits noise* — and the cross-enzyme
consistency test, the project's headline QC signal, is exactly what picks that up.
It is a sensitive instrument, and it was sensitive to the pipeline's own
over-fitting first.
