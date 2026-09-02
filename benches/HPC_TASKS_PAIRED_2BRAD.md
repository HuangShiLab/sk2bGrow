# HPC task list — paired WGS + 2bRAD, and the first measurement of σ_eff

**Priority: above everything else queued, including the PTR comparisons.**
Additive to `RESEARCH_PLAN.md`; ranks above `HPC_TASKS_SKETCH_MODE.md`. Nothing
here disturbs runs already in flight.

`docs/design/04-roadmap.md` lists σ_eff — per-site capture efficiency noise in a
real 2bRAD library — as **open question R1, "the single largest open risk"**, and
gives its resolution path as *"P2 technical replicates"*: a wet-lab experiment
that has not been run. Every route-B claim in the paper rests on assumed values
of 0.3–0.6 and **no published per-site depth CV for 2bRAD exists.**

Two paired datasets can measure it from data alone:

| | Sun et al. 2022, Genome Biol | Hou et al. 2025, Front Microbiol |
|---|---|---|
| paired WGS + 2bRAD | faecal, n = 3 | milk 34, maternal faecal 33, meconium 30 |
| same DNA extract | yes | yes |
| enzyme | BcgI | BcgI |
| 2bRAD depth | 13.5 M reads | 13 M raw / 10 M clean |
| WGS depth | 437 M reads (≈ 240 Gb, **verify**) | 200 M+ (faecal) |
| accession | PRJNA689204 | PRJCA030517 |

Same extract, same enzyme — so the wet-lab/in-silico difference is isolable.

---

## 1. Three arms, not two

| arm | reads | landmarks | estimator |
|---|---|---|---|
| **A** | WGS | FracMinHash | Pilea (reference) |
| **B** | WGS | *in-silico* BcgI digest | sk2bGrow V-fit |
| **C** | 2bRAD | *real* BcgI digest | sk2bGrow V-fit |

A vs C confounds two variables at once — different molecules **and** different
estimator. **B vs C is the clean contrast**: one DNA extract, one landmark set,
one estimator; the only difference is whether the tags were cut in a tube or in
software. That difference is σ_eff, and it is worth more than the PTR comparison
itself.

Arm A is the FracMinHash arm of `HPC_TASKS_SKETCH_MODE.md`. The two documents
meet here.

---

## 2. Task 1 — measure σ_eff directly. Do this before any PTR analysis. ★★★

No profiling run, no V-fit, no origin search. Per-anchor counts only.

For each species in each sample, over anchors unique in **both** pipelines:

```
e_i  =  (C_i / Σ_j C_j)  ÷  (B_i / Σ_j B_j)          # sums over that species' anchors
σ_eff = SD(log e_i), Poisson-corrected
```

### Four corrections, in order of how badly each one bites

**(a) Normalise within species, never by library total.** A species' share of a
2bRAD library is its abundance *weighted by its own BcgI site density*, and that
density varies **2.6× across GC** in our own measurement — 420.8 /Mb in
*B. subtilis* (GC 43.5%), 632.3 in *E. coli* (50.8%), 1087.4 in *P. putida*
(61.5%) (`data/density_3genomes.tsv`). Normalising by library totals folds that
systematic factor into the CV. Summing over the species' own anchors cancels it
exactly.

**(b) Subtract the Poisson term — it is not bookkeeping, it is most of the
variance.**

```
Var(log e_i) = σ_eff² + 1/m_C + 1/m_B          m = mean counts per anchor
```

| m/anchor | Poisson SD | observed SD if σ_eff = 0.4 | capture share of variance |
|---:|---:|---:|---:|
| 3 | 0.816 | 0.909 | **19%** |
| 5 | 0.632 | 0.748 | 29% |
| 10 | 0.447 | 0.600 | 44% |
| 20 | 0.316 | 0.510 | **62%** |
| 30 | 0.258 | 0.476 | 71% |
| 50 | 0.200 | 0.447 | 80% |

At 3 counts/anchor you would be inferring 19% of the variance as the difference
between two numbers of similar size — the answer would be dominated by how well
the Poisson model holds, not by capture. **Restrict this task to m ≥ 20.**

**(c) Unique in both pipelines.** Otherwise the estimate absorbs mapping
differences rather than capture differences. Note our exposure is at
`--max-mismatch 2` (see `HPC_TASKS_SKETCH_MODE.md` §1.4 / F3).

**(d) Report it as an upper bound.** Arm B carries its own per-site structure —
WGS mappability, GC bias, coverage waves. B-vs-C discordance is
σ_eff ⊕ (WGS site effects), so the measurement bounds σ_eff from above. That is
the right direction for a risk gate: **if the bound is below 0.8, route B is
safe regardless of how the bound splits.**

### Abundance needed

BcgI on a 3 Mb genome ≈ 1,900 anchors at *E. coli*-like GC (1,262 at GC 43%,
3,262 at GC 62%). For m ≥ 20:

| dataset | 2bRAD share needed | WGS abundance needed |
|---|---:|---:|
| Sun (13.5 M 2bRAD) | q ≥ 0.28% | p ≥ 0.03% |
| Hou (10 M clean 2bRAD) | q ≥ 0.38% | p ≥ 0.20% |

Species above ~0.3% in a faecal metagenome: typically 15–40 per sample. Across
Hou's 33 maternal faecal samples that is **500–1300 species × sample
observations** — amply powered, and the observation unit is species × sample,
not sample.

### Second product, for free

The same computation yields a **per-site capture-efficiency prior**, `e_i` per
anchor. That is the table needed to *correct* the noise it measures. Measure the
risk and produce its mitigation in one pass.

**Decision rule.** σ_eff bound > 0.8 → route-B gain estimates revise downward,
and we learn it on day one rather than after a full benchmark. Bound < 0.5 → the
assumed 0.3–0.6 range is vindicated and the roadmap's P2 wet-lab experiment can
be descoped.

---

## 3. Task 2 — three-arm PTR comparison, Sun as primary ★★★

**Sun's 3 samples are the primary validation set** despite n = 3, because they
are the only ones where the reference arm is credible across the abundance
range. Hou's 33 faecal samples are the scale-replication set, restricted to
where arm A can report.

### Where Pilea's gate actually bites — bracketed from our own runs

`data/pilea_gate_statistics.tsv`, and the fact that Pilea-with-defaults returned
n = 16 at 10× and **n = 0 at 0.5/1/2/5×**:

| nominal depth | Pilea's own coverage | ratio | dispersion | containment |
|---:|---:|---:|---:|---:|
| 2× | 1.97 | 0.98 | 0.541 | 0.682 |
| 5× | 3.97 | 0.79 | 0.867 | 0.929 |
| 8× | 6.41 | 0.80 | 0.819 | 0.884 |
| 10× | 7.50 | 0.75 | 1.039 | 0.977 |

The gate lies in **(3.97, 7.50]** — consistent with a `coverage ≥ 5` rule — so
arm A needs roughly **6.7× nominal** on that species:

| dataset | WGS | arm A floor | arms B/C floor | binding |
|---|---:|---:|---:|---|
| Sun | 240 Gb | p ≥ 0.008% | p ≥ 0.02% | **our arms** |
| Sun, if 65 Gb | 65 Gb | p ≥ 0.031% | p ≥ 0.02% | arm A |
| Hou faecal | 30 Gb | **p ≥ 0.067%** | p ≥ 0.03% | **arm A** |

Two caveats. This bracket comes from one clean genome; in a community, strain
heterogeneity degrades containment, so treat these as optimistic floors and
re-derive the gate empirically per dataset. And **verify Sun's depth from SRA
metadata before relying on it** — 437 M reads × 150 bp is 65 Gb, not 240 Gb; the
240 Gb figure implies ~550 bp/read, which needs an explanation (2×250 paired?).
The Sun-favoured conclusion survives either way, but the threshold moves 4×.

### Correction to the "≥ 3 reads per anchor" criterion

That criterion is roughly **4× stricter than the method needs**, because the
V-fit consumes *windows*, not anchors. Measured on the Zheng grid
(`data/exemplar_M1_2x_per_enzyme.tsv`), BcgI at 2× nominal:
**mean_rate = 1.64 counts per anchor**, 2,872 anchors in 29 windows ≈ 99 anchors
per window ≈ 162 counts per window. At 1× it is ~0.8 counts per anchor — and
BcgI alone still reaches r = 0.956 there.

So for **PTR** the floor is ~0.8–1.6 counts/anchor, i.e. q ≥ 0.02–0.03%, not
0.045%. For **σ_eff** the floor really is m ≥ 20, because σ_eff is a per-anchor
quantity and cannot be averaged into windows. Two tasks, two thresholds; they
must not share one.

Milk and meconium are still excluded: low biomass, few species, PTR not
estimable.

---

## 4. Task 3 — Hou, scale replication ★★

33 maternal faecal samples, species above the arm-A floor (p ≥ 0.067%, re-derived
empirically). Purpose is reproducibility of the Sun result at scale and across
subjects, not a second primary result.

---

## 5. Task 4 — report agreement, not correlation ★★

Faecal communities are mostly slow-growing or in transit; the true log₂PTR range
is plausibly 0–1, against 0.4–1.7 for the *E. coli* media panel. **Over a
compressed range Pearson r measures the range, not the agreement.**

Report **Bland–Altman** (bias + limits of agreement) and **CCC**, with the
observed truth range printed beside every coefficient. This is the same lesson as
"always report slope beside r", one step further: on the Zheng grid r stayed
0.86–0.97 while slope fell 0.95 → 0.21 under fragmentation, and r alone would
have missed it entirely.

---

## 6. The structural limit — and how much of it is real

Both libraries are **BcgI single-enzyme**. So on this data route B has **one
stratum**: no 16 strata, no cross-enzyme fusion, no Cochran's Q. D4 in
`README.md` — *"16 enzymes = 16 independent measurements"* — cannot be validated
here, and neither can the QC built on it.

That much stands. But the accuracy half of the concern does not survive our own
per-enzyme data. **BcgI alone against the 16-enzyme fusion** (`data/per_enzyme_zheng.tsv`):

| depth | r, BcgI alone | slope | mean r² | r, 16-enzyme fusion | gap |
|---:|---:|---:|---:|---:|---:|
| 0.5× | 0.651 | 0.536 | 0.118 | 0.913 | **+0.262** |
| 1× | **0.956** | 0.972 | 0.387 | 0.981 | +0.025 |
| 2× | 0.952 | 0.945 | 0.718 | 0.982 | +0.029 |
| 5× | 0.961 | 1.026 | 0.872 | 0.979 | +0.018 |
| 10× | 0.957 | 1.088 | 0.854 | 0.968 | +0.011 |

**Above 1×, BcgI alone loses 0.011–0.029 r.** The single-enzyme catastrophe is
confined to 0.5×.

The "mean per-enzyme r² = −0.009 at 0.5×, individually worthless" figure is an
average over all 16 enzymes, including sparse ones like AloI (350 anchors). BcgI
is the *densest* enzyme in the panel (2,872 anchors) and reaches r = 0.651,
mean r² = 0.118 at 0.5× — poor, but not the panel average.

So state the boundary precisely, before a reviewer states it for us:

> This data validates **σ_eff, and the transfer of the method to real
> communities**. It cannot validate the **multi-enzyme architecture** — a
> BcgI-only library has one stratum, so cross-enzyme fusion and the Q-based QC
> are untested here and still need the roadmap's M3 experiment (one BcgI library
> plus a 3–4 enzyme combination per sample). What the single-enzyme data does
> show is that above 1× per-species depth, one dense enzyme is within 0.03 r of
> the full panel; the panel's contribution is concentrated at 0.5× and in the QC.

That last clause is worth having independently: it tells a reader what the
16 enzymes are actually *for*, which the paper currently asserts rather than
bounds.

---

## 7. Execution order

1. **σ_eff from per-anchor B/C ratios** — cheapest, highest-risk, biggest payoff.
   No PTR machinery. Restricted to m ≥ 20. Ships a capture-efficiency prior.
2. **Sun n = 3, three arms**, species above the empirically re-derived arm-A gate.
3. **Hou n = 33 faecal**, scale replication, same gate.
4. **Bland–Altman + CCC throughout**, truth range reported beside every number.

Ahead of `HPC_TASKS_SKETCH_MODE.md` F1–F7, which stay queued behind these.

## 8. What not to do

- **Do not run PTR before σ_eff.** If σ_eff > 0.8 the PTR comparison is being
  read against a noise floor nobody has measured.
- **Do not use library-relative abundance** to normalise `e_i` — §2(a).
- **Do not report r alone** on compressed-range faecal data — §5.
- **Do not claim the enzyme-panel architecture is validated** by BcgI-only data.
- **Do not quote Sun's 240 Gb** until SRA metadata confirms it.
