# Zheng 2020 E. coli benchmark (P1 / M2 gate)

Head-to-head against [Pilea](https://github.com/xinehc/pilea) on real data with
an independent ground truth. This is the roadmap's decisive gate: if sk2bGrow
does not beat the Pilea baseline in the 1x band, the premise needs rethinking
rather than more engineering.

## Data

[Zheng et al. 2020, *Nature Microbiology*](https://doi.org/10.1038/s41564-020-0717-x),
"General quantitative relations linking cell growth and the cell cycle in
*Escherichia coli*" — BioProject **PRJNA615952**, 45 WGS runs of E. coli K-12
MG1655 across **16 growth media plus a stationary (RUN_OUT) control**, ~490x
each, most media in triplicate.

`fetch.sh` pulls only the first 600 000 reads per run (~19.4x) rather than the
full 55 GB, because everything is subsampled to <=10x anyway. SRA preserves
original flowcell read order, which is random with respect to genome position —
fine for a coverage profile, though not a substitute for true random sampling if
you care about tile-level bias. `fetch3.sh` retries from sibling replicates when
ENA throttles a particular path, which happens often and unpredictably.

## Ground truth

`growth_rates.tsv`, extracted from the paper's supplementary source data
(`41564_2020_717_MOESM3_ESM.xlsx`). Two different targets, and they are **not**
equally independent:

| target | what it is | independent? |
|---|---|---|
| `growth_rate` (lambda, /h) | steady-state growth rate, 0.40–1.72 across the 16 sequenced media | **yes** — measured by optical density / microscopy |
| `pred_log2ptr` = lambda*C/ln2 | theory-predicted log2(PTR) | **no** — the authors derived C from *this same sequencing* by marker-frequency analysis |

So correlation against `growth_rate` is the real test, and it is what Pilea's
reported r = 0.9764 refers to. RMSE against `pred_log2ptr` is a method
comparison against an independent analysis of the same reads — informative about
magnitude, but circular if quoted as accuracy.

The RUN_OUT control is a negative control: stationary phase, expected
log2(PTR) ~ 0.

## Arms

Four, so that the sketch and the estimator are separately attributable:

| arm | sketch | estimator |
|---|---|---|
| A | 16-enzyme anchors | adaptive windows + V-shape fit |
| B | 16-enzyme anchors | 25 kb windows + sorted/RANSAC (Pilea parity) |
| C default | FracMinHash (Pilea) | Pilea, shipped defaults |
| C relaxed | FracMinHash (Pilea) | Pilea, `-x 0 -z 0 -c 0` |

**A vs B** isolates the estimator. **B vs C** isolates the sketch. A vs C is the
end-to-end number.

Pilea is run twice because its default `--min-cove 5` refuses to report below
5x. Without the relaxed arm the low-coverage comparison degenerates into
"estimate vs no estimate", which is not a comparison. Both are reported; the
default arm is what a user actually gets.

## Running

```bash
./fetch.sh && ./fetch3.sh          # download (~2 GB, retries on throttling)
COVS="0.5 1 2 5 10" ./bench.sh     # arms A and B
./pilea_arm.sh                     # arm C
python3 analyze.py                 # comparison table
```

Requires `pilea` on PATH (`conda create -n pilea -c bioconda pilea`) and a
release build of `sk2bgrow`.

## Reporting

Report per coverage: Pearson r and Spearman rho against growth rate, RMSE
against the predicted log2(PTR), the fitted slope, and **the number of
conditions that yielded an estimate at all** — a tool that reports on 4 of 16
conditions is not comparable to one that reports on 16, and correlation alone
hides that.

Record the sk2bGrow commit, the Pilea version, and the subsampling read counts
with any published number.

## Result (2026-08-24)

sk2bGrow commit at time of run: see `git log`. Pilea v1.3.8 (bioconda).
16 media + stationary control, 600 000 reads/run subsampled to each level.
Full table in `RESULTS.txt`, per-run values in `results_raw.tsv`.

### Correlation with measured growth rate

| coverage | sk2bGrow | Pilea (defaults) | Pilea (gates off) |
|---|---:|---:|---:|
| 0.5x | **0.907** (n=16) | no estimate (n=0) | undefined — every sample returned PTR=1.0 |
| 1x | **0.954** (n=16) | no estimate (n=0) | 0.897 (n=15) |
| 2x | **0.975** | no estimate (n=0) | 0.955 |
| 5x | **0.979** | no estimate (n=0) | 0.953 |
| 10x | 0.974 | 0.972 (n=15) | 0.972 |

Pilea's published figure on this dataset is r = 0.9764 at full depth; it reaches
0.972 here at 10x, so the reimplementation-free comparison is consistent with its
paper.

**The gate is passed.** sk2bGrow is at or above Pilea at every coverage, and the
margin is where the design report predicted it would be — the 1–2x band. At its
shipped defaults (`--min-cove 5`) Pilea returns **no estimate at all** below 10x,
which is defect D1 exactly as described.

### What is actually responsible

| coverage | sk2bGrow (V-shape fit) | same anchors, Pilea-parity estimator |
|---|---:|---:|
| 0.5x | 0.907 | 0.445 |
| 1x | 0.954 | 0.605 |
| 5x | 0.979 | 0.912 |

Same sketch, same reads. **The estimator, not the sketch, produces most of the
gain.** Fitting the V-shape on real coordinates is doing the work — defect D3.
Comparing arm B against Pilea (0.605 vs 0.897 at 1x) shows the deterministic
sketch by itself is *behind* FracMinHash under a sorted-regression estimator;
the anchors only pay off once their coordinates are used. That is a sharper and
less flattering finding than "deterministic anchors are better", and it is the
one the data supports.

### Negative control — RUN_OUT (stationary phase, true log2(PTR) ~ 0)

| arm | 0.5x | 1x | 2x | 5x | 10x |
|---|---:|---:|---:|---:|---:|
| sk2bGrow | 0.113 | 0.045 | 0.100 | 0.081 | 0.054 |
| Pilea-parity estimator | 2.168 | 1.788 | 1.481 | 0.851 | 0.748 |
| Pilea | — | — | — | — | 0.168 |

The sorted-regression estimator reports **log2(PTR) = 2.17 for a non-growing
culture** at 0.5x — it manufactures a gradient out of rank-ordered noise. This is
the clearest single demonstration of D3 in the whole benchmark, and it is why the
coordinate fit matters beyond its correlation score.

### Caveats

1. **Pilea's "gates off" arm runs it outside its supported regime.** Its authors
   gate at 5x deliberately. The relaxed arm exists so the low-coverage
   comparison is not "estimate vs no estimate"; both columns are reported and
   the default column is what a user actually gets.
2. **Slope < 1 at low coverage** for sk2bGrow (0.52 at 0.5x, rising to 0.95 at
   10x). Conditions are ranked correctly but the PTR *range is compressed*.
   Correlation flatters this; RMSE against the predicted log2(PTR)
   (0.51 -> 0.03) shows it plainly. Do not quote r alone.
3. **RMSE is against a non-independent target.** See the ground-truth section.
4. **This is the easiest possible case**: one organism, one strain, a complete
   single-contig reference, no community. It says nothing yet about the
   metagenomic setting, which is where PTR estimation is actually hard.
5. n = 15 for Pilea vs 16 for sk2bGrow at some levels (a few runs were still
   pending when the table was generated).
