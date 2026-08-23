# Benchmarks

The A/B protocol against Pilea. Not a speed benchmark — the question is whether
deterministic anchors beat a random sketch at the same coverage, and by how much
in the 1–2× band where the argument lives.

## Principle

Change **one** thing at a time. Pilea and sk2bGrow differ in both the sketch and
the estimator, so comparing them head to head confounds the two. Every comparison
below fixes one and varies the other.

| comparison | sketch | estimator | isolates |
|---|---|---|---|
| A | random ↔ union | sorted regression (both) | the **sketch** |
| B | union (both) | sorted ↔ V-shape | the **estimator** |
| C | Pilea ↔ sk2bGrow defaults | — | the end-to-end delta |

Comparison A is the one the design report's central claim rests on. Run it first.

`sk2bgrow profile --windowing bp --window-bp 25000` plus
`python -m sk2bgrow.cli profile --method sorted --use-rust-windows` puts
sk2bGrow into full Pilea-parity mode for the estimator side.

## B0 — simulation (no data download)

```bash
python -m sk2bgrow.cli simulate a --reps 150 --output benches/work/route_a.tsv
python -m sk2bgrow.cli simulate b --reps 150 --output benches/work/route_b.tsv
```

Reproduces report §5. Expected: `union_16 < random_sketch < BcgI_single` in RMSE
at every coverage, single-enzyme collapse below 2×, and in route B the
same-budget union beating the deep single enzyme once σ_eff ≳ 0.3.

The built-in anchor sets are synthetic (uniform for the random sketch, even for
the union), which reproduces the *ordering* but not the absolute RMSE. For that,
digest a real genome and feed the coordinates in:

```bash
sk2bgrow index ecoli.fna -o benches/work/db --enzymes all --write-tgt
python - <<'PY'
from sk2bgrow import simulate
union = simulate.anchors_from_digest("benches/work/db/tgt/ecoli.tgt")
bcgi  = simulate.anchors_from_digest("benches/work/db/tgt/ecoli.tgt", enzyme="BcgI")
rand  = simulate.synthetic_anchors(18_600, simulate.ECOLI_LEN, "uniform")
df = simulate.route_a(anchor_sets={"BcgI_single": bcgi, "random_sketch": rand, "union_16": union},
                      n_reps=150)
df.to_csv("benches/work/route_a_real.tsv", sep="\t", index=False)
PY
```

## B1 — anchor density on real genomes (report §4.1)

```bash
sk2bgrow digest genomes/*.fna --enzymes all -o benches/work/density.tsv
```

Targets, from the report:

| genome | GC | union anchors | per 25 kb | max gap |
|---|---:|---:|---:|---:|
| E. coli K-12 MG1655 | 50.8 % | 28 381 | 153 | 1 447 bp |
| B. subtilis 168 | 43.5 % | 23 928 | 142 | ~1.5 kb |
| P. putida | 61.6 % | 39 451 | 160 | ~1.5 kb |

Per-enzyme reference values are in report table §4.1. Densities depend only on
the recognition patterns, so this is the direct check that the vendored enzyme
table matches upstream — see [`../docs/enzymes.md`](../docs/enzymes.md), and
check `BslFI` first.

## B2 — Zheng E. coli, the real gate (P1)

PRJNA615952: 16 growth conditions, >300× coverage, measured steady-state growth
rates 0.4–1.7 h⁻¹. Pilea reaches r = 0.9764.

```bash
# subsample to 0.5, 1, 2, 5, 10x and run both tools on identical inputs
for cov in 0.5 1 2 5 10; do
  seqtk sample reads.fq.gz $(fraction_for $cov) > sub_${cov}x.fq
  sk2bgrow profile sub_${cov}x.fq -d db -o out_${cov}x/
  pilea profile sub_${cov}x.fq -d pilea_db -o pilea_${cov}x/
done
```

Then correlate `log2(PTR)` against the measured growth rates at each coverage.

**Pass:** correlation ≥ Pilea at every coverage, and a clear win at 1×.
**Fail:** the premise needs rethinking, not more engineering. The entire argument
for a deterministic sketch is the low-coverage band; if the gain is not there, it
is not anywhere.

Report per coverage: Pearson r, RMSE against the growth-rate-derived expectation,
the number of conditions yielding an estimate at all, and — for sk2bGrow only —
median `enzyme_consistency`. That last column has no Pilea counterpart and is
worth reporting on its own: it says whether the 16 enzymes agreed, which is
information no single-sketch method can provide.

## B3 — simulated communities (P0)

Pilea's design: 400 samples, 4–32 strains, 4–32×, log2 PTR ~ U[0,2]. Score by L2
distance between estimated and true PTR vectors (Pilea: 10.681; iRep close;
GRiD 39.070).

## Recording results

Write results under `benches/work/` (git-ignored). For anything reported
publicly, record:

* commit hash of sk2bGrow, and the Pilea version;
* exact database build command, including `--enzymes`;
* subsampling seed;
* whether parity mode was used, and which parts of it.

An A/B number without the parity settings attached is not interpretable, because
the sketch and the estimator both changed.
