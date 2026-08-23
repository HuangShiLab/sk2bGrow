# sk2bGrow

Estimating bacterial replication rates (PTR) from **deterministic 2bRAD/TGT
anchors** instead of a random k-mer sketch.

Sixteen Type IIB restriction enzymes cut a genome at motif-defined positions. The
resulting tags are a sketch — but unlike FracMinHash, one whose loci are *known
in advance, identical across every sample, and grouped into sixteen independent
strata*. That difference is what this project is built on.

```
  genomes/*.fna ──digest──▶ TGT v2 ──index──▶ anchor db        (offline, once)
                                                  │
  reads *.fq.gz ─────────── count ──────────────▶ counts.tsv   (Rust, per sample)
                                                  │
                            python/sk2bgrow ─────▶ output.tsv  (PTR + QC)
```

## Why

[Pilea](https://doi.org/10.1186/s40168-025-02268-7) established alignment-free
PTR estimation and set the accuracy bar (r = 0.976 against measured growth
rates). Its remaining weaknesses are structural, and they all trace to one
choice — a random, position-less, unstratified sketch:

| | Pilea | sk2bGrow |
|---|---|---|
| **D1** low-coverage truncation bias | needs ≳5× | more observations per window; usable near 1× |
| **D2** no window guarantee | geometric spacing, unbounded gaps | motif-defined; worst gap auditable at build time |
| **D3** sorted regression rides on extreme order statistics | RANSAC + Tukey patches | anchors have coordinates → fit the V directly |
| **D4** no real replicates | bootstrap over mixture components | 16 enzymes = 16 independent measurements |
| **D6** GC correction is global and post hoc | one loess for everything | per-enzyme loess at anchor resolution |
| **D8** multi-fork profiles are non-linear | linear fit | segmented fit, selected by BIC |

The full argument, with the simulations behind it, is in
[`docs/design/`](docs/design/).

## Install

```bash
cargo build --release --workspace
pip install -e python/
```

The binary lands at `target/release/sk2bgrow`. The Rust layer runs standalone
with `--no-stats`; the Python layer runs standalone against any count table.

## Use

```bash
# 1. build the anchor database (once per reference set)
sk2bgrow index genomes/*.fna -a taxonomy.tsv -o db --enzymes all

# 2. check it before trusting it — density, blind spots, thin enzymes
sk2bgrow audit db/ -o audit.html

# 3. profile samples: route A, shotgun metagenome reads
sk2bgrow profile reads/*.fq.gz -d db -o out/

# 3b. route B, real 2bRAD reads — digestion already happened at the bench
sk2bgrow profile 2brad/*.fq.gz -d db -o out/ --mode 2brad

# 4. compare across samples
sk2bgrow dynamics out/*.tsv -o delta_ptr.tsv --baseline T0
```

Fragmented MAG? Give it a coordinate system first:

```bash
sk2bgrow scaffold mag.fna -d db -r close_relative -o mag.tgt
```

`output.tsv` keeps Pilea's column names (`coverage`, `dispersion`, `fraction`,
`containment`, `PTR`, `log2(PTR)`) so existing benchmark scripts work unchanged,
and appends `enzyme_consistency`, `n_anchors` and `ori_confidence`.

## Layout

```
crates/sk2bgrow-core/   digestion, TGT v2, anchor db, counting, EM, ori, scaffold
crates/sk2bgrow-cli/    the `sk2bgrow` binary
python/sk2bgrow/        ZTP/NB rates, GC correction, V-shape fit, fusion, dynamics
docs/                   design notes, data formats, CLI reference
benches/                A/B benchmark protocol against Pilea
tests/                  Rust integration + Python pytest suites
```

The two halves meet at a **file**, not an FFI boundary: Rust writes a per-anchor
count table, Python reads it. Either side can be rewritten without touching the
other, and every intermediate is inspectable with `head`.

## Reproduce the design report's simulations

```bash
python -m sk2bgrow.cli simulate a    # route A: union vs random sketch vs single enzyme
python -m sk2bgrow.cli simulate b    # route B: 16 enzymes at a fixed read budget
```

Route A reproduces the ordering `union < random sketch < single enzyme` in RMSE at
every coverage, and the single-enzyme collapse below 2×. Route B reproduces the
fixed-budget result: spreading reads over 16 enzymes beats concentrating them on
one, because window averaging quenches site-efficiency noise as √n.

Absolute RMSE depends on using real digested coordinates rather than the built-in
synthetic anchor sets — pass `sk2bgrow index --write-tgt` output through
`simulate.anchors_from_digest()` for that.

## Test

```bash
make test          # cargo test --workspace && pytest
```

## Status

Milestone **M1** of the plan in [`docs/design/04-roadmap.md`](docs/design/04-roadmap.md):
the full pipeline runs end to end and recovers planted PTR values in simulation.
82 Rust tests and 106 Python tests, plus `scripts/smoke.sh`, which builds a
1.5 Mb genome with a planted gradient and checks the whole stack recovers it —
currently log2(PTR) 0.98 against a planted 1.0, origin within 1.3 kb, all 15
enzymes agreeing (I² = 0).

M2 — the A/B benchmark against Pilea on the Zheng E. coli dataset — is the next
gate, and the honest one: if sk2bGrow does not beat the Pilea baseline at 1×
subsampling, the premise needs rethinking rather than more engineering. See
[`benches/README.md`](benches/README.md).

### Things the design documents did not anticipate

Found by running the pipeline, not by reasoning about it. Each is documented
where it lives; collected in [`docs/design/00-overview.md`](docs/design/00-overview.md).

* **Palindromic enzymes need two-orientation matching.** AlfI, BplI, FalI and
  HaeIV read identically on either strand, so a reverse-complemented read yields
  the reverse complement of the stored tag and the motif scan cannot tell.
  Without probing both, four of sixteen strata silently lose half their reads.
* **HaeIV's recognition site is a strict subset of Hin4I's** (`R = A/G ⊂ V = A/C/G`,
  identical flanks), so every HaeIV site is a byte-identical Hin4I tag. Those two
  strata are not independent — a caveat for the χ² test, and ~18 % of all tags
  were being discarded as "ambiguous" before it was handled.
* **Shared anchors must keep their counts** through the counter, or the EM
  reassignment has nothing to reassign.
* **A fixed window size drops the sparse enzymes.** At 100 anchors/window, PpiI,
  BplI and PsrI get 3–4 windows even on E. coli, which defeats a 16-enzyme design.
* **The origin belongs to the genome, not the enzyme.** Searching it per enzyme
  turns search divergence into apparent enzyme disagreement.

The last three share a root: in a stratified design, anything fitted *per
stratum* adds between-stratum variance when it fits noise — and the cross-enzyme
consistency test, this project's headline QC signal, is exactly what detects
that. It flagged the pipeline's own over-fitting before it ever saw real data.

## License

MIT.
