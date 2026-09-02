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

[Pilea](https://doi.org/10.1186/s40168-026-02374-0) established alignment-free
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
86 Rust tests and 106 Python tests, plus `scripts/smoke.sh`, which builds a
1.5 Mb genome with a planted gradient and checks the whole stack recovers it —
currently log2(PTR) 0.986 against a planted 1.0, origin within 417 bp of the
planted 150 kb, all 16 enzymes fitting and agreeing (I² = 0).

M2 — the A/B benchmark against Pilea on the Zheng E. coli dataset — is the next
gate, and the honest one: if sk2bGrow does not beat the Pilea baseline at 1×
subsampling, the premise needs rethinking rather than more engineering. See
[`benches/README.md`](benches/README.md).

### Findings, and corrections to earlier ones

The enzyme table is transcribed from
[Fast2bRAD-M](https://github.com/HuangShiLab/Fast2bRAD-M) and verified against
three real genomes spanning 43.5–61.5 % GC (E. coli K-12, B. subtilis 168,
P. putida KT2440): **46 of 48 enzyme × genome cells reproduce the design report's
Table §4.1 within 3 %**, and the union's worst-case gap matches to 1 bp (1 446 vs
1 447). A position-by-position diff of all 16 enzymes, 30 patterns and 442
constrained positions against both the Perl and Fast2bRAD-M found **zero**
discrepancies. See [`docs/enzymes.md`](docs/enzymes.md).

**Real, and confirmed on real data:**

* **Bsp24I ⊂ CjePI.** Both Bsp24I patterns are strict refinements of a CjePI
  pattern at the same tag length and offsets, so every Bsp24I tag is a
  byte-identical CjePI tag — 100 % on all three genomes (1 636/1 636, 891/891,
  2 910/2 910). A second partial relation, Bsp24I pattern 0 ⊂ CjeI pattern 1,
  covers a further ~48–51 %. The panel offers at most ~15 independent strata,
  not 16.
* **BslFI is Type IIS, not Type IIB** — `GGGAC(10/14)`, cut on one side only, so
  it excises nothing. Fine as an in-silico marker stratum; not executable as a
  bench 2bRAD protocol (a real digest gives 1.3–1.7 kb fragments, ~3 % under
  40 bp). The original Perl flags it `??some question??`.
* **The report's Hin4I density is not reproducible** from `2bRADExtraction.pl`
  under any convention tested, while 46 of 48 enzyme × genome cells reproduce
  within 3 %. Its ratio to the report is 0.741 / 0.687 / 0.794 — not even
  constant, so no scale factor explains it.
* **Uniqueness must be counted over loci, not anchor rows.** Two enzymes claiming
  one window is a locus seen twice, not a duplicated locus. Counting rows masked
  *100 % of Bsp24I* and 24 % of CjePI out of coverage modelling — 12.4 % of all
  anchors instead of the correct 3.4 %.
* **Shared anchors must keep their counts** through the counter, or the EM
  reassignment has nothing to reassign.
* **A fixed window size drops the sparse enzymes.** At 100 anchors/window, PpiI,
  BplI and PsrI get 3–4 windows even on E. coli.
* **The origin belongs to the genome, not the enzyme.** Searching it per enzyme
  turns search divergence into apparent enzyme disagreement.

**Withdrawn.** Two earlier claims did not survive checking against the reference
implementation, and both traced to the same root cause — five of my sixteen flank
definitions (CjeI, PpiI, HaeIV, CjePI, BslFI) were wrong, and all measurements
had been taken on synthetic random-sequence genomes rather than a real one:

* *"HaeIV's tags are byte-identical to Hin4I's."* The recognition sites do
  overlap — HaeIV's site is exactly the intersection of Hin4I's two patterns —
  but with correct flanks (HaeIV 7/9, Hin4I 8/8) the tag *windows* differ, so the
  tags are distinct sequences. The real containment is Bsp24I ⊂ CjePI.
* *"Palindromic enzymes silently lose half their reads."* This described a bug in
  an extraction model that reverse-complemented tags and compared them
  byte-exactly. The reference model never reverse-complements during extraction —
  the tag is the forward-strand window, and both orientations are covered by the
  enzyme having two patterns — so the situation does not arise. Matching is
  strand-canonical, as in Fast2bRAD-M.

## License

MIT.
# sk2bGrow
