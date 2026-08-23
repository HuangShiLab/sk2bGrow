# CLI reference

Two-stage shape, deliberately matching Pilea's (`index` once, `profile` per
batch) so a benchmark script can swap tools with minimal edits.

Global flags: `--threads N` (0 = one per core), `--quiet`.

---

## `sk2bgrow index`

Build an anchor database from reference genomes. Run once per reference set.

```bash
sk2bgrow index genomes/*.fna -o db --enzymes all -a taxonomy.tsv --write-tgt
```

| flag | default | |
|---|---|---|
| `-o, --output` | required | database directory |
| `-e, --enzymes` | `all` | `all`, or `BcgI,AlfI,…` |
| `-a, --taxonomy` | — | TSV `genome_name<TAB>lineage` |
| `--ori` | — | TSV `genome<TAB>position[<TAB>confidence][<TAB>source]` |
| `--min-contig-len` | 500 | skip shorter contigs |
| `--gc-flank` | 250 | half-width of the local GC window |
| `--write-tgt` | off | also write per-genome TGT text dumps |

Positional arguments may be files or directories. Genome ids follow the sorted
input order, so rebuilding the same inputs gives byte-identical ids.

Reports how many anchors were masked and why — a high masked fraction means the
reference set contains near-duplicate genomes.

## `sk2bgrow profile`

Count reads and estimate PTR.

```bash
sk2bgrow profile reads/*.fq.gz -d db -o out/
sk2bgrow profile 2brad/*.fq.gz -d db -o out/ --mode 2brad
```

| flag | default | |
|---|---|---|
| `-d, --db` | required | database directory |
| `-o, --output` | required | output directory |
| `--mode` | `wms` | `wms` or `2brad` |
| `-e, --enzymes` | all indexed | must be a subset of the database's |
| `--max-mismatch` | 2 | Hamming budget against the reference tag |
| `--per-file` | off | one sample per file instead of grouping by name |
| `--windowing` | `anchors` | `anchors` (TGT-native) or `bp` (Pilea parity) |
| `--window-anchors` | 100 | anchors per window |
| `--window-bp` | 25000 | bp per window |
| `--no-stats` | off | stop after the count tables |
| `--python` | `python3` | interpreter for the statistics layer |

Read files are grouped into samples by name, with `_R1`/`_R2`/`_1`/`_2` mate
markers stripped, so paired files land in one count table.

**Outputs**

```
out/
├── windows.tsv           union windows (Rust)
├── <sample>.counts.tsv   per-anchor counts        ← Rust/Python interface
├── <sample>.stats.json   counting + EM diagnostics
├── windows.rates.tsv     per-enzyme window rates  (Python)
├── per_enzyme.tsv        one PTR fit per enzyme   (Python)
└── output.tsv            fused PTR + QC           (Python)
```

If the Python layer cannot start, the run fails with the exact command it tried.
"Counts produced, no PTR" is the worst possible outcome for a batch job, so it is
never silent.

## `sk2bgrow dynamics`

```bash
sk2bgrow dynamics out/*.tsv -o delta_ptr.tsv --baseline T0 -m metadata.tsv
```

`--metadata` is a TSV of `sample<TAB>group[<TAB>timepoint]`; `--baseline` names a
sample or a group. With a `timepoint` column, a per-genome trend test is written
alongside.

## `sk2bgrow audit`

Anchor density and blind spots for a database — the build-quality gate.

```bash
sk2bgrow audit db/ -o audit.html --wide-gap 5000
```

`.html` renders a report; any other extension writes TSV. Reports per genome:
usable anchors, anchors per 25 kb, mean/median/p99 spacing, worst gap, count of
gaps beyond `--wide-gap`, and which enzymes are too thin
(`--min-anchors-per-enzyme`). Problems are warnings, not a non-zero exit, so
batch pipelines are not broken by an advisory.

This subcommand is the concrete form of the report's §4.2 argument: with a
deterministic sketch, sparse regions are knowable *before* you analyse anything.

## `sk2bgrow digest`

Digest genomes and print the density table — reproduces report §4.1. Needs no
database, so it is the fastest way to vet a new reference.

```bash
sk2bgrow digest genomes/*.fna --enzymes all -o density.tsv
```

## `sk2bgrow scaffold`

Order and orient draft MAG contigs against a reference, using shared tags.

```bash
sk2bgrow scaffold mag.fna -d db -r ecoli_k12 -o mag.tgt --min-tags 3
```

Writes the scaffolded TGT plus a `.scaffold.json` of placements. Contigs that
cannot be placed are parked past the placed region and reported — the statistics
layer drops their anchors, because a wrong coordinate is worse than a missing one
for a gradient fit.

---

## The Python layer standalone

```bash
python -m sk2bgrow.cli profile out/*.counts.tsv --db db --output out/ --figures
python -m sk2bgrow.cli dynamics out/output.tsv --output delta.tsv --baseline T0
python -m sk2bgrow.cli simulate a --reps 150 --output route_a.tsv
python -m sk2bgrow.cli manifest db/
```

### `profile` flags

| flag | default | |
|---|---|---|
| `--window-anchors` | `auto` | anchors per window; an integer, or `auto` to size each enzyme by its anchor count |
| `--window-cap` | 100 | upper bound on anchors per window in auto mode |
| `--use-rust-windows` | off | group by the Rust `window_id` instead (Pilea-parity windowing) |
| `--count-model` | `auto` | `auto`, `ztp` or `nb` |
| `--method` | `auto` | `auto`, `v_shape` or `sorted` |
| `--per-enzyme-ori` | off | let every enzyme search for its own origin instead of sharing one |
| `--no-gc-correct` | off | skip per-enzyme GC correction |
| `--min-coverage` | 1.0 | QC floor; pass 5 for Pilea-comparable strictness |
| `--min-enzyme-fit-rate` | 0.8 | flag a genome when fewer enzymes than this produced a fit |
| `--alpha` | 0.05 | cross-enzyme consistency threshold |
| `--figures` | off | write QC figures |

`--window-anchors auto` is the default for a reason: at a flat 100 anchors/window
the three sparsest enzymes get 3–4 windows even on E. coli and drop out of the
panel entirely. See [`design/02-algorithm.md`](design/02-algorithm.md).

### Two different "thin enzyme" reports

They mean different things and live in different places:

* **`sk2bgrow audit --min-anchors-per-enzyme`** — a *build-time* property: this
  genome is too small, or too GC-skewed, for that enzyme. Nothing about any
  sample.
* **`--min-enzyme-fit-rate`** (in `output.tsv` as `enzyme_fit_rate`) — a
  *per-sample* property: an enzyme that had enough anchors nonetheless produced
  no usable fit. That is a signal about the sample — a flat profile, a digestion
  failure, a mis-assembled region.

Conflating them would turn "PpiI is sparse on this genome" into a false alarm
about every sample profiled against it.
