# Architecture

## Layering

| layer | role | language | reuse source |
|---|---|---|---|
| anchor index | in silico digestion, TGT IO, anchor library | Rust | Syn2b (`bsyn`), Syn2bANI |
| counting | reads → anchor counts, shared-anchor EM | Rust | Syn2bANI `tag_matcher.rs`; Pilea `kmc.pyx`/`io.pyx` |
| statistics | ZTP/NB rates, GC correction, V-fit, fusion | Python | Pilea `ztp.py`, `profile.py`, `sketch.py` |
| CLI | subcommands, scheduling, reporting | Rust (clap) + Python (argparse) | Pilea `cli.py`, Syn2b `main.rs` |

The split is not aesthetic. Digestion and counting are throughput-bound and have
existing Rust implementations to mirror. The statistical layer is the part that
will change most often during validation, and it needs `scipy`, `statsmodels` and
a REPL. Putting the boundary between them keeps each side in the language that
makes it cheap.

## The interface is a file

Rust writes a per-anchor count table; Python reads it. No shared memory, no FFI,
no build coupling.

* Either half runs standalone — `--no-stats` stops after counting;
  `python -m sk2bgrow.cli profile` works on any count table.
* Every intermediate is inspectable with `head`.
* Re-running the statistics after a model change costs seconds, not a re-count.

The contract is specified in [`03-data-formats.md`](03-data-formats.md) and
enforced by `sk2bgrow_core::count::COUNT_TABLE_HEADER` on one side and
`sk2bgrow.io.COUNT_COLUMNS` on the other.

## Module map

### `crates/sk2bgrow-core`

| module | role | notes |
|---|---|---|
| `seq` | IUPAC, revcomp, 2-bit packing, bounded Hamming | vendored from `bsyn::seq` |
| `enzyme` | the 16-enzyme panel | `enzyme_idx` is persisted — never reorder |
| `digest` | in silico digestion, density report | palindrome dedup lives here |
| `tgt` | TGT v2 binary (48 B/record) and text | `bsyn::tgt` layout, plus `flags` and `local_gc` |
| `fasta` | FASTA/FASTQ streaming, transparent gzip | `MultiGzDecoder`, not `GzDecoder` |
| `anchor_db` | anchor library, uniqueness masking, persistence | masking computed once at build time |
| `window` | equal-anchor and fixed-bp windows | never spans a contig |
| `count` | motif-seeded matching with a mismatch budget | see below |
| `em` | shared-anchor reassignment | abundance from unique anchors only |
| `ori` | annotation table + coarse circular grid search | refined fit is in Python |
| `scaffold` | order and orient MAG contigs against a reference | `bsyn scaffold` |

### `python/sk2bgrow`

| module | role | replaces |
|---|---|---|
| `io` | interface files | — |
| `ztp` | ZTP mixtures (EM + BIC) and a zero-truncated NB branch | Pilea `ztp.py` |
| `gc_bias` | per-enzyme loess at anchor resolution, Tukey fences | Pilea `profile.py` (correction) |
| `fit` | V-shape MLE on coordinates; sorted+RANSAC for parity | Pilea `profile.py` (fitting) |
| `fusion` | inverse-variance fusion + Cochran's Q | **new** — the core contribution |
| `dynamics` | ΔPTR, anchor × sample matrix, trend tests | **new** |
| `report` | `output.tsv` and QC figures | Pilea output schema |
| `simulate` | Monte-Carlo harness for report §5 | **new** |

## Three implementation choices worth stating

### Counting scans for motifs instead of k-merising reads

A 2bRAD tag always contains its enzyme's recognition site at a known offset. So
rather than hashing every k-mer of a read, the counter scans the read for the 16
motifs and extracts the implied tag span at each hit. This makes matching
`O(read_len × n_enzymes)` motif tests plus one hash lookup per hit, and — the
real win — it makes route A and route B *the same code path*: in route B the read
simply is the tag, motif included.

A tag counts only when it lies wholly inside the read. For 150 bp reads and a
32 bp tag that keeps (150−32+1)/150 ≈ 0.79 of local depth — the same 0.8 factor
the design report uses for Pilea's k=31 sketch, so the two methods are compared
at matched effective depth rather than matched nominal coverage.

Mismatch tolerance uses pigeonhole seeding: with a budget of *m*, a tag splits
into *m+1* contiguous seeds, at least one of which must survive intact.
Candidates are then verified by full Hamming distance.

### Windows are cut per enzyme, in the Python layer

The Rust layer emits a `window_id` computed over the anchor union — that is the
Pilea-parity path, and it is what `--use-rust-windows` selects. The default path
re-cuts windows inside each `(genome, enzyme)` series, because that is what makes
each enzyme an independent channel for fusion. Windowing policy is a statistical
choice, so it lives with the statistics.

Either way a window never spans a contig boundary: two anchors on different
contigs have no defined genomic distance, and pooling them fabricates a
coordinate.

### Uniqueness masking is computed once, at build time

`UNIQUE_IN_GENOME` and `UNIQUE_ACROSS_DB` are properties of the whole database,
so `recompute_uniqueness` runs after assembly and after any mutation. Analysis
then filters for free. The two mask reasons are kept separate because they mean
different things: a multi-copy anchor is unusable for coverage, while a
cross-genome shared anchor is unusable for coverage *but is exactly the input the
EM needs*.

## Dependencies

**Rust** — `clap`, `rayon`, `serde` + `bincode`, `flate2`, `anyhow`, `thiserror`.
A `bsyn` cargo feature is reserved for swapping the vendored digest/TGT
implementation for upstream Syn2b.

**Python** — `numpy`, `scipy`, `pandas`, `statsmodels`, `pyarrow`; `matplotlib`
optional (missing it degrades to "no figures", never to an error).
