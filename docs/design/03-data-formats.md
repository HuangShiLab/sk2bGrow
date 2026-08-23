# Data formats

Three on-disk contracts. Changing any of them is a breaking change.

## TGT v2 — a digested genome

Ordered tags + inter-tag gaps + contig metadata. Reuse target: `bsyn::tgt`.

### Binary

```
magic      "TGT2"          4 B
version    u32 = 2         4 B
meta_len   u32             4 B
meta       JSON            meta_len B   (genome name, contig table)
n_records  u64             8 B
records    48 B each
```

Record layout (little-endian):

| offset | size | field | |
|---:|---:|---|---|
| 0 | 8 | `tag_hash` | strand-canonical hash of the tag |
| 8 | 8 | `position` | 0-based start of the **recognition site** |
| 16 | 4 | `gap` | bp since the previous tag on this contig; 0 for the first |
| 20 | 2 | `contig_id` | |
| 22 | 1 | `enzyme_idx` | index into `enzyme::PANEL` |
| 23 | 1 | `pattern` | which enzyme pattern matched: 0 as written, 1 its rc reading |
| 24 | 1 | `tag_len` | bp |
| 25 | 1 | `flags` | see below |
| 26 | 1 | `local_gc` | quantised ±250 bp GC; 255 = undefined |
| 27 | 1 | — | pad |
| 28 | 12 | `tag_2bit` | packed tag, 4 bases/byte, up to 48 bp |
| 40 | 8 | — | reserved |

The stored coordinate is the **tag window start**, matching Fast2bRAD-M, and the
tag is the forward-strand window at that position — nothing is
reverse-complemented at extraction. `pattern` records which strand reading
matched; it is not a genomic strand.

`gap` restarts at each contig boundary. A gap spanning a boundary is not a
genomic distance, and treating it as one inflates the blind-spot statistics that
`sk2bgrow audit` reports.

### Text

```
#TGT2	<genome>
#contig	<id>	<name>	<length>	<offset>	<kind>
#columns	tag	contig_id	position	strand	gap	flags	local_gc
BcgI:ACGTACGT…	0	4193	+	312	3	102
```

Grouped by enzyme, in the `Enzyme:SEQUENCE` form Syn2b uses. Round-trips
everything except `tag_hash`, which is recomputed from the sequence.

## Anchor database — a directory

```
db/
├── manifest.json     build parameters, genome table, anchor count
├── anchors.bin       bincode: (Vec<Anchor>, Vec<[u8; 12]>)
└── tgt/              optional per-genome text dumps (--write-tgt)
```

`manifest.json` is deliberately human-readable: it is the contract the Python
layer reads. `anchors.bin` holds the anchors and their packed tags as two
parallel vectors; a length mismatch between them is a load-time error rather than
a silent misalignment.

Anchors are sorted by `(genome_id, contig_id, position, enzyme_idx)`, so each
genome's slice is contiguous and `AnchorDb::genome_range` is a binary search.

### `Anchor` (32 B)

```rust
seq_hash:   u64   // strand-canonical tag hash
genome_id:  u32
contig_id:  u16
position:   u64   // recognition site start
enzyme_idx: u8    // index into PANEL
strand:     u8    // pattern index (0 as written, 1 rc reading)
flags:      u8
local_gc:   u8    // quantised, 255 = undefined
```

Extends the Syn2b tag layout with `flags` and `local_gc`. `tag_len` is *derived*
from `enzyme_idx` rather than stored.

### Flags

| bit | name | meaning |
|---:|---|---|
| 0 | `UNIQUE_IN_GENOME` | tag occurs once in its own genome |
| 1 | `UNIQUE_ACROSS_DB` | tag occurs in exactly one genome |
| 2 | `MASKED_MULTICOPY` | repeated within its genome |
| 3 | `MASKED_SHARED` | shared with another genome — **still counted**, for the EM |
| 4 | `NON_CHROMOSOMAL` | plasmid or similar: no ori-ter gradient |
| 5 | `GC_UNDEFINED` | the ±250 bp window was all N |

`USABLE_MASK = MASKED_MULTICOPY | MASKED_SHARED | NON_CHROMOSOMAL`. An anchor is
usable for coverage modelling when none of those are set. Bits 0–1 and 2–3 encode
the same facts from opposite directions on purpose: "unique" drives inclusion,
"masked" carries the *reason* for exclusion into the report.

## Count table — the Rust ↔ Python interface

TSV, one row per anchor per sample:

```
sample  genome_id  genome  contig_id  position  global_position
enzyme  strand  flags  local_gc  window_id  count
```

Defined by `sk2bgrow_core::count::COUNT_TABLE_HEADER` and
`sk2bgrow.io.COUNT_COLUMNS`; a Python-side test asserts they agree.

Notes:

* `global_position` = `contig.offset + position`. Meaningful only for a closed
  chromosome or a scaffolded MAG — it is the V-fit's x-axis.
* `local_gc` is quantised: 0–200 maps 0–100 % in 0.5 % steps, **255 means
  undefined**. `io.read_counts` converts it to a `gc` float column with `NaN` for
  the sentinel, so 255 can never be read as 127.5 % GC.
* `window_id` is the Rust union-based window. The Python layer re-windows per
  enzyme by default and uses this only with `--use-rust-windows`.
* Masked anchors are written too, with their flags, so the Python layer can audit
  what was excluded rather than inferring it from a row count.

### Sidecar: `<sample>.stats.json`

Counting diagnostics and the EM result. `containment` for each genome comes from
here.

`resolved_rate` (matched ÷ *extracted* tags) is the reference-distance
diagnostic. It deliberately does not divide by motif hits: with the 16-enzyme
union a 150 bp read spans several anchors and its edge ones are always truncated,
so a motif-hit denominator reports ~73 % for a perfectly healthy run.

## `output.tsv`

Column order is fixed by `report.OUTPUT_COLUMNS`, and a test asserts it.

**Pilea-compatible block** — `sample`, `genome`, `taxonomy`, `coverage`,
`dispersion`, `fraction`, `containment`, `PTR`, `log2(PTR)`.

**sk2bGrow block** — `enzyme_consistency` (Cochran's Q p-value), `n_anchors`,
`ori_confidence`, `se`, `ci_low`, `ci_high`, `n_enzymes`,
`n_enzymes_attempted`, `enzyme_fit_rate`, `enzymes_used`, `enzyme_i2`,
`fusion_model`, `ori`, `method`, `n_windows`, `pass_qc`, `qc_reason`,
`excluded`, `note`.

Two of these carry the cross-enzyme QC, and they cover different failure modes:
`enzyme_consistency` catches enzymes that fit to a *different slope*, while
`enzyme_fit_rate` catches enzymes that produced *no fit at all* — which never
reach the Q statistic. `fusion_model` reads `random` when Q rejected and the
interval was widened to account for real between-enzyme variance.

Written with `%.10g`, not `%.6g`: the `ori` column is a genome coordinate, and six
significant digits round 3 923 883 to 3 923 880.

## `bsyn` interop

The `bsyn` cargo feature is reserved for replacing the vendored `digest`/`tgt`
implementations with upstream Syn2b. The reconciliation checklist —
what must match exactly, what may differ harmlessly — is in
[`../enzymes.md`](../enzymes.md).
