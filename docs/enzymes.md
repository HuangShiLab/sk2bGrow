# The 16-enzyme Type IIB panel

Source of truth: [`crates/sk2bgrow-core/src/enzyme.rs`](../crates/sk2bgrow-core/src/enzyme.rs).
`enzyme_idx` is an index into `PANEL` and is persisted inside anchor databases
and count tables — **never reorder the table**.

## Panel

| idx | enzyme | recognition | tag (bp) | up/down flank | palindromic |
|----:|--------|-------------|---------:|---------------|:-----------:|
| 0 | BcgI | CGA-N6-TGC | 32 | 10 / 10 | |
| 1 | AlfI | GCA-N6-TGC | 32 | 10 / 10 | ● |
| 2 | AloI | GAAC-N6-TCC | 27 | 7 / 7 | |
| 3 | BaeI | AC-N4-GTAYC | 28 | 10 / 7 | |
| 4 | BplI | GAG-N5-CTC | 27 | 8 / 8 | ● |
| 5 | BsaXI | AC-N5-CTCC | 27 | 9 / 7 | |
| 6 | BslFI | GGGAC | 25 | 0 / 20 | |
| 7 | Bsp24I | GAC-N6-TGG | 27 | 8 / 7 | |
| 8 | CjeI | CCA-N6-GT | 28 | 9 / 8 | |
| 9 | CjePI | CCA-N7-TC | 27 | 8 / 7 | |
| 10 | CspCI | CAA-N5-GTGG | 33 | 11 / 10 | |
| 11 | FalI | AAG-N5-CTT | 27 | 8 / 8 | ● |
| 12 | HaeIV | GAY-N5-RTC | 27 | 8 / 8 | ● |
| 13 | Hin4I | GAY-N5-VTC | 27 | 8 / 8 | |
| 14 | PpiI | GAAC-N5-CTC | 28 | 8 / 8 | |
| 15 | PsrI | GAAC-N6-TAC | 27 | 7 / 7 | |

Invariants enforced by tests in `enzyme.rs`:

* `up_flank + len(pattern) + down_flank == tag_len` for every entry;
* `idx` equals the position in `PANEL`;
* the `display` string expands to exactly the stored `pattern`;
* palindromic entries have `up_flank == down_flank` (see below).

## Two panel properties that change how the code must behave

### 1. Four enzymes are palindromic

`AlfI`, `BplI`, `FalI` and `HaeIV` have recognition patterns equal to their own
reverse complement. Two consequences, both load-bearing:

**Digestion must not double-count them.** A palindromic motif is found by both
the forward and the reverse-complement scan at the same coordinate. `digest.rs`
scans only the forward strand for these, which is the deduplication the design
report calls for in §4.1.

**Counting must try both tag orientations.** The excised duplex is symmetric, but
the *flanks* are not palindromic, so the two strand readings of the tag are
different byte strings. A read arriving reverse-complemented yields the reverse
complement of the stored tag, and the motif scan has no way to tell — the motif
looks identical either way. `AnchorIndex::lookup` therefore probes both
orientations. Without that, these four enzymes lose roughly half their reads,
which would look like a coverage deficit and would also make the cross-enzyme
consistency test fire spuriously.

This is also why the flanks must be symmetric for a palindromic enzyme: an
asymmetric split would make the two strand readings span different intervals,
which is physically impossible for a symmetric cut. REBASE lists HaeIV as
`(7/13)GAY N5 RTC(14/9)`; those coordinates describe the staggered single-strand
nicks that create the overhangs, not the duplex core, so this table splits its
16 bp of flank evenly.

### 2. HaeIV's site is a strict subset of Hin4I's

```
HaeIV   GAY-N5-RTC      R = A/G
Hin4I   GAY-N5-VTC      V = A/C/G      ⊃ R
```

Both have an 11 bp pattern and 8/8 flanks, so **every HaeIV site is also a Hin4I
site, with a byte-identical tag**. The report's own density table is consistent
with this (E. coli: HaeIV 745/Mb < Hin4I 1650/Mb).

Consequences:

* A read tag at such a locus matches two anchors. `count.rs` classifies this as
  `tag_multi_enzyme` — one locus, two enzyme strata — and credits both, rather
  than discarding it as an ambiguous multimapper. Discarding cost ~18 % of all
  tags in the integration test before this was handled.
* **The HaeIV and Hin4I strata are not independent** over the HaeIV subset. The
  cross-enzyme χ² consistency test in `fusion.py` treats enzymes as independent
  channels, so agreement between these two is partly structural rather than
  evidential. With 16 enzymes one correlated pair is a minor effect, but a run
  restricted to `--enzymes HaeIV,Hin4I` would produce a consistency p-value that
  means nothing.

## Reconciliation checklist against `bsyn::enzyme`

This table is a vendored copy so the workspace builds standalone. When Syn2b is
brought in as a dependency, check in this order:

1. **Recognition patterns** — these drive anchor *density* and therefore every
   quantitative claim. They must match exactly.
2. **Flank offsets** — these drive tag *sequence*. Anchor coordinates and
   densities are unaffected (everything downstream keys on `site_start`), so a
   one-base difference shifts tags without moving statistics. Still worth
   aligning, because tags are what reads are matched against.
3. **`BslFI`** is the entry to check first. It is listed here with the panel's
   `GGGAC` / 25 bp tag, but `GGGAC(10/14)` is a Type IIS cut — one-sided, which
   is why its flanks are 0/20 rather than balanced. On a uniform-composition
   sequence `GGGAC` predicts ~2 000 sites/Mb, whereas the report measures 576/Mb
   in E. coli; some of that gap is genuine motif under-representation, but the
   size of it is worth confirming against the upstream registry.

## Reproducing the report's density table

```bash
sk2bgrow digest genomes/GCF_000005845.2_ASM584v2_genomic.fna --enzymes all
```

Prints sites and sites/Mb per enzyme, the merged 16-enzyme union, mean spacing
and the worst-case gap — the columns of report §4.1. Densities depend only on the
recognition patterns, so this is the direct check on item 1 above.
