# The 16-enzyme Type IIB panel

**Source of truth: [`Fast2bRAD-M/src/enzymes.rs`](https://github.com/HuangShiLab/Fast2bRAD-M/blob/main/src/enzymes.rs)**,
which derives its patterns from the `@site` regexes in `2bRADExtraction.pl`.
[`crates/sk2bgrow-core/src/enzyme.rs`](../crates/sk2bgrow-core/src/enzyme.rs) is
a transcription of that table, verified two ways: a reverse-complement closure
test in the unit tests, and measured densities on E. coli K-12 MG1655
(NC_000913.3) reproduced below.

`enzyme_idx` is an index into `PANEL` and is persisted inside anchor databases
and count tables — **never reorder the table**.

## The window model

A tag is **a fixed-length window of the forward strand** that satisfies one of
the enzyme's patterns. A pattern is a set of `(offset, IUPAC motif)` anchors
positioned inside that window:

```
BcgI, tag_len 32, pattern 0:  [ACGT]{10} CGA [ACGT]{6} TGC [ACGT]{10}
                               offset 10 ─┘        offset 19 ─┘
```

Most enzymes carry **two** patterns — the recognition motif as it reads on each
strand. Both are tested against the forward strand and the extracted tag is the
forward-strand window either way; **nothing is reverse-complemented during
extraction**. Comparison is instead strand-canonical (the lexicographically
smaller of a sequence and its reverse complement), because a read may be
sequenced from either strand.

Three enzymes need only one pattern, because their whole *window* pattern is its
own reverse complement: **AlfI, BplI, FalI**. HaeIV's recognition core is
palindromic too, but its flanks are 7/9, so its two readings occupy windows
offset by 2 bp and it needs two patterns.

## Panel

| idx | enzyme | tag | pattern 0 | pattern 1 |
|----:|--------|----:|-----------|-----------|
| 0 | BcgI | 32 | (10,CGA)(19,TGC) | (10,GCA)(19,TCG) |
| 1 | AlfI | 32 | (10,GCA)(19,TGC) | — |
| 2 | AloI | 27 | (7,GAAC)(17,TCC) | (7,GGA)(16,GTTC) |
| 3 | BaeI | 28 | (10,AC)(16,GTAYC) | (7,GRTAC)(16,GT) |
| 4 | BplI | 27 | (8,GAG)(16,CTC) | — |
| 5 | BsaXI | 27 | (9,AC)(16,CTCC) | (7,GGAG)(16,GT) |
| 6 | BslFI | 25 | (6,GGGAC) | (14,GTCCC) |
| 7 | Bsp24I | 27 | (8,GAC)(17,TGG) | (7,CCA)(16,GTC) |
| 8 | CjeI | 28 | (8,CCA)(17,GT) | (9,AC)(17,TGG) |
| 9 | CjePI | 27 | (7,CCA)(17,TC) | (8,GA)(17,TGG) |
| 10 | CspCI | 33 | (11,CAA)(19,GTGG) | (10,CCAC)(19,TTG) |
| 11 | FalI | 27 | (8,AAG)(16,CTT) | — |
| 12 | HaeIV | 27 | (7,GAY)(15,RTC) | (9,GAY)(17,RTC) |
| 13 | Hin4I | 27 | (8,GAY)(16,VTC) | (8,GAB)(16,RTC) |
| 14 | PpiI | 27 | (7,GAAC)(16,CTC) | (8,GAG)(16,GTTC) |
| 15 | PsrI | 27 | (7,GAAC)(17,TAC) | (7,GTA)(16,GTTC) |

Invariants enforced by tests:

* every motif fits inside the tag window;
* **reverse-complement closure** — the reverse complement of any pattern is
  itself a pattern of the same enzyme. This is what catches a mis-transcribed
  offset or motif, and it holds for all 16;
* `idx` equals the position in `PANEL`.

## Measured densities, three genomes

`sk2bgrow digest`, against the design report's Table §4.1. Counting distinct tag
**windows** (start positions), pure ACGT, all overlapping occurrences, union of
both patterns — exactly the Perl's semantics.

| enzyme | E. coli | rpt | B. subtilis | rpt | P. putida | rpt |
|---|---:|---:|---:|---:|---:|---:|
| CspCI | 132 | 132 | 66 | 65 | 142 | 142 |
| AloI | 113 | 112 | 104 | 104 | 144 | 143 |
| BsaXI | 212 | 212 | 337 | 336 | 273 | 272 |
| BaeI | 172 | 171 | 116 | 116 | 177 | 175 |
| BcgI | 632 | 632 | 421 | 421 | 1087 | 1085 |
| CjeI | 1962 | 1910 | 996 | 983 | 2065 | 2012 |
| PpiI | 74 | 73 | 74 | 74 | 138 | 137 |
| PsrI | 92 | 92 | 84 | 84 | 93 | 93 |
| BplI | 83 | 83 | 115 | 115 | 119 | 119 |
| FalI | 158 | 158 | 458 | 457 | 209 | 209 |
| Bsp24I | 352 | 351 | 211 | 211 | 471 | 469 |
| **HaeIV** | **1492** | 745 | **1515** | 755 | **841** | 420 |
| CjePI | 1712 | 1701 | 1470 | 1464 | 1665 | 1657 |
| **Hin4I** | **1223** | 1650 | **1380** | 2010 | **839** | 1057 |
| AlfI | 436 | 436 | 238 | 237 | 762 | 760 |
| BslFI | 577 | 576 | 784 | 784 | 690 | 690 |

**46 of 48 cells reproduce the report within 3 %**, most within 0.5 %. Genomes:
`GCF_000005845.2` (4 641 652 bp, GC 50.8 %), `GCF_000009045.1` (4 215 606 bp,
GC 43.5 %), `GCF_000007565.2` (6 181 873 bp, GC 61.5 %) — matching the report's
stated sizes and GC exactly.

Two rows differ, one explained and one not. A third row is worth noting: **CjeI
runs consistently 1.4–2.7 % high** on all three genomes, the only enzyme outside
±1 % besides these two.

### HaeIV — resolved: locus versus window

Ratio to the report is **2.002 / 2.006 / 2.003** — exactly 2.00 on all three
genomes, which is a convention difference, not noise.

HaeIV's core `GAY-N5-RTC` is its own reverse complement, so *both* its patterns
match at every locus, at window offsets `p−7` and `p−9`. Counting **loci**
instead of windows gives ratios of **1.001 / 1.003 / 1.001**. The report states
it deduplicated palindromic double-strand hits and names HaeIV explicitly.

An independent check: searching the 11 bp motif `GAY-N5-RTC` directly reproduces
745 / 755 / 420 at **0.3 % max error**, and a sweep of all 225 IUPAC combinations
at the two degenerate positions ranks `GAY-N5-RTC` first by a wide margin.

HaeIV is the *only* enzyme where the two patterns can hit one locus: it is the
only one whose recognition core is self-complementary while its flanks are not.
AlfI, BplI and FalI have symmetric flanks and so carry a single pattern.

### Hin4I — unresolved

Ratio to the report is **0.741 / 0.687 / 0.794** — *not* constant, which rules
out any single scale factor or dedup convention. Ruled out by direct measurement
on all three genomes:

| hypothesis | E. coli | B. sub | P. put | max err |
|---|---:|---:|---:|---:|
| union of both patterns (the Perl definition) | 1223 | 1380 | 839 | 31 % |
| sum of both patterns, no dedup | 1968 | 2137 | 1259 | 19 % |
| pattern 0 alone | 986 | 1072 | 631 | 47 % |
| best of all 225 IUPAC variants (`GAV-N5-DTC`) | 1541 | 2175 | 1001 | 8.2 % |
| best spacer-length variant | 1973 | 2114 | 1314 | 24 % |
| report Hin4I misaligned from another row | — | — | — | ≥ 25 % |
| **report target** | **1650** | **2010** | **1057** | |

Nothing fits. The report's Hin4I row is **not reproducible** from
`2bRADExtraction.pl` under any convention tested, while every other enzyme
reproduces. Since the definition is confirmed identical across the Perl,
Fast2bRAD-M and this codebase, the discrepancy is in the report's Hin4I figures
rather than in the enzyme definition.

One structural fact worth recording: **Hin4I's two patterns intersect in exactly
the HaeIV site set** — 3 462 / 3 193 / 2 600 loci, matching HaeIV's locus count
on all three genomes. `GA[C/T]-N5-[A/C/G]TC ∩ GA[C/G/T]-N5-[A/G]TC =
GA[C/T]-N5-[A/G]TC`. That is a real relation between the two enzymes, but it does
*not* produce identical tags: HaeIV's windows sit at `p−7`/`p−9` and Hin4I's at
`p−8`, so the extracted 27-mers differ.

### The union row

The report's union (28 381 / 23 928 / 39 451 loci) does not correspond to the
33 bp merge radius it states. Measured: 40 775 / 33 283 / 54 607 distinct
windows, falling to 26 299 / 21 544 / 37 351 when merged at 33 bp. The report's
numbers are best matched by a merge distance of about **17–26 bp**, not 33.

Union max gap is reproduced closely: **1 446 bp** on E. coli against the
report's 1 447.

## Enzyme containment: Bsp24I ⊂ CjePI

The two enzymes' regexes look unrelated because they are written in **opposite
strand orientations**. Bsp24I's first regex must be compared with CjePI's
*second*, not its first. Laid out over the 27 bp window they differ at exactly
one position:

```
window offset:              111111111122222222
                  0123456789012345678901234567

Bsp24I p0         ........GAC......TGG.......     N8 GAC N6 TGG N7
CjePI  p1         ........GA.......TGG.......     N8 GA  N7 TGG N7
                            ^
                            position 10: Bsp24I fixes C, CjePI leaves it free

Bsp24I p1         .......CCA......GTC........     N7 CCA N6 GTC N8
CjePI  p0         .......CCA.......TC........     N7 CCA N7 TC  N8
                                  ^
                                  position 16: Bsp24I fixes G, CjePI leaves it free
```

Or as recognition sites in one common orientation:

```
Bsp24I           GAC NNNNNN TGG     12 bp, 6 bases fixed
CjePI (revcomp)  GA  NNNNNNN TGG    12 bp, 5 bases fixed   [revcomp of CCA-N7-TC]
```

**CjePI is the less specific enzyme.** Bsp24I fixes one base that CjePI leaves
free, so CjePI's site set is strictly larger and contains Bsp24I's. This is an
ordinary nested specificity between two genuinely different enzymes — not a claim
that they are the same enzyme.

So **every Bsp24I tag is a byte-identical CjePI tag** — measured 1 636 / 1 636 on
E. coli K-12. Two consequences the code has to handle:

* **Counting.** A read tag at such a window matches two anchors. That is one
  locus in two enzyme strata, not an ambiguous multimapper, so it is credited to
  both (`tag_multi_enzyme`, 9.2 % of window hits on E. coli).
* **Uniqueness masking.** Within-genome uniqueness is counted over distinct
  **loci**, not anchor rows. Counting rows treated a co-located pair as a 2-copy
  repeat and masked *100 % of Bsp24I* and 24 % of CjePI out of coverage
  modelling; over the whole panel that was 12.4 % of anchors masked instead of
  the correct 3.4 %.
* **Statistics.** Bsp24I carries no information independent of CjePI. Treating
  the panel as 16 independent strata overstates its replication, and the
  cross-enzyme χ² consistency test in `fusion.py` will read the two as agreeing
  for structural rather than evidential reasons. `enzyme::CONTAINED_PAIRS`
  records the relation.

Measured on all three genomes: **100.0 % of Bsp24I windows are CjePI windows** —
1 636/1 636, 891/891, 2 910/2 910.

A second, partial relation: **Bsp24I pattern 0 ⊂ CjeI pattern 1** over their
shared 27-base prefix (CjeI's tag is one base longer, so this needs a
shift-aware comparison that an equal-length test misses). Empirically 48.4 % /
47.4 % / 50.9 % of Bsp24I windows are also CjeI windows — exactly its pattern-0
share.

Other measured co-located groups on E. coli: BsaXI+Hin4I 229, BsaXI+CjePI 76,
AloI+BsaXI+PpiI 55 (a genuine three-way core).

Both relations are recorded in `enzyme::CONTAINMENTS`, and
`enzyme::redundant_enzymes()` returns the enzymes that carry no independent
information. **The panel offers at most ~15 independent strata, not 16.**

## BslFI is Type IIS, not Type IIB

The Perl source comments this enzyme `??some question?? single enzyme`. Its
density is not the problem — it reproduces the report at 577/784/690 vs
576/784/690. The issue is biological.

REBASE gives **`GGGAC(10-11/14-15)`, "Type II restriction enzyme, subtype: S"**
(prototype FinI). A contiguous 5 bp site, cut **downstream only**, at variable
positions. It excises nothing.

The two Perl regexes are **not two cuts**. `GTCCC` is the reverse complement of
`GGGAC`, so `N14 GTCCC N6` is the *same window viewed from the other strand* as
`N6 GGGAC N14` — one window shape (6 bp on the motif's 5' side, 14 bp on its 3'
side), two strand views. That is exactly the same fwd/rev pattern-pair structure
every other enzyme in the panel has.

Where BslFI differs is whether the window's *edges* are cuts. Compare:

```
BcgI   (10/12) CGA N6 TGC (12/10)      panel window: N10 CGA N6 TGC N10 = 32 bp

       cut                                          cut
        v                                            v
  5' ---|----------CGAnnnnnnTGC----------|---  3'
  3' -----------|--GCTnnnnnnACG------------|-  5'
        |<-------------- 32 bp ------------->|
  Both edges are cuts; the enzyme excises this fragment (2 nt 3' overhangs).


BslFI  GGGAC(10-11/14-15)               panel window: N6 GGGAC N14 = 25 bp

       |<-- 6 bp -->|GGGAC|<------ 14 bp ------>|
       ^                              ^         ^
   NO CUT HERE                    top cut   bottom cut
  (arbitrary padding)              (+10)      (+14)
  Only the RIGHT edge is a cut. The 4 nt between the two cut positions is the
  5' overhang. There is no upstream cut in the notation at all.
```

So for the fifteen Type IIB enzymes the panel's window *is* the excised
fragment's core, both edges being real cuts. For BslFI the right edge is a real
cut and the left 6 bp is padding chosen to round the window to 25 bp.

* **In silico this is harmless.** The window is deterministic,
  reverse-complement closed, and applied identically to reference and reads — a
  perfectly good marker stratum.
* **At the bench it is not executable as 2bRAD.** A real BslFI digest yields
  kilobase fragments (mean 1.3–1.7 kb across the three genomes; only ~3 % below
  40 bp), so there is no short band to size-select and each fragment has one
  arbitrary end.

Consequence for this project: BslFI's route-B stratum can never be validated
against real 2bRAD data, and it should not appear in a density table a user reads
as "pick the densest enzyme". Recorded as `enzyme::BSLFI_IS_TYPE_IIS`.

Its density is also driven by taxon-specific motif avoidance rather than GC:
observed/iid-expected is 0.28 / 0.62 / 0.20, and P. putida (61.5 % GC) has *fewer*
BslFI sites per Mb (690) than B. subtilis (43.5 % GC, 784) despite a ~4× higher
expectation — the opposite of BcgI, which tracks GC monotonically. That makes it
a poor choice for cross-species normalisation independent of the Type IIS issue.

## Transcription check

A mechanical position-by-position comparison of `2bRADExtraction.pl`,
`Fast2bRAD-M/src/enzymes.rs` and this codebase's `PANEL` — all three expanded to
per-position allowed-base sets over the full tag window:

| check | result |
|---|---|
| enzymes compared | 16 / 16 |
| patterns compared | 30 |
| constrained positions compared | 442 |
| tag-length / pattern-count / offset / motif mismatches | 0 |
| within-enzyme pattern-order mismatches | 0 |
| IUPAC compressions (`[GAC]`→V, `[CTG]`→B, `[CT]`→Y, `[AG]`→R) | correct |

**Zero discrepancies.** The panel ordering differs (this codebase pins BcgI at
index 0, the rest alphabetical; the Perl numbers them 1–16) but that is
intentional and carries no collision risk — `parse_selection` accepts only `all`
or enzyme names, never a numeric index.

One behavioural divergence outside the definition table: this codebase digests
**soft-masked (lowercase) reference regions**, where Fast2bRAD-M silently skips
them. Masking is a repeat annotation rather than a quality call, and a masked
anchor is better handled by the multi-copy flag than by omission — but on a
lowercase-masked reference the two tools will report different anchor sets. See
`seq::normalize_in_place`.

## Reproducing this table

```bash
sk2bgrow digest ecoli.fna --enzymes all
sk2bgrow index ecoli.fna -o db --enzymes all && sk2bgrow audit db/ -o audit.html
```

Genome used: `GCF_000005845.2_ASM584v2` (NC_000913.3), 4 641 652 bp, GC 50.8 %.
