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

## Measured densities: E. coli K-12 MG1655

`sk2bgrow digest ecoli.fna --enzymes all`, against the design report's Table §4.1.

| enzyme | tags | per Mb | report | ratio |
|---|---:|---:|---:|---:|
| BcgI | 2 935 | 632 | 632 | 1.00 |
| AlfI | 2 023 | 436 | 436 | 1.00 |
| AloI | 523 | 113 | 112 | 1.01 |
| BaeI | 797 | 172 | 171 | 1.00 |
| BplI | 386 | 83 | 83 | 1.00 |
| BsaXI | 984 | 212 | 212 | 1.00 |
| BslFI | 2 679 | 577 | 576 | 1.00 |
| Bsp24I | 1 636 | 352 | 351 | 1.00 |
| CjeI | 9 108 | 1 962 | 1 910 | 1.03 |
| CjePI | 7 947 | 1 712 | 1 701 | 1.01 |
| CspCI | 613 | 132 | 132 | 1.00 |
| FalI | 735 | 158 | 158 | 1.00 |
| **HaeIV** | 6 924 | **1 492** | 745 | **2.00** |
| **Hin4I** | 5 675 | **1 223** | 1 650 | **0.74** |
| PpiI | 341 | 73 | 73 | 1.01 |
| PsrI | 429 | 92 | 92 | 1.00 |

Union: 40 775 tag windows → **27 054 loci** merged within 33 bp = 146 per 25 kb,
mean spacing 172 bp, **max gap 1 446 bp** (report: 28 381 / 153 / 1 447 bp).

Fourteen of sixteen match to within 3 %, and the worst-case gap matches to 1 bp.
The two that differ do so for understood reasons:

**HaeIV, ratio exactly 2.00.** Its core `GAY-N5-RTC` is palindromic, so *both*
its patterns match at every locus, at window offsets `p−7` and `p−9`. Counting
**windows** gives 1 492/Mb; counting **loci** gives 746/Mb, which is the report's
745. The report states it deduplicated palindromic double-strand hits, naming
HaeIV explicitly. Both conventions are self-consistent — they answer "how many
distinct tag sequences" versus "how many cut sites". This codebase counts
windows, since a window is what a read matches.

**Hin4I, ratio 0.74.** Its two patterns share offsets 8/16, so a locus matching
both collapses to one window; the union is 5 675. Their intersection is
`GA[C/T]-N5-[A/G]TC` — **exactly the HaeIV recognition site**, verified
computationally (3 462 loci, matching HaeIV's locus count). The report's
1 650/Mb corresponds to neither the union nor the sum, and is the one number here
not yet reconciled.

## Enzyme containment: Bsp24I ⊂ CjePI

Both Bsp24I patterns are strict refinements of a CjePI pattern, at the same tag
length and the same offsets:

```
Bsp24I p0  (8,GAC)(17,TGG)   ⊂   CjePI p1  (8,GA )(17,TGG)     GAC at 8 implies GA at 8
Bsp24I p1  (7,CCA)(16,GTC)   ⊂   CjePI p0  (7,CCA)(17,TC )
```

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

Measured co-located window groups on E. coli: Bsp24I+CjePI 844,
Bsp24I+CjeI+CjePI 759, BsaXI+Hin4I 229, BsaXI+CjePI 76, AloI+BsaXI+PpiI 55.

## Reproducing this table

```bash
sk2bgrow digest ecoli.fna --enzymes all
sk2bgrow index ecoli.fna -o db --enzymes all && sk2bgrow audit db/ -o audit.html
```

Genome used: `GCF_000005845.2_ASM584v2` (NC_000913.3), 4 641 652 bp, GC 50.8 %.
