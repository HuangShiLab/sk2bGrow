//! The 16-enzyme panel used by 2bRAD-M / Syn2b / Fast2bRAD-M.
//!
//! Fifteen are Type IIB: a bipartite site, cut on *both* flanks, excising a
//! short fixed-length fragment. **BslFI is not** — see [`BSLFI_IS_TYPE_IIS`].
//!
//! Transcribed from `Fast2bRAD-M/src/enzymes.rs`, which in turn derives its
//! patterns from the `@site` regexes in `2bRADExtraction.pl`. That is the
//! authoritative definition; this table is verified against it by the tests
//! below and by measured densities on E. coli K-12 MG1655 (see
//! `docs/enzymes.md`).
//!
//! ## The window model
//!
//! A tag is **a fixed-length window of the forward strand** that satisfies one of
//! the enzyme's patterns. A pattern is a set of `(offset, motif)` anchors
//! positioned inside that window:
//!
//! ```text
//!   BcgI, tag_len 32, pattern 0:   [ACGT]{10} CGA [ACGT]{6} TGC [ACGT]{10}
//!                                   ^offset 10 ─┘         ^offset 19 ─┘
//! ```
//!
//! Most enzymes have **two** patterns — the recognition motif as it reads on
//! each strand. Both are matched against the forward strand, and the extracted
//! tag is the forward-strand window either way; nothing is reverse-complemented
//! during extraction. Enzymes whose *whole window pattern* is its own reverse
//! complement (AlfI, BplI, FalI) need only one.
//!
//! This is why tag comparison is strand-canonical (see
//! [`crate::seq::canonical_hash`]): a read sequenced from the opposite strand
//! yields the reverse complement of the reference window.

use serde::{Deserialize, Serialize};

use crate::seq::iupac_matches;

/// Number of enzymes in the panel. `enzyme_idx` in [`crate::anchor_db::Anchor`]
/// is an index into [`PANEL`] and fits in `u8` precisely because this is 16.
pub const N_ENZYMES: usize = 16;

/// One anchor inside a tag window: a motif (IUPAC) at a fixed offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Motif {
    pub offset: u8,
    pub bases: &'static str,
}

impl Motif {
    #[inline]
    fn matches(&self, window: &[u8]) -> bool {
        let o = self.offset as usize;
        let b = self.bases.as_bytes();
        if o + b.len() > window.len() {
            return false;
        }
        b.iter()
            .zip(&window[o..o + b.len()])
            .all(|(&c, &x)| iupac_matches(c, x))
    }
}

/// One way an enzyme's recognition site can sit inside a tag window.
pub type Pattern = &'static [Motif];

/// One Type IIB enzyme. Not serialised: databases persist `enzyme_idx` into
/// [`PANEL`], never the definition itself, so the table can be corrected
/// without invalidating stored anchors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Enzyme {
    /// Index into [`PANEL`]; stored in anchors as `enzyme_idx`.
    pub idx: u8,
    pub name: &'static str,
    /// Excised tag length in bp.
    pub tag_len: u8,
    /// One entry per strand orientation. A single entry means the whole window
    /// pattern is its own reverse complement.
    pub patterns: &'static [Pattern],
    /// Human-readable recognition site, for reports.
    pub display: &'static str,
}

impl Enzyme {
    /// Does `window` (of length [`Self::tag_len`]) match pattern `p`?
    #[inline]
    pub fn pattern_matches(&self, window: &[u8], p: usize) -> bool {
        window.len() == self.tag_len as usize && self.patterns[p].iter().all(|m| m.matches(window))
    }

    /// Index of the first pattern `window` matches, if any.
    #[inline]
    pub fn match_window(&self, window: &[u8]) -> Option<usize> {
        (0..self.patterns.len()).find(|&p| self.pattern_matches(window, p))
    }

    /// True when one pattern suffices because the window pattern is its own
    /// reverse complement.
    pub fn is_self_complementary(&self) -> bool {
        self.patterns.len() == 1
    }
}

macro_rules! motifs {
    ($(($off:expr, $bases:expr)),+ $(,)?) => {
        &[$(Motif { offset: $off, bases: $bases }),+]
    };
}

// Patterns, transcribed from Fast2bRAD-M `src/enzymes.rs`. The second pattern of
// each pair is the reverse complement of the first *as a window pattern*; a test
// below asserts that closure holds for every entry, which is what catches a
// mis-transcribed offset.
static BCGI_P: [Pattern; 2] = [
    motifs![(10, "CGA"), (19, "TGC")],
    motifs![(10, "GCA"), (19, "TCG")],
];
static ALFI_P: [Pattern; 1] = [motifs![(10, "GCA"), (19, "TGC")]];
static ALOI_P: [Pattern; 2] = [
    motifs![(7, "GAAC"), (17, "TCC")],
    motifs![(7, "GGA"), (16, "GTTC")],
];
static BAEI_P: [Pattern; 2] = [
    motifs![(10, "AC"), (16, "GTAYC")],
    motifs![(7, "GRTAC"), (16, "GT")],
];
static BPLI_P: [Pattern; 1] = [motifs![(8, "GAG"), (16, "CTC")]];
static BSAXI_P: [Pattern; 2] = [
    motifs![(9, "AC"), (16, "CTCC")],
    motifs![(7, "GGAG"), (16, "GT")],
];
static BSLFI_P: [Pattern; 2] = [motifs![(6, "GGGAC")], motifs![(14, "GTCCC")]];
static BSP24I_P: [Pattern; 2] = [
    motifs![(8, "GAC"), (17, "TGG")],
    motifs![(7, "CCA"), (16, "GTC")],
];
static CJEI_P: [Pattern; 2] = [
    motifs![(8, "CCA"), (17, "GT")],
    motifs![(9, "AC"), (17, "TGG")],
];
static CJEPI_P: [Pattern; 2] = [
    motifs![(7, "CCA"), (17, "TC")],
    motifs![(8, "GA"), (17, "TGG")],
];
static CSPCI_P: [Pattern; 2] = [
    motifs![(11, "CAA"), (19, "GTGG")],
    motifs![(10, "CCAC"), (19, "TTG")],
];
static FALI_P: [Pattern; 1] = [motifs![(8, "AAG"), (16, "CTT")]];
static HAEIV_P: [Pattern; 2] = [
    motifs![(7, "GAY"), (15, "RTC")],
    motifs![(9, "GAY"), (17, "RTC")],
];
static HIN4I_P: [Pattern; 2] = [
    motifs![(8, "GAY"), (16, "VTC")],
    motifs![(8, "GAB"), (16, "RTC")],
];
static PPII_P: [Pattern; 2] = [
    motifs![(7, "GAAC"), (16, "CTC")],
    motifs![(8, "GAG"), (16, "GTTC")],
];
static PSRI_P: [Pattern; 2] = [
    motifs![(7, "GAAC"), (17, "TAC")],
    motifs![(7, "GTA"), (16, "GTTC")],
];

/// The panel, in the order used by `enzyme_idx`. Do not reorder: `enzyme_idx`
/// is persisted inside anchor databases and count tables.
pub static PANEL: [Enzyme; N_ENZYMES] = [
    Enzyme {
        idx: 0,
        name: "BcgI",
        tag_len: 32,
        patterns: &BCGI_P,
        display: "CGA-N6-TGC",
    },
    Enzyme {
        idx: 1,
        name: "AlfI",
        tag_len: 32,
        patterns: &ALFI_P,
        display: "GCA-N6-TGC",
    },
    Enzyme {
        idx: 2,
        name: "AloI",
        tag_len: 27,
        patterns: &ALOI_P,
        display: "GAAC-N6-TCC",
    },
    Enzyme {
        idx: 3,
        name: "BaeI",
        tag_len: 28,
        patterns: &BAEI_P,
        display: "AC-N4-GTAYC",
    },
    Enzyme {
        idx: 4,
        name: "BplI",
        tag_len: 27,
        patterns: &BPLI_P,
        display: "GAG-N5-CTC",
    },
    Enzyme {
        idx: 5,
        name: "BsaXI",
        tag_len: 27,
        patterns: &BSAXI_P,
        display: "AC-N5-CTCC",
    },
    Enzyme {
        idx: 6,
        name: "BslFI",
        tag_len: 25,
        patterns: &BSLFI_P,
        display: "GGGAC",
    },
    Enzyme {
        idx: 7,
        name: "Bsp24I",
        tag_len: 27,
        patterns: &BSP24I_P,
        display: "GAC-N6-TGG",
    },
    Enzyme {
        idx: 8,
        name: "CjeI",
        tag_len: 28,
        patterns: &CJEI_P,
        display: "CCA-N6-GT",
    },
    Enzyme {
        idx: 9,
        name: "CjePI",
        tag_len: 27,
        patterns: &CJEPI_P,
        display: "CCA-N7-TC",
    },
    Enzyme {
        idx: 10,
        name: "CspCI",
        tag_len: 33,
        patterns: &CSPCI_P,
        display: "CAA-N5-GTGG",
    },
    Enzyme {
        idx: 11,
        name: "FalI",
        tag_len: 27,
        patterns: &FALI_P,
        display: "AAG-N5-CTT",
    },
    Enzyme {
        idx: 12,
        name: "HaeIV",
        tag_len: 27,
        patterns: &HAEIV_P,
        display: "GAY-N5-RTC",
    },
    Enzyme {
        idx: 13,
        name: "Hin4I",
        tag_len: 27,
        patterns: &HIN4I_P,
        display: "GAY-N5-VTC",
    },
    Enzyme {
        idx: 14,
        name: "PpiI",
        tag_len: 27,
        patterns: &PPII_P,
        display: "GAAC-N5-CTC",
    },
    Enzyme {
        idx: 15,
        name: "PsrI",
        tag_len: 27,
        patterns: &PSRI_P,
        display: "GAAC-N6-TAC",
    },
];

/// Look an enzyme up by name (case-insensitive).
pub fn by_name(name: &str) -> Option<&'static Enzyme> {
    PANEL.iter().find(|e| e.name.eq_ignore_ascii_case(name))
}

/// Look an enzyme up by panel index.
pub fn by_idx(idx: u8) -> Option<&'static Enzyme> {
    PANEL.get(idx as usize)
}

/// Parse a `--enzymes` selection: `all`, or a comma-separated list of names.
///
/// The result is deduplicated and sorted by panel index, so enzyme order in a
/// database never depends on how the user typed the flag.
pub fn parse_selection(spec: &str) -> Result<Vec<&'static Enzyme>, EnzymeError> {
    let spec = spec.trim();
    if spec.eq_ignore_ascii_case("all") {
        return Ok(PANEL.iter().collect());
    }
    let mut picked: Vec<&'static Enzyme> = Vec::new();
    for token in spec.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        let e = by_name(token).ok_or_else(|| EnzymeError::Unknown(token.to_string()))?;
        if !picked.iter().any(|p| p.idx == e.idx) {
            picked.push(e);
        }
    }
    if picked.is_empty() {
        return Err(EnzymeError::EmptySelection);
    }
    picked.sort_by_key(|e| e.idx);
    Ok(picked)
}

/// Bitset over the panel; records which enzymes an anchor database holds.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnzymeSet(pub u16);

impl EnzymeSet {
    pub fn from_slice(enzymes: &[&'static Enzyme]) -> Self {
        EnzymeSet(enzymes.iter().fold(0u16, |acc, e| acc | (1 << e.idx)))
    }
    #[inline]
    pub fn contains(&self, idx: u8) -> bool {
        self.0 & (1 << idx) != 0
    }
    pub fn len(&self) -> usize {
        self.0.count_ones() as usize
    }
    pub fn is_empty(&self) -> bool {
        self.0 == 0
    }
    pub fn iter(&self) -> impl Iterator<Item = &'static Enzyme> + '_ {
        PANEL.iter().filter(move |e| self.contains(e.idx))
    }
}

/// BslFI is Type IIS, not Type IIB, unlike the other fifteen.
///
/// REBASE gives `GGGAC(10-11/14-15)`, "Type II restriction enzyme, subtype: S"
/// (prototype FinI): a contiguous 5 bp motif, cut downstream only and at
/// *variable* positions, excising nothing. The panel's `N6 GGGAC N14` window
/// therefore has a real (4 nt staggered) cut at its right edge and arbitrary
/// padding at its left, where the fifteen Type IIB windows have genuine cuts at
/// both edges.
///
/// This is *not* about the two patterns: `GTCCC` is the reverse complement of
/// `GGGAC`, so BslFI's second pattern is the same window seen from the other
/// strand, exactly as for every other enzyme in the panel.
///
/// In silico the window is still deterministic, reverse-complement closed and
/// applied identically to reference and reads — a perfectly good marker stratum.
/// At the bench it is not executable as 2bRAD: a real digest yields kilobase
/// fragments (mean 1.3–1.7 kb across E. coli, B. subtilis and P. putida; only
/// ~3 % below 40 bp), so there is no short band to size-select. Consequence for
/// this project: BslFI's route-B stratum can never be validated against real
/// 2bRAD data.
pub const BSLFI_IS_TYPE_IIS: bool = true;

/// One enzyme's tags being a subset of another's.
// No `Eq`: `measured_fraction` is f64.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Containment {
    pub subset: &'static str,
    pub superset: &'static str,
    /// `true` when *every* tag of `subset` is also a tag of `superset`.
    pub total: bool,
    /// Measured share of `subset`'s tags also claimed by `superset`, on
    /// E. coli K-12 / B. subtilis 168 / P. putida KT2440.
    pub measured_fraction: [f64; 3],
}

/// Containment relations within the panel, proved from the patterns and measured
/// on three genomes spanning 43.5–61.5 % GC.
///
/// These matter because `python/sk2bgrow/fusion.py` treats the enzymes as
/// independent measurement strata. They are not:
///
/// * **Bsp24I ⊂ CjePI, totally.** The two enzymes' patterns are written in
///   opposite strand orientations, so the comparison is Bsp24I pattern 0 against
///   CjePI pattern *1*. Aligned over the 27 bp window they differ at exactly one
///   position — Bsp24I fixes a `C` at offset 10 that CjePI leaves free — so
///   CjePI is the less specific enzyme and its site set contains Bsp24I's.
///   Measured 1 636/1 636, 891/891, 2 910/2 910.
/// * **Bsp24I pattern 0 ⊂ CjeI pattern 1**, over their shared 27-base prefix
///   (CjeI's tag is one base longer, so an equal-length test misses it). About
///   half of Bsp24I's tags — exactly its pattern-0 share.
///
/// So the panel offers at most ~15 independent strata, not 16.
pub static CONTAINMENTS: &[Containment] = &[
    Containment {
        subset: "Bsp24I",
        superset: "CjePI",
        total: true,
        measured_fraction: [1.0, 1.0, 1.0],
    },
    Containment {
        subset: "Bsp24I",
        superset: "CjeI",
        total: false,
        measured_fraction: [0.484, 0.474, 0.509],
    },
];

/// Enzymes carrying no information independent of another panel member. A caller
/// building a stratum set should drop these or declare them sub-strata.
pub fn redundant_enzymes() -> Vec<&'static str> {
    CONTAINMENTS
        .iter()
        .filter(|c| c.total)
        .map(|c| c.subset)
        .collect()
}

/// True when `a`'s tags are a subset of `b`'s.
pub fn is_contained_in(a: &str, b: &str) -> bool {
    CONTAINMENTS
        .iter()
        .any(|c| c.total && c.subset.eq_ignore_ascii_case(a) && c.superset.eq_ignore_ascii_case(b))
}

#[derive(Debug, thiserror::Error)]
pub enum EnzymeError {
    #[error("unknown enzyme '{0}'; the 16-enzyme panel is: {}", PANEL.iter().map(|e| e.name).collect::<Vec<_>>().join(", "))]
    Unknown(String),
    #[error("empty enzyme selection")]
    EmptySelection,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{iupac_mask, revcomp};

    /// Expand a pattern into a per-position set of allowed bases.
    fn constraints(e: &Enzyme, p: usize) -> Vec<u8> {
        let mut c = vec![0b1111u8; e.tag_len as usize];
        for m in e.patterns[p] {
            for (i, &b) in m.bases.as_bytes().iter().enumerate() {
                c[m.offset as usize + i] = iupac_mask(b);
            }
        }
        c
    }

    #[test]
    fn every_pattern_set_is_reverse_complement_closed() {
        // The reverse complement of a tag window is itself a valid tag of the
        // same enzyme, read from the other strand. So reverse-complementing a
        // pattern must yield another pattern of the same enzyme. This is the
        // check that catches a mis-transcribed offset or motif.
        for e in PANEL.iter() {
            let all: Vec<Vec<u8>> = (0..e.patterns.len()).map(|p| constraints(e, p)).collect();
            for c in &all {
                let rc: Vec<u8> = c
                    .iter()
                    .rev()
                    .map(|&mask| {
                        // complement a 4-bit ACGT mask: A<->T (bit0<->bit3), C<->G (bit1<->bit2)
                        ((mask & 0b0001) << 3)
                            | ((mask & 0b0010) << 1)
                            | ((mask & 0b0100) >> 1)
                            | ((mask & 0b1000) >> 3)
                    })
                    .collect();
                assert!(
                    all.contains(&rc),
                    "{}: reverse complement of a pattern is not in the panel",
                    e.name
                );
            }
        }
    }

    #[test]
    fn motif_offsets_fit_inside_the_tag() {
        for e in PANEL.iter() {
            for p in e.patterns {
                for m in *p {
                    assert!(
                        m.offset as usize + m.bases.len() <= e.tag_len as usize,
                        "{} motif {} at {} overruns the {} bp tag",
                        e.name,
                        m.bases,
                        m.offset,
                        e.tag_len
                    );
                }
            }
        }
    }

    #[test]
    fn panel_indices_match_position() {
        for (i, e) in PANEL.iter().enumerate() {
            assert_eq!(e.idx as usize, i, "{} has a stale idx", e.name);
        }
        assert_eq!(PANEL.len(), N_ENZYMES);
    }

    #[test]
    fn self_complementary_enzymes_are_exactly_these_three() {
        // AlfI, BplI and FalI have symmetric flanks around a palindromic core,
        // so one pattern covers both strands. HaeIV's core is palindromic too,
        // but its flanks are 7/9, so its two readings occupy different windows
        // and it needs two patterns.
        let one: Vec<&str> = PANEL
            .iter()
            .filter(|e| e.is_self_complementary())
            .map(|e| e.name)
            .collect();
        assert_eq!(one, vec!["AlfI", "BplI", "FalI"]);
        assert_eq!(by_name("HaeIV").unwrap().patterns.len(), 2);
    }

    #[test]
    fn matches_a_planted_bcgi_window() {
        let bcgi = by_name("BcgI").unwrap();
        // 10 filler, CGA, 6 filler, TGC, 10 filler = 32
        let w = b"AAAAAAAAAACGAACGTACTGCTTTTTTTTTT";
        assert_eq!(w.len(), 32);
        assert_eq!(bcgi.match_window(w), Some(0));
        // The reverse complement must match the other pattern.
        assert_eq!(bcgi.match_window(&revcomp(w)), Some(1));
        // A window of the wrong length never matches.
        assert_eq!(bcgi.match_window(&w[..31]), None);
    }

    #[test]
    fn degenerate_motifs_are_honoured() {
        let hae = by_name("HaeIV").unwrap();
        // pattern 0: N7 GAY N5 RTC N9
        let mut w = vec![b'A'; 27];
        w[7..10].copy_from_slice(b"GAT"); // Y = T
        w[15..18].copy_from_slice(b"ATC"); // R = A
        assert_eq!(hae.match_window(&w), Some(0));
        w[9] = b'A'; // Y cannot be A
        assert_eq!(hae.match_window(&w), None);
    }

    #[test]
    fn bsp24i_tags_are_all_cjepi_tags() {
        // A structural property of the panel, not an accident of one genome:
        // Bsp24I's patterns are strict refinements of CjePI's.
        let bsp = by_name("Bsp24I").unwrap();
        let cje = by_name("CjePI").unwrap();
        assert_eq!(bsp.tag_len, cje.tag_len);
        for p in 0..bsp.patterns.len() {
            let c_bsp = constraints(bsp, p);
            let ok = (0..cje.patterns.len()).any(|q| {
                let c_cje = constraints(cje, q);
                c_bsp.iter().zip(&c_cje).all(|(&a, &b)| a & b == a)
            });
            assert!(
                ok,
                "Bsp24I pattern {p} is not contained in any CjePI pattern"
            );
        }
        assert!(is_contained_in("Bsp24I", "CjePI"));
        assert!(!is_contained_in("CjePI", "Bsp24I"));
        assert_eq!(redundant_enzymes(), vec!["Bsp24I"]);
        // The CjeI relation is partial and must not be reported as total.
        assert!(!is_contained_in("Bsp24I", "CjeI"));
        assert!(CONTAINMENTS
            .iter()
            .any(|c| c.superset == "CjeI" && !c.total));
    }

    #[test]
    fn selection_is_order_independent() {
        let a = parse_selection("PsrI,BcgI").unwrap();
        let b = parse_selection("BcgI,PsrI,BcgI").unwrap();
        assert_eq!(a, b);
        assert_eq!(a[0].name, "BcgI");
        assert_eq!(parse_selection("all").unwrap().len(), N_ENZYMES);
        assert!(parse_selection("EcoRI").is_err());
    }

    #[test]
    fn enzyme_set_roundtrips() {
        let sel = parse_selection("BcgI,AlfI,PsrI").unwrap();
        let set = EnzymeSet::from_slice(&sel);
        assert_eq!(set.len(), 3);
        assert!(set.contains(0) && !set.contains(2));
        assert_eq!(set.iter().count(), 3);
    }
}
