//! The 16 Type IIB enzyme panel used by 2bRAD-M / Syn2b.
//!
//! Reuse target: `bsyn::enzyme`. The table below is the vendored copy so this
//! workspace builds standalone; recognition patterns, tag lengths and per-genome
//! densities are the ones tabulated in the design report §4.1.
//!
//! Anatomy of a tag. A Type IIB enzyme excises a fragment spanning its own
//! recognition site plus a fixed flank on each side:
//!
//! ```text
//!        up_flank            recognition (pattern)            down_flank
//!   |<--------------->|<--------------------------->|<--------------------->|
//!   ^ tag_start                                                   tag_end ^
//!                     ^ site_start (anchor coordinate)
//! ```
//!
//! `tag_len == up_flank + pattern.len() + down_flank` holds for every entry and
//! is asserted in the tests. The anchor coordinate reported downstream is
//! `site_start` — the recognition site, not the tag start — because that is the
//! quantity invariant to flank conventions and therefore comparable against
//! `bsyn`.

use serde::{Deserialize, Serialize};

use crate::seq::{is_palindromic, iupac_matches};

/// Number of enzymes in the panel. `enzyme_idx` in [`crate::anchor_db::Anchor`]
/// is an index into [`PANEL`] and fits in `u8` precisely because this is 16.
pub const N_ENZYMES: usize = 16;

/// One Type IIB enzyme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Enzyme {
    /// Index into [`PANEL`]; stored in anchors as `enzyme_idx`.
    pub idx: u8,
    pub name: &'static str,
    /// Recognition pattern in IUPAC, N-runs expanded (e.g. `CGANNNNNNTGC`).
    pub pattern: &'static str,
    /// Human-readable form used in reports (e.g. `CGA-N6-TGC`).
    pub display: &'static str,
    /// Excised tag length in bp, as tabulated in the report.
    pub tag_len: u8,
    /// Bases 5' of the recognition site that belong to the tag.
    pub up_flank: u8,
    /// Bases 3' of the recognition site that belong to the tag.
    pub down_flank: u8,
    /// REBASE cut-site notation, retained for provenance when reconciling this
    /// table against `bsyn::enzyme`.
    pub rebase: &'static str,
}

impl Enzyme {
    /// Recognition pattern as bytes.
    #[inline]
    pub fn pattern_bytes(&self) -> &'static [u8] {
        self.pattern.as_bytes()
    }

    /// Length of the recognition site in bp.
    #[inline]
    pub fn pattern_len(&self) -> usize {
        self.pattern.len()
    }

    /// A palindromic pattern hits the same physical locus from both strands;
    /// [`crate::digest`] deduplicates those so a locus is not double-counted.
    pub fn is_palindromic(&self) -> bool {
        is_palindromic(self.pattern_bytes())
    }

    /// Does the recognition site occur at `seq[pos..]` on the forward strand?
    pub fn matches_at(&self, seq: &[u8], pos: usize) -> bool {
        let pat = self.pattern_bytes();
        if pos + pat.len() > seq.len() {
            return false;
        }
        pat.iter()
            .zip(&seq[pos..pos + pat.len()])
            .all(|(&code, &base)| iupac_matches(code, base))
    }

    /// Half-open tag span `[start, end)` for a site starting at `site_start` on
    /// strand `fwd`. Returns `None` when the tag would run off the contig — a
    /// site within one tag length of a contig end yields no anchor.
    pub fn tag_span(
        &self,
        site_start: usize,
        contig_len: usize,
        fwd: bool,
    ) -> Option<(usize, usize)> {
        let (up, down) = if fwd {
            (self.up_flank as usize, self.down_flank as usize)
        } else {
            // On the reverse strand the flanks swap: what is "upstream of the
            // site" in enzyme coordinates lies downstream in contig coordinates.
            (self.down_flank as usize, self.up_flank as usize)
        };
        let start = site_start.checked_sub(up)?;
        let end = site_start + self.pattern_len() + down;
        if end > contig_len {
            return None;
        }
        Some((start, end))
    }
}

/// The panel, in the order used by `enzyme_idx`. Do not reorder: `enzyme_idx`
/// is persisted inside anchor databases and count tables.
///
/// `up_flank`/`down_flank` follow REBASE cut coordinates where those are
/// unambiguous, and split the remaining flank evenly otherwise. Because every
/// downstream statistic keys on `site_start`, a flank convention that differs
/// from `bsyn` by a base or two shifts tag sequences but not anchor positions or
/// densities. See `docs/enzymes.md` for the reconciliation checklist.
pub static PANEL: [Enzyme; N_ENZYMES] = [
    Enzyme {
        idx: 0,
        name: "BcgI",
        pattern: "CGANNNNNNTGC",
        display: "CGA-N6-TGC",
        tag_len: 32,
        up_flank: 10,
        down_flank: 10,
        rebase: "(10/12)CGA N6 TGC(12/10)",
    },
    Enzyme {
        idx: 1,
        name: "AlfI",
        pattern: "GCANNNNNNTGC",
        display: "GCA-N6-TGC",
        tag_len: 32,
        up_flank: 10,
        down_flank: 10,
        rebase: "(10/12)GCA N6 TGC(12/10)",
    },
    Enzyme {
        idx: 2,
        name: "AloI",
        pattern: "GAACNNNNNNTCC",
        display: "GAAC-N6-TCC",
        tag_len: 27,
        up_flank: 7,
        down_flank: 7,
        rebase: "(7/12)GAAC N6 TCC(12/7)",
    },
    Enzyme {
        idx: 3,
        name: "BaeI",
        pattern: "ACNNNNGTAYC",
        display: "AC-N4-GTAYC",
        tag_len: 28,
        up_flank: 10,
        down_flank: 7,
        rebase: "(10/15)AC N4 GTAYC(12/7)",
    },
    Enzyme {
        idx: 4,
        name: "BplI",
        pattern: "GAGNNNNNCTC",
        display: "GAG-N5-CTC",
        tag_len: 27,
        up_flank: 8,
        down_flank: 8,
        rebase: "(8/13)GAG N5 CTC(13/8)",
    },
    Enzyme {
        idx: 5,
        name: "BsaXI",
        pattern: "ACNNNNNCTCC",
        display: "AC-N5-CTCC",
        tag_len: 27,
        up_flank: 9,
        down_flank: 7,
        rebase: "(9/12)AC N5 CTCC(10/7)",
    },
    Enzyme {
        idx: 6,
        name: "BslFI",
        pattern: "GGGAC",
        display: "GGGAC",
        tag_len: 25,
        up_flank: 0,
        down_flank: 20,
        rebase: "GGGAC(10/14)",
    },
    Enzyme {
        idx: 7,
        name: "Bsp24I",
        pattern: "GACNNNNNNTGG",
        display: "GAC-N6-TGG",
        tag_len: 27,
        up_flank: 8,
        down_flank: 7,
        rebase: "(8/13)GAC N6 TGG(12/7)",
    },
    Enzyme {
        idx: 8,
        name: "CjeI",
        pattern: "CCANNNNNNGT",
        display: "CCA-N6-GT",
        tag_len: 28,
        up_flank: 9,
        down_flank: 8,
        rebase: "(8/14)CCA N6 GT(9/8)",
    },
    Enzyme {
        idx: 9,
        name: "CjePI",
        pattern: "CCANNNNNNNTC",
        display: "CCA-N7-TC",
        tag_len: 27,
        up_flank: 8,
        down_flank: 7,
        rebase: "(8/13)CCA N7 TC(13/8)",
    },
    Enzyme {
        idx: 10,
        name: "CspCI",
        pattern: "CAANNNNNGTGG",
        display: "CAA-N5-GTGG",
        tag_len: 33,
        up_flank: 11,
        down_flank: 10,
        rebase: "(11/13)CAA N5 GTGG(12/10)",
    },
    Enzyme {
        idx: 11,
        name: "FalI",
        pattern: "AAGNNNNNCTT",
        display: "AAG-N5-CTT",
        tag_len: 27,
        up_flank: 8,
        down_flank: 8,
        rebase: "(8/13)AAG N5 CTT(13/8)",
    },
    // HaeIV's recognition pattern is its own reverse complement, so the excised
    // duplex must be strand-symmetric: the flanks are split evenly rather than
    // following the asymmetric REBASE overhang coordinates.
    Enzyme {
        idx: 12,
        name: "HaeIV",
        pattern: "GAYNNNNNRTC",
        display: "GAY-N5-RTC",
        tag_len: 27,
        up_flank: 8,
        down_flank: 8,
        rebase: "(7/13)GAY N5 RTC(14/9)",
    },
    Enzyme {
        idx: 13,
        name: "Hin4I",
        pattern: "GAYNNNNNVTC",
        display: "GAY-N5-VTC",
        tag_len: 27,
        up_flank: 8,
        down_flank: 8,
        rebase: "(8/13)GAY N5 VTC(13/8)",
    },
    Enzyme {
        idx: 14,
        name: "PpiI",
        pattern: "GAACNNNNNCTC",
        display: "GAAC-N5-CTC",
        tag_len: 28,
        up_flank: 8,
        down_flank: 8,
        rebase: "(7/12)GAAC N5 CTC(13/8)",
    },
    Enzyme {
        idx: 15,
        name: "PsrI",
        pattern: "GAACNNNNNNTAC",
        display: "GAAC-N6-TAC",
        tag_len: 27,
        up_flank: 7,
        down_flank: 7,
        rebase: "(7/12)GAAC N6 TAC(12/7)",
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
/// The returned vector is deduplicated and sorted by panel index, so the enzyme
/// order in an anchor database never depends on how the user typed the flag.
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

/// Bitset over the panel; used to record which enzymes an anchor database holds.
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

    #[test]
    fn tag_length_is_consistent_with_flanks() {
        for e in PANEL.iter() {
            assert_eq!(
                e.up_flank as usize + e.pattern_len() + e.down_flank as usize,
                e.tag_len as usize,
                "{} flank/tag_len mismatch",
                e.name
            );
        }
    }

    #[test]
    fn palindromic_enzymes_have_symmetric_flanks() {
        // A palindromic recognition site is cut symmetrically on the duplex, so
        // the forward and reverse copies of a locus must excise the same span.
        // Without this the two strands would disagree on the tag sequence.
        let pal: Vec<&str> = PANEL
            .iter()
            .filter(|e| e.is_palindromic())
            .map(|e| e.name)
            .collect();
        assert_eq!(
            pal,
            vec!["AlfI", "BplI", "FalI", "HaeIV"],
            "palindrome set drifted"
        );
        for e in PANEL.iter().filter(|e| e.is_palindromic()) {
            assert_eq!(
                e.up_flank, e.down_flank,
                "{} palindrome with asymmetric flanks",
                e.name
            );
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
    fn display_expands_to_pattern() {
        // `CGA-N6-TGC` must expand to exactly the stored pattern.
        for e in PANEL.iter() {
            let expanded: String = e
                .display
                .split('-')
                .map(|part| {
                    if let Some(n) = part.strip_prefix('N') {
                        if let Ok(k) = n.parse::<usize>() {
                            return "N".repeat(k);
                        }
                    }
                    part.to_string()
                })
                .collect();
            assert_eq!(expanded, e.pattern, "{} display/pattern mismatch", e.name);
        }
    }

    #[test]
    fn recognises_sites_with_degenerate_bases() {
        let haeiv = by_name("HaeIV").unwrap();
        // GAY-N5-RTC with Y=T, R=A
        assert!(haeiv.matches_at(b"GATCCCCCATC", 0));
        // Y=C, R=G
        assert!(haeiv.matches_at(b"GACGGGGGGTC", 0));
        // Y=A is not allowed
        assert!(!haeiv.matches_at(b"GAACCCCCATC", 0));
    }

    #[test]
    fn tag_span_rejects_contig_edges() {
        let bcgi = by_name("BcgI").unwrap();
        // 10 bp of upstream flank required.
        assert!(bcgi.tag_span(5, 1_000, true).is_none());
        assert_eq!(bcgi.tag_span(10, 1_000, true), Some((0, 32)));
        // and 10 bp downstream of the 12 bp site.
        assert!(bcgi.tag_span(980, 1_000, true).is_none());
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
        assert!(set.contains(0));
        assert!(!set.contains(2));
        assert_eq!(set.iter().count(), 3);
    }
}
