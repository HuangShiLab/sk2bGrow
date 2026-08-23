//! In silico Type IIB digestion: reference sequence -> tags.
//!
//! Follows `Fast2bRAD-M/src/enzymes.rs::find_all_tags`. A tag is a fixed-length
//! **window of the forward strand** that satisfies one of the enzyme's patterns.
//! Both strand orientations are covered by the enzyme having two patterns (see
//! [`crate::enzyme`]), so nothing is reverse-complemented here: the tag is the
//! window, verbatim.
//!
//! Every occurrence is reported, including overlapping ones — the scan advances
//! one base at a time rather than skipping past a match. That mirrors
//! `2bRADExtraction.pl`, which rewinds its regex cursor to `match_start + 1`.

use serde::{Deserialize, Serialize};

use crate::enzyme::Enzyme;
use crate::seq::canonical_hash;

/// One digested tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub enzyme_idx: u8,
    pub contig_id: u16,
    /// 0-based start of the tag window on the forward strand.
    pub position: u64,
    /// Which of the enzyme's patterns matched: 0 is the motif as written,
    /// 1 its reverse-complement reading. Recorded as `strand` downstream.
    pub pattern: u8,
    /// The tag: `contig[position .. position + tag_len]`, forward strand.
    pub tag: Vec<u8>,
}

impl Site {
    /// Strand-canonical hash — a read off the opposite strand yields the reverse
    /// complement of this window, and must hash to the same value.
    pub fn tag_hash(&self) -> u64 {
        canonical_hash(&self.tag)
    }
    #[inline]
    pub fn end(&self) -> u64 {
        self.position + self.tag.len() as u64
    }
}

/// Cross-enzyme union merge radius, in bp: two tags whose windows start within
/// this distance are treated as one effective locus for spacing statistics
/// (report §4.1). 33 bp is the longest tag in the panel (CspCI).
pub const DEFAULT_MERGE_WINDOW: u64 = 33;

#[derive(Debug, Clone)]
pub struct DigestConfig {
    /// Skip windows containing a non-ACGT base. Reference `N` runs would
    /// otherwise produce tags no read can match.
    pub reject_ambiguous_tags: bool,
    /// Skip contigs shorter than this (bp).
    pub min_contig_len: usize,
}

impl Default for DigestConfig {
    fn default() -> Self {
        Self {
            reject_ambiguous_tags: true,
            min_contig_len: 500,
        }
    }
}

#[inline]
fn is_pure_acgt(window: &[u8]) -> bool {
    window
        .iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

/// Digest one contig with one enzyme, appending to `out`.
pub fn digest_contig_with(
    seq: &[u8],
    contig_id: u16,
    enzyme: &Enzyme,
    cfg: &DigestConfig,
    out: &mut Vec<Site>,
) {
    let len = enzyme.tag_len as usize;
    if seq.len() < cfg.min_contig_len || seq.len() < len {
        return;
    }
    for start in 0..=(seq.len() - len) {
        let window = &seq[start..start + len];
        if cfg.reject_ambiguous_tags && !is_pure_acgt(window) {
            continue;
        }
        // One entry per window: a window matching several of an enzyme's own
        // patterns is still one tag, credited to the first that matched.
        if let Some(p) = enzyme.match_window(window) {
            out.push(Site {
                enzyme_idx: enzyme.idx,
                contig_id,
                position: start as u64,
                pattern: p as u8,
                tag: window.to_vec(),
            });
        }
    }
}

/// Digest one contig with a whole enzyme selection, sorted by `(position, enzyme)`.
pub fn digest_contig(
    seq: &[u8],
    contig_id: u16,
    enzymes: &[&'static Enzyme],
    cfg: &DigestConfig,
) -> Vec<Site> {
    let mut sites = Vec::new();
    for e in enzymes {
        digest_contig_with(seq, contig_id, e, cfg, &mut sites);
    }
    sites.sort_by_key(|s| (s.position, s.enzyme_idx, s.pattern));
    sites
}

/// Per-enzyme counts and union spacing statistics.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DensityReport {
    pub genome_len: u64,
    /// Tag windows per enzyme.
    pub per_enzyme: Vec<(String, u64)>,
    /// Distinct tag windows across all enzymes (a window shared by two enzymes
    /// counts once).
    pub union_windows: u64,
    /// Independent loci after merging within [`DEFAULT_MERGE_WINDOW`].
    pub union_sites: u64,
    pub mean_spacing: f64,
    /// Largest gap between consecutive union loci — the "blind spot" statistic.
    pub max_gap: u64,
    /// Union loci per 25 kb, comparable to Pilea's ~100 k-mers per window.
    pub per_25kb: f64,
}

impl DensityReport {
    /// Tags per Mb for one enzyme, in the units of report table §4.1.
    pub fn density_per_mb(&self, enzyme_name: &str) -> Option<f64> {
        if self.genome_len == 0 {
            return None;
        }
        self.per_enzyme
            .iter()
            .find(|(n, _)| n == enzyme_name)
            .map(|(_, c)| *c as f64 / (self.genome_len as f64 / 1e6))
    }
}

/// Collapse tags into independent loci: sorted by position, merging any two
/// whose windows start within `merge_window` bp *on the same contig*.
pub fn merge_union(sites: &[Site], merge_window: u64) -> Vec<(u16, u64)> {
    let mut keys: Vec<(u16, u64)> = sites.iter().map(|s| (s.contig_id, s.position)).collect();
    keys.sort_unstable();
    keys.dedup();
    let mut merged: Vec<(u16, u64)> = Vec::with_capacity(keys.len());
    for (contig, pos) in keys {
        match merged.last() {
            Some(&(c, p)) if c == contig && pos.saturating_sub(p) <= merge_window => {}
            _ => merged.push((contig, pos)),
        }
    }
    merged
}

/// Build the density audit for a digested genome.
///
/// Spacing is computed within each contig; a gap across a contig boundary is not
/// a genomic gap and would inflate `max_gap` on fragmented MAGs.
pub fn density_report(
    sites: &[Site],
    contig_lens: &[u64],
    enzymes: &[&'static Enzyme],
    merge_window: u64,
) -> DensityReport {
    let genome_len: u64 = contig_lens.iter().sum();
    let mut per_enzyme: Vec<(String, u64)> = enzymes
        .iter()
        .map(|e| {
            (
                e.name.to_string(),
                sites.iter().filter(|s| s.enzyme_idx == e.idx).count() as u64,
            )
        })
        .collect();
    per_enzyme.sort_by(|a, b| a.0.cmp(&b.0));

    let mut windows: Vec<(u16, u64)> = sites.iter().map(|s| (s.contig_id, s.position)).collect();
    windows.sort_unstable();
    windows.dedup();

    let merged = merge_union(sites, merge_window);
    let mut max_gap = 0u64;
    let mut gap_sum = 0u64;
    let mut gap_n = 0u64;
    for w in merged.windows(2) {
        let ((c0, p0), (c1, p1)) = (w[0], w[1]);
        if c0 != c1 {
            continue;
        }
        let g = p1 - p0;
        max_gap = max_gap.max(g);
        gap_sum += g;
        gap_n += 1;
    }
    DensityReport {
        genome_len,
        per_enzyme,
        union_windows: windows.len() as u64,
        union_sites: merged.len() as u64,
        mean_spacing: if gap_n > 0 {
            gap_sum as f64 / gap_n as f64
        } else {
            f64::NAN
        },
        max_gap,
        per_25kb: if genome_len > 0 {
            merged.len() as f64 / (genome_len as f64 / 25_000.0)
        } else {
            0.0
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzyme::by_name;
    use crate::seq::revcomp;

    fn cfg() -> DigestConfig {
        DigestConfig {
            reject_ambiguous_tags: true,
            min_contig_len: 0,
        }
    }

    #[test]
    fn finds_a_planted_bcgi_tag() {
        let bcgi = by_name("BcgI").unwrap();
        let mut seq = vec![b'A'; 20];
        seq.extend_from_slice(b"AAAAAAAAAACGAACGTACTGCTTTTTTTTTT"); // 32 bp tag at offset 20
        seq.extend_from_slice(&[b'A'; 20]);
        let sites = digest_contig(&seq, 0, &[bcgi], &cfg());
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].position, 20);
        assert_eq!(sites[0].pattern, 0);
        assert_eq!(sites[0].tag.len(), 32);
        // The tag is the forward-strand window, verbatim.
        assert_eq!(sites[0].tag, seq[20..52].to_vec());
    }

    #[test]
    fn the_reverse_reading_is_found_by_the_second_pattern() {
        let bcgi = by_name("BcgI").unwrap();
        let fwd = b"AAAAAAAAAACGAACGTACTGCTTTTTTTTTT";
        let mut seq = vec![b'A'; 20];
        seq.extend_from_slice(&revcomp(fwd));
        seq.extend_from_slice(&[b'A'; 20]);
        let sites = digest_contig(&seq, 0, &[bcgi], &cfg());
        assert_eq!(sites.len(), 1);
        assert_eq!(
            sites[0].pattern, 1,
            "reverse reading should match pattern 1"
        );
        // Still stored forward-strand; canonical hashing is what links the two.
        assert_eq!(sites[0].tag, revcomp(fwd));
        assert_eq!(sites[0].tag_hash(), canonical_hash(fwd));
    }

    #[test]
    fn self_complementary_enzyme_yields_one_tag_per_window() {
        let bpli = by_name("BplI").unwrap();
        assert!(bpli.is_self_complementary());
        let mut seq = vec![b'A'; 20];
        seq.extend_from_slice(b"AAAAAAAAGAGACGTACTCAAAAAAAA"); // 27 bp
        seq.extend_from_slice(&[b'A'; 20]);
        let sites = digest_contig(&seq, 0, &[bpli], &cfg());
        assert_eq!(sites.len(), 1, "one window, one tag");
    }

    #[test]
    fn haeiv_yields_two_windows_per_palindromic_locus() {
        // HaeIV's core GAY-N5-RTC is its own reverse complement but its flanks
        // are 7/9, so the two readings occupy windows offset by 2 bp. Both are
        // real tags. (The design report deduplicates these to a locus count,
        // which is why its HaeIV density is exactly half the window count.)
        let hae = by_name("HaeIV").unwrap();
        let mut seq = vec![b'C'; 40];
        // Put GAT-N5-ATC (a valid GAY-N5-RTC) at position 40.
        seq.extend_from_slice(b"GATCCCCCATC");
        seq.extend_from_slice(&[b'C'; 40]);
        let sites = digest_contig(&seq, 0, &[hae], &cfg());
        assert_eq!(sites.len(), 2, "expected both readings of the locus");
        assert_eq!(sites[1].position - sites[0].position, 2);
        assert_eq!(sites[0].pattern, 1);
        assert_eq!(sites[1].pattern, 0);
    }

    #[test]
    fn windows_with_n_are_rejected() {
        let bcgi = by_name("BcgI").unwrap();
        let mut seq = vec![b'A'; 20];
        seq.extend_from_slice(b"AAAAANNNNNCGAACGTACTGCTTTTTTTTTT");
        seq.extend_from_slice(&[b'A'; 20]);
        assert!(digest_contig(&seq, 0, &[bcgi], &cfg()).is_empty());
        let lax = DigestConfig {
            reject_ambiguous_tags: false,
            min_contig_len: 0,
        };
        assert_eq!(digest_contig(&seq, 0, &[bcgi], &lax).len(), 1);
    }

    #[test]
    fn overlapping_occurrences_are_all_reported() {
        // BslFI: N6 GGGAC N14. Two GGGAC 1 bp apart would give two windows;
        // use a run that produces several overlapping hits.
        let bslfi = by_name("BslFI").unwrap();
        let mut seq = vec![b'A'; 10];
        seq.extend_from_slice(b"GGGACGGGAC");
        seq.extend_from_slice(&[b'A'; 30]);
        let sites = digest_contig(&seq, 0, &[bslfi], &cfg());
        assert_eq!(sites.len(), 2, "the scan must not skip past a match");
        assert_eq!(sites[1].position - sites[0].position, 5);
    }

    #[test]
    fn union_merges_nearby_windows_and_skips_contig_joins() {
        let mk = |e: u8, c: u16, p: u64| Site {
            enzyme_idx: e,
            contig_id: c,
            position: p,
            pattern: 0,
            tag: vec![b'A'; 27],
        };
        let sites = vec![mk(0, 0, 100), mk(1, 0, 110), mk(2, 0, 600), mk(3, 1, 105)];
        assert_eq!(
            merge_union(&sites, DEFAULT_MERGE_WINDOW),
            vec![(0, 100), (0, 600), (1, 105)]
        );
        let e = vec![by_name("BcgI").unwrap()];
        let rep = density_report(&sites, &[1_000, 1_000], &e, DEFAULT_MERGE_WINDOW);
        assert_eq!(rep.union_sites, 3);
        assert_eq!(rep.max_gap, 500, "a contig boundary leaked into max_gap");
        assert_eq!(rep.genome_len, 2_000);
    }
}
