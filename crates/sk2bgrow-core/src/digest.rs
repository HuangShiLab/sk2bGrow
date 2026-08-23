//! In silico Type IIB digestion: reference sequence -> anchor sites.
//!
//! Reuse target: `bsyn::digest`. Semantics follow the design report §4.1 —
//! double-stranded search with IUPAC degeneracy, palindromic enzymes
//! deduplicated at the locus level, and the cross-enzyme union merged within
//! [`DEFAULT_MERGE_WINDOW`] bp so overlapping tags count as one effective site.

use serde::{Deserialize, Serialize};

use crate::enzyme::Enzyme;
use crate::seq::{hash_tag, revcomp};

/// Strand of an anchor. Stored in [`crate::anchor_db::Anchor::strand`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[repr(u8)]
pub enum Strand {
    Fwd = 0,
    Rev = 1,
}

impl Strand {
    #[inline]
    pub fn as_u8(self) -> u8 {
        self as u8
    }
    #[inline]
    pub fn from_u8(v: u8) -> Self {
        if v == 0 {
            Strand::Fwd
        } else {
            Strand::Rev
        }
    }
    #[inline]
    pub fn symbol(self) -> char {
        match self {
            Strand::Fwd => '+',
            Strand::Rev => '-',
        }
    }
}

/// One digested locus.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub enzyme_idx: u8,
    pub contig_id: u16,
    /// 0-based coordinate of the first base of the recognition site.
    pub site_start: u64,
    pub strand: Strand,
    /// Half-open tag span in contig coordinates.
    pub tag_start: u64,
    pub tag_end: u64,
    /// The excised tag, already oriented 5'->3' along its own strand.
    pub tag: Vec<u8>,
}

impl Site {
    /// Canonical (strand-independent) hash of the tag sequence.
    pub fn tag_hash(&self) -> u64 {
        // The tag is already strand-oriented, so hashing it directly makes the
        // forward and reverse copy of a palindromic locus agree.
        let f = hash_tag(&self.tag);
        let r = hash_tag(&revcomp(&self.tag));
        f.min(r)
    }
}

/// Cross-enzyme union merge radius, in bp. Two sites from different enzymes
/// whose recognition starts fall within this distance are one effective locus
/// (report §4.1: "跨酶并集按 33 bp 内合并为同一有效位点"); 33 bp is the longest
/// tag in the panel (CspCI).
pub const DEFAULT_MERGE_WINDOW: u64 = 33;

/// Digestion options.
#[derive(Debug, Clone)]
pub struct DigestConfig {
    /// Skip tags containing an ambiguous base. Reference `N` runs would
    /// otherwise produce anchors that no read can ever match.
    pub reject_ambiguous_tags: bool,
    /// Skip contigs shorter than this (bp). Short contigs contribute anchors
    /// whose position carries no ori-ter information.
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

/// Digest one contig with one enzyme.
///
/// Both strands are searched. For a palindromic recognition pattern the reverse
/// hit is at the same locus as the forward hit, so only the forward strand is
/// scanned — this is the deduplication the report requires for AlfI, BplI, FalI,
/// HaeIV and friends.
pub fn digest_contig_with(
    seq: &[u8],
    contig_id: u16,
    enzyme: &Enzyme,
    cfg: &DigestConfig,
    out: &mut Vec<Site>,
) {
    if seq.len() < cfg.min_contig_len || seq.len() < enzyme.tag_len as usize {
        return;
    }
    let fwd_pat = enzyme.pattern_bytes();
    let rc_pat = revcomp(fwd_pat);
    let palindromic = enzyme.is_palindromic();
    let plen = fwd_pat.len();

    for pos in 0..=(seq.len() - plen) {
        // Forward-strand hit.
        if enzyme.matches_at(seq, pos) {
            push_site(seq, contig_id, enzyme, cfg, pos, Strand::Fwd, out);
        }
        // Reverse-strand hit: the reverse complement of the pattern occurring on
        // the forward strand. Skipped for palindromes, which would double-count.
        if !palindromic && matches_pattern(&rc_pat, seq, pos) {
            push_site(seq, contig_id, enzyme, cfg, pos, Strand::Rev, out);
        }
    }
}

#[inline]
fn matches_pattern(pat: &[u8], seq: &[u8], pos: usize) -> bool {
    if pos + pat.len() > seq.len() {
        return false;
    }
    pat.iter()
        .zip(&seq[pos..pos + pat.len()])
        .all(|(&code, &base)| crate::seq::iupac_matches(code, base))
}

fn push_site(
    seq: &[u8],
    contig_id: u16,
    enzyme: &Enzyme,
    cfg: &DigestConfig,
    site_start: usize,
    strand: Strand,
    out: &mut Vec<Site>,
) {
    let Some((ts, te)) = enzyme.tag_span(site_start, seq.len(), strand == Strand::Fwd) else {
        return; // tag would run off the contig end
    };
    let raw = &seq[ts..te];
    if cfg.reject_ambiguous_tags && raw.iter().any(|b| !matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return;
    }
    let tag = match strand {
        Strand::Fwd => raw.to_vec(),
        Strand::Rev => revcomp(raw),
    };
    out.push(Site {
        enzyme_idx: enzyme.idx,
        contig_id,
        site_start: site_start as u64,
        strand,
        tag_start: ts as u64,
        tag_end: te as u64,
        tag,
    });
}

/// Digest one contig with a whole enzyme selection. Sites come back sorted by
/// `(site_start, enzyme_idx, strand)` so downstream windowing can stream them.
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
    sites.sort_by_key(|s| (s.site_start, s.enzyme_idx, s.strand));
    sites
}

/// Per-enzyme site counts and the merged union count for one contig set.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DensityReport {
    /// Total sequence length considered, in bp.
    pub genome_len: u64,
    /// Site count per panel index.
    pub per_enzyme: Vec<(String, u64)>,
    /// Independent loci after merging within [`DEFAULT_MERGE_WINDOW`].
    pub union_sites: u64,
    /// Mean spacing between consecutive union loci, in bp.
    pub mean_spacing: f64,
    /// Largest gap between consecutive union loci, in bp — the "blind spot"
    /// statistic the report uses to justify the 16-enzyme overlay.
    pub max_gap: u64,
    /// Union loci per 25 kb, comparable to Pilea's ~100 k-mers per window.
    pub per_25kb: f64,
}

impl DensityReport {
    /// Sites per Mb for one enzyme, matching the units of report table §4.1.
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

/// Collapse sites into independent loci: sorted by position, merging any two
/// whose recognition starts are within `merge_window` bp *on the same contig*.
///
/// Returns the merged locus start coordinates per contig, which is what the
/// spacing/gap statistics are computed over.
pub fn merge_union(sites: &[Site], merge_window: u64) -> Vec<(u16, u64)> {
    let mut keys: Vec<(u16, u64)> = sites.iter().map(|s| (s.contig_id, s.site_start)).collect();
    keys.sort_unstable();
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
/// `contig_lens` must be parallel to contig ids. Spacing statistics are computed
/// within each contig; the gap across a contig boundary is not a real genomic
/// gap and would inflate `max_gap` on fragmented MAGs.
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

    let merged = merge_union(sites, merge_window);
    let mut max_gap = 0u64;
    let mut gap_sum = 0u64;
    let mut gap_n = 0u64;
    for w in merged.windows(2) {
        let ((c0, p0), (c1, p1)) = (w[0], w[1]);
        if c0 != c1 {
            continue; // contig boundary is not a genomic gap
        }
        let g = p1 - p0;
        max_gap = max_gap.max(g);
        gap_sum += g;
        gap_n += 1;
    }
    DensityReport {
        genome_len,
        per_enzyme,
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

    fn cfg() -> DigestConfig {
        DigestConfig {
            reject_ambiguous_tags: true,
            min_contig_len: 0,
        }
    }

    #[test]
    fn finds_a_planted_bcgi_site() {
        let bcgi = by_name("BcgI").unwrap();
        // 20 bp lead-in, then CGA-N6-TGC, then 20 bp tail.
        let seq = b"AAAAAAAAAAAAAAAAAAAACGAACGTACTGCTTTTTTTTTTTTTTTTTTTT".to_vec();
        let sites = digest_contig(&seq, 0, &[bcgi], &cfg());
        assert_eq!(sites.len(), 1);
        let s = &sites[0];
        assert_eq!(s.site_start, 20);
        assert_eq!(s.strand, Strand::Fwd);
        assert_eq!(s.tag.len(), 32);
        // 10 bp up-flank + 12 bp site + 10 bp down-flank
        assert_eq!((s.tag_start, s.tag_end), (10, 42));
    }

    #[test]
    fn finds_the_reverse_strand_copy() {
        let bcgi = by_name("BcgI").unwrap();
        // revcomp(CGA-N6-TGC) = GCA-N6-TCG planted on the forward strand.
        let seq = b"AAAAAAAAAAAAAAAAAAAAGCAACGTACTCGTTTTTTTTTTTTTTTTTTTT".to_vec();
        let sites = digest_contig(&seq, 0, &[bcgi], &cfg());
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].strand, Strand::Rev);
        // The stored tag is oriented along its own strand, so it reads as a
        // forward BcgI tag.
        let t = &sites[0].tag;
        assert!(by_name("BcgI").unwrap().matches_at(t, 10));
    }

    #[test]
    fn palindromic_enzyme_is_not_double_counted() {
        let plai = by_name("BplI").unwrap(); // GAG-N5-CTC is its own revcomp
        assert!(plai.is_palindromic());
        let seq = b"AAAAAAAAAAAAAAAAAAAAGAGACGTACTCTTTTTTTTTTTTTTTTTTTTT".to_vec();
        let sites = digest_contig(&seq, 0, &[plai], &cfg());
        assert_eq!(sites.len(), 1, "palindromic hit reported twice");
    }

    #[test]
    fn tags_with_n_are_rejected() {
        let bcgi = by_name("BcgI").unwrap();
        // N run at 10..15, i.e. inside the 10..42 tag span but outside the site.
        let seq = b"AAAAAAAAAANNNNNAAAAACGAACGTACTGCTTTTTTTTTTTTTTTTTTTT".to_vec();
        assert!(digest_contig(&seq, 0, &[bcgi], &cfg()).is_empty());
        let lax = DigestConfig {
            reject_ambiguous_tags: false,
            min_contig_len: 0,
        };
        assert_eq!(digest_contig(&seq, 0, &[bcgi], &lax).len(), 1);
    }

    #[test]
    fn sites_near_contig_ends_are_dropped() {
        let bcgi = by_name("BcgI").unwrap();
        // Site at position 2: no room for the 10 bp upstream flank.
        let seq = b"AACGAACGTACTGCTTTTTTTTTTTTTTTTTTTT".to_vec();
        assert!(digest_contig(&seq, 0, &[bcgi], &cfg()).is_empty());
    }

    #[test]
    fn union_merges_overlapping_loci() {
        let sites = vec![
            mk_site(0, 0, 100),
            mk_site(1, 0, 110), // within 33 bp of the previous -> same locus
            mk_site(2, 0, 400),
            mk_site(3, 1, 105), // different contig -> never merged
        ];
        let merged = merge_union(&sites, DEFAULT_MERGE_WINDOW);
        assert_eq!(merged, vec![(0, 100), (0, 400), (1, 105)]);
    }

    #[test]
    fn density_report_ignores_contig_boundaries_in_gaps() {
        let sites = vec![mk_site(0, 0, 100), mk_site(0, 0, 600), mk_site(0, 1, 50)];
        let e = vec![by_name("BcgI").unwrap()];
        let rep = density_report(&sites, &[1_000, 1_000], &e, DEFAULT_MERGE_WINDOW);
        assert_eq!(rep.union_sites, 3);
        assert_eq!(rep.max_gap, 500, "cross-contig gap leaked into max_gap");
        assert_eq!(rep.genome_len, 2_000);
        assert_eq!(rep.density_per_mb("BcgI"), Some(1_500.0));
    }

    fn mk_site(enzyme_idx: u8, contig_id: u16, site_start: u64) -> Site {
        Site {
            enzyme_idx,
            contig_id,
            site_start,
            strand: Strand::Fwd,
            tag_start: site_start,
            tag_end: site_start + 27,
            tag: b"ACGTACGTACGTACGTACGTACGTACG".to_vec(),
        }
    }
}
