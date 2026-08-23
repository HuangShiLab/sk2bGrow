//! Order and orient draft (MAG) contigs against a reference, using shared tags.
//!
//! Reuse target: Syn2b's `scaffold` subcommand. PTR estimation on a fragmented
//! MAG is hard for the obvious reason that contig order is unknown; the report
//! (§6.3) points out that shared 2bRAD tags already solve this, so scaffolding
//! first and fitting second turns the hardest case into a solved one — and does
//! it from a *single* sample, unlike the multi-sample reordering DEMIC and
//! CoPTR-Contig need.
//!
//! The placement is deliberately coarse. It only has to be good enough to give
//! every anchor a global coordinate; a few kb of positional error is negligible
//! against a megabase-scale ori-ter gradient.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::anchor_db::AnchorDb;
use crate::tgt::{ContigMeta, Tgt};

/// Orientation of a draft contig relative to the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Orientation {
    Forward,
    Reverse,
}

/// Where one draft contig landed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Placement {
    pub contig_id: u16,
    pub contig_len: u64,
    /// Reference coordinate the contig start maps to.
    pub ref_start: u64,
    pub orientation: Orientation,
    /// Tags shared with the reference.
    pub n_tags: usize,
    /// Share of shared tags whose ordering agrees with the chosen orientation.
    /// Below ~0.9 the contig is likely chimaeric or repeat-driven.
    pub concordance: f64,
}

#[derive(Debug, Clone)]
pub struct ScaffoldConfig {
    /// Minimum shared tags before a contig may be placed. Syn2b's default is 3.
    pub min_tags: usize,
    /// Minimum concordance for a placement to be accepted.
    pub min_concordance: f64,
}

impl Default for ScaffoldConfig {
    fn default() -> Self {
        ScaffoldConfig {
            min_tags: 3,
            min_concordance: 0.8,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ScaffoldResult {
    pub placements: Vec<Placement>,
    /// Contigs with too few shared tags, or too little agreement, to place.
    /// They keep their input order and are appended after the placed ones.
    pub unplaced: Vec<u16>,
}

impl ScaffoldResult {
    pub fn placed_bp(&self, draft: &Tgt) -> u64 {
        self.placements
            .iter()
            .filter_map(|p| draft.contigs.iter().find(|c| c.id == p.contig_id))
            .map(|c| c.length)
            .sum()
    }
    /// Fraction of draft sequence that received a reference coordinate.
    pub fn placed_fraction(&self, draft: &Tgt) -> f64 {
        let total: u64 = draft.contigs.iter().map(|c| c.length).sum();
        if total == 0 {
            0.0
        } else {
            self.placed_bp(draft) as f64 / total as f64
        }
    }
}

/// Scaffold `draft` against the anchors of `reference_genome_id` in `db`.
pub fn scaffold(
    draft: &Tgt,
    db: &AnchorDb,
    reference_genome_id: u32,
    cfg: &ScaffoldConfig,
) -> ScaffoldResult {
    // Reference tag hash -> global coordinate. Multi-copy tags are dropped: a
    // repeat would vote for several incompatible placements at once.
    let mut ref_pos: HashMap<u64, u64> = HashMap::new();
    let mut ambiguous: Vec<u64> = Vec::new();
    for a in &db.anchors[db.genome_range(reference_genome_id)] {
        let Some(gp) = db.global_position(a) else {
            continue;
        };
        if ref_pos.insert(a.seq_hash, gp).is_some() {
            ambiguous.push(a.seq_hash);
        }
    }
    for h in ambiguous {
        ref_pos.remove(&h);
    }

    let mut result = ScaffoldResult::default();
    for contig in &draft.contigs {
        let pairs: Vec<(u64, u64)> = draft
            .records
            .iter()
            .filter(|r| r.contig_id == contig.id)
            .filter_map(|r| ref_pos.get(&r.tag_hash).map(|&rp| (r.position, rp)))
            .collect();

        if pairs.len() < cfg.min_tags {
            result.unplaced.push(contig.id);
            continue;
        }
        let (orientation, concordance) = orient(&pairs);
        if concordance < cfg.min_concordance {
            result.unplaced.push(contig.id);
            continue;
        }
        // Project each shared tag back to where the contig's *start* would be,
        // then take the median. The median rather than the mean because a single
        // spurious tag match would drag a mean across the chromosome.
        let mut starts: Vec<i64> = pairs
            .iter()
            .map(|&(dp, rp)| match orientation {
                Orientation::Forward => rp as i64 - dp as i64,
                Orientation::Reverse => rp as i64 - (contig.length as i64 - 1 - dp as i64),
            })
            .collect();
        starts.sort_unstable();
        let ref_start = starts[starts.len() / 2].max(0) as u64;

        result.placements.push(Placement {
            contig_id: contig.id,
            contig_len: contig.length,
            ref_start,
            orientation,
            n_tags: pairs.len(),
            concordance,
        });
    }
    result
        .placements
        .sort_by_key(|p| (p.ref_start, p.contig_id));
    result
}

/// Decide orientation by counting concordant vs discordant ordered pairs — a
/// Kendall-style vote. Robust to a few bad matches, and needs no regression.
fn orient(pairs: &[(u64, u64)]) -> (Orientation, f64) {
    let mut concordant = 0usize;
    let mut discordant = 0usize;
    for i in 0..pairs.len() {
        for j in (i + 1)..pairs.len() {
            let (d0, r0) = pairs[i];
            let (d1, r1) = pairs[j];
            if d0 == d1 || r0 == r1 {
                continue;
            }
            if (d1 > d0) == (r1 > r0) {
                concordant += 1;
            } else {
                discordant += 1;
            }
        }
    }
    let total = concordant + discordant;
    if total == 0 {
        return (Orientation::Forward, 0.0);
    }
    if concordant >= discordant {
        (Orientation::Forward, concordant as f64 / total as f64)
    } else {
        (Orientation::Reverse, discordant as f64 / total as f64)
    }
}

/// Rewrite contig offsets from a scaffold result.
///
/// Placed contigs take their reference coordinate; unplaced contigs are parked
/// beyond the end of the placed region so their anchors keep distinct global
/// coordinates without being interleaved into the fit. The statistics layer
/// drops anchors whose contig was not placed — a wrong coordinate is worse than
/// a missing one for a gradient fit.
pub fn apply(draft: &mut Tgt, result: &ScaffoldResult) {
    let placed: HashMap<u16, &Placement> =
        result.placements.iter().map(|p| (p.contig_id, p)).collect();
    let mut park = result
        .placements
        .iter()
        .map(|p| p.ref_start + p.contig_len)
        .max()
        .unwrap_or(0);
    for c in draft.contigs.iter_mut() {
        match placed.get(&c.id) {
            Some(p) => c.offset = p.ref_start,
            None => {
                c.offset = park;
                park += c.length;
            }
        }
    }
    // A reverse-oriented contig has its internal coordinates mirrored so that
    // ascending position still means ascending reference coordinate.
    for p in &result.placements {
        if p.orientation != Orientation::Reverse {
            continue;
        }
        for r in draft
            .records
            .iter_mut()
            .filter(|r| r.contig_id == p.contig_id)
        {
            r.position = p.contig_len.saturating_sub(1).saturating_sub(r.position);
        }
    }
    draft.recompute_gaps();
}

/// Contigs the scaffold could not place, as a set of ids for masking.
pub fn unplaced_contigs(result: &ScaffoldResult) -> Vec<u16> {
    result.unplaced.clone()
}

/// Convenience: contig metadata sorted by final offset.
pub fn ordered_contigs(draft: &Tgt) -> Vec<&ContigMeta> {
    let mut v: Vec<&ContigMeta> = draft.contigs.iter().collect();
    v.sort_by_key(|c| c.offset);
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchor_db::{Anchor, BuildParams, GenomeMeta};
    use crate::tgt::{ContigKind, TgtRecord};

    fn ref_db(tags: &[(u64, u64)]) -> AnchorDb {
        // tags: (hash, position) on a single 1 Mb reference contig.
        let anchors: Vec<Anchor> = tags
            .iter()
            .map(|&(h, p)| Anchor {
                seq_hash: h,
                genome_id: 0,
                contig_id: 0,
                position: p,
                enzyme_idx: 0,
                strand: 0,
                flags: 0,
                local_gc: 100,
            })
            .collect();
        let n = anchors.len();
        let mut db = AnchorDb {
            params: BuildParams::default(),
            genomes: vec![GenomeMeta {
                id: 0,
                name: "ref".into(),
                taxonomy: None,
                contigs: vec![ContigMeta {
                    id: 0,
                    name: "chr".into(),
                    length: 1_000_000,
                    offset: 0,
                    kind: ContigKind::Chromosome,
                }],
                genome_len: 1_000_000,
                ori: None,
                ori_confidence: 0.0,
            }],
            anchors,
            tags: vec![[0u8; 12]; n],
        };
        db.recompute_uniqueness();
        db
    }

    fn draft_of(contigs: &[(u16, u64)], recs: &[(u16, u64, u64)]) -> Tgt {
        let mut t = Tgt::new(
            "draft",
            contigs
                .iter()
                .map(|&(id, len)| ContigMeta {
                    id,
                    name: format!("ctg{id}"),
                    length: len,
                    offset: 0,
                    kind: ContigKind::Chromosome,
                })
                .collect(),
        );
        t.records = recs
            .iter()
            .map(|&(contig_id, pos, hash)| TgtRecord {
                tag_hash: hash,
                position: pos,
                gap: 0,
                contig_id,
                enzyme_idx: 0,
                pattern: 0,
                tag_len: 32,
                flags: 0,
                local_gc: 100,
                tag_2bit: [0u8; 12],
            })
            .collect();
        t
    }

    #[test]
    fn places_a_forward_contig() {
        let db = ref_db(&[(1, 100_000), (2, 101_000), (3, 102_000)]);
        let draft = draft_of(&[(0, 5_000)], &[(0, 0, 1), (0, 1_000, 2), (0, 2_000, 3)]);
        let r = scaffold(&draft, &db, 0, &ScaffoldConfig::default());
        assert_eq!(r.placements.len(), 1);
        let p = &r.placements[0];
        assert_eq!(p.orientation, Orientation::Forward);
        assert_eq!(p.ref_start, 100_000);
        assert_eq!(p.concordance, 1.0);
    }

    #[test]
    fn detects_an_inverted_contig() {
        let db = ref_db(&[(1, 100_000), (2, 101_000), (3, 102_000)]);
        // Draft coordinates run the other way: tag 3 first, tag 1 last.
        let draft = draft_of(&[(0, 3_000)], &[(0, 0, 3), (0, 1_000, 2), (0, 2_000, 1)]);
        let r = scaffold(&draft, &db, 0, &ScaffoldConfig::default());
        let p = &r.placements[0];
        assert_eq!(p.orientation, Orientation::Reverse);
        assert_eq!(p.concordance, 1.0);
    }

    #[test]
    fn too_few_shared_tags_leaves_a_contig_unplaced() {
        let db = ref_db(&[(1, 100_000), (2, 101_000)]);
        let draft = draft_of(
            &[(0, 5_000), (1, 5_000)],
            &[(0, 0, 1), (0, 100, 2), (1, 0, 999)],
        );
        let r = scaffold(&draft, &db, 0, &ScaffoldConfig::default());
        assert_eq!(r.unplaced, vec![0, 1], "min_tags=3 should reject both");
        assert!(r.placements.is_empty());
    }

    #[test]
    fn repeated_reference_tags_are_ignored() {
        // Hash 1 occurs twice in the reference: it must not vote.
        let db = ref_db(&[
            (1, 100_000),
            (1, 900_000),
            (2, 101_000),
            (3, 102_000),
            (4, 103_000),
        ]);
        let draft = draft_of(
            &[(0, 5_000)],
            &[(0, 0, 1), (0, 1_000, 2), (0, 2_000, 3), (0, 3_000, 4)],
        );
        let r = scaffold(&draft, &db, 0, &ScaffoldConfig::default());
        let p = &r.placements[0];
        assert_eq!(p.n_tags, 3, "the repeated tag was counted");
        assert_eq!(p.ref_start, 100_000);
    }

    #[test]
    fn apply_sets_offsets_and_mirrors_reverse_contigs() {
        let db = ref_db(&[(1, 100_000), (2, 101_000), (3, 102_000)]);
        let mut draft = draft_of(
            &[(0, 3_000), (1, 500)],
            &[(0, 0, 3), (0, 1_000, 2), (0, 2_000, 1), (1, 10, 777)],
        );
        let r = scaffold(&draft, &db, 0, &ScaffoldConfig::default());
        apply(&mut draft, &r);
        let c0 = draft.contigs.iter().find(|c| c.id == 0).unwrap();
        assert_eq!(c0.offset, r.placements[0].ref_start);
        // Unplaced contig parked after the placed region.
        let c1 = draft.contigs.iter().find(|c| c.id == 1).unwrap();
        assert!(c1.offset >= c0.offset + c0.length);
        // The reverse contig's positions were mirrored: the tag that was at 0 is
        // now at len-1.
        let mirrored: Vec<u64> = draft
            .records
            .iter()
            .filter(|r| r.contig_id == 0)
            .map(|r| r.position)
            .collect();
        assert!(mirrored.contains(&2_999));
        assert_eq!(r.placed_fraction(&draft), 3_000.0 / 3_500.0);
    }
}
