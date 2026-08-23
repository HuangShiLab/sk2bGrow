//! EM reassignment of anchors shared between genomes.
//!
//! Reuse target: the containment-driven shared-k-mer reassignment Pilea borrows
//! from sylph (`sketch.py`), with anchors substituted for k-mers.
//!
//! The problem: a tag present in two references cannot be attributed to either
//! by sequence alone. Discarding such tags throws away signal in exactly the
//! situation where it is scarcest (closely related strains); assigning them to
//! the best-matching genome outright biases that genome's coverage upward. The
//! EM compromise splits each shared count in proportion to the current abundance
//! estimates, then re-estimates abundance, until the split stops moving.
//!
//! Genome abundance is estimated from **unique anchors only** on every
//! iteration. That keeps the fixed point identifiable: if shared anchors fed
//! back into the abundance they are then used to split, a genome with no unique
//! anchors could inflate itself without bound.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::anchor_db::AnchorDb;

#[derive(Debug, Clone)]
pub struct EmConfig {
    pub max_iter: usize,
    /// Stop when the largest relative change in any genome abundance is below
    /// this.
    pub tol: f64,
    /// Pseudo-count added to every genome's abundance, so a genome with zero
    /// observed unique anchors still receives a vanishing rather than undefined
    /// share.
    pub prior: f64,
}

impl Default for EmConfig {
    fn default() -> Self {
        EmConfig {
            max_iter: 100,
            tol: 1e-6,
            prior: 1e-9,
        }
    }
}

/// Per-genome summary produced alongside the reassignment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenomeAbundance {
    pub genome_id: u32,
    /// Mean count per unique anchor — the coverage proxy fed to the statistics
    /// layer as `coverage`.
    pub lambda: f64,
    /// Fraction of this genome's usable anchors with a non-zero count. This is
    /// the `containment` column, and Pilea's 95 % ANI cut-off is a threshold on
    /// exactly this quantity (report D5).
    pub containment: f64,
    pub n_unique_anchors: usize,
    pub n_detected_anchors: usize,
    /// Total count mass assigned to this genome after reassignment.
    pub assigned_mass: f64,
}

#[derive(Debug, Clone)]
pub struct EmResult {
    /// Fractional count per anchor, parallel to `db.anchors`. Unique anchors
    /// keep their raw count; shared anchors hold their assigned share.
    pub weights: Vec<f64>,
    pub genomes: Vec<GenomeAbundance>,
    pub iterations: usize,
    pub converged: bool,
}

/// Run the reassignment.
///
/// `counts` must be parallel to `db.anchors`.
pub fn reassign(db: &AnchorDb, counts: &[u32], cfg: &EmConfig) -> EmResult {
    assert_eq!(
        counts.len(),
        db.anchors.len(),
        "counts must be parallel to anchors"
    );

    let groups = db.shared_groups();
    let mut weights: Vec<f64> = counts.iter().map(|&c| c as f64).collect();
    // Shared anchors start unassigned; their mass is distributed each E step.
    for idxs in groups.values() {
        for &i in idxs {
            weights[i] = 0.0;
        }
    }

    let genome_ids: Vec<u32> = db.genomes.iter().map(|g| g.id).collect();
    let mut lambda: HashMap<u32, f64> = HashMap::new();
    let mut converged = false;
    let mut iterations = 0usize;

    // Initial abundance from unique anchors alone.
    for &g in &genome_ids {
        lambda.insert(g, unique_lambda(db, counts, g) + cfg.prior);
    }

    for it in 1..=cfg.max_iter {
        iterations = it;
        // --- E step: split each shared group by current abundance ------------
        for idxs in groups.values() {
            // Every member of a group carries the same tag, so the *observed*
            // count is the same number recorded once per member. Take it once.
            let observed = idxs
                .iter()
                .map(|&i| counts[i] as f64)
                .fold(0.0f64, f64::max);
            if observed == 0.0 {
                for &i in idxs {
                    weights[i] = 0.0;
                }
                continue;
            }
            let denom: f64 = idxs.iter().map(|&i| lambda[&db.anchors[i].genome_id]).sum();
            for &i in idxs {
                let share = if denom > 0.0 {
                    lambda[&db.anchors[i].genome_id] / denom
                } else {
                    1.0 / idxs.len() as f64
                };
                weights[i] = observed * share;
            }
        }

        // --- M step: re-estimate abundance from unique anchors ---------------
        let mut max_rel = 0.0f64;
        for &g in &genome_ids {
            let new = unique_lambda(db, counts, g) + cfg.prior;
            let old = lambda[&g];
            let rel = if old > 0.0 {
                (new - old).abs() / old
            } else {
                new.abs()
            };
            max_rel = max_rel.max(rel);
            lambda.insert(g, new);
        }
        if max_rel < cfg.tol {
            converged = true;
            break;
        }
    }

    let genomes = genome_ids
        .iter()
        .map(|&g| {
            let range = db.genome_range(g);
            let mut n_unique = 0usize;
            let mut n_detected = 0usize;
            let mut mass = 0.0f64;
            for i in range {
                let a = &db.anchors[i];
                mass += weights[i];
                if a.is_usable() {
                    n_unique += 1;
                    if counts[i] > 0 {
                        n_detected += 1;
                    }
                }
            }
            GenomeAbundance {
                genome_id: g,
                lambda: lambda[&g] - cfg.prior,
                containment: if n_unique > 0 {
                    n_detected as f64 / n_unique as f64
                } else {
                    0.0
                },
                n_unique_anchors: n_unique,
                n_detected_anchors: n_detected,
                assigned_mass: mass,
            }
        })
        .collect();

    EmResult {
        weights,
        genomes,
        iterations,
        converged,
    }
}

/// Mean count over a genome's usable (unique, chromosomal) anchors.
fn unique_lambda(db: &AnchorDb, counts: &[u32], genome_id: u32) -> f64 {
    let range = db.genome_range(genome_id);
    let mut n = 0usize;
    let mut s = 0u64;
    for i in range {
        if db.anchors[i].is_usable() {
            n += 1;
            s += counts[i] as u64;
        }
    }
    if n == 0 {
        0.0
    } else {
        s as f64 / n as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchor_db::{Anchor, BuildParams, GenomeMeta};
    use crate::tgt::{ContigKind, ContigMeta};

    /// `spec` is (seq_hash, genome_id).
    fn db_of(spec: &[(u64, u32)]) -> AnchorDb {
        let anchors: Vec<Anchor> = spec
            .iter()
            .enumerate()
            .map(|(i, &(h, g))| Anchor {
                seq_hash: h,
                genome_id: g,
                contig_id: 0,
                position: i as u64 * 100,
                enzyme_idx: 0,
                strand: 0,
                flags: 0,
                local_gc: 100,
            })
            .collect();
        let mut ids: Vec<u32> = spec.iter().map(|&(_, g)| g).collect();
        ids.sort_unstable();
        ids.dedup();
        let n = anchors.len();
        let mut db = AnchorDb {
            params: BuildParams::default(),
            genomes: ids
                .into_iter()
                .map(|id| GenomeMeta {
                    id,
                    name: format!("g{id}"),
                    taxonomy: None,
                    contigs: vec![ContigMeta {
                        id: 0,
                        name: "c0".into(),
                        length: 100_000,
                        offset: 0,
                        kind: ContigKind::Chromosome,
                    }],
                    genome_len: 100_000,
                    ori: None,
                    ori_confidence: 0.0,
                })
                .collect(),
            anchors,
            tags: vec![[0u8; 12]; n],
        };
        db.anchors.sort_by_key(|a| (a.genome_id, a.position));
        db.recompute_uniqueness();
        db
    }

    #[test]
    fn shared_mass_splits_by_abundance() {
        // g0: 3 unique anchors at depth 10; g1: 3 unique at depth 2.
        // One anchor (hash 99) is shared, observed 12 times.
        let db = db_of(&[
            (1, 0),
            (2, 0),
            (3, 0),
            (99, 0),
            (11, 1),
            (12, 1),
            (13, 1),
            (99, 1),
        ]);
        let counts: Vec<u32> = db
            .anchors
            .iter()
            .map(|a| {
                if a.seq_hash == 99 {
                    12
                } else if a.genome_id == 0 {
                    10
                } else {
                    2
                }
            })
            .collect();
        let r = reassign(&db, &counts, &EmConfig::default());
        assert!(
            r.converged,
            "EM failed to converge in {} iterations",
            r.iterations
        );

        let shared: Vec<(u32, f64)> = db
            .anchors
            .iter()
            .zip(&r.weights)
            .filter(|(a, _)| a.seq_hash == 99)
            .map(|(a, &w)| (a.genome_id, w))
            .collect();
        let total: f64 = shared.iter().map(|&(_, w)| w).sum();
        assert!(
            (total - 12.0).abs() < 1e-9,
            "shared mass was not conserved: {total}"
        );
        let g0 = shared.iter().find(|&&(g, _)| g == 0).unwrap().1;
        // 10 : 2 abundance ratio -> 10/12 of the shared count.
        assert!((g0 - 10.0).abs() < 1e-6, "g0 share was {g0}, expected 10");
    }

    #[test]
    fn unique_anchors_keep_their_raw_counts() {
        let db = db_of(&[(1, 0), (2, 0), (99, 0), (99, 1), (11, 1)]);
        let counts: Vec<u32> = db
            .anchors
            .iter()
            .map(|a| if a.seq_hash == 99 { 8 } else { 5 })
            .collect();
        let r = reassign(&db, &counts, &EmConfig::default());
        for (a, &w) in db.anchors.iter().zip(&r.weights) {
            if a.seq_hash != 99 {
                assert_eq!(w, 5.0);
            }
        }
    }

    #[test]
    fn containment_counts_only_usable_anchors() {
        let db = db_of(&[(1, 0), (2, 0), (3, 0), (4, 0)]);
        // Two of the four detected.
        let counts = vec![7u32, 0, 3, 0];
        let r = reassign(&db, &counts, &EmConfig::default());
        let g = &r.genomes[0];
        assert_eq!(g.n_unique_anchors, 4);
        assert_eq!(g.n_detected_anchors, 2);
        assert!((g.containment - 0.5).abs() < 1e-12);
        assert!((g.lambda - 2.5).abs() < 1e-12);
    }

    #[test]
    fn an_undetected_shared_anchor_stays_at_zero() {
        let db = db_of(&[(1, 0), (99, 0), (99, 1), (11, 1)]);
        let counts: Vec<u32> = db
            .anchors
            .iter()
            .map(|a| if a.seq_hash == 99 { 0 } else { 4 })
            .collect();
        let r = reassign(&db, &counts, &EmConfig::default());
        for (a, &w) in db.anchors.iter().zip(&r.weights) {
            if a.seq_hash == 99 {
                assert_eq!(w, 0.0);
            }
        }
    }

    #[test]
    fn a_genome_with_no_unique_anchors_does_not_explode() {
        // g1's only anchor is shared with g0.
        let db = db_of(&[(1, 0), (2, 0), (99, 0), (99, 1)]);
        // Uniform depth: the point is g1's *lack of unique anchors*, not any
        // depth difference between the two genomes.
        let counts: Vec<u32> = vec![6; db.anchors.len()];
        let r = reassign(&db, &counts, &EmConfig::default());
        let g1 = r.genomes.iter().find(|g| g.genome_id == 1).unwrap();
        assert_eq!(g1.n_unique_anchors, 0);
        assert_eq!(g1.lambda, 0.0);
        assert!(g1.assigned_mass.is_finite());
        // With zero abundance g1 takes a vanishing share, and g0 takes the rest.
        assert!(
            g1.assigned_mass < 1e-6,
            "g1 absorbed {} of the shared mass",
            g1.assigned_mass
        );
    }
}
