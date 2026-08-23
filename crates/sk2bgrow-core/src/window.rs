//! Window construction over anchors.
//!
//! Pilea slices a sketch into 25 kb fixed windows, which holds ~100 k-mers on
//! average but varies with the Bernoulli sampling and collapses on fragmented
//! references (report D2). Because TGT anchor positions are deterministic and
//! their spacing is known at build time, the window can instead hold a *fixed
//! number of anchors* — equalising statistical power across the genome
//! (report §6.3).
//!
//! Both policies are implemented: [`WindowPolicy::FixedBp`] exists so an A/B
//! benchmark against Pilea can hold windowing constant and vary only the sketch.

use serde::{Deserialize, Serialize};

use crate::anchor_db::{Anchor, AnchorDb};

/// Sentinel window id for anchors that were not assigned to any window.
pub const NO_WINDOW: u32 = u32::MAX;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum WindowPolicy {
    /// `n` anchors per window (default 100, matching Pilea's mean occupancy).
    EqualAnchors { n: usize },
    /// `bp` base pairs per window (Pilea parity: 25 000).
    FixedBp { bp: u64 },
}

impl Default for WindowPolicy {
    fn default() -> Self {
        WindowPolicy::EqualAnchors { n: 100 }
    }
}

/// A window's extent and membership.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Window {
    pub id: u32,
    pub genome_id: u32,
    pub contig_id: u16,
    /// Inclusive start / exclusive end in contig coordinates.
    pub start: u64,
    pub end: u64,
    /// Midpoint in *global* coordinates — the x-axis of the V-shape fit.
    pub global_mid: u64,
    pub n_anchors: usize,
}

impl Window {
    /// Span in bp. For an equal-anchor window this varies with local anchor
    /// density and is itself a useful QC signal.
    pub fn span(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }
}

/// Assign every anchor of one genome to a window.
///
/// Returns the window table and a per-anchor id vector aligned to
/// `db.anchors[db.genome_range(genome_id)]` offsets — i.e. index `i` of the
/// returned vector corresponds to global anchor index `range.start + i`.
///
/// Windows never straddle a contig boundary: two anchors on different contigs
/// have no defined genomic distance, so pooling them would fabricate a
/// coordinate. On a fragmented MAG this means short contigs produce short
/// windows, which the statistics layer down-weights rather than merges.
pub fn assign_windows(
    db: &AnchorDb,
    genome_id: u32,
    policy: WindowPolicy,
    usable_only: bool,
) -> (Vec<Window>, Vec<u32>) {
    let range = db.genome_range(genome_id);
    let anchors = &db.anchors[range.clone()];
    let mut ids = vec![NO_WINDOW; anchors.len()];
    let mut windows: Vec<Window> = Vec::new();

    let mut i = 0usize;
    while i < anchors.len() {
        let contig = anchors[i].contig_id;
        let contig_end = {
            let mut j = i;
            while j < anchors.len() && anchors[j].contig_id == contig {
                j += 1;
            }
            j
        };
        match policy {
            WindowPolicy::EqualAnchors { n } => assign_equal_anchors(
                db,
                genome_id,
                anchors,
                i,
                contig_end,
                n.max(1),
                usable_only,
                &mut windows,
                &mut ids,
            ),
            WindowPolicy::FixedBp { bp } => assign_fixed_bp(
                db,
                genome_id,
                anchors,
                i,
                contig_end,
                bp.max(1),
                usable_only,
                &mut windows,
                &mut ids,
            ),
        }
        i = contig_end;
    }
    (windows, ids)
}

#[allow(clippy::too_many_arguments)]
fn assign_equal_anchors(
    db: &AnchorDb,
    genome_id: u32,
    anchors: &[Anchor],
    lo: usize,
    hi: usize,
    n: usize,
    usable_only: bool,
    windows: &mut Vec<Window>,
    ids: &mut [u32],
) {
    let members: Vec<usize> = (lo..hi)
        .filter(|&k| !usable_only || anchors[k].is_usable())
        .collect();
    for chunk in members.chunks(n) {
        // A trailing chunk far below the target would have inflated variance;
        // fold it into the previous window rather than reporting a weak one.
        let wid = if chunk.len() * 2 < n
            && !windows.is_empty()
            && windows.last().unwrap().contig_id == anchors[chunk[0]].contig_id
        {
            windows.last().unwrap().id
        } else {
            windows.push(new_window(
                db,
                genome_id,
                anchors,
                chunk,
                windows.len() as u32,
            ));
            windows.last().unwrap().id
        };
        for &k in chunk {
            ids[k] = wid;
        }
        if let Some(w) = windows.iter_mut().find(|w| w.id == wid) {
            extend_window(w, anchors, chunk);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn assign_fixed_bp(
    db: &AnchorDb,
    genome_id: u32,
    anchors: &[Anchor],
    lo: usize,
    hi: usize,
    bp: u64,
    usable_only: bool,
    windows: &mut Vec<Window>,
    ids: &mut [u32],
) {
    let mut bucket: Vec<usize> = Vec::new();
    let mut bucket_idx = u64::MAX;
    for k in lo..hi {
        if usable_only && !anchors[k].is_usable() {
            continue;
        }
        let b = anchors[k].position / bp;
        if b != bucket_idx {
            if !bucket.is_empty() {
                flush_bp_bucket(
                    db, genome_id, anchors, &bucket, bucket_idx, bp, windows, ids,
                );
                bucket.clear();
            }
            bucket_idx = b;
        }
        bucket.push(k);
    }
    if !bucket.is_empty() {
        flush_bp_bucket(
            db, genome_id, anchors, &bucket, bucket_idx, bp, windows, ids,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn flush_bp_bucket(
    db: &AnchorDb,
    genome_id: u32,
    anchors: &[Anchor],
    bucket: &[usize],
    bucket_idx: u64,
    bp: u64,
    windows: &mut Vec<Window>,
    ids: &mut [u32],
) {
    let id = windows.len() as u32;
    let mut w = new_window(db, genome_id, anchors, bucket, id);
    // Fixed-bp windows report their nominal extent, not the extent of the
    // anchors that happen to fall in them — otherwise a sparse window would look
    // narrow and be over-weighted.
    w.start = bucket_idx * bp;
    w.end = w.start + bp;
    w.global_mid = global_of(db, genome_id, w.contig_id, w.start + bp / 2);
    windows.push(w);
    for &k in bucket {
        ids[k] = id;
    }
}

fn new_window(
    db: &AnchorDb,
    genome_id: u32,
    anchors: &[Anchor],
    members: &[usize],
    id: u32,
) -> Window {
    let first = &anchors[members[0]];
    let last = &anchors[*members.last().unwrap()];
    let start = first.position;
    let end = last.position + 1;
    Window {
        id,
        genome_id,
        contig_id: first.contig_id,
        start,
        end,
        global_mid: global_of(db, genome_id, first.contig_id, (start + end) / 2),
        n_anchors: members.len(),
    }
}

fn extend_window(w: &mut Window, anchors: &[Anchor], members: &[usize]) {
    for &k in members {
        w.start = w.start.min(anchors[k].position);
        w.end = w.end.max(anchors[k].position + 1);
    }
    w.n_anchors = w.n_anchors.max(members.len());
}

fn global_of(db: &AnchorDb, genome_id: u32, contig_id: u16, pos: u64) -> u64 {
    db.genome(genome_id)
        .and_then(|g| g.contigs.iter().find(|c| c.id == contig_id))
        .map(|c| c.offset + pos)
        .unwrap_or(pos)
}

/// Spacing diagnostics used by `sk2bgrow audit`: the "how blind is this genome"
/// numbers the design report quotes (mean spacing, max gap).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SpacingStats {
    pub n_anchors: usize,
    pub mean_spacing: f64,
    pub median_spacing: f64,
    pub p99_spacing: u64,
    pub max_gap: u64,
    /// Number of gaps wider than the audit threshold.
    pub n_wide_gaps: usize,
}

/// Compute spacing statistics for one genome, skipping contig boundaries.
pub fn spacing_stats(
    db: &AnchorDb,
    genome_id: u32,
    wide_gap_bp: u64,
    usable_only: bool,
) -> SpacingStats {
    let anchors: Vec<&Anchor> = db.anchors[db.genome_range(genome_id)]
        .iter()
        .filter(|a| !usable_only || a.is_usable())
        .collect();
    let mut gaps: Vec<u64> = Vec::new();
    for w in anchors.windows(2) {
        if w[0].contig_id != w[1].contig_id {
            continue;
        }
        gaps.push(w[1].position.saturating_sub(w[0].position));
    }
    if gaps.is_empty() {
        return SpacingStats {
            n_anchors: anchors.len(),
            ..Default::default()
        };
    }
    gaps.sort_unstable();
    let n = gaps.len();
    SpacingStats {
        n_anchors: anchors.len(),
        mean_spacing: gaps.iter().sum::<u64>() as f64 / n as f64,
        median_spacing: gaps[n / 2] as f64,
        p99_spacing: gaps[((n as f64 * 0.99) as usize).min(n - 1)],
        max_gap: *gaps.last().unwrap(),
        n_wide_gaps: gaps.iter().filter(|&&g| g > wide_gap_bp).count(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchor_db::{BuildParams, GenomeMeta};
    use crate::tgt::{ContigKind, ContigMeta};

    fn db_with(positions: &[(u16, u64)], contig_lens: &[u64]) -> AnchorDb {
        let anchors: Vec<Anchor> = positions
            .iter()
            .enumerate()
            .map(|(i, &(c, p))| Anchor {
                seq_hash: i as u64,
                genome_id: 0,
                contig_id: c,
                position: p,
                enzyme_idx: (i % 2) as u8,
                strand: 0,
                flags: crate::anchor_db::flags::UNIQUE_IN_GENOME
                    | crate::anchor_db::flags::UNIQUE_ACROSS_DB,
                local_gc: 100,
            })
            .collect();
        let mut offset = 0;
        let contigs: Vec<ContigMeta> = contig_lens
            .iter()
            .enumerate()
            .map(|(i, &l)| {
                let c = ContigMeta {
                    id: i as u16,
                    name: format!("c{i}"),
                    length: l,
                    offset,
                    kind: ContigKind::Chromosome,
                };
                offset += l;
                c
            })
            .collect();
        let n = anchors.len();
        AnchorDb {
            params: BuildParams::default(),
            genomes: vec![GenomeMeta {
                id: 0,
                name: "g0".into(),
                taxonomy: None,
                contigs,
                genome_len: contig_lens.iter().sum(),
                ori: None,
                ori_confidence: 0.0,
            }],
            anchors,
            tags: vec![[0u8; 12]; n],
        }
    }

    #[test]
    fn equal_anchor_windows_hold_the_target_count() {
        let pos: Vec<(u16, u64)> = (0..25).map(|i| (0u16, i as u64 * 1_000)).collect();
        let db = db_with(&pos, &[30_000]);
        let (windows, ids) = assign_windows(&db, 0, WindowPolicy::EqualAnchors { n: 10 }, true);
        // 25 anchors, target 10 -> 10 / 10 / 5; the trailing 5 is >= n/2 so it
        // stands on its own.
        assert_eq!(windows.len(), 3);
        assert_eq!(windows[0].n_anchors, 10);
        assert_eq!(windows[2].n_anchors, 5);
        assert!(ids.iter().all(|&w| w != NO_WINDOW));
    }

    #[test]
    fn a_tiny_trailing_chunk_is_folded_back() {
        let pos: Vec<(u16, u64)> = (0..21).map(|i| (0u16, i as u64 * 1_000)).collect();
        let db = db_with(&pos, &[30_000]);
        let (windows, ids) = assign_windows(&db, 0, WindowPolicy::EqualAnchors { n: 10 }, true);
        assert_eq!(windows.len(), 2, "a 1-anchor window was reported");
        assert_eq!(ids[20], windows[1].id);
    }

    #[test]
    fn windows_never_straddle_contigs() {
        let pos = vec![(0u16, 0u64), (0, 100), (0, 200), (1, 0), (1, 100)];
        let db = db_with(&pos, &[1_000, 1_000]);
        let (windows, ids) = assign_windows(&db, 0, WindowPolicy::EqualAnchors { n: 10 }, true);
        assert_eq!(windows.len(), 2, "one window spanned two contigs");
        assert_ne!(ids[2], ids[3]);
        assert_eq!(windows[1].contig_id, 1);
    }

    #[test]
    fn fixed_bp_windows_report_nominal_extent() {
        let pos = vec![(0u16, 10u64), (0, 20), (0, 25_010), (0, 60_000)];
        let db = db_with(&pos, &[100_000]);
        let (windows, ids) = assign_windows(&db, 0, WindowPolicy::FixedBp { bp: 25_000 }, true);
        assert_eq!(windows.len(), 3);
        assert_eq!((windows[0].start, windows[0].end), (0, 25_000));
        assert_eq!((windows[1].start, windows[1].end), (25_000, 50_000));
        // Window 2 (50k-75k) holds the 60 000 anchor; the empty 3rd bucket is
        // simply absent, which is what "fraction of covered windows" measures.
        assert_eq!(windows[2].start, 50_000);
        assert_eq!(ids[3], windows[2].id);
    }

    #[test]
    fn global_mid_uses_contig_offsets() {
        let pos = vec![(1u16, 100u64), (1, 200)];
        let db = db_with(&pos, &[1_000, 1_000]);
        let (windows, _) = assign_windows(&db, 0, WindowPolicy::EqualAnchors { n: 10 }, true);
        assert_eq!(windows[0].global_mid, 1_000 + 150);
    }

    #[test]
    fn spacing_stats_skip_contig_joins() {
        let pos = vec![(0u16, 0u64), (0, 500), (1, 0), (1, 100)];
        let db = db_with(&pos, &[1_000, 1_000]);
        let s = spacing_stats(&db, 0, 400, true);
        assert_eq!(s.n_anchors, 4);
        assert_eq!(s.max_gap, 500);
        assert_eq!(s.n_wide_gaps, 1);
    }
}
