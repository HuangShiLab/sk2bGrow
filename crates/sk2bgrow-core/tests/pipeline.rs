//! End-to-end: synthetic genome -> anchor database -> simulated reads -> counts,
//! EM reassignment and windows.
//!
//! The point is not to re-test the units but to check that the layers agree on
//! coordinates, strands and flags — the places where an integration bug hides.

use sk2bgrow_core::anchor_db::{assemble, build_genome, flags, AnchorDb, BuildParams, GC_FLANK};
use sk2bgrow_core::count::{AnchorIndex, CountMode, CountStats, MatchConfig};
use sk2bgrow_core::digest::DigestConfig;
use sk2bgrow_core::em::{reassign, EmConfig};
use sk2bgrow_core::enzyme::{parse_selection, EnzymeSet};
use sk2bgrow_core::seq::revcomp;
use sk2bgrow_core::window::{assign_windows, spacing_stats, WindowPolicy};

/// Deterministic xorshift64* — a fixed genome without a dependency on `rand`.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

fn random_genome(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = Rng(seed);
    const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|_| B[(rng.next_u64() >> 33) as usize % 4])
        .collect()
}

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let p = std::env::temp_dir().join(format!("sk2bgrow-it-{}-{}", std::process::id(), tag));
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn write_fasta(dir: &std::path::Path, name: &str, seq: &[u8]) -> std::path::PathBuf {
    let p = dir.join(format!("{name}.fna"));
    let mut body = format!(">{name}_chr\n");
    for chunk in seq.chunks(80) {
        body.push_str(std::str::from_utf8(chunk).unwrap());
        body.push('\n');
    }
    std::fs::write(&p, body).unwrap();
    p
}

fn build_db(dir: &std::path::Path, genomes: &[(&str, Vec<u8>)], enzymes: &str) -> AnchorDb {
    let sel = parse_selection(enzymes).unwrap();
    let cfg = DigestConfig::default();
    let parts: Vec<_> = genomes
        .iter()
        .enumerate()
        .map(|(i, (name, seq))| {
            let path = write_fasta(dir, name, seq);
            let (meta, anchors, tags, _) =
                build_genome(&path, i as u32, &sel, &cfg, GC_FLANK).unwrap();
            (meta, anchors, tags)
        })
        .collect();
    let params = BuildParams {
        enzymes: EnzymeSet::from_slice(&sel),
        ..BuildParams::default()
    };
    assemble(params, parts)
}

/// Emit `depth` synthetic 150 bp reads covering each anchor, half of them from
/// the reverse strand so the counter's strand handling is exercised.
fn simulate_reads(
    db: &AnchorDb,
    genome: &[u8],
    genome_id: u32,
    depth_at: impl Fn(u64) -> u32,
) -> Vec<Vec<u8>> {
    let mut reads = Vec::new();
    let mut rng = Rng(0xC0FFEE);
    for a in &db.anchors[db.genome_range(genome_id)] {
        let tag_len = a.tag_len() as u64;
        let d = depth_at(a.position);
        for _ in 0..d {
            // Place the read so it fully contains the tag, at a random offset.
            let slack = 150u64.saturating_sub(tag_len + 24);
            let jitter = if slack > 0 { rng.below(slack) } else { 0 };
            let start = a.position.saturating_sub(12 + jitter);
            let end = (start + 150).min(genome.len() as u64);
            if end <= start + tag_len + 24 {
                continue;
            }
            let read = genome[start as usize..end as usize].to_vec();
            reads.push(if rng.next_u64() % 2 == 0 {
                read
            } else {
                revcomp(&read)
            });
        }
    }
    reads
}

fn count_reads(
    index: &AnchorIndex<'_>,
    reads: &[Vec<u8>],
    cfg: &MatchConfig,
    n: usize,
) -> (Vec<u32>, CountStats) {
    let mut counts = vec![0u32; n];
    let mut stats = CountStats::default();
    for r in reads {
        stats.reads_total += 1;
        if index.count_read(r, cfg, &mut counts, &mut stats) > 0 {
            stats.reads_with_anchor += 1;
        }
    }
    (counts, stats)
}

#[test]
fn sixteen_enzymes_are_denser_than_one() {
    let dir = tmpdir("density");
    let g = random_genome(300_000, 42);
    let single = build_db(&dir, &[("g0", g.clone())], "BcgI");
    let union = build_db(&dir, &[("g0", g.clone())], "all");

    let per_mb = |db: &AnchorDb| db.n_anchors() as f64 / 0.3;
    // The report's headline density argument, on a synthetic genome: the union
    // is an order of magnitude denser than any single enzyme.
    assert!(
        per_mb(&union) > 10.0 * per_mb(&single),
        "union {:.0}/Mb vs BcgI {:.0}/Mb",
        per_mb(&union),
        per_mb(&single)
    );

    // And, more importantly, the worst-case blind spot shrinks by orders of
    // magnitude — the "single enzymes go blind locally" claim.
    let s_single = spacing_stats(&single, 0, 5_000, true);
    let s_union = spacing_stats(&union, 0, 5_000, true);
    // On a 300 kb synthetic genome the single-enzyme worst case cannot get as
    // extreme as it does on a real 4.6 Mb chromosome (the report measures 104 kb
    // for BplI on E. coli), so this is a floor, not the full effect size.
    assert!(
        s_union.max_gap * 5 < s_single.max_gap,
        "union max gap {} vs single {}",
        s_union.max_gap,
        s_single.max_gap
    );
    assert!(
        s_union.max_gap < 3_000,
        "union max gap {} exceeds the ~1.5 kb the report measures",
        s_union.max_gap
    );
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn counts_track_a_planted_replication_gradient() {
    let dir = tmpdir("gradient");
    let len = 300_000u64;
    let g = random_genome(len as usize, 7);
    let db = build_db(&dir, &[("g0", g.clone())], "all");
    assert!(db.n_anchors() > 1_000, "only {} anchors", db.n_anchors());

    // Depth falls linearly from position 0 (ori) to the midpoint (ter).
    let depth_at = |pos: u64| {
        let d = pos.min(len - pos) as f64 / (len as f64 / 2.0);
        (12.0 * (1.0 - 0.5 * d)).round() as u32
    };
    let reads = simulate_reads(&db, &g, 0, depth_at);
    let index = AnchorIndex::build(&db, 2);
    let (counts, stats) = count_reads(&index, &reads, &MatchConfig::default(), db.n_anchors());

    assert!(stats.tag_matched > 0);
    // Every tag that was fully contained in a read must find its anchor. The
    // denominator is extracted tags, not motif hits: with the 16-enzyme union a
    // 150 bp read spans several anchors and its edge ones are always truncated.
    assert!(
        stats.resolved_rate() > 0.98,
        "only {:.1}% of extracted tags resolved ({} unmatched of {} extracted)",
        100.0 * stats.resolved_rate(),
        stats.tag_unmatched,
        stats.tags_extracted()
    );
    assert_eq!(
        stats.mismatch_hist[1] + stats.mismatch_hist[2],
        0,
        "exact reads should match exactly"
    );

    // Anchors near the origin must be deeper than anchors near the terminus.
    let near: Vec<u32> = db
        .anchors
        .iter()
        .zip(&counts)
        .filter(|(a, _)| a.position < len / 10)
        .map(|(_, &c)| c)
        .collect();
    let far: Vec<u32> = db
        .anchors
        .iter()
        .zip(&counts)
        .filter(|(a, _)| (a.position as i64 - (len / 2) as i64).unsigned_abs() < len / 10)
        .map(|(_, &c)| c)
        .collect();
    let mean = |v: &[u32]| v.iter().map(|&c| c as f64).sum::<f64>() / v.len() as f64;
    assert!(!near.is_empty() && !far.is_empty());
    assert!(
        mean(&near) > 1.3 * mean(&far),
        "ori {:.2} vs ter {:.2}",
        mean(&near),
        mean(&far)
    );
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn shared_anchors_are_masked_and_reassigned_by_abundance() {
    let dir = tmpdir("shared");
    // Two genomes sharing a 30 kb block: those anchors are cross-genome shared.
    let base = random_genome(200_000, 11);
    let mut g1 = random_genome(200_000, 12);
    g1[50_000..80_000].copy_from_slice(&base[50_000..80_000]);
    let db = build_db(&dir, &[("g0", base.clone()), ("g1", g1.clone())], "all");

    let shared = db
        .anchors
        .iter()
        .filter(|a| a.flags & flags::MASKED_SHARED != 0)
        .count();
    assert!(
        shared > 50,
        "expected a shared block, found {shared} shared anchors"
    );
    assert!(
        db.anchors.iter().any(|a| a.is_usable()),
        "everything was masked"
    );

    // g0 is 5x more abundant, so it should absorb ~5/6 of each shared count.
    let counts: Vec<u32> = db
        .anchors
        .iter()
        .map(|a| if a.genome_id == 0 { 10 } else { 2 })
        .collect();
    let em = reassign(&db, &counts, &EmConfig::default());
    assert!(em.converged);

    let groups = db.shared_groups();
    let (_, idxs) = groups.iter().next().expect("no shared group");
    let total: f64 = idxs.iter().map(|&i| em.weights[i]).sum();
    let g0_share: f64 = idxs
        .iter()
        .filter(|&&i| db.anchors[i].genome_id == 0)
        .map(|&i| em.weights[i])
        .sum();
    assert!(total > 0.0);
    assert!(
        (g0_share / total - 5.0 / 6.0).abs() < 0.05,
        "g0 took {:.3} of the shared mass",
        g0_share / total
    );
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn database_survives_a_save_load_roundtrip() {
    let dir = tmpdir("roundtrip");
    let db = build_db(&dir, &[("g0", random_genome(120_000, 3))], "BcgI,AlfI,CjeI");
    let dbdir = dir.join("db");
    db.save(&dbdir).unwrap();
    let back = AnchorDb::load(&dbdir).unwrap();
    assert_eq!(back, db);
    assert_eq!(back.params.enzymes.len(), 3);
    // Tag sequences must survive packing, not just the metadata.
    for i in (0..db.n_anchors()).step_by(37) {
        assert_eq!(back.tag(i), db.tag(i));
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn equal_anchor_windows_beat_fixed_bp_on_occupancy_variance() {
    let dir = tmpdir("windows");
    let db = build_db(&dir, &[("g0", random_genome(400_000, 5))], "all");
    let (eq, _) = assign_windows(&db, 0, WindowPolicy::EqualAnchors { n: 100 }, true);
    let (bp, _) = assign_windows(&db, 0, WindowPolicy::FixedBp { bp: 25_000 }, true);

    let cv = |w: &[sk2bgrow_core::window::Window]| {
        let n: Vec<f64> = w.iter().map(|x| x.n_anchors as f64).collect();
        let m = n.iter().sum::<f64>() / n.len() as f64;
        (n.iter().map(|x| (x - m).powi(2)).sum::<f64>() / n.len() as f64).sqrt() / m
    };
    assert!(!eq.is_empty() && !bp.is_empty());
    // The whole point of the adaptive window: equal statistical power per window.
    assert!(
        cv(&eq) < cv(&bp),
        "equal-anchor CV {:.3} vs fixed-bp CV {:.3}",
        cv(&eq),
        cv(&bp)
    );
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn two_brad_mode_counts_one_tag_per_read() {
    let dir = tmpdir("2brad");
    let g = random_genome(200_000, 9);
    let db = build_db(&dir, &[("g0", g.clone())], "BcgI");
    let index = AnchorIndex::build(&db, 2);

    // Route B reads are the tags themselves.
    let reads: Vec<Vec<u8>> = db.anchors[db.genome_range(0)]
        .iter()
        .take(50)
        .filter_map(|a| {
            let e = sk2bgrow_core::enzyme::by_idx(a.enzyme_idx)?;
            let (s, t) = e.tag_span(a.position as usize, g.len(), a.strand == 0)?;
            Some(if a.strand == 0 {
                g[s..t].to_vec()
            } else {
                revcomp(&g[s..t])
            })
        })
        .collect();
    assert_eq!(reads.len(), 50);

    let cfg = MatchConfig {
        mode: CountMode::TwoBrad,
        ..MatchConfig::default()
    };
    let (counts, stats) = count_reads(&index, &reads, &cfg, db.n_anchors());
    assert_eq!(
        stats.reads_with_anchor, 50,
        "a bare tag should always match its own anchor"
    );
    assert_eq!(counts.iter().sum::<u32>(), 50);
    std::fs::remove_dir_all(dir).ok();
}
