//! Local benchmark for the M4 screen: full-index vs genome-screened counting,
//! and serial vs per-file parallel counting, on a database large enough that
//! the lookup CSR does not fit in cache — the regime the HPC benchmark hit.
//!
//! The database is 2 "present" genomes (3 Mb each, real digestion) padded
//! with fake genomes carrying random unique tags, so the exact/seed CSRs
//! span millions of keys like the real 24.1M-anchor database. Reads are
//! sampled from the present genomes with ~3% divergence (so most tag hits
//! travel the mismatch-tolerant seed path, as with real metagenomes) plus
//! foreign reads.
//!
//! ```sh
//! cargo run --release --example screen_bench -- serial    # main-style serial loop
//! cargo run --release --example screen_bench -- parallel  # count_sample over 2 files
//! cargo run --release --example screen_bench -- screen    # pass 2 with a screened index
//! ```
//!
//! RAYON_NUM_THREADS controls the pool for the parallel/screen modes.

use std::path::PathBuf;

use sk2bgrow_core::anchor_db::{assemble, build_genome, Anchor, BuildParams, GenomeMeta, GC_FLANK};
use sk2bgrow_core::count::{count_sample, AnchorIndex, MatchConfig};
use sk2bgrow_core::digest::DigestConfig;
use sk2bgrow_core::enzyme::{EnzymeSet, PANEL};
use sk2bgrow_core::tgt::{pack_bases, ContigKind, ContigMeta};

/// splitmix64-style LCG; deterministic so the benchmark is reproducible.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() & 3) as usize]
    }
}

const N_PRESENT: usize = 2;
const N_FAKE_GENOMES: usize = 150;
const ANCHORS_PER_FAKE: usize = 40_000;
const N_READS: usize = 1_000_000;

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "parallel".into());
    let n_fake = env_usize("BENCH_FAKE_GENOMES", N_FAKE_GENOMES);
    let n_reads = env_usize("BENCH_READS", N_READS);
    let dir = std::env::temp_dir().join(format!("sk2bgrow-screen-bench-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    let mut rng = Rng(0x5eed);
    // Present genomes: real digestion.
    let enzymes: Vec<&'static _> = PANEL.iter().collect();
    let mut parts = Vec::new();
    let mut present_seqs = Vec::new();
    for gid in 0..N_PRESENT {
        let seq: Vec<u8> = (0..3_000_000).map(|_| rng.base()).collect();
        let p = dir.join(format!("present{gid}.fna"));
        std::fs::write(&p, format!(">c0\n{}\n", String::from_utf8_lossy(&seq))).unwrap();
        parts.push(build_genome(&p, gid as u32, &enzymes, &DigestConfig::default(), GC_FLANK).unwrap());
        present_seqs.push(seq);
        std::fs::remove_file(&p).ok();
    }
    // Fake genomes: unique random hashes/tags so lookups miss and probe the
    // seed tables, modelling the hundreds of genomes a real sample lacks.
    let t_build = std::time::Instant::now();
    let mut fake_anchors = Vec::with_capacity(n_fake * ANCHORS_PER_FAKE);
    let mut fake_tags = Vec::with_capacity(n_fake * ANCHORS_PER_FAKE);
    for gid in N_PRESENT..N_PRESENT + n_fake {
        for i in 0..ANCHORS_PER_FAKE {
            let h = rng.next() | 1; // odd: cannot collide with a packed 2-bit tag hash
            let tag: Vec<u8> = (0..32).map(|_| rng.base()).collect();
            fake_anchors.push(Anchor {
                seq_hash: h,
                genome_id: gid as u32,
                contig_id: 0,
                position: i as u64,
                enzyme_idx: enzymes[0].idx,
                strand: 0,
                flags: 0,
                local_gc: 100,
            });
            fake_tags.push(pack_bases(&tag));
        }
    }
    let mut metas: Vec<GenomeMeta> = (0..N_PRESENT)
        .map(|i| parts[i].0.clone())
        .collect();
    metas.extend((N_PRESENT..N_PRESENT + n_fake).map(|gid| GenomeMeta {
        id: gid as u32,
        name: format!("fake{gid}"),
        taxonomy: None,
        contigs: vec![ContigMeta {
            id: 0,
            name: "c0".into(),
            length: 3_000_000,
            offset: 0,
            kind: ContigKind::Chromosome,
        }],
        genome_len: 3_000_000,
        ori: None,
        ori_confidence: 0.0,
    }));
    let mut db = assemble(
        BuildParams {
            enzymes: EnzymeSet::from_slice(&enzymes),
            ..BuildParams::default()
        },
        parts
            .into_iter()
            .map(|(m, a, t, _)| (m, a, t))
            .collect(),
    );
    db.genomes = metas;
    db.anchors.extend(fake_anchors);
    db.tags.extend(fake_tags);
    db.recompute_uniqueness();
    eprintln!(
        "db: {} anchors, {} genomes (assembled in {:.1}s)",
        db.n_anchors(),
        db.genomes.len(),
        t_build.elapsed().as_secs_f64()
    );

    let t = std::time::Instant::now();
    let index = AnchorIndex::build(&db, 2);
    eprintln!("full index built in {:.1}s", t.elapsed().as_secs_f64());

    // Reads: 70% from present genomes at 3% divergence, 30% foreign.
    let foreign: Vec<u8> = (0..3_000_000).map(|_| rng.base()).collect();
    let mut files: Vec<PathBuf> = Vec::new();
    for f in 0..2 {
        let mut body = String::with_capacity(n_reads / 2 * 160);
        for r in 0..n_reads / 2 {
            if (f * (n_reads / 2) + r) % 10 < 7 {
                let src = &present_seqs[(f * (n_reads / 2) + r) % N_PRESENT];
                let start = (rng.next() as usize) % (src.len() - 150);
                let mut read: Vec<u8> = src[start..start + 150].to_vec();
                for b in read.iter_mut() {
                    if rng.next() % 100 < 3 {
                        *b = b"ACGT"[((rng.next() & 3) as usize + 1) % 4];
                    }
                }
                body.push_str(&format!(">r{f}_{r}\n{}\n", String::from_utf8_lossy(&read)));
            } else {
                let start = (rng.next() as usize) % (foreign.len() - 150);
                body.push_str(&format!(
                    ">r{f}_{r}\n{}\n",
                    String::from_utf8_lossy(&foreign[start..start + 150])
                ));
            }
        }
        let p = dir.join(format!("reads{f}.fna"));
        std::fs::write(&p, body).unwrap();
        files.push(p);
    }

    let cfg = MatchConfig::default();
    match mode.as_str() {
        // main-branch serial loop, for a like-for-like baseline.
        "serial" => {
            let t = std::time::Instant::now();
            let mut counts = vec![0u32; db.n_anchors()];
            let mut stats = sk2bgrow_core::count::CountStats::default();
            for path in &files {
                sk2bgrow_core::fasta::for_each_read(path, |read| {
                    stats.reads_total += 1;
                    if index.count_read(read, &cfg, &mut counts, &mut stats) > 0 {
                        stats.reads_with_anchor += 1;
                    }
                })
                .unwrap();
            }
            eprintln!(
                "[bench] serial: {:.2}s, {} reads, {} matched tags",
                t.elapsed().as_secs_f64(),
                stats.reads_total,
                stats.tag_matched
            );
        }
        "parallel" => {
            let t = std::time::Instant::now();
            let (_, stats) = count_sample(&index, &files, &cfg).unwrap();
            eprintln!(
                "[bench] parallel: {:.2}s, {} reads, {} matched tags",
                t.elapsed().as_secs_f64(),
                stats.reads_total,
                stats.tag_matched
            );
        }
        "screen" => {
            let mut filter = vec![false; db.genomes.len()];
            for g in 0..N_PRESENT {
                filter[g] = true;
            }
            let t = std::time::Instant::now();
            let screened = AnchorIndex::build_screened(&db, 2, None, Some(&filter));
            let t_build = t.elapsed();
            let t = std::time::Instant::now();
            let (_, stats) = count_sample(&screened, &files, &cfg).unwrap();
            eprintln!(
                "[bench] screen pass2: {:.2}s (index build {:.2}s), {} reads, {} matched tags",
                t.elapsed().as_secs_f64(),
                t_build.as_secs_f64(),
                stats.reads_total,
                stats.tag_matched
            );
        }
        other => eprintln!("unknown mode {other}"),
    }
    std::fs::remove_dir_all(&dir).ok();
}
