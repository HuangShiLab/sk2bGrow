//! Local benchmark: motif scan vs hash lookup in the count loop.
//!
//! Reproduces the M4 phase-timing question on synthetic data: a ~3 Mb random
//! genome digested with the full 16-enzyme panel, then 100k reads (150 bp for
//! `wms`, 33 bp for `2brad`), half sampled from the indexed genome and half
//! from unrelated random sequence so both the matched and the pure-scan cases
//! are represented. Run with:
//!
//! ```sh
//! cargo run --release --example phase_timing
//! ```
//!
//! The example sets `SK2B_COUNT_TIMING=1` itself; the counters and the split
//! are printed on stderr by `count_sample`.
//!
//! Measured on 2026-09-03 (Apple Silicon, release build): see the git history
//! of this file's commit message for the numbers.

use sk2bgrow_core::anchor_db::{assemble, build_genome, BuildParams, GC_FLANK};
use sk2bgrow_core::count::{count_sample, AnchorIndex, CountMode, MatchConfig};
use sk2bgrow_core::digest::DigestConfig;
use sk2bgrow_core::enzyme::{EnzymeSet, PANEL};

/// splitmix64-style LCG; deterministic so the benchmark is reproducible.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() & 3) as usize]
    }
}

fn main() {
    std::env::set_var("SK2B_COUNT_TIMING", "1");
    let dir = std::env::temp_dir().join(format!("sk2bgrow-phase-timing-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    // One ~3 Mb random contig as the reference.
    let mut rng = Rng(0x5eed);
    let genome: Vec<u8> = (0..3_000_000).map(|_| rng.base()).collect();
    let gpath = dir.join("genome.fna");
    std::fs::write(
        &gpath,
        format!(">c0\n{}\n", String::from_utf8_lossy(&genome)),
    )
    .unwrap();

    // An unrelated random sequence: reads from it exercise the pure scan case
    // (every window misses every pattern) that dominates real metagenomes.
    let foreign: Vec<u8> = (0..3_000_000).map(|_| rng.base()).collect();

    let enzymes: Vec<&'static _> = PANEL.iter().collect();
    let (meta, anchors, tags, _) =
        build_genome(&gpath, 0, &enzymes, &DigestConfig::default(), GC_FLANK).unwrap();
    let db = assemble(
        BuildParams {
            enzymes: EnzymeSet::from_slice(&enzymes),
            ..BuildParams::default()
        },
        vec![(meta, anchors, tags)],
    );
    eprintln!("db: {} anchors on one 3 Mb genome", db.n_anchors());

    let index = AnchorIndex::build(&db, 2);

    for (label, mode, len) in [("wms/150bp", CountMode::Wms, 150), ("2brad/33bp", CountMode::TwoBrad, 33)] {
        let rpath = dir.join(format!("reads-{len}.fna"));
        // Half the reads come from the indexed genome, half from foreign DNA.
        let mut file_rng = Rng(0xbeef + len as u64);
        let mut genome_reads = String::new();
        let mut foreign_reads = String::new();
        for (i, src, buf) in [
            (0..50_000, &genome, &mut genome_reads),
            (50_000..100_000, &foreign, &mut foreign_reads),
        ] {
            for j in i {
                let start = (file_rng.next() as usize) % (src.len() - len);
                buf.push_str(&format!(
                    ">r{j}\n{}\n",
                    String::from_utf8_lossy(&src[start..start + len])
                ));
            }
        }
        std::fs::write(&rpath, genome_reads + &foreign_reads).unwrap();

        let cfg = MatchConfig {
            max_mismatch: 2,
            mode,
            keep_multimappers: true,
        };
        eprintln!("--- {label}: 100k reads (50k from the genome, 50k foreign) ---");
        let t = std::time::Instant::now();
        let (counts, stats) = count_sample(&index, &[rpath.clone()], &cfg).unwrap();
        eprintln!(
            "[phase-timing] {label}: count_sample wall {:.3} s, {} reads, {} matched tags, counts {}",
            t.elapsed().as_secs_f64(),
            stats.reads_total,
            stats.tag_matched,
            counts.iter().map(|&c| c as u64).sum::<u64>(),
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}
