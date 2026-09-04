//! Genome-level containment pre-screen (M4 tier 1).
//!
//! Before the anchor-level count (tier 2) decides *which loci* are present,
//! a FracMinHash-style containment screen decides *which genomes* are worth
//! counting at all. At GTDB scale (137 k representatives) the full
//! [`crate::count::AnchorIndex`] over every anchor is what will not fit in
//! memory; the screen reduces pass 2 to the genomes a sample actually
//! contains, usually a tiny fraction of the database.
//!
//! ## Sketch
//!
//! At index time (`sk2bgrow index --screen-scale S`) each reference genome is
//! sketched once: every canonical k-mer hash with `h < u64::MAX / S` is kept
//! (FracMinHash sampling; k defaults to 21). A global CSR maps kept hash ->
//! genome ids, so one lookup per read k-mer votes for every genome sharing
//! that k-mer. The sketch lives next to the anchor database as
//! `screen.meta` + `screen.csr`; databases built without it simply lack the
//! files, `load` returns `None`, and `profile --screen` says so instead of
//! guessing. The anchor files themselves are untouched, so old binaries read
//! new databases and vice versa.
//!
//! ## Correctness of the two-tier count
//!
//! Pass 2 builds a [`crate::count::AnchorIndex::build_screened`] restricted to
//! the selected genomes, which behaves exactly like a database built from
//! only those genomes:
//!
//! * An anchor whose tag is unique to a selected genome resolves to the same
//!   CSR entries as in the full index, so its counts are identical to an
//!   unscreened run.
//! * A tag shared between a selected and an unselected genome is credited
//!   only at the selected locus — the same semantics `--enzymes` restriction
//!   already has (a restricted index behaves like a database built with
//!   fewer enzymes). The EM then sees the counts a reduced database would
//!   have produced.
//! * Unselected genomes' anchors are absent from the lookup tables, so they
//!   are never credited; their count-table rows are zeros, which downstream
//!   layers already handle (`--enzymes` produces the same shape).
//!
//! ## What pass 2 costs — and why a screened pass 2 is fast
//!
//! Pass 2 scans **every read, every window** — it does not skip reads that
//! pass 1 missed. The motif scan is index-size-independent, so the screen's
//! speedup comes from the lookup step: on the full database the exact CSR
//! spans millions of keys (~200 MB at 24.1M anchors), every binary-search
//! probe is a cache miss, and a tag that does not exact-match — the common
//! case for reads from diverged organisms — then probes both orientations ×
//! (m+1) seed tables before per-candidate verification. On a 24.1M-anchor
//! database that lookup dominates the count (measured: 2M diverged reads,
//! 22M anchors, ~25 µs/read total with the scan only ~5 µs of it). The
//! screened index holds only the selected genomes' anchors and is
//! cache-resident, so the same lookups cost almost nothing; what remains is
//! the motif scan, which is why a screened pass 2 runs in seconds while the
//! full count runs in minutes.
//!
//! Because the genome subset shrinks the best-hit set, screened counts can
//! differ slightly from a full run even on selected genomes: a tag matching
//! genome A at distance 1 and genome B (unselected) at distance 0 is credited
//! to B alone in the full run and to A in the screened run. Measured on real
//! data: 296 of 1.43M anchors, all with counts ≤ 3.
//!
//! ## A tempting wrong optimisation: skipping pass-2 reads
//!
//! The naive argument — "a tag matching within m mismatches shares a
//! k = 21 mer with the reference, so reads with no sketch hit cannot
//! match" — is **false** for m = 2: two mismatches placed ~k apart inside a
//! 32 bp tag break every 21-mer of it (mismatches at offsets 11 and 21 break
//! all twelve windows), and the pigeonhole guarantee for (len, m) = (32, 2)
//! is only an exact run of ⌈(32−2)/3⌉ = 10 bp. A strictly sound read-skip
//! filter would need the sketch k at the matcher's seed length (~10 bp),
//! which is far too short for containment specificity. Per-genome skipping
//! at k = 21 is sound only for exact and 1-mismatch matches. The current
//! implementation therefore never skips reads; the screen filters *which
//! genomes enter the index*, not which reads are scanned. Count correctness
//! never depends on the sketch — only speed does — and the sketch threshold
//! plus the 50-hit floor bound how much signal a false-negative genome can
//! lose.

use std::path::{Path, PathBuf};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::anchor_db::AnchorDb;
use crate::error::{Result, Sk2bError};
use crate::fasta::for_each_read;
use crate::seq::two_bit;

/// Default k-mer length for the containment sketch. 21 is the sourmash/
/// mash default and leaves 13 usable k-mers on a 33 bp 2bRAD read; reads
/// shorter than k simply contribute no k-mers (the roller never fills).
pub const SCREEN_K: usize = 21;

// hash_tag's finaliser, factored out so the rolling sketch is bit-identical
// to hashing each window with seq::canonical_hash.
#[inline]
fn mix(packed: u64, k: usize) -> u64 {
    let mut z = packed
        .wrapping_add((k as u64) << 58)
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Rolling strand-canonical k-mer hash, O(1) per base and bit-identical to
/// `seq::canonical_hash(&window)` on every window. Maintains the 2-bit packed
/// forward window and its reverse complement simultaneously; the canonical
/// hash is the smaller of the two finalised values. Any ambiguous base resets
/// the window, so k-mers spanning an `N` are never emitted — the same set a
/// per-window `canonical_hash` over ACGTN-normalised sequence would produce.
struct KmerRoller {
    k: usize,
    mask: u64,
    /// 2-bit packed forward window (first base in the most significant bits,
    /// the same convention pack_2bit uses).
    fwd: u64,
    /// 2-bit packed reverse complement of the window.
    rev: u64,
    /// Consecutive ACGT bases seen since the last ambiguous base.
    valid: usize,
}

impl KmerRoller {
    fn new(k: usize) -> Self {
        assert!((1..=32).contains(&k), "roller supports k in 1..=32");
        KmerRoller {
            k,
            mask: (1u64 << (2 * k)) - 1,
            fwd: 0,
            rev: 0,
            valid: 0,
        }
    }

    /// Feed one base; returns the canonical hash once k valid bases fill the
    /// window, `None` otherwise.
    #[inline]
    fn feed(&mut self, b: u8) -> Option<u64> {
        let code = match two_bit(b) {
            Some(c) => c as u64,
            None => {
                self.valid = 0;
                return None;
            }
        };
        self.fwd = ((self.fwd << 2) | code) & self.mask;
        // complement in 2-bit code space is 3 - c (A<->T, C<->G); the new
        // base's complement enters the revcomp window at its left end.
        self.rev = (self.rev >> 2) | ((3 - code) << (2 * (self.k - 1)));
        self.valid += 1;
        if self.valid < self.k {
            return None;
        }
        Some(mix(self.fwd, self.k).min(mix(self.rev, self.k)))
    }
}

/// Per-genome containment sketches plus the shared hash -> genomes CSR.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScreenSketch {
    pub k: usize,
    /// FracMinHash scale: a k-mer hash is kept when `h < u64::MAX / scale`.
    pub scale: u64,
    /// Kept-hash count per genome, indexed by genome id.
    pub sizes: Vec<u32>,
    /// Sorted kept hashes.
    keys: Vec<u64>,
    /// `offsets[j]..offsets[j + 1]` of `vals` is `keys[j]`'s genome list.
    offsets: Vec<u32>,
    /// Genome ids per hash, sorted and deduplicated.
    vals: Vec<u32>,
}

/// JSON sidecar, the same manifest-plus-blob split the anchor database uses.
#[derive(Serialize, Deserialize)]
struct ScreenMeta {
    format_version: u32,
    k: usize,
    scale: u64,
    sizes: Vec<u32>,
}

impl ScreenSketch {
    /// Sketch every genome. `genome_paths[i]` becomes genome id `i`, which
    /// must match the id assignment of the anchor database built from the
    /// same paths (index assigns ids from sorted input order).
    ///
    /// Note the canonical hash is the *minimum* of two strand hashes, so the
    /// keep probability is ~2/scale, not 1/scale — the effective sampling is
    /// twice as dense as the flag suggests. Build and query share the same
    /// rule, so the sketch stays self-consistent; only the sizing expectation
    /// changes.
    pub fn build(genome_paths: &[PathBuf], k: usize, scale: u64) -> Result<Self> {
        if scale < 2 {
            return Err(Sk2bError::Config(format!(
                "screen scale must be >= 2, got {scale}"
            )));
        }
        let threshold = u64::MAX / scale;
        // Per-genome kept hashes, deduplicated. A genome is digested once
        // for anchors and once for the sketch; sharing the FASTA read would
        // couple the two builders for little gain (indexing is offline).
        let per_genome: Vec<Vec<u64>> = genome_paths
            .par_iter()
            .map(|path| {
                let mut kept: Vec<u64> = Vec::new();
                for_each_read(path, |seq| {
                    let mut roller = KmerRoller::new(k);
                    for &b in seq {
                        if let Some(h) = roller.feed(b) {
                            if h < threshold {
                                kept.push(h);
                            }
                        }
                    }
                })
                .map_err(|e| Sk2bError::Db(format!("sketching {}: {e}", path.display())))?;
                kept.sort_unstable();
                kept.dedup();
                Ok(kept)
            })
            .collect::<Result<Vec<_>>>()?;
        let sizes: Vec<u32> = per_genome
            .iter()
            .map(|v| v.len() as u32)
            .collect();

        let mut pairs: Vec<(u64, u32)> = Vec::new();
        for (gid, hashes) in per_genome.iter().enumerate() {
            pairs.reserve(hashes.len());
            for &h in hashes {
                pairs.push((h, gid as u32));
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        let mut sketch = ScreenSketch {
            k,
            scale,
            sizes,
            keys: Vec::with_capacity(pairs.len()),
            offsets: vec![0],
            vals: Vec::with_capacity(pairs.len()),
        };
        for &(h, g) in &pairs {
            if sketch.keys.last() != Some(&h) {
                sketch.keys.push(h);
                sketch.offsets.push(sketch.vals.len() as u32);
            }
            sketch.vals.push(g);
            *sketch.offsets.last_mut().unwrap() = sketch.vals.len() as u32;
        }
        Ok(sketch)
    }

    pub fn save(&self, dir: &Path) -> Result<()> {
        let meta = ScreenMeta {
            format_version: 1,
            k: self.k,
            scale: self.scale,
            sizes: self.sizes.clone(),
        };
        let mpath = dir.join("screen.meta");
        std::fs::write(&mpath, serde_json::to_vec_pretty(&meta)?).map_err(|e| {
            Sk2bError::Io {
                path: mpath,
                source: e,
            }
        })?;
        let cpath = dir.join("screen.csr");
        let f = std::fs::File::create(&cpath).map_err(|e| Sk2bError::Io {
            path: cpath.clone(),
            source: e,
        })?;
        let mut w = std::io::BufWriter::new(f);
        bincode::serialize_into(&mut w, &(&self.keys, &self.offsets, &self.vals))?;
        use std::io::Write;
        w.flush().map_err(|e| Sk2bError::Io {
            path: cpath,
            source: e,
        })?;
        Ok(())
    }

    /// Load the sketch next to an anchor database. `Ok(None)` when the
    /// database was built without `--screen-scale` — the caller turns that
    /// into a "rebuild the index" error rather than silently profiling
    /// without a screen.
    pub fn load(dir: &Path) -> Result<Option<Self>> {
        let mpath = dir.join("screen.meta");
        if !mpath.exists() {
            return Ok(None);
        }
        let bytes = std::fs::read(&mpath).map_err(|e| Sk2bError::Io {
            path: mpath,
            source: e,
        })?;
        let meta: ScreenMeta = serde_json::from_slice(&bytes)?;
        if meta.format_version != 1 {
            return Err(Sk2bError::Db(format!(
                "screen sketch format version {} is not supported (expected 1)",
                meta.format_version
            )));
        }
        let cpath = dir.join("screen.csr");
        let f = std::fs::File::open(&cpath).map_err(|e| Sk2bError::Io {
            path: cpath.clone(),
            source: e,
        })?;
        let (keys, offsets, vals): (Vec<u64>, Vec<u32>, Vec<u32>) =
            bincode::deserialize_from(std::io::BufReader::new(f))?;
        if keys.len() + 1 != offsets.len() {
            return Err(Sk2bError::Db(
                "corrupt screen sketch: offsets do not match keys".into(),
            ));
        }
        Ok(Some(ScreenSketch {
            k: meta.k,
            scale: meta.scale,
            sizes: meta.sizes,
            keys,
            offsets,
            vals,
        }))
    }

    /// Genomes whose sketch kept hash `h` (empty when `h` was not kept by
    /// any genome).
    #[inline]
    pub fn genomes_for(&self, h: u64) -> &[u32] {
        match self.keys.binary_search(&h) {
            Ok(j) => &self.vals[self.offsets[j] as usize..self.offsets[j + 1] as usize],
            Err(_) => &[],
        }
    }

    /// Select the genomes pass 2 should count.
    ///
    /// `hits_g` is the number of reads carrying at least one k-mer kept by
    /// genome g (see [`screen_sample`]). Abundance is estimated as
    /// `m_g = hits_g / sketch_size_g * genome_len_g` — the kept fraction of
    /// the genome's sketch that was observed, scaled back to genome size —
    /// and normalised against the total. A genome is kept when its estimated
    /// fraction reaches `min_frac`, or unconditionally when it has at least
    /// 50 hits: rare genomes below the relative threshold can still carry
    /// usable signal, and the cost of keeping one extra genome in pass 2 is
    /// small next to the cost of missing it.
    pub fn select(&self, db: &AnchorDb, hits: &[u32], min_frac: f64) -> Vec<u32> {
        let mut mass = vec![0.0f64; db.genomes.len()];
        let mut total = 0.0f64;
        for g in &db.genomes {
            let h = hits.get(g.id as usize).copied().unwrap_or(0);
            let s = self.sizes.get(g.id as usize).copied().unwrap_or(0);
            if h == 0 || s == 0 {
                continue;
            }
            let m = h as f64 / s as f64 * g.genome_len as f64;
            mass[g.id as usize] = m;
            total += m;
        }
        db.genomes
            .iter()
            .filter(|g| {
                let h = hits.get(g.id as usize).copied().unwrap_or(0);
                h >= 50 || (total > 0.0 && mass[g.id as usize] / total >= min_frac)
            })
            .map(|g| g.id)
            .collect()
    }
}

/// Pass 1: for every read, vote once per genome sharing any kept k-mer with
/// it. Returns per-genome read-hit counts.
///
/// Files are processed in parallel like [`crate::count::count_sample`]; the
/// per-file partials are summed, which is order-independent.
pub fn screen_sample(sketch: &ScreenSketch, reads: &[PathBuf]) -> Result<Vec<u32>> {
    let threshold = u64::MAX / sketch.scale;
    let parts: Vec<Result<Vec<u32>>> = reads
        .par_iter()
        .map(|path| {
            let mut hits = vec![0u32; sketch.sizes.len()];
            // One vote per genome per read: an epoch stamp per genome avoids
            // both double-counting a genome on repeated k-mers and a hash
            // set per read.
            let mut stamp = vec![u32::MAX; sketch.sizes.len()];
            let mut epoch = 0u32;
            for_each_read(path, |read| {
                epoch = epoch.wrapping_add(1);
                if epoch == u32::MAX {
                    stamp.iter_mut().for_each(|s| *s = u32::MAX);
                    epoch = 0;
                }
                let mut roller = KmerRoller::new(sketch.k);
                for &b in read {
                    if let Some(h) = roller.feed(b) {
                        if h >= threshold {
                            continue;
                        }
                        for &g in sketch.genomes_for(h) {
                            let g = g as usize;
                            if stamp[g] != epoch {
                                stamp[g] = epoch;
                                hits[g] += 1;
                            }
                        }
                    }
                }
            })?;
            Ok(hits)
        })
        .collect();
    let mut hits = vec![0u32; sketch.sizes.len()];
    for part in parts {
        for (dst, src) in hits.iter_mut().zip(part?.iter()) {
            *dst += src;
        }
    }
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::canonical_hash;

    /// splitmix64-style LCG; deterministic test sequence source.
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

    fn write_fasta(name: &str, seqs: &[Vec<u8>]) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-{name}", std::process::id()));
        let mut body = String::new();
        for (i, s) in seqs.iter().enumerate() {
            body.push_str(&format!(">c{i}\n{}\n", String::from_utf8_lossy(s)));
        }
        std::fs::write(&p, body).unwrap();
        p
    }

    #[test]
    fn roller_matches_window_by_window_hashing() {
        // The pass-1 lookup and the index-time sketch must agree bit for bit,
        // or genomes would vote for hashes they never kept.
        let mut rng = Rng(42);
        for case in 0..20 {
            let len = 30 + (rng.next() as usize) % 40;
            let seq: Vec<u8> = (0..len).map(|_| rng.base()).collect();
            let k = 21usize;
            let mut roller = KmerRoller::new(k);
            let mut emitted = Vec::new();
            for &b in &seq {
                if let Some(h) = roller.feed(b) {
                    emitted.push(h);
                }
            }
            let want: Vec<u64> = (0..=(seq.len() - k)).map(|i| canonical_hash(&seq[i..i + k])).collect();
            assert_eq!(emitted, want, "case {case}: roller diverged");
        }
        // Ambiguous bases break the window: nothing spanning an N is emitted.
        let mut seq: Vec<u8> = (0..50).map(|_| rng.base()).collect();
        seq[30] = b'N';
        let mut roller = KmerRoller::new(21);
        let mut emitted = Vec::new();
        for &b in &seq {
            if let Some(h) = roller.feed(b) {
                emitted.push(h);
            }
        }
        let want: Vec<u64> = (0..=seq.len() - 21)
            .filter(|&i| seq[i..i + 21].iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')))
            .map(|i| canonical_hash(&seq[i..i + 21]))
            .collect();
        assert_eq!(emitted, want, "N-containing windows leaked through");
    }

    /// Two genomes, 200 kb each: reads come from genome 0, so the screen
    /// must select exactly {0}.
    #[test]
    fn screen_selects_the_present_genome() {
        let mut rng = Rng(7);
        let g0: Vec<u8> = (0..200_000).map(|_| rng.base()).collect();
        let g1: Vec<u8> = (0..200_000).map(|_| rng.base()).collect();
        let p0 = write_fasta("screen-g0.fna", &[g0.clone()]);
        let p1 = write_fasta("screen-g1.fna", &[g1.clone()]);

        let sketch = ScreenSketch::build(&[p0.clone(), p1.clone()], SCREEN_K, 100).unwrap();
        assert_eq!(sketch.sizes.len(), 2);
        assert!(sketch.sizes[0] > 0 && sketch.sizes[1] > 0);

        // 5k reads sampled from genome 0.
        let mut reads = String::new();
        let mut rr = Rng(99);
        for i in 0..5000 {
            let start = (rr.next() as usize) % (g0.len() - 150);
            reads.push_str(&format!(
                ">r{i}\n{}\n",
                String::from_utf8_lossy(&g0[start..start + 150])
            ));
        }
        let rpath = {
            let mut p = std::env::temp_dir();
            p.push(format!("sk2bgrow-{}-screen-reads.fna", std::process::id()));
            std::fs::write(&p, reads).unwrap();
            p
        };
        let hits = screen_sample(&sketch, std::slice::from_ref(&rpath)).unwrap();
        assert!(hits[0] > 4000, "genome 0 reads barely hit: {}", hits[0]);
        assert_eq!(hits[1], 0, "foreign genome received hits");

        // A database stub with matching genome ids/lengths for `select`.
        let db = crate::anchor_db::AnchorDb {
            params: crate::anchor_db::BuildParams::default(),
            genomes: [0u32, 1]
                .iter()
                .map(|&id| crate::anchor_db::GenomeMeta {
                    id,
                    name: format!("g{id}"),
                    taxonomy: None,
                    contigs: vec![],
                    genome_len: 200_000,
                    ori: None,
                    ori_confidence: 0.0,
                })
                .collect(),
            anchors: vec![],
            tags: vec![],
        };
        let sel = sketch.select(&db, &hits, 5e-4);
        assert_eq!(sel, vec![0], "screen kept the wrong genomes: {sel:?}");

        std::fs::remove_file(p0).ok();
        std::fs::remove_file(p1).ok();
        std::fs::remove_file(rpath).ok();
    }

    /// A genome with 50 hits must survive even at a relative-abundance
    /// threshold it would otherwise fall under.
    #[test]
    fn fifty_hits_keep_a_genome_unconditionally() {
        let sketch = ScreenSketch {
            k: SCREEN_K,
            scale: 100,
            sizes: vec![10_000, 10_000],
            keys: vec![],
            offsets: vec![0],
            vals: vec![],
        };
        let db = crate::anchor_db::AnchorDb {
            params: crate::anchor_db::BuildParams::default(),
            genomes: [0u32, 1]
                .iter()
                .map(|&id| crate::anchor_db::GenomeMeta {
                    id,
                    name: format!("g{id}"),
                    taxonomy: None,
                    contigs: vec![],
                    genome_len: 200_000,
                    ori: None,
                    ori_confidence: 0.0,
                })
                .collect(),
            anchors: vec![],
            tags: vec![],
        };
        // Genome 0 dominates; genome 1 has 49 hits (dropped) ...
        let hits = vec![100_000, 49];
        assert_eq!(sketch.select(&db, &hits, 5e-4), vec![0]);
        // ... and 50 hits (kept, fraction ~1e-5 << 5e-4).
        let hits = vec![100_000, 50];
        assert_eq!(sketch.select(&db, &hits, 5e-4), vec![0, 1]);
    }

    #[test]
    fn sketch_roundtrips_and_absent_files_load_as_none() {
        let mut rng = Rng(3);
        let g: Vec<u8> = (0..50_000).map(|_| rng.base()).collect();
        let p = write_fasta("screen-rt.fna", &[g]);
        let sketch = ScreenSketch::build(std::slice::from_ref(&p), SCREEN_K, 50).unwrap();
        let dir = std::env::temp_dir().join(format!("sk2bgrow-{}-screendb", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        // No files yet: load reports "no sketch", the cue for the CLI error.
        assert!(ScreenSketch::load(&dir).unwrap().is_none());

        sketch.save(&dir).unwrap();
        let back = ScreenSketch::load(&dir).unwrap().expect("sketch vanished");
        assert_eq!(back, sketch);

        std::fs::remove_dir_all(&dir).ok();
        std::fs::remove_file(p).ok();
    }
}
