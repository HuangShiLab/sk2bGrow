//! Read -> anchor counting.
//!
//! Reuse target: Syn2bANI's `tag_matcher.rs` (O(1) hash lookup plus a mismatch
//! budget) and, on the Pilea side, `kmc.pyx`/`io.pyx` — except the counted unit
//! is an enzyme tag rather than a k-mer.
//!
//! ## The scan model
//!
//! Reads are scanned exactly as reference genomes are ([`crate::digest`]): slide
//! a window of each enzyme's tag length and test it against that enzyme's
//! patterns. A 2bRAD tag is not an arbitrary k-mer — it is a fixed-length window
//! satisfying a motif constraint — so this is both the correct model and cheaper
//! than hashing every k-mer of every read. It also makes route A (WMS reads) and
//! route B (real 2bRAD reads) one code path: in route B the read simply *is* the
//! tag.
//!
//! Reads may be sequenced from either strand, so a read window can be the
//! reverse complement of the reference window. Matching is therefore
//! strand-canonical throughout — the same `canonical_hash` convention
//! Fast2bRAD-M uses.
//!
//! A tag is only counted when it lies wholly inside the read. For 150 bp reads
//! and a 32 bp tag that retains (150−32+1)/150 ≈ 0.79 of the local depth — the
//! same 0.8 factor the design report uses for Pilea's k=31 sketch, so the two
//! methods are compared at matched effective depth rather than matched coverage.
//!
//! ## Mismatch tolerance
//!
//! With a budget of `m` mismatches, a tag is split into `m+1` contiguous seeds;
//! by the pigeonhole principle at least one seed survives intact, so probing all
//! `m+1` seed slots cannot miss a true match. Candidates are then verified by
//! full Hamming distance. `m = 0` skips the seed table entirely.

use std::collections::HashMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::anchor_db::AnchorDb;
use crate::enzyme::{Enzyme, PANEL};
use crate::error::Result;
use crate::seq::{canonical_hash, hamming_within, hash_tag, revcomp};

/// How reads are interpreted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CountMode {
    /// Shotgun metagenome reads: a read may span zero, one or several anchors.
    Wms,
    /// Real 2bRAD reads: each read carries exactly one tag, so counting stops at
    /// the first confirmed anchor. Extra motif hits within a 2bRAD read are
    /// artefacts (adapter chimaeras) and would inflate counts.
    TwoBrad,
}

#[derive(Debug, Clone)]
pub struct MatchConfig {
    /// Maximum Hamming distance between read tag and reference tag.
    pub max_mismatch: u32,
    pub mode: CountMode,
    /// Record a tag on every anchor it matches at the best distance.
    ///
    /// This must stay `true` for the pipeline to work. A tag shared between two
    /// genomes is precisely the input [`crate::em`] reassigns; discarding it
    /// here would leave every shared anchor at zero and give the EM nothing to
    /// split. Anchors that are shared or repeated are already flagged at build
    /// time and excluded from coverage modelling, so recording them costs
    /// nothing and enables the reassignment.
    ///
    /// Setting it to `false` records only tags that resolve to a single locus,
    /// which is useful for diagnosing how much signal the shared fraction
    /// carries.
    pub keep_multimappers: bool,
}

impl Default for MatchConfig {
    fn default() -> Self {
        MatchConfig {
            max_mismatch: 2,
            mode: CountMode::Wms,
            keep_multimappers: true,
        }
    }
}

/// Lookup structure over an [`AnchorDb`].
pub struct AnchorIndex<'a> {
    db: &'a AnchorDb,
    /// Canonical tag hash -> anchor indices (exact matching).
    exact: HashMap<u64, Vec<u32>>,
    /// (seed slot, seed hash) -> anchor indices (mismatch-tolerant matching).
    seeds: HashMap<(u8, u64), Vec<u32>>,
    n_seeds: u8,
    /// Panel entries actually present in this database.
    /// Active enzymes grouped by tag length. Two enzymes of the same tag length
    /// can satisfy their patterns on the *same physical read window* — the panel
    /// guarantees it, because every Bsp24I tag is byte-identical to a CjePI tag.
    /// Scanning per enzyme would then visit that window once per enzyme and
    /// credit every anchor at the locus once per visit; scanning per length
    /// visits it once, which is the correct model: one physical tag, one
    /// observation, credited to each stratum the locus belongs to.
    by_len: Vec<(usize, Vec<&'static Enzyme>)>,
}

impl<'a> AnchorIndex<'a> {
    /// Build the index. Cost is O(n_anchors × (m+1)).
    pub fn build(db: &'a AnchorDb, max_mismatch: u32) -> Self {
        let n_seeds = (max_mismatch + 1) as u8;
        let mut exact: HashMap<u64, Vec<u32>> = HashMap::with_capacity(db.anchors.len());
        let mut seeds: HashMap<(u8, u64), Vec<u32>> = HashMap::new();
        for (i, a) in db.anchors.iter().enumerate() {
            exact.entry(a.seq_hash).or_default().push(i as u32);
            if max_mismatch > 0 {
                let tag = db.tag(i);
                for (slot, (lo, hi)) in seed_ranges(tag.len(), n_seeds as usize)
                    .into_iter()
                    .enumerate()
                {
                    seeds
                        .entry((slot as u8, hash_tag(&tag[lo..hi])))
                        .or_default()
                        .push(i as u32);
                }
            }
        }
        let mut by_len: Vec<(usize, Vec<&'static Enzyme>)> = Vec::new();
        for e in PANEL.iter().filter(|e| db.params.enzymes.contains(e.idx)) {
            let len = e.tag_len as usize;
            match by_len.iter_mut().find(|(l, _)| *l == len) {
                Some((_, v)) => v.push(e),
                None => by_len.push((len, vec![e])),
            }
        }
        AnchorIndex {
            db,
            exact,
            seeds,
            n_seeds,
            by_len,
        }
    }

    pub fn db(&self) -> &'a AnchorDb {
        self.db
    }

    /// Anchors whose tag is within `budget` mismatches of `query`, together with
    /// the distance.
    ///
    /// Both orientations of `query` are tried, because a read may be sequenced
    /// from either strand: a read window can be the reverse complement of the
    /// reference window that produced the anchor. This is the same
    /// strand-canonical convention Fast2bRAD-M uses (`canonical_hash` takes the
    /// lexicographically smaller of a sequence and its reverse complement).
    pub fn lookup(&self, query: &[u8], budget: u32, out: &mut Vec<(u32, u32)>) {
        out.clear();
        // The exact hash is strand-canonical, so one bucket lookup serves both
        // orientations; only the verification has to consider each.
        let rc = revcomp(query);
        if let Some(hits) = self.exact.get(&canonical_hash(query)) {
            for &i in hits {
                // Verify against the stored bases: a hash collision would
                // otherwise become a phantom count.
                let tag = self.db.tag(i as usize);
                if tag == query || tag == rc {
                    out.push((i, 0));
                }
            }
            if !out.is_empty() {
                return;
            }
        }
        if budget == 0 {
            return;
        }
        let mut seen: Vec<u32> = Vec::new();
        for probe in [query, rc.as_slice()] {
            for (slot, (lo, hi)) in seed_ranges(probe.len(), self.n_seeds as usize)
                .into_iter()
                .enumerate()
            {
                let Some(cands) = self.seeds.get(&(slot as u8, hash_tag(&probe[lo..hi]))) else {
                    continue;
                };
                for &i in cands {
                    if seen.contains(&i) {
                        continue;
                    }
                    seen.push(i);
                    let tag = self.db.tag(i as usize);
                    if let Some(d) = hamming_within(&tag, probe, budget) {
                        out.push((i, d));
                    }
                }
            }
        }
    }

    /// Count anchors hit by one read. Returns the number of anchor observations
    /// recorded (0, or 1 in `TwoBrad` mode).
    pub fn count_read(
        &self,
        read: &[u8],
        cfg: &MatchConfig,
        counts: &mut [u32],
        stats: &mut CountStats,
    ) -> u32 {
        let mut recorded = 0u32;
        let mut hits: Vec<(u32, u32)> = Vec::new();
        for (len, enzymes) in self.by_len.iter() {
            let len = *len;
            if read.len() < len {
                continue;
            }
            for start in 0..=(read.len() - len) {
                let window = &read[start..start + len];
                // One physical window, tested against every enzyme of this tag
                // length. `any` rather than a per-enzyme loop is the whole point:
                // a window satisfying two patterns is still one observation.
                if !enzymes.iter().any(|e| e.match_window(window).is_some()) {
                    continue;
                }
                stats.motif_hits += 1;
                if window
                    .iter()
                    .any(|b| !matches!(b, b'A' | b'C' | b'G' | b'T'))
                {
                    stats.tag_ambiguous += 1;
                    continue;
                }
                self.lookup(window, cfg.max_mismatch, &mut hits);
                if hits.is_empty() {
                    stats.tag_unmatched += 1;
                    continue;
                }
                let best = hits.iter().map(|&(_, d)| d).min().unwrap();
                let best_hits: Vec<u32> = hits
                    .iter()
                    .filter(|&&(_, d)| d == best)
                    .map(|&(i, _)| i)
                    .collect();

                // Several anchors at one genomic locus is not the same thing as
                // several loci. The panel contains a containment relation —
                // every Bsp24I tag is byte-identical to a CjePI tag — so one
                // physical tag legitimately belongs to two enzyme strata.
                let mut loci: Vec<(u32, u16, u64)> = best_hits
                    .iter()
                    .map(|&i| {
                        let a = &self.db.anchors[i as usize];
                        (a.genome_id, a.contig_id, a.position)
                    })
                    .collect();
                loci.sort_unstable();
                loci.dedup();
                if loci.len() > 1 {
                    stats.tag_multi_locus += 1;
                    if !cfg.keep_multimappers {
                        continue;
                    }
                } else if best_hits.len() > 1 {
                    stats.tag_multi_enzyme += 1;
                }
                for i in &best_hits {
                    counts[*i as usize] += 1;
                }
                stats.tag_matched += 1;
                stats.mismatch_hist[best.min(3) as usize] += 1;
                recorded += 1;
                if cfg.mode == CountMode::TwoBrad {
                    return recorded;
                }
            }
        }
        recorded
    }
}

/// Split a tag of `len` bp into `n` contiguous seeds. Remainder bases go to the
/// trailing seeds, so every seed is `len/n` or `len/n + 1` bp.
pub fn seed_ranges(len: usize, n: usize) -> Vec<(usize, usize)> {
    let n = n.max(1);
    let base = len / n;
    let rem = len % n;
    let mut out = Vec::with_capacity(n);
    let mut lo = 0usize;
    for k in 0..n {
        let sz = base + usize::from(k >= n - rem);
        out.push((lo, lo + sz));
        lo += sz;
    }
    out
}

/// Diagnostics for one sample. Written next to the count table so QC can tell
/// "no reads matched" apart from "reads matched but coverage is low".
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CountStats {
    pub reads_total: u64,
    pub reads_with_anchor: u64,
    pub motif_hits: u64,
    /// Retained for format stability. The window scan cannot produce a
    /// truncated tag: a window shorter than `tag_len` never matches a pattern.
    pub tag_truncated: u64,
    /// Tag contained an ambiguous base.
    pub tag_ambiguous: u64,
    /// Tag matched no anchor within the mismatch budget.
    pub tag_unmatched: u64,
    /// Tag matched anchors at more than one genomic locus — a genome-shared tag
    /// or an intra-genome repeat. These are the counts [`crate::em`] reassigns.
    pub tag_multi_locus: u64,
    /// Tag matched several anchors at one locus, i.e. overlapping enzyme
    /// definitions. Not ambiguity in any harmful sense; recorded so the
    /// correlation it induces between those enzyme strata is visible.
    pub tag_multi_enzyme: u64,
    pub tag_matched: u64,
    /// Counts at 0, 1, 2, >=3 mismatches.
    pub mismatch_hist: [u64; 4],
}

impl CountStats {
    /// Tags fully extracted from reads, i.e. motif hits that were not truncated
    /// by a read end or rejected for an ambiguous base.
    pub fn tags_extracted(&self) -> u64 {
        self.motif_hits
            .saturating_sub(self.tag_truncated + self.tag_ambiguous)
    }

    /// Fraction of extracted tags that found an anchor. This is the reference-
    /// distance diagnostic: a low value with many motif hits means the reads come
    /// from something too diverged from the reference (report defect D5).
    ///
    /// Deliberately *not* divided by `motif_hits`: with the 16-enzyme union a
    /// 150 bp read spans several anchors and the ones at its edges are always
    /// truncated, so a motif-hit denominator would report a low rate for a
    /// perfectly healthy run.
    pub fn resolved_rate(&self) -> f64 {
        let n = self.tags_extracted();
        if n == 0 {
            f64::NAN
        } else {
            self.tag_matched as f64 / n as f64
        }
    }

    /// Fraction of motif hits that produced a counted tag. Lower than
    /// [`Self::resolved_rate`] by construction; useful for spotting read lengths
    /// too short for the panel's tags.
    pub fn match_rate(&self) -> f64 {
        if self.motif_hits == 0 {
            f64::NAN
        } else {
            self.tag_matched as f64 / self.motif_hits as f64
        }
    }
}

/// Count one sample's reads against the database.
pub fn count_sample(
    index: &AnchorIndex<'_>,
    reads: &[std::path::PathBuf],
    cfg: &MatchConfig,
) -> Result<(Vec<u32>, CountStats)> {
    let mut counts = vec![0u32; index.db().n_anchors()];
    let mut stats = CountStats::default();
    for path in reads {
        let n = crate::fasta::for_each_read(path, |read| {
            stats.reads_total += 1;
            let rec = index.count_read(read, cfg, &mut counts, &mut stats);
            if rec > 0 {
                stats.reads_with_anchor += 1;
            }
        })?;
        debug_assert!(n <= stats.reads_total);
    }
    Ok((counts, stats))
}

/// Column header of the count table — the Rust/Python interface contract.
pub const COUNT_TABLE_HEADER: &str = "sample\tgenome_id\tgenome\tcontig_id\tposition\tglobal_position\tenzyme\tstrand\tflags\tlocal_gc\twindow_id\tcount";

/// Write the per-anchor count table as TSV.
///
/// One row per anchor **that the database considers usable**, plus shared
/// anchors when `include_masked` is set (the EM step needs them). Masked-out
/// anchors are still written with their flags so the Python layer can audit what
/// was excluded rather than inferring it from a row count.
pub fn write_count_table(
    path: &Path,
    sample: &str,
    db: &AnchorDb,
    counts: &[u32],
    window_ids: &[u32],
    include_masked: bool,
) -> Result<()> {
    use std::io::Write;
    let f = std::fs::File::create(path).map_err(|e| crate::error::Sk2bError::Io {
        path: path.into(),
        source: e,
    })?;
    let mut w = std::io::BufWriter::with_capacity(1 << 20, f);
    let io = |e: std::io::Error| crate::error::Sk2bError::Io {
        path: path.to_path_buf(),
        source: e,
    };
    writeln!(w, "{COUNT_TABLE_HEADER}").map_err(io)?;
    for (i, a) in db.anchors.iter().enumerate() {
        if !include_masked && !a.is_usable() {
            continue;
        }
        let gname = db
            .genome(a.genome_id)
            .map(|g| g.name.as_str())
            .unwrap_or("?");
        let gpos = db.global_position(a).unwrap_or(a.position);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sample,
            a.genome_id,
            gname,
            a.contig_id,
            a.position,
            gpos,
            crate::enzyme::by_idx(a.enzyme_idx)
                .map(|e| e.name)
                .unwrap_or("?"),
            if a.strand == 0 { '+' } else { '-' },
            a.flags,
            a.local_gc,
            window_ids.get(i).copied().unwrap_or(u32::MAX),
            counts[i]
        )
        .map_err(io)?;
    }
    w.flush().map_err(io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchor_db::{assemble, build_genome, BuildParams, GC_FLANK};
    use crate::digest::DigestConfig;
    use crate::enzyme::{by_name, EnzymeSet};

    /// Unique per call: tests run in parallel and share a temp directory, so a
    /// PID-only filename lets one test delete the fixture another is reading.
    fn tiny_db() -> (AnchorDb, Vec<u8>) {
        use std::sync::atomic::{AtomicU32, Ordering};
        static N: AtomicU32 = AtomicU32::new(0);
        let mut p = std::env::temp_dir();
        p.push(format!(
            "sk2bgrow-{}-{}-count.fna",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        let mut seq = String::new();
        seq.push_str(&"A".repeat(60));
        seq.push_str("CGAACGTACTGC"); // BcgI site at 60
        seq.push_str(&"T".repeat(300));
        seq.push_str("CGAGGGTTCTGC"); // BcgI site at 372
        seq.push_str(&"C".repeat(160)); // clears the 500 bp min_contig_len default
        std::fs::write(&p, format!(">c0\n{seq}\n")).unwrap();
        let enzymes = vec![by_name("BcgI").unwrap()];
        let (meta, anchors, tags, _) =
            build_genome(&p, 0, &enzymes, &DigestConfig::default(), GC_FLANK).unwrap();
        let params = BuildParams {
            enzymes: EnzymeSet::from_slice(&enzymes),
            ..BuildParams::default()
        };
        std::fs::remove_file(&p).ok();
        (
            assemble(params, vec![(meta, anchors, tags)]),
            seq.into_bytes(),
        )
    }

    /// Regression: one physical read window that satisfies two enzymes' patterns
    /// is ONE observation, credited once to each stratum -- not one per enzyme.
    ///
    /// The panel makes this unavoidable: every Bsp24I site is also a CjePI site
    /// (their patterns, written in opposite strand orientations, differ at a
    /// single position). Scanning per enzyme visited such a window twice and
    /// incremented every anchor at the locus on each visit, inflating CjePI by
    /// 21% on E. coli and distorting its coverage profile -- exactly the enzyme
    /// carrying the second-largest weight in the fusion.
    #[test]
    fn a_shared_locus_is_counted_once_per_enzyme_not_once_per_pass() {
        use std::sync::atomic::{AtomicU32, Ordering};
        static N: AtomicU32 = AtomicU32::new(0);

        // 27 bp: GAC at offset 8 and TGG at 17 satisfies Bsp24I pattern 0, and
        // the same bases satisfy CjePI pattern 1 (GA at 8, TGG at 17).
        let tag = "TTTTTTTTGACTTTTTTTGGTTTTTTT";
        assert_eq!(tag.len(), 27);
        let seq = format!("{}{}{}", "A".repeat(300), tag, "A".repeat(300));

        let mut path = std::env::temp_dir();
        path.push(format!(
            "sk2bgrow-{}-{}-shared.fna",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, format!(">c0\n{seq}\n")).unwrap();
        let enzymes = vec![by_name("Bsp24I").unwrap(), by_name("CjePI").unwrap()];
        let (meta, anchors, tags, _) =
            build_genome(&path, 0, &enzymes, &DigestConfig::default(), GC_FLANK).unwrap();
        let params = BuildParams {
            enzymes: EnzymeSet::from_slice(&enzymes),
            ..BuildParams::default()
        };
        std::fs::remove_file(&path).ok();
        let db = assemble(params, vec![(meta, anchors, tags)]);

        // Both enzymes must have found the locus, or the test proves nothing.
        assert_eq!(db.anchors.len(), 2, "expected one anchor per enzyme");
        assert_eq!(db.tag(0), db.tag(1), "the two anchors must share a tag");

        let cfg = MatchConfig::default();
        let index = AnchorIndex::build(&db, cfg.max_mismatch);
        let mut counts = vec![0u32; db.anchors.len()];
        let mut stats = CountStats::default();
        let read = format!("GGGGG{tag}GGGGG");
        index.count_read(read.as_bytes(), &cfg, &mut counts, &mut stats);

        assert_eq!(
            counts,
            vec![1, 1],
            "one read over one shared locus must give each enzyme exactly one count"
        );
    }

    #[test]
    fn seed_ranges_cover_the_tag_exactly() {
        for len in [25usize, 27, 28, 32, 33] {
            for n in 1..=3 {
                let r = seed_ranges(len, n);
                assert_eq!(r.len(), n);
                assert_eq!(r[0].0, 0);
                assert_eq!(r.last().unwrap().1, len);
                for w in r.windows(2) {
                    assert_eq!(w[0].1, w[1].0, "seeds must be contiguous");
                }
                // Pigeonhole needs every seed non-empty.
                assert!(
                    r.iter().all(|&(lo, hi)| hi > lo),
                    "empty seed for len={len} n={n}"
                );
            }
        }
    }

    #[test]
    fn exact_read_matches_its_anchor() {
        let (db, genome) = tiny_db();
        let idx = AnchorIndex::build(&db, 2);
        // A 150 bp read fully containing the first anchor's 32 bp tag (50..82).
        let read = &genome[20..170];
        let mut counts = vec![0u32; db.n_anchors()];
        let mut stats = CountStats::default();
        let n = idx.count_read(read, &MatchConfig::default(), &mut counts, &mut stats);
        assert_eq!(n, 1, "the planted anchor was not counted");
        assert_eq!(counts.iter().sum::<u32>(), 1);
        assert_eq!(stats.mismatch_hist[0], 1);
    }

    #[test]
    fn reverse_complement_read_matches_too() {
        let (db, genome) = tiny_db();
        let idx = AnchorIndex::build(&db, 2);
        let read = revcomp(&genome[20..170]);
        let mut counts = vec![0u32; db.n_anchors()];
        let mut stats = CountStats::default();
        assert_eq!(
            idx.count_read(&read, &MatchConfig::default(), &mut counts, &mut stats),
            1
        );
    }

    #[test]
    fn flank_mismatches_are_tolerated_within_budget() {
        let (db, genome) = tiny_db();
        let idx = AnchorIndex::build(&db, 2);
        let mut read = genome[20..170].to_vec();
        // Mutate two bases inside the tag's flank (tag spans 50..82 in genome
        // coordinates, i.e. 30..62 in read coordinates); keep the motif intact.
        read[31] = if read[31] == b'A' { b'C' } else { b'A' };
        read[33] = if read[33] == b'A' { b'G' } else { b'A' };
        let mut counts = vec![0u32; db.n_anchors()];
        let mut stats = CountStats::default();
        assert_eq!(
            idx.count_read(&read, &MatchConfig::default(), &mut counts, &mut stats),
            1
        );
        assert_eq!(stats.mismatch_hist[2], 1, "expected a 2-mismatch hit");

        // A third mismatch pushes it past the budget.
        read[35] = if read[35] == b'A' { b'T' } else { b'A' };
        let mut counts2 = vec![0u32; db.n_anchors()];
        let mut stats2 = CountStats::default();
        assert_eq!(
            idx.count_read(&read, &MatchConfig::default(), &mut counts2, &mut stats2),
            0
        );
        assert_eq!(stats2.tag_unmatched, 1);
    }

    #[test]
    fn reads_match_regardless_of_sequencing_strand() {
        // A read may be sequenced from either strand, so its window can be the
        // reverse complement of the reference window. AlfI is the sharpest case:
        // its window pattern is its own reverse complement, so the scan matches
        // either way and only canonical comparison links the two readings.
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-pal.fna", std::process::id()));
        let mut seq = String::new();
        seq.push_str(&"A".repeat(60));
        seq.push_str("GCAACGTTCTGC"); // AlfI site at 60
        seq.push_str(&"TTCGA".repeat(120));
        std::fs::write(&p, format!(">c0\n{seq}\n")).unwrap();
        let enzymes = vec![by_name("AlfI").unwrap()];
        let (meta, anchors, tags, _) =
            build_genome(&p, 0, &enzymes, &DigestConfig::default(), GC_FLANK).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(!anchors.is_empty(), "fixture has no AlfI anchor");
        let params = BuildParams {
            enzymes: EnzymeSet::from_slice(&enzymes),
            ..BuildParams::default()
        };
        let db = assemble(params, vec![(meta, anchors, tags)]);
        assert!(by_name("AlfI").unwrap().is_self_complementary());

        let genome = seq.as_bytes();
        let idx = AnchorIndex::build(&db, 2);
        let read = &genome[20..170];
        for (label, r) in [("forward", read.to_vec()), ("reverse", revcomp(read))] {
            let mut counts = vec![0u32; db.n_anchors()];
            let mut stats = CountStats::default();
            let n = idx.count_read(&r, &MatchConfig::default(), &mut counts, &mut stats);
            assert_eq!(n, 1, "{label} read did not match the palindromic anchor");
        }
    }

    #[test]
    fn shared_anchors_keep_their_counts_for_the_em() {
        // Two genomes carrying the same tag: both anchors must be incremented,
        // or em::reassign has nothing to split.
        let (db1, genome) = tiny_db();
        let mut db = db1.clone();
        let mut dup: Vec<_> = db
            .anchors
            .iter()
            .map(|a| {
                let mut b = *a;
                b.genome_id = 1;
                b
            })
            .collect();
        let dup_tags = db.tags.clone();
        db.anchors.append(&mut dup);
        db.tags.extend(dup_tags);
        db.genomes.push(crate::anchor_db::GenomeMeta {
            id: 1,
            ..db.genomes[0].clone()
        });
        db.recompute_uniqueness();

        let idx = AnchorIndex::build(&db, 2);
        let read = &genome[20..170];
        let mut counts = vec![0u32; db.n_anchors()];
        let mut stats = CountStats::default();
        idx.count_read(read, &MatchConfig::default(), &mut counts, &mut stats);
        assert_eq!(stats.tag_multi_locus, 1);
        assert_eq!(
            counts.iter().sum::<u32>(),
            2,
            "the shared tag reached only one genome"
        );
        assert!(counts.iter().filter(|&&c| c == 1).count() == 2);
    }

    #[test]
    fn discarding_multimappers_is_opt_in() {
        let (db1, genome) = tiny_db();
        let mut db = db1.clone();
        let mut dup: Vec<_> = db
            .anchors
            .iter()
            .map(|a| {
                let mut b = *a;
                b.genome_id = 1;
                b
            })
            .collect();
        let dup_tags = db.tags.clone();
        db.anchors.append(&mut dup);
        db.tags.extend(dup_tags);
        db.genomes.push(crate::anchor_db::GenomeMeta {
            id: 1,
            ..db.genomes[0].clone()
        });
        db.recompute_uniqueness();

        let idx = AnchorIndex::build(&db, 2);
        let cfg = MatchConfig {
            keep_multimappers: false,
            ..MatchConfig::default()
        };
        let mut counts = vec![0u32; db.n_anchors()];
        let mut stats = CountStats::default();
        idx.count_read(&genome[20..170], &cfg, &mut counts, &mut stats);
        assert_eq!(counts.iter().sum::<u32>(), 0);
        assert_eq!(stats.tag_multi_locus, 1);
    }

    #[test]
    fn resolved_rate_ignores_edge_truncated_tags() {
        let s = CountStats {
            motif_hits: 100,
            tag_truncated: 30,
            tag_ambiguous: 0,
            tag_matched: 70,
            ..Default::default()
        };
        assert_eq!(s.tags_extracted(), 70);
        assert!((s.resolved_rate() - 1.0).abs() < 1e-12);
        assert!((s.match_rate() - 0.7).abs() < 1e-12);
    }

    #[test]
    fn zero_budget_rejects_any_mismatch() {
        let (db, genome) = tiny_db();
        let idx = AnchorIndex::build(&db, 0);
        let mut read = genome[20..170].to_vec();
        read[31] = if read[31] == b'A' { b'C' } else { b'A' };
        let cfg = MatchConfig {
            max_mismatch: 0,
            ..MatchConfig::default()
        };
        let mut counts = vec![0u32; db.n_anchors()];
        let mut stats = CountStats::default();
        assert_eq!(idx.count_read(&read, &cfg, &mut counts, &mut stats), 0);
    }

    #[test]
    fn a_partial_tag_at_a_read_edge_simply_does_not_match() {
        // In the window model there is no "truncated tag" to report: a window
        // shorter than tag_len can never satisfy a pattern, so a read carrying
        // only part of a tag produces no hit at all rather than a near-miss.
        let (db, genome) = tiny_db();
        let idx = AnchorIndex::build(&db, 2);
        // Anchor 0's window is genome[50..82]; start the read inside it.
        let read = &genome[58..120];
        let mut counts = vec![0u32; db.n_anchors()];
        let mut stats = CountStats::default();
        assert_eq!(
            idx.count_read(read, &MatchConfig::default(), &mut counts, &mut stats),
            0
        );
        assert_eq!(stats.motif_hits, 0);
        assert_eq!(stats.tag_unmatched, 0);
        assert_eq!(
            stats.tag_truncated, 0,
            "the field is vestigial in this model"
        );
    }

    #[test]
    fn two_brad_mode_counts_one_tag_per_read() {
        let (db, genome) = tiny_db();
        let idx = AnchorIndex::build(&db, 2);
        // A chimaeric read carrying both tags; WMS mode counts two, 2bRAD one.
        let mut read = genome[50..82].to_vec();
        read.extend_from_slice(&genome[362..394]);
        let mut c1 = vec![0u32; db.n_anchors()];
        let mut s1 = CountStats::default();
        assert_eq!(
            idx.count_read(&read, &MatchConfig::default(), &mut c1, &mut s1),
            2
        );

        let cfg = MatchConfig {
            mode: CountMode::TwoBrad,
            ..MatchConfig::default()
        };
        let mut c2 = vec![0u32; db.n_anchors()];
        let mut s2 = CountStats::default();
        assert_eq!(idx.count_read(&read, &cfg, &mut c2, &mut s2), 1);
        assert_eq!(c2.iter().sum::<u32>(), 1);
    }
}
