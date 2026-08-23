//! The anchor database: the build-time product that every analysis reads.
//!
//! Corresponds to Pilea's `sketch.py` (index side), with two substitutions the
//! design report argues for: random FracMinHash sampling becomes deterministic
//! 16-enzyme digestion, and the 25 kb fixed window becomes an equal-anchor-count
//! adaptive window (see [`crate::window`]).
//!
//! Uniqueness masking is computed **once at build time** so analysis-time
//! filtering is free — the same strategy Pilea uses for multi-copy k-mers.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::digest::{digest_contig, DigestConfig, Site};
use crate::enzyme::{by_idx, Enzyme, EnzymeSet};
use crate::error::{Result, Sk2bError};
use crate::fasta;
use crate::seq::{gc_fraction, quantize_gc, GC_UNDEFINED};
use crate::tgt::{pack_bases, unpack_bases, ContigKind, ContigMeta, Tgt, TgtRecord};

/// Flank half-width used for the `local_gc` field, in bp (report §7.1: ±250 bp).
pub const GC_FLANK: usize = 250;

/// Anchor flag bits.
pub mod flags {
    /// The tag sequence occurs exactly once in its own genome.
    pub const UNIQUE_IN_GENOME: u8 = 1 << 0;
    /// The tag sequence occurs in exactly one genome of the database.
    pub const UNIQUE_ACROSS_DB: u8 = 1 << 1;
    /// Masked: repeated within its own genome (multi-copy contamination).
    pub const MASKED_MULTICOPY: u8 = 1 << 2;
    /// Masked: shared with another genome. Such anchors still enter the EM
    /// reassignment step; they are only excluded from coverage modelling.
    pub const MASKED_SHARED: u8 = 1 << 3;
    /// Sits on a plasmid or other non-chromosomal contig — no ori-ter gradient
    /// (report §8.2).
    pub const NON_CHROMOSOMAL: u8 = 1 << 4;
    /// The ±250 bp GC window was undefined (all `N`).
    pub const GC_UNDEFINED: u8 = 1 << 5;

    /// Anchors that may be used directly for coverage modelling.
    pub const USABLE_MASK: u8 = MASKED_MULTICOPY | MASKED_SHARED | NON_CHROMOSOMAL;
}

/// One anchor. Layout follows the architecture doc §5, extending the Syn2b 48 B
/// tag with `flags` and `local_gc`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Anchor {
    /// Strand-canonical hash of the tag sequence.
    pub seq_hash: u64,
    pub genome_id: u32,
    pub contig_id: u16,
    pub position: u64,
    pub enzyme_idx: u8,
    pub strand: u8,
    pub flags: u8,
    /// Quantised ±[`GC_FLANK`] bp GC, 255 = undefined.
    pub local_gc: u8,
}

impl Anchor {
    /// True when this anchor can enter coverage modelling unaltered.
    #[inline]
    pub fn is_usable(&self) -> bool {
        self.flags & flags::USABLE_MASK == 0
    }
    /// Tag length in bp, derived from the enzyme rather than stored.
    #[inline]
    pub fn tag_len(&self) -> usize {
        by_idx(self.enzyme_idx)
            .map(|e| e.tag_len as usize)
            .unwrap_or(0)
    }
}

/// Per-genome metadata.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenomeMeta {
    pub id: u32,
    pub name: String,
    /// Free-form lineage string from `-a taxonomy.tsv`, if supplied.
    pub taxonomy: Option<String>,
    pub contigs: Vec<ContigMeta>,
    pub genome_len: u64,
    /// Global ori coordinate when known (DoriC/Ori-Finder annotation, or a
    /// previous grid search). `None` means the fit must search for it.
    pub ori: Option<u64>,
    /// 0..1 confidence in `ori`; 1.0 for a curated annotation.
    pub ori_confidence: f64,
}

impl GenomeMeta {
    /// Whether the reference is contiguous enough for real-coordinate fitting.
    /// Pilea's own sensitivity analysis puts the boundary near 100 contigs /
    /// N50 > 10 kbp (report §1.2).
    pub fn is_contiguous(&self) -> bool {
        self.contigs.len() <= 100
    }
}

/// Build-time parameters, persisted so an analysis can assert compatibility.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BuildParams {
    pub enzymes: EnzymeSet,
    pub gc_flank: usize,
    pub min_contig_len: usize,
    pub reject_ambiguous_tags: bool,
    pub sk2bgrow_version: String,
}

impl Default for BuildParams {
    fn default() -> Self {
        BuildParams {
            enzymes: EnzymeSet(0),
            gc_flank: GC_FLANK,
            min_contig_len: 500,
            reject_ambiguous_tags: true,
            sk2bgrow_version: env!("CARGO_PKG_VERSION").to_string(),
        }
    }
}

/// The database. Anchors are stored sorted by `(genome_id, contig_id, position)`
/// so per-genome slices are contiguous and windowing can stream them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AnchorDb {
    pub params: BuildParams,
    pub genomes: Vec<GenomeMeta>,
    pub anchors: Vec<Anchor>,
    /// Packed tag bases, parallel to `anchors`. Needed for mismatch-tolerant
    /// verification during counting.
    pub tags: Vec<[u8; 12]>,
}

impl AnchorDb {
    pub fn n_anchors(&self) -> usize {
        self.anchors.len()
    }

    /// Tag bases of anchor `i`.
    pub fn tag(&self, i: usize) -> Vec<u8> {
        unpack_bases(&self.tags[i], self.anchors[i].tag_len())
    }

    /// Half-open index range of one genome's anchors.
    pub fn genome_range(&self, genome_id: u32) -> std::ops::Range<usize> {
        let lo = self.anchors.partition_point(|a| a.genome_id < genome_id);
        let hi = self.anchors.partition_point(|a| a.genome_id <= genome_id);
        lo..hi
    }

    pub fn genome(&self, genome_id: u32) -> Option<&GenomeMeta> {
        self.genomes.iter().find(|g| g.id == genome_id)
    }

    /// Global coordinate of an anchor (contig offset + position).
    pub fn global_position(&self, a: &Anchor) -> Option<u64> {
        let g = self.genome(a.genome_id)?;
        g.contigs
            .iter()
            .find(|c| c.id == a.contig_id)
            .map(|c| c.offset + a.position)
    }

    /// Recompute `UNIQUE_IN_GENOME` / `UNIQUE_ACROSS_DB` and the corresponding
    /// mask bits from scratch.
    ///
    /// Must be called after every mutation that adds or removes anchors —
    /// uniqueness is a property of the whole database, not of one genome.
    pub fn recompute_uniqueness(&mut self) {
        // hash -> (occurrences, first genome, saw a second distinct genome)
        let mut seen: HashMap<u64, (u32, u32, bool)> = HashMap::with_capacity(self.anchors.len());
        for a in &self.anchors {
            match seen.get_mut(&a.seq_hash) {
                None => {
                    seen.insert(a.seq_hash, (1, a.genome_id, false));
                }
                Some((n, g0, multi)) => {
                    *n += 1;
                    if *g0 != a.genome_id {
                        *multi = true;
                    }
                }
            }
        }
        // Count per-genome occurrences separately: a tag can be unique across the
        // database yet repeated inside its own genome (a tandem duplication).
        let mut per_genome: HashMap<(u32, u64), u32> = HashMap::with_capacity(self.anchors.len());
        for a in &self.anchors {
            *per_genome.entry((a.genome_id, a.seq_hash)).or_insert(0) += 1;
        }

        for a in self.anchors.iter_mut() {
            let (_, _, multi_genome) = seen[&a.seq_hash];
            let in_genome = per_genome[&(a.genome_id, a.seq_hash)];
            a.flags &= !(flags::UNIQUE_IN_GENOME
                | flags::UNIQUE_ACROSS_DB
                | flags::MASKED_MULTICOPY
                | flags::MASKED_SHARED);
            if in_genome == 1 {
                a.flags |= flags::UNIQUE_IN_GENOME;
            } else {
                a.flags |= flags::MASKED_MULTICOPY;
            }
            if !multi_genome {
                a.flags |= flags::UNIQUE_ACROSS_DB;
            } else {
                a.flags |= flags::MASKED_SHARED;
            }
        }
    }

    /// Anchors sharing a tag with another genome, grouped by hash. This is the
    /// input to [`crate::em`].
    pub fn shared_groups(&self) -> HashMap<u64, Vec<usize>> {
        let mut by_hash: HashMap<u64, Vec<usize>> = HashMap::new();
        for (i, a) in self.anchors.iter().enumerate() {
            if a.flags & flags::MASKED_SHARED != 0 {
                by_hash.entry(a.seq_hash).or_default().push(i);
            }
        }
        by_hash.retain(|_, v| {
            let g0 = v.first().map(|&i| self.anchors[i].genome_id);
            v.iter().any(|&i| Some(self.anchors[i].genome_id) != g0)
        });
        by_hash
    }

    /// Per-enzyme usable-anchor counts for one genome. A genome with too few
    /// anchors under some enzyme must be degraded to an enzyme subset rather
    /// than reported (report §7.1 step 5).
    pub fn per_enzyme_counts(&self, genome_id: u32) -> [usize; crate::enzyme::N_ENZYMES] {
        let mut out = [0usize; crate::enzyme::N_ENZYMES];
        for a in &self.anchors[self.genome_range(genome_id)] {
            if a.is_usable() {
                out[a.enzyme_idx as usize] += 1;
            }
        }
        out
    }

    // ---- persistence -----------------------------------------------------

    /// Write the database as a directory: `manifest.json` + `anchors.bin`.
    ///
    /// The manifest is deliberately human-readable — it is the contract between
    /// the Rust counting layer and the Python statistics layer.
    pub fn save(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir).map_err(|e| Sk2bError::Io {
            path: dir.into(),
            source: e,
        })?;
        let manifest = Manifest {
            format_version: 1,
            params: self.params.clone(),
            genomes: self.genomes.clone(),
            n_anchors: self.anchors.len(),
        };
        let mpath = dir.join("manifest.json");
        std::fs::write(&mpath, serde_json::to_vec_pretty(&manifest)?).map_err(|e| {
            Sk2bError::Io {
                path: mpath,
                source: e,
            }
        })?;
        let apath = dir.join("anchors.bin");
        let f = std::fs::File::create(&apath).map_err(|e| Sk2bError::Io {
            path: apath.clone(),
            source: e,
        })?;
        let mut w = std::io::BufWriter::with_capacity(1 << 20, f);
        bincode::serialize_into(&mut w, &(&self.anchors, &self.tags))?;
        use std::io::Write;
        w.flush().map_err(|e| Sk2bError::Io {
            path: apath,
            source: e,
        })?;
        Ok(())
    }

    pub fn load(dir: &Path) -> Result<Self> {
        let mpath = dir.join("manifest.json");
        let mbytes = std::fs::read(&mpath).map_err(|e| Sk2bError::Io {
            path: mpath,
            source: e,
        })?;
        let manifest: Manifest = serde_json::from_slice(&mbytes)?;
        if manifest.format_version != 1 {
            return Err(Sk2bError::Db(format!(
                "anchor database format version {} is not supported (expected 1)",
                manifest.format_version
            )));
        }
        let apath = dir.join("anchors.bin");
        let f = std::fs::File::open(&apath).map_err(|e| Sk2bError::Io {
            path: apath.clone(),
            source: e,
        })?;
        let r = std::io::BufReader::with_capacity(1 << 20, f);
        let (anchors, tags): (Vec<Anchor>, Vec<[u8; 12]>) = bincode::deserialize_from(r)?;
        if anchors.len() != tags.len() {
            return Err(Sk2bError::Db(format!(
                "corrupt database: {} anchors but {} tags",
                anchors.len(),
                tags.len()
            )));
        }
        if anchors.len() != manifest.n_anchors {
            return Err(Sk2bError::Db(format!(
                "manifest claims {} anchors, anchors.bin holds {}",
                manifest.n_anchors,
                anchors.len()
            )));
        }
        Ok(AnchorDb {
            params: manifest.params,
            genomes: manifest.genomes,
            anchors,
            tags,
        })
    }
}

#[derive(Serialize, Deserialize)]
struct Manifest {
    format_version: u32,
    params: BuildParams,
    genomes: Vec<GenomeMeta>,
    n_anchors: usize,
}

/// What [`build_genome`] returns: genome metadata, its anchors, their packed
/// tags, and the TGT record set. Aliased so the signature stays readable.
pub type GenomeBuild = (GenomeMeta, Vec<Anchor>, Vec<[u8; 12]>, Tgt);

/// Digest one genome FASTA into anchors plus a TGT record set.
///
/// `genome_id` must be unique within the database. Contig ids are assigned in
/// input order and `offset` is the running concatenation offset, which is only a
/// meaningful coordinate for a closed chromosome or a scaffolded MAG.
pub fn build_genome(
    path: &Path,
    genome_id: u32,
    enzymes: &[&'static Enzyme],
    cfg: &DigestConfig,
    gc_flank: usize,
) -> Result<GenomeBuild> {
    let records = fasta::read_fasta(path)?;
    if records.len() > u16::MAX as usize {
        return Err(Sk2bError::Db(format!(
            "{} has {} contigs, more than the u16 contig_id field allows",
            path.display(),
            records.len()
        )));
    }
    let genome_name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| format!("genome{genome_id}"));

    let mut contigs = Vec::with_capacity(records.len());
    let mut anchors = Vec::new();
    let mut tags = Vec::new();
    let mut tgt_records = Vec::new();
    let mut offset = 0u64;

    for (i, rec) in records.iter().enumerate() {
        let contig_id = i as u16;
        let kind = classify_contig(&rec.name);
        contigs.push(ContigMeta {
            id: contig_id,
            name: rec.name.clone(),
            length: rec.seq.len() as u64,
            offset,
            kind,
        });
        offset += rec.seq.len() as u64;

        let sites = digest_contig(&rec.seq, contig_id, enzymes, cfg);
        for s in &sites {
            let (anchor, packed) = anchor_from_site(s, genome_id, &rec.seq, gc_flank, kind);
            anchors.push(anchor);
            tags.push(packed);
            tgt_records.push(TgtRecord::from_site(s)?);
        }
    }

    let genome_len = contigs.iter().map(|c| c.length).sum();
    let meta = GenomeMeta {
        id: genome_id,
        name: genome_name.clone(),
        taxonomy: None,
        contigs: contigs.clone(),
        genome_len,
        ori: None,
        ori_confidence: 0.0,
    };
    let mut tgt = Tgt::new(genome_name, contigs);
    tgt.records = tgt_records;
    tgt.recompute_gaps();
    // Carry the GC annotation into the TGT copy so the text dump is complete.
    for (r, a) in tgt.records.iter_mut().zip(anchors.iter()) {
        r.local_gc = a.local_gc;
        r.flags = a.flags;
    }
    Ok((meta, anchors, tags, tgt))
}

fn anchor_from_site(
    site: &Site,
    genome_id: u32,
    contig_seq: &[u8],
    gc_flank: usize,
    kind: ContigKind,
) -> (Anchor, [u8; 12]) {
    let centre = site.site_start as usize;
    let lo = centre.saturating_sub(gc_flank);
    let hi = (centre + gc_flank).min(contig_seq.len());
    let gcq = quantize_gc(gc_fraction(&contig_seq[lo..hi]));
    let mut f = 0u8;
    if gcq == GC_UNDEFINED {
        f |= flags::GC_UNDEFINED;
    }
    if kind != ContigKind::Chromosome {
        f |= flags::NON_CHROMOSOMAL;
    }
    (
        Anchor {
            seq_hash: site.tag_hash(),
            genome_id,
            contig_id: site.contig_id,
            position: site.site_start,
            enzyme_idx: site.enzyme_idx,
            strand: site.strand.as_u8(),
            flags: f,
            local_gc: gcq,
        },
        pack_bases(&site.tag),
    )
}

/// Heuristic contig classification from the FASTA header.
///
/// Deliberately conservative: only an explicit "plasmid" marker demotes a
/// contig, because misclassifying the chromosome would silently discard the
/// entire ori-ter signal.
pub fn classify_contig(name: &str) -> ContigKind {
    let n = name.to_ascii_lowercase();
    if n.contains("plasmid") || n.starts_with('p') && n.len() <= 4 {
        ContigKind::Plasmid
    } else {
        ContigKind::Chromosome
    }
}

/// Assemble a database from per-genome build results and finalise uniqueness.
pub fn assemble(
    params: BuildParams,
    parts: Vec<(GenomeMeta, Vec<Anchor>, Vec<[u8; 12]>)>,
) -> AnchorDb {
    let mut genomes = Vec::with_capacity(parts.len());
    let mut anchors = Vec::new();
    let mut tags = Vec::new();
    for (meta, a, t) in parts {
        genomes.push(meta);
        anchors.extend(a);
        tags.extend(t);
    }
    // Sort anchors and tags together, then finalise flags.
    let mut order: Vec<usize> = (0..anchors.len()).collect();
    order.sort_by_key(|&i| {
        let a = &anchors[i];
        (a.genome_id, a.contig_id, a.position, a.enzyme_idx)
    });
    let anchors: Vec<Anchor> = order.iter().map(|&i| anchors[i]).collect();
    let tags: Vec<[u8; 12]> = order.iter().map(|&i| tags[i]).collect();
    genomes.sort_by_key(|g| g.id);

    let mut db = AnchorDb {
        params,
        genomes,
        anchors,
        tags,
    };
    db.recompute_uniqueness();
    db
}

/// Attach lineages from a two-column `name<TAB>lineage` TSV (the `-a` flag).
/// Genomes absent from the file keep `taxonomy = None`; that is a warning, not
/// an error, because a lineage is metadata rather than an input to the model.
pub fn attach_taxonomy(db: &mut AnchorDb, path: &Path) -> Result<usize> {
    let text = std::fs::read_to_string(path).map_err(|e| Sk2bError::Io {
        path: path.into(),
        source: e,
    })?;
    let mut map: HashMap<&str, &str> = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some((k, v)) = line.split_once('\t') {
            map.insert(k.trim(), v.trim());
        }
    }
    let mut hit = 0;
    for g in db.genomes.iter_mut() {
        if let Some(t) = map.get(g.name.as_str()) {
            g.taxonomy = Some((*t).to_string());
            hit += 1;
        }
    }
    Ok(hit)
}

/// Where per-genome TGT dumps live inside a database directory.
pub fn tgt_path(dir: &Path, genome_name: &str) -> PathBuf {
    dir.join("tgt").join(format!("{genome_name}.tgt"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzyme::by_name;

    fn anchor(hash: u64, genome: u32, pos: u64) -> Anchor {
        Anchor {
            seq_hash: hash,
            genome_id: genome,
            contig_id: 0,
            position: pos,
            enzyme_idx: 0,
            strand: 0,
            flags: 0,
            local_gc: 100,
        }
    }

    fn db_of(anchors: Vec<Anchor>) -> AnchorDb {
        let n = anchors.len();
        let genome_ids: Vec<u32> = {
            let mut v: Vec<u32> = anchors.iter().map(|a| a.genome_id).collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        AnchorDb {
            params: BuildParams::default(),
            genomes: genome_ids
                .into_iter()
                .map(|id| GenomeMeta {
                    id,
                    name: format!("g{id}"),
                    taxonomy: None,
                    contigs: vec![ContigMeta {
                        id: 0,
                        name: "c0".into(),
                        length: 10_000,
                        offset: 0,
                        kind: ContigKind::Chromosome,
                    }],
                    genome_len: 10_000,
                    ori: None,
                    ori_confidence: 0.0,
                })
                .collect(),
            anchors,
            tags: vec![[0u8; 12]; n],
        }
    }

    #[test]
    fn uniqueness_separates_multicopy_from_shared() {
        let mut db = db_of(vec![
            anchor(1, 0, 100), // unique everywhere
            anchor(2, 0, 200), // repeated within genome 0
            anchor(2, 0, 300),
            anchor(3, 0, 400), // shared with genome 1
            anchor(3, 1, 500),
        ]);
        db.recompute_uniqueness();
        let a: Vec<u8> = db.anchors.iter().map(|a| a.flags).collect();
        assert!(a[0] & flags::UNIQUE_IN_GENOME != 0 && a[0] & flags::UNIQUE_ACROSS_DB != 0);
        assert!(db.anchors[0].is_usable());
        // Repeated inside genome 0, but still absent from genome 1.
        assert!(a[1] & flags::MASKED_MULTICOPY != 0);
        assert!(a[1] & flags::UNIQUE_ACROSS_DB != 0);
        assert!(!db.anchors[1].is_usable());
        // Present once per genome: genome-unique but cross-genome shared.
        let shared = db.anchors.iter().position(|x| x.seq_hash == 3).unwrap();
        assert!(a[shared] & flags::UNIQUE_IN_GENOME != 0);
        assert!(a[shared] & flags::MASKED_SHARED != 0);
        assert!(!db.anchors[shared].is_usable());
    }

    #[test]
    fn shared_groups_only_span_genomes() {
        let mut db = db_of(vec![
            anchor(2, 0, 200),
            anchor(2, 0, 300),
            anchor(3, 0, 400),
            anchor(3, 1, 500),
        ]);
        db.recompute_uniqueness();
        let groups = db.shared_groups();
        assert_eq!(
            groups.len(),
            1,
            "an intra-genome duplicate leaked into the EM input"
        );
        assert!(groups.contains_key(&3));
    }

    #[test]
    fn genome_range_slices_correctly() {
        let mut db = db_of(vec![
            anchor(1, 0, 1),
            anchor(2, 1, 1),
            anchor(3, 1, 2),
            anchor(4, 2, 1),
        ]);
        db.recompute_uniqueness();
        assert_eq!(db.genome_range(1), 1..3);
        assert_eq!(db.genome_range(0), 0..1);
        assert_eq!(db.genome_range(9), 4..4);
    }

    #[test]
    fn plasmid_contigs_are_flagged_not_dropped() {
        assert_eq!(
            classify_contig("NZ_CP009273.1_plasmid_pX"),
            ContigKind::Plasmid
        );
        assert_eq!(classify_contig("NC_000913.3"), ContigKind::Chromosome);
    }

    #[test]
    fn build_and_roundtrip_a_tiny_genome() {
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-mini.fna", std::process::id()));
        // Two BcgI sites, far enough apart to stay distinct loci.
        // Long enough to clear DigestConfig::default().min_contig_len (500 bp).
        let mut body = String::from(">c0 test contig\n");
        body.push_str(&"A".repeat(60));
        body.push_str("CGAACGTACTGC");
        body.push_str(&"T".repeat(400));
        body.push_str("CGAGGGGGGTGC");
        body.push_str(&"C".repeat(200));
        body.push('\n');
        std::fs::write(&p, body).unwrap();

        let enzymes = vec![by_name("BcgI").unwrap()];
        let cfg = DigestConfig::default();
        let (meta, anchors, tags, tgt) = build_genome(&p, 0, &enzymes, &cfg, GC_FLANK).unwrap();
        assert_eq!(anchors.len(), 2, "expected two planted BcgI anchors");
        assert_eq!(tgt.records.len(), 2);
        assert_ne!(anchors[0].local_gc, GC_UNDEFINED);

        let params = BuildParams {
            enzymes: EnzymeSet::from_slice(&enzymes),
            ..BuildParams::default()
        };
        let db = assemble(params, vec![(meta, anchors, tags)]);
        assert_eq!(db.per_enzyme_counts(0)[0], 2);

        let dir = std::env::temp_dir().join(format!("sk2bgrow-{}-db", std::process::id()));
        db.save(&dir).unwrap();
        let back = AnchorDb::load(&dir).unwrap();
        assert_eq!(back, db);
        std::fs::remove_dir_all(&dir).ok();
        std::fs::remove_file(&p).ok();
    }
}
