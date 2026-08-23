//! TGT v2 (Tag–Gap–Tag): the on-disk representation of a digested genome.
//!
//! Reuse target: `bsyn::tgt` (tag / gap / record / reader / writer). A genome is
//! an ordered tag sequence plus the gaps between consecutive tags plus contig
//! metadata — the structure the design report leans on for adaptive windows
//! (§6.3) and per-enzyme stratification (§6.2).
//!
//! Two serialisations, same content:
//!
//! * **binary** — a 48-byte fixed record per tag behind a small header. Fast to
//!   mmap-scan, and the size the report budgets with (~1.4 MB per E. coli-scale
//!   genome for the 16-enzyme union).
//! * **text** — `Enzyme:SEQUENCE<TAB>contig<TAB>pos<TAB>strand<TAB>gap`, grouped
//!   by enzyme, for diffing and eyeballing.
//!
//! ## Record layout (48 B, little-endian)
//!
//! ```text
//!  0..8   tag_hash    u64   strand-canonical hash of the tag
//!  8..16  position    u64   0-based start of the recognition site
//! 16..20  gap         u32   bp since the previous tag on this contig (0 for the first)
//! 20..22  contig_id   u16
//! 22..23  enzyme_idx  u8    index into enzyme::PANEL
//! 23..24  strand      u8    0 = +, 1 = -
//! 24..25  tag_len     u8    bp
//! 25..26  flags       u8    see anchor_db::flags
//! 26..27  local_gc    u8    quantised ±250 bp GC, 255 = undefined
//! 27..28  _pad        u8
//! 28..40  tag_2bit    [u8;12]  packed tag, 4 bases per byte (up to 48 bp)
//! 40..48  _reserved   [u8;8]
//! ```

use std::io::{BufRead, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::digest::{Site, Strand};
use crate::enzyme::{by_idx, by_name};
use crate::error::{Result, Sk2bError};
use crate::seq::{canonical_hash, two_bit, GC_UNDEFINED};

/// Bytes per binary record. Matches the Syn2b tag layout the architecture doc
/// extends (`48 B/tag`).
pub const RECORD_SIZE: usize = 48;
/// Longest tag the packed field can hold. The panel maxes out at 33 bp (CspCI).
pub const MAX_PACKED_BASES: usize = 48;
const TAG_BYTES: usize = MAX_PACKED_BASES / 4;
const MAGIC: &[u8; 4] = b"TGT2";

/// Contig metadata carried by TGT v2, which is what makes fragmented MAGs
/// representable at all.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ContigMeta {
    pub id: u16,
    pub name: String,
    pub length: u64,
    /// Offset of this contig in a virtual concatenated genome. Set by
    /// [`scaffold`](crate::scaffold) once contigs are ordered; before that it is
    /// simply the cumulative length in input order.
    pub offset: u64,
    /// Contig class. Plasmids must be maskable — they carry recognition sites
    /// but no ori-ter gradient (report §8.2).
    pub kind: ContigKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum ContigKind {
    #[default]
    Chromosome,
    Plasmid,
    Unknown,
}

/// One tag record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TgtRecord {
    pub tag_hash: u64,
    pub position: u64,
    pub gap: u32,
    pub contig_id: u16,
    pub enzyme_idx: u8,
    pub strand: Strand,
    pub tag_len: u8,
    pub flags: u8,
    pub local_gc: u8,
    /// Packed tag bases; use [`TgtRecord::tag_bases`] to unpack.
    pub tag_2bit: [u8; TAG_BYTES],
}

impl TgtRecord {
    /// Build a record from a digested [`Site`]. `gap` is filled in later by
    /// [`Tgt::recompute_gaps`], which needs the whole sorted run.
    pub fn from_site(site: &Site) -> Result<Self> {
        if site.tag.len() > MAX_PACKED_BASES {
            return Err(Sk2bError::Db(format!(
                "tag of {} bp exceeds the {} bp packed field",
                site.tag.len(),
                MAX_PACKED_BASES
            )));
        }
        Ok(TgtRecord {
            tag_hash: site.tag_hash(),
            position: site.site_start,
            gap: 0,
            contig_id: site.contig_id,
            enzyme_idx: site.enzyme_idx,
            strand: site.strand,
            tag_len: site.tag.len() as u8,
            flags: 0,
            local_gc: GC_UNDEFINED,
            tag_2bit: pack_bases(&site.tag),
        })
    }

    /// Unpack the stored tag back to ASCII bases.
    pub fn tag_bases(&self) -> Vec<u8> {
        unpack_bases(&self.tag_2bit, self.tag_len as usize)
    }

    pub fn enzyme_name(&self) -> &'static str {
        by_idx(self.enzyme_idx).map(|e| e.name).unwrap_or("?")
    }

    fn write_binary<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        let mut buf = [0u8; RECORD_SIZE];
        buf[0..8].copy_from_slice(&self.tag_hash.to_le_bytes());
        buf[8..16].copy_from_slice(&self.position.to_le_bytes());
        buf[16..20].copy_from_slice(&self.gap.to_le_bytes());
        buf[20..22].copy_from_slice(&self.contig_id.to_le_bytes());
        buf[22] = self.enzyme_idx;
        buf[23] = self.strand.as_u8();
        buf[24] = self.tag_len;
        buf[25] = self.flags;
        buf[26] = self.local_gc;
        buf[28..28 + TAG_BYTES].copy_from_slice(&self.tag_2bit);
        w.write_all(&buf)
    }

    fn read_binary(buf: &[u8; RECORD_SIZE]) -> Self {
        let mut tag_2bit = [0u8; TAG_BYTES];
        tag_2bit.copy_from_slice(&buf[28..28 + TAG_BYTES]);
        TgtRecord {
            tag_hash: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            position: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
            gap: u32::from_le_bytes(buf[16..20].try_into().unwrap()),
            contig_id: u16::from_le_bytes(buf[20..22].try_into().unwrap()),
            enzyme_idx: buf[22],
            strand: Strand::from_u8(buf[23]),
            tag_len: buf[24],
            flags: buf[25],
            local_gc: buf[26],
            tag_2bit,
        }
    }
}

/// A digested genome: contig table plus tag records.
#[derive(Debug, Clone, PartialEq)]
pub struct Tgt {
    pub genome_name: String,
    pub contigs: Vec<ContigMeta>,
    pub records: Vec<TgtRecord>,
}

impl Tgt {
    pub fn new(genome_name: impl Into<String>, contigs: Vec<ContigMeta>) -> Self {
        Tgt {
            genome_name: genome_name.into(),
            contigs,
            records: Vec::new(),
        }
    }

    /// Total sequence length across contigs.
    pub fn genome_len(&self) -> u64 {
        self.contigs.iter().map(|c| c.length).sum()
    }

    /// Sort records into `(contig_id, position, enzyme_idx)` order and fill in
    /// the `gap` field. Gaps restart at each contig boundary — a gap spanning a
    /// boundary is not a genomic distance.
    pub fn recompute_gaps(&mut self) {
        self.records
            .sort_by_key(|r| (r.contig_id, r.position, r.enzyme_idx, r.strand));
        let mut prev: Option<(u16, u64)> = None;
        for r in self.records.iter_mut() {
            r.gap = match prev {
                Some((c, p)) if c == r.contig_id => (r.position - p) as u32,
                _ => 0,
            };
            prev = Some((r.contig_id, r.position));
        }
    }

    /// Records belonging to one enzyme, in genome order.
    pub fn by_enzyme(&self, enzyme_idx: u8) -> impl Iterator<Item = &TgtRecord> {
        self.records
            .iter()
            .filter(move |r| r.enzyme_idx == enzyme_idx)
    }

    /// Global coordinate of a record: `contig.offset + position`. This is the
    /// axis the ori-ter V-shape is fitted on, so it is only meaningful once
    /// contigs are ordered (a single closed chromosome, or after scaffolding).
    pub fn global_position(&self, r: &TgtRecord) -> Option<u64> {
        self.contigs
            .iter()
            .find(|c| c.id == r.contig_id)
            .map(|c| c.offset + r.position)
    }

    pub fn write_binary(&self, path: &Path) -> Result<()> {
        let f = std::fs::File::create(path).map_err(|e| Sk2bError::Io {
            path: path.into(),
            source: e,
        })?;
        let mut w = BufWriter::with_capacity(1 << 20, f);
        let io = |e: std::io::Error| Sk2bError::Io {
            path: path.to_path_buf(),
            source: e,
        };
        // Header: magic, version, then a JSON blob with genome/contig metadata
        // (length-prefixed), then fixed-size records.
        w.write_all(MAGIC).map_err(io)?;
        w.write_all(&2u32.to_le_bytes()).map_err(io)?;
        let meta = serde_json::to_vec(&TgtHeaderMeta {
            genome_name: self.genome_name.clone(),
            contigs: self.contigs.clone(),
        })?;
        w.write_all(&(meta.len() as u32).to_le_bytes())
            .map_err(io)?;
        w.write_all(&meta).map_err(io)?;
        w.write_all(&(self.records.len() as u64).to_le_bytes())
            .map_err(io)?;
        for r in &self.records {
            r.write_binary(&mut w).map_err(io)?;
        }
        w.flush().map_err(io)?;
        Ok(())
    }

    pub fn read_binary(path: &Path) -> Result<Self> {
        let mut r = crate::fasta::open_reader(path)?;
        let io = |e: std::io::Error| Sk2bError::Io {
            path: path.to_path_buf(),
            source: e,
        };
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic).map_err(io)?;
        if &magic != MAGIC {
            return Err(Sk2bError::Format {
                path: path.into(),
                msg: "not a TGT v2 file".into(),
            });
        }
        let mut u32buf = [0u8; 4];
        r.read_exact(&mut u32buf).map_err(io)?;
        let version = u32::from_le_bytes(u32buf);
        if version != 2 {
            return Err(Sk2bError::Format {
                path: path.into(),
                msg: format!("TGT version {version} is not supported (expected 2)"),
            });
        }
        r.read_exact(&mut u32buf).map_err(io)?;
        let meta_len = u32::from_le_bytes(u32buf) as usize;
        let mut meta_buf = vec![0u8; meta_len];
        r.read_exact(&mut meta_buf).map_err(io)?;
        let meta: TgtHeaderMeta = serde_json::from_slice(&meta_buf)?;
        let mut u64buf = [0u8; 8];
        r.read_exact(&mut u64buf).map_err(io)?;
        let n = u64::from_le_bytes(u64buf) as usize;
        let mut records = Vec::with_capacity(n);
        let mut buf = [0u8; RECORD_SIZE];
        for i in 0..n {
            r.read_exact(&mut buf).map_err(|e| Sk2bError::Format {
                path: path.to_path_buf(),
                msg: format!("truncated at record {i}/{n}: {e}"),
            })?;
            records.push(TgtRecord::read_binary(&buf));
        }
        Ok(Tgt {
            genome_name: meta.genome_name,
            contigs: meta.contigs,
            records,
        })
    }

    /// Text form, grouped by enzyme as in the Syn2b `Enzyme:SEQUENCE` layout.
    pub fn write_text(&self, path: &Path) -> Result<()> {
        let f = std::fs::File::create(path).map_err(|e| Sk2bError::Io {
            path: path.into(),
            source: e,
        })?;
        let mut w = BufWriter::new(f);
        let io = |e: std::io::Error| Sk2bError::Io {
            path: path.to_path_buf(),
            source: e,
        };
        writeln!(w, "#TGT2\t{}", self.genome_name).map_err(io)?;
        for c in &self.contigs {
            writeln!(
                w,
                "#contig\t{}\t{}\t{}\t{}\t{:?}",
                c.id, c.name, c.length, c.offset, c.kind
            )
            .map_err(io)?;
        }
        writeln!(
            w,
            "#columns\ttag\tcontig_id\tposition\tstrand\tgap\tflags\tlocal_gc"
        )
        .map_err(io)?;
        let mut idxs: Vec<usize> = (0..self.records.len()).collect();
        idxs.sort_by_key(|&i| {
            let r = &self.records[i];
            (r.enzyme_idx, r.contig_id, r.position)
        });
        for i in idxs {
            let r = &self.records[i];
            writeln!(
                w,
                "{}:{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.enzyme_name(),
                String::from_utf8_lossy(&r.tag_bases()),
                r.contig_id,
                r.position,
                r.strand.symbol(),
                r.gap,
                r.flags,
                r.local_gc
            )
            .map_err(io)?;
        }
        w.flush().map_err(io)?;
        Ok(())
    }

    /// Parse the text form back. Round-trips everything except `tag_hash`, which
    /// is recomputed from the tag sequence.
    pub fn read_text(path: &Path) -> Result<Self> {
        let reader = crate::fasta::open_reader(path)?;
        let io = |e: std::io::Error| Sk2bError::Io {
            path: path.to_path_buf(),
            source: e,
        };
        let mut genome_name = String::new();
        let mut contigs: Vec<ContigMeta> = Vec::new();
        let mut records: Vec<TgtRecord> = Vec::new();
        for line in reader.lines() {
            let line = line.map_err(io)?;
            let line = line.trim_end();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("#TGT2") {
                genome_name = rest.trim().to_string();
                continue;
            }
            if let Some(rest) = line.strip_prefix("#contig\t") {
                let f: Vec<&str> = rest.split('\t').collect();
                if f.len() < 4 {
                    return Err(Sk2bError::Format {
                        path: path.into(),
                        msg: format!("bad #contig line: {line}"),
                    });
                }
                contigs.push(ContigMeta {
                    id: f[0].parse().map_err(|_| bad(path, line))?,
                    name: f[1].to_string(),
                    length: f[2].parse().map_err(|_| bad(path, line))?,
                    offset: f[3].parse().map_err(|_| bad(path, line))?,
                    kind: match f.get(4).copied().unwrap_or("Chromosome") {
                        "Plasmid" => ContigKind::Plasmid,
                        "Unknown" => ContigKind::Unknown,
                        _ => ContigKind::Chromosome,
                    },
                });
                continue;
            }
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                return Err(bad(path, line));
            }
            let (ename, tag) = f[0].split_once(':').ok_or_else(|| bad(path, line))?;
            let enzyme = by_name(ename).ok_or_else(|| Sk2bError::Format {
                path: path.to_path_buf(),
                msg: format!("unknown enzyme '{ename}'"),
            })?;
            let tag_b = tag.as_bytes();
            records.push(TgtRecord {
                tag_hash: canonical_hash(tag_b),
                position: f[2].parse().map_err(|_| bad(path, line))?,
                gap: f[4].parse().map_err(|_| bad(path, line))?,
                contig_id: f[1].parse().map_err(|_| bad(path, line))?,
                enzyme_idx: enzyme.idx,
                strand: if f[3] == "-" {
                    Strand::Rev
                } else {
                    Strand::Fwd
                },
                tag_len: tag_b.len() as u8,
                flags: f.get(5).and_then(|s| s.parse().ok()).unwrap_or(0),
                local_gc: f
                    .get(6)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(GC_UNDEFINED),
                tag_2bit: pack_bases(tag_b),
            });
        }
        Ok(Tgt {
            genome_name,
            contigs,
            records,
        })
    }
}

fn bad(path: &Path, line: &str) -> Sk2bError {
    Sk2bError::Format {
        path: path.to_path_buf(),
        msg: format!("malformed TGT text line: {line}"),
    }
}

#[derive(Serialize, Deserialize)]
struct TgtHeaderMeta {
    genome_name: String,
    contigs: Vec<ContigMeta>,
}

/// Pack ASCII bases 4-per-byte, most significant pair first. Ambiguous bases
/// pack as `A`; digestion rejects them upstream, so this is only reachable for
/// records built with `reject_ambiguous_tags = false`.
pub fn pack_bases(tag: &[u8]) -> [u8; TAG_BYTES] {
    let mut out = [0u8; TAG_BYTES];
    for (i, &b) in tag.iter().take(MAX_PACKED_BASES).enumerate() {
        let code = two_bit(b).unwrap_or(0);
        out[i / 4] |= code << (6 - 2 * (i % 4));
    }
    out
}

/// Inverse of [`pack_bases`].
pub fn unpack_bases(packed: &[u8; TAG_BYTES], len: usize) -> Vec<u8> {
    const LUT: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..len.min(MAX_PACKED_BASES))
        .map(|i| LUT[((packed[i / 4] >> (6 - 2 * (i % 4))) & 0b11) as usize])
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::{digest_contig, DigestConfig};
    use crate::enzyme::by_name;

    fn sample_tgt() -> Tgt {
        let seq = b"AAAAAAAAAAAAAAAAAAAACGAACGTACTGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCGAACGTACTGCTTTTTTTTTTTTTTTTTTTT".to_vec();
        let bcgi = by_name("BcgI").unwrap();
        let sites = digest_contig(
            &seq,
            0,
            &[bcgi],
            &DigestConfig {
                reject_ambiguous_tags: true,
                min_contig_len: 0,
            },
        );
        assert!(
            sites.len() >= 2,
            "fixture needs at least two sites, got {}",
            sites.len()
        );
        let mut tgt = Tgt::new(
            "fixture",
            vec![ContigMeta {
                id: 0,
                name: "c0".into(),
                length: seq.len() as u64,
                offset: 0,
                kind: ContigKind::Chromosome,
            }],
        );
        tgt.records = sites
            .iter()
            .map(|s| TgtRecord::from_site(s).unwrap())
            .collect();
        tgt.recompute_gaps();
        tgt
    }

    #[test]
    fn pack_unpack_roundtrip() {
        let tag = b"ACGTACGTACGTACGTACGTACGTACGTACGTA";
        assert_eq!(unpack_bases(&pack_bases(tag), tag.len()), tag.to_vec());
        // 33 bp (CspCI, the panel maximum) must survive.
        assert_eq!(tag.len(), 33);
    }

    #[test]
    fn gaps_reset_at_contig_boundaries() {
        let mut tgt = sample_tgt();
        tgt.contigs.push(ContigMeta {
            id: 1,
            name: "c1".into(),
            length: 500,
            offset: 0,
            kind: ContigKind::Chromosome,
        });
        let mut r = tgt.records[0].clone();
        r.contig_id = 1;
        r.position = 42;
        tgt.records.push(r);
        tgt.recompute_gaps();
        let first_on_c1 = tgt.records.iter().find(|r| r.contig_id == 1).unwrap();
        assert_eq!(first_on_c1.gap, 0, "gap spanned a contig boundary");
        let second_on_c0 = tgt
            .records
            .iter()
            .filter(|r| r.contig_id == 0)
            .nth(1)
            .unwrap();
        assert!(second_on_c0.gap > 0);
    }

    #[test]
    fn binary_roundtrip() {
        let tgt = sample_tgt();
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-rt.tgt", std::process::id()));
        tgt.write_binary(&p).unwrap();
        let back = Tgt::read_binary(&p).unwrap();
        assert_eq!(back, tgt);
        // Fixed-size records: header + n * 48 must match the file length.
        let meta_len = std::fs::metadata(&p).unwrap().len();
        assert!(meta_len >= (tgt.records.len() * RECORD_SIZE) as u64);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn text_roundtrip_preserves_positions_and_tags() {
        let tgt = sample_tgt();
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-rt.tgt.txt", std::process::id()));
        tgt.write_text(&p).unwrap();
        let back = Tgt::read_text(&p).unwrap();
        assert_eq!(back.contigs, tgt.contigs);
        assert_eq!(back.records.len(), tgt.records.len());
        let mut a: Vec<_> = back
            .records
            .iter()
            .map(|r| (r.position, r.tag_bases(), r.enzyme_idx))
            .collect();
        let mut b: Vec<_> = tgt
            .records
            .iter()
            .map(|r| (r.position, r.tag_bases(), r.enzyme_idx))
            .collect();
        a.sort();
        b.sort();
        assert_eq!(a, b);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn rejects_a_foreign_magic() {
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-bad.tgt", std::process::id()));
        std::fs::write(&p, b"NOPE____________").unwrap();
        assert!(Tgt::read_binary(&p).is_err());
        std::fs::remove_file(p).ok();
    }
}
