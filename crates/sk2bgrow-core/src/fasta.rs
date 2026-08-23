//! Minimal streaming FASTA/FASTQ reader with transparent gzip support.
//!
//! Deliberately dependency-light: the counter only needs sequence bytes, and a
//! hand-rolled reader keeps the anchor pipeline free of a parser dependency that
//! would have to be reconciled with Syn2b's.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::{Result, Sk2bError};

/// One record: name (up to the first whitespace) plus sequence.
#[derive(Debug, Clone)]
pub struct Record {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Open a path as a buffered reader, transparently decompressing `.gz`.
///
/// `MultiGzDecoder` rather than `GzDecoder`: concatenated gzip members are the
/// norm for bgzip-style read files, and the single-member decoder silently stops
/// at the first member boundary — a truncation bug that shows up as "half my
/// reads vanished".
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| Sk2bError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    let is_gz = path.extension().is_some_and(|e| e == "gz");
    let inner: Box<dyn Read> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::with_capacity(1 << 20, inner)))
}

/// Read an entire FASTA file into memory.
///
/// Reference genomes are read whole because digestion needs random access to the
/// flanks around every site; reads are streamed instead (see [`for_each_read`]).
pub fn read_fasta(path: &Path) -> Result<Vec<Record>> {
    let reader = open_reader(path)?;
    let mut records: Vec<Record> = Vec::new();
    let mut cur: Option<Record> = None;
    for line in reader.lines() {
        let line = line.map_err(|e| Sk2bError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(rec) = cur.take() {
                records.push(rec);
            }
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            cur = Some(Record {
                name,
                seq: Vec::new(),
            });
        } else if let Some(rec) = cur.as_mut() {
            rec.seq.extend_from_slice(line.as_bytes());
        } else {
            return Err(Sk2bError::Format {
                path: path.to_path_buf(),
                msg: "sequence data before the first '>' header".into(),
            });
        }
    }
    if let Some(rec) = cur.take() {
        records.push(rec);
    }
    for rec in records.iter_mut() {
        crate::seq::normalize_in_place(&mut rec.seq);
    }
    if records.is_empty() {
        return Err(Sk2bError::Format {
            path: path.to_path_buf(),
            msg: "no FASTA records".into(),
        });
    }
    Ok(records)
}

/// Stream sequences from a FASTA or FASTQ file, calling `f` once per record.
///
/// Format is detected from the first non-empty byte, so `--mode 2brad` inputs
/// (FASTQ) and in silico inputs (FASTA) share one code path. Sequences are
/// normalised to uppercase ACGTN before the callback sees them.
pub fn for_each_read<F>(path: &Path, mut f: F) -> Result<u64>
where
    F: FnMut(&[u8]),
{
    let mut reader = open_reader(path)?;
    let mut first = [0u8; 1];
    let mut n = 0u64;
    // Peek without consuming: fill_buf keeps the byte in the buffer.
    {
        let buf = reader.fill_buf().map_err(|e| Sk2bError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
        if buf.is_empty() {
            return Ok(0);
        }
        first[0] = buf[0];
    }
    let io = |e: std::io::Error| Sk2bError::Io {
        path: path.to_path_buf(),
        source: e,
    };

    match first[0] {
        b'>' => {
            let mut seq: Vec<u8> = Vec::new();
            let mut have = false;
            for line in reader.lines() {
                let line = line.map_err(io)?;
                let line = line.trim_end();
                if line.starts_with('>') {
                    if have && !seq.is_empty() {
                        crate::seq::normalize_in_place(&mut seq);
                        f(&seq);
                        n += 1;
                    }
                    seq.clear();
                    have = true;
                } else {
                    seq.extend_from_slice(line.as_bytes());
                }
            }
            if have && !seq.is_empty() {
                crate::seq::normalize_in_place(&mut seq);
                f(&seq);
                n += 1;
            }
        }
        b'@' => {
            let mut line = String::new();
            let mut which = 0u8;
            loop {
                line.clear();
                let got = reader.read_line(&mut line).map_err(io)?;
                if got == 0 {
                    break;
                }
                if which == 1 {
                    let mut seq = line.trim_end().as_bytes().to_vec();
                    crate::seq::normalize_in_place(&mut seq);
                    f(&seq);
                    n += 1;
                }
                which = (which + 1) % 4;
            }
            if which != 0 {
                return Err(Sk2bError::Format {
                    path: path.to_path_buf(),
                    msg: "truncated FASTQ: record count is not a multiple of 4 lines".into(),
                });
            }
        }
        other => {
            return Err(Sk2bError::Format {
                path: path.to_path_buf(),
                msg: format!("expected FASTA '>' or FASTQ '@', found {:?}", other as char),
            })
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp(name: &str, body: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-{}", std::process::id(), name));
        let mut f = File::create(&p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    #[test]
    fn reads_multiline_fasta() {
        let p = tmp("a.fna", ">c1 some description\nACGT\nacgt\n>c2\nTTTT\n");
        let recs = read_fasta(&p).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(
            recs[0].name, "c1",
            "description leaked into the contig name"
        );
        assert_eq!(recs[0].seq, b"ACGTACGT".to_vec(), "soft-masked bases lost");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn non_acgt_becomes_n() {
        let p = tmp("b.fna", ">c\nACGTRYKM\n");
        let recs = read_fasta(&p).unwrap();
        assert_eq!(recs[0].seq, b"ACGTNNNN".to_vec());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn streams_fastq_sequences_only() {
        let p = tmp("c.fq", "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n");
        let mut seen: Vec<Vec<u8>> = Vec::new();
        let n = for_each_read(&p, |s| seen.push(s.to_vec())).unwrap();
        assert_eq!(n, 2);
        assert_eq!(seen, vec![b"ACGT".to_vec(), b"TTTT".to_vec()]);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn truncated_fastq_is_an_error() {
        let p = tmp("d.fq", "@r1\nACGT\n+\n");
        assert!(for_each_read(&p, |_| {}).is_err());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn streams_fasta_reads_too() {
        let p = tmp("e.fa", ">r1\nACGT\n>r2\nTT\nTT\n");
        let mut seen = 0;
        let n = for_each_read(&p, |_| seen += 1).unwrap();
        assert_eq!((n, seen), (2, 2));
        std::fs::remove_file(p).ok();
    }
}
