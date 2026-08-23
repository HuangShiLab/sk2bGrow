//! Subcommand implementations.

pub mod audit;
pub mod digest;
pub mod dynamics;
pub mod index;
pub mod profile;
pub mod scaffold;

/// Shared run context.
pub struct Ctx {
    pub quiet: bool,
}

impl Ctx {
    /// Progress line on stderr, so stdout stays a clean data stream.
    pub fn say(&self, msg: impl AsRef<str>) {
        if !self.quiet {
            eprintln!("[sk2bgrow] {}", msg.as_ref());
        }
    }
}

/// Expand shell-style inputs: a directory becomes every FASTA/FASTQ it holds.
///
/// Shells normally expand globs, but a large `genomes/*.fna` can blow past
/// `ARG_MAX`, so passing the directory has to work too.
pub fn expand_inputs(
    paths: &[std::path::PathBuf],
    exts: &[&str],
) -> anyhow::Result<Vec<std::path::PathBuf>> {
    let mut out = Vec::new();
    for p in paths {
        if p.is_dir() {
            let mut kids: Vec<std::path::PathBuf> = std::fs::read_dir(p)?
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|c| c.is_file() && has_ext(c, exts))
                .collect();
            kids.sort();
            if kids.is_empty() {
                anyhow::bail!("{} contains no files matching {:?}", p.display(), exts);
            }
            out.extend(kids);
        } else {
            out.push(p.clone());
        }
    }
    if out.is_empty() {
        anyhow::bail!("no input files");
    }
    for p in &out {
        if !p.exists() {
            anyhow::bail!("input does not exist: {}", p.display());
        }
    }
    Ok(out)
}

fn has_ext(p: &std::path::Path, exts: &[&str]) -> bool {
    let name = p
        .file_name()
        .map(|n| n.to_string_lossy().to_lowercase())
        .unwrap_or_default();
    exts.iter().any(|e| name.ends_with(e))
}

/// Extensions accepted for reference genomes.
pub const GENOME_EXTS: &[&str] = &[".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"];
/// Extensions accepted for reads.
pub const READ_EXTS: &[&str] = &[
    ".fq",
    ".fastq",
    ".fq.gz",
    ".fastq.gz",
    ".fa",
    ".fasta",
    ".fa.gz",
    ".fasta.gz",
];

/// Strip read-file suffixes to get a sample name.
pub fn sample_name(path: &std::path::Path) -> String {
    let mut n = path
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();
    for suffix in [".gz", ".fastq", ".fq", ".fasta", ".fa", ".fna"] {
        if let Some(stripped) = n.strip_suffix(suffix) {
            n = stripped.to_string();
        }
    }
    // Drop a trailing _R1 / _1 mate marker so paired files share a sample name.
    for suffix in ["_R1", "_R2", "_1", "_2"] {
        if let Some(stripped) = n.strip_suffix(suffix) {
            n = stripped.to_string();
            break;
        }
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn sample_names_drop_suffixes_and_mate_markers() {
        assert_eq!(
            sample_name(&PathBuf::from("a/b/SRR123_R1.fastq.gz")),
            "SRR123"
        );
        assert_eq!(sample_name(&PathBuf::from("SRR123_2.fq")), "SRR123");
        assert_eq!(sample_name(&PathBuf::from("sample.fa")), "sample");
        // A name that merely ends in a digit must survive intact.
        assert_eq!(sample_name(&PathBuf::from("day3.fastq")), "day3");
    }
}
