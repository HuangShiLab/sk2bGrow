//! `sk2bgrow digest` — digest genomes and print the density table.
//!
//! This is the command that reproduces report §4.1: per-enzyme anchor density
//! per Mb, the 16-enzyme union count, mean spacing and the worst-case gap. It
//! needs no database, so it is the fastest way to check a new reference before
//! committing it to an index.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use rayon::prelude::*;
use sk2bgrow_core::digest::{density_report, digest_contig, DigestConfig, DEFAULT_MERGE_WINDOW};
use sk2bgrow_core::enzyme::parse_selection;
use sk2bgrow_core::fasta;

use super::{expand_inputs, Ctx, GENOME_EXTS};

#[derive(ClapArgs)]
pub struct Args {
    /// Genome FASTA files or directories.
    #[arg(required = true)]
    pub genomes: Vec<PathBuf>,

    #[arg(short = 'e', long, default_value = "all")]
    pub enzymes: String,

    /// Write TSV here instead of stdout.
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Merge radius for the cross-enzyme union, in bp.
    #[arg(long, default_value_t = DEFAULT_MERGE_WINDOW)]
    pub merge_window: u64,

    #[arg(long, default_value_t = 500)]
    pub min_contig_len: usize,
}

pub fn run(args: Args, ctx: &Ctx) -> Result<()> {
    let genomes = expand_inputs(&args.genomes, GENOME_EXTS)?;
    let enzymes = parse_selection(&args.enzymes)?;
    let cfg = DigestConfig {
        reject_ambiguous_tags: true,
        min_contig_len: args.min_contig_len,
    };

    let rows: Vec<(String, sk2bgrow_core::digest::DensityReport)> = genomes
        .par_iter()
        .map(|path| -> Result<_> {
            let records =
                fasta::read_fasta(path).with_context(|| format!("reading {}", path.display()))?;
            let mut sites = Vec::new();
            let mut lens = Vec::new();
            for (i, rec) in records.iter().enumerate() {
                lens.push(rec.seq.len() as u64);
                sites.extend(digest_contig(&rec.seq, i as u16, &enzymes, &cfg));
            }
            let rep = density_report(&sites, &lens, &enzymes, args.merge_window);
            let name = path
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_default();
            Ok((name, rep))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut out = String::new();
    out.push_str("genome\tenzyme\tsites\tdensity_per_mb\n");
    for (name, rep) in &rows {
        let mb = rep.genome_len as f64 / 1e6;
        for (enzyme, n) in &rep.per_enzyme {
            out.push_str(&format!("{name}\t{enzyme}\t{n}\t{:.1}\n", *n as f64 / mb));
        }
        out.push_str(&format!(
            "{name}\tUNION\t{}\t{:.1}\n",
            rep.union_sites,
            rep.union_sites as f64 / mb
        ));
    }
    out.push_str("\n# summary\ngenome\tgenome_len\tunion_sites\tper_25kb\tmean_spacing\tmax_gap\n");
    for (name, rep) in &rows {
        out.push_str(&format!(
            "{name}\t{}\t{}\t{:.1}\t{:.1}\t{}\n",
            rep.genome_len, rep.union_sites, rep.per_25kb, rep.mean_spacing, rep.max_gap
        ));
    }

    match &args.output {
        Some(p) => {
            std::fs::write(p, &out)?;
            ctx.say(format!("wrote {}", p.display()));
        }
        None => print!("{out}"),
    }
    Ok(())
}
