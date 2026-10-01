//! `sk2bgrow scaffold` — order and orient draft MAG contigs against a reference.
//!
//! Wraps [`sk2bgrow_core::scaffold`]. Run this before `profile` when the
//! reference is a fragmented MAG: without a contig order there is no genomic
//! coordinate, and without a coordinate the V-shape fit has no x-axis.

use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use sk2bgrow_core::anchor_db::{build_genome, AnchorDb, GC_FLANK};
use sk2bgrow_core::digest::DigestConfig;
use sk2bgrow_core::enzyme::parse_selection;
use sk2bgrow_core::fasta;
use sk2bgrow_core::scaffold::{apply, scaffold, ScaffoldConfig};
use sk2bgrow_core::seq::revcomp;

use super::Ctx;

#[derive(ClapArgs)]
pub struct Args {
    /// Draft MAG FASTA.
    pub draft: PathBuf,

    /// Database holding the reference genome.
    #[arg(short, long)]
    pub db: PathBuf,

    /// Reference genome name inside the database.
    #[arg(short, long)]
    pub reference: String,

    /// Output TGT file for the scaffolded draft.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Single-contig FASTA suitable for `sk2bgrow index` and the default
    /// coordinate fit. If omitted, it is written next to `output` as
    /// `<output-stem>.scaffolded.fna`.
    #[arg(long)]
    pub index_fasta: Option<PathBuf>,

    #[arg(short = 'e', long, default_value = "all")]
    pub enzymes: String,

    /// Minimum shared tags before a contig may be placed.
    #[arg(long, default_value_t = 3)]
    pub min_tags: usize,

    /// Minimum ordering agreement for a placement to be accepted.
    #[arg(long, default_value_t = 0.8)]
    pub min_concordance: f64,
}

pub fn run(args: Args, ctx: &Ctx) -> Result<()> {
    let db = AnchorDb::load(&args.db)
        .with_context(|| format!("loading database {}", args.db.display()))?;
    let Some(reference) = db.genomes.iter().find(|g| g.name == args.reference) else {
        bail!(
            "reference '{}' is not in the database; it holds: {}",
            args.reference,
            db.genomes
                .iter()
                .map(|g| g.name.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
    };
    let enzymes = parse_selection(&args.enzymes)?;
    let (_, _, _, mut draft) = build_genome(
        &args.draft,
        u32::MAX,
        &enzymes,
        &DigestConfig::default(),
        GC_FLANK,
    )
    .with_context(|| format!("digesting {}", args.draft.display()))?;
    ctx.say(format!(
        "draft: {} contigs, {} tags",
        draft.contigs.len(),
        draft.records.len()
    ));

    let cfg = ScaffoldConfig {
        min_tags: args.min_tags,
        min_concordance: args.min_concordance,
    };
    let result = scaffold(&draft, &db, reference.id, &cfg);
    apply(&mut draft, &result);

    let reversed = result
        .placements
        .iter()
        .filter(|p| p.orientation == sk2bgrow_core::scaffold::Orientation::Reverse)
        .count();
    ctx.say(format!(
        "placed {}/{} contigs ({:.1}% of draft bp), {} reversed, {} unplaced",
        result.placements.len(),
        draft.contigs.len(),
        100.0 * result.placed_fraction(&draft),
        reversed,
        result.unplaced.len()
    ));

    draft.write_text(&args.output)?;
    let fasta_out = args.index_fasta.clone().unwrap_or_else(|| {
        let stem = args
            .output
            .file_stem()
            .map(|s| format!("{}.scaffolded.fna", s.to_string_lossy()))
            .unwrap_or_else(|| "scaffolded.fna".to_string());
        args.output.with_file_name(stem)
    });
    write_index_fasta(&args.draft, &result, &fasta_out)?;
    let json = args.output.with_extension("scaffold.json");
    std::fs::write(
        &json,
        serde_json::to_vec_pretty(&serde_json::json!({
            "draft": args.draft,
            "reference": args.reference,
            "placements": result.placements,
            "unplaced": result.unplaced,
        }))?,
    )?;
    ctx.say(format!(
        "wrote {}, {}, and {}",
        args.output.display(),
        fasta_out.display(),
        json.display()
    ));
    Ok(())
}

/// Emit placed contigs as one pseudo-chromosome in placement order.
///
/// Concatenation intentionally uses order-preserving offsets rather than the
/// original reference starts: overlapping placements therefore cannot overwrite
/// one another. Rotation and stretch are harmless to the V-fit because it
/// searches for the origin; order and orientation are what it needs.
fn write_index_fasta(
    draft: &std::path::Path,
    result: &sk2bgrow_core::scaffold::ScaffoldResult,
    output: &std::path::Path,
) -> Result<()> {
    let records = fasta::read_fasta(draft)?;
    let sequences: Vec<Vec<u8>> = records.iter().map(|record| record.seq.clone()).collect();
    let mut chunks: Vec<Vec<u8>> =
        Vec::with_capacity(result.placements.len() + result.unplaced.len());
    for placement in &result.placements {
        let seq = sequences.get(placement.contig_id as usize).ok_or_else(|| {
            anyhow::anyhow!(
                "scaffold refers to missing contig id {}",
                placement.contig_id
            )
        })?;
        chunks.push(
            if placement.orientation == sk2bgrow_core::scaffold::Orientation::Reverse {
                revcomp(seq)
            } else {
                seq.clone()
            },
        );
    }
    for id in &result.unplaced {
        let seq = sequences
            .get(*id as usize)
            .ok_or_else(|| anyhow::anyhow!("scaffold refers to missing contig id {id}"))?;
        chunks.push(seq.clone());
    }
    let mut joined = Vec::with_capacity(chunks.iter().map(Vec::len).sum());
    for chunk in chunks {
        joined.extend_from_slice(&chunk);
    }
    let mut body = String::from(">scaffolded\n");
    for chunk in joined.chunks(70) {
        body.push_str(&String::from_utf8_lossy(chunk));
        body.push('\n');
    }
    std::fs::write(output, body)?;
    Ok(())
}
