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
use sk2bgrow_core::scaffold::{apply, scaffold, ScaffoldConfig};

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
        "wrote {} and {}",
        args.output.display(),
        json.display()
    ));
    Ok(())
}
