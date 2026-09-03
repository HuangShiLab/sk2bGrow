//! `sk2bgrow index` — build an anchor database from reference genomes.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use rayon::prelude::*;
use sk2bgrow_core::anchor_db::{
    assemble, attach_taxonomy, build_genome, tgt_path, BuildParams, GC_FLANK,
};
use sk2bgrow_core::digest::DigestConfig;
use sk2bgrow_core::enzyme::{parse_selection, EnzymeSet};
use sk2bgrow_core::ori;

use super::{expand_inputs, Ctx, GENOME_EXTS};

#[derive(ClapArgs)]
pub struct Args {
    /// Reference genome FASTA files, or directories containing them.
    #[arg(required = true)]
    pub genomes: Vec<PathBuf>,

    /// Output database directory.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Enzyme selection: `all`, or a comma-separated list (e.g. `BcgI,AlfI`).
    #[arg(short = 'e', long, default_value = "all")]
    pub enzymes: String,

    /// Two-column TSV of `genome_name<TAB>lineage`.
    #[arg(short = 'a', long)]
    pub taxonomy: Option<PathBuf>,

    /// Ori annotation TSV: `genome<TAB>position[<TAB>confidence][<TAB>source]`.
    #[arg(long)]
    pub ori: Option<PathBuf>,

    /// Skip contigs shorter than this many bp.
    #[arg(long, default_value_t = 500)]
    pub min_contig_len: usize,

    /// Half-width of the local GC window, in bp.
    #[arg(long, default_value_t = GC_FLANK)]
    pub gc_flank: usize,

    /// Also write a per-genome TGT text dump under `<output>/tgt/`.
    #[arg(long)]
    pub write_tgt: bool,

    /// Also build a FracMinHash containment sketch for `profile --screen`
    /// with this sampling scale (a k-mer hash is kept when
    /// `h < u64::MAX / scale`; larger is sparser). The sketch lands next to
    /// the anchor database as `screen.meta` + `screen.csr` and roughly
    /// doubles the disk footprint of the genome sequences.
    #[arg(long)]
    pub screen_scale: Option<u64>,
}

pub fn run(args: Args, ctx: &Ctx) -> Result<()> {
    let genomes = expand_inputs(&args.genomes, GENOME_EXTS)?;
    let enzymes = parse_selection(&args.enzymes)?;
    ctx.say(format!(
        "indexing {} genome(s) with {} enzyme(s)",
        genomes.len(),
        enzymes.len()
    ));

    let cfg = DigestConfig {
        reject_ambiguous_tags: true,
        min_contig_len: args.min_contig_len,
    };
    let gc_flank = args.gc_flank;
    let write_tgt = args.write_tgt;
    let out = args.output.clone();

    // Genomes are independent, so digestion parallelises cleanly. Genome ids are
    // assigned from the sorted input order, not from completion order, so a
    // rebuild of the same inputs yields byte-identical ids.
    let parts: Vec<_> = genomes
        .par_iter()
        .enumerate()
        .map(|(i, path)| -> Result<_> {
            let (meta, anchors, tags, tgt) = build_genome(path, i as u32, &enzymes, &cfg, gc_flank)
                .with_context(|| format!("digesting {}", path.display()))?;
            if write_tgt {
                let p = tgt_path(&out, &meta.name);
                if let Some(dir) = p.parent() {
                    std::fs::create_dir_all(dir)?;
                }
                tgt.write_text(&p)?;
            }
            Ok((meta, anchors, tags))
        })
        .collect::<Result<Vec<_>>>()?;

    for (meta, anchors, _) in &parts {
        if anchors.is_empty() {
            ctx.say(format!("  warning: {} yielded no anchors", meta.name));
        }
    }

    let params = BuildParams {
        enzymes: EnzymeSet::from_slice(&enzymes),
        gc_flank,
        min_contig_len: args.min_contig_len,
        reject_ambiguous_tags: true,
        sk2bgrow_version: sk2bgrow_core::VERSION.to_string(),
    };
    let mut db = assemble(params, parts);
    ctx.say(format!(
        "{} anchors across {} genomes",
        db.n_anchors(),
        db.genomes.len()
    ));

    if let Some(tax) = &args.taxonomy {
        let hit = attach_taxonomy(&mut db, tax)?;
        ctx.say(format!(
            "taxonomy attached to {}/{} genomes",
            hit,
            db.genomes.len()
        ));
    }
    if let Some(op) = &args.ori {
        let ann = ori::read_annotations(op)?;
        let hit = ori::attach(&mut db, &ann);
        ctx.say(format!(
            "ori annotated for {}/{} genomes",
            hit,
            db.genomes.len()
        ));
    }

    let usable = db.anchors.iter().filter(|a| a.is_usable()).count();
    ctx.say(format!(
        "{usable} usable anchors ({:.1}% masked as multi-copy, shared or non-chromosomal)",
        100.0 * (1.0 - usable as f64 / db.n_anchors().max(1) as f64)
    ));

    db.save(&args.output)
        .with_context(|| format!("writing database to {}", args.output.display()))?;
    ctx.say(format!("database written to {}", args.output.display()));

    if let Some(scale) = args.screen_scale {
        let t = std::time::Instant::now();
        // Genome ids were assigned from this same sorted path list, so the
        // sketch's per-genome indexing matches the anchor database's.
        let sketch = sk2bgrow_core::screen::ScreenSketch::build(
            &genomes,
            sk2bgrow_core::screen::SCREEN_K,
            scale,
        )?;
        sketch.save(&args.output)?;
        ctx.say(format!(
            "containment sketch written (scale {scale}, {} genomes, {:.1}s)",
            sketch.sizes.len(),
            t.elapsed().as_secs_f64()
        ));
    }
    Ok(())
}
