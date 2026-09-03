//! `sk2bgrow profile` — count reads against a database, then hand off to the
//! Python statistics layer.
//!
//! The Rust half stops at a per-anchor count table. Everything downstream —
//! ZTP/NB window rates, GC correction, the V-shape fit, cross-enzyme fusion — is
//! `python/sk2bgrow`, invoked as a subprocess. The interface is the file, not an
//! FFI boundary, which is what lets the two layers iterate independently.

use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Args as ClapArgs, ValueEnum};
use rayon::prelude::*;
use sk2bgrow_core::anchor_db::AnchorDb;
use sk2bgrow_core::count::{count_sample, write_count_table, AnchorIndex, CountMode, MatchConfig};
use sk2bgrow_core::em::{reassign, EmConfig};
use sk2bgrow_core::window::{assign_windows, WindowPolicy, NO_WINDOW};

use super::{expand_inputs, sample_name, Ctx, READ_EXTS};

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum Mode {
    /// Shotgun metagenome reads (route A).
    Wms,
    /// Real 2bRAD reads (route B): one tag per read, digestion already done at
    /// the bench.
    #[value(name = "2brad")]
    TwoBrad,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum Windowing {
    /// Fixed anchor count per window (the TGT-native policy).
    Anchors,
    /// Fixed base pairs per window (Pilea parity, for A/B benchmarks).
    Bp,
}

#[derive(ClapArgs)]
pub struct Args {
    /// Read files (FASTQ/FASTA, optionally gzipped), or directories of them.
    #[arg(required = true)]
    pub reads: Vec<PathBuf>,

    /// Anchor database directory built by `sk2bgrow index`.
    #[arg(short, long)]
    pub db: PathBuf,

    /// Output directory.
    #[arg(short, long)]
    pub output: PathBuf,

    #[arg(long, value_enum, default_value_t = Mode::Wms)]
    pub mode: Mode,

    /// Restrict counting to a subset of the database's enzymes.
    #[arg(short = 'e', long)]
    pub enzymes: Option<String>,

    /// Maximum Hamming distance between read tag and reference tag.
    #[arg(long, default_value_t = 2)]
    pub max_mismatch: u32,

    /// Treat every read file as its own sample instead of grouping by name.
    #[arg(long)]
    pub per_file: bool,

    #[arg(long, value_enum, default_value_t = Windowing::Anchors)]
    pub windowing: Windowing,

    /// Anchors per window (with `--windowing anchors`).
    #[arg(long, default_value_t = 100)]
    pub window_anchors: usize,

    /// Base pairs per window (with `--windowing bp`).
    #[arg(long, default_value_t = 25_000)]
    pub window_bp: u64,

    /// Stop after writing count tables; do not run the Python statistics layer.
    #[arg(long)]
    pub no_stats: bool,

    /// Two-tier counting (M4): pass 1 screens reads against the per-genome
    /// containment sketch built by `index --screen-scale`, pass 2 counts
    /// anchors only on genomes the screen selected. Requires a sketch in the
    /// database directory.
    #[arg(long)]
    pub screen: bool,

    /// Minimum estimated relative abundance for a genome to enter pass 2
    /// (with `--screen`). Genomes with at least 50 sketch hits are kept
    /// regardless.
    #[arg(long, default_value_t = 5e-4)]
    pub screen_min_frac: f64,

    /// Python interpreter used for the statistics layer.
    #[arg(long, default_value = "python3")]
    pub python: String,
}

pub fn run(args: Args, ctx: &Ctx) -> Result<()> {
    let read_files = expand_inputs(&args.reads, READ_EXTS)?;
    let db = AnchorDb::load(&args.db)
        .with_context(|| format!("loading database {}", args.db.display()))?;
    ctx.say(format!(
        "database: {} anchors, {} genomes",
        db.n_anchors(),
        db.genomes.len()
    ));

    // `--enzymes` restricts the *counting*, not just the report: excluded
    // anchors are left out of the lookup tables, so a restricted run costs what
    // a database built with only those enzymes would cost.
    let restrict = match &args.enzymes {
        None => None,
        Some(sel) => {
            let wanted = sk2bgrow_core::enzyme::parse_selection(sel)?;
            let missing: Vec<&str> = wanted
                .iter()
                .filter(|e| !db.params.enzymes.contains(e.idx))
                .map(|e| e.name)
                .collect();
            if !missing.is_empty() {
                bail!(
                    "the database was not built with {}; rebuild the index or drop them from --enzymes",
                    missing.join(", ")
                );
            }
            let set = sk2bgrow_core::enzyme::EnzymeSet::from_slice(&wanted);
            ctx.say(format!(
                "restricting to {} of {} enzymes",
                wanted.len(),
                db.params.enzymes.len()
            ));
            Some(set)
        }
    };

    std::fs::create_dir_all(&args.output)?;

    // Group read files into samples so paired mates land in one count table.
    let mut samples: Vec<(String, Vec<PathBuf>)> = Vec::new();
    for f in &read_files {
        let name = if args.per_file {
            f.file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_default()
        } else {
            sample_name(f)
        };
        match samples.iter_mut().find(|(n, _)| *n == name) {
            Some((_, v)) => v.push(f.clone()),
            None => samples.push((name, vec![f.clone()])),
        }
    }
    ctx.say(format!(
        "{} read file(s) -> {} sample(s)",
        read_files.len(),
        samples.len()
    ));

    let policy = match args.windowing {
        Windowing::Anchors => WindowPolicy::EqualAnchors {
            n: args.window_anchors,
        },
        Windowing::Bp => WindowPolicy::FixedBp { bp: args.window_bp },
    };

    // Window assignment depends only on the database, so compute it once and
    // reuse it for every sample — that also guarantees window ids are comparable
    // across samples, which `sk2bgrow dynamics` relies on.
    let mut window_ids = vec![NO_WINDOW; db.n_anchors()];
    let mut all_windows = Vec::new();
    for g in &db.genomes {
        let (windows, ids) = assign_windows(&db, g.id, policy, true);
        let range = db.genome_range(g.id);
        let base = all_windows.len() as u32;
        for (k, &w) in ids.iter().enumerate() {
            window_ids[range.start + k] = if w == NO_WINDOW { NO_WINDOW } else { w + base };
        }
        all_windows.extend(windows.into_iter().map(|mut w| {
            w.id += base;
            w
        }));
    }
    ctx.say(format!(
        "{} windows ({})",
        all_windows.len(),
        describe(policy)
    ));
    let wpath = args.output.join("windows.tsv");
    write_windows(&wpath, &db, &all_windows)?;

    let index = AnchorIndex::build_restricted(&db, args.max_mismatch, restrict);
    let cfg = MatchConfig {
        max_mismatch: args.max_mismatch,
        mode: match args.mode {
            Mode::Wms => CountMode::Wms,
            Mode::TwoBrad => CountMode::TwoBrad,
        },
        // Shared anchors must keep their counts: the EM below is what resolves them.
        keep_multimappers: true,
    };

    // Pass 1 of the two-tier architecture (M4): the containment sketch picks
    // each sample's genome subset. Without the sketch files --screen is an
    // error rather than a silent full scan — a screen the user asked for and
    // did not get would silently cost the full-index price.
    let screen_sketch = if args.screen {
        let sk = sk2bgrow_core::screen::ScreenSketch::load(&args.db)?.ok_or_else(|| {
            anyhow::anyhow!(
                "--screen needs a containment sketch, but {} has none; rebuild the index with `sk2bgrow index --screen-scale <N>`",
                args.db.display()
            )
        })?;
        ctx.say(format!(
            "screen: k={}, scale={}, {} genomes sketched",
            sk.k,
            sk.scale,
            sk.sizes.len()
        ));
        Some(sk)
    } else {
        None
    };

    let outputs: Vec<PathBuf> = samples
        .par_iter()
        .map(|(name, files)| -> Result<PathBuf> {
            // With --screen the genome subset is per-sample, so pass 2's
            // index is built per sample inside the loop; the merge and the
            // count tables stay identical to the unscreened path (the counts
            // vector is always full-length — unselected genomes are zeros).
            let (counts, stats, screen_json) = if let Some(sk) = &screen_sketch {
                let t1 = std::time::Instant::now();
                let hits = sk2bgrow_core::screen::screen_sample(sk, files)?;
                let selected = sk.select(&db, &hits, args.screen_min_frac);
                let pass1 = t1.elapsed();
                let mut filter = vec![false; db.genomes.len()];
                for &g in &selected {
                    filter[g as usize] = true;
                }
                let index =
                    AnchorIndex::build_screened(&db, args.max_mismatch, restrict, Some(&filter));
                let t2 = std::time::Instant::now();
                let (counts, stats) = count_sample(&index, files, &cfg)?;
                let pass2 = t2.elapsed();
                let genomes: Vec<_> = selected
                    .iter()
                    .map(|&g| {
                        let meta = db.genome(g);
                        serde_json::json!({
                            "genome_id": g,
                            "genome": meta.map(|m| m.name.as_str()).unwrap_or("?"),
                            "hits": hits.get(g as usize).copied().unwrap_or(0),
                            "sketch_size": sk.sizes.get(g as usize).copied().unwrap_or(0),
                        })
                    })
                    .collect();
                (
                    counts,
                    stats,
                    serde_json::json!({
                        "k": sk.k,
                        "scale": sk.scale,
                        "min_frac": args.screen_min_frac,
                        "selected_genomes": selected.len(),
                        "pass1_seconds": pass1.as_secs_f64(),
                        "pass2_seconds": pass2.as_secs_f64(),
                        "genomes": genomes,
                    }),
                )
            } else {
                let (counts, stats) = count_sample(&index, files, &cfg)?;
                (counts, stats, serde_json::Value::Null)
            };
            let em = reassign(&db, &counts, &EmConfig::default());
            let tsv = args.output.join(format!("{name}.counts.tsv"));
            write_count_table(&tsv, name, &db, &counts, &window_ids, true, restrict)?;
            let statj = args.output.join(format!("{name}.stats.json"));
            std::fs::write(
                &statj,
                serde_json::to_vec_pretty(&serde_json::json!({
                    "sample": name,
                    "files": files,
                    "mode": match args.mode { Mode::Wms => "wms", Mode::TwoBrad => "2brad" },
                    "max_mismatch": args.max_mismatch,
                    "counting": stats,
                    "screen": screen_json,
                    "em": {
                        "iterations": em.iterations,
                        "converged": em.converged,
                        "genomes": em.genomes,
                    },
                }))?,
            )?;
            Ok(tsv)
        })
        .collect::<Result<Vec<_>>>()?;

    for (name, _) in &samples {
        let statj = args.output.join(format!("{name}.stats.json"));
        ctx.say(format!("  {name}: counts + {}", statj.display()));
    }

    if args.no_stats {
        ctx.say("--no-stats: stopping after the count tables");
        return Ok(());
    }
    run_python_stats(&args, ctx, &outputs)
}

fn describe(p: WindowPolicy) -> String {
    match p {
        WindowPolicy::EqualAnchors { n } => format!("{n} anchors/window"),
        WindowPolicy::FixedBp { bp } => format!("{bp} bp/window"),
    }
}

fn write_windows(
    path: &std::path::Path,
    db: &AnchorDb,
    windows: &[sk2bgrow_core::window::Window],
) -> Result<()> {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "window_id\tgenome_id\tgenome\tcontig_id\tstart\tend\tglobal_mid\tspan\tn_anchors"
    )?;
    for win in windows {
        let gname = db
            .genome(win.genome_id)
            .map(|g| g.name.as_str())
            .unwrap_or("?");
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            win.id,
            win.genome_id,
            gname,
            win.contig_id,
            win.start,
            win.end,
            win.global_mid,
            win.span(),
            win.n_anchors
        )?;
    }
    Ok(())
}

/// Hand the count tables to `python -m sk2bgrow.cli profile`.
///
/// A missing or broken Python layer is reported as an error naming the exact
/// command, rather than silently producing counts and no PTR — "the run
/// succeeded but there is no output.tsv" is the worst possible failure mode for
/// a batch job.
fn run_python_stats(args: &Args, ctx: &Ctx, count_tables: &[PathBuf]) -> Result<()> {
    let mut cmd = std::process::Command::new(&args.python);
    cmd.arg("-m")
        .arg("sk2bgrow.cli")
        .arg("profile")
        .arg("--db")
        .arg(&args.db)
        .arg("--output")
        .arg(&args.output)
        .arg("--windows")
        .arg(args.output.join("windows.tsv"));
    for t in count_tables {
        cmd.arg(t);
    }
    ctx.say(format!("running statistics layer: {}", render(&cmd)));
    let status = cmd
        .status()
        .with_context(|| format!("could not launch '{}'; install the Python layer with `pip install -e python/`, or pass --no-stats", args.python))?;
    if !status.success() {
        bail!(
            "statistics layer failed ({status}); count tables are in {}",
            args.output.display()
        );
    }
    ctx.say(format!(
        "done: {}",
        args.output.join("output.tsv").display()
    ));
    Ok(())
}

fn render(cmd: &std::process::Command) -> String {
    let mut s = cmd.get_program().to_string_lossy().to_string();
    for a in cmd.get_args() {
        s.push(' ');
        s.push_str(&a.to_string_lossy());
    }
    s
}
