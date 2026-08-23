//! `sk2bgrow` — command line front end.
//!
//! The subcommand shape mirrors Pilea's two-stage design (`index` once,
//! `profile` per batch) so an existing Pilea benchmark script can swap tools
//! with minimal edits — see `docs/cli.md`.

use anyhow::Result;
use clap::{Parser, Subcommand};

mod cmd;

#[derive(Parser)]
#[command(
    name = "sk2bgrow",
    version,
    about = "PTR estimation from deterministic 2bRAD/TGT anchors",
    long_about = "Estimate bacterial peak-to-trough ratios (PTR) from WMS reads (route A) \
or real 2bRAD reads (route B), using the 16-enzyme Type IIB anchor panel as a \
deterministic sketch.\n\n\
Typical run:\n  \
sk2bgrow index genomes/*.fna -o db --enzymes all\n  \
sk2bgrow profile reads/*.fq.gz -d db -o out/"
)]
struct Cli {
    /// Worker threads. 0 uses one per available core.
    #[arg(long, global = true, default_value_t = 0)]
    threads: usize,

    /// Suppress progress messages on stderr.
    #[arg(short, long, global = true)]
    quiet: bool,

    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Build an anchor database from reference genomes (offline, once).
    Index(cmd::index::Args),
    /// Count reads against a database and estimate PTR.
    Profile(cmd::profile::Args),
    /// Compare PTR across samples (delta-PTR / time series).
    Dynamics(cmd::dynamics::Args),
    /// Report anchor density and blind spots for a database.
    Audit(cmd::audit::Args),
    /// Digest one genome and print the density table (reproduces report §4.1).
    Digest(cmd::digest::Args),
    /// Order and orient draft MAG contigs against a reference.
    Scaffold(cmd::scaffold::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()?;
    }
    let ctx = cmd::Ctx { quiet: cli.quiet };
    match cli.command {
        Command::Index(a) => cmd::index::run(a, &ctx),
        Command::Profile(a) => cmd::profile::run(a, &ctx),
        Command::Dynamics(a) => cmd::dynamics::run(a, &ctx),
        Command::Audit(a) => cmd::audit::run(a, &ctx),
        Command::Digest(a) => cmd::digest::run(a, &ctx),
        Command::Scaffold(a) => cmd::scaffold::run(a, &ctx),
    }
}
