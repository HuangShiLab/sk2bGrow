//! `sk2bgrow dynamics` — delegate multi-sample delta-PTR to the Python layer.
//!
//! The longitudinal model (anchor x sample count matrix, repeated-measures
//! mixed model, Poisson-PCA reordering check) lives in
//! `python/sk2bgrow/dynamics.py`. This subcommand exists so the whole workflow
//! is reachable from one binary.

use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;

use super::Ctx;

#[derive(ClapArgs)]
pub struct Args {
    /// Per-sample `output.tsv` files from `sk2bgrow profile`.
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output TSV of delta-PTR.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Sample metadata TSV: `sample<TAB>group[<TAB>timepoint]`.
    #[arg(short = 'm', long)]
    pub metadata: Option<PathBuf>,

    /// Reference sample or group that deltas are measured against.
    #[arg(long)]
    pub baseline: Option<String>,

    #[arg(long, default_value = "python3")]
    pub python: String,
}

pub fn run(args: Args, ctx: &Ctx) -> Result<()> {
    for p in &args.inputs {
        if !p.exists() {
            bail!("input does not exist: {}", p.display());
        }
    }
    let mut cmd = std::process::Command::new(&args.python);
    cmd.arg("-m")
        .arg("sk2bgrow.cli")
        .arg("dynamics")
        .arg("--output")
        .arg(&args.output);
    if let Some(m) = &args.metadata {
        cmd.arg("--metadata").arg(m);
    }
    if let Some(b) = &args.baseline {
        cmd.arg("--baseline").arg(b);
    }
    for p in &args.inputs {
        cmd.arg(p);
    }
    ctx.say(format!(
        "{} sample file(s) -> {}",
        args.inputs.len(),
        args.output.display()
    ));
    let status = cmd.status().with_context(|| {
        format!(
            "could not launch '{}'; install the Python layer with `pip install -e python/`",
            args.python
        )
    })?;
    if !status.success() {
        bail!("dynamics failed ({status})");
    }
    Ok(())
}
