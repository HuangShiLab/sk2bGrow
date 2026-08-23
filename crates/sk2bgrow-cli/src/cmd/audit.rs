//! `sk2bgrow audit` — anchor density and blind-spot report for a database.
//!
//! The report argues (§4.2) that the real advantage of a deterministic sketch is
//! that sparse regions are *knowable at build time*: a random sketch's blind
//! spots are hidden, an anchor library's are auditable. This subcommand is that
//! audit — the build-quality gate before a database is used for anything.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use sk2bgrow_core::anchor_db::AnchorDb;
use sk2bgrow_core::enzyme::PANEL;
use sk2bgrow_core::window::spacing_stats;

use super::Ctx;

#[derive(ClapArgs)]
pub struct Args {
    /// Database directory.
    pub db: PathBuf,

    /// Output file; `.html` renders a report, anything else writes TSV.
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Gaps wider than this many bp count as blind spots.
    #[arg(long, default_value_t = 5_000)]
    pub wide_gap: u64,

    /// Flag genomes with fewer than this many usable anchors per enzyme.
    #[arg(long, default_value_t = 50)]
    pub min_anchors_per_enzyme: usize,
}

pub fn run(args: Args, ctx: &Ctx) -> Result<()> {
    let db = AnchorDb::load(&args.db)
        .with_context(|| format!("loading database {}", args.db.display()))?;
    ctx.say(format!("auditing {} genome(s)", db.genomes.len()));

    let mut rows = Vec::new();
    for g in &db.genomes {
        let s = spacing_stats(&db, g.id, args.wide_gap, true);
        let per_enzyme = db.per_enzyme_counts(g.id);
        let thin: Vec<&str> = PANEL
            .iter()
            .filter(|e| {
                db.params.enzymes.contains(e.idx)
                    && per_enzyme[e.idx as usize] < args.min_anchors_per_enzyme
            })
            .map(|e| e.name)
            .collect();
        rows.push((g, s, per_enzyme, thin));
    }

    let mut tsv = String::from(
        "genome\tgenome_len\tcontigs\tusable_anchors\tper_25kb\tmean_spacing\tmedian_spacing\tp99_spacing\tmax_gap\twide_gaps\tthin_enzymes\tori\n",
    );
    for (g, s, _, thin) in &rows {
        tsv.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{}\n",
            g.name,
            g.genome_len,
            g.contigs.len(),
            s.n_anchors,
            if g.genome_len > 0 {
                s.n_anchors as f64 / (g.genome_len as f64 / 25_000.0)
            } else {
                0.0
            },
            s.mean_spacing,
            s.median_spacing,
            s.p99_spacing,
            s.max_gap,
            s.n_wide_gaps,
            if thin.is_empty() {
                "-".to_string()
            } else {
                thin.join(",")
            },
            g.ori.map(|o| o.to_string()).unwrap_or_else(|| "-".into()),
        ));
    }

    let mut per_enzyme_tsv = String::from("genome\tenzyme\tusable_anchors\tper_mb\n");
    for (g, _, per_enzyme, _) in &rows {
        let mb = (g.genome_len as f64 / 1e6).max(f64::MIN_POSITIVE);
        for e in PANEL.iter().filter(|e| db.params.enzymes.contains(e.idx)) {
            let n = per_enzyme[e.idx as usize];
            per_enzyme_tsv.push_str(&format!(
                "{}\t{}\t{}\t{:.1}\n",
                g.name,
                e.name,
                n,
                n as f64 / mb
            ));
        }
    }

    match &args.output {
        Some(p) if p.extension().is_some_and(|e| e == "html") => {
            std::fs::write(p, render_html(&tsv, &per_enzyme_tsv, &db))?;
            ctx.say(format!("wrote {}", p.display()));
        }
        Some(p) => {
            std::fs::write(p, format!("{tsv}\n# per enzyme\n{per_enzyme_tsv}"))?;
            ctx.say(format!("wrote {}", p.display()));
        }
        None => print!("{tsv}\n# per enzyme\n{per_enzyme_tsv}"),
    }

    // A non-zero exit would break batch pipelines, so problems are warnings.
    for (g, s, _, thin) in &rows {
        if s.n_wide_gaps > 0 {
            ctx.say(format!(
                "  {}: {} gaps > {} bp (worst {} bp)",
                g.name, s.n_wide_gaps, args.wide_gap, s.max_gap
            ));
        }
        if !thin.is_empty() {
            ctx.say(format!("  {}: thin enzymes {}", g.name, thin.join(",")));
        }
    }
    Ok(())
}

fn render_html(summary: &str, per_enzyme: &str, db: &AnchorDb) -> String {
    let table = |tsv: &str| -> String {
        let mut rows = tsv.lines();
        let head = rows.next().unwrap_or("");
        let mut h = String::from("<table><thead><tr>");
        for c in head.split('\t') {
            h.push_str(&format!("<th>{c}</th>"));
        }
        h.push_str("</tr></thead><tbody>");
        for r in rows.filter(|l| !l.trim().is_empty()) {
            h.push_str("<tr>");
            for c in r.split('\t') {
                h.push_str(&format!("<td>{c}</td>"));
            }
            h.push_str("</tr>");
        }
        h.push_str("</tbody></table>");
        h
    };
    format!(
        "<!doctype html><meta charset=\"utf-8\"><title>sk2bgrow audit</title>\
<style>body{{font:14px/1.5 system-ui,sans-serif;margin:2rem;max-width:70rem}}\
table{{border-collapse:collapse;margin:1rem 0;font-variant-numeric:tabular-nums}}\
th,td{{border:1px solid #ddd;padding:.3rem .6rem;text-align:right}}\
th:first-child,td:first-child{{text-align:left}}th{{background:#f5f5f5}}\
h1{{font-size:1.4rem}}h2{{font-size:1.1rem;margin-top:2rem}}</style>\
<h1>sk2bgrow anchor database audit</h1>\
<p>{} genomes, {} anchors, {} enzymes, built with sk2bgrow {}.</p>\
<h2>Per genome</h2>{}<h2>Per enzyme</h2>{}",
        db.genomes.len(),
        db.n_anchors(),
        db.params.enzymes.len(),
        db.params.sk2bgrow_version,
        table(summary),
        table(per_enzyme)
    )
}
