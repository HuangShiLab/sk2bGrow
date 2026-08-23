//! Drives the installed binary end to end: index -> profile -> audit.
//!
//! The Python statistics layer is skipped (`--no-stats`); it has its own suite in
//! `tests/python/`. What is checked here is the Rust half's contract with it:
//! the count table, the window table and the stats sidecar.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_sk2bgrow")
}

struct Rng(u64);
impl Rng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
}

fn workdir(tag: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("sk2bgrow-cli-{}-{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn write_genome(dir: &Path, name: &str, len: usize, seed: u64) -> (PathBuf, Vec<u8>) {
    let mut rng = Rng(seed);
    const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let seq: Vec<u8> = (0..len)
        .map(|_| B[(rng.next_u64() >> 33) as usize % 4])
        .collect();
    let p = dir.join(format!("{name}.fna"));
    let mut body = format!(">{name}_chr\n");
    for c in seq.chunks(70) {
        body.push_str(std::str::from_utf8(c).unwrap());
        body.push('\n');
    }
    std::fs::write(&p, body).unwrap();
    (p, seq)
}

/// A FASTQ of tiling 150 bp reads.
fn write_reads(dir: &Path, name: &str, genome: &[u8], step: usize) -> PathBuf {
    let p = dir.join(format!("{name}.fq"));
    let mut body = String::new();
    let mut i = 0usize;
    let mut n = 0usize;
    while i + 150 <= genome.len() {
        body.push_str(&format!(
            "@r{n}\n{}\n+\n{}\n",
            std::str::from_utf8(&genome[i..i + 150]).unwrap(),
            "I".repeat(150)
        ));
        i += step;
        n += 1;
    }
    std::fs::write(&p, body).unwrap();
    p
}

fn run(args: &[&str]) -> std::process::Output {
    let out = Command::new(bin())
        .args(args)
        .output()
        .expect("failed to run sk2bgrow");
    if !out.status.success() {
        panic!(
            "sk2bgrow {:?} failed ({})\nstdout:\n{}\nstderr:\n{}",
            args,
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    out
}

#[test]
fn index_profile_audit_round_trip() {
    let dir = workdir("e2e");
    let (g0, seq0) = write_genome(&dir, "g0", 120_000, 21);
    let (g1, _) = write_genome(&dir, "g1", 120_000, 22);
    let db = dir.join("db");

    // taxonomy is optional metadata; supply it so the manifest path is exercised.
    let tax = dir.join("tax.tsv");
    std::fs::write(
        &tax,
        "g0\td__Bacteria;s__Synthetic one\ng1\td__Bacteria;s__Synthetic two\n",
    )
    .unwrap();

    run(&[
        "index",
        g0.to_str().unwrap(),
        g1.to_str().unwrap(),
        "-o",
        db.to_str().unwrap(),
        "--enzymes",
        "BcgI,AlfI,CjeI,Hin4I",
        "-a",
        tax.to_str().unwrap(),
        "--write-tgt",
        "--quiet",
    ]);
    assert!(db.join("manifest.json").exists());
    assert!(db.join("anchors.bin").exists());
    assert!(
        db.join("tgt").join("g0.tgt").exists(),
        "--write-tgt produced no dump"
    );

    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(db.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["format_version"], 1);
    assert_eq!(manifest["genomes"].as_array().unwrap().len(), 2);
    assert_eq!(
        manifest["genomes"][0]["taxonomy"],
        "d__Bacteria;s__Synthetic one"
    );
    let n_anchors = manifest["n_anchors"].as_u64().unwrap();
    assert!(n_anchors > 500, "only {n_anchors} anchors");

    // --- profile --------------------------------------------------------
    let reads = write_reads(&dir, "S1", &seq0, 50);
    let out = dir.join("out");
    run(&[
        "profile",
        reads.to_str().unwrap(),
        "-d",
        db.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-stats",
        "--quiet",
    ]);

    let counts = out.join("S1.counts.tsv");
    let text = std::fs::read_to_string(&counts).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        sk2bgrow_core::count::COUNT_TABLE_HEADER
    );
    let rows: Vec<&str> = lines.collect();
    assert_eq!(rows.len() as u64, n_anchors, "one row per anchor expected");

    // Reads came only from g0, so g0 must carry essentially all the counts.
    let mut g0_total = 0u64;
    let mut g1_total = 0u64;
    for r in &rows {
        let f: Vec<&str> = r.split('\t').collect();
        let c: u64 = f[11].parse().unwrap();
        if f[2] == "g0" {
            g0_total += c
        } else {
            g1_total += c
        }
    }
    assert!(g0_total > 100, "g0 got {g0_total} counts");
    assert!(
        g1_total * 20 < g0_total,
        "g1 (absent from the sample) got {g1_total} vs g0 {g0_total}"
    );

    let stats: serde_json::Value =
        serde_json::from_slice(&std::fs::read(out.join("S1.stats.json")).unwrap()).unwrap();
    assert_eq!(stats["sample"], "S1");
    assert_eq!(stats["mode"], "wms");
    assert!(stats["counting"]["tag_matched"].as_u64().unwrap() > 0);
    assert_eq!(stats["em"]["converged"], true);

    let windows = std::fs::read_to_string(out.join("windows.tsv")).unwrap();
    assert!(windows.lines().next().unwrap().starts_with("window_id\t"));
    assert!(windows.lines().count() > 2);

    // --- audit ----------------------------------------------------------
    let report = dir.join("audit.html");
    run(&[
        "audit",
        db.to_str().unwrap(),
        "-o",
        report.to_str().unwrap(),
        "--quiet",
    ]);
    let html = std::fs::read_to_string(&report).unwrap();
    assert!(html.contains("<title>sk2bgrow audit</title>"));
    assert!(html.contains("g0") && html.contains("g1"));
    assert!(html.contains("max_gap"));

    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn digest_reproduces_the_density_table_shape() {
    let dir = workdir("digest");
    let (g0, _) = write_genome(&dir, "g0", 200_000, 31);
    let out = run(&[
        "digest",
        g0.to_str().unwrap(),
        "--enzymes",
        "all",
        "--quiet",
    ]);
    let text = String::from_utf8(out.stdout).unwrap();
    assert!(text.starts_with("genome\tenzyme\tsites\tdensity_per_mb"));
    assert!(
        text.contains("\tUNION\t"),
        "the union row is the point of the table"
    );
    assert!(text.contains("mean_spacing") && text.contains("max_gap"));
    // All 16 enzymes must appear, including the sparse ones.
    for e in ["BcgI", "PpiI", "PsrI", "CspCI", "BslFI"] {
        assert!(
            text.contains(&format!("\t{e}\t")),
            "{e} missing from the density table"
        );
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn unknown_enzyme_names_are_rejected_with_the_panel_listed() {
    let dir = workdir("badenzyme");
    let (g0, _) = write_genome(&dir, "g0", 20_000, 41);
    let out = Command::new(bin())
        .args([
            "index",
            g0.to_str().unwrap(),
            "-o",
            dir.join("db").to_str().unwrap(),
            "--enzymes",
            "EcoRI",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("EcoRI"));
    assert!(
        err.contains("BcgI"),
        "the error should list the panel; got: {err}"
    );
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn profile_rejects_enzymes_absent_from_the_database() {
    let dir = workdir("mismatch");
    let (g0, seq0) = write_genome(&dir, "g0", 60_000, 51);
    let db = dir.join("db");
    run(&[
        "index",
        g0.to_str().unwrap(),
        "-o",
        db.to_str().unwrap(),
        "--enzymes",
        "BcgI",
        "--quiet",
    ]);
    let reads = write_reads(&dir, "S1", &seq0, 400);
    let out = Command::new(bin())
        .args([
            "profile",
            reads.to_str().unwrap(),
            "-d",
            db.to_str().unwrap(),
            "-o",
            dir.join("out").to_str().unwrap(),
            "--enzymes",
            "AlfI",
            "--no-stats",
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "profiling with an unindexed enzyme should fail loudly"
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("AlfI") && err.contains("rebuild"),
        "unhelpful error: {err}"
    );
    std::fs::remove_dir_all(dir).ok();
}
