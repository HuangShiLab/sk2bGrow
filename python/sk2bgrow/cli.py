"""Command line entry point for the statistics layer.

Normally invoked by the Rust binary (``sk2bgrow profile`` shells out to
``python -m sk2bgrow.cli profile``), but it runs standalone against any count
table, which is what makes it possible to iterate on the statistics without
re-counting reads.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from . import __version__, dynamics, fit, fusion, gc_bias, io as sk_io, report, simulate, ztp


def _profile(args: argparse.Namespace) -> int:
    if args.window_anchors != "auto":
        args.window_anchors = int(args.window_anchors)
    manifest = sk_io.read_manifest(args.db)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    counts = sk_io.concat_counts(args.counts)
    print(f"[sk2bgrow.py] {len(counts):,} anchor rows, {counts['sample'].nunique()} sample(s), "
          f"{counts['genome'].nunique()} genome(s)", file=sys.stderr)

    if args.gc_correct:
        curves = gc_bias.fit_curves(counts, frac=args.loess_frac)
        counts = gc_bias.add_anchor_offsets(counts, curves)
        if curves:
            worst = max(curves.values(), key=lambda c: c.amplitude)
            print(f"[sk2bgrow.py] GC curves for {len(curves)} enzyme(s); "
                  f"largest amplitude {worst.amplitude:.2f} log2 ({worst.enzyme})", file=sys.stderr)
        else:
            print("[sk2bgrow.py] no enzyme had enough anchors for a GC curve; skipping correction", file=sys.stderr)

    windows = ztp.window_rates(
        counts,
        anchors_per_window=args.window_anchors,
        window_cap=args.window_cap,
        model=args.count_model,
        use_precomputed_windows=args.use_rust_windows,
    )
    if windows.empty:
        print("[sk2bgrow.py] no windows could be built; is the coverage zero?", file=sys.stderr)
        return 2
    if args.gc_correct:
        windows = gc_bias.apply_to_windows(windows)
    windows.to_csv(outdir / "windows.rates.tsv", sep="\t", index=False, float_format="%.10g")

    per_enzyme = fit.fit_windows(
        windows, manifest, method=args.method, min_windows=args.min_windows, shared_ori=args.shared_ori
    )
    per_enzyme.to_csv(outdir / "per_enzyme.tsv", sep="\t", index=False, float_format="%.10g")
    n_ok = int(per_enzyme["ok"].sum())
    print(f"[sk2bgrow.py] {n_ok}/{len(per_enzyme)} per-enzyme fits usable", file=sys.stderr)

    fused = fusion.fuse_table(per_enzyme, alpha=args.alpha, min_anchors=args.min_anchors_per_enzyme)

    stats = {}
    for c in args.counts:
        sp = sk_io.stats_path_for(c)
        if sp.exists():
            s = sk_io.read_stats(sp)
            stats[str(s.get("sample", sp.stem))] = s

    table = report.assemble(fused, manifest, stats)
    table = report.apply_qc(
        table,
        min_coverage=args.min_coverage,
        min_fraction=args.min_fraction,
        max_dispersion=args.max_dispersion,
        min_containment=args.min_containment,
        consistency_alpha=args.alpha,
        min_enzyme_fit_rate=args.min_enzyme_fit_rate,
    )
    out = report.write_output(table, outdir / "output.tsv")
    passed = int(table["pass_qc"].sum())
    print(f"[sk2bgrow.py] {passed}/{len(table)} genome-sample estimates passed QC -> {out}", file=sys.stderr)
    for _, r in table[~table["pass_qc"]].iterrows():
        print(f"[sk2bgrow.py]   {r['sample']}/{r['genome']}: {r['qc_reason']}", file=sys.stderr)

    if args.figures:
        written = report.write_qc_figures(windows, per_enzyme, table, outdir / "qc")
        print(f"[sk2bgrow.py] {len(written)} QC figure(s) in {outdir / 'qc'}", file=sys.stderr)
    return 0


def _dynamics(args: argparse.Namespace) -> int:
    outputs = dynamics.read_outputs(args.inputs)
    meta = pd.read_csv(args.metadata, sep="\t") if args.metadata else None
    deltas = dynamics.delta_ptr(outputs, baseline=args.baseline, metadata=meta, qc_only=not args.include_failed)
    if deltas.empty:
        print("[sk2bgrow.py] nothing to compare (check --baseline and QC status)", file=sys.stderr)
        return 2
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    deltas.to_csv(out, sep="\t", index=False, float_format="%.10g")
    print(f"[sk2bgrow.py] {len(deltas)} delta rows -> {out}", file=sys.stderr)

    if meta is not None and "timepoint" in meta.columns:
        trends = dynamics.trend_test(deltas)
        tp = out.with_suffix(".trend.tsv")
        trends.to_csv(tp, sep="\t", index=False, float_format="%.10g")
        print(f"[sk2bgrow.py] trend tests -> {tp}", file=sys.stderr)
    return 0


def _simulate(args: argparse.Namespace) -> int:
    if args.route == "a":
        df = simulate.route_a(n_reps=args.reps, estimator=args.estimator, seed=args.seed)
        cols = ["anchor_set", "n_anchors", "coverage", "true_log2_ptr", "bias", "rmse", "sd", "n_ok"]
    else:
        df = simulate.route_b(n_reps=args.reps, estimator=args.estimator, seed=args.seed)
        cols = ["anchor_set", "n_anchors", "per_site_depth", "sigma_eff", "true_log2_ptr", "bias", "rmse", "sd", "n_ok"]
    print(df[cols].to_string(index=False))
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, sep="\t", index=False, float_format="%.10g")
        print(f"\n[sk2bgrow.py] wrote {args.output}", file=sys.stderr)
    return 0


def _manifest(args: argparse.Namespace) -> int:
    m = sk_io.read_manifest(args.db)
    print(json.dumps(
        {
            "version": m.version,
            "enzymes": m.enzymes,
            "gc_flank": m.gc_flank,
            "genomes": [
                {"id": g.id, "name": g.name, "len": g.genome_len, "contigs": g.n_contigs,
                 "ori": g.ori, "contiguous": g.is_contiguous}
                for g in m.genomes.values()
            ],
        },
        indent=2,
    ))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="sk2bgrow.cli", description="sk2bGrow statistics layer")
    p.add_argument("--version", action="version", version=f"sk2bgrow {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("profile", help="count tables -> output.tsv")
    pr.add_argument("counts", nargs="+", help="per-sample *.counts.tsv from `sk2bgrow profile`")
    pr.add_argument("--db", required=True, help="anchor database directory")
    pr.add_argument("--output", required=True, help="output directory")
    pr.add_argument("--windows", help="window table from the Rust layer (used with --use-rust-windows)")
    pr.add_argument("--window-anchors", default="auto",
                    help="anchors per window: an integer, or 'auto' to size each enzyme by its anchor count")
    pr.add_argument("--window-cap", type=int, default=100, help="upper bound on anchors per window in auto mode")
    pr.add_argument("--per-enzyme-ori", dest="shared_ori", action="store_false",
                    help="let every enzyme search for its own origin instead of sharing one")
    pr.add_argument("--use-rust-windows", action="store_true",
                    help="group by the Rust window_id instead of re-windowing per enzyme (Pilea-parity path)")
    pr.add_argument("--count-model", choices=["auto", "ztp", "nb"], default="auto")
    pr.add_argument("--method", choices=["auto", "v_shape", "sorted"], default="auto")
    pr.add_argument("--min-windows", type=int, default=5)
    pr.add_argument("--no-gc-correct", dest="gc_correct", action="store_false", help="skip per-enzyme GC correction")
    pr.add_argument("--loess-frac", type=float, default=0.4)
    pr.add_argument("--alpha", type=float, default=0.05, help="cross-enzyme consistency threshold")
    pr.add_argument("--min-anchors-per-enzyme", type=int, default=fusion.MIN_ANCHORS_PER_ENZYME)
    pr.add_argument("--min-coverage", type=float, default=1.0, help="QC floor (Pilea uses 5)")
    pr.add_argument("--min-fraction", type=float, default=0.75)
    pr.add_argument("--max-dispersion", type=float, default=5.0)
    pr.add_argument("--min-containment", type=float, default=0.5)
    pr.add_argument("--min-enzyme-fit-rate", type=float, default=0.8,
                    help="flag a genome when fewer than this fraction of enzymes produced a usable fit")
    pr.add_argument("--figures", action="store_true", help="write QC figures")
    pr.set_defaults(func=_profile, gc_correct=True, shared_ori=True)

    dy = sub.add_parser("dynamics", help="compare PTR across samples")
    dy.add_argument("inputs", nargs="+", help="output.tsv files")
    dy.add_argument("--output", required=True)
    dy.add_argument("--metadata", help="TSV: sample<TAB>group[<TAB>timepoint]")
    dy.add_argument("--baseline", help="baseline sample or group")
    dy.add_argument("--include-failed", action="store_true", help="include rows that failed QC")
    dy.set_defaults(func=_dynamics)

    si = sub.add_parser("simulate", help="reproduce the design report's section 5 simulations")
    si.add_argument("route", choices=["a", "b"])
    si.add_argument("--reps", type=int, default=150)
    si.add_argument("--estimator", choices=["sorted", "v_shape"], default="sorted")
    si.add_argument("--seed", type=int, default=0)
    si.add_argument("--output")
    si.set_defaults(func=_simulate)

    mf = sub.add_parser("manifest", help="print a database manifest summary")
    mf.add_argument("db")
    mf.set_defaults(func=_manifest)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
