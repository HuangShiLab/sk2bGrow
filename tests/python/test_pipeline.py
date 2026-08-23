"""End-to-end: count table -> output.tsv, through the real CLI."""

import json

import numpy as np
import pandas as pd
import pytest

from conftest import ENZYMES, GENOME_LEN, ORI, make_counts
from sk2bgrow import cli, io as sk_io, report


def build_db(tmp_path, ori=None, n_contigs=1):
    db = tmp_path / "db"
    db.mkdir()
    bits = 0
    for name in ENZYMES:
        bits |= 1 << sk_io.ENZYME_PANEL.index(name)
    contigs = [
        {"id": i, "name": f"c{i}", "length": GENOME_LEN // n_contigs, "offset": i * (GENOME_LEN // n_contigs),
         "kind": "Chromosome"}
        for i in range(n_contigs)
    ]
    (db / "manifest.json").write_text(
        json.dumps(
            {
                "format_version": 1,
                "params": {"enzymes": bits, "gc_flank": 250, "min_contig_len": 500,
                           "reject_ambiguous_tags": True, "sk2bgrow_version": "0.1.0"},
                "genomes": [
                    {"id": 0, "name": "ecoli", "taxonomy": "d__Bacteria;s__Escherichia coli",
                     "contigs": contigs, "genome_len": GENOME_LEN,
                     "ori": ori, "ori_confidence": 1.0 if ori else 0.0}
                ],
                "n_anchors": 0,
            }
        )
    )
    return db


def write_sample(tmp_path, df, sample, containment=0.98):
    p = tmp_path / f"{sample}.counts.tsv"
    df[sk_io.COUNT_COLUMNS].to_csv(p, sep="\t", index=False)
    (tmp_path / f"{sample}.stats.json").write_text(
        json.dumps(
            {
                "sample": sample,
                "counting": {"reads_total": 1_000_000, "tag_matched": 900_000},
                "em": {"iterations": 3, "converged": True,
                       "genomes": [{"genome_id": 0, "containment": containment, "lambda": 8.0}]},
            }
        )
    )
    return p


@pytest.mark.parametrize("truth", [0.5, 1.0, 2.0])
def test_pipeline_recovers_the_planted_ptr(tmp_path, truth):
    db = build_db(tmp_path)
    counts = make_counts(log2_ptr=truth, per_anchor_depth=8.0, n_per_enzyme=1500, seed=11)
    p = write_sample(tmp_path, counts, "S1")
    out = tmp_path / "out"
    rc = cli.main(["profile", str(p), "--db", str(db), "--output", str(out)])
    assert rc == 0

    res = pd.read_csv(out / "output.tsv", sep="\t", na_values=["NA"])
    assert len(res) == 1
    row = res.iloc[0]
    assert row["log2(PTR)"] == pytest.approx(truth, abs=0.25)
    assert row["PTR"] == pytest.approx(2**truth, rel=0.25)
    assert row["n_enzymes"] == len(ENZYMES)
    assert row["containment"] == pytest.approx(0.98)
    assert row["taxonomy"].startswith("d__Bacteria")


def test_output_columns_stay_pilea_compatible(tmp_path):
    db = build_db(tmp_path)
    p = write_sample(tmp_path, make_counts(n_per_enzyme=800, seed=12), "S1")
    out = tmp_path / "out"
    cli.main(["profile", str(p), "--db", str(db), "--output", str(out)])
    cols = pd.read_csv(out / "output.tsv", sep="\t", nrows=0).columns.tolist()
    for c in ["coverage", "dispersion", "fraction", "containment", "PTR", "log2(PTR)"]:
        assert c in cols, f"Pilea-compatible column {c} disappeared"
    for c in ["enzyme_consistency", "n_anchors", "ori_confidence"]:
        assert c in cols, f"sk2bGrow column {c} missing"
    assert cols == report.OUTPUT_COLUMNS


def test_known_ori_is_taken_from_the_manifest(tmp_path):
    db = build_db(tmp_path, ori=ORI)
    p = write_sample(tmp_path, make_counts(log2_ptr=1.0, n_per_enzyme=1500, seed=13), "S1")
    out = tmp_path / "out"
    cli.main(["profile", str(p), "--db", str(db), "--output", str(out)])
    per_enzyme = pd.read_csv(out / "per_enzyme.tsv", sep="\t")
    assert (per_enzyme["ori"] == ORI).all()
    assert (per_enzyme["ori_confidence"] == 1.0).all()


def test_intermediate_tables_are_written(tmp_path):
    db = build_db(tmp_path)
    p = write_sample(tmp_path, make_counts(n_per_enzyme=800, seed=14), "S1")
    out = tmp_path / "out"
    cli.main(["profile", str(p), "--db", str(db), "--output", str(out)])
    for name in ["windows.rates.tsv", "per_enzyme.tsv", "output.tsv"]:
        assert (out / name).exists(), f"{name} was not written"
    w = pd.read_csv(out / "windows.rates.tsv", sep="\t")
    assert set(w["enzyme"]) == set(ENZYMES)
    assert w["gc_corrected"].all()


def test_two_samples_produce_two_rows_and_a_delta(tmp_path):
    db = build_db(tmp_path)
    slow = make_counts(log2_ptr=0.3, n_per_enzyme=1500, sample="T0", seed=15)
    fast = slow.copy()
    fast["sample"] = "T1"
    # Re-draw counts under a steeper gradient at the same anchor positions —
    # exactly what a deterministic sketch guarantees across samples.
    rng = np.random.default_rng(16)
    wrapped = np.mod(fast["position"].to_numpy() - ORI, GENOME_LEN)
    d = np.minimum(wrapped, GENOME_LEN - wrapped)
    fast["count"] = rng.poisson(8.0 * np.exp2(-1.8 * d / (GENOME_LEN / 2)))
    p0 = write_sample(tmp_path, slow, "T0")
    p1 = write_sample(tmp_path, fast, "T1")
    out = tmp_path / "out"
    assert cli.main(["profile", str(p0), str(p1), "--db", str(db), "--output", str(out)]) == 0

    res = pd.read_csv(out / "output.tsv", sep="\t", na_values=["NA"])
    assert len(res) == 2
    t0 = res[res["sample"] == "T0"]["log2(PTR)"].iloc[0]
    t1 = res[res["sample"] == "T1"]["log2(PTR)"].iloc[0]
    assert t1 > t0 + 1.0

    delta = tmp_path / "delta.tsv"
    assert cli.main(["dynamics", str(out / "output.tsv"), "--output", str(delta), "--baseline", "T0"]) == 0
    d = pd.read_csv(delta, sep="\t")
    assert d["delta_log2_ptr"].iloc[0] == pytest.approx(t1 - t0, abs=1e-6)
    assert d["p"].iloc[0] < 0.01


def test_a_rogue_enzyme_fails_the_consistency_gate(tmp_path):
    db = build_db(tmp_path)
    counts = make_counts(log2_ptr=1.0, n_per_enzyme=1500, seed=17)
    # AlfI sees a flat profile: a digestion failure, or a wrong reference region.
    rng = np.random.default_rng(18)
    sel = counts["enzyme"] == "AlfI"
    counts.loc[sel, "count"] = rng.poisson(8.0, size=int(sel.sum()))
    p = write_sample(tmp_path, counts, "S1")
    out = tmp_path / "out"
    cli.main(["profile", str(p), "--db", str(db), "--output", str(out)])
    row = pd.read_csv(out / "output.tsv", sep="\t", na_values=["NA"]).iloc[0]
    # An enzyme that sees no gradient fails to fit rather than reporting a
    # deviant number, so it never reaches the Q statistic. The fit-rate gate is
    # what catches it — without that, the surviving enzymes agree with each
    # other and a sample where a quarter of the panel saw nothing reads as clean.
    assert row["n_enzymes"] < row["n_enzymes_attempted"]
    assert row["enzyme_fit_rate"] < 0.8
    assert not bool(row["pass_qc"])
    assert "produced a fit" in row["qc_reason"]
    assert "AlfI" in row["excluded"]
    # The estimate is still reported, not deleted.
    assert np.isfinite(row["log2(PTR)"])


def test_a_deviant_enzyme_fails_the_consistency_gate(tmp_path):
    """The other discordance mode: an enzyme that fits, but to a different slope."""
    db = build_db(tmp_path)
    counts = make_counts(log2_ptr=1.0, n_per_enzyme=1500, seed=21)
    # CjeI sees a much steeper gradient — a mis-assembled region, say.
    rng = np.random.default_rng(22)
    sel = counts["enzyme"] == "CjeI"
    wrapped = np.mod(counts.loc[sel, "position"].to_numpy() - ORI, GENOME_LEN)
    d = np.minimum(wrapped, GENOME_LEN - wrapped)
    counts.loc[sel, "count"] = rng.poisson(8.0 * np.exp2(-2.6 * d / (GENOME_LEN / 2)))
    p = write_sample(tmp_path, counts, "S1")
    out = tmp_path / "out"
    cli.main(["profile", str(p), "--db", str(db), "--output", str(out)])
    row = pd.read_csv(out / "output.tsv", sep="\t", na_values=["NA"]).iloc[0]
    assert row["enzyme_consistency"] < 0.05
    assert not bool(row["pass_qc"])
    assert "disagree" in row["qc_reason"]
    # Heterogeneity detected -> the interval widens rather than staying tight.
    assert row["fusion_model"] == "random"
    assert np.isfinite(row["log2(PTR)"])


def test_qc_flags_low_coverage_without_dropping_the_row(tmp_path):
    db = build_db(tmp_path)
    p = write_sample(tmp_path, make_counts(per_anchor_depth=1.2, n_per_enzyme=1500, seed=19), "S1")
    out = tmp_path / "out"
    cli.main(["profile", str(p), "--db", str(db), "--output", str(out), "--min-coverage", "5"])
    row = pd.read_csv(out / "output.tsv", sep="\t", na_values=["NA"]).iloc[0]
    assert not bool(row["pass_qc"])
    assert "coverage" in row["qc_reason"]
    assert np.isfinite(row["log2(PTR)"]), "a failing row keeps its estimate"


def test_fragmented_reference_falls_back_to_sorted_regression(tmp_path):
    db = build_db(tmp_path, n_contigs=300)
    counts = make_counts(log2_ptr=1.0, n_per_enzyme=1500, seed=20)
    p = write_sample(tmp_path, counts, "S1")
    out = tmp_path / "out"
    cli.main(["profile", str(p), "--db", str(db), "--output", str(out)])
    per_enzyme = pd.read_csv(out / "per_enzyme.tsv", sep="\t")
    assert (per_enzyme["method"] == "sorted_ransac").all()
