"""Interface files between the Rust and Python layers."""

import json

import numpy as np
import pandas as pd
import pytest

from conftest import make_counts
from sk2bgrow import io as sk_io


def write_counts(tmp_path, df, name="S1.counts.tsv"):
    p = tmp_path / name
    df[sk_io.COUNT_COLUMNS].to_csv(p, sep="\t", index=False)
    return p


def test_local_gc_is_dequantised_and_the_sentinel_becomes_nan(tmp_path):
    df = make_counts(n_per_enzyme=50, enzymes=("BcgI",), seed=0)
    df.loc[df.index[0], "local_gc"] = sk_io.GC_UNDEFINED
    df.loc[df.index[1], "local_gc"] = 100
    p = write_counts(tmp_path, df)
    back = sk_io.read_counts(p)
    assert np.isnan(back["gc"].iloc[0]), "255 must not be read as 127.5% GC"
    assert back["gc"].iloc[1] == pytest.approx(0.5)


def test_usable_flag_matches_the_rust_mask(tmp_path):
    df = make_counts(n_per_enzyme=40, enzymes=("BcgI",), seed=1)
    df.loc[df.index[0], "flags"] = sk_io.FLAG_UNIQUE_IN_GENOME | sk_io.FLAG_MASKED_SHARED
    df.loc[df.index[1], "flags"] = sk_io.FLAG_UNIQUE_IN_GENOME | sk_io.FLAG_NON_CHROMOSOMAL
    df.loc[df.index[2], "flags"] = sk_io.FLAG_UNIQUE_IN_GENOME | sk_io.FLAG_UNIQUE_ACROSS_DB
    back = sk_io.read_counts(write_counts(tmp_path, df))
    assert not back["usable"].iloc[0]
    assert not back["usable"].iloc[1]
    assert back["usable"].iloc[2]


def test_missing_columns_are_an_error(tmp_path):
    p = tmp_path / "bad.tsv"
    pd.DataFrame({"sample": ["S1"], "count": [3]}).to_csv(p, sep="\t", index=False)
    with pytest.raises(sk_io.Sk2bIoError, match="missing columns"):
        sk_io.read_counts(p)


def test_manifest_decodes_the_enzyme_bitset(tmp_path):
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "format_version": 1,
                "params": {"enzymes": 0b11, "gc_flank": 250, "sk2bgrow_version": "0.1.0"},
                "genomes": [
                    {"id": 0, "name": "ecoli", "genome_len": 4_641_652,
                     "contigs": [{"id": 0, "name": "chr", "length": 4_641_652, "offset": 0, "kind": "Chromosome"}],
                     "taxonomy": "s__E. coli", "ori": 3_923_883, "ori_confidence": 1.0},
                ],
                "n_anchors": 100,
            }
        )
    )
    m = sk_io.read_manifest(tmp_path)
    assert m.enzymes == ["BcgI", "AlfI"], "bit i must map to PANEL[i]"
    g = m.genome(0)
    assert g.ori == 3_923_883
    assert g.is_contiguous


def test_manifest_flags_a_fragmented_reference(tmp_path):
    contigs = [{"id": i, "name": f"c{i}", "length": 1000, "offset": i * 1000, "kind": "Chromosome"} for i in range(250)]
    (tmp_path / "manifest.json").write_text(
        json.dumps({"format_version": 1, "params": {"enzymes": 1},
                    "genomes": [{"id": 0, "name": "mag", "genome_len": 250_000, "contigs": contigs}],
                    "n_anchors": 10})
    )
    assert not sk_io.read_manifest(tmp_path).genome(0).is_contiguous


def test_unknown_genome_id_is_an_error(tmp_path):
    (tmp_path / "manifest.json").write_text(
        json.dumps({"format_version": 1, "params": {"enzymes": 1}, "genomes": [], "n_anchors": 0})
    )
    with pytest.raises(sk_io.Sk2bIoError, match="not in the database manifest"):
        sk_io.read_manifest(tmp_path).genome(7)


def test_missing_manifest_points_at_the_fix(tmp_path):
    with pytest.raises(sk_io.Sk2bIoError, match="sk2bgrow index"):
        sk_io.read_manifest(tmp_path)


def test_concat_rejects_tables_from_different_databases(tmp_path):
    a = make_counts(n_per_enzyme=50, enzymes=("BcgI",), sample="S1", seed=0)
    b = make_counts(n_per_enzyme=50, enzymes=("BcgI",), sample="S2", seed=99)  # different anchor positions
    pa = write_counts(tmp_path, a, "S1.counts.tsv")
    pb = write_counts(tmp_path, b, "S2.counts.tsv")
    with pytest.raises(sk_io.Sk2bIoError, match="different anchor database"):
        sk_io.concat_counts([pa, pb])


def test_concat_accepts_matching_tables(tmp_path):
    a = make_counts(n_per_enzyme=50, enzymes=("BcgI",), sample="S1", seed=0)
    b = a.copy()
    b["sample"] = "S2"
    pa = write_counts(tmp_path, a, "S1.counts.tsv")
    pb = write_counts(tmp_path, b, "S2.counts.tsv")
    joined = sk_io.concat_counts([pa, pb])
    assert set(joined["sample"]) == {"S1", "S2"}
    assert len(joined) == 2 * len(a)


def test_stats_path_is_derived_from_the_count_path(tmp_path):
    p = tmp_path / "SRR123.counts.tsv"
    assert sk_io.stats_path_for(p).name == "SRR123.stats.json"
