"""Adaptive per-enzyme window sizing and the shared origin estimate.

Both exist because a fixed policy quietly discards the sparse half of the panel,
which would defeat the point of a 16-enzyme design.
"""

import numpy as np
import pandas as pd
import pytest

from conftest import ENZYMES, GENOME_LEN, ORI, make_counts
from sk2bgrow import fit, ztp

# Anchors per Mb on E. coli K-12, from design report table 4.1.
ECOLI_DENSITY = {
    "PpiI": 73, "BplI": 83, "PsrI": 92, "AloI": 112, "CspCI": 132, "FalI": 158,
    "BaeI": 171, "BsaXI": 212, "Bsp24I": 351, "AlfI": 436, "BslFI": 576,
    "BcgI": 632, "HaeIV": 745, "CjePI": 1701, "CjeI": 1910, "Hin4I": 1650,
}
ECOLI_MB = 4.64


def test_auto_sizing_keeps_every_enzyme_above_the_fit_minimum():
    """At a flat 100 anchors/window the three sparsest enzymes drop out entirely."""
    fixed_losses, auto_losses = [], []
    for enzyme, per_mb in ECOLI_DENSITY.items():
        n = int(per_mb * ECOLI_MB)
        if n // 100 < 5:
            fixed_losses.append(enzyme)
        if n // ztp.auto_window_size(n) < 5:
            auto_losses.append(enzyme)
    assert set(fixed_losses) == {"PpiI", "BplI", "PsrI"}
    assert auto_losses == [], f"auto sizing still loses {auto_losses}"


def test_auto_sizing_is_bounded_at_both_ends():
    assert ztp.auto_window_size(10) == ztp.MIN_ANCHORS_PER_WINDOW
    assert ztp.auto_window_size(1_000_000) == 100
    assert ztp.auto_window_size(1_000_000, cap=250) == 250
    # Dense enzymes stay at the cap; sparse ones shrink rather than disappear.
    assert ztp.auto_window_size(int(1910 * ECOLI_MB)) == 100
    assert ztp.auto_window_size(int(73 * ECOLI_MB)) == ztp.MIN_ANCHORS_PER_WINDOW


def test_auto_sizing_is_monotone_in_anchor_count():
    sizes = [ztp.auto_window_size(n) for n in range(50, 5_000, 50)]
    assert all(b >= a for a, b in zip(sizes, sizes[1:]))


def test_sparse_and_dense_enzymes_both_get_windows():
    dense = make_counts(n_per_enzyme=2000, enzymes=("CjeI",), seed=1)
    sparse = make_counts(n_per_enzyme=300, enzymes=("PpiI",), seed=2)
    df = pd.concat([dense, sparse], ignore_index=True)
    w = ztp.window_rates(df, anchors_per_window="auto")
    per_enzyme = w.groupby("enzyme").size()
    assert per_enzyme["CjeI"] >= 15
    assert per_enzyme["PpiI"] >= 8, "the sparse enzyme was squeezed out"
    # ...at different window sizes.
    sizes = w.groupby("enzyme")["anchors_per_window"].first()
    assert sizes["CjeI"] > sizes["PpiI"]


def test_explicit_window_size_is_respected():
    df = make_counts(n_per_enzyme=1000, enzymes=("BcgI",), seed=3)
    w = ztp.window_rates(df, anchors_per_window=50)
    assert (w["anchors_per_window"] == 50).all()
    assert 15 <= len(w) <= 25


def test_shared_ori_beats_per_enzyme_search_when_windows_are_few():
    """Searching per enzyme injects variance that looks like enzyme disagreement."""
    df = make_counts(log2_ptr=1.0, per_anchor_depth=8.0, n_per_enzyme=600, seed=4)
    windows = ztp.window_rates(df, anchors_per_window=50)
    ori, conf = fit.find_shared_ori(windows, GENOME_LEN)
    assert ori is not None
    assert fit.circular_distance(np.array([ori]), ORI, GENOME_LEN)[0] < 250_000
    assert 0.0 <= conf <= 1.0


def test_shared_ori_pins_every_enzyme_to_one_coordinate(manifest):
    df = make_counts(log2_ptr=1.0, n_per_enzyme=900, seed=5)
    windows = ztp.window_rates(df, anchors_per_window=60)
    shared = fit.fit_windows(windows, manifest, shared_ori=True)
    assert shared["ori"].nunique() == 1, "enzymes were fitted at different origins"

    per_enzyme = fit.fit_windows(windows, manifest, shared_ori=False)
    assert per_enzyme["ori"].nunique() > 1, "the per-enzyme path should search independently"


def test_pooling_locates_the_origin_better_than_any_single_enzyme(manifest):
    """The direct claim: more windows in one search beats N smaller searches.

    Asserted on origin *accuracy* rather than on the spread of the resulting PTR
    estimates — with enough windows per enzyme the two agree closely and a spread
    comparison is a coin flip, but the pooled search sees 4x the data and locates
    the origin better regardless.
    """
    df = make_counts(log2_ptr=0.6, per_anchor_depth=5.0, n_per_enzyme=500, seed=6)
    windows = ztp.window_rates(df, anchors_per_window=50)

    pooled_ori, _ = fit.find_shared_ori(windows, GENOME_LEN)
    pooled_err = fit.circular_distance(np.array([pooled_ori]), ORI, GENOME_LEN)[0]

    independent = fit.fit_windows(windows, manifest, shared_ori=False)
    per_enzyme_err = np.median(
        [fit.circular_distance(np.array([o]), ORI, GENOME_LEN)[0] for o in independent["ori"].dropna()]
    )
    assert pooled_err <= per_enzyme_err, f"pooled {pooled_err/1e3:.0f} kb vs per-enzyme median {per_enzyme_err/1e3:.0f} kb"


def test_a_curated_ori_overrides_the_search(manifest_with_ori):
    df = make_counts(log2_ptr=1.0, n_per_enzyme=900, seed=7)
    windows = ztp.window_rates(df, anchors_per_window=60)
    out = fit.fit_windows(windows, manifest_with_ori, shared_ori=True)
    assert (out["ori"] == ORI).all()
    assert (out["ori_confidence"] == 1.0).all()


def test_find_shared_ori_centres_each_enzyme_first():
    """A per-enzyme efficiency offset must not tilt the pooled origin search."""
    df = make_counts(log2_ptr=1.0, n_per_enzyme=900, enzymes=tuple(ENZYMES), seed=8)
    # Give one enzyme 8x the depth: a pooled fit that ignored the offset would
    # be dragged toward that enzyme's windows.
    df.loc[df["enzyme"] == "AlfI", "count"] = df.loc[df["enzyme"] == "AlfI", "count"] * 8
    windows = ztp.window_rates(df, anchors_per_window=60)
    ori, _ = fit.find_shared_ori(windows, GENOME_LEN)
    assert fit.circular_distance(np.array([ori]), ORI, GENOME_LEN)[0] < 300_000


def test_find_shared_ori_returns_none_without_enough_data():
    windows = pd.DataFrame(
        {"enzyme": ["BcgI"] * 2, "global_mid": [1.0, 2.0], "log2_rate": [1.0, 2.0]}
    )
    ori, conf = fit.find_shared_ori(windows, GENOME_LEN)
    assert ori is None
    assert np.isnan(conf)
