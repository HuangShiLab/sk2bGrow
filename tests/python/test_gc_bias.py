"""Per-enzyme GC correction."""

import numpy as np
import pandas as pd
import pytest

from conftest import make_counts
from sk2bgrow import gc_bias, ztp


def test_curve_recovers_a_planted_gc_slope():
    rng = np.random.default_rng(0)
    n = 3000
    gc = rng.uniform(0.35, 0.65, n)
    # log2 count rises 2.0 per unit GC -> 0.6 log2 across the 0.3-wide range.
    df = pd.DataFrame(
        {"enzyme": "BcgI", "gc": gc, "count": rng.poisson(8 * 2 ** (2.0 * (gc - 0.5))), "usable": True}
    )
    curve = gc_bias.fit_curves(df)["BcgI"]
    assert curve.amplitude == pytest.approx(0.6, abs=0.15)
    assert curve(np.array([0.36]))[0] < curve(np.array([0.64]))[0]


def test_curves_are_fitted_independently_per_enzyme():
    rng = np.random.default_rng(1)
    n = 2000
    frames = []
    for enzyme, slope in [("BcgI", 2.0), ("AlfI", -2.0)]:
        gc = rng.uniform(0.35, 0.65, n)
        frames.append(
            pd.DataFrame(
                {"enzyme": enzyme, "gc": gc, "count": rng.poisson(8 * 2 ** (slope * (gc - 0.5))), "usable": True}
            )
        )
    curves = gc_bias.fit_curves(pd.concat(frames, ignore_index=True))
    # Opposite slopes: one global curve would cancel them into nothing.
    assert curves["BcgI"](np.array([0.64]))[0] > 0
    assert curves["AlfI"](np.array([0.64]))[0] < 0


def test_thin_enzymes_get_no_curve():
    df = pd.DataFrame({"enzyme": "PpiI", "gc": [0.4, 0.5, 0.6], "count": [3, 4, 5], "usable": True})
    assert gc_bias.fit_curves(df, min_anchors=50) == {}


def test_a_curve_fitted_to_noise_shrinks_to_nothing():
    """Without this, each enzyme's loess traces its own scatter, and the
    cross-enzyme consistency test fires on an artefact of the correction."""
    rng = np.random.default_rng(10)
    for n in (200, 3000):
        gc = rng.uniform(0.35, 0.65, n)
        df = pd.DataFrame({"enzyme": "BcgI", "gc": gc, "count": rng.poisson(8, n), "usable": True})
        c = gc_bias.fit_curves(df)["BcgI"]
        assert c.shrinkage == 0.0, f"n={n}: a noise curve kept {c.shrinkage:.2f} of its amplitude"
        assert c.amplitude == 0.0


def test_a_real_gc_gradient_survives_shrinkage():
    rng = np.random.default_rng(11)
    n = 3000
    gc = rng.uniform(0.35, 0.65, n)
    df = pd.DataFrame(
        {"enzyme": "BcgI", "gc": gc, "count": rng.poisson(8 * 2 ** (2.0 * (gc - 0.5))), "usable": True}
    )
    c = gc_bias.fit_curves(df)["BcgI"]
    assert c.shrinkage > 0.9
    assert c.amplitude == pytest.approx(0.6, abs=0.15)


def test_curve_is_clamped_outside_the_fitted_range():
    rng = np.random.default_rng(2)
    gc = rng.uniform(0.45, 0.55, 500)
    df = pd.DataFrame({"enzyme": "BcgI", "gc": gc, "count": rng.poisson(8, 500), "usable": True})
    c = gc_bias.fit_curves(df)["BcgI"]
    assert c(np.array([0.05]))[0] == c.log2_offset[0]
    assert c(np.array([0.95]))[0] == c.log2_offset[-1]


def test_correction_flattens_a_biased_profile():
    """End to end: a GC-biased count table should fit better after correction."""
    df = make_counts(log2_ptr=0.0, gc_slope=3.0, n_per_enzyme=2000, enzymes=("BcgI",), sigma_eff=0.0, seed=3)
    curves = gc_bias.fit_curves(df)
    df = gc_bias.add_anchor_offsets(df, curves)
    w = ztp.window_rates(df, anchors_per_window=100)
    corrected = gc_bias.apply_to_windows(w)
    assert corrected["gc_corrected"].all()
    # With no replication gradient, the corrected rates must be flatter.
    assert corrected["log2_rate"].std() < corrected["log2_rate_raw"].std()


def test_missing_offsets_are_marked_not_silently_skipped():
    w = pd.DataFrame({"log2_rate": [1.0, 2.0]})
    out = gc_bias.apply_to_windows(w)
    assert not out["gc_corrected"].any()
    assert (out["log2_rate"] == out["log2_rate_raw"]).all()


def test_anchors_of_uncorrected_enzymes_get_zero_offset():
    # A real GC gradient, so BcgI's curve survives shrinkage and the contrast
    # with the uncorrected enzyme is meaningful.
    df = make_counts(n_per_enzyme=2000, enzymes=("BcgI", "PpiI"), gc_slope=3.0, sigma_eff=0.0, seed=4)
    curves = {k: v for k, v in gc_bias.fit_curves(df).items() if k == "BcgI"}
    assert curves["BcgI"].shrinkage > 0.5, "fixture no longer carries a detectable GC bias"
    out = gc_bias.add_anchor_offsets(df, curves)
    assert (out.loc[out["enzyme"] == "PpiI", "gc_offset"] == 0).all()
    assert out.loc[out["enzyme"] == "BcgI", "gc_offset"].abs().sum() > 0


def test_tukey_mask_fences_outliers():
    v = np.array([1.0, 2.0, 3.0, 4.0, 100.0])
    m = gc_bias.tukey_mask(v)
    assert m.sum() == 4
    assert not m[-1]
    # A degenerate distribution keeps everything rather than rejecting all of it.
    assert gc_bias.tukey_mask(np.full(10, 5.0)).all()
    assert gc_bias.tukey_mask(np.array([1.0, np.nan, 3.0])).sum() == 2


def test_enzyme_efficiency_is_relative_to_the_panel():
    df = make_counts(n_per_enzyme=500, enzymes=("BcgI", "AlfI", "CjeI"), seed=5)
    df.loc[df["enzyme"] == "AlfI", "count"] = (df.loc[df["enzyme"] == "AlfI", "count"] * 0.25).astype(int)
    eff = gc_bias.enzyme_efficiency(df)
    alfi = float(eff.loc[eff["enzyme"] == "AlfI", "rel_efficiency"].iloc[0])
    assert alfi < 0.4
    assert float(eff.loc[eff["enzyme"] == "BcgI", "rel_efficiency"].iloc[0]) == pytest.approx(1.0, abs=0.15)
