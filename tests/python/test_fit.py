"""PTR estimators: coordinate V-shape fit and the Pilea-parity sorted regression."""

import numpy as np
import pytest

from conftest import GENOME_LEN, ORI
from sk2bgrow import fit


def v_profile(log2_ptr, n=300, noise=0.05, seed=0, ori=ORI, kink=None, second_slope=None):
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.uniform(0, GENOME_LEN, n))
    d = fit.circular_distance(pos, ori, GENOME_LEN)
    half = GENOME_LEN / 2
    if kink is None:
        y = 6.0 - log2_ptr * d / half
    else:
        b1 = log2_ptr / half
        b2 = second_slope if second_slope is not None else b1 / 3
        y = 6.0 - b1 * np.minimum(d, kink) - b2 * np.maximum(d - kink, 0)
    return pos, y + rng.normal(0, noise, n)


def test_circular_distance_wraps_and_peaks_at_the_terminus():
    g = 1000.0
    assert fit.circular_distance(np.array([10.0]), 990.0, g)[0] == 20.0
    assert fit.circular_distance(np.array([500.0]), 0.0, g)[0] == 500.0
    assert fit.circular_distance(np.array([0.0]), 0.0, g)[0] == 0.0


@pytest.mark.parametrize("truth", [0.3, 1.0, 2.0])
def test_v_shape_recovers_ptr_and_ori(truth):
    pos, y = v_profile(truth, seed=1)
    f = fit.fit_v_shape(pos, y, GENOME_LEN, se=np.full(pos.size, 0.05))
    assert f.log2_ptr == pytest.approx(truth, abs=0.08)
    assert fit.circular_distance(np.array([f.ori]), ORI, GENOME_LEN)[0] < 120_000
    assert f.ok
    assert f.ptr == pytest.approx(2**truth, rel=0.1)


def test_known_ori_is_used_verbatim():
    pos, y = v_profile(1.0, seed=2)
    f = fit.fit_v_shape(pos, y, GENOME_LEN, se=np.full(pos.size, 0.05), ori=ORI)
    assert f.ori == ORI
    assert f.ori_confidence == 1.0
    assert f.log2_ptr == pytest.approx(1.0, abs=0.06)


def test_flat_profile_gives_low_ori_confidence():
    rng = np.random.default_rng(3)
    pos = np.sort(rng.uniform(0, GENOME_LEN, 200))
    y = 6.0 + rng.normal(0, 0.05, 200)
    f = fit.fit_v_shape(pos, y, GENOME_LEN, se=np.full(200, 0.05))
    assert f.log2_ptr < 0.15, "a flat profile must not imply growth"
    assert f.ori_confidence < 0.95


def test_uphill_slope_is_rejected_rather_than_reported_as_the_ori():
    """The terminus fits the same line read backwards; it must never win."""
    pos, y = v_profile(1.5, seed=4)
    f = fit.fit_v_shape(pos, y, GENOME_LEN, se=np.full(pos.size, 0.05))
    ter = (ORI + GENOME_LEN // 2) % GENOME_LEN
    assert fit.circular_distance(np.array([f.ori]), ter, GENOME_LEN)[0] > GENOME_LEN / 4
    assert all(s >= 0 for s in f.slopes)


def test_segmented_model_is_chosen_only_when_there_is_a_kink():
    pos, y = v_profile(1.0, n=400, noise=0.03, seed=5)
    plain = fit.fit_v_shape(pos, y, GENOME_LEN, se=np.full(pos.size, 0.03), ori=ORI)
    assert not plain.segmented, "BIC picked the 4-parameter model on a plain V"

    kpos, ky = v_profile(2.0, n=400, noise=0.03, seed=6, kink=GENOME_LEN * 0.25, second_slope=0.0)
    kinked = fit.fit_v_shape(kpos, ky, GENOME_LEN, se=np.full(kpos.size, 0.03), ori=ORI)
    assert kinked.segmented, "a hard kink should select the segmented model"


def test_standard_error_reflects_scatter_not_just_the_input_errors():
    """Window errors describe counting noise only; residual scatter must widen them."""
    pos, clean = v_profile(1.0, n=300, noise=0.01, seed=7)
    pos2, noisy = v_profile(1.0, n=300, noise=0.30, seed=7)
    se_in = np.full(300, 0.01)
    a = fit.fit_v_shape(pos, clean, GENOME_LEN, se=se_in, ori=ORI)
    b = fit.fit_v_shape(pos2, noisy, GENOME_LEN, se=se_in, ori=ORI)
    assert b.se > 5 * a.se
    assert b.reduced_chi2 > 10 * a.reduced_chi2


def test_too_few_windows_is_reported_not_guessed():
    f = fit.fit_v_shape(np.array([0.0, 1.0, 2.0]), np.array([1.0, 2.0, 3.0]), GENOME_LEN)
    assert not f.ok
    assert np.isnan(f.log2_ptr)
    assert "5" in f.note


def test_sorted_ransac_runs_without_coordinates():
    _, y = v_profile(1.0, n=200, noise=0.05, seed=8)
    f = fit.fit_sorted_ransac(y)
    assert f.method == "sorted_ransac"
    # The trimmed sorted estimator is known to be biased low; the report's D3
    # critique is exactly this. Bound it rather than pretend it is unbiased.
    assert 0.6 < f.log2_ptr < 1.2


def test_sorted_ransac_survives_extreme_outliers():
    _, y = v_profile(1.0, n=200, noise=0.05, seed=9)
    clean = fit.fit_sorted_ransac(y)
    dirty = y.copy()
    dirty[:3] = -20.0
    dirty[-3:] = 20.0
    robust = fit.fit_sorted_ransac(dirty)
    assert abs(robust.log2_ptr - clean.log2_ptr) < 0.6


def test_fit_windows_falls_back_on_fragmented_references(manifest, fragmented_manifest):
    import pandas as pd

    pos, y = v_profile(1.0, n=120, seed=10)
    windows = pd.DataFrame(
        {
            "sample": "S1", "genome_id": 0, "genome": "ecoli", "enzyme": "BcgI",
            "global_mid": pos, "log2_rate": y, "log2_se": 0.05, "rate": 2**y,
            "detected_fraction": 1.0, "dispersion": 1.0, "n_anchors": 100,
        }
    )
    good = fit.fit_windows(windows, manifest, method="auto")
    assert good["method"].iloc[0].startswith("v_shape")
    frag = fit.fit_windows(windows, fragmented_manifest, method="auto")
    assert frag["method"].iloc[0] == "sorted_ransac"
    assert "fragmented" in frag["note"].iloc[0]
