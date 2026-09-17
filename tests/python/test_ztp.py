"""Window rate model: zero-truncated Poisson mixtures and the NB branch."""

import numpy as np
import pytest

from conftest import make_counts
from sk2bgrow import ztp


def test_ztp_mean_and_inverse_are_consistent():
    for lam in [0.05, 0.5, 1.0, 3.0, 12.0, 60.0]:
        assert ztp.solve_ztp_lambda(ztp.ztp_mean(lam)) == pytest.approx(lam, rel=1e-6)


def test_ztp_mean_tends_to_one_at_zero_rate():
    # A truncated Poisson with a vanishing rate almost surely yields exactly 1.
    assert ztp.ztp_mean(1e-9) == pytest.approx(1.0, abs=1e-6)
    # A sample mean at or below 1 carries no information about lambda.
    assert ztp.solve_ztp_lambda(1.0) == 0.0
    assert ztp.solve_ztp_lambda(0.4) == 0.0


def test_ztp_variance_matches_simulation():
    rng = np.random.default_rng(0)
    for lam in [0.8, 4.0]:
        x = rng.poisson(lam, 400_000)
        x = x[x >= 1]
        assert x.mean() == pytest.approx(ztp.ztp_mean(lam), rel=0.01)
        assert x.var() == pytest.approx(ztp.ztp_var(lam), rel=0.03)


@pytest.mark.parametrize("true_lam", [0.8, 3.0, 12.0])
def test_rate_recovered_from_truncated_data(true_lam):
    """The point of truncation handling: a plain mean of positives is biased up."""
    rng = np.random.default_rng(1)
    x = rng.poisson(true_lam, 200_000)
    x = x[x >= 1][:20_000]
    est = ztp.estimate_window_rate(x, model="ztp")
    assert est.rate == pytest.approx(true_lam, rel=0.05)
    assert est.model == "ztp"
    if true_lam < 2:
        assert x.mean() > true_lam * 1.15, "fixture should be visibly truncation-biased"


def test_standard_error_shrinks_as_sqrt_n():
    rng = np.random.default_rng(2)
    x = rng.poisson(5.0, 100_000)
    x = x[x >= 1]
    small = ztp.estimate_window_rate(x[:200], model="ztp")
    large = ztp.estimate_window_rate(x[:20_000], model="ztp")
    assert small.se / large.se == pytest.approx(10.0, rel=0.35)


def test_log2_standard_error_uses_the_delta_method():
    est = ztp.estimate_window_rate(np.full(500, 5), model="ztp")
    assert est.log2_se == pytest.approx(est.se / (est.rate * np.log(2)), rel=1e-9)


def test_negative_binomial_wins_on_overdispersed_data():
    rng = np.random.default_rng(3)
    lam = rng.gamma(shape=1.5, scale=6.0, size=6_000)  # gamma-Poisson
    x = rng.poisson(lam)
    x = x[x >= 1]
    est = ztp.estimate_window_rate(x, model="auto")
    assert est.model == "ztnb", "auto should switch to NB when dispersion is real"
    assert est.dispersion > 2.0


def test_poisson_data_stays_on_ztp():
    rng = np.random.default_rng(4)
    x = rng.poisson(6.0, 3_000)
    assert ztp.estimate_window_rate(x[x >= 1], model="auto").model == "ztp"


def test_empty_window_is_reported_not_crashed():
    est = ztp.estimate_window_rate(np.zeros(50), model="auto")
    assert est.model == "empty"
    assert np.isnan(est.rate)
    assert est.detected_fraction == 0.0


def test_all_ones_window_is_uninformative_not_a_garbage_rate():
    """C7 regression: an all-ones window carries no rate information.

    The ZTP likelihood has a finite supremum at the boundary (loglik -> 0 as
    lam -> 0+), so BIC must keep the ZTP branch and the window must come back
    NaN. Before the fix the ZTP loglik was hard-coded to -inf at lam == 0, BIC
    handed the window to ZTNB, and ZTNB returned a finite but meaningless
    rate of ~1e-6 (log2 ~ -20) that the downstream Tukey fence had to catch.
    """
    counts = np.ones(31)
    m = ztp.fit_ztp_mixture(counts)
    assert m.lambdas[0] == 0.0
    assert m.loglik == pytest.approx(0.0, abs=1e-12)
    assert np.isfinite(m.bic)
    est = ztp.estimate_window_rate(counts, model="auto")
    assert est.model == "ztp", "BIC must not hand an all-ones window to ZTNB"
    assert np.isnan(est.rate)


def test_near_boundary_window_with_a_single_two_still_fits():
    # 29 ones and one 2: the truncated mean is above 1, so the ordinary ZTP
    # path applies and must yield a small but finite rate.
    counts = np.ones(30)
    counts[0] = 2
    est = ztp.estimate_window_rate(counts, model="auto")
    assert est.model == "ztp"
    assert 0 < est.rate < 0.2


def test_ztp_logpmf_is_accurate_at_the_zero_boundary():
    """logpmf(1) must tend to 0 as lam -> 0+, never positive (a probability
    above one). The naive log1p(-exp(-lam)) form lost the mantissa here and
    handed EM phantom likelihood on mostly-ones windows (C7)."""
    for lam in [1e-12, 1e-9, 1e-6]:
        assert ztp.ztp_logpmf(np.array([1.0]), lam)[0] == pytest.approx(0.0, abs=1e-6)
        assert ztp.ztp_logpmf(np.array([2.0]), lam)[0] < -5.0
    # Away from the boundary the expm1 form matches the plain one.
    from scipy import special

    k = np.arange(1, 8, dtype=float)
    for lam in [0.05, 0.5, 3.0]:
        plain = k * np.log(lam) - lam - special.gammaln(k + 1.0) - np.log1p(-np.exp(-lam))
        assert ztp.ztp_logpmf(k, lam) == pytest.approx(plain, rel=1e-12)


def test_mixture_does_not_collapse_a_component_to_zero_on_mostly_ones():
    """34 ones and one 2: a phantom logpmf at lam ~ 0 used to reward a
    two-component fit whose dominant component sat at lam ~ 1e-12, returning
    a garbage rate. With accurate boundary numerics the single component wins.
    """
    counts = np.concatenate([np.ones(34), [2.0]])
    est = ztp.estimate_window_rate(counts, model="auto")
    assert est.n_components == 1
    assert est.rate == pytest.approx(ztp.solve_ztp_lambda(36 / 35), rel=1e-6)


def test_ztnb_ridge_fit_falls_back_to_ztp():
    """C7 regression: mostly-one counts with a heavy tail send the ZTNB MLE
    onto the non-identifiable (mu -> 0, alpha -> inf) ridge, whose implied
    detection probability contradicts the observed detected fraction. The
    window must not report the ridge's garbage mu; it falls back to ZTP.
    """
    # 12 zeros, then 1x8, 2x2, 3x2, 6: ZTNB's MLE is mu ~ 0.002 at alpha ~ 1e3,
    # implying P(X >= 1) ~ 0.001 while 13 of 25 anchors were observed.
    counts = np.array([0] * 12 + [1] * 8 + [2] * 2 + [3] * 2 + [6], dtype=float)
    est = ztp.estimate_window_rate(counts, model="auto")
    assert est.model == "ztp"
    assert est.rate == pytest.approx(ztp.solve_ztp_lambda(24 / 13), rel=1e-6)


def test_mixture_separates_two_components():
    rng = np.random.default_rng(5)
    x = np.concatenate([rng.poisson(2.0, 3_000), rng.poisson(20.0, 3_000)])
    x = x[x >= 1]
    m = ztp.fit_ztp_mixture(x, n_components=2)
    lo, hi = np.sort(m.lambdas)
    assert lo == pytest.approx(2.0, rel=0.2)
    assert hi == pytest.approx(20.0, rel=0.2)
    assert m.bic < ztp.fit_ztp_mixture(x, n_components=1).bic


def test_windows_never_span_contigs():
    df = make_counts(n_per_enzyme=400, enzymes=("BcgI",), seed=7)
    df.loc[df.index[200:], "contig_id"] = 1
    w = ztp.window_rates(df, anchors_per_window=50)
    assert w["contig_id"].nunique() == 2
    assert len(w) >= 8


def test_window_rates_are_per_enzyme():
    df = make_counts(n_per_enzyme=400, enzymes=("BcgI", "AlfI"), seed=8)
    w = ztp.window_rates(df, anchors_per_window=100)
    assert set(w["enzyme"]) == {"BcgI", "AlfI"}
    # Each enzyme is windowed inside its own anchor series, so both get windows.
    assert (w.groupby("enzyme").size() >= 3).all()
