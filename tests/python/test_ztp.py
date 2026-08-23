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
