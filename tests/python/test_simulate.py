"""Monte-Carlo harness: does it reproduce the report's section 5 conclusions?"""

import numpy as np
import pytest

from sk2bgrow import simulate


def test_replication_profile_is_flat_at_zero_ptr():
    pos = np.linspace(0, simulate.ECOLI_LEN, 100)
    y = simulate.replication_log2_copy(pos, simulate.ECOLI_ORI, simulate.ECOLI_LEN, 0.0)
    assert np.allclose(y, 0.0)


def test_replication_profile_drops_by_log2_ptr_at_the_terminus():
    g, ori = simulate.ECOLI_LEN, simulate.ECOLI_ORI
    ter = (ori + g // 2) % g
    y = simulate.replication_log2_copy(np.array([ori, ter]), ori, g, 1.5)
    assert y[0] == pytest.approx(0.0)
    assert y[1] == pytest.approx(-1.5, abs=1e-3)


def test_efficiency_noise_keeps_the_mean_depth():
    """Changing sigma_eff must change the noise, not the depth."""
    pos = simulate.synthetic_anchors(20_000, simulate.ECOLI_LEN, "even")
    quiet = simulate.SimConfig(log2_ptr=0.0, per_site_depth=10.0, sigma_eff=0.0)
    noisy = simulate.SimConfig(log2_ptr=0.0, per_site_depth=10.0, sigma_eff=0.6)
    a = simulate.simulate_counts(pos, quiet, np.random.default_rng(0))
    b = simulate.simulate_counts(pos, noisy, np.random.default_rng(0))
    assert a.mean() == pytest.approx(10.0, rel=0.03)
    assert b.mean() == pytest.approx(10.0, rel=0.05)
    assert b.var() > 2 * a.var()


def test_estimator_recovers_the_truth_at_high_depth():
    pos = simulate.synthetic_anchors(28_381, simulate.ECOLI_LEN, "even")
    cfg = simulate.SimConfig(log2_ptr=1.0, per_site_depth=20.0, sigma_eff=0.05)
    counts = simulate.simulate_counts(pos, cfg, np.random.default_rng(1))
    assert simulate.estimate_sorted(pos, counts, cfg) == pytest.approx(1.0, abs=0.15)
    assert simulate.estimate(pos, counts, cfg, "v_shape") == pytest.approx(1.0, abs=0.15)


@pytest.mark.parametrize("coverage", [1.0, 5.0])
def test_route_a_ordering_matches_the_report(coverage):
    """Union beats a random sketch beats a single enzyme, at every coverage.

    This is the report's central route-A claim (section 5.2). Absolute RMSE
    depends on using the real digested coordinates; the *ordering* is the
    reproducible part and is what this asserts.
    """
    df = simulate.route_a(coverages=[coverage], log2_ptrs=[1.0], n_reps=30, seed=0)
    rmse = df.set_index("anchor_set")["rmse"]
    assert rmse["union_16"] < rmse["random_sketch"] < rmse["BcgI_single"]


def test_route_a_improves_monotonically_with_coverage():
    df = simulate.route_a(coverages=[1.0, 2.0, 5.0, 10.0], log2_ptrs=[1.0], n_reps=25, seed=1)
    for name, grp in df.groupby("anchor_set"):
        r = grp.sort_values("coverage")["rmse"].to_numpy()
        assert r[-1] < r[0], f"{name} did not improve with coverage"


def test_single_enzyme_collapses_at_low_coverage():
    df = simulate.route_a(coverages=[1.0], log2_ptrs=[1.0], n_reps=30, seed=2)
    single = df[df["anchor_set"] == "BcgI_single"].iloc[0]
    union = df[df["anchor_set"] == "union_16"].iloc[0]
    assert single["rmse"] > 4 * union["rmse"], "the single-enzyme failure mode did not reproduce"


def test_route_b_union_wins_at_equal_read_budget():
    """The wet-lab claim: spreading a fixed budget over 16 enzymes is a net win.

    Averaging within a window quenches independent site-efficiency noise at rate
    sqrt(n), and the union has ~10x more anchors per window, so trading per-site
    depth for window replication pays off whenever efficiency noise dominates.
    """
    df = simulate.route_b(sigma_effs=[0.6], log2_ptrs=[1.0], n_reps=30, seed=3)
    rmse = df.set_index("anchor_set")["rmse"]
    assert rmse["union_16_same_budget"] < rmse["BcgI_single"]
    assert rmse["union_16_full_budget"] <= rmse["union_16_same_budget"] * 1.2


def test_route_b_gap_widens_with_efficiency_noise():
    df = simulate.route_b(sigma_effs=[0.3, 0.6], log2_ptrs=[1.0], n_reps=30, seed=4)
    ratios = {}
    for sig, grp in df.groupby("sigma_eff"):
        r = grp.set_index("anchor_set")["rmse"]
        ratios[sig] = r["BcgI_single"] / r["union_16_same_budget"]
    assert ratios[0.6] > ratios[0.3], "the union's advantage should grow with sigma_eff"


def test_anchors_from_digest_reads_a_tgt_dump(tmp_path):
    p = tmp_path / "g.tgt"
    p.write_text(
        "#TGT2\tg\n"
        "#contig\t0\tc0\t1000\t0\tChromosome\n"
        "#columns\ttag\tcontig_id\tposition\tstrand\tgap\n"
        "BcgI:ACGT\t0\t100\t+\t0\t3\t100\n"
        "BcgI:TTTT\t0\t400\t+\t300\t3\t100\n"
        "AlfI:GGGG\t0\t250\t-\t0\t3\t100\n"
    )
    assert list(simulate.anchors_from_digest(p)) == [100, 250, 400]
    assert list(simulate.anchors_from_digest(p, enzyme="BcgI")) == [100, 400]


def test_synthetic_anchor_kinds_differ_in_spacing():
    n, g = 5_000, 1_000_000
    even = simulate.synthetic_anchors(n, g, "even")
    uniform = simulate.synthetic_anchors(n, g, "uniform", seed=0)
    # Bernoulli sampling gives geometric spacings with a long tail; regular
    # spacing has none. This is defect D2 in one line.
    assert np.diff(uniform).max() > 5 * np.diff(even).max()
    with pytest.raises(ValueError):
        simulate.synthetic_anchors(10, g, "spiral")
