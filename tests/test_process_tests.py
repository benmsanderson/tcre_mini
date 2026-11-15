"""process_tests on the CMIP6 example tables: input handling, the least-squares tests, and (slow) one fit."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from tcre_evidence import process_tests as pt
from tcre_evidence.paths import DATA

VALUES = DATA / "examples" / "cmip6_land_feedbacks.csv"
FLAGS = DATA / "examples" / "cmip6_process_flags.csv"


@pytest.fixture(scope="module")
def gamma():
    return pt.load_ensemble(VALUES, FLAGS, "gamma_land_coubgc")


def test_load_ensemble(gamma):
    y, flags, groups = gamma
    assert len(y) == 13 and list(flags.columns) == ["nitrogen", "permafrost", "fire", "veg"]
    assert set(flags.values.ravel()) <= {0, 1}
    assert groups["CESM2"] == groups["NorESM2-LM"] == "NCAR"


def test_load_rejects_non_binary_flags(tmp_path):
    F = pd.read_csv(FLAGS)
    F.loc[0, "fire"] = 2
    F.to_csv(tmp_path / "flags.csv", index=False)
    with pytest.raises(ValueError, match="0/1"):
        pt.load_ensemble(VALUES, tmp_path / "flags.csv", "gamma_land_coubgc")


def test_flag_budget_matches_note(gamma):
    y, flags, _ = gamma
    b = pt.flag_budget(y, flags)
    assert b["r2"] == pytest.approx(0.202, abs=1e-3)
    assert b["loo_q2"] == pytest.approx(-0.697, abs=1e-3)
    assert b["resid_sd"] == pytest.approx(28.70, abs=0.01)


def test_detection_power_matches_m8(gamma):
    y, flags, _ = gamma
    # the balanced-design search draws random designs, so pass the effects in m8's order
    P = pt.detection_power(y, flags, {"nitrogen": 10.0, "permafrost": 30.0, "fire": 5.0, "veg": 5.0}).iloc[1]
    assert P.se_now == pytest.approx(34.52, abs=0.01)
    assert P.n_needed_same_mix == 78 and P.n_needed_balanced == 20


def test_default_priors_new_process():
    pr = pt.default_priors(pt.ADDITIVE, ["nitrogen", "demography"])
    assert pr["delta"]["demography"] == dict(dist="Normal", mu=0.0, sigma=10.0)
    assert pr["delta"]["nitrogen"]["dist"] == "LogNormal"   # the development code's own prior
    pm_ = pt.default_priors(pt.MULTIPLICATIVE, ["demography"])
    assert pm_["eta"]["demography"] == dict(mu=0.0, sigma=0.3)


def test_synthetic_ensemble_known_effects(gamma):
    _, flags, _ = gamma
    y = pt.synthetic_ensemble(flags, pt.ADDITIVE, -15.0, {"permafrost": -30.0}, scatter_sd=0.0)
    assert np.allclose(y, -15.0 - 30.0 * flags.permafrost)
    b = pt.synthetic_ensemble(flags, pt.MULTIPLICATIVE, 1.2, {"nitrogen": 0.7}, scatter_sd=0.0)
    assert np.allclose(b, 1.2 * 0.7 ** flags.nitrogen)


@pytest.mark.slow
def test_fit_and_own_flags_predictive(gamma):
    _, flags, _ = gamma
    priors = pt.default_priors(pt.ADDITIVE, list(flags.columns))
    y = pt.synthetic_ensemble(flags, pt.ADDITIVE, -15.0, {"nitrogen": 10.0}, scatter_sd=5.0,
                              rng=np.random.default_rng(1))
    idata, diag = pt.fit(pt.build(y, flags, priors), sampler=dict(draws=300, tune=300, chains=2))
    assert diag["rhat_max"] < 1.05
    pred = pt.own_flags_predictive(idata, flags.loc["CESM2"], priors)
    assert pred.shape == (600,) and np.isfinite(pred).all()


@pytest.mark.parametrize("kind", [pt.ADDITIVE, pt.MULTIPLICATIVE])
def test_build_paired_spread_modes(gamma, kind):
    y, flags, _ = gamma
    pr = pt.default_priors(kind, list(flags.columns))
    vals = {"fire": [5.0, -3.0]} if kind == pt.ADDITIVE else {"fire": [0.9, 1.05]}
    m0 = pt.build(y, flags, pr, paired=vals)
    assert f"{pt.VAR}_fire_1" not in m0.named_vars          # upstream label: delta_fire_1 / eta_fire_1
    m1 = pt.build(y, flags, pr, paired=vals, paired_spread="struct")
    assert f"{pt.VAR}_paired_fire" in m1.named_vars and f"{pt.VAR}_paired_tau_fire" not in m1.named_vars
    m2 = pt.build(y, flags, pr, paired=vals, paired_spread="estimate")
    assert f"{pt.VAR}_paired_tau_fire" in m2.named_vars
    obs = m2[f"{pt.VAR}_paired_fire"]
    expect = vals["fire"] if kind == pt.ADDITIVE else np.log(vals["fire"])
    assert np.allclose(m2.rvs_to_values[obs].data, expect)
    with pytest.raises(ValueError, match="paired_spread"):
        pt.build(y, flags, pr, paired=vals, paired_spread="fixed")
    with pytest.raises(ValueError, match="not a flag column"):
        pt.build(y, flags, pr, paired={"demography": [1.0]}, paired_spread="struct")


def test_paired_struct_observation_sd(gamma):
    """In "struct" mode the paired likelihood sd is sqrt(sigma_process_struct^2 + paired_sigma^2)."""
    y, flags, _ = gamma
    pr = pt.default_priors(pt.ADDITIVE, list(flags.columns))
    m = pt.build(y, flags, pr, paired={"fire": [4.0]}, paired_spread="struct", paired_sigma=2.0)
    pt_ = m.compile_logp(vars=[m[f"{pt.VAR}_paired_fire"]], sum=True)
    point = m.initial_point()
    d = float(point[f"{pt.VAR}_fire"])                    # fire's prior is Normal, so its value is untransformed
    sd = np.hypot(pr["sigma_process_struct"]["fire"], 2.0)
    assert float(pt_(point)) == pytest.approx(stats.norm.logpdf(4.0, d, sd), rel=1e-6)
