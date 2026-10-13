"""Regression tests for HANDOVER §2.3 and the upstream alpha alignment."""
import numpy as np
import pytest

from tcre_evidence import covariance as cv
from tcre_evidence import upstream_data as ud


@pytest.fixture(scope="module")
def hardcoded():
    cp = ud.carbon_parameters("2xCO2")
    return {"alpha": ud.HARDCODED_ALPHA, **{c: cp[c].values for c in cv.COMPONENTS[1:]}}


def test_common_models_order():
    cp = ud.carbon_parameters("2xCO2")
    assert list(cp.model.values) == ["ACCESS-ESM1-5", "CESM2", "CMCC-ESM2", "CanESM5", "EC-Earth3-CC",
                                     "GFDL-ESM4", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR",
                                     "NorESM2-LM", "UKESM1-0-LL"]


def test_section_2_3_correlations(hardcoded):
    c = cv.pair_correlations(hardcoded)
    target = {"beta_L,gamma_L": -0.38, "beta,alpha*gamma": -0.29, "alpha,alpha*gamma": -0.53,
              "beta_O,gamma_O": -0.56, "alpha,beta_O": 0.52}
    for k, v in target.items():
        assert round(c[k], 2) == v, k


@pytest.mark.parametrize("scheme,sd_af,mean,sd,p05,p95", [
    ("paired", 0.036, 1.84, 0.475, 1.31, 2.85),
    ("upstream", 0.044, 1.84, 0.443, 1.26, 2.59),
    ("independent", 0.044, 1.84, 0.442, 1.26, 2.59),
])
def test_section_2_3_table(hardcoded, scheme, sd_af, mean, sd, p05, p95):
    s = cv.summary(cv.resample(hardcoded, scheme, n=400_000, seed=0))
    assert s["sd_af"] == pytest.approx(sd_af, abs=1.5e-3)
    assert s["tcre_mean"] == pytest.approx(mean, abs=0.01)
    assert s["tcre_sd"] == pytest.approx(sd, abs=3e-3)
    assert s["tcre_p05"] == pytest.approx(p05, abs=0.015)
    assert s["tcre_p95"] == pytest.approx(p95, abs=0.015)


def test_hardcoded_alpha_not_aligned_with_common_models():
    """Upstream issue 8: the hard-coded alpha array does not follow intersect1d order."""
    cp = ud.carbon_parameters("2xCO2")
    by_name = ud.cmip_sheet_dataset().alpha.reindex(model=cp.model).values
    ok = ~np.isnan(by_name)
    assert np.isnan(by_name).sum() == 2  # CMCC-ESM2 and EC-Earth3-CC have no alpha in the sheet
    mismatch = np.abs(ud.HARDCODED_ALPHA[ok] - by_name[ok]) > 2e-4
    assert mismatch.sum() == 6  # CESM2, CanESM5, GFDL, IPSL, MIROC, MPI


def test_identity_af_tcre():
    af, t = cv.tcre_from_components(1.7 / 280, 1.0, -40.0, 0.8, -15.0)
    assert t == pytest.approx(1e3 * (1.7 / 280) * af / ud.K)
