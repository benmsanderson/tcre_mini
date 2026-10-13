"""The pinned upstream builders fail to construct; the patched copy constructs with the notebook's names."""
import sys

import pytest

from tcre_evidence import upstream
from tcre_evidence.paths import UPSTREAM_REPO


@pytest.fixture(scope="module")
def inp():
    return upstream.inputs()


def test_pinned_upstream_beta_land_fails(inp):
    sys.path.insert(0, str(UPSTREAM_REPO))
    try:
        from bayes import model_builders as mb
        d = inp["beta_land_data"]
        with pytest.raises(TypeError, match="var_name"):
            mb.build_beta_land_model(d["beta_land"], d["lookup_table"], upstream.priors.beta_land)
    finally:
        sys.path.remove(str(UPSTREAM_REPO))


def test_patched_builders_construct(inp):
    from tcre_evidence.upstream_patched import model_builders as mb
    d = inp["gamma_land_data"]
    m = mb.build_gamma_land_model(d["gamma_land_tropics"], d["gamma_land_extratropics"], d["lookup_table"],
                                  inp["iav_df"], inp["iav_obs_mean"], inp["iav_obs_sigma"],
                                  upstream.priors.gamma_land)
    names = {v.name for v in m.free_RVs}
    # names listed in the notebook's saved sampler output
    for n in ["gamma_land_tropics_unscaled", "gamma_land_tropics_sigma_struct", "gamma_land_tropics_rho",
              "gamma_land_tropics_bias_m_raw", "gamma_land_tropics_nitrogen", "gamma_land_tropics_eps_add",
              "negative_gamma_land_extratropics_permafrost"]:
        assert n in names, n
    b = inp["beta_land_data"]
    m = mb.build_beta_land_model(b["beta_land"], b["lookup_table"], upstream.priors.beta_land)
    assert {"beta_land_log_eta_offset", "beta_land_eps_mult"} <= {v.name for v in m.free_RVs}
