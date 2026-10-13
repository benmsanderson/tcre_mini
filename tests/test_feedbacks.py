"""C4MIP reader, feedback conventions and closure (A1), and HANDOVER §2.4 reproduction."""
import numpy as np
import pytest

from tcre_evidence import c4mip, feedbacks as fb, identity as idt


@pytest.fixture(scope="module")
def t_fs():
    return fb.table("2x", "forced_start")


@pytest.fixture(scope="module")
def t_pi():
    return fb.table("2x", "picontrol")


@pytest.mark.parametrize("model", list(c4mip.MODELS))
def test_signs_and_units(model):
    runs = c4mip.read_model(model)
    c4mip.assert_signs(runs, model)
    for r in runs.values():
        assert 5 < r.tas.mean() < 25


def test_coubgc_closure_exact(t_pi):
    r = idt.reconstruct(t_pi, "coubgc")
    assert np.abs(r.af_err).max() < 1e-10
    assert np.abs(r.tcre_err).max() < 1e-10


def test_section_2_4(t_fs):
    """HANDOVER §2.4 quick computation (first-10-year baseline): TCRE 1.70 ± 0.36, RAD bias +0.008/+0.03."""
    assert t_fs.tcre.mean() == pytest.approx(1.70, abs=0.01)
    assert t_fs.tcre.std() == pytest.approx(0.36, abs=0.01)
    r = idt.reconstruct(t_fs, "bgcrad")
    assert r.af_err.mean() == pytest.approx(0.008, abs=0.002)
    assert r.tcre_err.mean() == pytest.approx(0.03, abs=0.005)
    assert r.tcre_err.idxmax() == "ACCESS-ESM1-5"
    assert t_fs.loc["ACCESS-ESM1-5", "eps_land"] == pytest.approx(116, abs=2)
    assert t_fs.loc["CNRM-ESM2-1", "eps_ocean"] == pytest.approx(49, abs=2)


def test_eps_definition(t_pi):
    for x in ("land", "ocean"):
        e = t_pi[f"dC_{x}_cou"] - t_pi[f"dC_{x}_bgc"] - t_pi[f"dC_{x}_rad"]
        assert np.allclose(e, t_pi[f"eps_{x}"])
        # exact: gamma_coubgc = (eps + gamma_rad*(dT_cou - dT_rad) + ...) ; check via C_cou - C_bgc
        lhs = t_pi[f"gamma_{x}_coubgc"] * (t_pi.dT_cou - t_pi.dT_bgc)
        assert np.allclose(lhs, t_pi[f"dC_{x}_rad"] + t_pi[f"eps_{x}"])
