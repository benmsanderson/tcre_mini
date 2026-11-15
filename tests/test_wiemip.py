"""WIEMIP feedbacks on a synthetic ensemble with known beta, gamma and process effects, and on the real
extraction when data/wiemip/ is present (it exists only on the WIEMIP hub)."""
import numpy as np
import pandas as pd
import pytest

from tcre_evidence import wiemip as wm

YEARS = np.arange(1850, 2001)


def _forcing():
    co2 = 280 * 1.01 ** (YEARS - 1850)
    return pd.DataFrame(dict(co2_ppm=co2, gmst_k=0.006 * (YEARS - 1850)), index=pd.Index(YEARS, name="year"))


def _ensemble(fc):
    """cland = 1000 + beta*dCa + gamma*dT (cou), beta*dCa (bgc), gamma_rad*dT (rad); ctrl drifts linearly.

    Model A: baseline beta 1.0, gamma -40; noNitrogen beta 2.0, gamma -60 (eta 0.5, delta +20).
    """
    dCa, dT = fc.co2_ppm - 280, fc.gmst_k
    drift = 0.1 * (YEARS - 1850)
    spec = {("A", "baseline"): (1.0, -40.0, -50.0), ("A", "noNitrogen"): (2.0, -60.0, -70.0),
            ("B", "baseline"): (0.5, -10.0, -12.0)}
    rows = []
    for (m, fac), (b, g, gr) in spec.items():
        runs = {("ctrl", "stable"): 0 * dCa, ("bgc", "stable"): b * dCa, ("cou", "ukesm"): b * dCa + g * dT,
                ("cou", "gfdl"): b * dCa + (g - 10) * dT, ("rad", "ukesm"): gr * dT}
        for (sim, frc), c in runs.items():
            rows += [(m, fac, sim, frc, y, 1000 + d + v) for y, d, v in zip(YEARS, drift, c)]
    return pd.DataFrame(rows, columns=["model", "factorial", "sim", "forcing", "year", "cland_pgc"])


@pytest.fixture(scope="module")
def synth():
    fc = _forcing()
    cl = _ensemble(fc)
    return cl, fc, wm.feedbacks(cl, fc, "2x")


def test_processes():
    assert wm.processes("baseline") == frozenset()
    assert wm.processes("noFire_noNitrogen") == {"fire", "nitrogen"}
    assert wm.processes("noBVOC") == {"bvoc"}
    with pytest.raises(ValueError):
        wm.processes("fire_off")


def test_recovers_beta_gamma(synth):
    _, _, fb = synth
    r = fb.set_index(["model", "factorial", "pattern"])
    assert r.loc[("A", "baseline", "ukesm")].beta_land_coubgc == pytest.approx(1.0)
    assert r.loc[("A", "baseline", "ukesm")].gamma_land_coubgc == pytest.approx(-40.0)
    assert r.loc[("A", "baseline", "gfdl")].gamma_land_coubgc == pytest.approx(-50.0)
    assert r.loc[("A", "baseline", "ukesm")].gamma_land_rad == pytest.approx(-50.0)
    assert r.loc[("A", "baseline", "ukesm")].eps_land == pytest.approx(10.0 * r.loc[("A", "baseline", "ukesm")].dT)
    assert np.isnan(r.loc[("A", "baseline", "gfdl")].gamma_land_rad)


def test_values_and_paired(synth):
    _, _, fb = synth
    v = wm.values_table(fb, "mean").set_index("model")
    assert v.loc["A", "gamma_land_coubgc"] == pytest.approx(-45.0)
    assert v.loc["A", "n_patterns"] == 2
    pe = wm.paired_effects(fb).query("pattern == 'ukesm'").iloc[0]
    assert (pe.model, pe.process, pe.background) == ("A", "nitrogen", "")
    assert pe.eta_beta_land == pytest.approx(0.5)
    assert pe.delta_gamma_land_coubgc == pytest.approx(20.0)


def test_screen_drops_nonfinite(synth):
    cl, fc, _ = synth
    cl = cl.copy()
    cl.loc[(cl.model == "B") & (cl.sim == "bgc") & (cl.year == 1915), "cland_pgc"] = np.inf
    with pytest.warns(UserWarning, match="B baseline bgc/stable dropped"):
        fb = wm.feedbacks(cl, fc, "2x")
    assert "B" not in set(fb.model)


def test_sign_assertion(synth):
    cl, fc, _ = synth
    cl = cl.copy()
    m = (cl.model == "B") & (cl.sim == "bgc")
    cl.loc[m, "cland_pgc"] = 2000 - cl.loc[m, "cland_pgc"]
    with pytest.raises(AssertionError, match="BGC land uptake"):
        wm.feedbacks(cl, fc, "2x")


@pytest.mark.skipif(not (wm.WIEMIP_DATA / "cland_annual.csv").exists(), reason="WIEMIP extraction not present")
def test_real_extraction():
    cl, fc = wm.load()
    fb = wm.feedbacks(cl, fc, "2x")
    assert fb.dT.iloc[0] == pytest.approx(1.742, abs=1e-3)
    assert fb.dCa.iloc[0] == pytest.approx(280.0, abs=0.5)
    assert set(wm.values_table(fb, "ukesm").model) >= {"CLM", "CLM_FATES", "JULES", "LPJ_GUESS"}


def test_factorial_flags(synth):
    cl, _, _ = synth
    f = wm.factorial_flags(cl)
    assert f.loc["A", "nitrogen"] == 1 and np.isnan(f.loc["B", "nitrogen"])


def test_design_se_factor_balanced():
    from tcre_evidence import process_tests as pt
    # one flag, half the models flagged: se factor = sqrt(1/k + 1/(n-k)) = 2/sqrt(n)
    assert pt.design_se_factor(16, 8, 0) == pytest.approx(0.5)
    assert pt.design_se_factor(18, 3, 3) > pt.design_se_factor(18, 9, 3)


def test_paired_by_model_one_value_per_model():
    pe = pd.DataFrame(dict(model=["A", "A", "A", "B"], pattern=["ukesm", "gfdl", "ipsl", "ukesm"],
                           process="fire", background="", eta_beta_land=[0.5, 0.5, 2.0, 0.8],
                           delta_gamma_land_coubgc=[1.0, 2.0, 3.0, np.nan]))
    assert wm.paired_by_model(pe, "delta_gamma_land_coubgc") == {"fire": [2.0]}
    eta = wm.paired_by_model(pe, "eta_beta_land")["fire"]
    assert eta == pytest.approx([0.5 ** (2 / 3) * 2.0 ** (1 / 3), 0.8])
