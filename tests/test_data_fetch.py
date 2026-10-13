"""Sanity checks on fetched CMIP global means and observations (fast; skip if absent)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tcre_evidence.paths import DATA

CMIP = DATA / "derived" / "cmip_globalmean"
OBS = DATA / "obs"


def _csvs(exp):
    d = CMIP / exp
    return sorted(d.glob("*.csv")) if d.exists() else []


def _read(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, index_col="year")


def _control(model, exp="piControl"):
    p = CMIP / exp / f"{model}.csv"
    return _read(p) if p.exists() else None


def test_availability_table():
    p = CMIP / "availability.csv"
    if not p.exists():
        pytest.skip("availability.csv not built")
    df = pd.read_csv(p)
    assert {"source_id", "experiment_id", "variable", "source"} <= set(df.columns)
    assert set(df.source) <= {"pangeo", "esgf", "missing", "qc_fail"}


@pytest.mark.parametrize("p", _csvs("1pctCO2"), ids=lambda p: p.stem)
def test_1pct_warming_year70(p):
    df = _read(p)
    ctl = _control(p.stem)
    if "tas" not in df or ctl is None or "tas" not in ctl:
        pytest.skip("tas or piControl tas missing")
    t70 = df.tas.iloc[60:80].mean()  # 20-yr mean centred on year 70
    dT = t70 - ctl.tas.mean()
    assert 1.0 < dT < 3.5, f"{p.stem}: dT70={dT:.2f}"


@pytest.mark.parametrize("p", _csvs("historical") + _csvs("esm-hist"),
                         ids=lambda p: f"{p.parent.name}/{p.stem}")
def test_hist_cumulative_fgco2(p):
    df = _read(p)
    if "fgco2" not in df or df.fgco2.isna().all():
        pytest.skip("fgco2 missing")
    f = df.fgco2.loc[1850:2014]
    ctl = _control(p.stem, "esm-piControl" if p.parent.name == "esm-hist" else "piControl")
    if ctl is not None and "fgco2" in ctl and ctl.fgco2.notna().sum() > 50:
        # remove the pre-industrial flux (e.g. river-induced outgassing in CNRM, BCC),
        # but only when the first 200 control years are representative of the
        # branch state (MIROC-ES2L's early control is not)
        if abs(ctl.fgco2.mean() - f.iloc[:20].mean()) < 0.3:
            f = f - ctl.fgco2.mean()
    cum = f.sum()
    # esm-hist: models with low simulated CO2 (CNRM-ESM2-1 348 ppm, BCC-CSM2-MR 368 ppm
    # in 2014) take up less, so the lower bound is relaxed
    lo = 50 if p.parent.name == "esm-hist" else 80
    assert lo < cum < 200, f"{p.stem}: cumulative fgco2 {cum:.0f} PgC"


@pytest.mark.parametrize("p", _csvs("historical") + _csvs("esm-hist"),
                         ids=lambda p: f"{p.parent.name}/{p.stem}")
def test_hist_tas_and_co2_ranges(p):
    df = _read(p)
    assert df.index.min() == 1850 and df.index.max() == 2014
    if "tas" in df:
        assert df.tas.between(270, 300).all()
        warm = df.tas.loc[1995:2014].mean() - df.tas.loc[1850:1900].mean()
        assert -0.2 < warm < 2.0
    if "co2" in df and df.co2.notna().any():
        assert 260 < df.co2.loc[1850] < 300
        assert 340 < df.co2.loc[2014] < 460
    for v in ("nbp", "fgco2"):
        if v in df:
            assert df[v].abs().max() < 30


def test_gcb():
    p = OBS / "gcb_global.csv"
    if not p.exists():
        pytest.skip("gcb_global.csv missing")
    from tcre_evidence.obs import gcb
    g = gcb()
    assert 9.0 < g.loc[2020, "fossil"] < 10.5
    assert 2.0 < g.loc[2014:2023, "ocean_sink"].mean() < 3.5
    # budget closes: E_FF + E_LUC = G_ATM + S_OCEAN + S_LAND + S_CEMENT + B_IM
    lhs = g.fossil + g.luc
    rhs = g.atm_growth + g.ocean_sink + g.land_sink + g.cement_carbonation_sink + g.budget_imbalance
    assert np.allclose(lhs, rhs, atol=0.02)


def test_noaa_and_hadcrut():
    if not (OBS / "noaa_co2.csv").exists() or not (OBS / "hadcrut5_global.csv").exists():
        pytest.skip("obs missing")
    from tcre_evidence.obs import hadcrut5, noaa_co2
    n = noaa_co2()
    assert 335 < n.loc[1980, "co2_gl"] < 342
    assert 1.5 < n.loc[2010:2019, "growth_gl"].mean() < 3.0
    h = hadcrut5()
    assert h.index.min() == 1850
    assert 0.8 < h.loc[2015:2024, "anomaly"].mean() - h.loc[1850:1900, "anomaly"].mean() < 1.6
