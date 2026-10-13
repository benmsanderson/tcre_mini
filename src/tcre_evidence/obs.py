"""Observational datasets for Workstream B (see data/obs/README.md).

Loaders return pandas DataFrames indexed by integer year:

- :func:`gcb`: Global Carbon Budget 2025 v1.0 global budget, PgC/yr (= GtC/yr).
- :func:`noaa_co2`: NOAA GML global annual mean CO2 (ppm), global growth rate
  (ppm/yr) and Mauna Loa annual mean (ppm).
- :func:`hadcrut5`: HadCRUT5 analysis global annual mean anomaly (K, relative
  to 1961-1990) with 95% limits.

Run ``python -m tcre_evidence.obs`` to re-download the raw files into
``data/obs/raw/`` and rebuild the tidy CSVs in ``data/obs/``.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .paths import DATA

OBS = DATA / "obs"
RAW = OBS / "raw"

GCB_VERSION = "Global_Carbon_Budget_2025_v1.0"
GCB_ICOS_ID = "qSjPBsV1drZnYdH-yCJMmkGn"  # https://doi.org/10.18160/GCP-2025
HADCRUT_VERSION = "HadCRUT.5.2.0.0"
NOAA_BASE = "https://gml.noaa.gov/webdata/ccgg/trends/co2/"
NOAA_FILES = ["co2_annmean_gl.txt", "co2_gr_gl.txt", "co2_annmean_mlo.txt"]
HADCRUT_URL = ("https://www.metoffice.gov.uk/hadobs/hadcrut5/data/{v}/analysis/diagnostics/"
               "{v}.analysis.summary_series.global.annual.csv")

GCB_COLUMNS = {
    "fossil emissions excluding carbonation": "fossil",
    "land-use change emissions": "luc",
    "atmospheric growth": "atm_growth",
    "ocean sink": "ocean_sink",
    "land sink": "land_sink",
    "cement carbonation sink": "cement_carbonation_sink",
    "budget imbalance": "budget_imbalance",
}
# 1 ppm = 2.124 PgC (GCB convention)
PPM_TO_PGC = 2.124


# ------------------------------------------------------------------ download
def download_all(raw: Path = RAW) -> None:
    import requests

    raw.mkdir(parents=True, exist_ok=True)
    r = requests.get(f"https://data.icos-cp.eu/objects/{GCB_ICOS_ID}",
                     cookies={"CpLicenseAcceptedFor": GCB_ICOS_ID}, timeout=120)
    r.raise_for_status()
    if not r.content.startswith(b"PK"):
        raise RuntimeError("GCB download did not return an xlsx file")
    (raw / f"{GCB_VERSION}.xlsx").write_bytes(r.content)
    for f in NOAA_FILES:
        r = requests.get(NOAA_BASE + f, timeout=60)
        r.raise_for_status()
        (raw / f).write_text(r.text)
    url = HADCRUT_URL.format(v=HADCRUT_VERSION)
    r = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=60)
    r.raise_for_status()
    (raw / url.rsplit("/", 1)[1]).write_text(r.text)


# ------------------------------------------------------------------ parsing
def _read_gcb_sheet(path: Path, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet, header=None)
    hdr = df.index[df.iloc[:, 0].astype(str).str.strip() == "Year"][0]
    out = df.iloc[hdr + 1:].copy()
    out.columns = [str(c).strip() for c in df.iloc[hdr]]
    out = out[pd.to_numeric(out["Year"], errors="coerce").notna()]
    out["Year"] = out["Year"].astype(int)
    keep = {c: GCB_COLUMNS[c] for c in out.columns if c in GCB_COLUMNS}
    if len(keep) != len(GCB_COLUMNS):
        raise RuntimeError(f"unexpected GCB columns in {sheet}: {list(out.columns)}")
    out = out.rename(columns=keep).set_index("Year")[list(keep.values())].astype(float)
    out.index.name = "year"
    return out


def _read_noaa(path: Path, cols: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", sep=r"\s+", header=None, names=["year", *cols])
    return df.set_index("year")


def build_tidy(raw: Path = RAW, out: Path = OBS) -> None:
    """Write gcb_global.csv, gcb_historical.csv, noaa_co2.csv, hadcrut5_global.csv."""
    x = raw / f"{GCB_VERSION}.xlsx"
    _read_gcb_sheet(x, "Global Carbon Budget").to_csv(out / "gcb_global.csv", float_format="%.6g")
    _read_gcb_sheet(x, "Historical Budget").to_csv(out / "gcb_historical.csv", float_format="%.6g")
    gl = _read_noaa(raw / "co2_annmean_gl.txt", ["co2_gl", "co2_gl_unc"])
    gr = _read_noaa(raw / "co2_gr_gl.txt", ["growth_gl", "growth_gl_unc"])
    mlo = _read_noaa(raw / "co2_annmean_mlo.txt", ["co2_mlo", "co2_mlo_unc"])
    noaa = pd.concat([gl, gr, mlo], axis=1).sort_index()
    noaa.index.name = "year"
    noaa.to_csv(out / "noaa_co2.csv", float_format="%.6g")
    hc = pd.read_csv(raw / f"{HADCRUT_VERSION}.analysis.summary_series.global.annual.csv")
    hc.columns = ["year", "anomaly", "lower_2p5", "upper_97p5"]
    hc.set_index("year").to_csv(out / "hadcrut5_global.csv", float_format="%.6g")


# ------------------------------------------------------------------ loaders
def gcb(historical: bool = False) -> pd.DataFrame:
    """GCB 2025 global budget, PgC/yr; columns fossil, luc, atm_growth,
    ocean_sink, land_sink, cement_carbonation_sink, budget_imbalance.

    ``historical=True`` returns the 1750-2024 'Historical Budget' sheet;
    otherwise the 1959-2024 'Global Carbon Budget' sheet.
    """
    f = OBS / ("gcb_historical.csv" if historical else "gcb_global.csv")
    return pd.read_csv(f, index_col="year")


def noaa_co2() -> pd.DataFrame:
    """NOAA GML annual CO2: co2_gl (ppm, 1979-), growth_gl (ppm/yr, 1959-),
    co2_mlo (ppm, 1959-), each with a 1-sigma ``*_unc`` column."""
    return pd.read_csv(OBS / "noaa_co2.csv", index_col="year")


def hadcrut5(complete_years_only: bool = True) -> pd.DataFrame:
    """HadCRUT5.2.0.0 analysis global annual anomaly (K, vs 1961-1990).

    The final row is a partial year when the file is retrieved mid-year;
    ``complete_years_only`` drops years after 2025 (retrieval 2026-10-13).
    """
    df = pd.read_csv(OBS / "hadcrut5_global.csv", index_col="year")
    return df.loc[:2025] if complete_years_only else df


if __name__ == "__main__":
    download_all()
    build_tidy()
