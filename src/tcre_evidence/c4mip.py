"""Reader for the esm-deck C4MIP CMIP6 global-mean ascii files (1pctCO2 BGC/COU/RAD, piControl).

Files: external/esm-deck/C4MIP_data/c4mip_cmip6_ascii_files_11_models/<model>/<model>_{1pc_BGC,1pc_COU,1pc_RAD,PI_CONTROL}.txt

Harmonisation (HANDOVER §6.0), all explicit in MODELS below:
- column names normalised (quoted headers, HETERES, UKESM 'Year' and lower-case piControl header);
- land and ocean fluxes converted to positive = uptake by a per-model sign;
- temperature in K (UKESM forced runs) converted to degC;
- rows indexed by simulation year 0..N-1 (calendar stamps differ: 1850, 1850.0, 1850.5, 1851.0);
- piControl CO2 is not used (missing for CNRM, a 1 %/yr ramp for JAMSTEC); the CO2 baseline is the
  forced run's first-year value.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .paths import ESM_DECK

C4MIP_DIR = ESM_DECK / "C4MIP_data" / "c4mip_cmip6_ascii_files_11_models"
RUNS = {"bgc": "1pc_BGC", "cou": "1pc_COU", "rad": "1pc_RAD", "ctrl": "PI_CONTROL"}


@dataclass(frozen=True)
class ModelInfo:
    cmip_name: str      # CMIP6 source_id, as used upstream
    uptake_sign: int    # multiply raw fluxes by this to get positive = uptake
    family: str         # shared-code family for leave-one-family-out


# Sign table from HANDOVER §6.0, verified by assert_signs() (cumulative BGC uptake at year 70 > 0).
MODELS = {
    "ACCESS-ESM1.5": ModelInfo("ACCESS-ESM1-5", -1, "UM-ACCESS"),
    "BCC":           ModelInfo("BCC-CSM2-MR", -1, "BCC"),
    "CESM2":         ModelInfo("CESM2", +1, "NCAR"),
    "CNRM-ESM2.1":   ModelInfo("CNRM-ESM2-1", +1, "CNRM-IPSL"),
    "CanESM5":       ModelInfo("CanESM5", +1, "CCCma"),
    "GFDL":          ModelInfo("GFDL-ESM4", +1, "GFDL"),
    "IPSL":          ModelInfo("IPSL-CM6A-LR", +1, "CNRM-IPSL"),
    "JAMSTEC":       ModelInfo("MIROC-ES2L", +1, "MIROC"),
    "MPIESM":        ModelInfo("MPI-ESM1-2-LR", -1, "MPI"),
    "NorESM2LM":     ModelInfo("NorESM2-LM", -1, "NCAR"),
    "UKESM":         ModelInfo("UKESM1-0-LL", -1, "UM-ACCESS"),
}
CMIP_TO_DECK = {v.cmip_name: k for k, v in MODELS.items()}

_RENAME = {
    "YEAR": "year", "ST": "tas", "CO2": "co2", "LA_CO2_FLUX": "fland", "OA_CO2_FLUX": "focean",
    "GPP": "gpp", "NPP": "npp", "HETRES": "rh", "HETERES": "rh", "CLIVE": "cveg", "CDEAD": "csoil",
    "TAS_(K)": "tas", "CO2_(PPM)": "co2", "LANDFLX_(PGC/YR)": "fland", "OCNFLX_(PGC/YR)": "focean",
    "GPP_(PGC/YR)": "gpp", "NPP_(PGC/YR)": "npp", "RH_(PGC/YR)": "rh", "CVEG_(PGC)": "cveg",
    "CSOIL_(PGC)": "csoil",
}


def read_run(model: str, run: str) -> pd.DataFrame:
    """One run as a DataFrame indexed by simulation year 0..N-1; fluxes positive = uptake, tas in degC."""
    f = C4MIP_DIR / model / f"{model}_{RUNS[run]}.txt"
    df = pd.read_csv(f, sep=r"\s+", engine="python")
    df.columns = [_RENAME[c.strip('"').upper()] for c in df.columns]
    df = df.dropna(subset=["tas", "fland", "focean"]).reset_index(drop=True)
    if df.tas.mean() > 200:  # K -> degC (UKESM forced runs)
        df["tas"] = df.tas - 273.15
    s = MODELS[model].uptake_sign
    df["fland"] *= s
    df["focean"] *= s
    df.index.name = "t"
    return df


def read_model(model: str) -> dict[str, pd.DataFrame]:
    """All C4MIP runs (COU, BGC, RAD, piControl) for one esm-deck model, keyed by run name."""
    return {r: read_run(model, r) for r in RUNS}


def assert_signs(runs: dict[str, pd.DataFrame], model: str = "") -> None:
    """Uptake-positive check: cumulative BGC land and ocean uptake at year 70 must be positive."""
    for x in ("fland", "focean"):
        c = runs["bgc"][x].iloc[:70].sum()
        assert c > 0, f"{model}: cumulative BGC {x} at year 70 = {c:.1f} PgC, expected uptake > 0"
    assert abs(runs["ctrl"].tas.mean() - runs["cou"].tas.iloc[:10].mean()) < 1.0, f"{model}: control/forced tas units"
