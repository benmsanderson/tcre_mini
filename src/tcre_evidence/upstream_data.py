"""Loaders for the upstream TCRE_Assessment inputs, read from data/upstream_snapshot.

These mirror the reads in upstream notebooks/TCRE_total.ipynb and
utils/data_utils.py at commit 98c74cb, but read the frozen snapshot rather than
the live Google Sheet. Upstream's convention for β/γ is undocumented (assumed
COU−BGC, to be verified in A1), so variables here keep upstream names.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from .paths import SNAPSHOT, SNAPSHOT_DATA

K = 2.12  # PgC ppm^-1, as upstream config.priors.K

# Upstream notebook cell 13: per-model alpha (K ppm^-1), "ordered to match common models".
# HARDCODED_ALPHA_ASSUMED_ORDER is the order the notebook implicitly assigns
# (np.intersect1d of land and ocean model names, a case-sensitive sort).
HARDCODED_ALPHA = np.array([0.0074, 0.0059, 0.0094, 0.0078, 0.0064,
                            0.0087, 0.0052, 0.0063, 0.0054, 0.0050, 0.0088])


def sheet(name: str) -> pd.DataFrame:
    """One tab of the snapshot Google Sheet: 'beta_gamma', 'gamma_iav' or 'ocean_predictors'."""
    return pd.read_csv(SNAPSHOT / f"gsheet_{name}.csv")


def lookup_table() -> xr.DataArray:
    """Process flags (process: nitrogen, permafrost, fire, veg) x 13 models."""
    return xr.open_dataarray(SNAPSHOT_DATA / "lookup_table.nc").load()


def cmip_sheet_dataset() -> xr.Dataset:
    """Notebook cells 4-6: per-model sheet values plus lookup table, model order as in sheet."""
    df = sheet("beta_gamma")
    co = dict(model=df.Model.values)
    ds = xr.Dataset({v: xr.DataArray(df[v].values.astype(float), coords=co)
                     for v in ["alpha", "beta_land", "gamma_land", "beta_land_tropics",
                               "gamma_land_tropics", "beta_ocean", "gamma_ocean"]})
    ds["gamma_land_extratropics"] = ds.gamma_land - ds.gamma_land_tropics
    ds["lookup_table"] = lookup_table().reindex(model=ds.model)
    return ds


def iav() -> tuple[pd.DataFrame, float, float]:
    """Notebook cells 8-9: IAV table (with OBS row, as upstream passes it) and observed mean, sigma."""
    df = sheet("gamma_iav")
    obs = df.loc[df.Model == "OBS"].iloc[0]
    return df, float(obs.gamma_IAV), float(obs.sigma_iav)


def terhaar_ocean(time_index: int = 69) -> xr.Dataset:
    """Ocean beta, gamma from beta_gamma_terhaar.nc at a time index (69 = year 70, 2xCO2)."""
    ds = xr.load_dataset(SNAPSHOT_DATA / "beta_gamma_terhaar.nc", decode_times=False)
    return xr.Dataset({"beta_ocean": ds.beta.isel(time=time_index),
                       "gamma_ocean": ds.gamma.isel(time=time_index)})


def terhaar_predictors() -> xr.Dataset:
    """AMOC, SSS, CUC per model from Terhaar.txt (as upstream load_ocean_ec_predictors)."""
    t = pd.read_csv(SNAPSHOT_DATA / "Terhaar.txt", header=None, sep=" ")
    return xr.Dataset({"AMOC": ("model", t[1].values), "SSS": ("model", t[2].values),
                       "CUC": ("model", t[3].values)}, coords={"model": t[0].values})


OCEAN_EC_OBS = {"AMOC": (17.02, 0.35), "SSS": (34.04, 0.04), "CUC": (0.546, 0.003)}


def standardize_ocean_ec(ec: xr.Dataset, obs: dict = OCEAN_EC_OBS) -> dict:
    """As upstream data_utils.standardize_ocean_ec (population std via xarray .std())."""
    vs = ["AMOC", "SSS", "CUC"]
    mean = {v: float(ec[v].mean()) for v in vs}
    std = {v: float(ec[v].std()) for v in vs}
    return dict(
        O_m=np.column_stack([(ec[v].values - mean[v]) / std[v] for v in vs]),
        O_obs=np.array([(obs[v][0] - mean[v]) / std[v] for v in vs]),
        O_obs_unc=np.array([obs[v][1] / std[v] for v in vs]),
        mean=mean, std=std,
    )


def land_csv(kind: str = "2xCO2", sources=("CMIP6", "CMIP6+")) -> xr.Dataset:
    """Land beta_L, gamma_L from TCREsource_betagamma.csv (as upstream load_carbon_parameters)."""
    df = pd.read_csv(SNAPSHOT_DATA / "TCREsource_betagamma.csv")
    df = df[df.source.isin(sources)]
    return xr.Dataset({"beta_land": ("model", df[f"beta_L_{kind}"].values),
                       "gamma_land": ("model", df[f"gamma_L_{kind}"].values)},
                      coords={"model": df.model.values})


def carbon_parameters(kind: str = "2xCO2") -> xr.Dataset:
    """Re-implementation of upstream data_utils.load_carbon_parameters.

    Land from TCREsource_betagamma.csv, ocean from beta_gamma_terhaar.nc (time index
    69 for 2xCO2, last for 4xCO2) labelled by Terhaar.txt row order, restricted to
    np.intersect1d of the two model lists (the netCDF model coordinate equals
    Terhaar.txt row order, which upstream uses as labels) (sorted, case-sensitive).
    """
    land = land_csv(kind)
    ocean = terhaar_ocean(69 if kind == "2xCO2" else -1)
    common = np.intersect1d(land.model.values, ocean.model.values)
    return xr.merge([land.sel(model=common), ocean.sel(model=common)])
