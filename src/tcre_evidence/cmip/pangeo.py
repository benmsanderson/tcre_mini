"""Pangeo CMIP6 (Google Cloud zarr) access and global-mean reduction.

Reduction conventions
---------------------
- Spatial reduction is done lazily (dask) on the monthly field; only the
  monthly global series is loaded. Annual values use month-length weights
  from the dataset's own calendar (``time.dt.days_in_month``); only complete
  years (12 months) are kept.
- ``tas`` (Amon, K): global mean weighted by ``areacella``; if no cell area is
  found, cos(lat) weights are used and recorded.
- ``nbp`` (Lmon, kg m-2 s-1): CMIP6 ``cell_methods`` is ``area: mean where
  land``, i.e. the value is per unit land area of the cell. The global total is
  sum(nbp * areacella * sftlf/100) * seconds-in-month, summed over the year,
  in PgC/yr. If ``cell_methods`` does not contain ``where land`` the sftlf
  factor is not applied. Positive = into land (CMIP6 standard name
  ``surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes``).
- ``fgco2`` (Omon, kg m-2 s-1): global total sum(fgco2 * areacello [* sftof/100])
  in PgC/yr, positive = into ocean. ``sftof`` is applied when available on
  the same grid (it matters for regridded ``gr`` ocean output and for
  partially wet cells); otherwise areacello alone, recorded.
- ``co2`` (ppm): from ``co2mass`` (Amon, kg CO2) as
  co2mass * 12/44 / 1e12 / 2.124; else the lowest model or pressure level of
  ``co2`` (mol mol-1), area-weighted mean ignoring missing values, * 1e6.

Cell areas are taken from the catalogue (``fx``/``Ofx`` tables, any
experiment and member of the same ``source_id`` and ``grid_label``) and must
match the data's horizontal shape.
"""
from __future__ import annotations

import json
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

CATALOGUE_CSV_URL = "https://storage.googleapis.com/cmip6/pangeo-cmip6.csv"
CACHE = Path(os.environ.get("TCRE_CACHE", Path.home() / ".cache" / "tcre_evidence"))

EXPERIMENTS = [
    "historical", "esm-hist", "piControl", "esm-piControl",
    "1pctCO2", "1pctCO2-bgc", "1pctCO2-rad",
]
TARGET_MODELS = [
    "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CESM2", "CMCC-ESM2",
    "CNRM-ESM2-1", "EC-Earth3-CC", "GFDL-ESM4", "IPSL-CM6A-LR", "MIROC-ES2L",
    "MPI-ESM1-2-LR", "NorESM2-LM", "UKESM1-0-LL", "MRI-ESM2-0", "GISS-E2-1-G",
]
VAR_TABLE = {"tas": "Amon", "nbp": "Lmon", "fgco2": "Omon", "co2": "Amon", "co2mass": "Amon"}
UNITS = {"tas": "K", "nbp": "PgC/yr (positive into land)",
         "fgco2": "PgC/yr (positive into ocean)", "co2": "ppm"}
PPM_TO_PGC = 2.124
SECONDS_PER_DAY = 86400.0
KG_PER_PG = 1e12

# Sibling source_ids on the identical atmosphere grid, used only as a last
# resort for missing fx fields (shape must still match; use is recorded).
FX_SIBLINGS = {
    "EC-Earth3-CC": ["EC-Earth3", "EC-Earth3-Veg"],
}


# ---------------------------------------------------------------- catalogue
_CAT: pd.DataFrame | None = None


def load_catalogue(refresh: bool = False) -> pd.DataFrame:
    """Return the Pangeo CMIP6 catalogue as a DataFrame (cached locally as parquet)."""
    global _CAT
    if _CAT is not None and not refresh:
        return _CAT
    CACHE.mkdir(parents=True, exist_ok=True)
    pq = CACHE / "pangeo-cmip6.parquet"
    if refresh or not pq.exists():
        cols = ["activity_id", "institution_id", "source_id", "experiment_id", "member_id",
                "table_id", "variable_id", "grid_label", "zstore", "version"]
        df = pd.read_csv(CATALOGUE_CSV_URL, usecols=cols)
        tmp = pq.with_suffix(f".tmp{os.getpid()}")
        df.to_parquet(tmp)
        os.replace(tmp, pq)
    _CAT = pd.read_parquet(pq)
    return _CAT


def query(**facets) -> pd.DataFrame:
    """Filter the catalogue on facets; list values are allowed (``isin``)."""
    df = load_catalogue()
    m = np.ones(len(df), bool)
    for k, v in facets.items():
        if v is None:
            continue
        m &= df[k].isin(v).to_numpy() if isinstance(v, (list, tuple, set)) else (df[k] == v).to_numpy()
    return df[m]


def _member_key(member: str):
    nums = [int(x) for x in re.findall(r"\d+", member)]
    return (member != "r1i1p1f1", tuple(nums))


def choose_member(source_id: str, experiment_id: str, variables: list[str]) -> str | None:
    """Pick one member for (model, experiment) maximising variable coverage.

    Ties are broken by preferring r1i1p1f1 and then the lowest (r, i, p, f).
    """
    rows = query(source_id=source_id, experiment_id=experiment_id,
                 variable_id=list(variables))
    rows = rows[[VAR_TABLE.get(v) == t for v, t in zip(rows.variable_id, rows.table_id)]]
    if rows.empty:
        return None
    cov = rows.groupby("member_id").variable_id.nunique()
    return sorted(cov.index, key=lambda m: (-cov[m], _member_key(m)))[0]


def find_store(source_id, experiment_id, variable_id, member_id=None):
    """Return the catalogue row (Series) for a variable, preferring ``member_id``."""
    rows = query(source_id=source_id, experiment_id=experiment_id,
                 variable_id=variable_id, table_id=VAR_TABLE[variable_id])
    if rows.empty:
        return None
    if member_id is not None and (rows.member_id == member_id).any():
        rows = rows[rows.member_id == member_id]
    else:
        rows = rows.loc[sorted(rows.index, key=lambda i: _member_key(rows.at[i, "member_id"]))]
        rows = rows[rows.member_id == rows.member_id.iloc[0]]
    # prefer native grid, then latest version
    rows = rows.assign(_g=(rows.grid_label != "gn").astype(int)).sort_values(
        ["_g", "version"], ascending=[True, False])
    return rows.iloc[0]


def open_store(zstore: str) -> xr.Dataset:
    return xr.open_zarr(zstore, storage_options={"token": "anon"}, consolidated=True,
                        use_cftime=True) if _old_xr() else xr.open_zarr(
        zstore, storage_options={"token": "anon"}, consolidated=True,
        decode_times=xr.coders.CFDatetimeCoder(use_cftime=True))


def _old_xr():
    return not hasattr(xr, "coders")


# ---------------------------------------------------------------- geometry
def spatial_dims(da: xr.DataArray) -> list[str]:
    skip = {"time", "plev", "lev", "level", "alevel", "olevel", "bnds", "nbnd", "vertices"}
    return [d for d in da.dims if d not in skip]


def _fx_candidates(source_id, variable_id, grid_label):
    table = "Ofx" if variable_id in ("areacello", "sftof") else "fx"
    rows = query(source_id=source_id, variable_id=variable_id, table_id=table)
    if rows.empty:
        return rows
    rows = rows.assign(_g=(rows.grid_label != grid_label).astype(int),
                       _m=[_member_key(m) for m in rows.member_id])
    return rows.sort_values(["_g", "_m", "version"], ascending=[True, True, False])


def _fx_cache_path(source_id, variable_id, shape):
    return CACHE / "fx" / f"{source_id}_{variable_id}_{'x'.join(map(str, shape))}.npy"


def get_fx(source_id: str, variable_id: str, grid_label: str, shape: tuple,
           allow_esgf: bool = True) -> tuple[np.ndarray | None, str]:
    """Return (field, provenance) for an fx variable matching ``shape``.

    Search order: local cache, Pangeo same source_id, ESGF same source_id,
    Pangeo sibling source_id. Returns (None, 'missing') if nothing matches.
    """
    shape = tuple(shape)
    cp = _fx_cache_path(source_id, variable_id, shape)
    if cp.exists():
        meta = cp.with_suffix(".json")
        prov = json.loads(meta.read_text())["provenance"] if meta.exists() else "cache"
        return np.load(cp), prov

    def _store(arr, prov):
        cp.parent.mkdir(parents=True, exist_ok=True)
        tmp = cp.with_suffix(f".tmp{os.getpid()}.npy")
        np.save(tmp, arr)
        os.replace(tmp, cp)
        cp.with_suffix(".json").write_text(json.dumps({"provenance": prov}))
        return arr, prov

    def _try_pangeo(src):
        for _, r in _fx_candidates(src, variable_id, grid_label).iterrows():
            try:
                ds = xr.open_zarr(r.zstore, storage_options={"token": "anon"}, consolidated=True)
                arr = ds[variable_id].squeeze(drop=True).values.astype("f8")
            except Exception:  # noqa: BLE001
                continue
            if arr.shape == shape:
                return arr, f"pangeo:{r.zstore}"
        return None

    got = _try_pangeo(source_id)
    if got:
        return _store(*got)
    if allow_esgf:
        from . import esgf
        got = esgf.get_fx(source_id, variable_id, shape)
        if got:
            return _store(*got)
    for sib in FX_SIBLINGS.get(source_id, []):
        got = _try_pangeo(sib)
        if got:
            return _store(got[0], f"sibling {sib}: {got[1]}")
    return None, "missing"


def _bounds_1d(c: np.ndarray) -> np.ndarray:
    mid = 0.5 * (c[1:] + c[:-1])
    return np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def regular_grid_area(da: xr.DataArray) -> np.ndarray | None:
    """Exact spherical cell area (m2) for a regular lat-lon grid, else None."""
    if "lat" not in da.dims or "lon" not in da.dims:
        return None
    lat = da["lat"].values.astype("f8")
    lon = da["lon"].values.astype("f8")
    R = 6371000.0
    lb = np.clip(_bounds_1d(lat), -90, 90)
    ob = _bounds_1d(lon)
    dlat = np.abs(np.sin(np.deg2rad(lb[1:])) - np.sin(np.deg2rad(lb[:-1])))
    dlon = np.abs(np.deg2rad(ob[1:] - ob[:-1]))
    area = R**2 * np.outer(dlat, dlon)
    order = [da.dims.index("lat"), da.dims.index("lon")]
    return area if order[0] < order[1] else area.T


def coslat_weights(da: xr.DataArray) -> np.ndarray | None:
    if "lat" in da.dims and "lon" in da.dims:
        w = np.cos(np.deg2rad(da["lat"].values))[:, None] * np.ones(da.sizes["lon"])[None, :]
        return w if da.dims.index("lat") < da.dims.index("lon") else w.T
    if "lat" in da.coords and da["lat"].ndim == 2:
        return np.cos(np.deg2rad(da["lat"].values))
    return None


# ---------------------------------------------------------------- reduction
def annual_from_monthly(series: xr.DataArray, how: str) -> pd.Series:
    """Monthly global series -> annual, month-length weighted.

    how='mean': sum(x*d)/sum(d). how='flux_total': sum(x*d*86400) (x in kg/s ->
    kg per year) / 1e12 -> PgC/yr. Only years with 12 monthly values are kept.
    """
    s = series.load()
    t = s["time"]
    years = t.dt.year.values
    days = t.dt.days_in_month.values.astype("f8")
    vals = s.values.astype("f8")
    df = pd.DataFrame({"year": years, "x": vals, "d": days})
    df = df[np.isfinite(df.x)]
    g = df.groupby("year")
    n = g.size()
    if how == "mean":
        out = g.apply(lambda q: (q.x * q.d).sum() / q.d.sum(), include_groups=False)
    elif how == "flux_total":
        out = g.apply(lambda q: (q.x * q.d).sum() * SECONDS_PER_DAY / KG_PER_PG,
                      include_groups=False)
    else:
        raise ValueError(how)
    return out[n == 12]


def _prep_time(ds: xr.Dataset, max_years: int | None) -> xr.Dataset:
    # drop duplicated time stamps occasionally present in the cloud stores
    _, idx = np.unique(ds["time"].values, return_index=True)
    if len(idx) != ds.sizes["time"]:
        ds = ds.isel(time=np.sort(idx))
    if max_years is not None:
        ds = ds.isel(time=slice(0, 12 * max_years))
    return ds


def reduce_variable(ds: xr.Dataset, variable_id: str, source_id: str, grid_label: str,
                    max_years: int | None = None, allow_esgf_fx: bool = True
                    ) -> tuple[pd.Series, dict]:
    """Reduce one CMIP6 variable dataset to an annual global series.

    Returns (series indexed by year, info dict describing the weighting).
    ``variable_id`` may be 'tas', 'nbp', 'fgco2', 'co2mass' or 'co2'; the
    latter two both return ppm.
    """
    ds = _prep_time(ds, max_years)
    da = ds[variable_id]
    info = {"units_in": da.attrs.get("units"), "cell_methods": da.attrs.get("cell_methods"),
            "cell_measures": da.attrs.get("cell_measures")}

    if variable_id == "co2mass":
        da = da.squeeze(drop=True)
        extra = [d for d in da.dims if d != "time"]
        if extra:
            da = da.sum(extra)
        ppm = da * 12.0 / 44.0 / KG_PER_PG / PPM_TO_PGC
        info["method"] = "co2mass kg CO2 * 12/44 / 1e12 / 2.124"
        return annual_from_monthly(ppm, "mean"), info

    if variable_id == "co2":
        vdim = next((d for d in ("plev", "lev", "alevel", "level") if d in da.dims), None)
        if vdim is not None:
            v = da[vdim]
            if vdim == "plev":
                k = int(np.nanargmax(v.values))
            else:
                pos = v.attrs.get("positive", "")
                k = int(np.nanargmin(v.values)) if pos == "up" else int(np.nanargmax(v.values))
                # hybrid sigma 'lev' in [0,1]: surface is max; in Pa: max
            da = da.isel({vdim: k})
            info["level"] = f"{vdim}={float(v.values[k])}"

    sd = spatial_dims(da)
    shape = tuple(da.sizes[d] for d in sd)
    kind = "areacello" if variable_id == "fgco2" else "areacella"
    area, prov = get_fx(source_id, kind, grid_label, shape, allow_esgf=allow_esgf_fx)
    info["area"] = prov
    if area is None:
        reg = regular_grid_area(da)
        if variable_id in ("tas", "co2"):
            w = coslat_weights(da)
            if w is None:
                raise RuntimeError(f"no area weights for {source_id} {variable_id}")
            area = w
            info["area"] = "cos(lat) weights (no cell area found)"
        elif reg is not None:
            area = reg
            info["area"] = "spherical cell area from regular lat/lon (no cell area found)"
        else:
            raise RuntimeError(f"no cell area for {source_id} {variable_id} {shape}")
    weights = area.copy()

    if variable_id == "nbp":
        cm = (info["cell_methods"] or "")
        if "where land" in cm:
            frac, fprov = get_fx(source_id, "sftlf", grid_label, shape, allow_esgf=allow_esgf_fx)
            if frac is None:
                info["land_fraction"] = "MISSING sftlf: areacella only (overestimates coastal cells)"
            else:
                frac = np.nan_to_num(frac, nan=0.0)
                if np.nanmax(frac) <= 1.0 + 1e-6:
                    frac = frac * 100.0
                weights = weights * frac / 100.0
                info["land_fraction"] = f"sftlf/100 ({fprov})"
        else:
            info["land_fraction"] = "not applied (cell_methods lacks 'where land')"
    if variable_id == "fgco2":
        frac, fprov = get_fx(source_id, "sftof", grid_label, shape, allow_esgf=False)
        if frac is not None:
            frac = np.nan_to_num(frac, nan=0.0)
            if np.nanmax(frac) <= 1.0 + 1e-6:
                frac = frac * 100.0
            ratio = float(np.nansum(area * frac / 100.0) / np.nansum(area))
            weights = weights * frac / 100.0
            info["sea_fraction"] = f"sftof/100 ({fprov}); sum(area*sftof)/sum(area)={ratio:.4f}"
        else:
            info["sea_fraction"] = "not applied (sftof not found)"

    w = xr.DataArray(np.nan_to_num(weights, nan=0.0), dims=sd)
    da = da.drop_vars([c for c in da.coords if c not in da.dims and c != "time"], errors="ignore")
    for d in sd:
        if d in da.coords:
            da = da.drop_vars(d)

    if variable_id in ("tas", "co2"):
        num = (da * w).sum(sd, skipna=True)
        den = (w * da.notnull()).sum(sd)
        s = num / den
        out = annual_from_monthly(s, "mean")
        if variable_id == "co2":
            # mol mol-1 -> ppm; some models (e.g. CanESM5) store ppm with units '1e-06'
            raw = float(np.nanmean(out.values))
            if 1e-5 < raw < 1e-2:
                out = out * 1e6
                info["scaling"] = "mol mol-1 * 1e6"
            elif 100 < raw < 5000:
                info["scaling"] = f"values already in ppm (units attr {info['units_in']!r})"
            else:
                raise RuntimeError(f"unrecognised co2 magnitude {raw} ({info['units_in']})")
    else:  # flux totals (kg m-2 s-1 * m2 -> kg/s)
        s = (da * w).sum(sd, skipna=True)
        out = annual_from_monthly(s, "flux_total")
    info["total_weight_m2"] = float(np.nansum(weights)) if variable_id != "tas" or "cos" not in str(info["area"]) else None
    return out, info


def dataset_attrs(ds: xr.Dataset) -> dict:
    keys = ["parent_experiment_id", "parent_variant_label", "branch_time_in_parent",
            "branch_time_in_child", "parent_time_units", "branch_method", "variant_label",
            "grid_label", "tracking_id", "further_info_url"]
    out = {k: ds.attrs.get(k) for k in keys if k in ds.attrs}
    try:
        out["calendar"] = ds["time"].encoding.get("calendar") or str(ds["time"].dt.calendar)
    except Exception:  # noqa: BLE001
        pass
    return {k: (v.item() if hasattr(v, "item") else v) for k, v in out.items()}


def fetch_variable(source_id: str, experiment_id: str, variable_id: str,
                   member_id: str | None = None, max_years: int | None = None
                   ) -> tuple[pd.Series, dict] | None:
    """Fetch and reduce one variable from Pangeo. Returns None if not in the catalogue.

    For ``variable_id='co2'`` tries ``co2mass`` first.
    """
    cands = ["co2mass", "co2"] if variable_id == "co2" else [variable_id]
    for v in cands:
        row = find_store(source_id, experiment_id, v, member_id)
        if row is None:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds = open_store(row.zstore)
            s, info = reduce_variable(ds, v, source_id, row.grid_label, max_years=max_years)
        info.update({"source": "pangeo", "variable_id": v, "member_id": row.member_id,
                     "grid_label": row.grid_label, "zstore": row.zstore,
                     "version": str(row.version), "attrs": dataset_attrs(ds)})
        return s, info
    return None
