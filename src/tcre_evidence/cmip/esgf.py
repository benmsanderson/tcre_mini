"""ESGF fallback (Solr search API + HTTPS file download).

Used only where Pangeo lacks a variable. Files are downloaded to a temporary
cache directory, reduced with :func:`tcre_evidence.cmip.pangeo.reduce_variable`,
and deleted.
"""
from __future__ import annotations

import os
import re
import shutil
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
import xarray as xr

from . import pangeo

INDEX_NODES = [
    "https://esgf-data.dkrz.de/esg-search/search",
    "https://esgf.ceda.ac.uk/esg-search/search",
]
# data nodes that have been reliable from this environment, tried first
PREFERRED_NODES = ["esgf.ceda.ac.uk", "esgf3.dkrz.de", "esgf-data1.llnl.gov",
                   "eagle.alcf.anl.gov", "esgf-node.ornl.gov"]
MAX_DOWNLOAD_BYTES = 4e9


def search_files(**facets) -> list[dict]:
    """Return ESGF File records (all replicas) for the given facets."""
    params = {"project": "CMIP6", "type": "File", "format": "application/solr+json",
              "latest": "true", "distrib": "true", "limit": 2000,
              "fields": "id,instance_id,dataset_id,data_node,member_id,grid_label,size,url,title,version"}
    params.update({k: v for k, v in facets.items() if v is not None})
    last = None
    for node in INDEX_NODES:
        for attempt in range(2):
            try:
                r = requests.get(node, params=params, timeout=120)
                r.raise_for_status()
                return r.json()["response"]["docs"]
            except Exception as e:  # noqa: BLE001
                last = e
                time.sleep(2)
    raise RuntimeError(f"ESGF search failed: {last}")


def _member_of(doc):
    m = doc.get("member_id")
    return m[0] if isinstance(m, list) else m


def _grid_of(doc):
    g = doc.get("grid_label")
    return g[0] if isinstance(g, list) else g


def _http_url(doc):
    for u in doc.get("url", []):
        url, _, kind = u.split("|")
        if kind == "HTTPServer":
            return url
    return None


def pick_dataset(docs: list[dict], member_id: str | None = None):
    """Group File docs into datasets; choose member (preference as Pangeo).

    Returns (member, grid_label, version_tag, {filename: [urls by replica]}).
    """
    if not docs:
        return None
    by_member = defaultdict(list)
    for d in docs:
        by_member[_member_of(d)].append(d)
    members = sorted(by_member, key=pangeo._member_key)
    mem = member_id if member_id in by_member else members[0]
    ds_docs = by_member[mem]
    # prefer gn grid, latest version tag (from dataset_id)
    def vtag(d):
        m = re.search(r"\.v(\d{8})\|", d["dataset_id"])
        return m.group(1) if m else "0"
    grids = sorted({_grid_of(d) for d in ds_docs}, key=lambda g: g != "gn")
    g = grids[0]
    ds_docs = [d for d in ds_docs if _grid_of(d) == g]
    v = max(vtag(d) for d in ds_docs)
    ds_docs = [d for d in ds_docs if vtag(d) == v]
    files = defaultdict(list)
    for d in ds_docs:
        u = _http_url(d)
        if u:
            files[d["title"]].append((d["data_node"], u, d.get("size", 0)))
    for k in files:
        files[k].sort(key=lambda t: PREFERRED_NODES.index(t[0]) if t[0] in PREFERRED_NODES else 99)
    return mem, g, v, dict(files)


def _download(urls, dest: Path) -> str:
    last = None
    for node, url, _ in urls:
        for u in ([url.replace("http://", "https://"), url] if url.startswith("http://") else [url]):
            try:
                with requests.get(u, stream=True, timeout=(30, 300)) as r:
                    r.raise_for_status()
                    tmp = dest.with_suffix(".part")
                    with open(tmp, "wb") as f:
                        for chunk in r.iter_content(1 << 20):
                            f.write(chunk)
                os.replace(tmp, dest)
                # quick integrity check
                with xr.open_dataset(dest, decode_times=False):
                    pass
                return u
            except Exception as e:  # noqa: BLE001
                last = e
                dest.unlink(missing_ok=True)
    raise RuntimeError(f"download failed for {dest.name}: {last}")


def _workdir(tag: str) -> Path:
    d = pangeo.CACHE / "esgf_tmp" / f"{tag}_{os.getpid()}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch_variable(source_id, experiment_id, variable_id, member_id=None, max_years=None):
    """ESGF analogue of :func:`pangeo.fetch_variable`. Returns None if not found."""
    cands = ["co2mass", "co2"] if variable_id == "co2" else [variable_id]
    for v in cands:
        docs = search_files(source_id=source_id, experiment_id=experiment_id,
                            variable_id=v, table_id=pangeo.VAR_TABLE[v])
        picked = pick_dataset(docs, member_id)
        if not picked:
            continue
        mem, grid, ver, files = picked
        names = sorted(files)
        if max_years is not None:
            names = names[: max(1, int(np.ceil(max_years / 50)) + 1)]  # files are time-sorted
        size = sum(files[n][0][2] or 0 for n in names)
        if size > MAX_DOWNLOAD_BYTES:
            warnings.warn(f"{source_id} {experiment_id} {v}: {size/1e9:.1f} GB exceeds cap; skipped")
            continue
        wd = _workdir(f"{source_id}_{experiment_id}_{v}")
        try:
            used = []
            for n in names:
                used.append(_download(files[n], wd / n))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                kw = ({"use_cftime": True} if pangeo._old_xr() else
                      {"decode_times": xr.coders.CFDatetimeCoder(use_cftime=True)})
                ds = xr.open_mfdataset(sorted(wd.glob("*.nc")), combine="by_coords",
                                       data_vars="minimal", coords="minimal",
                                       compat="override", chunks={"time": 120}, **kw)
                s, info = pangeo.reduce_variable(ds, v, source_id, grid, max_years=max_years)
                attrs = pangeo.dataset_attrs(ds)
                ds.close()
        finally:
            shutil.rmtree(wd, ignore_errors=True)
        info.update({"source": "esgf", "variable_id": v, "member_id": mem, "grid_label": grid,
                     "version": f"v{ver}", "urls": used, "bytes": size, "attrs": attrs})
        return s, info
    return None


def get_fx(source_id: str, variable_id: str, shape: tuple):
    """Return (array, provenance) for an fx field from ESGF matching shape, else None."""
    table = "Ofx" if variable_id in ("areacello", "sftof") else "fx"
    try:
        docs = search_files(source_id=source_id, variable_id=variable_id, table_id=table)
    except Exception:  # noqa: BLE001
        return None
    seen = set()
    docs.sort(key=lambda d: (_grid_of(d) != "gn", pangeo._member_key(_member_of(d))))
    for d in docs:
        if d["title"] in seen:
            continue
        seen.add(d["title"])
        same = [(x["data_node"], _http_url(x), 0) for x in docs if x["title"] == d["title"] and _http_url(x)]
        wd = _workdir(f"fx_{source_id}_{variable_id}")
        try:
            p = wd / d["title"]
            u = _download(same, p)
            with xr.open_dataset(p) as ds:
                arr = ds[variable_id].squeeze(drop=True).values.astype("f8")
        except Exception:  # noqa: BLE001
            continue
        finally:
            shutil.rmtree(wd, ignore_errors=True)
        if arr.shape == tuple(shape):
            return arr, f"esgf:{u}"
        if len(seen) > 6:
            break
    return None
