"""Fetch CMIP6 annual global means / totals for Workstream B.

Writes data/derived/cmip_globalmean/<experiment>/<source_id>.csv with columns
year, tas [K], nbp [PgC/yr, +into land], fgco2 [PgC/yr, +into ocean],
co2 [ppm] (where available), and a sidecar <source_id>_meta.json recording
member, grid, store/URLs, version and weighting per variable. Then rebuilds
data/derived/cmip_globalmean/availability.csv from all sidecars.

Source order: Pangeo (gs://cmip6); ESGF only for historical/esm-hist.
piControl/esm-piControl: first 200 years only.

Resumable (existing CSVs are skipped; per-variable results are cached under
$TCRE_CACHE/parts) and parallel-safe (per-output lock files, atomic writes).

Usage: OMP_NUM_THREADS=1 python scripts/fetch_cmip_globalmeans.py [--workers 3]
       [--experiments historical esm-hist] [--models CESM2 ...] [--no-esgf]
       [--availability-only]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from tcre_evidence.cmip import esgf, pangeo  # noqa: E402
from tcre_evidence.paths import DERIVED  # noqa: E402

OUT = DERIVED / "cmip_globalmean"
PARTS = pangeo.CACHE / "parts"
PRIORITY = ["historical", "esm-hist", "1pctCO2", "1pctCO2-bgc", "1pctCO2-rad",
            "piControl", "esm-piControl"]
ESGF_EXPS = {"historical", "esm-hist"}
CONTROL_YEARS = 200
LOCK_STALE_S = 3 * 3600


def variables_for(exp: str) -> list[str]:
    v = ["tas", "nbp", "fgco2"]
    if exp.startswith("esm-"):
        v.append("co2")
    return v


def _lock(path: Path):
    lock = path.with_suffix(".lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, str(os.getpid()).encode())
        os.close(fd)
        return lock
    except FileExistsError:
        if time.time() - lock.stat().st_mtime > LOCK_STALE_S:
            lock.unlink(missing_ok=True)
            return _lock(path)
        return None


def _atomic_write(path: Path, text: str):
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, path)


def fetch_one_var(model, exp, var, member, use_esgf):
    part = PARTS / exp / f"{model}_{var}.json"
    if part.exists():
        d = json.loads(part.read_text())
        return pd.Series(d["values"], index=d["years"], dtype=float), d["info"]
    maxy = CONTROL_YEARS if "piControl" in exp else None
    res, err = None, None
    try:
        res = pangeo.fetch_variable(model, exp, var, member_id=member, max_years=maxy)
    except Exception as e:  # noqa: BLE001
        err = f"pangeo error: {e!r}"
        traceback.print_exc()
    if res is None and use_esgf and exp in ESGF_EXPS:
        try:
            res = esgf.fetch_variable(model, exp, var, member_id=member, max_years=maxy)
        except Exception as e:  # noqa: BLE001
            err = (err or "") + f"; esgf error: {e!r}"
            traceback.print_exc()
    if res is None:
        return None, {"source": "missing", "error": err}
    s, info = res
    info["n_years"] = int(len(s))
    part.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write(part, json.dumps({"years": [int(y) for y in s.index],
                                    "values": [float(x) for x in s.values], "info": info},
                                   default=str))
    return s, info


def sanity_flags(df: pd.DataFrame, exp: str) -> list[str]:
    flags = []
    if "tas" in df and df.tas.notna().any() and not (180 < df.tas.mean() < 320):
        flags.append(f"tas mean {df.tas.mean():.1f} K out of range")
    if "co2" in df and df.co2.notna().any() and not (150 < df.co2.mean() < 2000):
        flags.append(f"co2 mean {df.co2.mean():.1f} ppm out of range")
    if exp in ("historical", "esm-hist", "1pctCO2", "1pctCO2-bgc") and "fgco2" in df:
        f = df.fgco2.dropna()
        if len(f) > 60:
            d = f.iloc[-20:].mean() - f.iloc[:20].mean()
            if d <= 0:
                flags.append(f"fgco2 last-20 minus first-20 yr mean {d:.2f} PgC/yr <= 0 (sign?)")
            if f.iloc[:20].mean() < -0.3:
                flags.append(f"fgco2 initial 20-yr mean {f.iloc[:20].mean():.2f} PgC/yr "
                             "(pre-industrial outgassing or drift; correct with control)")
    return flags


QC_MAX_ABS = {"nbp": 30.0, "fgco2": 30.0}  # PgC/yr, any annual value
QC_MAX_DIFF_STD = {"nbp": 8.0, "fgco2": 1.5}  # PgC/yr; CMIP6 range ~0.6-4.4 and ~0.1-0.3


def qc_series(var: str, s: pd.Series, exp: str = "") -> str | None:
    """Return a failure message if the series is physically implausible."""
    if var == "fgco2" and exp in ("1pctCO2", "1pctCO2-bgc") and len(s) >= 80:
        up = float(np.nanmean(s.values[60:80]))
        if up < 1.0:  # CMIP6 range at year 70 is about 3-6 PgC/yr (Arora et al. 2020)
            return f"fgco2 years 61-80 mean {up:.2f} PgC/yr < 1: no uptake response (data error)"
    lim = QC_MAX_ABS.get(var)
    if lim is not None and np.nanmax(np.abs(s.values)) > lim:
        return f"max |annual| {np.nanmax(np.abs(s.values)):.1f} PgC/yr exceeds {lim} (units/data error)"
    lim = QC_MAX_DIFF_STD.get(var)
    if lim is not None and len(s) > 10 and np.nanstd(np.diff(s.values)) > lim:
        return f"std of year-to-year change {np.nanstd(np.diff(s.values)):.2f} PgC/yr exceeds {lim}"
    return None


def process(task):
    model, exp, use_esgf = task
    out = OUT / exp / f"{model}.csv"
    if out.exists():
        return model, exp, "exists"
    lock = _lock(out)
    if lock is None:
        return model, exp, "locked"
    try:
        t0 = time.time()
        vars_ = variables_for(exp)
        member = pangeo.choose_member(model, exp, vars_ + (["co2mass"] if "co2" in vars_ else []))
        if member is None and not (use_esgf and exp in ESGF_EXPS):
            meta = {"source_id": model, "experiment_id": exp, "variables": {
                v: {"source": "missing"} for v in vars_}}
            _atomic_write(out.with_name(f"{model}_meta.json"), json.dumps(meta, indent=1))
            return model, exp, "no data"
        cols, meta_vars = {}, {}
        for v in vars_:
            s, info = fetch_one_var(model, exp, v, member, use_esgf)
            if s is not None:
                bad = qc_series(v, s, exp)
                if bad:
                    info = {**info, "source_retrieved": info.get("source"), "source": "qc_fail",
                            "qc": bad}
                    s = None
            meta_vars[v] = info
            if s is not None:
                cols[v] = s
            print(f"  {model} {exp} {v}: {info.get('source')} "
                  f"{info.get('member_id', '')} {info.get('n_years', '')}", flush=True)
        meta = {"source_id": model, "experiment_id": exp, "preferred_member": member,
                "units": {v: pangeo.UNITS[v] for v in vars_},
                "annual_method": "month-length-weighted from monthly; complete years only",
                "control_years_kept": CONTROL_YEARS if "piControl" in exp else None,
                "retrieved": time.strftime("%Y-%m-%d"), "variables": meta_vars}
        if not cols:
            _atomic_write(out.with_name(f"{model}_meta.json"), json.dumps(meta, indent=1, default=str))
            return model, exp, "no data"
        df = pd.DataFrame(cols)
        df.index.name = "year"
        df = df.sort_index()
        meta["sanity_flags"] = sanity_flags(df, exp)
        _atomic_write(out.with_name(f"{model}_meta.json"), json.dumps(meta, indent=1, default=str))
        _atomic_write(out, df.to_csv(float_format="%.6g"))
        return model, exp, f"ok {sorted(cols)} {time.time()-t0:.0f}s {meta['sanity_flags']}"
    finally:
        lock.unlink(missing_ok=True)


def build_availability():
    rows = []
    for mf in sorted(OUT.glob("*/*_meta.json")):
        m = json.loads(mf.read_text())
        for v, info in m["variables"].items():
            rows.append({"source_id": m["source_id"], "experiment_id": m["experiment_id"],
                         "variable": v, "source": info.get("source", "missing"),
                         "variable_id_used": info.get("variable_id"),
                         "member_id": info.get("member_id"), "grid_label": info.get("grid_label"),
                         "n_years": info.get("n_years"),
                         "area": info.get("area"), "land_fraction": info.get("land_fraction"),
                         "sea_fraction": info.get("sea_fraction")})
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["source_id", "experiment_id", "variable"])
        _atomic_write(OUT / "availability.csv", df.to_csv(index=False))
    return df


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiments", nargs="+", default=PRIORITY)
    ap.add_argument("--models", nargs="+", default=pangeo.TARGET_MODELS)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--no-esgf", action="store_true")
    ap.add_argument("--availability-only", action="store_true")
    a = ap.parse_args(argv)
    if not a.availability_only:
        pangeo.load_catalogue()
        tasks = [(m, e, not a.no_esgf) for e in a.experiments for m in a.models]
        if a.workers > 1:
            with Pool(a.workers, maxtasksperchild=4) as p:
                for r in p.imap_unordered(process, tasks):
                    print(*r, flush=True)
        else:
            for t in tasks:
                print(*process(t), flush=True)
    df = build_availability()
    print(df.pivot_table(index="source_id", columns=["experiment_id", "variable"],
                         values="source", aggfunc="first").to_string() if not df.empty else "empty")


if __name__ == "__main__":
    sys.exit(main())
