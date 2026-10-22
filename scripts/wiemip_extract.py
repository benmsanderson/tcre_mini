"""Extract WIEMIP 1pctCO2 global land carbon and forcing into data/wiemip/ (the only step that reads WIEMIP files).

Run with the JupyterHub Python, which has the maintained `wiemip_registry` installed:

    WIEMIP_DATA_ROOT=/srv/wiemip /opt/tljh/user/bin/python scripts/wiemip_extract.py

For every registry model, factorial and 1pctCO2 run that `exists()`, takes the registry's
`land_carbon_stock` (PgC, sum of the adapter's land_carbon_variables, from the registry CSV cache) and
annual-means it. ctrl and bgc are read with forcing "stable" (registry rule); cou and rad with each GCM
pattern. The forcing is the WIEMIP CO2 file and the FaIR GMST config the drivers were built from (the
climate inputs are named `clim3_50perc`, i.e. FaIR config 1051936, the 50th percentile, for all three
patterns).

Writes data/wiemip/cland_annual.csv (model, factorial, sim, forcing, year, cland_pgc),
data/wiemip/forcing.csv (year, co2_ppm, gmst_k) and data/wiemip/provenance.json.
This script does not depend on tcre_evidence, so it runs without PyMC.
"""
import datetime
import hashlib
import json
import os
import warnings
from pathlib import Path

import pandas as pd
import xarray as xr

import wiemip_registry as wr
from wiemip_registry import const

OUT = Path(__file__).resolve().parents[1] / "data" / "wiemip"
RUNS = [("ctrl", "stable"), ("bgc", "stable")] + [(s, p) for s in ("cou", "rad") for p in ("ukesm", "gfdl", "ipsl")]
FAIR_CONFIG = 1051936  # 50th percentile (attribute of FaIR_temperature.nc); matches 'clim3_50perc' driver names
PROBE_VAR = "cVeg"     # existence probe for a run; every adapter lists cVeg in its land carbon pools


def _tree_hash(root: Path) -> str:
    h = hashlib.sha256()
    for f in sorted(root.rglob("*.py")):
        h.update(f.relative_to(root).as_posix().encode())
        h.update(f.read_bytes())
    return h.hexdigest()[:16]


def extract() -> pd.DataFrame:
    rows = []
    for m in wr.models:
        for fac in wr.factorials:
            for sim, frc in RUNS:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    try:
                        f = wr.retrieve_one_pct_variable(model=m, forcing=frc, simulation=sim, factorial=fac,
                                                         variable=PROBE_VAR)
                        if not f.exists():
                            continue
                    except Exception:
                        continue
                    s = wr.land_carbon_stock("1pctCO2", m, frc, sim, fac)
                a = s.groupby(s.index.year).mean()
                rows += [(m, fac, sim, frc, int(y), float(v)) for y, v in a.items()]
                print(f"{m:12s} {fac:20s} {sim:5s} {frc:7s} {len(a)} yr", flush=True)
    return pd.DataFrame(rows, columns=["model", "factorial", "sim", "forcing", "year", "cland_pgc"])


def forcing() -> pd.DataFrame:
    root = Path(const.DATA_ROOT) / "1pctCO2" / "input"
    (co2_file,) = sorted((root / "co2").glob("*.txt"))
    co2 = pd.read_csv(co2_file, sep=r"\s+", header=None, names=["year", "co2_ppm"])
    with xr.open_dataset(root / "FaIR_temperature.nc") as d:
        t = d.temperature.sel(config=FAIR_CONFIG, scenario="1pctCO2").to_series().rename("gmst_k")
    return co2.merge(t.reset_index()[["year", "gmst_k"]], on="year", how="inner"), co2_file


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    cl = extract()
    cl.to_csv(OUT / "cland_annual.csv", index=False)
    fc, co2_file = forcing()
    fc.to_csv(OUT / "forcing.csv", index=False)
    reg = Path(wr.__file__).parent
    prov = dict(
        created=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        data_root=str(const.DATA_ROOT), csv_cache=str(const.CSV_ROOT),
        registry_path=str(reg), registry_py_sha256_16=_tree_hash(reg),
        co2_file=str(co2_file), fair_config=FAIR_CONFIG, probe_variable=PROBE_VAR,
        n_runs=int(cl.groupby(["model", "factorial", "sim", "forcing"]).ngroups),
        note="Global land carbon from wiemip_registry.land_carbon_stock, annual means. "
             "Unpublished WIEMIP output: share only with WIEMIP authors.",
    )
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=1) + "\n")
    print(f"wrote {OUT} ({prov['n_runs']} runs)")
