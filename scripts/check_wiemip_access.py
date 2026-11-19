"""Check that `wiemip_registry` can see the WIEMIP data before a long run (on or off the JupyterHub).

Reports the registry version and data paths, checks the 1pctCO2 forcing inputs that wiemip_extract.py
reads, and probes each registry model's baseline BGC run with the same existence test as the extraction.
With --read, also times one land_carbon_stock read (the slow step when the CSV cache is cold).

    WIEMIP_DATA_ROOT=~/wiemip WIEMIP_CSV_PATH=~/wiemip-csv python scripts/check_wiemip_access.py [--read]

Exits with status 1 if a required input is missing or no model is reachable. Like wiemip_extract.py, it
does not depend on tcre_evidence.
"""
import argparse
import os
import sys
import time
import warnings
from importlib import metadata
from pathlib import Path

import wiemip_registry as wr
from wiemip_registry import const

PROBE_VAR = "cVeg"  # as in wiemip_extract.py


def version() -> str:
    try:
        d = metadata.distribution("wiemip_data_processing")
        url = d.read_text("direct_url.json") or ""
        commit = url.split('"commit_id": "')[1][:7] if '"commit_id"' in url else "?"
        return f"{d.version} (commit {commit})"
    except metadata.PackageNotFoundError:
        return "unknown (not installed as wiemip_data_processing)"


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--read", action="store_true", help="also time one land_carbon_stock read")
    args = a.parse_args()
    ok = True

    print(f"registry       {version()} at {Path(wr.__file__).parent}")
    for var, path in (("WIEMIP_DATA_ROOT", const.DATA_ROOT), ("WIEMIP_CSV_PATH", const.CSV_ROOT)):
        src = "set" if var in os.environ else "default"
        print(f"{var:16s} {path} ({src}; {'exists' if Path(path).is_dir() else 'MISSING'})")
    if not Path(const.DATA_ROOT).is_dir():
        print("data root missing: mount or sync the bucket and set WIEMIP_DATA_ROOT")
        sys.exit(1)
    Path(const.CSV_ROOT).mkdir(parents=True, exist_ok=True)  # the registry fills it on first read

    inp = Path(const.DATA_ROOT) / "1pctCO2" / "input"
    co2 = sorted((inp / "co2").glob("*.txt"))
    for name, found in (("CO2 file", co2[0] if len(co2) == 1 else None),
                        ("FaIR file", inp / "FaIR_temperature.nc" if (inp / "FaIR_temperature.nc").is_file() else None)):
        print(f"{name:16s} {found or 'MISSING'}")
        ok &= found is not None

    reachable = []
    for m in wr.models:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                f = wr.retrieve_one_pct_variable(model=m, forcing="stable", simulation="bgc",
                                                 factorial="baseline", variable=PROBE_VAR)
                hit = f.exists()
            except Exception as e:  # report and carry on; one adapter failing should not hide the others
                hit = f"error: {type(e).__name__}"
        print(f"  {m:24s} {'ok' if hit is True else hit or 'baseline bgc not found'}")
        if hit is True:
            reachable.append(m)
    print(f"{len(reachable)} of {len(wr.models)} registry models reachable")
    ok &= bool(reachable)

    if args.read and reachable:
        t = time.time()
        s = wr.land_carbon_stock("1pctCO2", reachable[0], "stable", "bgc", "baseline")
        print(f"read {reachable[0]} baseline bgc: {len(s)} steps in {time.time() - t:.1f} s")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
