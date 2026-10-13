"""CMIP6 data access and global-mean reduction (Pangeo first, ESGF fallback).

Conventions used throughout (see ``pangeo`` module docstring for details):

- ``tas``: area-weighted global annual mean, K.
- ``nbp``: global annual total land-atmosphere flux, PgC/yr, positive = into land.
- ``fgco2``: global annual total air-sea flux, PgC/yr, positive = into ocean.
- ``co2``: global annual mean atmospheric CO2, ppm (from ``co2mass`` where
  available, else the lowest level of the 3D ``co2`` field).
"""
from .pangeo import (  # noqa: F401
    EXPERIMENTS,
    TARGET_MODELS,
    annual_from_monthly,
    load_catalogue,
    query,
    reduce_variable,
)
