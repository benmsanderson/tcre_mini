"""Upstream default priors (external/TCRE_Assessment/config/priors.py), loaded read-only."""
import importlib.util

from .paths import UPSTREAM_REPO

_spec = importlib.util.spec_from_file_location("_upstream_priors", UPSTREAM_REPO / "config" / "priors.py")
priors = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(priors)
