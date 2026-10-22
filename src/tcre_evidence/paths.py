"""Repository paths."""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = REPO_ROOT / "external"
UPSTREAM_REPO = EXTERNAL / "TCRE_Assessment"
DATA = REPO_ROOT / "data"
SNAPSHOT = DATA / "upstream_snapshot"
SNAPSHOT_DATA = SNAPSHOT / "DATA"
DERIVED = DATA / "derived"
C4MIP_DATA = DATA / "c4mip" / "c4mip_cmip6_ascii_files_11_models"  # from esm-deck eef04dc, see data/c4mip/README.md
FLAT10_DATA = DATA / "flat10"
RESULTS = REPO_ROOT / "results"
