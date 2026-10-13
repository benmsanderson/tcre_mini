"""Repository paths."""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = REPO_ROOT / "external"
UPSTREAM_REPO = EXTERNAL / "TCRE_Assessment"
ESM_DECK = EXTERNAL / "esm-deck"
DATA = REPO_ROOT / "data"
SNAPSHOT = DATA / "upstream_snapshot"
SNAPSHOT_DATA = SNAPSHOT / "DATA"
DERIVED = DATA / "derived"
RESULTS = REPO_ROOT / "results"
