# hermpymod/paths.py
from pathlib import Path
import os

def find_repo_root(marker: str = ".git") -> Path:
    """Walk up from this file until we find the repo root."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / marker).exists():
            return parent
    raise RuntimeError(f"Could not locate repo root (no {marker} found)")

REPO_ROOT = find_repo_root()
DATA_DIR = str(REPO_ROOT / "code" / ".cache/")
os.makedirs(DATA_DIR, exist_ok=True)

