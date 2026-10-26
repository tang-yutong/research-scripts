"""Put the repo root and every module directory on sys.path for imports."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
for sub in ("md-analysis", "docking-screen", "multiwfn-batch", "job-radar"):
    sys.path.insert(0, str(ROOT / sub))
