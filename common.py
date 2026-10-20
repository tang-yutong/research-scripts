"""Shared runtime helpers: exit codes, logging, CLI plumbing, run manifests.

Exit-code contract for every entry point:
  0 success | 2 bad arguments | 3 missing/malformed input | 4 external tool failed
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from dataclasses import fields
from pathlib import Path

EXIT_OK = 0
EXIT_ARGS = 2
EXIT_INPUT = 3
EXIT_TOOL = 4

ROOT = Path(__file__).resolve().parent


def get_logger(name: str, level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s",
                                         "%H:%M:%S"))
        logger.addHandler(h)
    logger.setLevel(getattr(logging, level, logging.INFO))
    return logger


def add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", type=Path, default=None, metavar="JSON",
                   help="override config fields from a JSON file")
    p.add_argument("--dry-run", action="store_true",
                   help="print the planned commands/outputs and exit")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])


def apply_overrides(cfg, path: Path) -> None:
    """Apply a {field_name: value} JSON mapping onto the Config dataclass."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    known = {f.name for f in fields(cfg)}
    for k, v in data.items():
        if k not in known:
            fail(get_logger("config"), f"unknown config field {k!r}", EXIT_ARGS)
        cur = getattr(cfg, k)
        if isinstance(cur, Path):
            v = Path(v)
        elif isinstance(cur, tuple) and isinstance(v, list):
            v = tuple(v)
        setattr(cfg, k, v)
    resolve(cfg)


def resolve(cfg) -> None:
    """Make every Path field absolute (relative ones against the repo root)."""
    for f in fields(cfg):
        v = getattr(cfg, f.name)
        if isinstance(v, Path) and not v.is_absolute():
            setattr(cfg, f.name, ROOT / v)


def fail(logger: logging.Logger, msg: str, code: int = EXIT_INPUT):
    logger.error(msg)
    sys.exit(code)


def require_inputs(logger: logging.Logger, *paths) -> None:
    missing = [str(p) for p in paths if not Path(p).exists()]
    if missing:
        fail(logger, "missing input(s): " + ", ".join(missing), EXIT_INPUT)


def run_tool(logger: logging.Logger, cmd, timeout: int | None = None, **kw):
    """Run an external tool; on non-zero exit, stop with EXIT_TOOL."""
    try:
        p = subprocess.run(cmd, timeout=timeout, **kw)
    except FileNotFoundError:
        fail(logger, f"tool not found on PATH: {cmd[0]}", EXIT_TOOL)
    except subprocess.TimeoutExpired:
        fail(logger, f"tool timed out after {timeout}s: {cmd[0]}", EXIT_TOOL)
    if p.returncode != 0:
        fail(logger, f"tool failed (rc={p.returncode}): {' '.join(map(str, cmd))}",
             EXIT_TOOL)
    return p


def tool_version(cmd) -> str:
    try:
        out = subprocess.run(list(cmd) + ["--version"], capture_output=True,
                             text=True, timeout=20)
        return ((out.stdout or out.stderr).strip().splitlines() or ["unknown"])[0][:80]
    except Exception:
        return "unknown"


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=10, cwd=ROOT)
        return out.stdout.strip() or "n/a"
    except Exception:
        return "n/a"


def write_manifest(out_dir, script: str, argv, cfg, inputs=()) -> Path:
    """Write manifest.json (timestamp, argv, config snapshot, input sizes,
    tool/git versions) next to the outputs of a run -- run provenance."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    snap = {f.name: (str(getattr(cfg, f.name)) if isinstance(getattr(cfg, f.name), Path)
                     else getattr(cfg, f.name))
            for f in fields(cfg)}
    manifest = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "script": script,
        "argv": [str(a) for a in argv],
        "config": snap,
        "inputs": {str(p): Path(p).stat().st_size for p in inputs if Path(p).exists()},
        "tools": {"gmx": tool_version("gmx"), "vina": tool_version("vina"),
                  "obabel": tool_version("obabel")},
        "git_commit": git_commit(),
    }
    path = out_dir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=1, ensure_ascii=False, default=str),
                    encoding="utf-8")
    return path
