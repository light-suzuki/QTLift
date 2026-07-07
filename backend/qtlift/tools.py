from __future__ import annotations

import shutil
import subprocess
import time
from pathlib import Path


_CACHE: tuple[float, dict] | None = None


def detect_tools(config: dict | None = None) -> dict:
    global _CACHE
    config = config or {}
    if not config and _CACHE and time.monotonic() - _CACHE[0] < 60:
        return _CACHE[1]
    result = {}
    for name in ("blastn", "makeblastdb", "minimap2", "paftools.js"):
        configured = config.get(name)
        path = str(Path(configured).resolve()) if configured and Path(configured).exists() else shutil.which(name)
        runtime = "windows" if path else None
        if not path and shutil.which("wsl.exe"):
            probe = subprocess.run(["wsl.exe", "bash", "-lc", f"source ~/.profile >/dev/null 2>&1; command -v '{name}' || true"], text=True, capture_output=True, timeout=10).stdout.strip()
            if probe: path, runtime = probe, "wsl"
        result[name] = {"available": bool(path), "path": path, "runtime": runtime, "configured": bool(configured)}
    if not config:
        _CACHE = (time.monotonic(), result)
    return result
