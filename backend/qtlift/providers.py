from __future__ import annotations

import shutil


PROVIDERS = {
    "auto": {"label": "Auto", "scope": "Use local BLAST when available, otherwise exact fallback"},
    "windows": {"label": "Windows BLAST+", "scope": "Native blastn against target A"},
    "wsl": {"label": "WSL BLAST+", "scope": "blastn inside the selected WSL distribution (recommended on Windows)"},
    "exact": {"label": "Exact fallback", "scope": "Perfect-match validation and artificial samples"},
}


def provider_status() -> dict:
    wsl = shutil.which("wsl.exe") or shutil.which("wsl")
    return {
        "auto": {"available": True},
        "windows": {"available": bool(shutil.which("blastn")), "path": shutil.which("blastn")},
        "wsl": {"available": bool(wsl), "path": wsl, "blast_check": "Run `command -v blastn` in the chosen distribution"},
        "exact": {"available": True},
    }


def validate_provider(name: str, options: dict | None = None) -> list[str]:
    options = options or {}; status = provider_status()
    if name not in PROVIDERS: raise ValueError(f"Unknown mapping backend: {name}")
    warnings=[]
    if name == "windows" and not status[name]["available"]: warnings.append("Windows BLAST+ selected but blastn is unavailable.")
    if name == "wsl" and not status[name]["available"]: warnings.append("WSL BLAST+ selected but wsl.exe is unavailable.")
    return warnings
