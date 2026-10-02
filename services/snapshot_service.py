"""Resolve legacy flat caches and immutable manifest-driven generations."""
import hashlib
from pathlib import Path


def artifact_path(base_path: Path, manifest: dict, filename: str) -> Path:
    entry = next((f for f in manifest.get("files", []) if Path(f["path"]).name == filename), None)
    path = (base_path / entry["path"] if entry else base_path / "data" / "cache" / filename).resolve()
    if not path.is_relative_to((base_path / "data" / "cache").resolve()):
        raise ValueError("Cache asset path escapes its cache directory")
    return path


def validate_asset(path: Path, manifest: dict):
    entry = next((f for f in manifest.get("files", []) if Path(f["path"]).name == path.name), None)
    if entry and (not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]):
        raise ValueError(f"Cache checksum mismatch: {path.name}")
