"""Selective FreeGenes snapshots, public backend metadata, and lazy GenBank bytes.

Snapshot publication is content-addressed and manifest-last. An interrupted refresh
cannot make readers observe a half-written index. Runtime caches are disposable.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import re
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote, urlparse

import requests

from .data_processing import normalize_id
from .genbank_service import InvalidGenBank, parse_genbank

REPO = "freegenes/freegenes.github.io"
SHEET_ID = "1LZCXzBtgey9xv5OH7YGYgp8UMJ27Eyj1aF9IhAW6M6o"
SHEET_URL = f"https://docs.google.com/spreadsheets/d/{SHEET_ID}/gviz/tq"
BBF_PATTERN = re.compile(r"^BBF10K_\d{6}$")
OMIT_FIELDS = {"", "production_sequence", "insert_sequence", "internal_notes_1", "internal_notes_2", "changed"}
# Public descriptive columns only. No customer, shipping, credentials, or internal sheets.
GENES_QUERY = "select " + ",".join(
    ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P",
     "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AE", "AG",
     "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT", "AW"])


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_bytes(raw)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def json_bytes(value) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def read_snapshot(directory: Path) -> tuple[dict, dict]:
    manifest_path = directory / "manifest.json"
    if not manifest_path.exists():
        return {}, {}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    filename = manifest["index_file"]
    if not re.fullmatch(r"index-[a-f0-9]{64}\.json(?:\.gz)?", filename):
        raise ValueError("Invalid FreeGenes index filename")
    raw = (directory / filename).read_bytes()
    if hashlib.sha256(raw).hexdigest() != manifest["index_sha256"]:
        raise ValueError("FreeGenes index checksum mismatch")
    payload = json.loads(gzip.decompress(raw) if filename.endswith(".gz") else raw)
    if manifest.get("schema_version") != 1 or not isinstance(payload.get("records"), list):
        raise ValueError("Unsupported FreeGenes snapshot schema")
    return manifest, payload


def publish_snapshot(directory: Path, payload: dict, manifest: dict) -> dict:
    # Repeated provenance compresses well; hash the stored bytes.
    raw = gzip.compress(json_bytes(payload), mtime=0)
    sha = hashlib.sha256(raw).hexdigest()
    result = {**manifest, "schema_version": 1, "index_file": f"index-{sha}.json.gz",
              "index_sha256": sha, "index_size_bytes": len(raw)}
    atomic_write(directory / result["index_file"], raw)
    atomic_write(directory / "manifest.json", json_bytes(result))
    return result


def parse_metadata_csv(raw: bytes, path: str, provenance: dict) -> list[dict]:
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if not reader.fieldnames or not {"id", "gene_name_short"}.issubset(reader.fieldnames):
        raise ValueError("FreeGenes metadata CSV schema mismatch")
    records = []
    for row_number, row in enumerate(reader, 2):
        part_id = normalize_id(row.get("id"))
        # The live sheet includes a descriptive row and empty reserved rows.
        if not part_id or not BBF_PATTERN.fullmatch(part_id):
            continue
        fields = {k: v.strip() for k, v in row.items() if k not in OMIT_FIELDS and k and isinstance(v, str) and v.strip()}
        fields["id"] = part_id
        records.append({"part_id": part_id, "fields": fields,
                        "provenance": {**provenance, "source_path": path, "row": row_number,
                                       "sha256": hashlib.sha256(raw).hexdigest()}})
    if not records:
        raise ValueError("FreeGenes CSV has no usable part identifiers")
    return records


class FreeGenesClient:
    """A bounded public-source transport; credentials never reach arbitrary URLs."""
    def __init__(self, token: str | None = None):
        self.token = token or os.getenv("GITHUB_TOKEN")

    def get(self, url: str, *, params=None, headers=None, limit=15_000_000):
        parsed = urlparse(url)
        if parsed.scheme != "https" or parsed.netloc not in {"api.github.com", "raw.githubusercontent.com", "docs.google.com"}:
            raise ValueError("Unconfigured FreeGenes source URL")
        request_headers = {"User-Agent": "open-dna-collections", **(headers or {})}
        if parsed.netloc == "api.github.com" and self.token:
            request_headers["Authorization"] = f"Bearer {self.token}"
        # At most one retry on a transient transport failure; do not retry rate limits.
        for attempt in range(2):
            response = None
            deadline = time.monotonic() + 30
            try:
                response = requests.get(url, params=params, headers=request_headers,
                                        timeout=(3, 7), stream=True, allow_redirects=False)
                if response.status_code == 304:
                    response.close()
                    return b"", dict(response.headers), 304
                response.raise_for_status()
                if response.status_code != 200:
                    raise ValueError("Unexpected redirect or response from FreeGenes source")
                chunks, size = [], 0
                try:
                    for chunk in response.iter_content(64 * 1024):
                        if time.monotonic() > deadline:
                            raise requests.Timeout("FreeGenes response time limit exceeded")
                        size += len(chunk)
                        if size > limit:
                            raise ValueError("FreeGenes source exceeds allowed response size")
                        chunks.append(chunk)
                finally:
                    response.close()
                return b"".join(chunks), dict(response.headers), 200
            except (requests.ConnectionError, requests.Timeout):
                if attempt:
                    raise
            finally:
                if response is not None:
                    response.close()
        raise RuntimeError("Unreachable")

    def head(self, etag: str | None = None):
        raw, headers, status = self.get(f"https://api.github.com/repos/{REPO}/commits/master",
                                       headers={"If-None-Match": etag} if etag else {})
        return (json.loads(raw) if raw else None), headers.get("ETag"), status

    def tree(self, commit: str) -> list[str]:
        raw, _, _ = self.get(f"https://api.github.com/repos/{REPO}/git/trees/{commit}", params={"recursive": 1})
        tree = json.loads(raw)
        if tree.get("truncated"):
            raise ValueError("FreeGenes tree is truncated; preserving the previous snapshot")
        return sorted(n["path"] for n in tree["tree"] if n["type"] == "blob")

    def file(self, commit: str, path: str) -> bytes:
        if not re.fullmatch(r"[a-f0-9]{40}", commit) or not re.fullmatch(r"(?:genbank/BBF10K_\d{6}\.gb|product-csvs/\d+\.csv)", path):
            raise ValueError("Invalid FreeGenes file identity")
        return self.get(f"https://raw.githubusercontent.com/{REPO}/{commit}/{quote(path, safe='/')}", limit=5_000_000)[0]

    def genes(self) -> bytes:
        raw, headers, _ = self.get(SHEET_URL, params={"tqx": "out:csv", "sheet": "Genes", "headers": 1, "tq": GENES_QUERY})
        if "text/csv" not in headers.get("Content-Type", ""):
            raise ValueError("FreeGenes backend did not return CSV")
        return raw


class FreeGenesService:
    def __init__(self, base_path: str | Path, client=None, ttl_seconds: int | None = None):
        self.base_path = Path(base_path)
        self.directory = self.base_path / "data" / "freegenes"
        self.runtime = self.base_path / ".runtime" / "freegenes"
        self.client = client or FreeGenesClient()
        self.ttl_seconds = ttl_seconds if ttl_seconds is not None else int(os.getenv("FREEGENES_TTL_SECONDS", "86400"))
        self._lock = threading.RLock()
        self.load_error = None
        try:
            self.manifest, self.payload = read_snapshot(self.directory)
        except (ValueError, KeyError, OSError) as exc:
            self.manifest, self.payload = {}, {}
            self.load_error = str(exc)

    @property
    def revision(self):
        return self.manifest.get("index_sha256", "missing")

    @property
    def records(self) -> list[dict]:
        return self.payload.get("records", [])

    @property
    def backend_records(self) -> list[dict]:
        return self.payload.get("backend_records", [])

    def status(self) -> dict:
        attempt = {}
        try:
            attempt = json.loads((self.runtime / "last_refresh.json").read_text())
        except (OSError, ValueError):
            pass
        return {"manifest": self.manifest, "load_error": self.load_error,
                "last_refresh": attempt,
                "location_status": self.payload.get("location_status", {"status": "unavailable",
                    "message": "No verified FreeGenes plate/well feed configured"})}

    def ensure_fresh(self, force: bool = False) -> dict:
        if os.getenv("FREEGENES_OFFLINE") == "1":
            return {"status": "offline", "message": "Using packaged/cache data; network access disabled"}
        attempt = self.status()["last_refresh"]
        stamp = attempt.get("checked_at") or self.manifest.get("checked_at")
        if not force and stamp:
            age = (datetime.now(timezone.utc) - datetime.fromisoformat(stamp)).total_seconds()
            # Retry failed refreshes after five minutes, rather than on every widget interaction.
            ttl = min(self.ttl_seconds, 300) if attempt.get("status") == "unavailable" else self.ttl_seconds
            if age < ttl:
                return attempt or {"status": "cached", "checked_at": stamp}
        with self._lock:
            try:
                self.refresh()
                result = {"status": "ok", "checked_at": utc_now(),
                          "message": "FreeGenes metadata refresh completed"}
            except (requests.RequestException, ValueError, KeyError, OSError) as exc:
                result = {"status": "unavailable", "checked_at": utc_now(),
                          "message": f"Refresh unavailable ({type(exc).__name__}); previous snapshot retained"}
            try:
                atomic_write(self.runtime / "last_refresh.json", json_bytes(result))
            except OSError:
                result["message"] += "; runtime status cache could not be written"
            return result

    def refresh(self):
        """Optional backend failure preserves its old data; GitHub batch failure promotes nothing."""
        now = utc_now()
        deadline = time.monotonic() + 120
        head, etag, status = self.client.head(self.manifest.get("head_etag"))
        commit = head["sha"] if head else self.manifest["source_commit"]
        if status == 304 or (commit == self.manifest.get("source_commit") and self.records):
            records = self.records
            genbank_paths = self.payload.get("genbank_paths", [])
            commit_date = self.manifest["commit_date"]
        else:
            paths = self.client.tree(commit)
            csv_paths = [p for p in paths if re.fullmatch(r"product-csvs/\d+\.csv", p)]
            genbank_paths = [p for p in paths if re.fullmatch(r"genbank/BBF10K_\d{6}\.gb", p)]
            if not csv_paths or not genbank_paths:
                raise ValueError("Incomplete FreeGenes source tree")
            commit_date = head["commit"]["committer"]["date"]
            def fetch(path):
                if time.monotonic() > deadline:
                    raise requests.Timeout("FreeGenes metadata refresh time limit exceeded")
                return parse_metadata_csv(self.client.file(commit, path), path,
                    {"source": "FreeGenes GitHub", "repo": REPO, "revision": commit,
                     "retrieved_at": now, "url": f"https://raw.githubusercontent.com/{REPO}/{commit}/{path}"})
            with ThreadPoolExecutor(max_workers=4) as pool:
                batches = list(pool.map(fetch, csv_paths))
            records = [row for batch in batches for row in batch]
        backend_records = self.backend_records
        try:
            if time.monotonic() > deadline:
                raise requests.Timeout("FreeGenes metadata refresh time limit exceeded")
            backend_records = parse_metadata_csv(self.client.genes(), "Genes",
                {"source": "FreeGenes backend", "url": SHEET_URL, "retrieved_at": now,
                 "revision": "content-sha256"})
            backend_status = {"status": "ok", "retrieved_at": now}
        except (requests.RequestException, ValueError, OSError):
            backend_status = {"status": "unavailable", "message": "Backend unavailable; keeping its last good cache"}
        payload = {"records": records, "backend_records": backend_records, "genbank_paths": genbank_paths,
                   "backend_status": backend_status, "locations": self.payload.get("locations", []),
                   "location_status": self.payload.get("location_status", {"status": "unavailable",
                       "message": "Genes has no plate/well columns. Packaging/Collections locations are unspecified; showing Reclone locations."})}
        manifest = {"source_repo": REPO, "source_branch": "master", "source_commit": commit,
                    "commit_date": commit_date, "checked_at": now, "head_etag": etag,
                    "counts": {"github_records": len(records), "backend_records": len(backend_records),
                               "genbank_files": len(genbank_paths)}}
        manifest = publish_snapshot(self.directory, payload, manifest)
        self.manifest, self.payload = manifest, payload

    def _previous_genbank(self, part_id: str) -> dict | None:
        try:
            meta = json.loads((self.runtime / "last_successful" / f"{part_id}.json").read_text())
            revision = meta["revision"]
            if not re.fullmatch(r"[a-f0-9]{40}", revision):
                return None
            raw = (self.runtime / "genbank" / revision / f"{part_id}.gb").read_bytes()
            parse_genbank(raw)
            if hashlib.sha256(raw).hexdigest() != meta["sha256"]:
                return None
            return {"status": "stale cached", "raw": raw,
                    "provenance": {**meta, "cache_status": "stale cached",
                                   "requested_revision": self.manifest.get("source_commit")}}
        except (OSError, ValueError, KeyError):
            return None

    def fetch_genbank(self, part_id: str, force: bool = False) -> dict:
        """Return a precise source status, using cached bytes only after validating them."""
        if not BBF_PATTERN.fullmatch(part_id):
            return {"status": "absent", "message": "No FreeGenes BBF identifier"}
        path = f"genbank/{part_id}.gb"
        links = {r["fields"].get("genbank_file_link", "") for r in self.backend_records if r["part_id"] == part_id}
        links.discard("")
        if links and links != {f"https://freegenes.github.io/{path}"}:
            return {"status": "invalid", "message": "FreeGenes backend GenBank links conflict with the exact indexed identity; no DNA selected automatically"}
        previous = self._previous_genbank(part_id)
        if not self.manifest:
            return {"status": "unavailable", "message": "FreeGenes index unavailable"}
        if path not in self.payload.get("genbank_paths", []):
            if previous:
                return {"status": "withdrawn", "message": "GenBank entry was removed from the current FreeGenes index; previous DNA was not substituted"}
            return {"status": "absent", "message": "No GenBank entry in the indexed FreeGenes revision"}
        commit = self.manifest["source_commit"]
        provenance = {"source": "FreeGenes GitHub", "repo": REPO, "revision": commit,
                      "snapshot_date": self.manifest["commit_date"],
                      "url": f"https://raw.githubusercontent.com/{REPO}/{commit}/{path}"}
        cache = self.runtime / "genbank" / commit / f"{part_id}.gb"
        with self._lock:
            if cache.exists():
                try:
                    raw = cache.read_bytes()
                    parse_genbank(raw)
                    meta = json.loads(cache.with_suffix(".json").read_text())
                    if hashlib.sha256(raw).hexdigest() != meta["sha256"]:
                        raise ValueError("Cached GenBank checksum mismatch")
                    return {"status": "cached", "raw": raw,
                            "provenance": {**provenance, **meta, "cache_status": "cached"}}
                except (OSError, ValueError, KeyError):
                    pass  # Corrupt cache is never used as a fallback.
            negative = cache.with_suffix(".negative.json")
            try:
                stored = json.loads(negative.read_text())
                if not force and time.time() - stored["timestamp"] < 120:
                    if stored["status"] == "unavailable" and previous:
                        return previous
                    return {k: v for k, v in stored.items() if k != "timestamp"}
            except (OSError, ValueError, KeyError):
                pass
            try:
                if os.getenv("FREEGENES_OFFLINE") == "1":
                    raise requests.ConnectionError("offline")
                raw = self.client.file(commit, path)
                parse_genbank(raw)
                meta = {"retrieved_at": utc_now(), "sha256": hashlib.sha256(raw).hexdigest()}
                atomic_write(cache, raw)
                atomic_write(cache.with_suffix(".json"), json_bytes(meta))
                atomic_write(self.runtime / "last_successful" / f"{part_id}.json", json_bytes({**provenance, **meta}))
                return {"status": "direct", "raw": raw,
                        "provenance": {**provenance, **meta, "cache_status": "direct"}}
            except InvalidGenBank:
                result = {"status": "invalid", "message": "FreeGenes returned an invalid GenBank record"}
            except requests.HTTPError as exc:
                result = {"status": "absent" if exc.response is not None and exc.response.status_code == 404 else "unavailable",
                          "message": "FreeGenes GenBank request failed"}
            except (requests.RequestException, OSError, ValueError):
                result = {"status": "unavailable", "message": "FreeGenes GenBank request unavailable"}
            try:
                atomic_write(negative, json_bytes({**result, "timestamp": time.time()}))
            except OSError:
                result["message"] += "; runtime cache could not be written"
            if result["status"] == "unavailable" and previous:
                return previous
            return result
