import json
from pathlib import Path

import pytest
import requests

from conftest import COMMIT, Client, make_freegenes, metadata_record
from services.freegenes_service import (FreeGenesClient, FreeGenesService, parse_metadata_csv,
                                       publish_snapshot, read_snapshot)

CSV = b"id,gene_name_short,product,internal_notes_1,production_sequence\nBBF10K_000001,One,Kit,private,ATGC\n"


def test_metadata_filters_internal_fields_and_keeps_duplicates():
    raw = CSV + b"BBF10K_000001,Duplicate,Kit,secret,TTAA\n"
    records = parse_metadata_csv(raw, "Genes", {"source": "test"})
    assert len(records) == 2
    assert "internal_notes_1" not in records[0]["fields"]
    assert "production_sequence" not in records[0]["fields"]
    assert records[0]["provenance"]["row"] == 2


def test_refresh_pins_content_and_reads_backend(tmp_path):
    client = Client()
    service = FreeGenesService(tmp_path, client=client)
    service.refresh()
    assert client.calls == [(COMMIT, "product-csvs/123.csv")]
    manifest, payload = read_snapshot(service.directory)
    assert manifest["source_commit"] == COMMIT
    assert payload["backend_records"][0]["provenance"]["source"] == "FreeGenes backend"
    assert payload["location_status"]["status"] == "unavailable"


def test_partial_batch_failure_preserves_manifest(tmp_path):
    client = Client()
    service = make_freegenes(tmp_path, records=[metadata_record()], client=client)
    before = (service.directory / "manifest.json").read_bytes()
    client.head = lambda etag=None: ({"sha": "b" * 40, "commit": {"committer": {"date": "2026-10-01T00:00:00Z"}}}, None, 200)
    client.fail = requests.Timeout()
    result = service.ensure_fresh(force=True)
    assert result["status"] == "unavailable"
    assert (service.directory / "manifest.json").read_bytes() == before
    assert service.status()["last_refresh"]["status"] == "unavailable"


def test_backend_failure_retains_old_backend_rows(tmp_path):
    client = Client()
    old = metadata_record(name="Previous backend", source="FreeGenes backend")
    service = make_freegenes(tmp_path, records=[metadata_record()], backend_records=[old], client=client)
    client.backend_fail = True
    service.refresh()
    assert service.backend_records == [old]
    assert service.payload["backend_status"]["status"] == "unavailable"


def test_snapshot_hash_corruption_detected(tmp_path):
    service = make_freegenes(tmp_path, records=[metadata_record()])
    (service.directory / service.manifest["index_file"]).write_text("{}")
    loaded = FreeGenesService(tmp_path)
    assert loaded.records == []
    assert "checksum" in loaded.load_error


def test_lazy_byte_cache_and_retrieval_provenance(tmp_path, gb_bytes):
    client = Client(gb_bytes)
    service = make_freegenes(tmp_path, paths=["genbank/BBF10K_000001.gb"], client=client)
    first = service.fetch_genbank("BBF10K_000001")
    assert first["status"] == "direct"
    assert first["raw"] == gb_bytes
    client.fail = requests.Timeout()
    second = service.fetch_genbank("BBF10K_000001")
    assert second["status"] == "cached"
    assert first["provenance"]["retrieved_at"] == second["provenance"]["retrieved_at"]
    assert len(client.calls) == 1


def test_previous_revision_cache_is_labeled_and_withdrawal_not_hidden(tmp_path, gb_bytes):
    client = Client(gb_bytes)
    service = make_freegenes(tmp_path, paths=["genbank/BBF10K_000001.gb"], client=client)
    service.fetch_genbank("BBF10K_000001")
    service.manifest["source_commit"] = "b" * 40
    client.fail = requests.Timeout()
    result = service.fetch_genbank("BBF10K_000001")
    assert result["status"] == "stale cached"
    assert result["provenance"]["revision"] == COMMIT
    assert result["provenance"]["requested_revision"] == "b" * 40
    service.payload["genbank_paths"] = []
    assert service.fetch_genbank("BBF10K_000001")["status"] == "withdrawn"


def test_invalid_genbank_and_short_negative_cache(tmp_path):
    client = Client(b"<html>not DNA</html>")
    service = make_freegenes(tmp_path, paths=["genbank/BBF10K_000001.gb"], client=client)
    assert service.fetch_genbank("BBF10K_000001")["status"] == "invalid"
    assert service.fetch_genbank("BBF10K_000001")["status"] == "invalid"
    assert len(client.calls) == 1


def test_transient_timeout_uses_negative_cache(tmp_path):
    client = Client()
    client.fail = requests.Timeout()
    service = make_freegenes(tmp_path, paths=["genbank/BBF10K_000001.gb"], client=client)
    assert service.fetch_genbank("BBF10K_000001")["status"] == "unavailable"
    service.fetch_genbank("BBF10K_000001")
    assert len(client.calls) == 1


def test_bad_source_url_is_not_requested(monkeypatch):
    monkeypatch.setattr(requests, "get", lambda *a, **k: pytest.fail("Network called"))
    with pytest.raises(ValueError):
        FreeGenesClient().get("https://example.invalid/credentials")


def test_http_304_keeps_revision_and_uses_conditional_request(tmp_path):
    client = Client()
    service = make_freegenes(tmp_path, records=[metadata_record()], paths=["genbank/BBF10K_000001.gb"], client=client)
    client.head = lambda etag=None: (None, "new-etag", 304)
    service.refresh()
    assert service.manifest["source_commit"] == COMMIT
    assert service.manifest["head_etag"] == "new-etag"
    assert not client.calls


def test_truncated_tree_rejected(monkeypatch):
    client = FreeGenesClient()
    monkeypatch.setattr(client, "get", lambda *a, **k: (b'{"truncated":true,"tree":[]}', {}, 200))
    with pytest.raises(ValueError, match="truncated"):
        client.tree(COMMIT)


def test_manifest_published_last(tmp_path, monkeypatch):
    service = make_freegenes(tmp_path, records=[metadata_record()])
    before = (service.directory / "manifest.json").read_bytes()
    from services import freegenes_service as module
    original = module.atomic_write
    def fail_manifest(path, raw):
        if path.name == "manifest.json":
            raise OSError("simulated interrupted promotion")
        original(path, raw)
    monkeypatch.setattr(module, "atomic_write", fail_manifest)
    with pytest.raises(OSError):
        publish_snapshot(service.directory, {"records": [metadata_record(name="New")]}, service.manifest)
    assert (service.directory / "manifest.json").read_bytes() == before
    assert read_snapshot(service.directory)[1]["records"][0]["fields"]["gene_name_short"] == "FreeGenes name"


def test_compressed_snapshot_roundtrip_and_legacy_compatibility(tmp_path):
    import hashlib
    from services.freegenes_service import json_bytes
    payload = {"records": [metadata_record()]}
    manifest = publish_snapshot(tmp_path, payload, {})
    assert manifest["index_file"].endswith(".json.gz")
    assert read_snapshot(tmp_path)[1] == payload
    raw = json_bytes(payload)
    sha = hashlib.sha256(raw).hexdigest()
    (tmp_path / f"index-{sha}.json").write_bytes(raw)
    (tmp_path / "manifest.json").write_bytes(json_bytes({"schema_version": 1,
        "index_file": f"index-{sha}.json", "index_sha256": sha}))
    assert read_snapshot(tmp_path)[1] == payload


@pytest.mark.parametrize("status", [403, 429])
def test_rate_limit_keeps_previous_snapshot(tmp_path, status):
    client = Client()
    service = make_freegenes(tmp_path, records=[metadata_record()], client=client)
    before = (service.directory / "manifest.json").read_bytes()
    response = requests.Response(); response.status_code = status
    def blocked(etag=None):
        raise requests.HTTPError(response=response)
    client.head = blocked
    assert service.ensure_fresh(force=True)["status"] == "unavailable"
    assert (service.directory / "manifest.json").read_bytes() == before


def test_conflicting_backend_genbank_link_blocks_sequence(tmp_path):
    row = metadata_record(source="FreeGenes backend")
    row["fields"]["genbank_file_link"] = "https://freegenes.github.io/genbank/BBF10K_000002.gb"
    client = Client()
    service = make_freegenes(tmp_path, backend_records=[row], paths=["genbank/BBF10K_000001.gb"], client=client)
    assert service.fetch_genbank("BBF10K_000001")["status"] == "invalid"
    assert not client.calls


def test_unwritable_runtime_cache_reports_failure_without_crashing(tmp_path, monkeypatch):
    from services import freegenes_service as module
    client = Client(); client.fail = requests.Timeout()
    service = make_freegenes(tmp_path, paths=["genbank/BBF10K_000001.gb"], client=client)
    def unwritable(*args):
        raise OSError("read-only runtime")
    monkeypatch.setattr(module, "atomic_write", unwritable)
    assert service.fetch_genbank("BBF10K_000001")["status"] == "unavailable"
    client.head = lambda etag=None: (_ for _ in ()).throw(requests.Timeout())
    assert service.ensure_fresh(force=True)["status"] == "unavailable"
