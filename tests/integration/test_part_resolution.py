import io
import json
from pathlib import Path

import pandas as pd
import pytest
import requests

from conftest import Client, make_freegenes, make_reclone, metadata_record
from services.genbank_service import export_part
from services.part_service import PartService


def test_search_restricts_inventory_but_enriches_matching_aliases(tmp_path):
    local = make_reclone(tmp_path)
    upstream = make_freegenes(tmp_path, records=[metadata_record(), metadata_record("BBF10K_000002", "Unique external enzyme")])
    parts = PartService(local, upstream)
    assert len(parts.search_parts("")) == 1
    assert parts.search_parts("odc1").iloc[0]["Part Key"] == "BBF10K_000001"
    assert parts.search_parts("Unique external").empty
    assert parts.search_parts("BBF10K_000002").empty
    assert parts.search_parts("FreeGenes name").iloc[0]["Part Key"] == "BBF10K_000001"
    assert "BBF10K_000002" not in parts.parts
    assert len(parts.search_parts(collection="Local Collection")) == 1
    assert parts.search_parts(collection="Upstream Collection").empty
    assert parts.collections == ["Local Collection"]


def test_empty_reclone_returns_no_inventory(tmp_path):
    local = make_reclone(tmp_path, rows=[])
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()]))
    assert parts.search_parts("FreeGenes").empty


def test_backend_precedence_preserves_conflicting_original_rows(tmp_path):
    local = make_reclone(tmp_path)
    backend = metadata_record(name="Backend name", source="FreeGenes backend", description="Current")
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()], backend_records=[backend]))
    d = parts.get_part_details("ODC_0001")
    assert d["name"] == "Backend name"
    assert d["metadata"]["description"] == "Current"
    assert d["conflicts"]["gene_name_short"] == ["Backend name", "FreeGenes name"]
    assert len(d["source_records"]) == 3


def test_multiple_provider_locations_and_same_record_filters(tmp_path):
    local = make_reclone(tmp_path, locations=[
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000001", "Well Location": "A1", "Bacterial Resistance": "amp", "Growth Strain": "DH5a"},
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000001", "Well Location": "B2", "Bacterial Resistance": "kan", "Growth Strain": "BL21"}])
    upstream = make_freegenes(tmp_path, records=[metadata_record()], locations=[
        {"part_id": "BBF10K_000001", "plate_name": "FG Plate", "plate_number": "7", "well": "C3", "revision": "feed-1"}])
    parts = PartService(local, upstream)
    d = parts.get_part_details("ODC_0001")
    assert len(d["locations"]) == 3
    assert d["locations"][2]["plate_number"] == "7"
    assert d["locations"][0]["plate_number"] is None
    assert len(parts.search_parts(platemap_filters={"resistance": "amp", "strain": "BL21"})) == 0
    assert len(parts.search_parts(platemap_filters={"well_pattern": "C*", "provider": "FreeGenes"})) == 1


def test_upstream_bytes_preferred_to_different_local_file(tmp_path, gb_bytes):
    local = make_reclone(tmp_path)
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "BBF10K_000001.gb").write_bytes(gb_bytes.replace(b"Annotated test", b"Local alternate"))
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()], paths=["genbank/BBF10K_000001.gb"], client=Client(gb_bytes)))
    details = parts.get_part_details("ODC_0001")
    assert details["sequence_status"] == "direct"
    assert export_part(details)["gb"] == gb_bytes
    assert details["provenance"]["genbank"]["source"] == "FreeGenes GitHub"


def test_invalid_upstream_not_silently_replaced_and_metadata_export_works(tmp_path, gb_bytes):
    local = make_reclone(tmp_path)
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "BBF10K_000001.gb").write_bytes(gb_bytes)
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()], paths=["genbank/BBF10K_000001.gb"], client=Client(b"bad")))
    details = parts.get_part_details("ODC_0001")
    assert details["sequence_status"] == "invalid"
    assert details["genbank"] is None
    assert set(export_part(details)) == {"csv", "txt"}


def test_no_manifest_local_genbank_fallback_is_labeled(tmp_path, gb_bytes):
    local = make_reclone(tmp_path)
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "ODC_0001.gb").write_bytes(gb_bytes)
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()]))
    details = parts.get_part_details("ODC_0001")
    assert details["sequence_status"] == "local fallback"
    assert details["genbank"]["raw"] == gb_bytes
    assert details["provenance"]["genbank"]["upstream_status"] == "absent"


def test_duplicate_alias_conflict_requires_exact_bbf_selection(tmp_path, gb_bytes):
    local = make_reclone(tmp_path, rows=[
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000001", "Name": "One", "Collection": "Kit"},
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000002", "Name": "Two", "Collection": "Kit"}])
    client = Client(gb_bytes)
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()], paths=["genbank/BBF10K_000001.gb"], client=client))
    d = parts.get_part_details("ODC_0001")
    assert d["sequence_status"] == "ambiguous identity"
    assert not client.calls
    assert parts.get_part_details("ODC_0001", "BBF10K_000001")["genbank"]["raw"] == gb_bytes


def test_duplicate_preferred_source_fields_remain_ambiguous(tmp_path):
    parts = PartService(make_reclone(tmp_path), make_freegenes(tmp_path, backend_records=[
        metadata_record(name="A", source="FreeGenes backend"), metadata_record(name="B", source="FreeGenes backend")]))
    d = parts.get_part_details("ODC_0001")
    assert d["metadata"]["gene_name_short"] == ["A", "B"]
    assert "gene_name_short" not in d["field_sources"]


def test_local_genbank_candidates_disagree_without_first_match(tmp_path, gb_bytes):
    local = make_reclone(tmp_path)
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "ODC_0001.gb").write_bytes(gb_bytes)
    (folder / "BBF10K_000001.gb").write_bytes(gb_bytes.replace(b"Annotated test", b"Alternate local"))
    parts = PartService(local, make_freegenes(tmp_path))
    d = parts.get_part_details("ODC_0001")
    assert d["genbank"] is None
    assert any("disagree" in w for w in d["warnings"])
