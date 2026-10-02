import io
import json
from pathlib import Path

import pandas as pd
import pytest
import requests

from conftest import Client, make_freegenes, make_reclone, metadata_record
from services.genbank_service import export_part
from services.part_service import PartService, part_display_name


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


def test_polymerases_keep_collection_names_and_search_all_source_descriptions(tmp_path):
    local = make_reclone(tmp_path, rows=[
        {"ODC ID": "ODC_0016", "BBF ID": "BBF10K_003257", "Name": "Taq DNA Polymerase", "Collection": "Enzymes"},
        {"ODC ID": "ODC_0019", "BBF ID": "BBF10K_003261", "Name": "Bst DNA Polymerase, Full Length", "Collection": "Enzymes"},
        {"ODC ID": "ODC_0020", "BBF ID": "BBF10K_003262", "Name": "Bst DNA Polymerase, Large Fragment", "Collection": "Enzymes"}])
    upstream = make_freegenes(tmp_path,
        records=[metadata_record("BBF10K_003257", "THEAQpolA", description="Thermostable enzyme")],
        backend_records=[metadata_record("BBF10K_003257", "THEAQpolA", source="FreeGenes backend", description="Current annotation"),
                         metadata_record("BBF10K_003261", "Bstpol"), metadata_record("BBF10K_003262", "BstpolLF"),
                         metadata_record("BBF10K_009999", "External polymerase")])
    parts = PartService(local, upstream)
    result = parts.search_parts("POLYMERASE").set_index("Part Key")
    assert result["Name"].to_dict() == {r["BBF ID"]: r["Name"] for r in local.main_df.to_dict("records")}
    assert parts.search_parts("THEAQpolA").iloc[0]["Name"] == "Taq DNA Polymerase"
    assert parts.search_parts("thermostable").iloc[0]["Part Key"] == "BBF10K_003257"
    assert parts.parts["BBF10K_003257"]["metadata"]["description"] == "Current annotation"
    assert parts.parts["BBF10K_003257"]["name"] == "THEAQpolA"
    assert parts.search_parts("External polymerase").empty


def test_plate_only_reclone_inventory_is_searchable_without_freegenes(tmp_path):
    local = make_reclone(tmp_path, rows=[], locations=[
        {"BBF ID": "BBF10K_000483", "Name": "Plate-only reporter", "Well Location": "G4"},
        {"BBF ID": "BBF10K_000483", "Name": "Alternate reporter name", "Well Location": "H4"}])
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record("BBF10K_009999", "External only")]))
    result = parts.search_parts("reporter")
    assert result["Part Key"].tolist() == ["BBF10K_000483"]
    assert result.iloc[0]["Locations"] == 2
    assert parts.collections == ["Local Collection"]
    assert parts.search_parts("Alternate reporter").iloc[0]["Part Key"] == "BBF10K_000483"
    assert parts.search_parts("External only").empty
    assert parts.parts["BBF10K_000483"]["source_records"][0]["provenance"]["source_path"].endswith("v1.csv")


def test_plate_conflict_does_not_expand_or_merge_master_identity(tmp_path):
    local = make_reclone(tmp_path, locations=[
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000002", "Name": "Wrong identity", "Well Location": "A1"}])
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record("BBF10K_000002", "External only")]))
    assert set(parts.parts) == {"BBF10K_000001"}
    assert parts.parts["BBF10K_000001"]["locations"][0]["identity_conflict"]
    assert parts.search_parts("Wrong identity").empty


def test_matching_plate_names_are_search_aliases(tmp_path):
    local = make_reclone(tmp_path, locations=[
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000001", "Name": "Well annotation alias", "Well Location": "A1"}])
    parts = PartService(local, make_freegenes(tmp_path))
    assert parts.search_parts("Well annotation alias").iloc[0]["Name"] == "Local name"


def test_legacy_label_preserves_reclone_names_and_handles_missing_fields():
    legacy = {"name": "THEAQpolA", "source_records": [
        {"source": "FreeGenes backend", "fields": {"gene_name_short": "THEAQpolA"}},
        {"source": "Reclone", "fields": {"Name": "Taq DNA Polymerase"}}]}
    assert part_display_name(legacy) == "Taq DNA Polymerase"
    assert part_display_name({"name": "Legacy name"}) == "Legacy name"
    assert part_display_name({}, "BBF10K_003257") == "BBF10K_003257"


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


def test_invalid_upstream_uses_labeled_local_genbank(tmp_path, gb_bytes):
    local = make_reclone(tmp_path)
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "BBF10K_000001.gb").write_bytes(gb_bytes)
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()], paths=["genbank/BBF10K_000001.gb"], client=Client(b"bad")))
    details = parts.get_part_details("ODC_0001")
    assert details["sequence_status"] == "local fallback"
    assert details["genbank"]["raw"] == gb_bytes
    assert details["provenance"]["genbank"]["upstream_status"] == "invalid"
    assert details["provenance"]["genbank"]["source_paths"] == ["genbank/BBF10K_000001.gb"]
    assert export_part(details)["gb"] == gb_bytes
    assert any("local Open DNA collection" in w for w in details["warnings"])


def test_no_manifest_local_genbank_fallback_is_labeled(tmp_path, gb_bytes):
    local = make_reclone(tmp_path)
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "ODC_0001.gb").write_bytes(gb_bytes)
    parts = PartService(local, make_freegenes(tmp_path, records=[metadata_record()]))
    details = parts.get_part_details("ODC_0001")
    assert details["sequence_status"] == "local fallback"
    assert details["genbank"]["raw"] == gb_bytes
    assert details["provenance"]["genbank"]["upstream_status"] == "absent"


@pytest.mark.parametrize("failure,status", [
    ("missing_index", "absent"), ("http404", "absent"),
    ("timeout", "unavailable"), ("withdrawn", "withdrawn")])
def test_freegenes_failures_resolve_local_odc_alias_and_gbk_file(tmp_path, gb_bytes, failure, status):
    local = make_reclone(tmp_path)
    folder = tmp_path / "Local Collection" / "genbank_seq"
    folder.mkdir(parents=True)
    (folder / "ODC_0001.gbk").write_bytes(gb_bytes)
    client = Client(gb_bytes)
    upstream = make_freegenes(tmp_path, records=[metadata_record()], client=client,
                             paths=[] if failure == "missing_index" else ["genbank/BBF10K_000001.gb"])
    if failure == "http404":
        response = requests.Response(); response.status_code = 404
        client.fail = requests.HTTPError(response=response)
    elif failure == "timeout":
        client.fail = requests.Timeout()
    elif failure == "withdrawn":
        assert upstream.fetch_genbank("BBF10K_000001")["status"] == "direct"
        upstream.payload["genbank_paths"] = []
    details = PartService(local, upstream).get_part_details("ODC_0001")
    assert details["sequence_status"] == "local fallback"
    assert details["genbank"]["raw"] == gb_bytes
    assert export_part(details)["gb"] == gb_bytes
    assert details["provenance"]["genbank"]["upstream_status"] == status
    assert details["provenance"]["genbank"]["source_paths"] == ["Local Collection/genbank_seq/ODC_0001.gbk"]


def test_local_file_preferred_to_stale_freegenes_but_stale_cache_remains_last_resort(tmp_path, gb_bytes):
    local = make_reclone(tmp_path)
    client = Client(gb_bytes)
    upstream = make_freegenes(tmp_path, records=[metadata_record()],
                             paths=["genbank/BBF10K_000001.gb"], client=client)
    assert upstream.fetch_genbank("BBF10K_000001")["status"] == "direct"
    upstream.manifest["source_commit"] = "b" * 40
    client.fail = requests.Timeout()
    folder = tmp_path / "genbank"; folder.mkdir()
    path = folder / "ODC_0001.gb"
    local_bytes = gb_bytes.replace(b"Annotated test", b"Local alternate")
    path.write_bytes(local_bytes)
    parts = PartService(local, upstream)
    details = parts.get_part_details("ODC_0001")
    assert details["sequence_status"] == "local fallback"
    assert details["genbank"]["raw"] == local_bytes
    assert details["provenance"]["genbank"]["upstream_status"] == "stale cached"
    path.unlink()
    details = parts.get_part_details("ODC_0001")
    assert details["sequence_status"] == "stale cached"
    assert details["genbank"]["raw"] == gb_bytes


@pytest.mark.parametrize("status", ["absent", "unavailable", "invalid", "withdrawn"])
def test_no_valid_matching_local_file_keeps_failure_and_metadata_exports(tmp_path, gb_bytes, monkeypatch, status):
    local = make_reclone(tmp_path)
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "ODC_0001.gb").write_bytes(b"invalid local file")
    (folder / "ODC_0002.gb").write_bytes(gb_bytes)
    upstream = make_freegenes(tmp_path)
    monkeypatch.setattr(upstream, "fetch_genbank", lambda part_id: {"status": status, "message": "Unavailable upstream file"})
    details = PartService(local, upstream).get_part_details("ODC_0001")
    assert details["sequence_status"] == status
    assert details["genbank"] is None
    assert set(export_part(details)) == {"csv", "txt"}
    assert any("local GenBank candidate" in w for w in details["warnings"])


def test_local_fallback_respects_explicit_bbf_choice_in_conflicting_aliases(tmp_path, gb_bytes):
    local = make_reclone(tmp_path, rows=[
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000001", "Name": "One", "Collection": "Kit"},
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000002", "Name": "Two", "Collection": "Kit"}])
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "BBF10K_000001.gb").write_bytes(gb_bytes)
    (folder / "ODC_0001.gb").write_bytes(gb_bytes)
    parts = PartService(local, make_freegenes(tmp_path))
    assert parts.get_part_details("ODC_0001")["sequence_status"] == "ambiguous identity"
    assert parts.get_part_details("ODC_0001", "BBF10K_000001")["sequence_status"] == "local fallback"
    assert parts.get_part_details("ODC_0001", "BBF10K_000002")["genbank"] is None


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
