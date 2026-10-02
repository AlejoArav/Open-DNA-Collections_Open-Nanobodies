from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from conftest import make_freegenes, make_reclone, metadata_record

APP = Path(__file__).resolve().parents[2] / "streamlit_app.py"


def button(app, label):
    return next(b for b in app.button if b.label == label)


@pytest.fixture
def app(tmp_path, monkeypatch, gb_bytes):
    make_reclone(tmp_path, rows=[
        {"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000001", "Name": "Local name", "Collection": "Local Collection"},
        {"ODC ID": "ODC_0002", "BBF ID": "BBF10K_000002", "Name": "Missing sequence", "Collection": "Local Collection"}])
    make_freegenes(tmp_path, records=[metadata_record(), metadata_record("BBF10K_000002", "Missing sequence"),
                                    metadata_record("BBF10K_000003", "External only")])
    folder = tmp_path / "genbank"; folder.mkdir()
    (folder / "ODC_0001.gb").write_bytes(gb_bytes)
    monkeypatch.setenv("OPEN_DNA_BASE_PATH", str(tmp_path))
    monkeypatch.setenv("FREEGENES_OFFLINE", "1")
    # Isolate cached services between independent deployment fixtures.
    import streamlit as st
    st.cache_resource.clear()
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    assert not app.exception
    return app


def test_navigation_home_and_debug_merge(app):
    labels = [b.label for b in app.sidebar.button]
    assert labels == ["Home", "Search & Browse", "BLAST Search", "Interactive Builder", "Debug"]
    assert not app.json
    assert all("Freshness" not in m.value and "Diagnostics" not in m.value for m in app.markdown)
    button(app, "Debug").click().run()
    assert not app.exception
    assert any("Diagnostics" in m.value for m in app.markdown)
    assert app.json
    assert {"Reclone main", "Reclone platemaps", "Local GenBank index"}.issubset(next(s for s in app.selectbox if s.label == "Dataset").options)


def test_search_open_close_reopen_and_missing_sequence(app):
    button(app, "Search & Browse").click().run()
    button(app, "Search").click().run()
    assert not app.exception
    assert len(app.session_state["search_results"]) == 2
    assert app.session_state["submitted_search"]["query"] == ""
    button(app, "View details").click().run()
    assert not app.exception
    assert app.session_state["details_open"]
    assert any(m.label == "Features" and m.value == "2" for m in app.metric)
    assert len(app.get("download_button")) == 6  # Search + five part formats.
    button(app, "Close details").click().run()
    assert not app.session_state["details_open"]
    button(app, "View details").click().run()
    assert app.session_state["details_open"]
    button(app, "Close details").click().run()
    part = next(s for s in app.selectbox if s.label == "Part on this page")
    part.select("BBF10K_000002").run()
    button(app, "View details").click().run()
    assert not app.exception
    assert any("No usable sequence" in i.value for i in app.info)
    assert len(app.get("download_button")) == 3  # Search + metadata CSV/TXT.


def test_search_state_survives_navigation_and_literal_filters(app):
    button(app, "Search & Browse").click().run()
    next(t for t in app.text_input if t.label == "Search parts").set_value("ODC-1")
    button(app, "Search").click().run()
    assert len(app.session_state["search_results"]) == 1
    button(app, "Debug").click().run()
    button(app, "Search & Browse").click().run()
    assert not app.exception
    assert app.session_state["submitted_search"]["query"] == "ODC-1"
    assert len(app.session_state["search_results"]) == 1


def test_stale_page_and_blast_page(app):
    app.session_state["current_page"] = "Analytics"
    app.run()
    assert app.session_state["current_page"] == "Search & Browse"
    button(app, "BLAST Search").click().run()
    assert not app.exception
    assert any("separate from" in m.value for m in app.markdown)


def test_dialog_hides_requested_diagnostics_and_provenance(app):
    button(app, "Search & Browse").click().run()
    button(app, "Search").click().run()
    button(app, "View details").click().run()
    assert not app.exception
    assert not app.json
    assert not app.expander or all("provenance" not in e.label.lower() for e in app.expander)
    visible = "\n".join(e.value for kind in (app.info, app.warning, app.caption) for e in kind)
    for hidden in ("Metadata pulled from", "Attempting to parse malformed locus", "LOCUS topology",
                   "No physical location records found", "distribution locations are shown separately"):
        assert hidden not in visible


def test_row_click_event_opens_exact_key_and_does_not_reopen_on_close(app):
    button(app, "Search & Browse").click().run()
    button(app, "Search").click().run()
    table_key = app.get("component_instance")[0].proto.id.rsplit("-", 1)[-1]
    app.session_state[table_key] = {"part_key": "BBF10K_000002", "event_id": "click-1"}
    app.run()
    assert app.session_state["selected_part_key"] == "BBF10K_000002"
    assert app.session_state["details_open"]
    button(app, "Close details").click().run()
    assert not app.session_state["details_open"]
    assert not app.exception


def select(app, label, value):
    return next(s for s in app.selectbox if s.label == label).select(value).run()


def test_builder_explicit_selection_generate_downloads_and_stale_inputs(app):
    button(app, "Interactive Builder").click().run()
    app.radio(key="builder_mode").set_value("Paste DNA / FASTA").run()
    app.text_area(key="builder_dna").set_value("GGTCTCAGGAGTTTTTTGCTTAGAGACC").run()
    button(app, "Analyze sequence").click().run()
    assert any("topology" in e.value for e in app.error)
    select(app, "Analysis topology", "linear")
    button(app, "Analyze sequence").click().run()
    assert not app.exception
    assert button(app, "Generate selected part").disabled
    select(app, "Digestion fragment", 1)
    assert not button(app, "Generate selected part").disabled
    button(app, "Generate selected part").click().run()
    assert not app.exception
    assert any("10 bp" in e.value for e in app.success)
    assert len(app.get("download_button")) == 8  # Three analysis reports, five generated formats.
    select(app, "Enzyme", "SapI")
    assert any("Inputs changed" in e.value for e in app.info)
    assert not app.get("download_button")


def test_builder_sapi_unmapped_and_invalid_fragment_blocks_generation(app):
    button(app, "Interactive Builder").click().run()
    app.radio(key="builder_mode").set_value("Paste DNA / FASTA").run()
    app.text_area(key="builder_dna").set_value("GCTCTTCAGCATTTTTTTGCAGAAGAGC").run()
    select(app, "Analysis topology", "linear")
    select(app, "Enzyme", "SapI")
    button(app, "Analyze sequence").click().run()
    assert not app.exception
    analysis = app.session_state["builder_analysis"]["analysis"]
    assert all("Unmapped" in site["junction"]["status"] for site in analysis["sites"])
    app.text_area(key="builder_dna").set_value("GCTCTTCAGCATTNTTTTGCAGAAGAGC").run()
    button(app, "Analyze sequence").click().run()
    select(app, "Digestion fragment", 1)
    assert button(app, "Generate selected part").disabled
    assert any("ambiguous" in e.value for e in app.error)


def test_builder_database_inventory_and_missing_sequence(app):
    button(app, "Interactive Builder").click().run()
    selector = next(s for s in app.selectbox if s.label == "Reclone part")
    assert len(selector.options) == 3  # Placeholder plus two Reclone identities; no upstream-only part.
    selector.select("BBF10K_000002").run()
    select(app, "Analysis topology", "linear")
    button(app, "Analyze sequence").click().run()
    assert not app.exception
    assert any("No usable sequence" in e.value for e in app.error)


def test_builder_database_revision_change_hides_previous_analysis(app, tmp_path):
    import json
    button(app, "Interactive Builder").click().run()
    select(app, "Reclone part", "BBF10K_000001")
    select(app, "Analysis topology", "circular")
    button(app, "Analyze sequence").click().run()
    assert not app.exception and app.session_state["builder_analysis"]
    manifest = tmp_path / "data/freegenes/manifest.json"
    value = json.loads(manifest.read_text())
    value["revision_test"] = "changed_generation"
    manifest.write_text(json.dumps(value))
    app.run()
    assert not app.exception
    assert any("Inputs changed" in e.value for e in app.info)
    assert not app.get("download_button")
