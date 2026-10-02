"""Inclusive-coordinate adapter and pinned local asset integrity."""
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from services.builder_input_service import input_record, with_topology
from services.viewer_service import viewer_data

ROOT = Path(__file__).resolve().parents[2]


def feature(spans, strand=1, index=0):
    return {"spans": spans, "index": index, "label": "joined", "type": "CDS", "strand": strand}


def span(a, b, strand=1):
    return {"start": a, "end": b, "strand": strand, "remote_ref": None}


def test_compound_reverse_origin_ranges_and_linear_topology():
    gb = input_record("ACGT" * 25, "circular")["genbank"]
    gb = {**gb, "features": [feature([span(0, 8, -1), span(90, 100, -1)], -1)]}
    data = viewer_data(gb)
    mapped = data["sequenceData"]["features"][0]
    assert mapped["start"] == 90 and mapped["end"] == 7 and not mapped["forward"]
    assert mapped["locations"] == [{"start": 90, "end": 99}, {"start": 0, "end": 7}]
    linear = viewer_data(with_topology(gb, "linear"))
    assert not linear["sequenceData"]["circular"]
    assert linear["sequenceData"]["features"][0]["start"] == 0
    assert linear["content_hash"] != data["content_hash"]
    assert data["sequenceData"]["sequence"] == gb["sequence"]


def test_mixed_strands_and_unmappable_features_keep_full_list():
    gb = input_record("ACGT" * 25)["genbank"]
    features = [feature([span(10, 20, 1), span(40, 50, -1)], None),
                feature([span(None, 25)], index=1),
                feature([{**span(0, 10), "remote_ref": "external"}], index=2)]
    gb = {**gb, "features": features}
    data = viewer_data(gb)
    assert [f["forward"] for f in data["sequenceData"]["features"]] == [True, False]
    assert data["unmapped_features"] == [1, 2]
    assert len(gb["features"]) == 3


def test_unknown_strand_is_not_drawn_as_a_directional_arrow():
    gb = input_record("ACGT" * 25)["genbank"]
    gb = {**gb, "features": [feature([span(10, 20, None)], None)]}
    mapped = viewer_data(gb)["sequenceData"]["features"][0]
    assert mapped["arrowheadType"] == "NONE" and mapped["strand"] is None
    assert "strand unknown" in mapped["name"]


def test_pinned_assets_and_license_match_manifest():
    folder = ROOT / "ui/sequence_viewer/vendor"
    manifest = json.loads((folder / "manifest.json").read_text())
    assert manifest["version"] == "0.8.42" and manifest["license"] == "MIT"
    for name, sha256 in manifest["files"].items():
        assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == sha256
    html = (folder.parent / "index.html").read_text()
    assert 'src="vendor/index.umd.js"' in html and 'href="vendor/ove.css"' in html
    assert 'src="http' not in html and 'href="http' not in html
    assert "Permission is hereby granted" in (folder / "LICENSE-TeselaGen.txt").read_text()


def test_viewer_readiness_and_measurement_order():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js required for viewer protocol regression")
    result = subprocess.run([node, str(ROOT / "tests/frontend/viewer_handshake.cjs"),
        str(ROOT / "ui/sequence_viewer/index.html")], capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
