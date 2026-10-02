import io
import json

import pandas as pd
import pytest
from Bio import SeqIO

from services.genbank_service import InvalidGenBank, export_part, feature_table, parse_genbank


def test_compound_reverse_and_export_consistency(gb_bytes):
    gb = parse_genbank(gb_bytes)
    assert gb["raw"] == gb_bytes
    assert gb["topology"] == "circular"
    assert gb["features"][0]["spans"] == [{"start": 2, "end": 12, "strand": -1,
                                         "start_expression": "2", "end_expression": "12", "remote_ref": None}]
    assert gb["features"][1]["operator"] == "join"
    assert [(p["start"], p["end"]) for p in gb["features"][1]["spans"]] == [(32, 40), (0, 4)]
    assert feature_table(gb["features"]).iloc[0]["Bases (1-based, inclusive)"] == "3..12"
    details = {"part_key": "BBF10K_000001", "aliases": ["ODC_0001", "BBF10K_000001"], "name": "Test",
               "collections": [], "metadata": {}, "locations": [], "source_records": [], "warnings": [],
               "provenance": {}, "sequence_status": "direct", "genbank": gb}
    exports = export_part(details)
    assert exports["gb"] == gb_bytes
    fasta = SeqIO.read(io.StringIO(exports["fasta"].decode()), "fasta")
    assert str(fasta.seq) == gb["sequence"]
    row = pd.read_csv(io.BytesIO(exports["csv"])).iloc[0]
    assert row["sequence"] == gb["sequence"]
    assert json.loads(row["features"]) == gb["features"]
    assert row["genbank_sha256"] == gb["sha256"]
    assert gb["sha256"] in exports["txt"].decode()
    assert gb["sequence"] in exports["txt"].decode()


@pytest.mark.parametrize("raw", [b"<html>error</html>", b"LOCUS invalid", b"LOCUS x\nORIGIN\n//"])
def test_invalid_file_is_rejected(raw):
    with pytest.raises(InvalidGenBank):
        parse_genbank(raw)
