from pathlib import Path
import io

import pandas as pd
import pytest
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from services.cache_service import DNACollectionDataService
from services.freegenes_service import FreeGenesService, publish_snapshot

COMMIT = "a" * 40


@pytest.fixture
def gb_bytes():
    record = SeqRecord(Seq("ATGCGCAT" * 5), id="example", name="example", description="Annotated test plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(SimpleLocation(2, 12, strand=-1), type="CDS", qualifiers={"label": ["Reverse CDS"], "note": ["one", "two"]}),
        SeqFeature(CompoundLocation([SimpleLocation(32, 40, strand=1), SimpleLocation(0, 4, strand=1)]),
                   type="misc_feature", qualifiers={"label": ["Origin crossing"]})]
    output = io.StringIO()
    SeqIO.write(record, output, "genbank")
    return output.getvalue().encode()


def metadata_record(part_id="BBF10K_000001", name="FreeGenes name", source="FreeGenes GitHub", **fields):
    return {"part_id": part_id, "fields": {"id": part_id, "gene_name_short": name, "product": "Upstream Collection", **fields},
            "provenance": {"source": source, "revision": COMMIT, "retrieved_at": "2026-10-01T00:00:00+00:00",
                           "url": "https://example.invalid/reference", "row": 2}}


def make_freegenes(tmp_path, records=None, backend_records=None, paths=None, locations=None, client=None):
    payload = {"records": records or [], "backend_records": backend_records or [],
               "genbank_paths": paths or [], "locations": locations or [],
               "backend_status": {"status": "ok", "retrieved_at": "2026-10-01T00:00:00+00:00"},
               "location_status": {"status": "unavailable"}}
    publish_snapshot(tmp_path / "data" / "freegenes", payload,
                     {"source_commit": COMMIT, "commit_date": "2023-09-15T00:00:00Z",
                      "checked_at": "2026-10-01T00:00:00+00:00", "source_repo": "freegenes/freegenes.github.io"})
    return FreeGenesService(tmp_path, client=client)


def make_reclone(tmp_path, rows=None, locations=None):
    rows = rows if rows is not None else [{"ODC ID": "ODC_0001", "BBF ID": "BBF10K_000001", "Name": "Local name", "Collection": "Local Collection"}]
    pd.DataFrame(rows, columns=["ODC ID", "BBF ID", "Name", "Collection"]).to_csv(tmp_path / "odc_plasmids.csv", index=False)
    if locations:
        folder = tmp_path / "Local Collection" / "Platemaps"
        folder.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(locations).to_csv(folder / "v1.csv", index=False)
    return DNACollectionDataService(str(tmp_path))

CSV = b"id,gene_name_short,product,internal_notes_1,production_sequence\nBBF10K_000001,One,Kit,private,ATGC\n"

class Client:
    def __init__(self, gb=None):
        self.gb = gb
        self.calls = []
        self.fail = None
        self.backend_fail = False

    def head(self, etag=None):
        return {"sha": COMMIT, "commit": {"committer": {"date": "2023-09-15T00:00:00Z"}}}, "etag", 200

    def tree(self, commit):
        return ["product-csvs/123.csv", "genbank/BBF10K_000001.gb"]

    def genes(self):
        if self.backend_fail:
            raise requests.Timeout()
        return CSV

    def file(self, commit, path):
        self.calls.append((commit, path))
        if self.fail:
            raise self.fail
        return self.gb if path.endswith(".gb") else CSV
