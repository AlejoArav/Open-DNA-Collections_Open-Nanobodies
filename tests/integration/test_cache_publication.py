import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.sync_upstream_data import GitHubUpstreamClient, SourceMetadata, write_artifacts
from services.blast_service import BlastService
from services.cache_service import DNACollectionDataService


def test_reclone_immutable_publication_and_blast_manifest_path(tmp_path, gb_bytes):
    (tmp_path / "genbank").mkdir()
    (tmp_path / "genbank" / "ODC_0001.gb").write_bytes(gb_bytes)
    main = pd.DataFrame({"BBF ID": ["BBF10K_000001"], "ODC ID": ["ODC_0001"], "Name": ["One"],
                         "Collection": ["Kit"], "ODC_ID_NORM": ["ODC_0001"], "BBF_ID_NORM": ["BBF10K_000001"]})
    source = SourceMetadata("repo", "branch", "a" * 40)
    manifest = write_artifacts(tmp_path, tmp_path / "data" / "cache", main, pd.DataFrame(), source)
    assert all("/snapshots/" in f["path"].replace("\\", "/") for f in manifest["files"])
    assert DNACollectionDataService(str(tmp_path)).main_df.iloc[0]["Name"] == "One"
    blast = BlastService(str(tmp_path))
    assert blast.local_fasta_path.exists()
    assert "snapshots" in blast.local_fasta_path.parts


def test_reclone_failed_generation_preserves_previous_pointer(tmp_path, monkeypatch):
    main = pd.DataFrame({"BBF ID": ["BBF10K_000001"], "ODC ID": ["ODC_0001"], "Name": ["One"],
                         "Collection": ["Kit"], "ODC_ID_NORM": ["ODC_0001"], "BBF_ID_NORM": ["BBF10K_000001"]})
    source = SourceMetadata("repo", "branch", "a" * 40)
    cache = tmp_path / "data" / "cache"
    write_artifacts(tmp_path, cache, main, pd.DataFrame(), source)
    before = (cache / "manifest.json").read_bytes()
    def fail(*args, **kwargs):
        raise OSError("simulated generation failure")
    monkeypatch.setattr("scripts.sync_upstream_data.build_genbank_assets", fail)
    with pytest.raises(OSError):
        write_artifacts(tmp_path, cache, main.assign(Name="New"), pd.DataFrame(), source)
    assert (cache / "manifest.json").read_bytes() == before
    assert DNACollectionDataService(str(tmp_path)).main_df.iloc[0]["Name"] == "One"


def test_upstream_files_pinned_to_tree_commit(monkeypatch):
    client = GitHubUpstreamClient("org/repo", "main")
    monkeypatch.setattr(client, "_get_json", lambda *a, **k: {"tree": [{"type": "blob", "path": "odc_plasmids.csv"}], "truncated": False})
    client.list_target_files("a" * 40)
    class Response:
        text = "BBF ID,Name\nBBF10K_000001,One\n"
        def raise_for_status(self):
            pass
    observed = []
    monkeypatch.setattr(client.session, "get", lambda url, **kwargs: observed.append(url) or Response())
    client.fetch_text_file("odc_plasmids.csv")
    assert "/" + "a" * 40 + "/" in observed[0]
    assert "/main/" not in observed[0]


def test_upstream_truncated_tree_fails(monkeypatch):
    client = GitHubUpstreamClient("org/repo", "main")
    monkeypatch.setattr(client, "_get_json", lambda *a, **k: {"truncated": True, "tree": []})
    with pytest.raises(ValueError, match="truncated"):
        client.list_target_files("a" * 40)
