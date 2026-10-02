from pathlib import Path

import pandas as pd

from services.cache_service import DNACollectionDataService
from services.freegenes_service import FreeGenesService
from services.part_service import PartService, identifier


def test_search_enrichment_works_with_local_fallback(tmp_path):
    # Build minimal local repo-like structure without cache manifest.
    main_df = pd.DataFrame(
        {
            "BBF ID": ["BBF10K_003247"],
            "ODC ID": ["ODC_0007"],
            "Name": ["9°N-7 DNA polymerase"],
            "Collection": ["Open Enzyme Collection"],
        }
    )
    main_df.to_csv(tmp_path / "odc_plasmids.csv", index=False)

    platemap_dir = tmp_path / "Open Enzyme Collection" / "Platemaps"
    platemap_dir.mkdir(parents=True)
    platemap_df = pd.DataFrame(
        {
            "Well Location": ["A1"],
            "ODC ID": ["ODC_0007"],
            "BBF ID": ["BBF10K_003247"],
            "Bacterial Resistance": ["ampicillin"],
            "Growth Strain": ["DH5a"],
            "Growth Conditions": ["LB"],
            "Name": ["9°N-7 DNA polymerase"],
        }
    )
    platemap_df.to_csv(platemap_dir / "OEC-v1_1.csv", index=False)

    service = DNACollectionDataService(base_path=str(tmp_path))
    result = service.search_parts(query="ODC_0007")

    assert len(result) == 1
    assert result.iloc[0]["Well_Location"] == "A1"
    assert result.iloc[0]["Bacterial_Resistance"] == "ampicillin"
    assert result.iloc[0]["Growth_Strain"] == "DH5a"


def test_checked_in_inventory_has_complete_name_and_polymerase_search_coverage():
    # This is a read-only audit of shipped caches; no sequence fetch or upstream request.
    root = Path(__file__).resolve().parents[2]
    local = DNACollectionDataService(str(root))
    parts = PartService(local, FreeGenesService(root))
    original = pd.read_csv(root / "odc_plasmids.csv").fillna("")
    pd.testing.assert_frame_equal(original, local.main_df[original.columns].fillna(""))
    for row in local.main_df.to_dict("records"):
        key = parts.aliases[identifier(row["BBF ID"]) or identifier(row["ODC ID"])]
        result = parts.search_parts(row["Name"]).set_index("Part Key")
        assert key in result.index, row["Name"]
        assert row["Name"] in result.loc[key, "Name"]
    expected = local.main_df[local.main_df["Name"].str.contains("polymerase", case=False)]
    found = parts.search_parts("polymerase").set_index("Part Key")
    for row in expected.to_dict("records"):
        assert found.loc[parts.aliases[identifier(row["BBF ID"])], "Name"] == row["Name"]
    for part_id in ("BBF10K_000483", "BBF10K_003338"):
        assert parts.search_parts(part_id).iloc[0]["Sources"].startswith("Reclone")
