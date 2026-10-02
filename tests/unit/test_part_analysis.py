"""Hand-calculated cuts, duplex ends, ambiguity, and source-backed labels."""
import pytest

from services.nomenclature_service import junction_label, load_scheme
from services.restriction_service import analyze_restriction, candidate_between

BSA = "GGTCTC" + "A" + "GGAG" + "TTTTTT" + "GCTT" + "A" + "GAGACC"
SAP = "GCTCTTC" + "A" + "GCA" + "TTTTTT" + "TGC" + "A" + "GAAGAGC"


@pytest.mark.parametrize("enzyme,dna,expected,fragment,bottom", [
    ("BsaI", BSA, [(0, 1, 7, 11, "GGAG"), (22, -1, 17, 21, "GCTT")], "GGAGTTTTTT", "AAGCAAAAAA"),
    ("SapI", SAP, [(0, 1, 8, 11, "GCA"), (21, -1, 17, 20, "TGC")], "GCATTTTTT", "GCAAAAAAA"),
])
def test_hand_calculated_both_strands_and_physical_ends(enzyme, dna, expected, fragment, bottom):
    result = analyze_restriction(dna, enzyme, "linear", load_scheme())
    assert [(s["recognition_start"], s["strand"], s["top_cut"], s["bottom_cut"], s["reference_overhang_5to3"])
            for s in result["sites"]] == expected
    assert len(result["products"]) == 3
    insert = result["products"][1]
    assert insert["sequence"] == fragment
    assert insert["bottom_strand_5to3"] == bottom
    assert insert["core_length"] == 6
    assert insert["valid"] and insert["retained_recognition_sites"] == 0
    assert insert["left_end"]["physical_5prime"] == expected[0][-1]
    assert insert["right_end"]["physical_5prime"] == bottom[:result["overhang_length"]]
    assert result["products"][0]["left_end"]["kind"] == "blunt"
    assert result["products"][2]["right_end"]["kind"] == "blunt"


@pytest.mark.parametrize("enzyme,dna,offset", [("BsaI", BSA, 10), ("SapI", SAP, 10)])
def test_origin_crossing_overhang_and_fragment(enzyme, dna, offset):
    rotated = dna[offset:] + dna[:offset]
    result = analyze_restriction(rotated, enzyme, "circular")
    forward = next(s for s in result["sites"] if s["strand"] == 1)
    assert forward["bottom_cut"] < forward["top_cut"]
    insert = next(p for p in result["products"] if p["retained_recognition_sites"] == 0)
    assert insert["wraps_origin"] and insert["valid"]
    assert insert["sequence"] == ("GGAG" if enzyme == "BsaI" else "GCA") + "TTTTTT"


def test_origin_crossing_recognition_and_single_circle_cut():
    dna = BSA[3:] + BSA[:3]
    result = analyze_restriction(dna, "BsaI", "circular")
    assert next(s for s in result["sites"] if s["strand"] == 1)["recognition_wraps"]
    single = analyze_restriction("GGTCTCAACGT", "BsaI", "circular")
    assert len(single["products"]) == 1
    assert single["products"][0]["top_length"] == 11
    assert single["products"][0]["bottom_strand_5to3"] == "ACGTTGAGACC"


def test_linear_end_incomplete_and_unknown_topology():
    for dna in ("GGTCTCAA", "GAGACCAAAAAAAA"):
        result = analyze_restriction(dna, "BsaI", "linear")
        assert result["sites"][0]["status"] == "incomplete linear end"
        assert not result["cuts"] and not result["products"]
    with pytest.raises(ValueError, match="topology"):
        analyze_restriction(BSA, "BsaI", "unknown")


def test_ambiguous_recognition_overhang_and_fragment_are_not_validated():
    uncertain = analyze_restriction("GGTNTCAGGAGTTTTTT", "BsaI", "linear")
    assert uncertain["sites"][0]["status"] == "ambiguous recognition/overhang"
    assert not uncertain["cuts"]
    overhang = analyze_restriction(BSA.replace("GGAG", "GNAG"), "BsaI", "linear")
    assert any(s["status"] == "ambiguous recognition/overhang" for s in overhang["sites"])
    middle = analyze_restriction(BSA.replace("TTTTTT", "TTNTTT"), "BsaI", "linear")
    assert not middle["products"][1]["valid"]
    assert "ambiguous" in middle["products"][1]["errors"][0]


def test_custom_boundaries_reject_internal_complete_and_possible_sites():
    for inner in ("GGTCTCAGGAGTTTTTT", "GGTNTCAGGAGTTTTTT"):
        result = analyze_restriction(BSA.replace("TTTTTT", inner), "BsaI", "linear")
        candidate = candidate_between(result, result["cuts"][0]["id"], result["cuts"][-1]["id"])
        assert candidate["internal_cut_sites"] and not candidate["valid"]
    with pytest.raises(ValueError, match="follow"):
        candidate_between(result, result["cuts"][-1]["id"], result["cuts"][0]["id"])


def test_workbook_provenance_exact_labels_and_unmapped_sapi():
    scheme = load_scheme()
    assert scheme["source"]["sha256"] == "f45eac44b2858fad7e74d20739ffd179a046aa32d4e9dc663ca0aa08eaa7786d"
    a = junction_label("GGAG", scheme)
    assert a["label"] == "A" and a["source_cells"]["sequence"] == "B6"
    assert junction_label("CCAT", scheme)["aliases"] == ["N1"]
    for word in ("GCA", "CTCC", "AAAA"):
        assert junction_label(word, scheme)["label"] is None
        assert "Unmapped" in junction_label(word, scheme)["status"]


def test_cut_at_linear_origin_reports_degenerate_product_without_crashing():
    result = analyze_restriction("GGAGAGAGACC", "BsaI", "linear")
    assert result["cuts"][0]["top_cut"] == 0
    assert result["products"][0]["top_length"] == 0
    assert not result["products"][0]["valid"]
    assert result["products"][1]["valid"]
    unchanged = candidate_between(result, None, None)
    assert not unchanged["valid"] and "At least one" in unchanged["errors"][0]


def test_ambiguous_match_explosion_is_bounded():
    with pytest.raises(ValueError, match="Too many possible"):
        analyze_restriction("N" * 10_020, "BsaI", "linear")


def test_overlapping_staggered_cuts_and_coincident_site_pair():
    overlapping = analyze_restriction("GGTCTCAGGAGAAGAGACC", "BsaI", "linear")
    assert [(c["top_cut"], c["bottom_cut"]) for c in overlapping["cuts"]] == [(7, 11), (8, 12)]
    middle = overlapping["products"][1]
    assert middle["top_length"] == 1 and not middle["valid"]
    assert "no positive double-stranded core" in middle["errors"][0]
    coincident = analyze_restriction("GGTCTCAGGAGAGAGACC", "BsaI", "circular")
    assert len(coincident["sites"]) == 2 and len(coincident["cuts"]) == 1
    assert len(coincident["cuts"][0]["site_ids"]) == 2
    assert len(coincident["products"]) == 1
