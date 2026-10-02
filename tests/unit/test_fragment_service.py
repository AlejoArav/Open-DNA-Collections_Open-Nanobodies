import base64
import hashlib
import io
import json
from copy import deepcopy

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from services.builder_input_service import input_record, with_topology
from services.fragment_service import analysis_report, generate_fragment, remap_features
from services.genbank_service import parse_genbank
from services.restriction_service import analyze_restriction
from test_part_analysis import BSA


def parsed(record):
    record.annotations["molecule_type"] = "DNA"
    output = io.StringIO()
    SeqIO.write(record, output, "genbank")
    return parse_genbank(output.getvalue().encode())


def test_single_record_inputs_and_topology_override_are_immutable():
    gb = input_record(">sample\nacgtry\n")["genbank"]
    assert gb["sequence"] == "ACGTRY" and gb["topology"] == "unknown"
    copy = with_topology(gb, "circular")
    assert copy["topology"] == "circular" and gb["topology"] == "unknown"
    assert copy["raw"] == gb["raw"] and "topology" not in gb["record"].annotations
    for text in (">a\nACGT\n>b\nAAAA", "ACGT123", "", "ACGU"):
        with pytest.raises(ValueError):
            input_record(text)
    with pytest.raises(ValueError, match="500,000"):
        input_record("A" * 500_001)


def test_complete_reverse_and_compound_features_preserve_extraction():
    record = SeqRecord(Seq(BSA), id="source", annotations={"topology": "linear"})
    record.features = [SeqFeature(SimpleLocation(8, 16, strand=-1), type="CDS"),
        SeqFeature(CompoundLocation([SimpleLocation(8, 10, strand=-1), SimpleLocation(13, 16, strand=-1)]), type="misc_feature"),
        SeqFeature(SimpleLocation(0, 28, strand=1), type="source"),
        SeqFeature(SimpleLocation(BeforePosition(9), 14, strand=1), type="misc_feature")]
    original = deepcopy(record)
    gb = parsed(record)
    source_raw = gb["raw"]
    analysis = analyze_restriction(BSA, "BsaI", "linear")
    generated = generate_fragment(gb, analysis, analysis["products"][1])
    output = generated["details"]["genbank"]
    assert output["sequence"] == "GGAGTTTTTT"
    for index in (0, 1, 3):
        assert output["record"].features[index].extract(output["record"].seq) == original.features[index].extract(original.seq)
    assert isinstance(output["record"].features[3].location.start, BeforePosition)
    assert gb["raw"] == source_raw and gb["sequence"] == BSA
    assert str(record.seq) == BSA and record.features[2].location == original.features[2].location
    assert pd.read_csv(io.BytesIO(generated["exports"]["csv"])).iloc[0]["sequence"] == output["sequence"]
    assert str(SeqIO.read(io.StringIO(generated["exports"]["fasta"].decode()), "fasta").seq) == output["sequence"]
    assert output["sequence"] in generated["exports"]["txt"].decode()
    assert output["sha256"] == hashlib.sha256(generated["exports"]["gb"]).hexdigest()
    comment = output["record"].annotations["comment"]
    provenance = json.loads(base64.b64decode("".join(comment.splitlines()[1:])))
    assert provenance["parent_genbank_sha256"] == gb["sha256"]
    assert generated["details"]["locations"] == []


@pytest.mark.parametrize("strand", [1, -1])
def test_clipped_cds_removes_translation_and_adjusts_frame(strand):
    record = SeqRecord(Seq("ACGT" * 10), id="s")
    record.features = [SeqFeature(SimpleLocation(5, 25, strand=strand), type="CDS",
                                 qualifiers={"translation": ["STALE"], "codon_start": ["1"]})]
    mapped, notices = remap_features(record, 7, 16, "linear")
    assert str(mapped[0].location) == f"[0:16]({'+' if strand == 1 else '-'})"
    assert "translation" not in mapped[0].qualifiers
    assert mapped[0].qualifiers["codon_start"] == ["2"]  # Two 5prime bases removed, either strand.
    assert notices and record.features[0].qualifiers["translation"] == ["STALE"]


@pytest.mark.parametrize("strand", [1, -1])
def test_origin_crossing_join_keeps_biological_order(strand):
    record = SeqRecord(Seq("ACGT" * 25), id="circle")
    parts = [SimpleLocation(90, 98, strand=strand), SimpleLocation(2, 8, strand=strand)]
    if strand == -1:
        parts.reverse()
    record.features = [SeqFeature(CompoundLocation(parts), type="CDS")]
    mapped, notices = remap_features(record, 85, 30, "circular")
    fragment = record.seq[85:] + record.seq[:15]
    assert mapped[0].extract(fragment) == record.features[0].extract(record.seq)
    assert not notices


def test_unmapped_remote_features_and_stale_invalid_analysis_are_rejected():
    record = SeqRecord(Seq(BSA), id="s")
    record.features = [SeqFeature(SimpleLocation(1, 5, ref="remote")), SeqFeature(None)]
    mapped, notices = remap_features(record, 7, 10, "linear")
    assert not mapped and len(notices) == 2
    gb = input_record(BSA)["genbank"]
    analysis = analyze_restriction(BSA.replace("TTTTTT", "TTNTTT"), "BsaI", "linear")
    with pytest.raises(ValueError, match="match"):
        generate_fragment(gb, analysis, analysis["products"][1])
    gb = input_record(analysis["sequence"])["genbank"]
    with pytest.raises(ValueError, match="ambiguous"):
        generate_fragment(gb, analysis, analysis["products"][1])


def test_fuzzy_clipped_cds_does_not_invent_a_reading_frame():
    record = SeqRecord(Seq("ACGT" * 10), id="s")
    record.features = [SeqFeature(SimpleLocation(BeforePosition(5), 25, strand=1), type="CDS",
                                 qualifiers={"translation": ["STALE"], "codon_start": ["1"]})]
    features, notices = remap_features(record, 7, 16, "linear")
    assert notices and "translation" not in features[0].qualifiers
    assert "codon_start" not in features[0].qualifiers


def test_invalid_candidate_report_keeps_boundaries_errors_and_hashes():
    analysis = analyze_restriction(BSA.replace("TTTTTT", "TTNTTT"), "BsaI", "linear")
    selected = analysis["products"][1]
    report = json.loads(analysis_report(analysis, {"input_sha256": "test"}, selected))
    assert not report["selected_fragment"]["valid"]
    assert report["selected_fragment"]["left_cut_id"] == "cut_7_11"
    assert report["selected_fragment"]["errors"] == selected["errors"]
    assert report["sequence_sha256"] == analysis["sequence_sha256"]


def test_checked_in_origin_fixture_reparse_preserves_retained_join():
    from pathlib import Path
    gb = parse_genbank((Path(__file__).resolve().parents[1] / "fixtures/goal2_circle.gb").read_bytes())
    analysis = analyze_restriction(gb["sequence"], "BsaI", "circular")
    candidate = next(p for p in analysis["products"] if not p["retained_recognition_sites"])
    generated = generate_fragment(gb, analysis, candidate)
    record = generated["details"]["genbank"]["record"]
    joined = next(f for f in record.features if f.qualifiers["label"] == ["Retained reverse join"])
    assert str(joined.extract(record.seq)) == "AAATC"
    assert len(joined.location.parts) == 2 and joined.location.strand == -1
