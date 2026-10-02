"""Project an explicitly selected reference-strand fragment and reparse its exports."""
from __future__ import annotations

import base64
import hashlib
import io
import json
import textwrap
from copy import deepcopy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, ExactPosition, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .genbank_service import export_part, parse_genbank
from .restriction_service import candidate_between


def remap_features(record, start, length, topology):
    """Retain biological part order; clipping endpoints become exact cut boundaries."""
    n = len(record.seq)
    segments = [(start, min(start + length, n), 0)]
    if start + length > n:
        if topology != "circular":
            raise ValueError("A linear fragment cannot cross the origin")
        segments.append((0, start + length - n, n - start))
    remapped, notices = [], []
    for index, feature in enumerate(record.features):
        location = feature.location
        if location is None or any(p.ref or p.ref_db for p in location.parts):
            notices.append(f"Feature {index} ({feature.type}) has an unknown/remote location and was omitted")
            continue
        try:
            if any(not 0 <= int(p.start) <= int(p.end) <= n for p in location.parts):
                raise ValueError("out-of-range feature")
            original_length = len(location)
        except (TypeError, ValueError):
            notices.append(f"Feature {index} ({feature.type}) has an unmappable location and was omitted")
            continue
        parts, retained, removed_prefix, started = [], 0, 0, False
        for part in location.parts:
            a, b = int(part.start), int(part.end)
            intersections = [(max(a, lo), min(b, hi), lo, offset)
                             for lo, hi, offset in segments if max(a, lo) < min(b, hi)]
            intersections.sort(key=lambda v: v[0], reverse=part.strand == -1)
            if not started:
                if intersections:
                    first = intersections[0]
                    removed_prefix += b - first[1] if part.strand == -1 else first[0] - a
                    started = True
                else:
                    removed_prefix += b - a
            for x, y, lo, offset in intersections:
                delta = offset - lo
                new_start = part.start + delta if x == a else ExactPosition(x + delta)
                new_end = part.end + delta if y == b else ExactPosition(y + delta)
                parts.append(SimpleLocation(new_start, new_end, strand=part.strand))
                retained += y - x
        if not parts:
            continue
        mapped = deepcopy(feature)
        mapped.location = parts[0] if len(parts) == 1 else CompoundLocation(
            parts, operator=getattr(location, "operator", "join"))
        if retained != original_length:
            mapped.qualifiers.setdefault("note", []).append(
                f"Partial feature retained after digestion; source location {location}")
            if mapped.type == "CDS":
                mapped.qualifiers.pop("translation", None)
                try:
                    old_frame = int(mapped.qualifiers.get("codon_start", ["1"])[0])
                    exact_source = all(isinstance(p.start, ExactPosition) and isinstance(p.end, ExactPosition)
                                       for p in location.parts)
                    if exact_source and location.strand in (-1, 1) and old_frame in (1, 2, 3):
                        mapped.qualifiers["codon_start"] = [str((old_frame - 1 - removed_prefix) % 3 + 1)]
                    else:
                        mapped.qualifiers.pop("codon_start", None)
                except (ValueError, TypeError):
                    mapped.qualifiers.pop("codon_start", None)
            notices.append(f"Feature {index} ({feature.type}) was clipped; stale CDS translation removed when present")
        remapped.append(mapped)
    return remapped, notices


def generate_fragment(gb, analysis, candidate, source=None):
    """Revalidate boundaries and source identity; never trust a stale UI candidate."""
    if str(gb["record"].seq).upper() != analysis["sequence"]:
        raise ValueError("Analysis does not match the source sequence")
    selected = candidate_between(analysis, candidate["left_cut_id"], candidate["right_cut_id"])
    if not selected["valid"]:
        raise ValueError("; ".join(selected["errors"]))
    features, notices = remap_features(gb["record"], selected["top_start"], selected["top_length"], analysis["topology"])
    identity = hashlib.sha256((gb["sha256"] + analysis["enzyme"] + analysis["topology"] +
                              selected["id"]).encode()).hexdigest()[:12]
    part_id = "fragment_" + identity
    provenance = {"parent_id": gb["record_id"], "parent_genbank_sha256": gb["sha256"],
        "parent_sequence_sha256": analysis["sequence_sha256"], "source": source or {},
        "analysis_topology": analysis["topology"], "enzyme": analysis["enzyme"],
        "enzyme_reference": analysis["enzyme_spec"]["reference"],
        "scheme_id": (analysis["scheme"] or {}).get("id"),
        "scheme_config_sha256": (analysis["scheme"] or {}).get("config_sha256"),
        "fragment": {k: v for k, v in selected.items() if k != "sequence"},
        "strand_convention": "Reference top strand, 5prime to 3prime: includes left overhang; excludes right bottom-strand overhang"}
    # GenBank COMMENT wrapping inserts newlines into long JSON strings. Encode
    # and pre-wrap it so provenance survives independent GenBank round-trips.
    encoded = base64.b64encode(json.dumps(provenance, ensure_ascii=True,
                                        separators=(",", ":")).encode()).decode()
    comment = "OpenDNA provenance (base64 UTF-8 JSON):\n" + "\n".join(textwrap.wrap(encoded, 60))
    record = SeqRecord(Seq(selected["sequence"]), id=part_id, name=part_id[:16],
                       description=f"Selected {analysis['enzyme']} digestion fragment of {gb['record_id']}",
                       features=features, annotations={"molecule_type": "DNA", "topology": "linear",
                           "comment": comment})
    buffer = io.StringIO()
    SeqIO.write(record, buffer, "genbank")
    parsed = parse_genbank(buffer.getvalue().encode("utf-8"))
    if parsed["sequence"] != selected["sequence"] or parsed["topology"] != "linear":
        raise ValueError("Generated GenBank did not round-trip consistently")
    # Independently check that serialization retained mapped biological locations.
    if len(parsed["record"].features) != len(features) or any(
            a.extract(record.seq) != b.extract(parsed["record"].seq)
            for a, b in zip(features, parsed["record"].features)):
        raise ValueError("Generated feature locations did not round-trip consistently")
    details = {"part_key": part_id, "aliases": [part_id], "name": record.description,
               "collections": [], "metadata": {"enzyme": analysis["enzyme"],
                   "left_end": selected["left_end"], "right_end": selected["right_end"]},
               "locations": [], "source_records": [], "warnings": notices,
               "provenance": provenance, "sequence_status": "Generated validated digestion fragment",
               "genbank": parsed}
    return {"details": details, "candidate": selected, "exports": export_part(details)}


def analysis_report(analysis, source=None, selected=None):
    """Portable facts and candidate reports, also available when generation is invalid."""
    report = {k: v for k, v in analysis.items() if k != "sequence"}
    report["source"] = source or {}
    if selected is not None:
        report["selected_fragment"] = selected
    return json.dumps(report, indent=2, ensure_ascii=False).encode("utf-8")
