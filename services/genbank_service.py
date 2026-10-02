"""One parsed record and structured feature model for every part export/view."""
from __future__ import annotations

import hashlib
import io
import json
import re
import textwrap
import warnings
from functools import lru_cache

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class InvalidGenBank(ValueError):
    pass


def _position(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


@lru_cache(maxsize=128)
def parse_genbank(raw: bytes) -> dict:
    """Preserve original bytes; locations use BioPython's 0-based half-open convention."""
    if len(raw) > 5_000_000 or b"LOCUS" not in raw[:300] or not raw.rstrip().endswith(b"//"):
        raise InvalidGenBank("Incomplete or non-GenBank response")
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            record = SeqIO.read(io.StringIO(raw.decode("utf-8-sig")), "genbank")
            sequence = str(record.seq)
        if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence):
            raise InvalidGenBank("Missing or invalid DNA sequence")
        if any("Expected sequence length" in str(w.message) for w in captured):
            raise InvalidGenBank("Sequence length differs from the LOCUS header")
    except (ValueError, UnicodeError) as exc:
        raise InvalidGenBank(str(exc)) from exc
    features = []
    for i, feature in enumerate(record.features):
        location = feature.location
        spans = [] if location is None else [
            {"start": _position(p.start), "end": _position(p.end), "strand": p.strand,
             "start_expression": str(p.start), "end_expression": str(p.end),
             "remote_ref": p.ref} for p in location.parts
        ]
        qualifiers = {k: [str(v) for v in vals] for k, vals in feature.qualifiers.items()}
        label = next((qualifiers[k][0] for k in ("label", "gene", "product", "note")
                      if qualifiers.get(k)), feature.type)
        features.append({"index": i, "type": feature.type, "label": label,
                         "location": str(location) if location is not None else "unknown",
                         "operator": getattr(location, "operator", None),
                         "strand": location.strand if location is not None else None,
                         "spans": spans, "qualifiers": qualifiers})
    return {"record": record, "raw": raw, "sha256": hashlib.sha256(raw).hexdigest(),
            "sequence": sequence, "length": len(sequence), "gc_content": gc_fraction(sequence) * 100,
            "record_id": record.id, "description": record.description,
            "topology": record.annotations.get("topology", "unknown"),
            "annotations": record.annotations, "features": features,
            "parse_warnings": list(dict.fromkeys(str(w.message) for w in captured))}


def feature_table(features: list[dict]) -> pd.DataFrame:
    rows = []
    for f in features:
        coordinates = []
        for span in f["spans"]:
            start, end = span["start"], span["end"]
            coordinates.append(f"{start + 1}..{end}" if start is not None and end is not None else "unknown")
        rows.append({"Type": f["type"], "Label": f["label"],
                     "Bases (1-based, inclusive)": "; ".join(coordinates) or "unknown",
                     "Strand": {1: "+", -1: "-"}.get(f["strand"], "unknown/mixed"),
                     "Location (0-based, half-open)": f["location"],
                     "Qualifiers": json.dumps(f["qualifiers"], ensure_ascii=False)})
    return pd.DataFrame(rows)


def export_part(details: dict) -> dict[str, bytes]:
    """CSV/TXT work even without usable DNA; never expose absolute server paths."""
    gb = details.get("genbank") or {}
    summary = {"part_id": details["part_key"], "aliases": details["aliases"],
               "name": details["name"], "collections": details["collections"],
               "metadata": details["metadata"], "locations": details["locations"],
               "source_records": details["source_records"], "warnings": details["warnings"],
               "provenance": details["provenance"], "sequence_status": details["sequence_status"],
               "topology": gb.get("topology"), "sequence": gb.get("sequence", ""),
               "length": gb.get("length"), "description": gb.get("description"),
               "features": gb.get("features", []), "annotations": gb.get("annotations", {}),
               "genbank_sha256": gb.get("sha256")}
    row = {k: json.dumps(v, ensure_ascii=False, default=str) if isinstance(v, (dict, list)) else v
           for k, v in summary.items()}
    report = [f"Part: {details['part_key']}", f"Name: {details['name']}",
              f"Sequence status: {details['sequence_status']}", "",
              "Metadata and provenance", json.dumps({k: v for k, v in summary.items()
                  if k not in ("sequence", "features")}, indent=2, ensure_ascii=False, default=str),
              "", "Feature list (internal coordinates: 0-based, half-open)",
              json.dumps(summary["features"], indent=2, ensure_ascii=False),
              "", "DNA sequence (5' to 3')", textwrap.fill(summary["sequence"], width=80)]
    exports = {"csv": pd.DataFrame([row]).to_csv(index=False).encode("utf-8"),
               "txt": "\n".join(report).encode("utf-8")}
    if gb:
        safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", details["part_key"])
        exports.update({"gb": gb["raw"], "fasta":
            (f">{safe_id} sha256={gb['sha256']}\n" + textwrap.fill(gb["sequence"], 80) + "\n").encode(),
            "features.csv": feature_table(gb["features"]).to_csv(index=False).encode("utf-8")})
    return exports
