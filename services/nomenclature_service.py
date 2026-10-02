"""Versioned workbook labels, independent of enzyme cleavage calculations."""
import hashlib
import json
import re
from pathlib import Path

SCHEME_PATH = Path(__file__).resolve().parents[1] / "data/nomenclature/user_workbook_four_base_v1.json"


def load_scheme(path=SCHEME_PATH):
    raw = Path(path).read_bytes()
    scheme = json.loads(raw)
    if scheme["schema_version"] != 1:
        raise ValueError("Unsupported nomenclature schema")
    sequences, labels = set(), set()
    for entry in scheme["junctions"]:
        if not re.fullmatch(f"[ACGT]{{{scheme['overhang_length']}}}", entry["sequence"]):
            raise ValueError("Invalid scheme junction sequence")
        names = [entry["label"], *entry["aliases"]]
        if entry["sequence"] in sequences or any(n in labels for n in names):
            raise ValueError("Duplicate scheme sequence or label")
        sequences.add(entry["sequence"]); labels.update(names)
    return {**scheme, "config_sha256": hashlib.sha256(raw).hexdigest()}


def junction_label(word, scheme):
    if not scheme:
        return {"label": None, "aliases": [], "status": "No scheme selected"}
    if len(word) != scheme["overhang_length"]:
        return {"label": None, "aliases": [], "status": "Unmapped in selected scheme (overhang length differs)"}
    for entry in scheme["junctions"]:
        if word == entry["sequence"]:
            return {"label": entry["label"], "aliases": entry["aliases"], "status": "Mapped",
                    "source_cells": {"labels": entry["label_cells"], "sequence": entry["sequence_cell"]}}
    return {"label": None, "aliases": [], "status": "Unmapped in selected scheme"}
