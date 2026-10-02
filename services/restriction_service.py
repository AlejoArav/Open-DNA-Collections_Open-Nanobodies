"""Deterministic Type IIS cuts and explicit digestion-fragment boundaries.

All coordinates are 0-based boundaries on the reference strand. For both enzymes,
top_cut < bottom_cut in unwrapped coordinates: the right fragment has the reference
overhang as its physical 5' end; the left fragment has its reverse complement.
"""
from __future__ import annotations

import hashlib
import re

from Bio.Seq import Seq

from .nomenclature_service import junction_label

ENZYMES = {
    "BsaI": {"motif": "GGTCTC", "top_offset": 1, "bottom_offset": 5,
             "reference": "https://www.neb.com/en/tools-and-resources/selection-charts/enzymes-with-nonpalindromic-sequences"},
    "SapI": {"motif": "GCTCTTC", "top_offset": 1, "bottom_offset": 4,
             "reference": "https://www.neb.com/en-ca/products/r0569-sapi"},
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def sequence_slice(sequence, start, length, topology):
    if length < 0 or length > len(sequence):
        raise ValueError("Invalid sequence span")
    if topology == "circular":
        start %= len(sequence)
        return (sequence + sequence)[start:start+length]
    if start < 0 or start + length > len(sequence):
        raise ValueError("Sequence span crosses a linear end")
    return sequence[start:start+length]


def _matches(sequence, motif, topology, max_matches=10_000):
    n, m = len(sequence), len(motif)
    if n < m:
        return []
    scan = sequence + (sequence[:m-1] if topology == "circular" else "")
    exact, cursor = set(), 0
    while True:
        start = scan.find(motif, cursor)
        if start == -1 or start >= n:
            break
        exact.add(start); cursor = start + 1
        if len(exact) > max_matches:
            raise ValueError("Too many recognition matches (limit 10,000 per orientation); use a shorter, more specific sequence")
    # Only windows touching ambiguous bases need the more expensive IUPAC check.
    possible = set()
    for i, base in enumerate(sequence):
        if base not in "ACGT":
            for offset in range(m):
                start = (i-offset) % n if topology == "circular" else i-offset
                if 0 <= start < n and start+m <= len(scan):
                    possible.add(start)
    matches = [(s, True) for s in exact]
    for start in possible - exact:
        if all(expected in IUPAC[actual] for actual, expected in zip(scan[start:start+m], motif)):
            matches.append((start, False))
            if len(matches) > max_matches:
                raise ValueError("Too many possible recognition matches (limit 10,000 per orientation); resolve ambiguous DNA or shorten the input")
    return matches


def analyze_restriction(sequence, enzyme, topology, scheme=None):
    sequence = re.sub(r"\s+", "", sequence).upper()
    if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
        raise ValueError("Enter a nonempty DNA sequence using IUPAC bases")
    if len(sequence) > 500_000:
        raise ValueError("Analysis is limited to 500,000 bases per record")
    if topology not in ("linear", "circular"):
        raise ValueError("Choose linear or circular topology explicitly")
    if enzyme not in ENZYMES:
        raise ValueError("Choose BsaI or SapI")
    spec, n = ENZYMES[enzyme], len(sequence)
    m, width = len(spec["motif"]), spec["bottom_offset"] - spec["top_offset"]
    sites = []
    for strand, motif in ((1, spec["motif"]), (-1, reverse_complement(spec["motif"]))):
        for start, certain in _matches(sequence, motif, topology):
            # Reverse sites cut upstream on the reference strand, with offsets swapped.
            top = start + m + spec["top_offset"] if strand == 1 else start - spec["bottom_offset"]
            bottom = start + m + spec["bottom_offset"] if strand == 1 else start - spec["top_offset"]
            geometry_ok = n > m + width if topology == "circular" else 0 <= top <= bottom <= n
            word = sequence_slice(sequence, top, width, topology) if geometry_ok else ""
            status = "complete" if certain and geometry_ok and set(word) <= set("ACGT") else (
                "incomplete linear end" if not geometry_ok and topology == "linear" else
                "record too short" if not geometry_ok else "ambiguous recognition/overhang")
            sites.append({"id": f"site_{start}_{strand}", "recognition_start": start,
                "recognition_length": m, "recognition_sequence": motif, "strand": strand,
                "recognition_wraps": start+m > n,
                "top_cut": top % n if topology == "circular" else top,
                "bottom_cut": bottom % n if topology == "circular" else bottom,
                "top_cut_unwrapped": top, "bottom_cut_unwrapped": bottom,
                "reference_overhang_5to3": word, "recognition_oriented_overhang_5to3": word if strand == 1 else reverse_complement(word),
                "right_fragment_physical_5prime": word, "left_fragment_physical_5prime": reverse_complement(word),
                "junction": junction_label(word, scheme), "status": status, "certain_recognition": certain})
    sites.sort(key=lambda s: (s["recognition_start"], -s["strand"]))
    groups = {}
    for site in sites:
        if site["status"] == "complete":
            identity = (site["top_cut"], site["bottom_cut"])
            group = groups.setdefault(identity, {"id": f"cut_{identity[0]}_{identity[1]}",
                "top_cut": identity[0], "bottom_cut": identity[1],
                "reference_overhang_5to3": site["reference_overhang_5to3"], "site_ids": []})
            group["site_ids"].append(site["id"])
    result = {"sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "length": n, "enzyme": enzyme, "enzyme_spec": spec, "topology": topology,
        "overhang_length": width, "sites": sites, "cuts": sorted(groups.values(), key=lambda g: g["top_cut"]),
        "scheme": scheme, "warnings": [], "products": []}
    if any(s["status"] != "complete" for s in sites):
        result["warnings"].append("Some recognition sites have incomplete or ambiguous cuts; inspect the site table.")
    cuts = result["cuts"]
    if cuts:
        pairs = list(zip([None, *cuts], [*cuts, None])) if topology == "linear" else list(zip(cuts, [*cuts[1:], cuts[0]]))
        result["products"] = [candidate_between(result, a["id"] if a else None, b["id"] if b else None)
                              for a, b in pairs]
    return result


def _end(word, side, scheme):
    physical = word if side == "left" else reverse_complement(word)
    return {"kind": "5prime overhang" if word else "blunt", "reference_junction_5to3": word,
            "physical_5prime": physical, "physical_reverse_complement": reverse_complement(physical),
            "junction": junction_label(word, scheme) if word else {"label": None, "status": "Blunt end", "aliases": []},
            "physical_word_label": junction_label(physical, scheme) if word else None}


def candidate_between(analysis, left_id, right_id):
    cuts = {cut["id"]: cut for cut in analysis["cuts"]}
    if (left_id is not None and left_id not in cuts) or (right_id is not None and right_id not in cuts):
        raise ValueError("Unknown cut boundary")
    left, right = cuts.get(left_id), cuts.get(right_id)
    n, circular = analysis["length"], analysis["topology"] == "circular"
    if circular and (not left or not right):
        raise ValueError("Circular fragments require two cut boundaries (or one full-length linearizing cut)")
    start, end = left["top_cut"] if left else 0, right["top_cut"] if right else n
    length = (end-start) % n if circular else end-start
    if circular and length == 0:
        length = n
    if length < 0:
        raise ValueError("The right boundary must follow the left boundary on linear DNA")
    sequence = sequence_slice(analysis["sequence"], start, length, analysis["topology"])
    left_word = left["reference_overhang_5to3"] if left else ""
    right_word = right["reference_overhang_5to3"] if right else ""
    chosen_sites = set((left or {}).get("site_ids", []) + (right or {}).get("site_ids", []))
    internal = []
    for site in analysis["sites"]:
        if site["id"] in chosen_sites or site["status"] in ("incomplete linear end", "record too short"):
            continue
        distance = (site["top_cut"] - start) % n if circular else site["top_cut"] - start
        if 0 < distance < length:
            internal.append(site["id"])
    errors = []
    if not left and not right:
        errors.append("At least one enzyme cut must bound the selected fragment")
    if length <= len(left_word):
        errors.append("Staggered cuts overlap or coincide: no positive double-stranded core")
    if internal:
        errors.append("Additional complete or possible cut sites lie inside the selected fragment")
    if set(sequence + left_word + right_word) - set("ACGT"):
        errors.append("The fragment or its physical ends contain ambiguous DNA bases")
    retained_motifs = _matches(sequence, analysis["enzyme_spec"]["motif"], "linear") + _matches(
        sequence, reverse_complement(analysis["enzyme_spec"]["motif"]), "linear")
    # The exported reference strand includes its left overhang, excludes the right
    # bottom-strand overhang. Physical ends and both-strand coordinates travel with it.
    bottom_start = start + len(left_word)
    bottom_length = length + len(right_word) - len(left_word)
    bottom = reverse_complement(sequence_slice(analysis["sequence"], bottom_start, bottom_length,
                                             analysis["topology"])) if bottom_length >= 0 and not errors else ""
    return {"id": f"fragment_{start}_{length}", "left_cut_id": left_id, "right_cut_id": right_id,
        "top_start": start, "top_length": length, "top_end": end,
        "wraps_origin": circular and start+length > n,
        "sequence": sequence, "bottom_strand_5to3": bottom,
        "core_start_in_fragment": len(left_word), "core_length": length-len(left_word),
        "left_end": _end(left_word, "left", analysis["scheme"]),
        "right_end": _end(right_word, "right", analysis["scheme"]),
        "internal_cut_sites": internal, "retained_recognition_sites": len(retained_motifs),
        "category": "Recognition-free fragment" if not retained_motifs else "Contains recognition sites / donor backbone candidate",
        "valid": not errors, "errors": errors}
