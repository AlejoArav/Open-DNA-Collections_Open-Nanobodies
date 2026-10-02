"""Adapt parsed half-open GenBank spans to OVE's inclusive annotation ranges."""
from __future__ import annotations

COLORS = {"CDS": "#2e86ab", "promoter": "#e79a32", "terminator": "#9967bc",
          "rep_origin": "#45a372", "misc_feature": "#7188a2", "source": "#b6b6b6"}


def viewer_data(gb: dict) -> dict:
    length = gb["length"]
    features, unmapped = [], []
    for feature in gb["features"]:
        spans = feature["spans"]
        if not spans or any(s.get("remote_ref") or s["start"] is None or s["end"] is None
                            or not 0 <= s["start"] < s["end"] <= length for s in spans):
            unmapped.append(feature["index"])
            continue
        ordered = sorted(spans, key=lambda s: s["start"])
        # An origin-crossing join begins with its terminal span, then wraps to zero.
        if gb["topology"] == "circular" and len(ordered) > 1 and ordered[0]["start"] == 0 and ordered[-1]["end"] == length:
            ordered = [ordered[-1], *ordered[:-1]]
        locations = [{"start": s["start"], "end": s["end"] - 1} for s in ordered]
        mixed = len({s["strand"] for s in spans}) > 1
        if mixed:
            # Keep biological identity while drawing each span's actual direction.
            for i, (span, location) in enumerate(zip(ordered, locations)):
                features.append({"id": f"feature-{feature['index']}-span-{i}",
                    "name": f"{feature['label']} (segment {i+1})" + (" (strand unknown)" if span["strand"] not in (-1, 1) else ""), "type": feature["type"],
                    **location, "forward": span["strand"] != -1,
                    "strand": span["strand"],
                    **({"arrowheadType": "NONE"} if span["strand"] not in (-1, 1) else {}),
                    "color": COLORS.get(feature["type"], "#7188a2")})
        else:
            features.append({"id": f"feature-{feature['index']}", "name": feature["label"] +
                    (" (strand unknown)" if feature["strand"] not in (-1, 1) else ""),
                "type": feature["type"], "start": locations[0]["start"], "end": locations[-1]["end"],
                "locations": locations if len(locations) > 1 else [],
                "forward": feature["strand"] != -1,
                "strand": feature["strand"],
                **({"arrowheadType": "NONE"} if feature["strand"] not in (-1, 1) else {}),
                "color": COLORS.get(feature["type"], "#7188a2")})
    return {"sequenceData": {"name": gb["record_id"], "sequence": gb["sequence"],
                "circular": gb["topology"] == "circular", "features": features},
            "content_hash": gb.get("viewer_revision", gb["sha256"]), "topology": gb["topology"],
            "unmapped_features": unmapped}
