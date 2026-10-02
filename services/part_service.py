"""Federated part identities, metadata conflicts, and source-specific locations."""
from __future__ import annotations

import copy
import fnmatch
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from .data_processing import normalize_id
from .freegenes_service import BBF_PATTERN, FreeGenesService
from .genbank_service import InvalidGenBank, parse_genbank


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if value is None or (not isinstance(value, (dict, list)) and pd.isna(value)):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def identifier(value):
    value = normalize_id(value)
    return value if value and re.fullmatch(r"(?:ODC_\d+|BBF10K_\d{6})", value) else None


def _unique(values):
    return list(dict.fromkeys(v for v in values if v not in (None, "")))


class PartService:
    def __init__(self, reclone, freegenes: FreeGenesService):
        self.reclone = reclone
        self.freegenes = freegenes
        self.parts = {}
        self.aliases = {}
        self._build_index()

    def _build_index(self):
        parent = {}
        def find(value):
            parent.setdefault(value, value)
            if parent[value] != value:
                parent[value] = find(parent[value])
            return parent[value]
        def union(a, b):
            parent[find(a)] = find(b)
        rows = []
        for row in self.reclone.main_df.to_dict("records"):
            row = clean(row)
            ids = _unique([identifier(row.get("BBF ID")), identifier(row.get("ODC ID"))])
            if not ids:
                ids = ["LOCAL_" + hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()[:16]]
            for value in ids:
                find(value)
                union(ids[0], value)
            rows.append((ids[0], {"source": "Reclone", "fields": row,
                        "provenance": {"source": "Reclone", "revision": self.reclone.manifest.get("source_commit", "local-checkout"),
                                       "source_path": "odc_plasmids.csv"}}))
        location_rows = [clean(r) for r in self.reclone.platemaps_df.to_dict("records")]
        main_aliases = set(parent)
        # Reclone plate records are inventory too, even when absent from the master CSV.
        # Do not use a plate row to change an identity established by the master CSV.
        for row in location_rows:
            ids = _unique([identifier(row.get("BBF ID")), identifier(row.get("ODC ID"))])
            if not ids or main_aliases.intersection(ids):
                continue
            for value in ids:
                find(value)
                union(ids[0], value)
            rows.append((ids[0], {"source": "Reclone", "fields": {
                "BBF ID": row.get("BBF ID"), "ODC ID": row.get("ODC ID"),
                "Name": row.get("Name"), "Collection": row.get("Toolkit")},
                "provenance": {"source": "Reclone", "revision": self.reclone.manifest.get("source_commit", "local-checkout"),
                               "source_path": row.get("Source_Path")}}))
        for record in [*self.freegenes.backend_records, *self.freegenes.records]:
            find(record["part_id"])
            rows.append((record["part_id"], {"source": record["provenance"]["source"],
                        "fields": record["fields"], "provenance": record["provenance"]}))
        # A GB-only record is still discoverable by exact ID even without metadata.
        for path in self.freegenes.payload.get("genbank_paths", []):
            find(Path(path).stem)
        grouped = {}
        for alias in parent:
            grouped.setdefault(find(alias), []).append(alias)
        root_keys = {}
        for root, aliases in grouped.items():
            aliases = sorted(aliases)
            bbf = [a for a in aliases if BBF_PATTERN.fullmatch(a)]
            key = bbf[0] if len(bbf) == 1 else aliases[0]
            root_keys[root] = key
            self.parts[key] = {"part_key": key, "aliases": aliases, "bbf_ids": bbf,
                               "source_records": [], "locations": [], "warnings": []}
            self.aliases.update({a: key for a in aliases})
        for alias, row in rows:
            self.parts[root_keys[find(alias)]]["source_records"].append(row)
        # FreeGenes enriches the accessible Reclone inventory; it does not expand it.
        self.parts = {key: part for key, part in self.parts.items()
                      if any(row["source"] == "Reclone" for row in part["source_records"])}
        self.aliases = {alias: key for alias, key in self.aliases.items() if key in self.parts}
        locations_by_alias = {}
        for i, row in enumerate(location_rows):
            for alias in _unique([identifier(row.get("ODC ID")), identifier(row.get("BBF ID"))]):
                locations_by_alias.setdefault(alias, []).append(i)
        for part in self.parts.values():
            records = part["source_records"]
            backend = [r for r in records if r["source"] == "FreeGenes backend"]
            github = [r for r in records if r["source"] == "FreeGenes GitHub"]
            local = [r for r in records if r["source"] == "Reclone"]
            fields = sorted({k for r in [*backend, *github] for k in r["fields"]})
            metadata, field_sources, conflicts = {}, {}, {}
            for field in fields:
                values = _unique(r["fields"].get(field) for r in [*backend, *github])
                tier = next((tier for tier in (backend, github) if any(r["fields"].get(field) for r in tier)), [])
                preferred = _unique(r["fields"].get(field) for r in tier)
                if len(preferred) == 1:
                    metadata[field] = preferred[0]
                    field_sources[field] = tier[0]["source"]
                elif preferred:
                    metadata[field] = preferred
                if len(values) > 1:
                    conflicts[field] = values
            local_names = _unique(r["fields"].get("Name") for r in local)
            all_names = _unique([*(r["fields"].get("gene_name_short") for r in [*backend, *github]), *local_names])
            preferred_name = metadata.get("gene_name_short")
            part["name"] = preferred_name if isinstance(preferred_name, str) else " / ".join(all_names) or part["part_key"]
            # Keep familiar collection names visible; upstream names remain metadata/search aliases.
            part["display_name"] = " / ".join(local_names) or part["name"]
            part["collections"] = _unique(r["fields"].get("Collection") for r in local)
            part["metadata"], part["field_sources"], part["conflicts"] = metadata, field_sources, conflicts
            if len(part["bbf_ids"]) > 1:
                part["warnings"].append("Conflicting ODC-to-BBF aliases: choose an exact BBF record for the sequence.")
            if any(len(_unique(r["fields"].get(k) for r in backend or github)) > 1 for k in fields):
                part["warnings"].append("Duplicate FreeGenes records disagree; all source rows are retained.")
            if len(local_names) > 1:
                part["warnings"].append("Reclone rows have differing names; all values are retained.")
            if conflicts:
                part["warnings"].append("Some source metadata differs. See source records for each original value.")
            # Add complete source rows; no first-match join or cross-provider field splicing.
            location_indexes = sorted({i for alias in part["aliases"] for i in locations_by_alias.get(alias, [])})
            location_names = []
            for i in location_indexes:
                row = location_rows[i]
                ids = _unique([identifier(row.get("ODC ID")), identifier(row.get("BBF ID"))])
                if not set(ids).intersection(part["aliases"]):
                    continue
                exact = all(a in part["aliases"] for a in ids)
                if exact and row.get("Name"):
                    location_names.append(row["Name"])
                part["locations"].append({"provider": "Reclone", "part_ids": ids,
                    "plate_name": row.get("Platemap_Key"), "plate_number": None,
                    "well": row.get("Well_Location"), "distribution": row.get("Toolkit"),
                    "version": row.get("Platemap_Version"), "resistance": row.get("Bacterial_Resistance"),
                    "strain": row.get("Growth_Strain"), "growth_conditions": row.get("Growth_Conditions"),
                    "source_path": row.get("Source_Path"), "identity_conflict": not exact,
                    "revision": self.reclone.manifest.get("source_commit", "local-checkout")})
                if not exact:
                    part["warnings"].append("A Reclone location contains incompatible identifiers; excluded from location filters.")
            for location in self.freegenes.payload.get("locations", []):
                if location.get("part_id") in part["bbf_ids"]:
                    part["locations"].append({**location, "provider": "FreeGenes"})
            part["warnings"] = _unique(part["warnings"])
            part["sources"] = _unique(r["source"] for r in records)
            if not part["sources"]:
                part["sources"] = ["FreeGenes GitHub"]
            search_fields = ("Name", "Description", "description", "gene_name_short", "gene_name_long", "product")
            part["search_text"] = " ".join([*part["aliases"], *all_names, *part["collections"], *location_names,
                *(str(r["fields"].get(field) or "") for r in records for field in search_fields)]).casefold()

    @property
    def collections(self):
        return sorted({c for p in self.parts.values() for c in p["collections"]})

    def search_parts(self, query="", collection=None, platemap_filters=None) -> pd.DataFrame:
        query = query.strip().casefold()
        normalized_query = (normalize_id(query) or "").casefold()
        rows = []
        for p in self.parts.values():
            if query and query not in p["search_text"] and normalized_query not in p["search_text"]:
                continue
            if collection and collection != "All Collections" and collection not in p["collections"]:
                continue
            locations = [loc for loc in p["locations"] if not loc.get("identity_conflict")]
            # A single physical location must satisfy all chosen filters.
            def matches(loc):
                f = platemap_filters or {}
                return (not f.get("resistance") or f["resistance"].casefold() in str(loc.get("resistance") or "").casefold()) and (
                    not f.get("strain") or f["strain"].casefold() in str(loc.get("strain") or "").casefold()) and (
                    not f.get("well_pattern") or fnmatch.fnmatchcase(str(loc.get("well") or "").upper(), f["well_pattern"].upper())) and (
                    not f.get("provider") or loc.get("provider") == f["provider"])
            if platemap_filters and not any(matches(loc) for loc in locations):
                continue
            rows.append({"Part Key": p["part_key"], "BBF ID": "; ".join(p["bbf_ids"]),
                         "ODC ID": "; ".join(a for a in p["aliases"] if a.startswith("ODC_")),
                         "Name": p["display_name"], "Collection": "; ".join(p["collections"]),
                         "Sources": "; ".join(p["sources"]), "Locations": len(locations),
                         "Well_Location": "; ".join(_unique(f"{loc['provider']}: {loc.get('plate_name') or 'unknown plate'} / {loc.get('well') or 'unknown well'}" for loc in locations)),
                         "Bacterial_Resistance": "; ".join(_unique(loc.get("resistance") for loc in locations)),
                         "Growth_Strain": "; ".join(_unique(loc.get("strain") for loc in locations)),
                         "Source Conflicts": len(p["conflicts"]) + (len(p["bbf_ids"]) > 1)})
        columns = ["Part Key", "BBF ID", "ODC ID", "Name", "Collection", "Sources", "Locations",
                   "Well_Location", "Bacterial_Resistance", "Growth_Strain", "Source Conflicts"]
        return pd.DataFrame(rows, columns=columns).sort_values("Part Key").reset_index(drop=True)

    def _local_genbank(self, part, selected_bbf=None):
        candidates = []
        allowed = [selected_bbf] if selected_bbf else part["aliases"]
        index = self.reclone.genbank_index_df
        if not index.empty:
            for row in index.to_dict("records"):
                if normalize_id(row.get("part_id")) in allowed:
                    path = (self.reclone.base_path / str(row.get("file_path", ""))).resolve()
                    if path.exists() and path.is_relative_to(self.reclone.base_path.resolve()):
                        candidates.append(path)
        # Manifest-free and incomplete-index fallback; discover only this exact identity.
        for alias in allowed:
            if not identifier(alias):
                continue
            for pattern in (f"genbank/{alias}.gb", f"*/Plasmids_Genbank/{alias}.gb", f"*/genbank_seq/{alias}.gb"):
                candidates.extend(self.reclone.base_path.glob(pattern))
        parsed, failures = {}, []
        for path in dict.fromkeys(candidates):
            try:
                gb = parse_genbank(path.read_bytes())
                parsed.setdefault(gb["sha256"], gb)
            except (InvalidGenBank, OSError):
                failures.append("A local GenBank candidate could not be parsed.")
        if len(parsed) > 1:
            return None, ["Local GenBank candidates disagree; no file selected automatically."]
        return next(iter(parsed.values()), None), failures

    def get_part_details(self, key: str, selected_bbf: str | None = None) -> dict:
        key = self.aliases.get(normalize_id(key), key)
        if key not in self.parts:
            return {}
        details = copy.deepcopy(self.parts[key])
        details.pop("search_text", None)
        details["provenance"] = {"metadata_fields": details["field_sources"],
                                 "freegenes_index": self.freegenes.revision,
                                 "reclone_revision": self.reclone.manifest.get("source_commit", "local-checkout")}
        details["genbank"] = None
        if len(details["bbf_ids"]) > 1 and selected_bbf not in details["bbf_ids"]:
            details["sequence_status"] = "ambiguous identity"
            return details
        bbf = selected_bbf or next(iter(details["bbf_ids"]), None)
        upstream = self.freegenes.fetch_genbank(bbf) if bbf else {"status": "absent", "message": "No BBF cross-reference"}
        details["sequence_status"] = upstream["status"]
        if upstream.get("raw"):
            details["genbank"] = parse_genbank(upstream["raw"])
            details["provenance"]["genbank"] = upstream["provenance"]
            if upstream["status"] == "stale cached":
                details["warnings"].append("Using a previous FreeGenes GenBank revision because the current request is unavailable.")
        elif upstream["status"] in ("invalid", "withdrawn"):
            details["warnings"].append(upstream.get("message", "Invalid upstream GenBank") + "; local DNA was not substituted.")
        else:
            gb, warnings = self._local_genbank(details, selected_bbf)
            details["warnings"].extend(warnings)
            if gb:
                details["genbank"] = gb
                details["sequence_status"] = "local fallback"
                details["provenance"]["genbank"] = {"source": "Reclone local checkout", "sha256": gb["sha256"],
                    "cache_status": "local fallback", "upstream_status": upstream["status"],
                    "revision": "checkout file hash", "reason": upstream.get("message")}
                details["warnings"].append(f"Using a local GenBank fallback: FreeGenes status is {upstream['status']}.")
            else:
                details["warnings"].append(upstream.get("message", "No usable GenBank record found"))
        gb = details["genbank"]
        if gb:
            details["warnings"].extend(gb["parse_warnings"])
            if gb["topology"] == "circular" and "linear" in gb["description"].casefold():
                details["warnings"].append("LOCUS topology is circular but the description mentions linear DNA.")
        details["warnings"] = _unique(details["warnings"])
        return details
