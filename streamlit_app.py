#!/usr/bin/env python3
"""Open DNA Collections Interactive Database (Streamlit, cache-first architecture)."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import time
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import streamlit as st

from services import BlastService, DNACollectionDataService
from services.data_processing import normalize_id
from services.freegenes_service import FreeGenesService
from services.part_service import PartService
from ui.part_details import open_details, part_details_dialog
from ui.results_table import render_results_table
from ui.debug import show_debug_page

APP_BASE = Path(os.getenv("OPEN_DNA_BASE_PATH", str(Path(__file__).parent))).resolve()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


st.set_page_config(
    page_title="Open DNA Collections Database",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
<style>
    .main-header {
        font-size: 2.8rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 1.2rem;
    }
    .sequence-box {
        background-color: #f8f9fa;
        color: #0f172a !important;
        padding: 1rem;
        border-radius: 0.5rem;
        font-family: 'Courier New', monospace;
        font-size: 0.9rem;
        word-break: break-all;
        border: 1px solid #dee2e6;
        white-space: pre-wrap;
    }
    .sequence-box * {
        color: #0f172a !important;
    }
    .freshness-card {
        background: #f1f5f9;
        border: 1px solid #cbd5e1;
        color: #0f172a !important;
        border-radius: 8px;
        padding: 12px;
        margin-bottom: 12px;
        line-height: 1.5;
    }
    .freshness-card * {
        color: #0f172a !important;
    }
</style>
""",
    unsafe_allow_html=True,
)


@st.cache_resource
def load_data_service(revision: str) -> DNACollectionDataService:
    base_path = str(APP_BASE)
    return DNACollectionDataService(base_path)


@st.cache_resource
def load_blast_service() -> BlastService:
    base_path = str(APP_BASE)
    ncbi_email = None

    # Prefer Streamlit secrets, then environment variable.
    try:
        if "NCBI_EMAIL" in st.secrets:
            ncbi_email = st.secrets["NCBI_EMAIL"]
    except Exception:
        ncbi_email = None

    if not ncbi_email:
        ncbi_email = st.session_state.get("ncbi_email") or None

    return BlastService(base_path=base_path, ncbi_email=ncbi_email)


def manifest_revision(path: Path) -> str:
    if path.exists():
        return hashlib.sha256(path.read_bytes()).hexdigest()
    local = APP_BASE / "odc_plasmids.csv"
    return f"missing:{local.stat().st_mtime_ns if local.exists() else 0}"


@st.cache_resource(max_entries=3)
def load_freegenes_service(revision: str) -> FreeGenesService:
    return FreeGenesService(APP_BASE)


@st.cache_resource(max_entries=3)
def load_part_service(reclone_revision: str, freegenes_revision: str) -> PartService:
    return PartService(load_data_service(reclone_revision), load_freegenes_service(freegenes_revision))


def extract_subject_candidates(subject_id: object) -> list[str]:
    if subject_id is None:
        return []

    text = str(subject_id).strip()
    if not text:
        return []

    candidates: list[str] = []

    def _add(value: object) -> None:
        if value is None:
            return
        val = str(value).strip()
        if val and val not in candidates:
            candidates.append(val)

    _add(text)
    _add(text.upper())

    for token in re.split(r"[|,\s;]+", text):
        token = token.strip()
        if token:
            _add(token)
            _add(token.upper())

    for pattern in (r"ODC[_-]?\d+", r"BBF10K[_-]?\d+"):
        for match in re.finditer(pattern, text, flags=re.IGNORECASE):
            _add(match.group(0))

    for raw in list(candidates):
        norm = normalize_id(raw)
        if norm:
            _add(norm)
            _add(norm.upper())
            _add(f"lcl|{norm}")

    return candidates


def build_part_metadata_lookup(
    main_df: pd.DataFrame, genbank_index_df: pd.DataFrame
) -> Dict[str, Dict[str, str]]:
    lookup: Dict[str, Dict[str, str]] = {}
    by_norm: Dict[str, Dict[str, str]] = {}

    def add_key(key: object, meta: Dict[str, str]) -> None:
        if key is None:
            return
        value = str(key).strip()
        if not value:
            return
        if value not in lookup:
            lookup[value] = meta
        upper = value.upper()
        if upper not in lookup:
            lookup[upper] = meta

    if not main_df.empty:
        for _, row in main_df.iterrows():
            row_meta = {
                "Part Name": row.get("Name"),
                "Collection": row.get("Collection"),
                "ODC ID": row.get("ODC ID"),
                "BBF ID": row.get("BBF ID"),
            }

            odc_norm = row.get("ODC_ID_NORM")
            bbf_norm = row.get("BBF_ID_NORM")

            if pd.notna(odc_norm) and odc_norm:
                by_norm[str(odc_norm)] = row_meta
            if pd.notna(bbf_norm) and bbf_norm:
                by_norm[str(bbf_norm)] = row_meta

            add_key(row.get("ODC ID"), row_meta)
            add_key(row.get("BBF ID"), row_meta)
            add_key(odc_norm, row_meta)
            add_key(bbf_norm, row_meta)

    if not genbank_index_df.empty:
        for _, row in genbank_index_df.iterrows():
            part_norm = normalize_id(row.get("part_id"))
            record_norm = normalize_id(row.get("record_id"))
            base_meta = by_norm.get(str(part_norm)) or by_norm.get(str(record_norm))

            row_meta = base_meta or {
                "Part Name": row.get("description") or row.get("record_id") or row.get("part_id"),
                "Collection": None,
                "ODC ID": part_norm if str(part_norm).startswith("ODC_") else None,
                "BBF ID": part_norm if str(part_norm).startswith("BBF10K_") else None,
            }

            add_key(row.get("part_id"), row_meta)
            add_key(row.get("record_id"), row_meta)
            add_key(part_norm, row_meta)
            add_key(record_norm, row_meta)
            if part_norm:
                add_key(f"lcl|{part_norm}", row_meta)
            if record_norm:
                add_key(f"lcl|{record_norm}", row_meta)

    return lookup


def resolve_subject_metadata(
    subject_id: object, lookup: Dict[str, Dict[str, str]]
) -> tuple[Optional[str], Dict[str, str]]:
    for candidate in extract_subject_candidates(subject_id):
        meta = lookup.get(candidate)
        if meta:
            return candidate, meta
    return None, {}


def show_home_page(service: DNACollectionDataService) -> None:
    st.markdown("## Welcome to the Open DNA Collections Database")
    st.markdown("Search Reclone collections and FreeGenes, explore annotated DNA parts, and download sequences and metadata.")
    summary = service.get_collections_summary()
    c1, c2, c3 = st.columns(3)
    c1.metric("Reclone parts", sum(info["count"] for key, info in summary.items() if key != "GenBank Files"))
    c2.metric("Reclone collections", len([k for k in summary if k != "GenBank Files"]))
    c3.metric("Indexed local GenBank files", summary.get("GenBank Files", {}).get("count", 0))
    st.markdown("## Collections overview")
    for collection, info in summary.items():
        if collection != "GenBank Files":
            with st.expander(f"{collection} ({info['count']} parts)"):
                st.write(f"BBF identifiers: {info['with_bbf_id']} · ODC identifiers: {info['with_odc_id']}")


def show_search_page(service, freegenes, reclone_revision) -> None:
    st.markdown("## Search & Browse Collections")
    parts = load_part_service(reclone_revision, manifest_revision(freegenes.directory / "manifest.json"))
    st.caption("Search Reclone parts with matching FreeGenes information. Submit an empty query to browse the Reclone inventory. Click anywhere on a row to open details.")
    with st.form("part_search"):
        c1, c2 = st.columns([2, 1])
        query = c1.text_input("Search parts", placeholder="Name, description, BBF ID, ODC ID, or collection", key="search_query")
        collection = c2.selectbox("Collection", ["All Collections", *parts.collections], key="search_collection")
        with st.expander("Location filters"):
            c1, c2, c3, c4 = st.columns(4)
            resistance = c1.text_input("Resistance contains", key="search_resistance")
            strain = c2.text_input("Strain contains", key="search_strain")
            well = c3.text_input("Well pattern", placeholder="A1, A*, *1", key="search_well")
            provider = c4.selectbox("Location provider", ["All providers", "Reclone", "FreeGenes"], key="search_provider")
        submitted = st.form_submit_button("Search", type="primary")
    if submitted:
        with st.spinner("Checking FreeGenes sources…"):
            refresh = freegenes.ensure_fresh()
        fg_revision = manifest_revision(freegenes.directory / "manifest.json")
        parts = load_part_service(reclone_revision, fg_revision)
        filters = {k: v.strip() for k, v in {"resistance": resistance, "strain": strain, "well_pattern": well}.items() if v.strip()}
        if provider != "All providers":
            filters["provider"] = provider
        st.session_state["submitted_search"] = {"query": query, "collection": collection, "platemap_filters": filters}
        st.session_state.pop("results_page", None)
        st.session_state["details_open"] = False
        st.session_state["selection_epoch"] = st.session_state.get("selection_epoch", 0) + 1
        st.session_state.pop("search_results", None)
        st.session_state["search_refresh"] = refresh
    last = st.session_state.get("submitted_search")
    if not last:
        return
    source_revision = "reclone-inventory-v2:" + reclone_revision + parts.freegenes.revision
    if "search_results" not in st.session_state or st.session_state.get("search_revision") != source_revision:
        st.session_state["search_results"] = parts.search_parts(**last)
        st.session_state["search_revision"] = source_revision
        st.session_state["selection_epoch"] = st.session_state.get("selection_epoch", 0) + 1
    results = st.session_state["search_results"]
    refresh = st.session_state.get("search_refresh", {})
    if refresh.get("status") in ("unavailable", "offline"):
        st.warning(refresh.get("message", "FreeGenes refresh unavailable; using cached sources."))
    fg = parts.freegenes
    if fg.backend_records:
        st.info("FreeGenes metadata is pulled from its backend database and cached between searches.")
        backend = fg.payload.get("backend_status", {})
        st.caption("Backend retrieval: " + str(backend.get("retrieved_at", "previous cache")))
        if backend.get("status") != "ok":
            st.warning("FreeGenes backend refresh unavailable; using its last successful cached metadata.")
    else:
        st.info("Using the FreeGenes GitHub metadata snapshot; backend metadata unavailable.")
    if fg.manifest:
        st.caption("GenBank source snapshot: " + fg.manifest["commit_date"][:10])
    else:
        st.warning("FreeGenes index unavailable. Reclone results remain searchable.")
    st.caption("FreeGenes plate/well data is unavailable. Reclone distribution locations are labeled by provider.")
    st.markdown(f"### Found {len(results)} results")
    st.caption("Submitted query: " + (last["query"] or "All parts") + " · Collection: " + last["collection"])
    if results.empty:
        st.info("No results found. Try adjusting your search terms or location filters.")
        return
    c1, c2, c3 = st.columns(3)
    show_all = c1.checkbox("Show all columns", key="search_full_table")
    page_size = c2.selectbox("Items per page", [10, 25, 50, 100], index=1, key="search_page_size")
    pages = (len(results) - 1) // page_size + 1
    if st.session_state.get("results_page", 1) > pages:
        st.session_state["results_page"] = 1
    page_number = c3.selectbox("Page", range(1, pages + 1), key="results_page")
    page = results.iloc[(page_number - 1) * page_size:page_number * page_size].reset_index(drop=True)
    page_identity = (source_revision, json.dumps(last, sort_keys=True), page_size, page_number, show_all)
    if st.session_state.get("displayed_page_identity") != page_identity:
        st.session_state["displayed_page_identity"] = page_identity
        st.session_state["details_open"] = False
        st.session_state["selection_epoch"] = st.session_state.get("selection_epoch", 0) + 1
    # Keep the same table mounted during dialog open/close; reset only with the page.
    table_key = "part_results_" + hashlib.sha256(json.dumps(page_identity).encode()).hexdigest()[:16]
    keys = page["Part Key"].tolist()
    display = page if show_all else page[["BBF ID", "ODC ID", "Name", "Collection", "Sources", "Locations", "Source Conflicts"]]
    render_results_table(display, keys, table_key)
    c1, c2 = st.columns([3, 1])
    chosen = c1.selectbox("Part on this page", keys, format_func=lambda k: k + " — " + parts.parts[k]["name"],
                          key=f"part_choice_{page_number}_{page_size}")
    if c2.button("View details", key="view_details"):
        open_details(chosen)
    st.download_button("Download search results CSV", results.to_csv(index=False), "search_results.csv", "text/csv", on_click="ignore")
    if st.session_state.get("details_open"):
        part_details_dialog(parts)


def show_blast_page(service: DNACollectionDataService) -> None:
    st.markdown("## BLAST Search")
    st.markdown(
        """
Run sequence similarity search against local Open DNA collection sequences.
The local BLAST index is separate from the FreeGenes records shown in part details.
Optional NCBI fallback can be used when local BLAST has no hits.
"""
    )

    if "ncbi_email" not in st.session_state:
        st.session_state["ncbi_email"] = ""

    with st.expander("NCBI Configuration (for remote fallback)"):
        st.session_state["ncbi_email"] = st.text_input(
            "NCBI Email",
            value=st.session_state["ncbi_email"],
            help="Required by NCBI for remote BLAST fallback.",
        )

    blast_service = BlastService(
        base_path=str(APP_BASE),
        ncbi_email=st.session_state["ncbi_email"] or None,
    )

    c1, c2, c3 = st.columns(3)
    with c1:
        mode = st.selectbox(
            "Execution Mode",
            ["hybrid", "local", "ncbi"],
            index=0,
            help="hybrid = local first, fallback to NCBI when local has no hits",
        )
    with c2:
        target = st.selectbox("Target", ["auto", "local", "ncbi"], index=0)
    with c3:
        max_hits = st.slider("Max Hits", min_value=5, max_value=100, value=25, step=5)

    p1, p2 = st.columns(2)
    with p1:
        evalue = st.number_input("E-value threshold", min_value=0.0, value=1e-5, format="%.1e")
    with p2:
        poll_interval = st.number_input(
            "NCBI poll interval (sec)",
            min_value=60,
            max_value=300,
            value=60,
            step=30,
            help="NCBI recommends no more than 1 status poll per minute for a RID.",
        )

    query_sequence = st.text_area(
        "Query Sequence (DNA or Protein)",
        height=180,
        placeholder="Paste sequence here",
    )

    if st.button("Run BLAST"):
        if not query_sequence.strip():
            st.error("Please provide a sequence.")
            return

        params = {
            "max_hits": int(max_hits),
            "evalue": float(evalue),
            "poll_interval_sec": int(poll_interval),
            "timeout_sec": 300,
        }

        with st.spinner("Running BLAST query..."):
            start = time.perf_counter()
            result = blast_service.run_blast(
                mode=mode,
                sequence=query_sequence,
                target=target,
                params=params,
            )
            elapsed = time.perf_counter() - start

        st.markdown("### BLAST Result")
        st.caption(f"Elapsed: {elapsed:.2f} s")

        status = result.get("status", "unknown")
        source = result.get("source", "n/a")
        program = result.get("program", "n/a")
        query_type = result.get("query_type", "n/a")
        job_id = result.get("job_id", "n/a")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Status", status)
        c2.metric("Source", source)
        c3.metric("Program", program)
        c4.metric("Query Type", query_type)
        st.text(f"Job ID: {job_id}")

        if result.get("error"):
            st.error(result["error"])

        hits = result.get("hits", [])
        if hits:
            hits_df = pd.DataFrame(hits)
            part_lookup = build_part_metadata_lookup(service.main_df, service.genbank_index_df)
            resolved = hits_df["subject_id"].apply(lambda sid: resolve_subject_metadata(sid, part_lookup))

            hits_df["Matched Key"] = resolved.apply(lambda item: item[0])
            hits_df["Part Name"] = resolved.apply(lambda item: item[1].get("Part Name"))
            hits_df["Collection"] = resolved.apply(lambda item: item[1].get("Collection"))
            hits_df["ODC ID"] = resolved.apply(lambda item: item[1].get("ODC ID"))
            hits_df["BBF ID"] = resolved.apply(lambda item: item[1].get("BBF ID"))

            preferred_order = [
                "subject_id",
                "Part Name",
                "Collection",
                "Matched Key",
                "ODC ID",
                "BBF ID",
                "identity",
                "alignment_length",
                "evalue",
                "bitscore",
                "qcov",
                "sstart",
                "send",
            ]
            ordered_cols = [c for c in preferred_order if c in hits_df.columns] + [
                c for c in hits_df.columns if c not in preferred_order
            ]
            hits_df = hits_df[ordered_cols]

            st.dataframe(hits_df, use_container_width=True)
            st.download_button(
                "Download Hits (CSV)",
                hits_df.to_csv(index=False),
                f"blast_hits_{job_id}.csv",
                "text/csv",
            )
        else:
            st.info("No hits found.")

        raw_ref = result.get("raw_output_ref")
        if raw_ref:
            raw_path = Path(str(raw_ref))
            if raw_path.exists():
                st.download_button(
                    "Download Raw Output",
                    raw_path.read_text(encoding="utf-8", errors="ignore"),
                    raw_path.name,
                    "text/plain",
                )



def main() -> None:
    st.markdown('<h1 class="main-header">Open DNA Collections Database</h1>', unsafe_allow_html=True)
    st.caption("Reclone collections and FreeGenes part information")

    try:
        reclone_revision = manifest_revision(APP_BASE / "data" / "cache" / "manifest.json")
        service = load_data_service(reclone_revision)
        freegenes = load_freegenes_service(manifest_revision(APP_BASE / "data" / "freegenes" / "manifest.json"))
    except Exception as exc:
        st.error(f"Failed to load data service: {exc}")
        st.stop()

    st.sidebar.title("Navigation")
    if "current_page" not in st.session_state:
        st.session_state.current_page = "Home"

    pages = [
        "Home",
        "Search & Browse",
        "BLAST Search",
        "Debug",
    ]

    if st.session_state.current_page not in pages:
        st.session_state.current_page = "Debug" if st.session_state.current_page == "Data Management" else "Search & Browse"

    for page_name in pages:
        if st.sidebar.button(
            page_name,
            use_container_width=True,
            type="primary" if st.session_state.current_page == page_name else "secondary",
        ):
            st.session_state.current_page = page_name
            st.session_state["details_open"] = False
            st.rerun()

    page = st.session_state.current_page
    if page == "Home":
        show_home_page(service)
    elif page == "Search & Browse":
        show_search_page(service, freegenes, reclone_revision)
    elif page == "BLAST Search":
        show_blast_page(service)
    elif page == "Debug":
        show_debug_page(service, freegenes)


if __name__ == "__main__":
    main()
