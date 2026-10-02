"""Read-only details dialog backed by one resolved record."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from services.genbank_service import export_part, feature_table
from .sequence_viewer import render_sequence_viewer


def dismiss_details():
    st.session_state["details_open"] = False


def open_details(key):
    st.session_state["selected_part_key"] = key
    st.session_state["details_open"] = True
    # New table widget identity clears the selection; the same row can be reopened.
    st.session_state["selection_epoch"] = st.session_state.get("selection_epoch", 0) + 1


def render_details(details, location_status):
    st.subheader(details["name"])
    st.text("Identifiers: " + ", ".join(details["aliases"]))
    st.text("Collections: " + "; ".join(details["collections"]))
    gb = details.get("genbank")
    hidden_warnings = set((gb or {}).get("parse_warnings", []))
    for warning in details["warnings"]:
        if warning not in hidden_warnings and not warning.startswith((
                "LOCUS topology is circular but", "Some source metadata differs.",
                "Duplicate FreeGenes records disagree;")):
            st.warning(warning)
    if details["locations"]:
        st.markdown("### Physical locations")
        locations = pd.DataFrame(details["locations"])
        labels = {"provider": "Provider", "plate_name": "Plate name / map", "plate_number": "Plate number",
                  "well": "Well", "distribution": "Distribution", "version": "Version",
                  "resistance": "Resistance", "strain": "Strain", "identity_conflict": "Identity conflict"}
        columns = [c for c in labels if c in locations.columns]
        st.dataframe(locations[columns].rename(columns=labels), hide_index=True, width="stretch")
    if details["metadata"]:
        st.markdown("### Part metadata")
        rows = [{"Field": k, "Value": str(v)}
                for k, v in details["metadata"].items()]
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    if gb:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Length", f"{gb['length']} bp")
        c2.metric("GC Content", f"{gb['gc_content']:.1f}%")
        c3.metric("Features", len(gb["features"]))
        c4.metric("Topology", gb["topology"])
        st.text(gb["description"])
        render_sequence_viewer(gb, key=f"details_viewer_{details['part_key']}")
        st.markdown("### DNA sequence (5′ to 3′)")
        # Plain code/text, never interpolated untrusted HTML.
        st.code("\n".join(gb["sequence"][i:i+80] for i in range(0, gb["length"], 80)), language=None)
        st.markdown("### GenBank feature list")
        st.caption("Displayed base ranges are 1-based and inclusive. Internal locations retain 0-based, half-open spans, strand, and compound operators.")
        if gb["features"]:
            st.dataframe(feature_table(gb["features"]), hide_index=True, width="stretch")
        else:
            st.info("This GenBank record contains no annotated features.")
    else:
        st.info("No usable sequence: " + details["sequence_status"] + ". Metadata and location exports are available.")
    st.markdown("### Downloads")
    exports = export_part(details)
    labels = {"gb": "Download GenBank", "csv": "Download CSV", "fasta": "Download FASTA",
              "txt": "Download TXT", "features.csv": "Download feature CSV"}
    for column, extension in zip(st.columns(len(exports)), exports):
        column.download_button(labels[extension], exports[extension], f"{details['part_key']}.{extension}",
                               "text/csv" if extension.endswith("csv") else "text/plain",
                               on_click="ignore", key=f"part_download_{extension}")


@st.dialog("Part details", width="large", on_dismiss=dismiss_details)
def part_details_dialog(parts):
    key = st.session_state.get("selected_part_key", "")
    part = parts.parts.get(key)
    if not part:
        st.info("This result is no longer in the current index. Submit the search again.")
        return
    selected_bbf = None
    if len(part["bbf_ids"]) > 1:
        st.warning("This record has incompatible BBF aliases. Select an exact sequence identity.")
        choice = st.selectbox("Sequence identity", ["Choose an exact BBF ID", *part["bbf_ids"]], key=f"identity_{key}")
        if choice in part["bbf_ids"]:
            selected_bbf = choice
    with st.spinner("Resolving part and GenBank source…"):
        details = parts.get_part_details(key, selected_bbf)
    render_details(details, parts.freegenes.status()["location_status"])
    if st.button("Close details", key="close_part_details"):
        dismiss_details()
        st.rerun()
