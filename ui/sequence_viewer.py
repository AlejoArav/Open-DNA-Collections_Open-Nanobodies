"""Read-only, locally bundled TeselaGen viewer sharing the resolved GenBank."""
from pathlib import Path

import streamlit as st
from streamlit.components.v1 import declare_component

from services.viewer_service import viewer_data

_viewer = declare_component("sequence_viewer", path=str(Path(__file__).parent / "sequence_viewer"))


def render_sequence_viewer(gb, key="sequence_viewer"):
    st.markdown("### Interactive sequence viewer")
    data = viewer_data(gb)
    if gb["topology"] not in ("linear", "circular"):
        st.caption("Topology is unknown; the viewer uses a linear display without assuming circular DNA.")
    if data["unmapped_features"]:
        st.caption(f"{len(data['unmapped_features'])} features have unmappable locations; see the full feature list.")
    st.caption("Select a feature to highlight its sequence. Use the map controls to zoom/rotate; scroll the Sequence Map to explore bases.")
    if any(len(f["spans"]) > 1 for f in gb["features"]):
        st.caption("Joined feature segments are drawn separately. Viewer selection/size spans their endpoints and can include gaps; the full feature list preserves exact locations.")
    if any(f["strand"] not in (-1, 1) for f in gb["features"]):
        st.caption("Features with unknown strand have no directional arrow; mixed strands are shown as separate segments.")
    _viewer(payload=data, default=None, key=key)
