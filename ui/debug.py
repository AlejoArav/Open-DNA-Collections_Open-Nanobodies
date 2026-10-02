"""Source status, manifests, and data exports in one debugging view."""
import pandas as pd
import streamlit as st


def show_debug_page(service, freegenes):
    st.markdown("## Debug")
    st.markdown("### Reclone data freshness")
    freshness = service.get_data_freshness()
    st.json(vars(freshness))
    st.markdown("### Diagnostics")
    st.json(service.get_diagnostics())
    st.markdown("### FreeGenes sources")
    st.info("GitHub provides versioned GenBank files and snapshot metadata. The public Genes sheet provides backend metadata. Physical FreeGenes plate/well information remains unavailable.")
    status = freegenes.status()
    st.json(status)
    if st.button("Refresh FreeGenes metadata", key="refresh_freegenes"):
        with st.spinner("Refreshing the FreeGenes index…"):
            result = freegenes.ensure_fresh(force=True)
        st.session_state["freegenes_refresh_result"] = result
        st.session_state.pop("search_results", None)
        st.session_state["details_open"] = False
        st.rerun()
    if st.session_state.get("freegenes_refresh_result"):
        result = st.session_state["freegenes_refresh_result"]
        (st.success if result["status"] == "ok" else st.warning)(result.get("message", result["status"]))
    st.caption("Refresh replaces the validated metadata snapshot and reloads the app index. It does not change original collection files or the local BLAST database.")
    st.markdown("### Export data")
    datasets = {"Reclone main": service.main_df, "Reclone platemaps": service.platemaps_df,
                "Local GenBank index": service.genbank_index_df,
                "FreeGenes GitHub metadata": pd.DataFrame([r["fields"] for r in freegenes.records]),
                "FreeGenes backend metadata": pd.DataFrame([r["fields"] for r in freegenes.backend_records])}
    choice = st.selectbox("Dataset", list(datasets), key="debug_export")
    # Raw internal file paths are appropriate only for this explicit Debug dataset export.
    st.download_button("Download dataset CSV", datasets[choice].to_csv(index=False),
                       choice.replace(" ", "_").lower() + ".csv", "text/csv", on_click="ignore")
    st.markdown("### Reclone cache manifest")
    st.json(service.manifest)
    st.markdown("### Available Reclone platemaps")
    summary = service.get_platemap_summary()
    if summary.empty:
        st.info("No Reclone platemap records.")
    else:
        st.dataframe(summary, hide_index=True, width="stretch")
        selected = st.selectbox("Platemap", summary["Platemap_Key"].tolist(), key="debug_platemap")
        data = service.get_platemap(selected)
        st.dataframe(data, hide_index=True, width="stretch")
        st.download_button("Download platemap CSV", data.to_csv(index=False), "reclone_platemap.csv", "text/csv", on_click="ignore")
