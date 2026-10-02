"""Analyze / Generate Part workflow; assembly belongs to the next goal."""
from __future__ import annotations

import hashlib
import json

import pandas as pd
import streamlit as st

from services.builder_input_service import input_record, with_topology
from services.fragment_service import analysis_report, generate_fragment
from services.genbank_service import feature_table
from services.nomenclature_service import load_scheme
from services.restriction_service import analyze_restriction, candidate_between
from .sequence_viewer import render_sequence_viewer


def _label(end):
    junction = end["junction"]
    return junction["label"] or junction["status"]


def _candidate_label(candidate):
    status = "valid" if candidate["valid"] else "unresolved"
    return (f"{candidate['id']} · {candidate['top_length']} bp · "
            f"{_label(candidate['left_end'])} → {_label(candidate['right_end'])} · "
            f"{candidate['category']} · {status}")


def _show_features(gb):
    if gb["features"]:
        st.dataframe(feature_table(gb["features"]), hide_index=True, width="stretch")
    else:
        st.caption("This sequence has no annotated GenBank features.")


def show_builder_page(parts, revision=""):
    st.header("Interactive Builder")
    st.subheader("Analyze / Generate Part")
    st.write("Analyze a source sequence, choose a digestion fragment, and export its DNA and retained annotations.")
    mode = st.radio("Sequence source", ["Database part", "Paste DNA / FASTA", "Upload GenBank / FASTA / DNA"],
                    horizontal=True, key="builder_mode")
    raw, part_key, selected_bbf = b"", None, None
    if mode == "Database part":
        part_key = st.selectbox("Reclone part", [None, *sorted(parts.parts)],
            format_func=lambda key: "Choose a part" if key is None else f"{key} · {parts.parts[key]['name']}",
            key="builder_part")
        if part_key and len(parts.parts[part_key]["bbf_ids"]) > 1:
            st.warning("Multiple sequence identities exist. Choose an exact BBF ID.")
            selected_bbf = st.selectbox("Exact BBF identity", [None, *parts.parts[part_key]["bbf_ids"]],
                                       key=f"builder_bbf_{part_key}")
    elif mode == "Paste DNA / FASTA":
        raw = st.text_area("DNA or single-record FASTA", height=150, key="builder_dna").encode("utf-8")
    else:
        uploaded = st.file_uploader("Sequence file", type=["gb", "gbk", "genbank", "fa", "fasta", "fna", "txt"],
                                    key="builder_upload")
        if uploaded:
            raw = uploaded.getvalue()
    c1, c2, c3 = st.columns(3)
    topology = c1.selectbox("Analysis topology", [None, "linear", "circular"],
                           format_func=lambda value: "Choose topology" if value is None else value.title(),
                           key="builder_topology")
    enzyme = c2.selectbox("Enzyme", ["BsaI", "SapI"], key="builder_enzyme")
    scheme_name = c3.selectbox("Nomenclature", ["Workbook four-base scheme", "Physical ends only"],
                              key="builder_scheme")
    st.caption("Choose topology explicitly. The source GenBank is preserved; an override applies only to this analysis. Limits: 500,000 bases, 5 MB.")
    signature = hashlib.sha256(json.dumps([revision if mode == "Database part" else "", mode, part_key, selected_bbf, hashlib.sha256(raw).hexdigest(),
        topology, enzyme, scheme_name], ensure_ascii=True).encode()).hexdigest()
    if st.button("Analyze sequence", type="primary", key="builder_analyze"):
        st.session_state.pop("builder_analysis", None)
        st.session_state.pop("builder_generated", None)
        try:
            if topology is None:
                raise ValueError("Choose linear or circular topology explicitly")
            if mode == "Database part":
                if not part_key:
                    raise ValueError("Choose a Reclone part")
                with st.spinner("Resolving source sequence…"):
                    details = parts.get_part_details(part_key, selected_bbf)
                gb = details.get("genbank")
                if not gb:
                    raise ValueError("No usable sequence: " + details["sequence_status"])
                source = {"part_key": part_key, "aliases": details["aliases"], "provenance": details["provenance"]}
            else:
                value = input_record(raw)
                gb = value["genbank"]
                source = {"type": "User input", "input_sha256": value["input_sha256"]}
            scheme = load_scheme() if scheme_name == "Workbook four-base scheme" else None
            analysis = analyze_restriction(gb["sequence"], enzyme, topology, scheme)
            st.session_state["builder_analysis"] = {"signature": signature, "genbank": gb,
                "source": source, "analysis": analysis}
        except (ValueError, UnicodeError) as exc:
            st.error(str(exc))
    saved = st.session_state.get("builder_analysis")
    if not saved or saved["signature"] != signature:
        if saved:
            st.info("Inputs changed. Analyze again before selecting or downloading a fragment.")
        return
    analysis, gb, source = saved["analysis"], saved["genbank"], saved["source"]
    st.caption(f"Source topology: {gb['topology']}; analyzed as {topology}. Source: {gb['record_id']}.")
    with st.expander("Source sequence and annotations"):
        render_sequence_viewer(with_topology(gb, topology), key=f"builder_source_{signature}")
        _show_features(gb)
    st.markdown("### Restriction analysis")
    a, b, c = st.columns(3)
    a.metric("Source length", f"{analysis['length']} bp")
    b.metric("Recognition sites", len(analysis["sites"]))
    c.metric("Complete cut boundaries", len(analysis["cuts"]))
    st.caption(f"{enzyme}: {analysis['enzyme_spec']['motif']}({analysis['enzyme_spec']['top_offset']}/{analysis['enzyme_spec']['bottom_offset']}); "
               f"{analysis['overhang_length']}-base 5′ overhangs. [Cleavage reference]({analysis['enzyme_spec']['reference']}).")
    st.caption("Recognition bases use 1-based coordinates. Cut boundaries are cuts after base X; circular boundary 0 lies between the last base and base 1. All DNA words read 5′ to 3′.")
    rows = [{"Site": site["id"], "Recognition first base": site["recognition_start"] + 1,
        "Recognition last base": (site["recognition_start"] + site["recognition_length"] - 1) % analysis["length"] + 1,
        "Across origin": site["recognition_wraps"], "Strand": "+" if site["strand"] == 1 else "-",
        "Top cut after base": site["top_cut"], "Bottom cut after base": site["bottom_cut"],
        "Reference junction": site["reference_overhang_5to3"],
        "Right fragment physical 5′ end": site["right_fragment_physical_5prime"],
        "Left fragment physical 5′ end": site["left_fragment_physical_5prime"],
        "Scheme label": site["junction"]["label"] or site["junction"]["status"], "Status": site["status"]}
        for site in analysis["sites"]]
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    else:
        st.info("No recognition sites found for this enzyme. A digestion fragment cannot be generated.")
    for warning in analysis["warnings"]:
        st.warning(warning)
    st.download_button("Download analysis JSON", analysis_report(analysis, source), "restriction_analysis.json",
                       "application/json", on_click="ignore", key=f"analysis_json_{signature}")
    if rows:
        st.download_button("Download cut-site CSV", pd.DataFrame(rows).to_csv(index=False).encode("utf-8"),
                           "restriction_sites.csv", "text/csv", on_click="ignore", key=f"analysis_csv_{signature}")
    if not analysis["cuts"]:
        return
    st.markdown("### Choose a fragment")
    selection = st.radio("Boundary selection", ["Complete digestion product", "Choose cut boundaries"],
                         horizontal=True, key=f"boundary_mode_{signature}")
    candidate = None
    if selection == "Complete digestion product":
        index = st.selectbox("Digestion fragment", [None, *range(len(analysis["products"]))],
            format_func=lambda i: "Choose a fragment" if i is None else _candidate_label(analysis["products"][i]),
            key=f"fragment_choice_{signature}")
        if index is not None:
            candidate = analysis["products"][index]
    else:
        cuts = {cut["id"]: cut for cut in analysis["cuts"]}
        options = [None, *cuts]
        if topology == "linear":
            options.append("linear_end")
        def boundary_label(value):
            if value is None:
                return "Choose a boundary"
            if value == "linear_end":
                return "Original linear end"
            return f"{value}: reference cut after base {cuts[value]['top_cut']}"
        left, right = st.columns(2)
        left_id = left.selectbox("Left boundary", options, format_func=boundary_label, key=f"left_cut_{signature}")
        right_id = right.selectbox("Right boundary", options, format_func=boundary_label, key=f"right_cut_{signature}")
        if left_id is not None and right_id is not None:
            try:
                candidate = candidate_between(analysis, None if left_id == "linear_end" else left_id,
                                               None if right_id == "linear_end" else right_id)
            except ValueError as exc:
                st.error(str(exc))
    if candidate:
        st.write(f"Selected reference strand: {candidate['top_length']} bp; double-stranded core: {candidate['core_length']} bp. {candidate['category']}.")
        st.dataframe(pd.DataFrame([{"End": name, "Reference junction (5′→3′)": end["reference_junction_5to3"],
            "Physical 5′ overhang (5′→3′)": end["physical_5prime"], "Kind": end["kind"], "Scheme label": _label(end)}
            for name, end in (("Left", candidate["left_end"]), ("Right", candidate["right_end"]))]),
            hide_index=True, width="stretch")
        st.caption("Exports represent the reference top strand: the left overhang is included; the right overhang belongs to the complementary strand. End data and both strand sequences are recorded in the report. Recognition-site retention alone does not establish insert/backbone identity.")
        for error in candidate["errors"]:
            st.error(error)
        if candidate["internal_cut_sites"]:
            st.write("Internal cut sites: " + ", ".join(candidate["internal_cut_sites"]))
        st.download_button("Download selected fragment report", analysis_report(analysis, source, candidate),
            "selected_fragment_report.json", "application/json", on_click="ignore",
            key=f"selected_report_{signature}_{candidate['id']}")
    if st.button("Generate selected part", disabled=not candidate or not candidate["valid"],
                 type="primary", key=f"generate_{signature}"):
        try:
            generated = generate_fragment(gb, analysis, candidate, source)
            st.session_state["builder_generated"] = {"signature": signature, "candidate_id": candidate["id"],
                                                    "value": generated}
        except ValueError as exc:
            st.error(str(exc))
    generated = st.session_state.get("builder_generated")
    if not candidate or not generated or generated["signature"] != signature or generated["candidate_id"] != candidate["id"]:
        return
    value = generated["value"]
    output = value["details"]["genbank"]
    st.success(f"Generated {output['record_id']}: {output['length']} bp, linear digestion fragment.")
    for warning in value["details"]["warnings"]:
        st.warning(warning)
    render_sequence_viewer(output, key=f"generated_viewer_{signature}_{candidate['id']}")
    st.code(output["sequence"], language=None)
    _show_features(output)
    labels = {"gb": "GenBank", "fasta": "FASTA", "csv": "CSV", "txt": "TXT", "features.csv": "feature CSV"}
    for extension, data in value["exports"].items():
        st.download_button("Download generated " + labels[extension], data, f"{output['record_id']}.{extension}",
            "text/csv" if extension.endswith("csv") else "text/plain", on_click="ignore",
            key=f"generated_download_{extension}_{signature}_{candidate['id']}")
