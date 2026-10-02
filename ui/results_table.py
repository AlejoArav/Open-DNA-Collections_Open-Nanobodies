"""Locally bundled table with whole-row click and keyboard activation."""
from pathlib import Path

import streamlit as st
from streamlit.components.v1 import declare_component

from .part_details import open_details

_table = declare_component("part_results_table", path=str(Path(__file__).parent / "results_table"))


def render_results_table(display, part_keys, key):
    rows = [{"part_key": part_key, "cells": [str(value) for value in row]}
            for part_key, row in zip(part_keys, display.itertuples(index=False, name=None))]
    event = _table(columns=list(display.columns), rows=rows, default=None, key=key)
    if not isinstance(event, dict) or event.get("part_key") not in part_keys:
        return
    token = (key, event.get("event_id"))
    if event.get("event_id") and token != st.session_state.get("last_part_row_event"):
        st.session_state["last_part_row_event"] = token
        open_details(event["part_key"])
