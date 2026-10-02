"""Exercise the actual bundled JavaScript with a late-registering host."""
import shutil
import subprocess
from pathlib import Path

import pytest


def test_frontend_recovers_after_early_ready_message_is_dropped():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is required for the frontend protocol regression")
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run([node, str(root / "tests/frontend/results_table_handshake.cjs"),
                             str(root / "ui/results_table/index.html")],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
    assert "whole-row click succeeds" in result.stdout
