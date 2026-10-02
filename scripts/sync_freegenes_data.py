"""Refresh the small FreeGenes metadata index; GenBank files stay lazy."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from services.freegenes_service import FreeGenesService


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default=str(ROOT))
    args = parser.parse_args()
    service = FreeGenesService(args.base_dir)
    result = service.ensure_fresh(force=True)
    print(json.dumps({"refresh": result, "counts": service.manifest.get("counts", {}),
                      "backend": service.payload.get("backend_status"),
                      "locations": service.status()["location_status"]}, indent=2))
    return 0 if result["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
