"""Bounded DNA/FASTA/GenBank input using the shared parsed-record contract."""
from __future__ import annotations

import hashlib
import io
import re
from copy import deepcopy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .genbank_service import parse_genbank

MAX_BASES = 500_000
MAX_BYTES = 5_000_000


def input_record(raw: bytes | str, topology: str | None = None) -> dict:
    raw = raw.encode("utf-8") if isinstance(raw, str) else raw
    if len(raw) > MAX_BYTES:
        raise ValueError("Input is limited to 5 MB")
    text = raw.decode("utf-8-sig").strip()
    if text.startswith("LOCUS"):
        gb = parse_genbank(raw)
    else:
        if text.startswith(">"):
            records = list(SeqIO.parse(io.StringIO(text), "fasta"))
            if len(records) != 1:
                raise ValueError("Provide exactly one FASTA record")
            record = records[0]
        else:
            record = SeqRecord(Seq(re.sub(r"\s+", "", text).upper()), id="user_DNA",
                               name="user_DNA", description="User-provided DNA")
        sequence = str(record.seq).upper()
        if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
            raise ValueError("Enter DNA using IUPAC bases only; remove numbering and punctuation")
        record.seq = Seq(sequence)
        record.name = re.sub(r"[^A-Za-z0-9_.-]", "_", record.id)[:16] or "user_DNA"
        record.annotations["molecule_type"] = "DNA"
        if topology in ("linear", "circular"):
            record.annotations["topology"] = topology
        output = io.StringIO()
        SeqIO.write(record, output, "genbank")
        gb = parse_genbank(output.getvalue().encode("utf-8"))
    if gb["length"] > MAX_BASES:
        raise ValueError("Analysis is limited to 500,000 bases per record")
    return {"genbank": gb, "input_sha256": hashlib.sha256(raw).hexdigest(), "source": "User input"}


def with_topology(gb: dict, topology: str) -> dict:
    """Analysis/view copy; cached original bytes and annotations stay untouched."""
    if topology not in ("linear", "circular"):
        raise ValueError("Choose linear or circular topology explicitly")
    copy = {**gb, "record": deepcopy(gb["record"]), "annotations": deepcopy(gb["annotations"])}
    copy["topology"] = topology
    copy["annotations"]["topology"] = topology
    copy["record"].annotations["topology"] = topology
    # Viewer identity includes the explicit override without claiming different GB bytes.
    copy["viewer_revision"] = f"{gb['sha256']}:{topology}"
    return copy
