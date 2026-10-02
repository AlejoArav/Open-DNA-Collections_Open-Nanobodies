# Reclone Open DNA Collections

> A global collaboration for equitable access to biotechnology

The [Reclone Reagent Collaboration Network](https://reclone.org) has worked for the creation and disponibilization of open DNA collections. While these collections were originally distributed by [Freegenes](https://stanford.freegenes.org)
), now some parts of this project are separating into their own organization structure. The parts in this repository will receive a new ID, and this is where we will centralize the issues board for any changes proposed for the DNA or documentation on these collections.

Here you can find the genbank files and metadata for the plasmids in different collections:
- Open Enzymes Collection
- Open Reporters Collection
- Open Yeast Collection
- Open Plasmids
- E.coli Protein Expression Toolkit
- Research in Diagnostics Toolkit (formerly Molecular Diagnostics Toolkit, and Research in Diagnostics Collection; branch for RiDT is retained as MDT for now, but description is referred to as RiDT).

## Interactive Database

🧬 **[Launch the Interactive Database](https://opendnacollections.streamlit.app/)** (deployed on Streamlit Community Cloud)

The local implementation provides:

- **Home**: collection overview.
- **Search & Browse**: Reclone inventory search enriched with matching FreeGenes records. FreeGenes-only parts are excluded. Submit an empty query to browse Reclone parts. Click any cell in a row to open details; column headers sort the displayed page. Keyboard users can focus a row and press Enter or Space.
- **Part details dialog**: downloads beside the name/description, followed by sequence metrics, the interactive viewer, source-specific physical locations, metadata and complete GenBank features. The separate raw DNA block is omitted; DNA remains available in the viewer and GenBank/CSV/FASTA/TXT downloads. Retrieval/parser notices and source/provenance sections are hidden in this dialog; provenance remains in exports. A separate feature CSV is also available.
- **BLAST Search**: the existing local sequence index and optional NCBI fallback.
- **Interactive Builder**: Analyze / Generate Part from a Reclone record, pasted DNA/FASTA, or an uploaded GenBank/FASTA/DNA file. Analyze BsaI/SapI cuts, select a digestion fragment explicitly, and export its sequence, physical ends, and retained features.
- **Debug**: freshness, diagnostics, manifests, FreeGenes refresh, dataset exports, and Reclone platemaps.

Analytics and the standalone Part Details/Data Management routes were removed. The deployed URL above may still run the previous version; these changes have not been deployed. Full multi-part plasmid assembly remains Goal 3 in [the implementation plan](docs/IMPLEMENTATION_PLAN.md).

### Viewing and generating parts

Open a search result to explore its annotated circular or linear map. The viewer always uses a white background and dark text, regardless of the app theme. Click a feature to select its sequence, use the plus/minus controls to zoom, rotate circular maps with the arrow controls, and scroll the Sequence Map. Assets for MIT-licensed TeselaGen OVE 0.8.42 are bundled locally with a version/hash manifest. The viewer uses the same resolved GenBank as downloads and preserves source DNA. Unknown topology uses a linear display. Exact source feature locations remain in the full feature table: OVE selects joined features by their bounding endpoints, including gaps, and its size column describes that selection span.

In **Interactive Builder**, choose an input and select topology explicitly, then choose an enzyme and nomenclature before clicking **Analyze sequence**. A topology override applies to this analysis/view only. Select a complete digestion product or choose cut boundaries yourself. The app lists recognition orientation, cuts on both strands, end sequences, mapped/unmapped labels, and fragment validity. **Generate selected part** stays disabled until a usable fragment is selected; ambiguous DNA, additional internal cuts, and overlapping cuts block validated exports. Changing input or analysis settings hides stale results until reanalysis.

GenBank/FASTA/CSV/TXT exports describe the selected reference strand in 5′→3′ orientation: its left sticky overhang is included, while its right sticky overhang belongs to the complementary strand. Reports preserve both physical end words, complementary-strand sequence, cut boundaries, source hashes, and scheme provenance. Generated GenBank records are linear and reparsed before export. Complete retained features preserve biological part order; clipped features receive a partial-feature note, and clipped CDS translations are removed. Original records and physical inventory locations are never rewritten or assigned to virtual fragments.

Analysis JSON, cut-site CSV, and a selected-fragment JSON report are downloadable even when the selected candidate cannot be generated. Limits are 5 MB input, 500,000 bases, and 10,000 recognition matches per orientation. Input accepts IUPAC ambiguity for inspection but does not silently resolve it or mutate DNA. The workbook provides four-base junction labels with cell/hash provenance; SapI three-base labels and additional uLoop level rules remain unmapped. No assembly, domestication, primer design, or ligation-efficiency prediction is performed. See [the nomenclature reference](docs/NOMENCLATURE_REFERENCE.md) and [Goal 2 validation](docs/GOAL2_VALIDATION.md).

### Local development

Tested with Python 3.12 and Streamlit 1.49.1. On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-dev.txt
.\.venv\Scripts\python.exe -m streamlit run streamlit_app.py
```

The checked-in caches support cold starts without cloning FreeGenes. If the Reclone cache is absent, its service falls back to local CSVs. To rebuild the local cache:

```powershell
.\.venv\Scripts\python.exe scripts/sync_upstream_data.py --base-dir . --source local
```

Validation:

```powershell
.\.venv\Scripts\python.exe -m pytest tests/unit tests/integration -q
.\.venv\Scripts\python.exe -X utf8 test_app.py
```

The frontend handshake regression uses Node.js when available; its test is skipped without Node. Browser verification remains required for UI changes.

Local BLAST requires the external `blastn`, `blastp`, and `makeblastdb` executables on PATH. Python packages alone do not install those binaries. NCBI requests are optional and require the existing email configuration.

### Data sources and refresh

Reclone metadata and physical distribution locations retain their original aliases, collection associations, and provider labels. FreeGenes descriptive metadata comes from the public **Genes** worksheet in [its backend spreadsheet](https://docs.google.com/spreadsheets/d/1LZCXzBtgey9xv5OH7YGYgp8UMJ27Eyj1aF9IhAW6M6o/edit#gid=954552604), with cached retrieval timestamps and hashes. A commit-pinned [FreeGenes GitHub snapshot](https://github.com/freegenes/freegenes.github.io) supplies fallback metadata and GenBank files. The packaged GitHub revision dates to **2023-09-15**; refreshing checks the available upstream revision and does not make that snapshot newer.

Search enriches matching Reclone identities from both indexes on every submitted query and never expands the accessible Reclone inventory. Metadata refresh has a default 24-hour TTL; **Debug → Refresh FreeGenes metadata** forces a check. Individual GenBank files are fetched only when details are opened and cached by commit/part ID. Downloads preserve the exact resolved GenBank bytes; sequence, features, and other exports derive from those bytes. Invalid upstream files or conflicting backend file links are surfaced instead of silently selecting local DNA. Unavailable upstream files may use a visibly labeled previous cache or local fallback.

**FreeGenes plate/well data remains unavailable.** Genes has no plate/well columns, and the public Packaging/Collections location values are unspecified. The app displays Reclone distribution locations separately and does not invent FreeGenes plates or infer wells from collection order. An authoritative location feed is still needed.

Metadata snapshots are content-addressed, checksummed, and published manifest-last. Reclone artifacts live under `data/cache/`; FreeGenes metadata uses a compressed index under `data/freegenes/`. Disposable GenBank bytes and refresh status live in ignored `.runtime/freegenes/`. Original collection CSVs/GenBank files and the local BLAST source sequences are preserved. Community Cloud runtime storage is ephemeral, so the small metadata index is packaged in the checkout.

Configuration (environment variables):

| Variable | Purpose |
| --- | --- |
| `FREEGENES_TTL_SECONDS` | Metadata freshness interval; default `86400` |
| `FREEGENES_OFFLINE=1` | Disable FreeGenes network calls and use available caches/fallbacks |
| `GITHUB_TOKEN` | Optional GitHub API authentication; never sent to other sources |
| `OPEN_DNA_BASE_PATH` | Optional alternate data root, used by isolated tests |

Transport uses explicit timeouts, response limits, at most one transient transport retry, and no rate-limit retries. Refresh work stops starting new source requests after its 120-second deadline; requests already in flight finish under their own bounded transport limits. Failed GitHub batches retain the previous snapshot; an unavailable Genes backend retains its last good records with a visible status.

Refresh commands:

```powershell
.\.venv\Scripts\python.exe scripts/sync_upstream_data.py --base-dir . --source upstream --repo Reclone-org/Open-DNA-Collections --branch main
.\.venv\Scripts\python.exe scripts/sync_freegenes_data.py --base-dir .
```

The scheduled workflow `.github/workflows/sync_upstream_cache.yml` refreshes the providers independently and publishes valid snapshots. Automated tests mock external network calls; bounded public-source/browser checks are recorded in [Goal 1 validation](docs/GOAL1_VALIDATION.md).
