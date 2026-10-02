# Goal 1 validation — 2026-10-02

## Scope and result

Goal 1 is complete within the staged prompt's explicit allowance for an unavailable authoritative location feed. UI/source work is implemented; **FreeGenes physical plate/well sourcing remains unmet**. No original collection CSV/GenBank changes, deployment, viewer, restriction analysis or assembly implementation were made.

## Automated checks

Environment: isolated `.venv`, Python 3.12.10, Streamlit 1.49.1, pandas 2.3.3, BioPython 1.88, pytest 8.4.2.

```powershell
.\.venv\Scripts\python.exe -m pytest tests/unit tests/integration -q
# 50 passed in 2.87s
.\.venv\Scripts\python.exe -X utf8 test_app.py
# Smoke tests passed
 git diff --check
# Passed
```

The pre-change baseline had 10 tests. Added offline fixtures cover metadata precedence, exact aliases and conflicts, duplicate records, empty Reclone/upstream-only search, multiple source-specific locations, combined location filters, upstream byte priority, missing/invalid/withdrawn files, previous revision caches, local disagreements, rate limits, timeout/negative caches, malformed source responses, interrupted manifest publication, compression/legacy loading, source-pinned requests, and dialog state. External calls are mocked in automated tests.

The smoke script validates Reclone data loading (337 main rows, 770 platemap rows), `ODC_0007` search and the BLAST response contract. Its BLAST response is **failed** because BLAST executables are missing, not an alignment success. Its emoji output needs `-X utf8` in this Windows shell; the plain command raises a console encoding error before checks run.

## Browser evidence

Local server: `http://127.0.0.1:8501`, run with `.venv` and telemetry disabled.

- Home shows the collection overview and the four retained navigation buttons; freshness/diagnostics are in Debug.
- Empty submitted search returns 2,349 merged identities from Reclone and FreeGenes. A cold server start loads the packaged compressed metadata successfully.
- `ODC_0007` opens `BBF10K_003247` / `9N7polA`. The dialog shows backend retrieval provenance and commit-pinned FreeGenes GenBank priority, all three Reclone location records, 4,451 bp, 12 full features and the circular/linear annotation disagreement.
- Page 2, Name descending: selecting the displayed `BBF10K_000046` row opens **Sc-pCYC1**, confirming sorting/page positions resolve to the correct stable identity.
- Native X and explicit Close details both dismiss. Page 2 persists; the same part can reopen. CSV/FASTA/GB/TXT downloads leave the dialog and submitted query intact.
- Debug displays merged freshness/diagnostics/manifests/exports/plate maps. Explicit metadata refresh completed successfully; source/index revision reloaded.
- BLAST page preserves its controls/local-source label and rejects an empty query. No remote NCBI request was sent and no local alignment was claimed.
- Browser console errors observed during the intentional server restart were connection/init-ping errors; the subsequent cold start and complete flows succeeded. No application exception was observed.

Preview: `C:/Users/Alejandro/.codex/visualizations/2026/10/02/01a0fa70-3e49-7b30-b05f-a2f42899b810/goal1-details.png`.

Downloaded sample GenBank: 9,130 bytes, SHA-256 `5a93497f1f4b445a37c8714a271165422bd9077fb0a617af3ae8d3fc8c80b469`. Downloaded FASTA, metadata CSV and TXT (after removing sequence wrapping) each reproduce the same 4,451-base sequence. CSV/TXT exports do not expose the server base path. GB preserves upstream bytes exactly.

## Public source findings and unresolved inputs

GitHub revision: `a477c45baac46483f03bc61ccde92867920fabe7`, dated 2023-09-15. Counts: 2,334 CSV rows, 2,346 usable public Genes backend records, 2,338 indexed GenBank files. Current compressed package: 620,214 bytes. Selective metadata/byte fetches avoid a full repository clone.

The public Genes worksheet is accessible via Google Visualization CSV with explicit `headers=1`. Its exact part IDs and public GenBank links are used; links that disagree with indexed identity block automatic DNA selection. Production/insert sequence and internal-note fields are omitted from the metadata request/cache.

Genes contains no plate/well columns. The bounded public Packaging/Collections check found unspecified locations and no verified physical plate/well identities. **Needed input:** an authoritative feed/worksheet/API with part ID, plate name/number, well and source/version semantics. Collection composition order is not treated as a well map. Reclone locations are labeled fallback distribution records.

The daily GitHub Actions workflow was updated but has not been run on GitHub in this task. External BLAST executables are still required for a real local alignment. Goal 2/3 retain the separate viewer/analysis/assembly validation obligations in the plan.


## Subsequent correction — 2026-10-02

Current behavior supersedes the original 2,349-part browse evidence above: only **337 Reclone identities** are listed. FreeGenes-only IDs are excluded even from exact-ID searches; matched upstream metadata and GenBank precedence are preserved. Automated fixtures verify Reclone-only, matched and excluded upstream-only cases, plus empty Reclone inventory.

Whole-row activation now uses `ui/results_table/index.html` and its Python adapter, bundled locally with no new dependencies. The [Streamlit custom component protocol](https://docs.streamlit.io/develop/concepts/custom-components/components-v1/intro) provides the event bridge. Names and other cells are escaped through DOM `textContent`; events contain stable keys rather than sorted positions. Header sorting and Enter/Space activation are supported.

The listed Part Details messages and provenance sections/captions are hidden; exports retain original diagnostic/source information. Updated suite: **52 passed in 3.28s**. Browser directly clicked the MERScontrol name cell after Name sorting, verified BBF10K_000004/ODC_0218, sequence/features and five downloads, confirmed none of the requested messages/provenance text appeared, dismissed using X, then reopened by clicking its ODC-ID cell. No checkbox was used. Preview: `corrections-details.png` in the existing visualization output folder.


## Component loading fix — 2026-10-02

Reproduced the reported blank table in a fresh browser session. The console showed `Received component message for unregistered ComponentInstance!`: the cached frontend emitted its one-shot ready message before Streamlit registered the iframe listener. The local asset was present; the failure was the startup handshake.

The bundled frontend now retries `streamlit:componentReady` every 100 ms until the first `streamlit:render` acknowledgement, then stops the timer. Unloading the iframe also clears it. Added a Node VM regression exercising the actual bundled script with an initially unregistered host: the initial message is dropped, a later retry recovers, rendering succeeds, and clicking the row returns its stable key.

Validation: **53 tests passed in 3.94s**, including the new frontend regression (Node is available locally; that test skips when Node is absent). Three consecutive browser reload/search cycles rendered the 337-part table with no component-loading message. Sorted row clicking, the details dialog, download controls, and native dismissal still worked. Screenshot: `table-loading-fixed.png` in the visualization output folder.
