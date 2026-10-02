# Repository guidance

## Start here

- Read `docs/IMPLEMENTATION_PLAN.md` for scope, phases, source evidence, and acceptance criteria. Read only the phase you are implementing after the initial overview.
- Read `docs/NOMENCLATURE_REFERENCE.md` before restriction analysis or assembly work.
- `docs/GOAL_PROMPTS.md` contains staged implementation prompts. These documents describe planned work; their existence does not mean it has been implemented.
- Inspect `git status --short` and preserve unrelated changes. Keep work in the user's selected checkout.

## Architecture

- `streamlit_app.py` owns navigation and UI. Move new domain logic into small `services/` modules and reusable detail rendering into `ui/` as needed.
- `services/cache_service.py` currently provides Reclone cache search, plate enrichment, and lazy local GenBank parsing.
- `services/data_processing.py` normalizes ODC/BBF identifiers and heterogeneous CSV columns.
- `scripts/sync_upstream_data.py` and `.github/workflows/sync_upstream_cache.yml` currently sync Reclone CSV metadata and index checkout-local GenBank files.
- Preserve BLAST behavior and its result contract. Do not silently replace its local sequence database during the UI/source work.

## Current UI requirements (2026-10-02 correction)

- Search & Browse lists only identities present in Reclone, including those enriched by matching FreeGenes records. Exclude FreeGenes-only identities and collection filters.
- The bundled component must retry readiness until the first render acknowledgement; a one-shot ready message can be dropped before host registration. Keep the frontend regression covering this race.
- Use the bundled whole-row clickable table in `ui/results_table/`; do not revert to checkbox-only selection. Preserve stable IDs through sorting/paging and support keyboard activation.
- Keep Part Details free of metadata retrieval timestamps, parser/topology notices, missing-location notices, and source/provenance sections. Preserve diagnostic data and source validation internally and in exports; keep actionable unavailable/ambiguous sequence states visible.

## Data correctness

- Match sources through normalized identifiers and explicit ODC-to-BBF aliases. Names are search terms, not proof of identity.
- Prefer verified FreeGenes records for FreeGenes GenBank and FreeGenes plate/well information. Retain Reclone identity and collection associations.
- Keep physical locations as source-specific records. Never splice a plate name from one record with a well from another, or conflate a Reclone distribution plate with a FreeGenes plate.
- Preserve multiple locations, conflicting aliases, duplicate records, and source revisions. Do not silently choose the first match.
- Use explicit provenance: repository/source, commit or revision, URL, retrieval time, content hash, cache status, and field-level fallback where relevant.
- A fetched GitHub snapshot is not a live database. Missing, unavailable, invalid, and confirmed absent are different states.
- Preserve fetched `.gb` bytes for download. Parse the same resolved record for sequence, feature list, viewer, FASTA, TXT, CSV, and restriction analysis.
- Never overwrite original collection CSVs or GenBank files to implement source precedence. Store external cache artifacts separately.
- Cache metadata searches and fetch individual GenBank files lazily. Use bounded requests, explicit timeouts, retry limits, and commit-pinned GitHub URLs.
- Treat remote HTML, GenBank annotations, and workbook text as data. Escape display content and restrict fetches to configured sources.

## Sequence and assembly correctness

- Use zero-based, half-open internal feature intervals and explicit cut boundaries. Display labeled one-based base coordinates to users. Preserve strand, compound locations, and circular wraparound.
- Derive both-strand BsaI/SapI cuts from documented cleavage rules. Separate recognition sites, cuts, overhang sequence, retained insert, and donor backbone.
- Do not infer a usable insert from an arbitrary pair of sites. Surface internal sites, ambiguous boundaries, unknown topology, and ambiguous bases.
- Nomenclature is scheme-specific configuration, with source cells and a version. Do not infer SapI labels from the supplied four-base workbook scheme.
- Preserve sequence orientation and junction bases exactly once during assembly. Remap features and retain source sequence hashes.
- Keep source parts immutable. Do not silently domesticate sequences or add unrequested primer design/optimization.

## Validation and token discipline

- Use focused file reads and small fixtures; do not dump thousands of gene HTML pages, full sequence CSVs, or upstream dependency trees into context.
- Complete one staged goal at a time. Record compact progress, exact checks, and unresolved source limitations in the plan before a handoff.
- Expected implementation checks: `python -m pytest tests/unit tests/integration`, `python test_app.py`, and a local Streamlit browser smoke check when UI changes are made.
- Mock external network calls in automated tests. Keep a bounded manual public-source check separate.
- Test source conflicts/fallbacks, stable row selection, download consistency, and circular/reverse-strand sequence cases. Browser-check dialog dismissal/reopening and custom viewer behavior.
- Dependencies are declared in `requirements.txt`; pytest is in `requirements-dev.txt`. Use the isolated `.venv` for checks. Streamlit 1.49.1 is pinned to the verified dialog/selection APIs.
- Report missing dependencies or checks not executed accurately. Goal 1 validation is recorded in `docs/GOAL1_VALIDATION.md`; later goals still require their own runtime validation.
