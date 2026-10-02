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

- Search & Browse lists identities from Reclone's master CSV and identified platemap records, including those enriched by matching FreeGenes records. Exclude FreeGenes-only identities and collection filters. Plate rows must not change identities established by the master CSV.
- Use Reclone collection names as visible labels in results, details and builder choices; keep FreeGenes names as searchable aliases and source metadata. Search every source's public names/descriptions and exact matching platemap names, not just the preferred metadata value. Preserve index-version invalidation for saved results and cached services.
- Pass `PART_INDEX_VERSION` explicitly to the cached PartService loader and retain its schema validator; a comment/docstring change alone is not an upgrade contract. Use the displayed page's labels for its dropdown and `part_display_name` for record labels. Never require `display_name` in legacy cached records.
- The bundled component must retry readiness until the first render acknowledgement; a one-shot ready message can be dropped before host registration. Keep the frontend regression covering this race.
- Use the bundled whole-row clickable table in `ui/results_table/`; do not revert to checkbox-only selection. Preserve stable IDs through sorting/paging and support keyboard activation.
- Keep Part Details free of metadata retrieval timestamps, parser/topology notices, missing-location notices, and source/provenance sections. Preserve diagnostic data and source validation internally and in exports; keep actionable unavailable/ambiguous sequence states visible.
- Part Details places downloads beside the name/description and the interactive viewer below the summary/metrics, before metadata/location tables. Omit the separate raw DNA code block; preserve full DNA in the viewer and exports.

## Data correctness

- Match sources through normalized identifiers and explicit ODC-to-BBF aliases. Names are search terms, not proof of identity.
- Prefer verified FreeGenes records for FreeGenes GenBank and FreeGenes plate/well information. Retain Reclone identity and collection associations.
- If FreeGenes cannot provide a current valid GenBank (missing, unavailable, invalid or withdrawn), use a valid, unambiguous local Open DNA file matched by exact BBF/ODC identity. Try local files before a stale FreeGenes cache; preserve the fallback label, original upstream status, relative source paths and hash. Current verified FreeGenes files still take precedence.
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

## Goal 2 contracts to preserve

- Viewer assets are pinned in `ui/sequence_viewer/vendor/` with a hash/license manifest. Never print/minify/rebuild these bundles during routine exploration. Preserve the readiness retry, positive-width guard, two-frame mount delay, and responsive observer. Register the component at module import rather than in a render function.
- `restriction_service.py` uses reference-strand cut boundaries. A generated top strand includes its left overhang and excludes the right complementary-strand overhang. Physical end words, bottom strand, cuts and provenance accompany it; Goal 3 must account for each junction exactly once.
- `fragment_service.py` revalidates selected boundaries, projects features in biological part order, and reparses GB before export. Partial CDS translations must not survive clipping; reading-frame adjustment requires exact source positions/strand. GB COMMENT provenance is base64 UTF-8 JSON with a labeled prefix and wrapped lines.
- Keep explicit candidate selection and stale-output guards for inputs and source revisions. Invalid fragments remain downloadable reports and cannot produce validated-part exports. Never attach original physical inventory locations to a generated virtual fragment.
- Goal 2 checks and limits are recorded in `docs/GOAL2_VALIDATION.md`. Full multi-part assembly remains Goal 3.

## Validation and token discipline

- Use focused file reads and small fixtures; do not dump thousands of gene HTML pages, full sequence CSVs, or upstream dependency trees into context.
- Complete one staged goal at a time. Record compact progress, exact checks, and unresolved source limitations in the plan before a handoff.
- Expected implementation checks: `python -m pytest tests/unit tests/integration`, `python test_app.py`, and a local Streamlit browser smoke check when UI changes are made.
- Mock external network calls in automated tests. Keep a bounded manual public-source check separate.
- Test source conflicts/fallbacks, stable row selection, download consistency, and circular/reverse-strand sequence cases. Browser-check dialog dismissal/reopening and custom viewer behavior.
- Dependencies are declared in `requirements.txt`; pytest is in `requirements-dev.txt`. Use the isolated `.venv` for checks. Streamlit 1.49.1 is pinned to the verified dialog/selection APIs.
- Report missing dependencies or checks not executed accurately. Goal 1 validation is recorded in `docs/GOAL1_VALIDATION.md`; later goals still require their own runtime validation.
