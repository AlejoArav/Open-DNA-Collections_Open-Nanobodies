# Streamlit implementation plan

Planning date: 2026-10-01 (America/Santiago). Status: analysis and documentation complete; application changes not started.

## Recommended product

Navigation: **Home, Search & Browse, BLAST Search, Interactive Builder, Debug**. Add Interactive Builder when its first working part-analysis mode is ready.

Clicking a search row opens a wide details dialog. Remove the standalone Part Details navigation item once that interaction works. Reuse the detail renderer inside the builder and, where useful, BLAST hits. Home keeps the collection overview. Debug combines freshness, diagnostics, cache manifests, dataset exports, and source status. Remove Analytics and its navigation branch.

Interactive Builder grows in two steps: **Analyze / Generate Part**, then **Assemble Plasmid**. The first step inspects real sequences and produces a validated part representation. The second assembles those parts into a selected backbone.

## Current implementation and concrete gaps

| Area | Current evidence | Required change |
| --- | --- | --- |
| Navigation | `streamlit_app.py:763` uses sidebar buttons and `current_page`, not actual `st.tabs` | Retain this simple routing; rename/remove pages and handle stale session page names |
| Home | `show_home_page`, line 248, calls `show_data_freshness` | Move freshness and diagnostics to Debug |
| Search | `show_search_page`, line 291; passive tables at lines 370/387 | Single-row selection, stable part identity, persistent search/filter/page state |
| Details | `show_part_details_page`, line 487; manual ID input | Shared dialog renderer with locations, full features, all downloads, later viewer |
| Features | `cache_service.py:_parse_genbank_details` retains a BioPython record but exposes only a count | Serialize/render full feature type, labels, coordinates, strand, qualifiers, and compound spans |
| Exports | Details currently offer GB and FASTA | Add resolved metadata/sequence/features CSV and readable TXT |
| Sources | Reclone metadata sync; GenBank indexing still scans local files | Separate FreeGenes provider/index and explicit resolution policy |
| Locations | `build_platemap_lookups` drops duplicate keys, keeping first | Keep source-specific one-to-many locations; remove ambiguous first-match behavior from the new resolver |
| Sync integrity | Tree is listed at a commit, but `fetch_text_file` fetches the branch | Pin content to the same SHA; handle truncated trees and failed/partial refreshes |
| Cache lifetime | Data service uses `st.cache_resource` without a revision argument | Invalidate/reload by manifest revision; never mutate shared resource state with session selections |
| GenBank fallback | No-manifest fallback loads CSVs/plates, not a GenBank index | Support explicit local sequence fallback without requiring an unnoticed prebuilt index |
| BLAST | Separate local/NCBI service; local FASTA may differ from future details | Preserve current contract and label local database provenance |

The checked-in manifest claims 337 main rows, 770 plate rows, and 2,387 indexed GenBank records, generated 2026-09-25. These are manifest values, not independently validated Parquet row counts. This checkout has 2,003 entries in root `genbank/` and 23 product CSVs, plus collection-level sequence directories and 2,347 `genes/` entries. Product CSVs and gene cards are not currently searched by the cache service.

The lower bound `streamlit>=1.28.0` does not guarantee the required APIs. Set a tested supported version with a floor of at least 1.37 for stable dialog and row-selection APIs, and verify signatures against that version. The default and bundled Python environments inspected here lack Streamlit, BioPython, PyArrow, and pytest. No app/test baseline was executed; establish an isolated implementation environment first.

## FreeGenes findings and source decision

Public API/file probes during planning established:

- `freegenes/freegenes.github.io` default branch is `master`. HEAD was `a477c45baac46483f03bc61ccde92867920fabe7`, committed `2023-09-15T18:31:57Z`.
- Its non-truncated recursive tree has 12,670 paths, including 2,338 root GenBank files and 23 product CSVs. GitHub reports repository size of 15,565,304 KB; this is a repository-size indicator, not a measured clone transfer size.
- All 23 product CSVs were inspected: 2,334 rows and 2,295 distinct nonempty IDs. None has a column whose name contains `plate` or `well`. This checks column names, not every possible nested description or old file.
- `genes/BBF10K_003247.html` and `genbank/BBF10K_003247.gb` both returned HTTP 200. The GB LOCUS header declares a 4,451 bp circular record; its description mentions linear DNA, illustrating why conflicting annotations must remain visible.
- The [upstream README](https://github.com/freegenes/freegenes.github.io#readme) identifies Google Sheets as its backend. A direct anonymous CSV export of its linked sheet/GID returned HTTP 400. This does not establish a permission denial or that all other export endpoints fail.

**Decision:** selectively cache FreeGenes metadata and fetch chosen GB files by pinned commit. Search this index together with Reclone on every submitted query. Avoid cloning the entire repository or using GitHub code search per query. GitHub is an upstream snapshot; it is not evidence of live plate inventory.

**Unresolved prerequisite:** locate a usable authoritative FreeGenes plate/well feed, worksheet, or API and verify its schema and identity mapping. Make one bounded investigation in Goal 1. If unavailable, implement the provider interface and truthful source status, use separately labeled Reclone locations, and report this specific unmet data requirement. Do not spend the whole goal probing endpoints or claim the plate sourcing request is fully satisfied.

### Resolution rules

| Field / record | Priority and fallback |
| --- | --- |
| ODC ID and Reclone collection | Preserve Reclone metadata and explicit BBF cross-reference |
| FreeGenes name/description/cloning metadata | Verified available backend record, otherwise pinned FreeGenes CSV snapshot; keep conflicting source values accessible |
| GenBank bytes and derived sequence/features | Usable verified FreeGenes resource, then its last successful cache; local file only with visible fallback status |
| FreeGenes plate/well | Verified location records from an accessible authoritative feed, then cached records from that feed |
| Reclone plate/well | Retain separately as Reclone distribution locations; do not relabel them as FreeGenes |
| Missing / invalid / ambiguous match | Preserve explicit status and conflict details; never name-match silently |

Treat plate name, plate number, well row/column/address, distribution/version, and provider as a single location record. Preserve every valid location. A confirmed absent field can fall back with field provenance; a withdrawn record or invalid upstream GB must not silently be replaced by an apparently current local record.

Use user-facing messages that reflect reality: `GenBank and metadata retrieved from FreeGenes GitHub (snapshot: 2023-09-15)`, `Using cached FreeGenes data; refresh unavailable`, or `FreeGenes plate/well data unavailable; showing Reclone locations`. Use `Information pulled directly from the FreeGenes database` only for a verified direct backend response.

## Service contracts and storage

Keep existing Reclone caches. Suggested additions, created only when needed:

- `services/freegenes_service.py`: source transport, metadata index, exact ID lookup, lazy byte fetch, bounded caches.
- `services/part_service.py`: federated search, aliases, conflict handling, resolved details.
- `services/genbank_service.py`: parsing, structured features, export payloads, viewer adapter.
- `services/restriction_service.py`: deterministic both-strand cut analysis and candidate insert extraction.
- `services/assembly_service.py`: selected fragments/backbone, end compatibility, assembly and feature remapping.
- `ui/part_details.py`, `ui/debug.py`, `ui/builder.py`; a small viewer component only after the viewer spike.
- `scripts/sync_freegenes_data.py`, `data/freegenes/` for curated metadata/manifests, an ignored runtime byte cache, and `data/nomenclature/*.json` for reviewed scheme configuration.

Resolved details should carry a stable key, ODC/BBF aliases, source-specific metadata, all locations, raw GB bytes, parsed record, structured features, topology, provenance, and warnings. Include source revision and content SHA-256. Parse/resolve once per content hash and reuse for all views/exports. Do not put BioPython objects in persisted JSON or expose server filesystem paths in user downloads.

Search result rows carry the same stable key and source badges. Include FreeGenes-only records; deduplicate by verified identity, retain aliases, and surface incompatible aliases. Apply collection/location filters after enrichment. A submitted empty query browses the combined index. Persist the submitted query so opening/downloading from a dialog does not trigger a new upstream search. Do not return early just because the Reclone dataset is empty.

Use configurable metadata refresh TTL (initial proposal: 24 hours), explicit refresh in Debug, a short negative-cache TTL, conditional requests where supported, and timeouts/retries with a bounded total wait. Key raw bytes by source revision and part ID. Build artifacts in staging, validate schemas/hashes/rows, then atomically promote a complete manifest. Failed refreshes retain the previous good snapshot and report the failure. Record separate provenance for local indexed sequences and externally synced metadata.

Scheduled sync may extend the existing workflow with a separate FreeGenes step and clear failure reporting. Community Cloud runtime storage is ephemeral: commit or package the small validated index for cold starts, and lazily rebuild byte caches. Measure index sizes before choosing tracked artifacts. Keep upstream credentials in environment/secrets, never logs, exports, or Debug JSON.

## Detail dialog and viewer

- Use the version-appropriate `st.dataframe` single-row selection API and `st.dialog(width="large")`. Resolve selected positions against the displayed page data, then retain the stable key. Reset selection on query/filter/page changes; handle sorting and repeated opening of the same part.
- Show metadata, source message, physical locations, description, topology, sequence, and full feature list. Missing sequences still allow metadata/location viewing and applicable exports.
- CSV: one resolved part row, with full sequence and explicit JSON columns for features/locations/provenance. Optionally provide a second tidy feature CSV. TXT: metadata, locations, provenance, sequence, and readable feature list. FASTA: parsed sequence with a sanitized ID/header and wrapped lines. GB: original resolved bytes. All formats derive from the same source/hash.
- Browser-check downloads without lost filters or unexpected modal reopening. Provide an explicit `View details` control if table selection is awkward on small screens. Delete the old Part Details route after parity is demonstrated.
- Start with a short viewer spike using the current MIT-licensed [`@teselagen/ove` in tg-oss](https://github.com/TeselaGen/tg-oss), rather than the deprecated `openVectorEditor` repository. Verify the selected package version, license, JSON schema, and Streamlit embedding before committing to it.
- Bundle pinned frontend assets locally. A read-only iframe/component may be sufficient for viewing. Use a real bidirectional component only if selection must update Python state.
- Viewer acceptance: circular and linear display, pan/zoom/sequence scrolling, feature click/hover, readable annotations, strand direction, origin-crossing features, and full feature list fallback. Never draw unknown or linear topology as circular by assumption. Use the same parsed record as exports.
- If the component cannot be made reliable within a bounded spike, use a documented interactive Plotly map as an interim viewer, with feature selection and a synchronized sequence window; report any rotation/scroll limitation.

## Analyze / Generate Part and assembly builder

Use the workbook-derived scheme in `NOMENCLATURE_REFERENCE.md`. Do not auto-fetch its Benchling examples or treat names as DNA sequences.

Analyze mode accepts a chosen database record or user-provided sequence/GenBank plus explicit topology. Compute BsaI and SapI recognition orientation, recognition span, top/bottom cut boundaries, physical 5' overhangs, candidate retained fragments, internal restriction sites, and mapped/unmapped labels. Account for reverse complements, boundary cuts, circular origin-spanning motifs, incomplete linear ends, compound features, and ambiguous bases. Separate facts calculated from DNA from labels inferred through a selected scheme.

`Generate part` initially means export a user-selected, validated digestion fragment with mapped ends and remapped features. If there are multiple plausible site pairs/fragments, require a choice in the product UI. No arbitrary first-pair extraction, automatic mutations, or inferred compatibility from names. Overhang analysis must remain useful even when SapI nomenclature is unmapped.

Assembly mode starts with a ordered list/reorder controls rather than drag-and-drop. User selects backbone/destination vector, scheme, enzyme/assembly level, oriented fragments, and candidate cut boundaries. Check physical compatible ends, scheme slots, internal sites, direction, duplicate/ambiguous ends, CDS frame/scars where relevant, and final circular closure. Distinguish donor plasmids from retained inserts and backbone. Report extra/unresolved fragments and nonunique assemblies.

Construct the predicted sequence with junction bases represented exactly once, remap features across cuts/reverse complements/origin, and verify final length and closure independently. Show a junction table and viewer. Export GB, FASTA, TXT/CSV report, and a portable JSON design including source IDs/hashes, orientation, cuts, scheme version, and enzyme. A restored design must not silently pick newer source sequences. Export failed/ambiguous designs as reports only, not a validated assembled plasmid.

Basic in-silico compatibility is not a ligation-efficiency prediction. Additional uLoop levels/enzymes require an explicitly sourced scheme. The workbook's four-base scheme does not supply SapI level rules.

## Staged implementation and acceptance

| Goal | Work | Completion evidence | Planning token allowance |
| --- | --- | --- | --- |
| 1: UI and FreeGenes | Establish environment; bounded location-source investigation; selective metadata/cache provider; resolved GB priority; federated search; clickable details/downloads; Debug merge; Analytics/Part Details removal; README refresh | Mocked source/conflict/fallback tests; search-click-feature-download flow; paginated/sorted selection; absent GB behavior; Debug/BLAST regressions; bounded public source probe | About 25,000-35,000 |
| 2: Viewer and part analysis | Bounded viewer spike; shared record adapter; BsaI/SapI analysis; reviewed workbook config; candidate extraction and part exports; working Interactive Builder Analyze mode | Browser viewer behavior; hand-calculated cuts for both enzymes/orientations; linear/circular edge fixtures; compound features; ambiguous/internal sites; scheme/unmapped labels; export reparse | About 15,000-25,000 |
| 3: Assemble Plasmid | Ordered selection, explicit backbone and scheme/level, end/slot validation, predicted sequence, feature remapping, design save/restore and exports | Known valid sequence/junction fixtures; failed/nonunique assembly cases; final closure/length; reverse-complement feature tests; browser end-to-end build/save/restore/download | About 25,000-40,000 |

These are suggested working budgets, not measured usage or guaranteed completion limits. Do not start all three under one small hard budget. If approaching a hard cap, preserve work and record remaining criteria without claiming completion. Viewer integration and upstream location access are the main uncertainty drivers.

Tests should use small offline fixtures for exact aliases, duplicate IDs, multiple locations, absent/unavailable/invalid sources, timeout/rate limit, partial cache refresh, source-specific fallback, identical export content hashes, reverse cuts and circular boundaries. Existing tests cover normalization, search enrichment, sync helpers/mocks, and BLAST contracts; extend them without live network dependence. Use browser checks for selection/dialog/components that unit tests cannot establish.

## Progress and open decisions

- [x] Inspect repository, existing service boundaries, sync workflow, tests, and dependency declarations.
- [x] Verify public FreeGenes branch/revision, sample GB/card, and all product CSV headers.
- [x] Extract workbook nomenclature with exact cell provenance.
- [x] Save plan, repository instructions, and staged goal prompts.
- [ ] Find an accessible authoritative FreeGenes plate/well source and document identity/version semantics.
- [x] Establish isolated `.venv` dependencies and verify the 10-test pre-change baseline.
- [x] Implement Goal 1 within its documented unavailable-location allowance; see `GOAL1_VALIDATION.md`.
- [x] Implement Goal 2 (viewer and part analysis); see `GOAL2_VALIDATION.md`.
- [ ] Implement Goal 3 (plasmid assembly).

Defaults chosen: dialog for details; selective external cache; merged search includes upstream-only records; ordered builder controls; raw sequence analysis remains available without a nomenclature match. Only exact authoritative location endpoint/access and additional SapI/uLoop scheme rules remain necessary source inputs. Viewer choice remains conditional on its spike.

## Primary references

- [FreeGenes repository / README](https://github.com/freegenes/freegenes.github.io) and [checked revision](https://github.com/freegenes/freegenes.github.io/commit/a477c45baac46483f03bc61ccde92867920fabe7).
- [Product CSV directory at checked revision](https://github.com/freegenes/freegenes.github.io/tree/a477c45baac46483f03bc61ccde92867920fabe7/product-csvs) and [sample GenBank](https://github.com/freegenes/freegenes.github.io/blob/a477c45baac46483f03bc61ccde92867920fabe7/genbank/BBF10K_003247.gb).
- [Backend sheet linked by FreeGenes](https://docs.google.com/spreadsheets/d/1LZCXzBtgey9xv5OH7YGYgp8UMJ27Eyj1aF9IhAW6M6o/edit#gid=954552604). Public Genes metadata is accessible; authoritative physical plate/well schema remains unavailable.
- [Streamlit 1.49 dialog](https://docs.streamlit.io/1.49.0/develop/api-reference/execution-flow/st.dialog) and [dataframe selection](https://docs.streamlit.io/1.49.0/develop/api-reference/data/st.dataframe).
- [TeselaGen tg-oss](https://github.com/TeselaGen/tg-oss), [NEB cleavage reference](https://www.neb.com/en/tools-and-resources/selection-charts/enzymes-with-nonpalindromic-sequences), and the user workbook documented in `NOMENCLATURE_REFERENCE.md`.


## Goal 1 completion checkpoint — 2026-10-02

Implemented Home/Search & Browse/BLAST Search/Debug navigation, merged diagnostics/freshness/data exports, federated search, and the reusable details dialog with source-specific locations, all GenBank features, sequence and five export formats. Removed Analytics, Part Details, and Data Management routes. Original collection data remains unchanged; no deployment was performed.

The bounded backend investigation found that public Google Visualization CSV queries **do** work for the Genes worksheet when `headers=1` is explicit. The provider caches 2,346 usable backend records alongside 2,334 GitHub metadata rows and 2,338 indexed GenBank paths. The merged index currently has 2,349 identities. Public metadata omits raw sequence/internal-note columns. Backend file links are checked against exact BBF identity, then fetched through the commit-pinned GitHub resource. Conflicting links block sequence selection.

Genes has no physical plate/well columns. Packaging location entries are N/A; Collections locations are unspecified. Neither establishes an authoritative physical plate/well feed. This remains an unmet user data requirement: all independent Goal 1 work is complete, but FreeGenes location sourcing cannot be claimed complete. The provider exposes an unavailable state and supports complete sourced location records when a verified feed becomes available. Reclone distribution plate maps remain separate.

The final packaged FreeGenes index is gzip-compressed to 620,214 bytes and validated before manifest promotion. Debug refresh succeeded in the browser at 2026-10-02T12:15:22Z; the upstream GitHub commit is still `a477c45baac46483f03bc61ccde92867920fabe7` (2023-09-15). Runtime caches are ignored; original local BLAST sequences remain separate. Both providers publish atomic metadata generations, and existing flat Reclone caches remain supported.

Validation: **50 focused/regression tests passed**, the existing smoke script passed with Windows UTF-8 mode, and browser checks passed for cold start, merged search, second-page sorting/selection, native/explicit dialog dismissal, reopening, sequence/features/source locations, consistent downloads, Debug refresh and BLAST page input validation. Actual local BLAST alignments remain unverified because external executables are absent. See `GOAL1_VALIDATION.md` for exact commands, source findings, download hash and browser evidence. Viewer and builder work is unstarted.


## Post-Goal 1 correction — 2026-10-02 (current behavior)

The user's correction supersedes the earlier upstream-only search requirement: only Reclone identities are searchable/browsable; matching FreeGenes records enrich them. FreeGenes-only identities and collection filters are excluded. Current browse count: 337.

Replaced checkbox-only dataframe selection with a locally bundled Streamlit v1 table component. Clicking any cell in a row sends the stable part key; sorting preserves identity. Rows also support Enter/Space. Paging, the explicit View details fallback, and search CSV downloads remain.

Part Details hides metadata retrieval banners/timestamps, parser warnings, circular/linear-description warnings, missing-location notices, and all source/provenance sections/captions. Actual locations, metadata values, sequence/features and exports remain. Source diagnostics are preserved internally/in exports; ambiguous or unavailable DNA is still handled explicitly.

Validation: 52 tests passed; smoke script and diff check passed. Browser verified the 337-part inventory, sorted name-cell click, a different cell reopening the same part after native dismissal, hidden messages, and all five detail download controls. No viewer/builder work started.

## Goal 2 completion checkpoint — 2026-10-02

Added locally pinned TeselaGen OVE 0.8.42 to shared Part Details and Builder views. The successful spike established the inclusive-coordinate adapter and required iframe sizing delay. Readiness retries, deferred mounting until visible width, responsive remounting, and module-level registration make the component work in dialogs and collapsed expanders. Actual browser checks covered circular/linear zoom, rotation, sequence scrolling, forward/reverse feature selection, and an origin-crossing reverse join. OVE's joined-feature selection/size covers the bounding span; exact segments are preserved and the difference is explained beside the viewer.

Interactive Builder now provides Analyze / Generate Part from Reclone records, pasted DNA/FASTA, or uploaded GenBank/FASTA/DNA. Explicit topology, enzyme and scheme settings produce both-strand BsaI/SapI recognition/cuts, physical overhangs, candidates and internal-site reports. The workbook configuration retains cell/hash/version provenance; SapI and missing uLoop rules stay unmapped. Fragment choice is explicit, invalid fragments remain reports only, and source changes invalidate previous analysis. Generated linear GB/FASTA/CSV/TXT and feature CSV share one validated reparsed record, with correctly projected features and source/end provenance. Original sequences and physical inventory locations remain unchanged.

Validation: **84 tests passed**, legacy smoke script and diff check passed. Browser downloads agreed across formats; an annotated generated GB preserved reverse joined extraction independently. Real database DNA matched the Part Details source hash. Ambiguous DNA and a custom interval containing an additional cut blocked generation. Final cold-start builder flow and console/server checks passed. Exact commands, source references, hashes, artifact names, scientific conventions, limits and the requirement-by-requirement audit are in `GOAL2_VALIDATION.md`.

No required Goal 2 acceptance item remains open. Full multi-part assembly, compatibility/coding-junction checks, reverse-oriented assembly products, and design save/restore remain Goal 3. FreeGenes physical location sourcing remains the unchanged Goal 1 limitation. No commit, push, deployment, or original dataset rewrite was performed.

## Part Details layout correction — 2026-10-02

The Search & Browse row-click dialog now places all applicable downloads to the right of the part name, identifiers, collection and description. Metrics and the white-background interactive viewer follow immediately, before location and metadata tables. Removed the separate raw DNA code block; DNA remains available inside the viewer and in exports. Full GenBank features and missing-sequence metadata exports remain available.

Validation: 84 tests passed, UTF-8 smoke checks passed, and application diff whitespace checks passed. Browser verified name-cell activation, five right-side download buttons, no raw DNA block, white viewer and reverse `ori` selection (589 bases, 2623–3211), a successful FASTA download without dismissing the dialog, and closing/reopening the same row without console errors. Screenshot: `part-details-downloads-right.png` in the task visualization directory. Changes remain local; the temporary verification server was stopped afterward.

## Search inventory and recognizable names — 2026-10-02

The shipped master cache matches all 337 rows of `odc_plasmids.csv`. All 24 records named polymerase already matched the query, but display names were replaced by FreeGenes codes (Taq → THEAQpolA; Bst → Bstpol/BstpolLF). The full query has 36 matches across two default pages. Reclone names now label results, details and builder choices; FreeGenes source names and descriptions remain searchable without changing preferred metadata or GenBank precedence. Added an explicit result range/page caption and upgraded the search index version to invalidate saved results; the cached constructor also changes version.

The inventory audit additionally found Reclone plate-only records `BBF10K_000483` (beta-galactosidasealphadomain, G4) and `BBF10K_003338` (meffCP, H4), both in `Open Reporters Collection/Platemaps/ORC-v1_0.csv`. Those identified plate records now participate in the inventory with source-path/version provenance, yielding 339 parts. Plate records cannot rewrite master-list identities; conflicting location identities remain excluded from filters and search aliases. FreeGenes-only records remain excluded. Home counts use the same inventory. Source CSVs, GenBank files and BLAST cache remain unchanged.

Validation: 90 unit/integration tests passed, including a read-only audit of every shipped master-list name, all named polymerases, plate-only membership, source-description aliases, identity conflicts and saved-result invalidation. UTF-8 smoke and diff checks passed. Browser verified full Taq/Bst names, Taq row-click details with viewer/five downloads, and page-two polymerase matches. The FreeGenes Taq sequence request was unavailable during this browser check; the existing labeled local fallback worked. Screenshot: `polymerase-search-fixed.png` in the task visualization directory. Changes remain local for review.

## Cloud cached-record KeyError correction — 2026-10-02

Reproduced `KeyError: 'display_name'` at the Search & Browse dropdown by retaining submitted results and removing the derived field from a cached PartService record. The dropdown previously indexed the cached record directly. It now uses labels from its displayed result page. Shared `part_display_name` derives Reclone names for legacy records and falls back to stored names/identifiers in details, builder and fresh search results.

`PART_INDEX_VERSION` v4 is now an explicit hashed loader argument and part of saved-result revisions. Each PartService carries its schema version; a cache validator rejects instances with a missing/outdated version. Regression checks prove reconstruction when only the schema changes or an incompatible instance survives, without changing dataset manifests. Legacy-record tests cover an existing dropdown, details, builder and a new search. All 94 tests, UTF-8 smoke and diff checks pass. These checks reproduce the reported failure locally; no direct access to Cloud runtime logs was used.

Browser validation: Taq search and dropdown selection survived reruns; View details opened Taq with the viewer and all five downloads. Closing details and navigating to Builder rendered its database dropdown without console errors. Screenshot: `search-cache-keyerror-fixed.png` in the task visualization directory. The temporary server was stopped; changes remain local and have not yet been pushed to Cloud's source branch.
