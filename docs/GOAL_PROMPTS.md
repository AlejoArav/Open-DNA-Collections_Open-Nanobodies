# Staged goal prompts

These prompts start implementation only when submitted by the user. No goal was started during planning. Suggested budgets are provisional; configure a hard budget explicitly if wanted. Read `AGENTS.md` and the relevant plan section once, then work from compact checkpoints rather than repeatedly rescanning the repository.

## Goal 1: UI and FreeGenes integration

```text
/goal Implement Goal 1 in docs/IMPLEMENTATION_PLAN.md in this checkout. Read AGENTS.md and the plan's source findings/resolution rules. Aim for approximately 30,000 tokens; prioritize the complete UI/data flow and focused validation.

Move Home freshness/diagnostics and Data Management into Debug. Remove Analytics. Make Search & Browse rows open a reusable wide part-details dialog and remove the standalone Part Details route after parity. Preserve query/filter/page state and resolve clicked rows by stable identity. Show metadata, all source-specific plate/well locations, sequence, full GenBank features, provenance, and GB/CSV/FASTA/TXT downloads from one resolved record.

Add a separate selective FreeGenes metadata index and lazy GenBank byte cache. Search both Reclone and FreeGenes on each submitted query, including upstream-only results. Prefer verified FreeGenes GenBank and FreeGenes location records over local equivalents. Preserve Reclone collection/ODC aliases and separate Reclone distribution locations. Handle conflicts, multiple locations, invalid files, rate limits, missing records, offline cache, and partial refreshes explicitly. Pin GitHub fetches to the recorded commit, retain exact GB bytes, and invalidate app caches by revision.

Make one bounded attempt to identify an accessible authoritative FreeGenes plate/well feed. The checked GitHub snapshot is from 2023 and its product CSVs have no plate/well columns. If the feed remains unavailable, finish all independent work, provide an explicit unavailable state and labeled local fallback, and report the missing feed as an unmet data requirement. Do not claim a GitHub snapshot is a live database or invent locations.

Establish an isolated development environment, update the tested Streamlit dependency/API usage and README, and preserve BLAST behavior. Run focused mocked source/search/export tests, existing regression tests, and a browser smoke check of search-click-dialog-download-Debug. Record results and outstanding acceptance items in the plan. Stop at Goal 1's boundary; viewer and builder work have separate goals. Do not rewrite original datasets or deploy as part of this goal.
```

Before Goals 2–3, follow the latest correction checkpoint in the implementation plan: Reclone-only inventory enriched by FreeGenes, whole-row activation, and simplified Part Details. The original Goal 1 prompt above is historical.

## Goal 2: Viewer and Analyze / Generate Part

```text
/goal Implement Goal 2 in docs/IMPLEMENTATION_PLAN.md after checking Goal 1's actual status. Read AGENTS.md and docs/NOMENCLATURE_REFERENCE.md. Aim for approximately 20,000 tokens. Add a reliable interactive plasmid/sequence viewer to the shared details renderer, using a bounded @teselagen/ove integration spike or the documented interactive fallback. Verify zoom/scroll, feature selection, linear/circular topology, compound/reverse/origin-crossing features, and export consistency in the browser.

Add Interactive Builder with a working Analyze / Generate Part mode. Calculate BsaI/SapI both-strand recognition and cleavage boundaries, physical overhangs, internal sites, candidate retained inserts, and scheme-specific labels. Convert the workbook reference to versioned configuration with cell/hash provenance. Keep unmapped SapI labels explicit. Require a choice when candidate boundaries are ambiguous, preserve original sequences, remap retained features, and export validated parts and analysis reports. Test hand-calculated forward/reverse/circular/linear edge fixtures and GB reparse. Update plan checkpoints and report exact limitations. Leave full plasmid assembly for Goal 3.
```

## Goal 3: Assemble Plasmid

```text
/goal Implement Goal 3 in docs/IMPLEMENTATION_PLAN.md using the validated part-analysis services from Goal 2. Read AGENTS.md and the nomenclature reference. Aim for approximately 30,000 tokens. Add ordered database-part selection, explicit backbone/fragment/orientation choices, sourced scheme/assembly-level settings, physical end and slot validation, internal-site and coding-junction checks, and final circular closure validation.

Build the predicted plasmid sequence with every junction represented exactly once and correctly remapped features. Show the junction report and interactive viewer. Export GenBank, FASTA, CSV/TXT reports, and portable JSON designs with source IDs/revisions/hashes. Save/restore must preserve the selected sequence versions. Reject ambiguous or incompatible assemblies visibly without silently mutating DNA or presenting them as validated products. Test known exact-sequence assemblies, reverse fragments, feature remapping, failed/nonunique ends, closure, and browser save/restore/download flows. Update plan checkpoints and report the implemented scope and validation evidence.
```
