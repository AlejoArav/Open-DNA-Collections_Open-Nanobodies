# Goal 2 validation — 2026-10-02

The viewer, Analyze / Generate Part services, and builder UI are implemented locally. Full assembly remains Goal 3. No original dataset was modified, and no commit, push, or deployment was performed.

## Implementation and scientific conventions

- Viewer: `ui/sequence_viewer/`, with `services/viewer_service.py` adapting shared parsed GenBank spans from half-open to OVE inclusive coordinates. Pinned `@teselagen/ove` 0.8.42 at Git commit `77702b6e8ebedb554af14ee0ea90cf04fc101954`; npm integrity, asset SHA-256 hashes, and the pinned MIT license are in `vendor/manifest.json` and `LICENSE-TeselaGen.txt`. JS is 7,398,558 bytes; CSS is 1,962,246 bytes, including embedded font resources. Frontend assets load locally. No Node build or added Python runtime dependency is needed.
- The bounded spike succeeded. OVE needed initialization after Streamlit assigned frame height, followed by two animation frames. Readiness retries until first host render. A resize observer defers mounting in collapsed expanders until width is positive and remounts on changed width. The component is registered once at module import, avoiding dialog-rerun registration failures. Source viewers are read-only; menu/toolbar/import/drop are disabled. Linear/unknown topology offers only a linear map; circular records also offer a linear view.
- Enzyme geometry was checked against [NEB's nonpalindromic cleavage chart](https://www.neb.com/en/tools-and-resources/selection-charts/enzymes-with-nonpalindromic-sequences): BsaI `GGTCTC(1/5)`, SapI `GCTCTTC(1/4)`. For forward recognition beginning at zero-based `s` with motif length `m`, top/bottom boundaries are `s+m+1` and `s+m+5` (SapI: `s+m+4`). For reverse recognition, they are `s-5` and `s-1` (SapI: `s-4` and `s-1`). Circular coordinates are reduced modulo record length while unwrapped values remain in reports.
- Both enzymes produce 5′ overhangs. The reference word between top and bottom boundaries is the right fragment's physical left end; its reverse complement is the left fragment's physical right end. Workbook labels refer to reference junction words; physical words and their separate lookup results are retained. Nomenclature never changes DNA.
- The checked-in `user_workbook_four_base_v1.json` contains eleven junctions, four source-backed module relationships, worksheet/cell provenance, workbook SHA-256 `f45eac44b2858fad7e74d20739ffd179a046aa32d4e9dc663ca0aa08eaa7786d`, and explicit mapping limitations. There is no dependency on the user's Downloads path at runtime.
- Generated sequence is the source/reference top strand from selected top-cut boundary to the next selected boundary, including the left overhang and excluding the right complementary-strand overhang. Reports also contain the physical ends and bottom strand (5′→3′). Double-stranded core length must be positive. Invalid candidates remain reports, not validated parts.
- Feature projection preserves original biological part order, strand, join/order operators, and retained fuzzy endpoints. Clipped endpoints become exact cut positions. Clipped features are noted, stale CDS translation removed, and codon start adjusted only when original positions/strand establish the reading frame. Unknown/remote/out-of-range features are omitted with a report notice. Generated virtual fragments have no physical inventory locations. GenBank sequence and mapped feature extraction are checked independently after serialization/reparse. Full source/cut/end provenance is retained in CSV/TXT; GB COMMENT carries a line-wrapped base64 UTF-8 JSON record to survive GenBank wrapping.

## Automated checks

Latest executed full suite:

```powershell
.\.venv\Scripts\python.exe -m pytest tests/unit tests/integration -q
# 84 passed in 4.49s

.\.venv\Scripts\python.exe -X utf8 test_app.py
# Smoke tests passed: 337 main rows, 770 platemap rows, ODC_0007 lookup, BLAST contract
```

New coverage includes hand-derived BsaI/SapI forward/reverse cuts and duplex ends; circular recognition/overhang/fragment wraparound; single-cut circle; incomplete linear ends; a cut at linear boundary zero; overlapping staggered cuts and coincident recognition pairs; ambiguous recognition/ends/interior; internal complete and possible cuts; bounded ambiguous-match expansion; explicit topology and single-record input; exact workbook labels/cells and unmapped SapI; source immutability; clipping/reverse/compound/fuzzy/remote feature cases; CDS frame handling; generated GB reparse; invalid-candidate reports; export DNA/hash agreement; viewer adapter/asset hashes; and actual frontend readiness/layout JavaScript executed under Node.

Streamlit AppTest additionally covers builder inventory membership, missing DNA, explicit topology, mandatory fragment choice, disabled invalid generation, five generated formats, and stale-output invalidation on changed inputs or source manifests. Existing Goal 1 search, source precedence, provenance hiding, whole-row handshake, cache, and BLAST contract tests remain included.

## Browser evidence

Real Streamlit app, Python 3.12.10 / Streamlit 1.49.1, at `http://127.0.0.1:8501/`:

- Reclone inventory count remained 337. Searching `ODC_0007` and clicking its name cell opened `BBF10K_003247` / 9N7polA. Part Details showed 4,451 bp and 12 source features with no requested diagnostic/provenance notices. Close/reopen worked after the component registration fix.
- Circular and sequence panels each measured 594 px in the wide dialog. Clicking `ori` selected 589 bases, 2623–3211, matching the downloaded GenBank's reverse-strand `[2622:3211]` location. Circular plus/minus zoom and rotation arrow changed slider positions; sequence scrolling advanced displayed base ranges. Downloaded source GB SHA-256: `5a93497f1f4b445a37c8714a271165422bd9077fb0a617af3ae8d3fc8c80b469`.
- BsaI hand fixture (28 bp) displayed forward boundaries 7/11 and reverse boundaries 17/21. Generation was disabled until explicit selection. The 10 bp insert exported `GGAGTTTTTT`, with reference labels A→E and physical ends GGAG/AAGC. Browser-downloaded GB/FASTA/CSV/TXT agreed; CSV recorded the exact downloaded GB hash `a827cd1dc8bff18ec9ee6a8f1432a86e7fce64e70605d63cd8495c41cdd76ff0`.
- Uploaded `tests/fixtures/goal2_circle.gb` (1,000 bp) exercised forward CDS and reverse joins. Expanding the source viewer produced two 409 px panels; selecting the origin-crossing reverse join selected 140 bases, 901→40. The feature panel retained both `(901–1000)(1–40)` and `(29–30)(34–36)` spans with strand −1. With an explicit linear override, no Circular Map tab appeared; linear zoom changed to 5%, and selecting Forward CDS selected 610 bases, 91–700.
- Generated fixture fragment `fragment_bfc6f3e8fc76` had 10 bases and two retained reverse features. Browser-downloaded GB reparsed as linear; the retained compound location was `join{[6:9](-), [1:3](-)}` and independently extracted `AAATC`, exactly matching the original feature. Clipped origin annotation carried an explicit partial-feature note.
- SapI hand fixture displayed forward boundaries 8/11 and reverse boundaries 17/20; physical ends and three-base junctions were reported with `Unmapped in selected scheme`. Explicit selection generated 9 bp `GCATTTTTT`; the GB download completed.
- Downloaded SapI GB SHA-256: `c736af9e28c2b07d5a870aadae8a9e672691e7aa8c30de67fc7b0a7e55fb9748`. It reparsed as linear with the expected nine bases. A selected ambiguous SapI candidate showed its ambiguity error and disabled generation; its downloaded selected-fragment report preserved the selected boundaries and validation errors.
- Actual Builder database input `BBF10K_003247` analyzed 4,451 bases with complete BsaI boundaries 27/31 (`AATG`) and 2359/2363 (`GCTT`). Its downloaded analysis sequence SHA-256 matched the Part Details GB DNA, with the same source identity/provenance. Selecting boundary `cut_27_31` at both ends attempted a full-length linearization; the second internal site `site_2364_-1` was reported and generation disabled.
- After a final cold server start with current imports, the analyze/select/generate flow produced the expected ten-base fragment, both nonzero-width viewer panels, the empty-feature fallback, all five generated download controls, and no console errors. The final server handle remained live with only normal cache-load messages. `git diff --check` passed. HEAD remained `47cf0398a94439b2ac554920d2998525e5146499`.
- Browser images were saved as `goal2-part-details.png` and `goal2-builder.png` in the task's visualization output directory. No required Goal 2 browser acceptance item remains open.

## Completion audit

Post-completion viewer correction: the shared Part Details/Builder iframe now forces a light color scheme, white background and dark text without inheriting the app's dark class. The frontend protocol regression sends a dark host theme. All 84 tests and the UTF-8 smoke script passed again. A real Search & Browse name-cell click opened the viewer with computed white background (`rgb(255, 255, 255)`) and dark text; selecting `ori` still selected 589 bases, 2623–3211. Screenshot: `viewer-white-background.png` in the task's visualization output directory.

Commit preparation: all three staged vendor assets match their pinned SHA-256 hashes. `.gitattributes` preserves their bytes across Windows/Linux checkouts. Staged whitespace validation excludes those unchanged upstream assets, which contain vendor whitespace; application code and documentation pass the check.

| Requirement | Evidence |
| --- | --- |
| Preserve Reclone membership, row activation, hidden notices, precedence, exports and BLAST | Existing tests included in the 84-test suite; real name-cell click, hidden-notice check, close/reopen, unchanged 337-part inventory and source GB; tracked app diff adds only viewer/builder hooks |
| Shared-record interactive viewer; topology, zoom, scrolling, selection and strands | Adapter tests; actual circular and linear interactions on 9N7polA/fixture; selected ranges independently checked against downloaded GB |
| Compound and origin-crossing features | Adapter and projection tests; reverse origin selection 901→40; separate viewer segments; generated compound GB reparse/extraction `AAATC` |
| Pinned local frontend assets and license; bounded spike | Manifest hash tests; exact vendor bytes matched the integrity-verified npm tar; MIT license at pinned source revision; real layout and readiness fixes |
| Database and user input with explicit topology | Input tests, AppTest, actual 9N7polA database analysis, pasted fixtures and uploaded GenBank; override preserved original source |
| BsaI/SapI recognition/cuts/physical ends/internal sites/candidates | NEB chart; hand-derived fixtures for both strands and topologies; real BsaI/SapI tables and rejected internal/ambiguous selections |
| Versioned source-backed nomenclature without invented mappings | Eleven exact four-base entries with cells/hash/version; current supplied workbook hash rechecked; explicit unmapped SapI in tests and browser |
| Explicit selection, feature remapping, validated GB/FASTA/CSV/TXT and reports | Disabled placeholder/invalid controls; generation tests; four actual browser downloads compared; annotated GB independently reparsed; invalid-candidate report downloaded |
| Regressions, browser checks, documentation and goal boundary | 84 passing tests, UTF-8 smoke script, diff check, current screenshots/logs; README/plan/reference guidance updated; original datasets/HEAD preserved; no assembly, commit, push or deployment |

## Deliberate limits

OVE represents compound segments accurately but its selection/size column uses their bounding span (including gaps). This is explained beside the viewer; exact original segments remain in the source/generated feature table and exports. Unknown strands have no directional arrow; unmappable positions retain a full-list fallback. Physical sticky-end geometry is in reports, while the viewer visualizes the exported reference strand.

Input limits are 5 MB, 500,000 bases, and 10,000 exact/possible recognition matches per orientation. Ambiguous/incomplete sites are reported and excluded from certain cut boundaries; possible internal cuts and ambiguous retained DNA prevent validation. Recognition retention only marks a donor-backbone candidate, not a proven biological identity. No digestion efficiency, domestication, primer design, protein-function guarantee, full assembly, or inferred SapI/uLoop level rules are claimed. FreeGenes plate/well availability and missing external BLAST executables remain unchanged Goal 1 limitations.
