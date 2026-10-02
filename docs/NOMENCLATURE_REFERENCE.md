# Nomenclature reference

Inspected on 2026-10-01. Source: `C:/Users/Alejandro/Downloads/DNA_NOMENCLATURE_FOR_GOLDENGATE_ASSEMBLY.xlsx`, worksheet `DNA modules scheme`.

Workbook SHA-256: `f45eac44b2858fad7e74d20739ffd179a046aa32d4e9dc663ca0aa08eaa7786d`.

The workbook was read without modification. It contains one schematic worksheet, with content in rows 2-52 and formatting extending to row 1001. Its 60 formulas are hyperlinks, not sequence calculations. Source text and links are reference data, not implementation instructions.

## Extracted four-base junctions

| Junction label | Sequence as written | Label cells | Sequence cell |
| --- | --- | --- | --- |
| A | GGAG | B4:B5 | B6 |
| B | TACT | D5:D6 | D7 |
| C | AATG | L9:L10 | L11 |
| D | AGGT | N9 | N11 |
| E | GCTT | T6:T7 | T8 |
| F | CGCT | V5:V6 | V7 |
| T1 (N1) | CCAT | F6:F7 | F8 |
| T2 (N2) | GTCA | H7:H8 | H9 |
| T3 (N3) | TCCA | J8:J9 | J10 |
| T4 (N4) | TTCG | P8:P9 | P10 |
| T5 (N5) | CGGC | R7:R8 | R9 |

The schematic also uses N3 at K9, C1 at O9, C2 at Q8, and C3 at S7 to describe domains. Treat possible C1/C2/C3 junction aliases as a review item rather than silently adding them as equivalent lookup keys.

## Module relationships visible in the schematic

- A-B: promoter / 5' UTR (C5); example links at C6:C11.
- B-C: RBS modules, including BC-prefixed examples at E34:E41.
- C-D: CDS (M10), with a listed set of example constructs at M11:M30.
- E-F: terminator (U6); examples at U7:U9.
- Intermediate T/N junctions permit tags, cleavage domains, reporters, and linkers. Composite modules span more than one junction, such as N1C, BN2, CE, and DE. Do not force every part into one immediate neighboring slot.
- Destination vectors and selection cassettes occupy separate schematic sections. Their actual digestion boundaries must come from their sequences.

These are interpretations of the sheet's labeled layout, not validated sequence identities. Benchling links identify examples; they were not fetched and cannot substitute for actual GenBank sequences.

## Limits and implementation contract

1. A4/W5 group the scheme under MoClo/CIDAR/uLoop. This does not establish that every variant of those systems uses identical rules.
2. Every listed junction is four bases. The workbook does not explicitly identify BsaI, SapI, their cleavage geometry, SapI three-base mappings, or alternating uLoop level rules.
3. Store the extracted table as a versioned scheme configuration during implementation, with the workbook hash and cell provenance. Suggested initial ID: `user_workbook_four_base_v1`.
4. Enzyme analysis is independent of labels. BsaI recognition is `GGTCTC(1/5)` and SapI is `GCTCTTC(1/4)`, giving four- and three-base 5' overhangs respectively. Verify against [NEB's cleavage reference](https://www.neb.com/en/tools-and-resources/selection-charts/enzymes-with-nonpalindromic-sequences) and test both orientations.
5. Normalize computed physical ends into a documented 5'-to-3' convention before scheme lookup. Show forward and reverse-complement interpretations explicitly when orientation matters.
6. For SapI or an unrecognized four-base junction, display the measured overhang with `Unmapped in selected scheme`. Do not invent labels.
7. Selecting a nomenclature label must never change the underlying DNA. Coding junctions require explicit reading-frame/scar checks; the workbook's gly/ser/met notes are context, not universal guarantees.
8. Convert this small curated reference into application configuration once. Do not make Streamlit depend on an absolute Downloads path or reparse the workbook on each rerun.
