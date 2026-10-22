# Report Schema

## Verdicts

| Verdict | Meaning |
|---|---|
| `HIGH_CONFIDENCE_DUPLICATE` | At least two strong range matches and global guard passed. Report as a high-confidence spectral duplicate. |
| `SUSPECTED_DUPLICATE` | At least two loose range matches, but below the high-confidence rule. Report as a review candidate, not a confirmed duplicate. |
| `NOT_DUPLICATE` | Did not meet either duplicate rule. |
| `INSUFFICIENT_INFORMATIVE_RANGES` | Fewer than two ranges had enough signal/peaks to support a decision. |
| `NOT_COMPARABLE_DIFFERENT_SOLVENT` | Different solvents. Do not interpret the numerical score as duplicate/non-duplicate evidence. |

## CSV outputs

### `nmr_duplicate_pairs.csv`

One row per high-confidence duplicate pair. Important columns:

- `name_a`, `name_b`: dataset labels.
- `doc_a`, `doc_b`: source document/root.
- `weighted_range_score`: weighted mean of informative range scores.
- `informative_ranges`: number of ranges that contributed.
- `strong_ranges`: number of ranges with score >= `strong_score`.
- `loose_ranges`: number of ranges with score >= `loose_score`.
- `bad_ranges`: informative ranges below 0.65.
- `global_cosine`, `global_bin_cos_020`, `global_deriv_cos`, `global_peak_coverage_020`.
- `verdict`, `range_summary`.

### `nmr_suspected_pairs.csv`

Same schema as duplicate pairs, but only `SUSPECTED_DUPLICATE` rows.

### `nmr_all_pairs.csv`

All pairwise comparisons, including `NOT_DUPLICATE`, `NOT_COMPARABLE_DIFFERENT_SOLVENT`, and `INSUFFICIENT_INFORMATIVE_RANGES`.

### `nmr_range_evidence.csv`

One row per compared pair and informative range:

- `range`
- `score`
- `coverage`
- `cosine`
- `energy_a`, `energy_b`
- `peaks_a`, `peaks_b`
- `weight`

### `nmr_dataset_qc.csv`

One row per parsed 1H dataset:

- source document/root
- dataset label
- nucleus, solvent, pulse program, field
- `source` (`1r` or `fid-magnitude-fallback`)
- `parse_ok`
- SHA-256 of `fid` and `pdata/1/1r`
- dataset path

### `nmr_excluded.csv`

Datasets excluded from the 1H comparison, including:

- 13C/DEPT/2D experiments
- parse failures
- datasets with insufficient readable spectral data

## XLSX workbook

The XLSX contains:

- `Summary`
- `High_Confidence_Pairs`
- `Suspected_Pairs`
- `Dataset_QC`
- `All_Pairs`
- `Method`

## Scope reporting for multiple sources

For a combined run, report counts separately for:

- internal pairs in source A;
- internal pairs in source B;
- internal pairs in source C;
- each pair of sources (A–B, A–C, B–C).

Do not collapse internal and cross-source pairs into one undifferentiated count.

## Minimum reporting sentence

Use wording close to:

> “X high-confidence spectral duplicate pairs and Y suspected pairs were found among N 1H spectra. These are spectral matches; compound identity requires 13C/2D NMR, MS, and sample-record confirmation.”

For cross-solvent or cross-nucleus pairs, state that they were marked non-comparable rather than “not duplicate”.

## Scope breakdown

`nmr_scope_breakdown.csv` reports pair counts separately for each source pair and comparison type:

- `NEW_NEW`
- `NEW_HISTORY`
- `HISTORY_HISTORY`

This is the verification table for incremental uploads: a newly added file should produce `NEW_HISTORY` rows against every previously registered source.
