# Method: Range-Based 1H NMR Duplicate Detection

## Goal

Find pairs of 1H spectra whose peak positions and local spectral shapes agree across multiple chemical-shift ranges. The method is designed for batch screening of Bruker datasets, including datasets referenced by MestReNova documents.

It is a **spectral duplicate** test, not a structure-elucidation or identity test.

## Data handling

- Prefer `pdata/1/1r`; use `pdata/1/procs` for the ppm axis (`OFFSET`, `SW_p`, `SF`).
- If `1r` is absent, the script falls back to a raw-FID magnitude spectrum. QC records `source=fid-magnitude-fallback`.
- Read nucleus and solvent from `acqus` (`NUC1`, `SOLVENT`, `PULPROG`, `SFO1`).
- Exclude 13C, DEPT, and 2D datasets from the 1H comparison; list them in `nmr_excluded.csv`.
- Compute SHA-256 for `fid` and `pdata/1/1r` and record it in QC.

## Preprocessing

- Interpolate to `0.0–12.0 ppm` on a `0.002 ppm` grid.
- Subtract the 5th percentile as a flat baseline and clip negative values to zero.
- Mask common solvent/water regions:
  - CD3OD: `3.24–3.38 ppm`
  - CDCl3: `7.20–7.32 ppm`
  - DMSO: `2.44–2.56 ppm`
  - D2O: `4.68–4.90 ppm`
  - Water: `4.60–5.00 ppm`
- Estimate a single global shift per pair by derivative cross-correlation within `±0.035 ppm`.

Use the same preprocessing for every spectrum. Do not tune parameters per spectrum.

## Ranges

Default ranges:

| Range | ppm |
|---|---|
| R1 | 0.50–2.00 |
| R2 | 2.00–4.00 |
| R3 | 4.00–6.50 |
| R4 | 6.50–8.50 |
| R5 | 8.50–12.00 |

A range is **informative** only when both spectra have:

- at least 2 detected peaks;
- at least 2% of the total masked spectral energy in that range;
- a range maximum of at least 1% of the global maximum.

R5 is often uninformative and will not contribute to a decision. This is intentional.

## Range score

For each informative range:

```text
range_score = 0.5 * weighted_peak_coverage + 0.5 * local_spectral_cosine
```

`weighted_peak_coverage` is one-to-many matching at `0.020 ppm` tolerance, weighted by local peak area. One-to-many matching is used because one spectrum may resolve more multiplet fine structure than another.

`local_spectral_cosine` is computed on the range-local normalized positive spectrum.

## Decision rule

A pair must have at least **two informative ranges** before it can be labeled a duplicate.

### HIGH_CONFIDENCE_DUPLICATE

All conditions:

- at least 2 informative ranges;
- at least 2 ranges with `range_score >= 0.94`;
- weighted range score `>= 0.90`;
- no informative range below `0.65`;
- global full-spectrum cosine `>= 0.97`.

### SUSPECTED_DUPLICATE

All conditions:

- at least 2 informative ranges;
- at least 2 ranges with `range_score >= 0.85`;
- weighted range score `>= 0.85`;
- at most one range below `0.65`;
- global full-spectrum cosine `>= 0.85`.

### Other verdicts

- `NOT_COMPARABLE_DIFFERENT_SOLVENT`: different solvents; do not interpret as duplicate/non-duplicate.
- `INSUFFICIENT_INFORMATIVE_RANGES`: fewer than two informative ranges.
- `NOT_DUPLICATE`: did not meet the above thresholds.

## Scope reporting

When multiple source documents are analyzed, always report:

- internal pairs within each source;
- cross-source pairs;
- excluded datasets;
- non-comparable pairs.

A high-confidence pair inside one source is not evidence of a cross-source duplicate.

## Tuning

The defaults are calibrated for routine small-molecule 1H NMR at approximately 400 MHz. For different fields or highly congested spectra, adjust only after reviewing known positive and negative pairs:

- `--peak-tol`: peak matching tolerance, default `0.020 ppm`.
- `--max-shift`: maximum global alignment, default `0.035 ppm`.
- `--strong-score`: strong range threshold, default `0.94`.
- `--loose-score`: loose range threshold, default `0.85`.
- `--high-global-cosine`: high-confidence guard, default `0.97`.
- `--suspicious-global-cosine`: suspected guard, default `0.85`.
- `--ranges`: custom chemical-shift windows.

Do not lower thresholds simply to increase the number of hits. Validate changes against known duplicates and known non-duplicates.

## Limitations

- Different solvents, concentrations, pH, temperature, and field strengths change spectra.
- Structural isomers may have similar 1H spectra.
- Mixtures and impurities can change local integrals.
- A 1H match does not prove compound identity.
- The raw-FID fallback is less reliable than the instrument-processed `1r`.
- The `.mnova` reader is a path-recovery helper, not a full Mnova parser; verify the extracted inventory.

## Recommended interpretation

- **HIGH_CONFIDENCE_DUPLICATE**: same 1H spectral pattern; review if the scientific question requires compound identity.
- **SUSPECTED_DUPLICATE**: related or similar spectra; preserve the range evidence and do not overclaim.
- **NOT_COMPARABLE**: metadata or acquisition differences prevent a valid comparison.

## Incremental history

Use the same `--history-dir` on every run for a project. The script stores a cumulative `nmr_inventory_history.json` keyed by resolved `fid` path. A new run merges its inventory with the history, computes `NEW_NEW` and `NEW_HISTORY` pairs, and reuses cached `HISTORY_HISTORY` results when dataset hashes and parameters are unchanged.

Each pair is labelled:

- `NEW_NEW`: both datasets came from the current input;
- `NEW_HISTORY`: one dataset is new and the other is from the history;
- `HISTORY_HISTORY`: both datasets came from earlier runs.

Use `nmr_scope_breakdown.csv` to verify that the new file was actually compared with the full history.
