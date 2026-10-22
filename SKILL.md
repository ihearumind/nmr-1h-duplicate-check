---
name: nmr-1h-duplicate-check
description: Analyze batches of Bruker 1H NMR spectra or MestReNova-referenced data for duplicate and near-duplicate spectra using range-wise peak comparison. Use when asked to batch-check 1H NMR spectra, compare spectra across Mnova documents, or produce a list of duplicate pairs; do not use for 13C-only/2D-only datasets or de novo structure elucidation.
---

# NMR 1H Duplicate Check

## Outcome

Produce an auditable batch report that answers: **which 1H spectra are duplicate or near-duplicate, within each source and across sources?** Report high-confidence duplicates, suspected duplicates, excluded/non-comparable data, pairwise metrics, and range-level evidence.

This skill detects **spectral duplication**, not compound identity. Never claim that a 1H match alone proves two samples are the same compound.

## Inputs

Preferred input is a set of Bruker dataset directories containing:

```text
<sample>/<expno>/
|-- fid
|-- acqus
`-- pdata/1/
    |-- 1r
    `-- procs
```

`pdata/1/1r` is the primary comparison trace. If it is absent, the helper falls back to a raw-FID magnitude spectrum and records that in QC; treat that fallback as lower-confidence evidence.

For `.mnova` input, extract referenced Bruker paths first. The `.mnova` container is proprietary, so verify the extracted inventory before analysis.

## Workflow

1. **Extract or enumerate datasets**
   - For one or more `.mnova` files:

     ```bash
     python scripts/extract_mnova_paths.py \
       --mnova "Document 3.mnova" "Document 4.mnova" "zfs.mnova" \
       --search-root "F:/课题步骤" \
       --out inventory.json
     ```

   - Or let the comparison script discover Bruker directories:

     ```bash
     python scripts/nmr_range_duplicate.py \
       --root "F:/path/to/NMR" \
       --out nmr_results --xlsx
     ```

2. **Check inventory and QC**
   - Verify every intended source is represented.
   - Filter to `NUC1=1H`; exclude 13C, DEPT, and 2D experiments from the 1H comparison.
   - Keep unparsed or excluded files in a separate table. Do not silently drop them.
   - Record nucleus, solvent, field, pulse program, source, and SHA-256 of `fid` and `pdata/1/1r`.

3. **Compare all pairs**
   - Compare every 1H spectrum against every other 1H spectrum, including within each source and across sources.
   - Mark pairs with different solvents as `NOT_COMPARABLE_DIFFERENT_SOLVENT`; do not turn a cross-solvent score into a duplicate verdict.
   - Use the same preprocessing and thresholds for every spectrum.

4. **Apply range-wise decisions**
   - Ranges: `0.50-2.00`, `2.00-4.00`, `4.00-6.50`, `6.50-8.50`, and `8.50-12.00` ppm.
   - Exclude common solvent/water regions.
   - A range is informative only when both spectra have at least two peaks, at least 2% of total spectral energy, and a range maximum of at least 1% of the global maximum.
   - Range score = `0.5 * weighted 0.020-ppm peak coverage + 0.5 * local spectral cosine`.
   - Require at least two informative ranges before making a duplicate decision.
   - HIGH confidence requires at least two ranges scoring >= 0.94, weighted range score >= 0.90, no range below 0.65, and global cosine >= 0.97.
   - SUSPECTED requires at least two ranges scoring >= 0.85, weighted range score >= 0.85, at most one range below 0.65, and global cosine >= 0.85.
   - See `references/method.md` for formulas and tuning guidance.

5. **Report**
   - Always include: high-confidence pairs, suspected pairs, source/scope breakdown, dataset QC, and range evidence.
   - For a batch containing multiple files, separately report:
     - internal pairs within each source;
     - cross-source pairs.
   - State separately when a pair is `NOT_COMPARABLE` because of different solvent or nucleus.
   - Use `references/report_schema.md` for output columns and sheet definitions.

## Commands

```bash
# CSV only
python scripts/nmr_range_duplicate.py \
  --inventory inventory.json \
  --out nmr_results

# CSV + XLSX
python scripts/nmr_range_duplicate.py \
  --inventory inventory.json \
  --out nmr_results --xlsx

# Custom ranges or thresholds
python scripts/nmr_range_duplicate.py \
  --root "F:/path/to/NMR" \
  --out nmr_results \
  --ranges "0.5-2,2-4,4-6.5,6.5-8.5,8.5-12" \
  --peak-tol 0.020 \
  --high-global-cosine 0.97 \
  --xlsx
```

Dependencies: `numpy`, `scipy`; `openpyxl` is optional but required for XLSX output.

## Reliability Rules

- Prefer the instrument-processed `pdata/1/1r` over a raw-FID fallback.
- Never compare 13C/DEPT/2D spectra as if they were 1H spectra.
- Do not treat different-solvent pairs as duplicates.
- Do not decide from a single chemical-shift range; require at least two informative ranges.
- Do not hide excluded files, parse failures, or non-comparable pairs.
- Do not tune processing parameters separately for each spectrum.
- Treat labels such as `HIGH_CONFIDENCE_DUPLICATE` as “same 1H spectral pattern”, not “same molecule”.
- For borderline cases, preserve the numerical evidence rather than forcing a binary answer.

## Supporting Files

- `scripts/nmr_range_duplicate.py` — main batch comparator and report generator.
- `scripts/extract_mnova_paths.py` — extract Bruker paths from `.mnova` binary containers.
- `references/method.md` — detailed method, formulas, thresholds, and limitations.
- `references/report_schema.md` — CSV/XLSX schema and verdict definitions.
