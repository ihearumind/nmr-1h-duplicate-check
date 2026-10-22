# nmr-1h-duplicate-check

A Codex/Agent Skill for batch duplicate screening of Bruker 1H NMR spectra.

## Install

Copy the `nmr-1h-duplicate-check` folder into your Codex skills directory:

```text
$CODEX_HOME/skills/nmr-1h-duplicate-check
```

or, on a default Windows install:

```text
C:/Users/<you>/.codex/skills/nmr-1h-duplicate-check
```

Restart Codex after installation.

## Dependencies

```bash
python -m pip install numpy scipy openpyxl
```

`openpyxl` is optional for CSV-only workflows but required for XLSX output.

## Quick start

From a directory containing Bruker datasets:

```bash
python scripts/nmr_range_duplicate.py --root "F:/path/to/NMR" --out nmr_results --xlsx
```

From one or more `.mnova` files:

```bash
python scripts/extract_mnova_paths.py   --mnova "Document 3.mnova" "Document 4.mnova" "zfs.mnova"   --search-root "F:/path/to/NMR"   --out inventory.json

python scripts/nmr_range_duplicate.py   --inventory inventory.json   --out nmr_results --xlsx
```

## What it outputs

- high-confidence duplicate pairs
- suspected duplicate pairs
- all pairwise comparisons
- range-level evidence
- dataset QC and hashes
- excluded/non-comparable datasets
- optional XLSX workbook

## Important note

High 1H spectral similarity is not proof of compound identity. Use 13C/2D NMR, MS, and sample records when identity matters.
