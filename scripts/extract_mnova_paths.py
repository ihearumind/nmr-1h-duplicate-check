#!/usr/bin/env python3
"""Extract Bruker dataset references from MestReNova .mnova files.

The .mnova format is proprietary. This helper scans the binary for UTF-16LE
and UTF-16BE path-like strings that end in /fid, then resolves them against
the filesystem or a search root. Always verify the resulting inventory.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def printable(cp: int) -> bool:
    return (
        0x20 <= cp <= 0x7E
        or 0x4E00 <= cp <= 0x9FFF
        or 0x3000 <= cp <= 0x303F
        or 0xFF00 <= cp <= 0xFFEF
        or cp in (9, 10, 13)
    )


def extract_strings(data: bytes, big_endian: bool) -> list[str]:
    strings: list[str] = []
    i = 0
    n = len(data)
    while i + 1 < n:
        cp = ((data[i] << 8) | data[i + 1]) if big_endian else (data[i] | (data[i + 1] << 8))
        if printable(cp):
            chars = []
            while i + 1 < n:
                cp = ((data[i] << 8) | data[i + 1]) if big_endian else (data[i] | (data[i + 1] << 8))
                if not printable(cp):
                    break
                chars.append(chr(cp))
                i += 2
            if len(chars) >= 4:
                strings.append("".join(chars))
        else:
            i += 2
    return strings


def candidate_paths(text: str) -> list[str]:
    out = []
    # Full Windows path, possibly prefixed by a control character.
    for m in re.finditer(r"([A-Za-z]:[/\\].*?/(?:fid|acqus))", text, re.I):
        p = m.group(1).replace("\\", "/")
        if p not in out:
            out.append(p)
    # Relative path /something/.../fid
    for m in re.finditer(r"(/[^\x00\r\n]*?/(?:fid|acqus))", text):
        p = m.group(1).replace("\\", "/")
        if p not in out:
            out.append(p)
    return out


def dataset_name(fid: Path) -> str:
    ds = fid.parent
    if ds.name.isdigit() and ds.parent.name:
        return f"{ds.parent.name} exp{ds.name}"
    return ds.name


def resolve_inventory(strings, mnova: Path, search_root, index, include_acqus_only, seen):
    cands = []
    evidence_dirs = set()
    for text in strings:
        for c in candidate_paths(text):
            if c not in cands:
                cands.append(c)
        for m in re.finditer(r"([A-Za-z]:[/\\].*?/pdata/1/(?:title|procs|1r)|/[^\x00\r\n]*?/pdata/1/(?:title|procs|1r))", text, re.I):
            d = re.sub(r"/pdata/1/(?:title|procs|1r)$", "", m.group(1).replace("\\", "/"), flags=re.I)
            evidence_dirs.add(d.lower())
    resolved_items = []
    for cand in cands:
        if cand.lower().endswith("/acqus") and not include_acqus_only:
            continue
        resolved = None
        reason = ""
        target = cand[:-6] + "/fid" if cand.lower().endswith("/acqus") else cand
        p = Path(target.replace("/", "\\")) if re.match(r"^[A-Za-z]:", target) else None
        if p and p.exists() and p.name.lower() == "fid":
            ds_candidate = p.parent.resolve()
            if cand.lower().endswith("/acqus") and str(ds_candidate).replace("\\", "/").lower() not in evidence_dirs:
                continue
            resolved = p.resolve(); reason = "absolute-path"
        elif search_root:
            needle = target.replace("\\", "/").lower()
            matches = []
            for fid in index:
                f = str(fid).replace("\\", "/").lower()
                if f.endswith(needle) or needle.endswith(f):
                    matches.append(fid)
            if len(matches) == 1:
                ds_candidate = matches[0].parent.resolve()
                if cand.lower().endswith("/acqus") and str(ds_candidate).replace("\\", "/").lower() not in evidence_dirs:
                    continue
                resolved = matches[0]; reason = "search-root-suffix"
        if resolved is None:
            continue
        key = str(resolved).lower()
        if key in seen:
            continue
        seen.add(key)
        ds = resolved.parent
        resolved_items.append({
            "doc": mnova.stem,
            "name": dataset_name(resolved),
            "fid": str(ds / "fid"),
            "acqus": str(ds / "acqus"),
            "pdata1r": str(ds / "pdata" / "1" / "1r"),
            "procs": str(ds / "pdata" / "1" / "procs"),
            "resolved_by": reason,
        })
    return resolved_items

def main():
    ap = argparse.ArgumentParser(description="Extract Bruker dataset references from .mnova files")
    ap.add_argument("--mnova", nargs="+", required=True)
    ap.add_argument("--search-root", help="Root directory containing Bruker datasets, used to resolve relative paths.")
    ap.add_argument("--include-acqus-only", action="store_true", help="Also include datasets referenced only by acqus + pdata/title/procs paths. Higher recall, lower precision; verify the inventory.")
    ap.add_argument("--out", required=True, help="Output inventory JSON.")
    args = ap.parse_args()

    search_root = Path(args.search_root).resolve() if args.search_root else None
    index = {}
    if search_root:
        for fid in search_root.rglob("fid"):
            try:
                index[fid.resolve()] = str(fid.resolve()).replace("\\", "/")
            except Exception:
                pass

    inventory = []
    seen = set()
    for mnova_name in args.mnova:
        mnova = Path(mnova_name)
        data = mnova.read_bytes()
        le_strings = extract_strings(data, False)
        be_strings = extract_strings(data, True)
        # Some .mnova files contain both encodings. Prefer the encoding whose
        # candidate paths can actually be resolved, then fall back to the other.
        le_items = resolve_inventory(le_strings, mnova, search_root, index, args.include_acqus_only, seen)
        be_items = resolve_inventory(be_strings, mnova, search_root, index, args.include_acqus_only, seen) if not le_items else []
        inventory.extend(le_items if le_items else be_items)

    out = Path(args.out)
    out.write_text(json.dumps(inventory, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Found {len(inventory)} Bruker datasets. Inventory: {out.resolve()}")


if __name__ == "__main__":
    main()
