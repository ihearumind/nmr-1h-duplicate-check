#!/usr/bin/env python3
"""Range-based duplicate detection for Bruker 1H NMR spectra.

Reads Bruker dataset directories or an inventory JSON, compares every 1H
spectrum against every other 1H spectrum in fixed chemical-shift ranges, and
writes CSV reports. If openpyxl is installed, also writes an XLSX workbook.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from collections import Counter

import numpy as np

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
except Exception as exc:
    raise SystemExit("This script requires numpy and scipy: " + str(exc))

STEP = 0.002
GRID = np.arange(0.0, 12.0001, STEP)
DEFAULT_RANGES = [
    ("R1 0.50-2.00", 0.50, 2.00),
    ("R2 2.00-4.00", 2.00, 4.00),
    ("R3 4.00-6.50", 4.00, 6.50),
    ("R4 6.50-8.50", 6.50, 8.50),
    ("R5 8.50-12.00", 8.50, 12.00),
]


def read_params(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    text = path.read_text(errors="ignore")
    for line in text.splitlines():
        m = re.match(r"##\$([^=]+)=\s*(.*)$", line)
        if m:
            k = m.group(1).strip()
            v = m.group(2).strip()
            if v.startswith("<") and v.endswith(">"):
                v = v[1:-1]
            out[k] = v
    return out


def fnum(p: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(str(p.get(key, default)).strip())
    except Exception:
        return float(default)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def dataset_name(ds: Path, root: Path | None) -> str:
    if ds.name.isdigit() and ds.parent.name:
        return f"{ds.parent.name} exp{ds.name}"
    return ds.name


def discover_datasets(roots: list[Path]) -> list[dict]:
    found = []
    seen = set()
    for root in roots:
        root = root.resolve()
        for acqus in root.rglob("acqus"):
            ds = acqus.parent
            if (ds / "fid").exists() or (ds / "pdata" / "1" / "1r").exists():
                key = str(ds.resolve()).lower()
                if key in seen:
                    continue
                seen.add(key)
                found.append({
                    "doc": root.name,
                    "name": dataset_name(ds, root),
                    "fid": str(ds / "fid"),
                    "acqus": str(ds / "acqus"),
                    "pdata1r": str(ds / "pdata" / "1" / "1r"),
                    "procs": str(ds / "pdata" / "1" / "procs"),
                })
    return found


def read_inventory(path: Path) -> list[dict]:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if isinstance(obj, dict) and "datasets" in obj:
        obj = obj["datasets"]
    if not isinstance(obj, list):
        raise SystemExit("Inventory JSON must be a list or {'datasets': [...]}.")
    return obj


def read_1r(ds: Path):
    procs = read_params(ds / "pdata" / "1" / "procs")
    p = ds / "pdata" / "1" / "1r"
    if not p.exists():
        return None
    endian = "<" if int(fnum(procs, "BYTORDP", 0)) == 0 else ">"
    y = np.fromfile(p, dtype=np.dtype(endian + "i4")).astype(float)
    si = int(fnum(procs, "SI", len(y)))
    if si > 0:
        y = y[:si]
    offset = fnum(procs, "OFFSET", 0.0)
    sw_p = fnum(procs, "SW_p", fnum(read_params(ds / "acqus"), "SW_h", 0.0))
    sf = fnum(procs, "SF", fnum(read_params(ds / "acqus"), "SFO1", 0.0))
    if len(y) < 100 or not sf or not sw_p:
        return None
    ppm = offset - np.arange(len(y)) * (sw_p / sf) / max(len(y) - 1, 1)
    return ppm[::-1], y[::-1], "1r"


def read_fid_magnitude(ds: Path):
    acq = read_params(ds / "acqus")
    p = ds / "fid"
    if not p.exists():
        return None
    td = int(fnum(acq, "TD", 0))
    sw = fnum(acq, "SW_h", 0.0)
    sf = fnum(acq, "SFO1", 0.0)
    o1 = fnum(acq, "O1", 0.0)
    if td <= 0 or sw <= 0 or sf <= 0:
        return None
    b = np.fromfile(p, dtype="<i4")
    n = min(len(b) // 2, td)
    if n < 100:
        return None
    z = b[:2 * n].reshape(n, 2)
    fid = z[:, 0].astype(float) + 1j * z[:, 1].astype(float)
    fid *= np.exp(-np.arange(n) * 2e-5)
    spec = np.abs(np.fft.fftshift(np.fft.fft(fid)))
    carrier = o1 / sf
    ppm = carrier + sw / sf / 2 - np.arange(n) * (sw / sf) / max(n - 1, 1)
    return ppm[::-1], spec[::-1], "fid-magnitude-fallback"


def clean_spectrum(ds: Path):
    r = read_1r(ds)
    if r is None:
        r = read_fid_magnitude(ds)
    if r is None:
        return None
    x, v, source = r
    m = (x >= 0.0) & (x <= 12.0)
    if m.sum() < 100:
        return None
    y = np.interp(GRID, x[m], v[m])
    y = np.maximum(y - np.percentile(y, 5), 0.0)
    acq = read_params(ds / "acqus")
    solvent = acq.get("SOLVENT", "")
    mask = np.ones(GRID.shape, dtype=bool)
    if "CD3OD" in solvent:
        mask &= ~((GRID >= 3.24) & (GRID <= 3.38))
    if "CDCl3" in solvent:
        mask &= ~((GRID >= 7.20) & (GRID <= 7.32))
    if "DMSO" in solvent:
        mask &= ~((GRID >= 2.44) & (GRID <= 2.56))
    if "D2O" in solvent:
        mask &= ~((GRID >= 4.68) & (GRID <= 4.90))
    mask &= ~((GRID >= 4.60) & (GRID <= 5.00))
    return y, mask, acq, source


def cosine(x: np.ndarray, y: np.ndarray) -> float:
    nx = float(np.linalg.norm(x))
    ny = float(np.linalg.norm(y))
    if not nx or not ny:
        return 0.0
    return float(np.dot(x, y) / (nx * ny))


def best_shift(a: np.ndarray, b: np.ndarray, mask: np.ndarray, max_shift: float) -> int:
    from scipy.signal import correlate
    da = gaussian_filter1d(a, 1.0, order=1)
    db = gaussian_filter1d(b, 1.0, order=1)
    c = correlate(da[mask], db[mask], mode="full", method="fft")
    mid = int(mask.sum()) - 1
    maxlag = int(round(max_shift / STEP))
    lo = max(0, mid - maxlag)
    hi = min(len(c), mid + maxlag + 1)
    return lo + int(np.argmax(c[lo:hi])) - mid


def shift_array(b: np.ndarray, shift: int) -> np.ndarray:
    if shift == 0:
        return b.copy()
    out = np.zeros_like(b)
    if shift > 0:
        out[shift:] = b[:-shift]
    else:
        out[:shift] = b[-shift:]
    return out


def peak_table(v: np.ndarray, mask: np.ndarray, max_peaks: int = 100):
    s = gaussian_filter1d(v, 1.4)
    vmax = float(np.max(s[mask])) if np.any(mask) else 0.0
    if vmax <= 0:
        return np.array([], dtype=int), np.array([])
    peaks, props = find_peaks(s, prominence=max(vmax * 0.004, 1e-12), height=max(vmax * 0.004, 1e-12), distance=max(2, int(0.004 / STEP)))
    keep = mask[peaks]
    peaks = peaks[keep]
    prom = props["prominences"][keep]
    if len(peaks) > max_peaks:
        order = np.argsort(prom)[-max_peaks:]
        peaks = peaks[order]
        prom = prom[order]
    areas = []
    for i in peaks:
        lo = max(0, i - int(0.012 / STEP))
        hi = min(len(v), i + int(0.012 / STEP) + 1)
        areas.append(float(np.trapz(np.maximum(v[lo:hi], 0.0), dx=STEP)))
    areas = np.asarray(areas)
    if len(areas) and areas.sum() > 0:
        weights = areas / areas.sum()
    elif len(prom) and prom.sum() > 0:
        weights = prom / prom.sum()
    else:
        weights = np.array([])
    return np.asarray(peaks), weights


def peak_coverage(p1, w1, p2, w2, tol_ppm: float = 0.020) -> float:
    if len(p1) == 0 or len(p2) == 0:
        return 0.0
    a = sum(float(w) for p, w in zip(p1, w1) if np.any(np.abs(GRID[p2] - GRID[p]) <= tol_ppm))
    b = sum(float(w) for p, w in zip(p2, w2) if np.any(np.abs(GRID[p1] - GRID[p]) <= tol_ppm))
    return 0.5 * (a + b)


def binned_cosine(x: np.ndarray, y: np.ndarray, bin_points: int) -> float:
    n = (len(x) // bin_points) * bin_points
    xx = x[:n].reshape(-1, bin_points).sum(axis=1)
    yy = y[:n].reshape(-1, bin_points).sum(axis=1)
    return cosine(xx, yy)


def parse_ranges(text: str | None):
    if not text:
        return DEFAULT_RANGES
    ranges = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^([0-9.]+)\s*-\s*([0-9.]+)$", part)
        if not m:
            raise SystemExit(f"Bad range '{part}'. Use e.g. 0.5-2,2-4.")
        lo, hi = float(m.group(1)), float(m.group(2))
        ranges.append((f"R{len(ranges)+1} {lo:.2f}-{hi:.2f}", lo, hi))
    return ranges


def build_xlsx(path: Path, rows, evidence, qc, high, susp, ranges, args):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
        from openpyxl.utils import get_column_letter
    except Exception:
        print("openpyxl not available; skipped XLSX.", file=sys.stderr)
        return
    hf = PatternFill("solid", fgColor="1F4E78")
    hfont = Font(color="FFFFFF", bold=True)
    green = PatternFill("solid", fgColor="C6EFCE")
    yellow = PatternFill("solid", fgColor="FFEB9C")
    grey = PatternFill("solid", fgColor="E7E6E6")
    thin = Side(style="thin", color="A6A6A6")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    def table(ws, headers, data, widths=None):
        ws.append(headers)
        for c in ws[1]:
            c.fill = hf; c.font = hfont; c.alignment = Alignment(horizontal="center", wrap_text=True)
        for r in data:
            ws.append([r.get(h, "") for h in headers])
        ws.freeze_panes = "A2"; ws.auto_filter.ref = ws.dimensions
        for i, h in enumerate(headers, 1):
            ws.column_dimensions[get_column_letter(i)].width = (widths or {}).get(h, 16)
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.border = border; c.alignment = Alignment(vertical="top", wrap_text=True)

    wb = Workbook(); ws = wb.active; ws.title = "Summary"
    ws["A1"] = "NMR 1H Range-Based Duplicate Report"; ws["A1"].font = Font(bold=True, size=16)
    ws["A2"] = "Method: compare peaks and local spectra in chemical-shift ranges; require at least two informative ranges."
    ws["A3"] = f"Datasets: {len(qc)} 1H; pairwise comparisons: {len(rows)}; high: {len(high)}; suspected: {len(susp)}"
    ws["A5"] = "High-confidence duplicate pairs"
    for i, r in enumerate(high, 1):
        ws.cell(6 + i, 1, f"HD-{i:02d}"); ws.cell(6 + i, 2, r["name_a"]); ws.cell(6 + i, 3, r["name_b"]); ws.cell(6 + i, 4, float(r["weighted_range_score"]))
    ws.column_dimensions["A"].width = 12; ws.column_dimensions["B"].width = 34; ws.column_dimensions["C"].width = 34; ws.column_dimensions["D"].width = 18

    ws = wb.create_sheet("High_Confidence_Pairs")
    table(ws, ["Pair ID","File A","File B","Document A","Document B","Weighted range score","Strong ranges","Global cosine","Range summary"], [{
        "Pair ID": f"HD-{i:02d}", "File A": r["name_a"], "File B": r["name_b"], "Document A": r["doc_a"], "Document B": r["doc_b"],
        "Weighted range score": float(r["weighted_range_score"]), "Strong ranges": int(float(r["strong_ranges"])),
        "Global cosine": float(r["global_cosine"]), "Range summary": r["range_summary"]} for i, r in enumerate(high, 1)], {"File A":34,"File B":34,"Range summary":90})
    for row in ws.iter_rows(min_row=2):
        for c in row: c.fill = green

    ws = wb.create_sheet("Suspected_Pairs")
    table(ws, ["File A","File B","Weighted range score","Informative ranges","Strong ranges","Loose ranges","Global cosine","Verdict","Range summary"], [{
        "File A": r["name_a"], "File B": r["name_b"], "Weighted range score": float(r["weighted_range_score"]),
        "Informative ranges": int(float(r["informative_ranges"])), "Strong ranges": int(float(r["strong_ranges"])),
        "Loose ranges": int(float(r["loose_ranges"])), "Global cosine": float(r["global_cosine"]),
        "Verdict": r["verdict"], "Range summary": r["range_summary"]} for r in susp], {"File A":34,"File B":34,"Range summary":90})
    for row in ws.iter_rows(min_row=2):
        for c in row: c.fill = yellow

    ws = wb.create_sheet("Dataset_QC")
    table(ws, ["doc","name","nucleus","solvent","pulprog","sfo1_mhz","source","parse_ok","fid_sha256","pdata_sha256","dataset"], qc, {"name":36,"dataset":90,"fid_sha256":66,"pdata_sha256":66})

    ws = wb.create_sheet("All_Pairs")
    table(ws, list(rows[0].keys()) if rows else ["status"], rows, {"name_a":34,"name_b":34,"range_summary":90,"verdict":30})
    for row in ws.iter_rows(min_row=2):
        v = row[-2].value
        if v == "HIGH_CONFIDENCE_DUPLICATE":
            for c in row: c.fill = green
        elif v == "SUSPECTED_DUPLICATE":
            for c in row: c.fill = yellow
        elif str(v).startswith("NOT_COMPARABLE"):
            for c in row: c.fill = grey

    ws = wb.create_sheet("Method")
    lines = [
        "Ranges: " + ", ".join(f"{n} ({lo}-{hi} ppm)" for n, lo, hi in ranges),
        "A range is informative only when both spectra have >=2 peaks, >=min_energy fraction and a range maximum >=min_maxfrac of the global maximum.",
        "Range score = 0.5 * weighted peak coverage + 0.5 * local spectral cosine.",
        "HIGH confidence: >=2 informative ranges, >=2 ranges with score >= strong_score, weighted score >=0.90, no range below 0.65, global cosine >= high_global_cosine.",
        "SUSPECTED: >=2 informative ranges, >=2 ranges with score >= loose_score, weighted score >=0.85, at most one range below 0.65, global cosine >= suspicious_global_cosine.",
        "Different solvents are not comparable. A 1H spectral match does not prove compound identity.",
    ]
    for i, text in enumerate(lines, 1):
        ws.cell(i, 1, text).alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 140
    wb.save(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Range-based duplicate detection for Bruker 1H NMR spectra")
    ap.add_argument("--root", action="append", default=[], help="Root directory containing Bruker datasets. Repeat for multiple roots.")
    ap.add_argument("--inventory", help="JSON inventory with name/fid/acqus/pdata1r/procs/doc entries.")
    ap.add_argument("--history-dir", help="Persistent directory for the cumulative inventory. New runs merge with all previous datasets and compare all of them.")
    ap.add_argument("--out", required=True, help="Output directory.")
    ap.add_argument("--ranges", help="Comma-separated ranges, e.g. 0.5-2,2-4,4-6.5,6.5-8.5,8.5-12")
    ap.add_argument("--min-energy", type=float, default=0.02)
    ap.add_argument("--min-maxfrac", type=float, default=0.01)
    ap.add_argument("--peak-tol", type=float, default=0.020)
    ap.add_argument("--max-shift", type=float, default=0.035)
    ap.add_argument("--strong-score", type=float, default=0.94)
    ap.add_argument("--loose-score", type=float, default=0.85)
    ap.add_argument("--high-global-cosine", type=float, default=0.97)
    ap.add_argument("--suspicious-global-cosine", type=float, default=0.85)
    ap.add_argument("--xlsx", action="store_true", help="Write XLSX if openpyxl is available.")
    args = ap.parse_args(argv)

    new_inventory = read_inventory(Path(args.inventory)) if args.inventory else (discover_datasets([Path(p) for p in args.root]) if args.root else None)
    if not new_inventory:
        raise SystemExit("Provide --root or --inventory and ensure datasets exist.")
    new_keys = {str(Path(x.get("fid", "")).resolve()).lower() for x in new_inventory if x.get("fid")}
    inventory = new_inventory
    if args.history_dir:
        history_dir = Path(args.history_dir)
        history_dir.mkdir(parents=True, exist_ok=True)
        history_file = history_dir / "nmr_inventory_history.json"
        previous = []
        if history_file.exists():
            previous = json.loads(history_file.read_text(encoding="utf-8-sig"))
        merged = {}
        for item in previous:
            if item.get("fid"):
                merged[str(Path(item["fid"]).resolve()).lower()] = item
        for item in new_inventory:
            if item.get("fid"):
                merged[str(Path(item["fid"]).resolve()).lower()] = item
        inventory = list(merged.values())
        history_file.write_text(json.dumps(inventory, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"History inventory: {history_file.resolve()} ({len(inventory)} datasets)")
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True); ranges = parse_ranges(args.ranges)
    datasets = []; excluded = []
    for item in inventory:
        ds = Path(item["fid"]).parent; acq = read_params(ds / "acqus")
        rec = {"doc": item.get("doc",""), "name": item.get("name", dataset_name(ds,None)), "fid": str(item.get("fid",ds/"fid")), "acqus": str(item.get("acqus",ds/"acqus")), "pdata1r": str(item.get("pdata1r",ds/"pdata"/"1"/"1r")), "procs": str(item.get("procs",ds/"pdata"/"1"/"procs")), "nucleus": acq.get("NUC1",""), "solvent": acq.get("SOLVENT",""), "pulprog": acq.get("PULPROG",""), "sfo1_mhz": fnum(acq,"SFO1",float("nan")), "fid_sha256": sha256(Path(item["fid"])) if Path(item["fid"]).exists() else "", "pdata_sha256": sha256(Path(item["pdata1r"])) if Path(item["pdata1r"]).exists() else "", "dataset": str(ds), "is_new": str(Path(item["fid"]).resolve()).lower() in new_keys}
        if rec["nucleus"] != "1H":
            excluded.append(rec); continue
        try:
            c = clean_spectrum(ds)
            if c is None:
                rec["parse_ok"] = False; rec["source"] = ""; excluded.append(rec); continue
            rec["parse_ok"] = True; rec["source"] = c[3]; datasets.append((rec,c))
        except Exception as exc:
            rec["parse_ok"] = False; rec["source"] = repr(exc); excluded.append(rec)
    if len(datasets) < 2: raise SystemExit("Fewer than two 1H datasets were parsed.")
    data = {rec["name"]: (y,mask,acq,source) for rec,(y,mask,acq,source) in datasets}
    rec_by_name = {rec["name"]: rec for rec, _ in datasets}
    names = list(data.keys()); feats = {}
    for name,(v,mask,acq,source) in data.items():
        feats[name] = {}; vv=np.maximum(v,0.0); total=float(vv[mask].sum()); gmax=float(vv[mask].max()) if np.any(mask) else 0.0
        for rn,lo,hi in ranges:
            sel=mask&(GRID>=lo)&(GRID<=hi)
            if not np.any(sel):
                feats[name][rn]={"info":False,"p":np.array([],int),"w":np.array([]),"energy":0.0}; continue
            p,w=peak_table(v,sel); energy=float(vv[sel].sum()/max(total,1e-12)); maxfrac=float(vv[sel].max()/max(gmax,1e-12))
            feats[name][rn]={"info":energy>=args.min_energy and maxfrac>=args.min_maxfrac and len(p)>=2,"p":p,"w":w,"energy":energy}
    history_cache = {}
    history_cache_file = None
    cache_config_signature = None
    if args.history_dir:
        history_dir = Path(args.history_dir)
        history_dir.mkdir(parents=True, exist_ok=True)
        history_cache_file = history_dir / "nmr_pair_history.json"
        config_obj = {
            "ranges": ranges, "min_energy": args.min_energy, "min_maxfrac": args.min_maxfrac,
            "peak_tol": args.peak_tol, "max_shift": args.max_shift, "strong_score": args.strong_score,
            "loose_score": args.loose_score, "high_global_cosine": args.high_global_cosine,
            "suspicious_global_cosine": args.suspicious_global_cosine,
        }
        cache_config_signature = hashlib.sha256(json.dumps(config_obj, sort_keys=True).encode("utf-8")).hexdigest()
        if history_cache_file.exists():
            try:
                cached = json.loads(history_cache_file.read_text(encoding="utf-8"))
                if cached.get("config_signature") == cache_config_signature:
                    history_cache = {tuple(sorted([item["row"]["name_a"], item["row"]["name_b"]])): item for item in cached.get("pairs", [])}
                    print(f"Loaded cached pair results: {len(history_cache)}")
                else:
                    print("Comparison parameters changed; cached pairs will not be reused.")
            except Exception as exc:
                print(f"Ignoring invalid pair cache: {exc}", file=sys.stderr)

    def pair_signature(ra, rb):
        payload = [cache_config_signature, sorted([ra.get("fid_sha256", ""), rb.get("fid_sha256", "")]), sorted([ra.get("pdata_sha256", ""), rb.get("pdata_sha256", "")])]
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()

    rows=[]; evidence=[]; cache_entries=[]
    for i,na in enumerate(names):
        for nb in names[i+1:]:
            va,ma,aca,_=data[na]; vb,mb,acb,_=data[nb]; mask=ma&mb
            ra, rb = rec_by_name[na], rec_by_name[nb]
            if ra["is_new"] and rb["is_new"]:
                comparison_type = "NEW_NEW"
            elif ra["is_new"] or rb["is_new"]:
                comparison_type = "NEW_HISTORY"
            else:
                comparison_type = "HISTORY_HISTORY"
            pkey = tuple(sorted([na, nb]))
            sig = pair_signature(ra, rb) if cache_config_signature else ""
            cached_entry = history_cache.get(pkey)
            if comparison_type == "HISTORY_HISTORY" and cached_entry and cached_entry.get("signature") == sig:
                rows.append(cached_entry["row"])
                cache_entries.append(cached_entry)
                continue
            shift=best_shift(va,vb,mask,args.max_shift); aa=va.copy(); bb=shift_array(vb,shift)
            same_solvent=aca.get("SOLVENT","")==acb.get("SOLVENT",""); range_rows=[]
            for rn,lo,hi in ranges:
                fa,fb=feats[na][rn],feats[nb][rn]; sel=mask&(GRID>=lo)&(GRID<=hi)
                if fa["info"] and fb["info"] and np.any(sel):
                    x=np.maximum(aa[sel],0); y=np.maximum(bb[sel],0); cov=peak_coverage(fa["p"],fa["w"],fb["p"],fb["w"],args.peak_tol); score=.5*(cov+cosine(x,y))
                    range_rows.append({"range":rn,"score":score,"coverage":cov,"cosine":cosine(x,y),"energy_a":fa["energy"],"energy_b":fb["energy"],"peaks_a":len(fa["p"]),"peaks_b":len(fb["p"]),"weight":max(.02,min(fa["energy"],fb["energy"]))})
            scores=[x["score"] for x in range_rows]; weights=[x["weight"] for x in range_rows]; weighted=float(np.average(scores,weights=weights)) if scores else 0.0
            n_inf=len(scores); n_strong=sum(1 for x in scores if x>=args.strong_score); n_loose=sum(1 for x in scores if x>=args.loose_score); n_bad=sum(1 for x in scores if x<.65)
            x=np.maximum(aa[mask],0); y=np.maximum(bb[mask],0); gcos=cosine(x,y); gbin=binned_cosine(x,y,10); gder=cosine(gaussian_filter1d(x,1,order=1),gaussian_filter1d(y,1,order=1)); gpa,gwa=peak_table(aa,mask); gpb,gwb=peak_table(bb,mask); gcov=peak_coverage(gpa,gwa,gpb,gwb,args.peak_tol)
            if not same_solvent: verdict="NOT_COMPARABLE_DIFFERENT_SOLVENT"
            elif n_inf>=2 and n_strong>=2 and weighted>=.90 and n_bad==0 and gcos>=args.high_global_cosine: verdict="HIGH_CONFIDENCE_DUPLICATE"
            elif n_inf>=2 and n_loose>=2 and weighted>=.85 and n_bad<=1 and gcos>=args.suspicious_global_cosine: verdict="SUSPECTED_DUPLICATE"
            elif n_inf<2: verdict="INSUFFICIENT_INFORMATIVE_RANGES"
            else: verdict="NOT_DUPLICATE"
            row = {"doc_a":ra["doc"],"doc_b":rb["doc"],"name_a":na,"name_b":nb,"is_new_a":ra["is_new"],"is_new_b":rb["is_new"],"comparison_type":comparison_type,"same_solvent":same_solvent,"shift_ppm":shift*STEP,"informative_ranges":n_inf,"strong_ranges":n_strong,"loose_ranges":n_loose,"bad_ranges":n_bad,"weighted_range_score":weighted,"global_cosine":gcos,"global_bin_cos_020":gbin,"global_deriv_cos":gder,"global_peak_coverage_020":gcov,"global_unique_ratio":1-gcov,"verdict":verdict,"range_summary":"; ".join(f"{x['range']}={x['score']:.3f}" for x in range_rows)}
            rows.append(row)
            for rr in range_rows: evidence.append({"name_a":na,"name_b":nb,**rr})
            if cache_config_signature:
                cache_entries.append({"signature": sig, "row": row})

    if history_cache_file is not None and cache_config_signature is not None:
        history_cache_file.write_text(json.dumps({"config_signature": cache_config_signature, "pairs": cache_entries}, ensure_ascii=False, indent=2), encoding="utf-8")
    rows.sort(key=lambda r:(-float(r["weighted_range_score"]),-int(r["strong_ranges"]))); high=[r for r in rows if r["verdict"]=="HIGH_CONFIDENCE_DUPLICATE"]; susp=[r for r in rows if r["verdict"]=="SUSPECTED_DUPLICATE"]; qc=[r for r,_ in datasets]
    scope_counts = Counter()
    for r in rows:
        scope = (f"{r['doc_a']} internal" if r['doc_a'] == r['doc_b'] else " + ".join(sorted([r['doc_a'], r['doc_b']])))
        scope_counts[(scope, r['comparison_type'])] += 1
    scope_rows = []
    for (scope, ctype), n in sorted(scope_counts.items()):
        subset = [r for r in rows if (r['doc_a'] == r['doc_b'] and f"{r['doc_a']} internal" == scope or r['doc_a'] != r['doc_b'] and " + ".join(sorted([r['doc_a'], r['doc_b']])) == scope) and r['comparison_type'] == ctype]
        h = sum(1 for r in subset if r['verdict'] == 'HIGH_CONFIDENCE_DUPLICATE')
        su = sum(1 for r in subset if r['verdict'] == 'SUSPECTED_DUPLICATE')
        scope_rows.append({'scope': scope, 'comparison_type': ctype, 'pair_count': n, 'high_confidence': h, 'suspected': su, 'other': n - h - su})
    def write_csv(name,fields,data):
        with open(out/name,"w",newline="",encoding="utf-8-sig") as f:
            w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore"); w.writeheader(); w.writerows(data)
    write_csv("nmr_duplicate_pairs.csv",list(rows[0].keys()),high); write_csv("nmr_suspected_pairs.csv",list(rows[0].keys()),susp); write_csv("nmr_all_pairs.csv",list(rows[0].keys()),rows); write_csv("nmr_dataset_qc.csv",["doc","name","nucleus","solvent","pulprog","sfo1_mhz","source","parse_ok","fid_sha256","pdata_sha256","dataset"],qc); write_csv("nmr_range_evidence.csv",["name_a","name_b","range","score","coverage","cosine","energy_a","energy_b","peaks_a","peaks_b","weight"],evidence)
    write_csv("nmr_scope_breakdown.csv", ["scope","comparison_type","pair_count","high_confidence","suspected","other"], scope_rows)
    with open(out/"nmr_excluded.csv","w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=["doc","name","nucleus","solvent","source","parse_ok","fid","dataset"],extrasaction="ignore"); w.writeheader(); w.writerows(excluded)
    if args.xlsx: build_xlsx(out/"nmr_duplicate_report.xlsx",rows,evidence,qc,high,susp,ranges,args)
    print(f"Datasets parsed: {len(datasets)}; pairs: {len(rows)}; high: {len(high)}; suspected: {len(susp)}"); print(f"Outputs written to: {out.resolve()}")

if __name__ == "__main__": main()
