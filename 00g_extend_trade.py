#!/usr/bin/env python3
"""
00g_extend_trade.py — fill missing months for core trade series via Comtrade
preview (curl-equivalent urllib). Does NOT overwrite existing months.
Merges into data/raw/{aus_exports,imports_jpn,imports_kor,imports_ind}.csv
"""
from __future__ import annotations
import csv, json, time, sys, urllib.request
from pathlib import Path

RAW = Path(__file__).resolve().parent / "data" / "raw"
CMD = "2701"
TARGETS = [
    ("aus_exports.csv", "36", "X"),
    ("imports_jpn.csv", "392", "M"),
    ("imports_kor.csv", "410", "M"),
    ("imports_ind.csv", "699", "M"),
]


def fetch(reporter, flow, period, tries=4):
    # Do NOT pass customsCode/motCode — preview often returns empty with them
    # for older years; without them we get the aggregate World row.
    url = (
        "https://comtradeapi.un.org/public/v1/preview/C/M/HS?"
        f"reporterCode={reporter}&partnerCode=0&partner2Code=0&"
        f"period={period}&flowCode={flow}&cmdCode={CMD}"
    )
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "coalEl/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                j = json.loads(r.read().decode())
            tons = 0.0
            for rec in j.get("data") or []:
                cc = str(rec.get("customsCode") or "C00")
                if cc not in ("C00", ""):
                    continue
                w = rec.get("netWgt")
                if w is None:
                    continue
                try:
                    wv = float(w)
                except (TypeError, ValueError):
                    continue
                if wv > 0:
                    tons += wv / 1000.0  # kg → tonnes
            return tons if tons > 0 else None
        except Exception as e:
            last = e
            time.sleep(2 * (k + 1))
    print(f"  FAIL {period}: {last}")
    return None


def load_existing(path):
    have = {}
    if path.exists():
        with path.open() as f:
            for r in csv.DictReader(f):
                have[r["date"][:7]] = float(r["tonnes"])
    return have


def save(path, have):
    rows = [{"date": d, "tonnes": have[d]} for d in sorted(have)]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["date", "tonnes"])
        w.writeheader()
        w.writerows(rows)


def main():
    y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 2008
    y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2011
    for fname, code, flow in TARGETS:
        path = RAW / fname
        have = load_existing(path)
        print(f"\n{fname} existing={len(have)}")
        n_new = 0
        for y in range(y0, y1 + 1):
            for m in range(1, 13):
                ds = f"{y}-{m:02d}"
                if ds in have:
                    continue
                period = f"{y}{m:02d}"
                tons = fetch(code, flow, period)
                time.sleep(0.7)
                if tons is None:
                    print(f"  {ds}  -- empty")
                    continue
                have[ds] = tons
                n_new += 1
                print(f"  {ds}  {tons/1e6:.2f} Mt")
        save(path, have)
        # annual QA
        by = {}
        for d, t in have.items():
            by.setdefault(d[:4], 0.0)
            by[d[:4]] += t
        print(f"  -> +{n_new} months; annual Mt:")
        for y in sorted(by):
            if int(y) >= y0 - 1:
                nmo = sum(1 for d in have if d.startswith(y))
                print(f"     {y}  {by[y]/1e6:6.1f} Mt  ({nmo} mo)")


if __name__ == "__main__":
    main()
