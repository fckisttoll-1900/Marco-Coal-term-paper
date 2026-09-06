#!/usr/bin/env python3
"""
04c_pull_bilateral_curl.py
-------------------------
Monthly bilateral coal imports (HS 2701) via Comtrade preview + curl.
Resume-safe: skips (date, partner) already on disk.

Default: LOW Aus-share clean-margin importers + Taiwan (HIGH border).
Japan/Korea/India already exist — left alone unless listed.

Usage:
  PYTHONUNBUFFERED=1 python3 04c_pull_bilateral_curl.py
  PYTHONUNBUFFERED=1 python3 04c_pull_bilateral_curl.py 2012 2024 deu,tur,phl
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

RAW = Path(__file__).resolve().parent / "data" / "raw"
CMD = "2701"

# stub -> (M49, label)
IMPORTERS = {
    "twn": ("490", "Taiwan"),
    "tur": ("792", "Turkey"),
    "deu": ("276", "Germany"),
    "nld": ("528", "Netherlands"),
    "gbr": ("826", "United Kingdom"),
    "phl": ("608", "Philippines"),
    "mar": ("504", "Morocco"),
    "vnm": ("704", "Vietnam"),
    "bra": ("76", "Brazil"),
    "jpn": ("392", "Japan"),
    "kor": ("410", "Korea"),
    "ind": ("699", "India"),
}

PARTNERS = {
    "0": "World",
    "36": "Australia",
    "360": "Indonesia",
    "643": "Russia",
    "842": "USA",
    "124": "Canada",
    "710": "South Africa",
    "170": "Colombia",
    "496": "Mongolia",
}
PARTNER_CODES = ",".join(PARTNERS.keys())


def curl_json(url: str, timeout: int = 60):
    tmp = Path("/tmp") / f"ct_{abs(hash(url)) % 10**8}.json"
    try:
        r = subprocess.run(
            ["curl", "-sS", "-m", str(timeout), "-o", str(tmp), url],
            capture_output=True,
            text=True,
        )
        if r.returncode != 0:
            raise RuntimeError(r.stderr.strip() or f"curl exit {r.returncode}")
        return json.loads(tmp.read_text())
    finally:
        tmp.exists() and tmp.unlink()


def fetch_month(reporter: str, period: str, tries: int = 4):
    # no customsCode/motCode — older years often empty with those filters
    url = (
        "https://comtradeapi.un.org/public/v1/preview/C/M/HS?"
        f"reporterCode={reporter}&partnerCode={PARTNER_CODES}&partner2Code=0&"
        f"period={period}&flowCode=M&cmdCode={CMD}"
    )
    last = None
    for k in range(tries):
        try:
            j = curl_json(url)
            rows = []
            for rec in j.get("data") or []:
                w = rec.get("netWgt")
                if w is None:
                    continue
                try:
                    wv = float(w)
                except (TypeError, ValueError):
                    continue
                if wv <= 0:
                    continue
                pc = str(rec.get("partnerCode"))
                if pc not in PARTNERS:
                    continue
                rows.append((PARTNERS[pc], wv / 1000.0))  # kg → tonnes
            # collapse duplicates: max
            by = {}
            for name, t in rows:
                by[name] = max(by.get(name, 0.0), t)
            return by
        except Exception as e:
            last = e
            time.sleep(2.0 * (k + 1))
    print(f"    FAIL {period}: {last}")
    return {}


def load_existing(path: Path):
    have = defaultdict(dict)  # date -> partner -> tonnes
    if path.exists():
        with path.open() as f:
            for r in csv.DictReader(f):
                have[r["date"][:7]][r["partner"]] = float(r["tonnes"])
    return have


def save(path: Path, have):
    rows = []
    for d in sorted(have):
        for p, t in sorted(have[d].items()):
            rows.append({"date": d, "partner": p, "tonnes": t})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["date", "partner", "tonnes"])
        w.writeheader()
        w.writerows(rows)


def main():
    y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 2012
    y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2024
    stubs = (
        sys.argv[3].split(",")
        if len(sys.argv) > 3
        else ["twn", "tur", "deu", "nld", "gbr", "phl", "mar"]
    )
    RAW.mkdir(parents=True, exist_ok=True)

    for stub in stubs:
        stub = stub.strip().lower()
        if stub not in IMPORTERS:
            print(f"unknown stub {stub}, skip")
            continue
        code, label = IMPORTERS[stub]
        path = RAW / f"bilateral_{stub}.csv"
        have = load_existing(path)
        print(f"\n{label} ({stub})  existing months={len(have)}")
        n_new = 0
        for y in range(y0, y1 + 1):
            yr_new = 0
            for m in range(1, 13):
                ds = f"{y}-{m:02d}"
                # re-pull if World missing (incomplete month)
                if "World" in have.get(ds, {}):
                    continue
                period = f"{y}{m:02d}"
                got = fetch_month(code, period)
                time.sleep(0.55)
                if not got:
                    print(f"  {ds}  -- empty")
                    continue
                have[ds].update(got)
                n_new += 1
                yr_new += 1
                aus = got.get("Australia", 0.0)
                world = got.get("World", 0.0)
                share = 100 * aus / world if world > 0 else float("nan")
                print(f"  {ds}  World {world/1e6:.2f} Mt  Aus {share:5.1f}%")
            if yr_new:
                save(path, have)
            print(f"  {y}: +{yr_new} months")
        save(path, have)
        # summary Aus share
        w = sum(v.get("World", 0) for v in have.values())
        a = sum(v.get("Australia", 0) for v in have.values())
        print(f"  -> {path.name}: {len(have)} months, overall Aus share "
              f"{(100*a/w if w else float('nan')):.1f}%  (+{n_new} new)")


if __name__ == "__main__":
    main()
