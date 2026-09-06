#!/usr/bin/env python3
"""
04b_pull_aus_shares_annual.py
-----------------------------
Annual Comtrade Aus-share pull (World + Australia only) for the advisor
HIGH vs LOW dependence contrast. Saves progress after every successful year.

Output:
  data/raw/aus_share_annual.csv
  out/import_shares_annual.csv
  out/import_shares_summary.csv
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "out"
CMD = "2701"

IMPORTERS = [
    ("jpn", "392", "Japan", True),
    ("kor", "410", "Korea", True),
    ("twn", "490", "Taiwan", True),
    ("chn", "156", "China", False),
    ("ind", "699", "India", False),
    ("tur", "792", "Turkey", True),
    ("vnm", "704", "Vietnam", True),
    ("deu", "276", "Germany", True),
    ("nld", "528", "Netherlands", True),
    ("bra", "76", "Brazil", True),
    ("tha", "764", "Thailand", True),
    ("phl", "608", "Philippines", True),
    ("esp", "724", "Spain", True),
    ("gbr", "826", "United Kingdom", True),
    ("mar", "504", "Morocco", True),
]

HIGH_CUT = 0.50
LOW_CUT = 0.25


def classify(s: float) -> str:
    if s >= HIGH_CUT:
        return "HIGH"
    if s < LOW_CUT:
        return "LOW"
    return "MID"


def fetch_year(reporter: str, year: int, tries: int = 4):
    url = (
        "https://comtradeapi.un.org/public/v1/preview/C/A/HS?"
        f"reporterCode={reporter}&partnerCode=0,36&partner2Code=0&"
        f"period={year}&flowCode=M&cmdCode={CMD}&customsCode=C00&motCode=0"
    )
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "coalEl/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                j = json.loads(r.read().decode())
            data = j.get("data") or []
            world = aus = 0.0
            for rec in data:
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
                if wv <= 0:
                    continue
                pc = str(rec.get("partnerCode"))
                if pc == "0":
                    world = max(world, wv)
                elif pc == "36":
                    aus = max(aus, wv)
            return world, aus
        except Exception as e:
            last = e
            time.sleep(2.0 * (k + 1))
    raise RuntimeError(f"failed {reporter} {year}: {last}")


def main():
    y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 2015
    y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2023
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    out_path = RAW / "aus_share_annual.csv"
    # resume if present
    have = set()
    rows = []
    if out_path.exists():
        with out_path.open() as f:
            for r in csv.DictReader(f):
                rows.append(r)
                have.add((r["importer"], int(r["year"])))
        print(f"resuming with {len(rows)} existing rows")

    print(f"Annual Comtrade Aus shares, {y0}–{y1}\n")
    for stub, code, label, clean in IMPORTERS:
        print(f"{label}{'' if clean else '  (domestic production)'}")
        for y in range(y0, y1 + 1):
            if (label, y) in have:
                print(f"  {y}  (cached)")
                continue
            try:
                world, aus = fetch_year(code, y)
            except Exception as e:
                print(f"  {y}  FAIL {e}")
                time.sleep(2.0)
                continue
            time.sleep(0.8)
            if world <= 0:
                print(f"  {y}  -- no World total")
                continue
            share = aus / world
            print(f"  {y}  Aus share {100*share:5.1f}%  (World {world/1e9:.1f} Mt)")
            rows.append(
                {
                    "importer": label,
                    "stub": stub,
                    "year": str(y),
                    "world_mt": f"{world/1e9:.6f}",
                    "aus_mt": f"{aus/1e9:.6f}",
                    "aus_share": f"{share:.6f}",
                    "clean_margin": str(clean).lower(),
                }
            )
            have.add((label, y))
            # incremental save
            with out_path.open("w", newline="") as f:
                w = csv.DictWriter(
                    f,
                    fieldnames=[
                        "importer",
                        "stub",
                        "year",
                        "world_mt",
                        "aus_mt",
                        "aus_share",
                        "clean_margin",
                    ],
                )
                w.writeheader()
                w.writerows(rows)
        print()

    # summaries
    by = defaultdict(list)
    for r in rows:
        by[r["importer"]].append(r)

    summary = []
    for label, rs in by.items():
        aus = sum(float(r["aus_mt"]) for r in rs)
        world = sum(float(r["world_mt"]) for r in rs)
        share = aus / world if world else 0.0
        years = sorted(int(r["year"]) for r in rs)
        clean = rs[0]["clean_margin"]
        summary.append(
            {
                "importer": label,
                "stub": rs[0]["stub"],
                "clean_margin": clean,
                "years": len(rs),
                "span": f"{years[0]}-{years[-1]}",
                "world_mt_yr": f"{world/len(rs):.3f}",
                "aus_share": f"{share:.6f}",
                "group": classify(share),
            }
        )
    summary.sort(key=lambda r: float(r["aus_share"]), reverse=True)

    with (OUT / "import_shares_annual.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with (OUT / "import_shares_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    print("=" * 66)
    print("GROUPING")
    for g in ("HIGH", "MID", "LOW"):
        names = [r["importer"] for r in summary if r["group"] == g]
        print(f"  {g}: {', '.join(names) if names else '(none)'}")
    print("\nwrote", out_path)
    print("wrote", OUT / "import_shares_summary.csv")


if __name__ == "__main__":
    main()
