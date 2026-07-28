# Coal demand elasticities from Australian supply disruptions

Minimal pipeline: download seven CSVs, run two Julia scripts, get elasticities.

## 0. Setup

```bash
julia --project=. -e 'using Pkg; Pkg.add(["CSV","DataFrames","Statistics","Random"])'
```

## 1. Download data into `data/raw/`

Each file needs a `date` column (`YYYY-MM` or any parseable date) and one
value column. Rename to the filenames below. Start with the first three —
that alone gives you a result.

| File | Series | Where |
|---|---|---|
| `coal_price.csv` | Newcastle thermal, USD/t, monthly from 1990 | World Bank Pink Sheet, "Monthly Prices" sheet → <https://www.worldbank.org/en/research/commodity-markets> |
| `imports_kor.csv` | Korea coal imports, tonnes, monthly from 1995 | Korea Customs Service / KOSIS trade statistics |
| `aus_exports.csv` | Australian coal exports, tonnes | ABS *International Trade in Goods* (AHECC 27011213 coking, 27011299 thermal), or DFAT monthly pivot tables → <https://www.dfat.gov.au/trade/trade-and-investment-data-information-and-publications/trade-statistics/merchandise-trade-statistical-pivot-tables> |
| `imports_jpn.csv` | Japan coal imports, tonnes | Ministry of Finance trade statistics via e-Stat → <https://www.e-stat.go.jp/en> |
| `imports_ind.csv` | India coal imports, tonnes | Ministry of Commerce tradestat, or UN Comtrade HS 2701 |
| `rea.csv` | Kilian global real activity index | Dallas Fed → <https://www.dallasfed.org/research/igrea> |
| `uscpi.csv` | US CPI (CPIAUCSL), for deflating | FRED |

Missing importer files are skipped automatically, so Korea-only works.

**Use tonnage, not value.** Export *value* moves with price, so a supply
disruption can raise it — feeding value into the VAR flips the sign of your
shock. If only value is available, divide by tonnage to get a realised unit
price (useful), never the reverse.

## 2. Run

```bash
julia --project=. 01_build_panel.jl    # → data/panel.csv
julia --project=. 02_estimate.jl       # → out/elasticities.csv, out/irf_svar.csv
```

## 3. What it does

`01_build_panel.jl` merges the raw files onto a monthly grid, deflates the
price by US CPI, and attaches the narrative shock series built from
`data/events.csv` (binary `z`, plus `z_major` for major events only).

`02_estimate.jl` runs:

1. **First stage** — does the shock actually cut Australian exports and raise
   the price? If not, stop; nothing downstream is valid.
2. **Naive OLS** — the uninstrumented elasticity, reported deliberately as the
   simultaneity-biased benchmark for the motivation section.
3. **LP-IV** — elasticity at horizons 0, 3, 6, 12 months by importer, with
   wild-bootstrap 90% bands.
4. **Proxy SVAR** — impact elasticity as a cross-check on the LP-IV.

## 4. Interpreting the output

- **Korea, Japan** (negligible domestic production): import demand ≈ total
  demand, so the estimate maps directly onto a COALMOD-World consumption
  elasticity.
- **India** (large domestic production): this is an **import** demand
  elasticity, more elastic than consumption. Do **not** feed it to COALMOD as
  a consumption parameter — the model already solves the domestic/import
  margin endogenously, so you would double-count the substitution.

## 5. Robustness switches

In `02_estimate.jl`:

- `SHOCKVAR = :z_major` — major events only.
- `NLAGS` — 6 is the default; try 3 and 12.
- `H` — IRF horizon.

## 6. Precision, before you start

Simulated with a true elasticity of −0.25, ~300 months, and a 12% supply loss
per event:

| events | sd | 90% interval |
|---|---|---|
| 6 | 0.13 | [−0.52, −0.11] |
| 11 | 0.10 | [−0.45, −0.13] |
| 20 | 0.06 | [−0.37, −0.17] |

Holding events at 11 and varying disruption size:

| supply loss | sd |
|---|---|
| 5% | 1.48 (useless) |
| 12% | 0.10 |
| 25% | 0.04 |

**Event size matters far more than event count.** Effort is better spent
documenting the tonnage lost in the few large events (2010–11 floods, Cyclone
Debbie) than adding marginal ones. Expect to distinguish "inelastic" from
"zero", not −0.1 from −0.3.
