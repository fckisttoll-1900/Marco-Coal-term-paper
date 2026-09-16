# Australian coal supply disruptions — replication

Replication package for the Macro II term paper on weather-driven Australian
coal supply shocks, replacement by other exporters, and Japanese import
adjustment.

## What lives here

| Path | Role |
|------|------|
| `01_build_panel.jl` | Builds `data/panel.csv` (incl. severity-weighted `z_w`) |
| `data/events.csv` | Narrative disruption list (severity labels) |
| `data/raw/` | Monthly inputs (exports, price, REA, CPI, bilateral) |
| `diagnostics/` | Harmonised local projections used in the paper |
| `diagnostics/out/` | Canonical paper tables / randomisation outputs |

Scratch pull scripts, draft reports, and one-off checks are **not** included.

## Setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Reproduce the paper baseline

```bash
# 1. Rebuild the panel (optional if data/panel.csv is already current)
julia --project=. 01_build_panel.jl

# 2. Baseline paths, leakage, and MBB bands (Tables 2–4 / Fig. 1 material)
julia --project=. diagnostics/task_UVX.jl

# 3. Headline keep-2022 randomisation p-values (valid-row Nov–Apr pool)
julia --project=. diagnostics/task_headline_rand.jl

# 4. Paper figures (optional)
julia --project=. diagnostics/make_paper_figures.jl
```

Baseline design: severity-weighted `z_w`, early-2022 retained, six outcome lags,
contemporaneous REA + six lags, month-of-year FE, no price lags.

## Canonical outputs (do not delete)

- `diagnostics/out/X_baseline_keep2022.csv` — horizon paths
- `diagnostics/out/U1_per_horizon_leakage_keep2022.csv` — impact leakage
- `diagnostics/out/U2_cum_leakage_keep2022.csv` — cumulative leakage
- `diagnostics/out/headline_rand_keep2022.csv` — headline randomisation *p*-values
- `diagnostics/out/V_randomisation.csv` — additional keep-2022 randomisation cells
- `diagnostics/out/C3_leave_one_event_qaus.csv` — leave-one-event (price-lag FS)
- `diagnostics/out/G_unpredictability.*` — instrument unpredictability checks
- `diagnostics/out/fig_data_F1_keep2022.csv`, `fig_data_F3_price_keep2022.csv`

## Notes

- The Mt Australian-export response is a **direct levels** local projection, not
  a rescaled log coefficient.
- Headline randomisation draws placebos from the November–April pool restricted
  to valid regression rows, preserving the effective treated-month weight multiset.
