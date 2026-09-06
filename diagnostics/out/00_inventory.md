# Pipeline inventory

Generated before any diagnostic estimation. Sources: scripts under `coalEl/`,
paper bodies in `Proposal&Paper/section2_empirical_strategy.tex` and
`section3_results_discussion.tex`, and current files in `out/`.

## Estimation sample (as coded)

| Object | Sample start | Instrument | Drop 2022 M2–M3? | Controls | Script |
|--------|--------------|------------|------------------|----------|--------|
| First stage / LP-IV elasticities | `SAMPLE_START = "2012-01"` | `SHOCKVAR = :z_w` | Yes (`DROP_EVENTS`) | Outcome & price lags, REA, month FE | `02_estimate.jl` |
| Exporter IRFs + leakage | `"2012-01"` on common exporter span | **Binary** `z` (all events = 1) | **No** | Own lags + month FE only (no REA) | `03_leakage.jl` |
| Importer IRFs | `"2012-01"` | **Binary** `z` | **No** | Own lags + month FE | `05_importer_substitution.jl` |
| HIGH/LOW group LP-IV | `"2012-01"` | `z_w` + drop 2022 | Yes | As in `02_estimate.jl` | `08_estimate_by_group.jl` |

**Important inconsistency:** the paper text prefers \(z_t^w\) with early-2022
excluded, but Figures 1–2 / exporter and leakage tables are produced by
`03_leakage.jl`, which uses **unweighted binary** \(z_t\) and **does not** zero
out 2022-02 / 2022-03.

Panel construction: `01_build_panel.jl` → `data/panel.csv` (shock series
`z`, `z_major`, `z_w` from `data/events.csv`).

## Paper figures → script and input

| Paper figure | File | Produced by | Input CSV |
|--------------|------|-------------|-----------|
| Fig 1 exporters | `out/fig1_exporters.png` | `06_figures.jl` | `out/exporter_responses.csv` |
| Fig 2 leakage | `out/fig2_leakage.png` | `06_figures.jl` | `out/leakage.csv` |
| Fig 3 importers | `out/fig3_importers.png` | `06_figures.jl` | `out/importer_responses.csv` |
| Fig 4 sources | `out/fig4_sources.png` | `06_figures.jl` | `out/importer_responses.csv` |

Fig 2 y-limits are hard-coded `ylims!(ax, -2, 3)` and band values are
`clamp`ed to \([-2, 3]\) in `06_figures.jl` (clipping).

Fig 3/4 plot **Japan, Korea, and India** (all importers present in
`importer_responses.csv`), not Japan alone.

Fig 1: shaded bands drawn **only** for Australia and Indonesia
(`nm in ("Australia", "Indonesia")` in `06_figures.jl`).

## Paper tables → script and input

| Paper table (section 3 labels) | Source |
|--------------------------------|--------|
| `tab:firststage` (Aus log exports, Newcastle price) | `out/first_stage.txt` ← `02_estimate.jl` |
| `tab:exporters` (Mt impact) | `out/exporter_responses.csv` ← `03_leakage.jl` |
| `tab:leakage` | `out/leakage.csv` ← `03_leakage.jl` |
| `tab:japan` | `out/importer_responses.csv` ← `05_importer_substitution.jl` |
| `tab:shares` | `out/import_shares_summary.csv` ← `07_import_shares.jl` |
| `tab:events` | `data/events.csv` (hand-coded) |

Also written but not primary paper tables: `out/elasticities.csv`,
`out/irf_svar.csv`, `out/loo_events.csv`, `out/elasticities_by_group.csv`.

## Bootstrap procedures currently used

| Script | Method | Draws | Block | Seed | What is bootstrapped |
|--------|--------|------:|------:|------|----------------------|
| `03_leakage.jl` | Moving-block on calendar indices; re-runs full `leakage_path` | 1000 | 12 | `MersenneTwister(20260728)` | **Leakage ratio recomputed inside each draw** (`BS[b,:] = lk`). Exporter Mt bands = quantiles of boot betas × fixed `qbar`. |
| `05_importer_substitution.jl` | Moving-block | 1000 | 12 | `20260728` | Log-response coefficients; bands on β |
| `02_estimate.jl` | Wild bootstrap of LP-IV **ratio** (`boot_lpiv`) | 2000 | n/a | `Random.seed!(20260905)` | Elasticity \(\beta^q/\beta^p\) |
| `08_estimate_by_group.jl` | Same family as `02_estimate.jl` | (see script) | | | Group elasticities |

**Answer to Task F5 (preview):** Table 5 / `leakage.csv` bands are percentiles of
the **ratio recomputed per bootstrap draw**, not intervals formed by dividing
marginal confidence bands (`03_leakage.jl` lines 221–230, 297–300).

## Quantity field (Comtrade)

`00_pull_comtrade.jl`: reads `netWgt` (kg), divides by 1000 → tonnes
(line ~173: `wv / 1000.0`). Filters `customsCode=C00`, `motCode=0`,
`partnerCode=0` (World). HS `cmdCode=2701`.

## Mean Aus tonnage used for Mt conversion (Table 3 / Fig 1)

`03_leakage.jl` uses exporter-specific `qbar[j] = mean(finite tonnes)`.
`06_figures.jl` hard-codes `QBAR["Australia"] => 30.69` for display only;
CSV Mt columns use the leakage script’s `qbar`.

## Proxy SVAR

`02_estimate.jl` writes `out/irf_svar.csv` (reduced-form / proxy cross-check
exists in code).
