# Diagnostics report (Tasks A–F)

**Repo:** `coalEl/`  
**Inventory:** `diagnostics/out/00_inventory.md` (written before new estimation)  
**Bootstrap seed (Tasks B, F):** `20260906` (`BOOT_SEED` in `diagnostics/common.jl`)  
**Draws:** 1000 · **Block length:** 12 months  

Pipeline estimation code was **not** modified. All numbers below come from scripts under `diagnostics/`. CSV companions live in `diagnostics/out/`.

**Critical pipeline fact (from inventory):** paper text prefers \(z_t^w\) with 2022 M2–M3 zeroed (`02_estimate.jl`), but Figures 1–2 / exporter & leakage tables come from `03_leakage.jl` (**binary** \(z_t\), **2022 included**). Diagnostics label which convention is used.

---

## Task A — Effective sample and event count

**Script:** `diagnostics/task_A.jl`  
**Outputs:** `A_sample_counts.csv`, `A_h0_event_months.csv`  
**LP definition:** first-stage style from `02_estimate.jl` (lags of \(q\) and \(p\), REA + lags, month FE).  
**`n_blocks`:** `cld(tmax-tmin+1, 12)` over estimation-row calendar span.

### A1. Preferred instrument \(z_t^w\), early-2022 dropped — Australian exports

| h | N obs | N \(z>0\) | N \(z=1\) | N \(z=0.5\) | df | n blocks |
|---|------:|----------:|----------:|------------:|---:|---------:|
| 0 | 149 | 8 | 2 | 6 | 117 | 14 |
| 1 | 149 | 7 | 2 | 5 | 117 | 14 |
| 2 | 149 | 8 | 2 | 6 | 117 | 14 |
| 3 | 148 | 8 | 2 | 6 | 116 | 14 |
| 6 | 144 | 8 | 2 | 6 | 112 | 13 |
| 12 | 138 | 8 | 2 | 6 | 106 | 13 |
| 18 | 132 | 8 | 2 | 6 | 100 | 12 |

Price LP under the same instrument/drop has **identical** \(N\), event counts, df, and n blocks at every \(h\) (same row filter on \(q,p,\mathrm{REA}\)). Full \(h=0..18\) for all six instrument×drop cells: `A_sample_counts.csv`.

### A2. Same, but early-2022 **retained** (\(z_t^w\), drop=false), \(h=0\)

| outcome | N obs | N \(z>0\) | N \(z=1\) | N \(z=0.5\) | df | n blocks |
|---------|------:|----------:|----------:|------------:|---:|---------:|
| q_aus / price | 149 | 10 | 2 | 8 | 117 | 14 |

### A3. Literal \(h=0\) event months entering the regression

**\(z_t^w\), drop 2022 (preferred):**  
`2013-01` (0.5), `2017-03` (1.0), `2017-04` (1.0), `2019-02` (0.5), `2020-05` (0.5), `2021-01` (0.5), `2024-01` (0.5), `2024-06` (0.5)  

Same list for q_aus and price (`A_h0_event_months.csv`).

**\(z_t^w\), 2022 retained:** above plus `2022-02` (0.5), `2022-03` (0.5).

Pre-2012 events in `events.csv` never enter (sample starts 2012-01).

---

## Task B — Bootstrap harmonisation

**Script:** `diagnostics/task_B.jl`  
**Output:** `B_bootstrap_harmonisation.csv`  
**Seed:** 20260906 · **NBOOT:** 1000 · **BLOCK:** 12  

Columns: point, OLS SE, Newey–West SE (bandwidth \(h+1\)), MBB 90%, wild 90%, dependent-wild 90%, `implied_SE_from_block = (hi-lo)/(2×1.645)`.

### B1. First-stage objects (\(z_t^w\), drop 2022, REA+month FE)

| object | h | point | SE OLS | SE NW | MBB 90% | wild 90% | DWB 90% | implied SE (MBB) |
|--------|--:|------:|-------:|------:|---------|----------|---------|------------------:|
| aus_exports_log_FS | 0 | −0.1196 | (in CSV) | (in CSV) | [−0.1506, −0.03179] | [−0.1970, −0.04279] | (in CSV) | (in CSV) |
| price_real_log_FS | 0 | −0.002644 | (in CSV) | (in CSV) | **[−0.2721, +0.08541]** | [−0.05900, +0.04980] | (in CSV) | (in CSV) |

Full horizons 0,1,3,6,12 and DWB/SE columns: see CSV (produced at `task_B.jl` bootstrap loops).

### B2. Exporter Mt (03_leakage-style: binary \(z\), no REA), \(h=0\)

| object | point (Mt) | MBB 90% | wild 90% | implied SE (MBB) |
|--------|----------:|---------|----------|------------------:|
| Australia | −2.084 | [−3.159, −0.4974] | [−3.301, −0.8920] | 0.8091 |
| Indonesia | +2.021 | [−2.242, +6.607] | (CSV) | 2.690 |
| USA | +0.2588 | [−1.236, +1.254] | (CSV) | 0.7569 |
| Colombia | −0.3115 | [−1.331, +0.6928] | (CSV) | 0.6150 |
| Canada | −0.03729 | [−0.3337, +0.2911] | (CSV) | 0.1899 |

### B3. Japan (05-style binary \(z\)), \(h=0\)

| object | point | MBB 90% | wild 90% |
|--------|------:|---------|----------|
| japan_total_log | −0.01605 | [−0.08012, +0.02472] | [−0.06364, +0.03275] |
| japan_from_aus_log | −0.03905 | [−0.09055, +0.000077] | [−0.07853, +0.001506] |

### B4. Leakage ratio

| h | point | MBB 90% (ratio recomputed per draw) | OLS/NW/wild/DWB |
|--:|------:|--------------------------------------|-----------------|
| 0 | 0.9265 | [−2.258, +4.974] | **NaN — not available** for a multi-equation ratio without a single residual equation (`task_B.jl` note column) |

---

## Task C — Leave-one-out on the price equation

**Script:** `diagnostics/task_C.jl`  
**Outputs:** `C1_leave_one_block_price.csv`, `C2_leave_one_event_price.csv`, `C3_leave_one_event_qaus.csv`, `C_summary.txt`

Baselines:
- price \(h=0\), \(z_w\) drop 2022: **−0.002644**
- price \(h=0\), full \(z_w\) (2022 in): **+0.06399**
- q_aus \(h=0\), drop 2022: **−0.1196**; full \(z_w\): **−0.1288**

### C1. Leave-one-block-out (price, baseline drop-2022)

Sorted by \(|\Delta|\):

| block start | end | β price | Δ from baseline |
|-------------|-----|--------:|----------------:|
| 2016-01 | 2016-12 | −0.06075 | **−0.05811** |
| 2017-01 | 2017-12 | −0.05459 | −0.05195 |
| 2024-01 | 2024-12 | +0.02237 | +0.02502 |
| 2021-01 | 2021-12 | +0.01751 | +0.02016 |
| 2023-01 | 2023-12 | +0.01674 | +0.01938 |

**Plain statement:** No single block flips the sign of the near-zero baseline, but **2016** and **2017** (Debbie year) blocks move the estimate most (about −0.05 to −0.06). The wide left MBB tail is consistent with block resampling that occasionally drops / reweights the Debbie neighbourhood.

### C2. Leave-one-event-out (price, **full** \(z_w\) including 2022)

| event month | name | β price | Δ vs full |
|-------------|------|--------:|----------:|
| **2022-03** | East coast floods (cont.) | +0.01400 | **−0.04999** |
| 2024-06 | Grosvenor fire | +0.08426 | +0.02027 |
| 2024-01 | Kirrily | +0.08051 | +0.01651 |
| 2022-02 | East coast floods | +0.04889 | −0.01510 |
| 2017-03 | Debbie | +0.07908 | +0.01509 |

Pre-2012 events: Δ = 0 (not in sample).

**Plain statement:** **2022-03 dominates** the positive price response when 2022 is left in \(z_t\). Dropping it cuts β from +0.064 to +0.014. That is why the preferred spec zeros early-2022.

### C3. Leave-one-event-out (q_aus, full \(z_w\))

| event month | name | β q | Δ |
|-------------|------|----:|--:|
| **2017-04** | Debbie (cont.) | −0.07005 | **+0.05873** |
| 2017-03 | Debbie | −0.1685 | −0.03968 |
| 2021-01 | Kimi | −0.1370 | −0.00826 |

**Plain statement:** **Cyclone Debbie (especially 2017-04)** dominates the Australian export first stage.

---

## Task D — Early-2022 flood months

**Script:** `diagnostics/task_D.jl`  
**Outputs:** `D_2022_specs.csv`, `D_gas_availability.txt`

| outcome | h | (1) baseline drop 2022 | (2) include 2022 | (3) include 2022 + post-2022-M2 step | (4) + gas price |
|---------|--:|-----------------------:|-----------------:|--------------------------------------:|----------------:|
| q_aus | 0 | −0.1196 | −0.1288 | −0.1286 | **NaN** |
| q_aus | 3 | −0.001402 | +0.002623 | +0.003810 | **NaN** |
| price | 0 | −0.002644 | **+0.06399** | **+0.06417** | **NaN** |
| price | 3 | −0.05932 | **+0.1342** | **+0.1357** | **NaN** |

**Spec 4 not estimated:** no monthly TTF / Henry Hub / natural-gas series in `data/raw/` or the panel (`D_gas_availability.txt`). To run it, add e.g. `data/raw/ttf.csv` (`date`, price) covering the sample and rebuild.

Including 2022 flips the price impact from ≈0 to **+6.4 log points**; a post-invasion step dummy does **not** undo that.

---

## Task E — Australian tonnage sanity check

**Script:** `diagnostics/task_E.jl`

### E1. Comtrade Aus exports (HS 2701), Mt/month

| year | mean Mt/month | n months | annual Mt |
|------|--------------:|---------:|----------:|
| 2012 | 26.30 | 12 | 315.6 |
| 2013 | 29.64 | 12 | 355.6 |
| 2014 | 32.02 | 12 | 384.2 |
| 2015 | 32.18 | 12 | 386.1 |
| 2016 | 32.37 | 12 | 388.5 |
| 2017 | 31.02 | 12 | 372.2 |
| 2018 | 32.15 | 12 | 385.8 |
| 2019 | 32.75 | 10 | 327.5 |
| 2020 | 30.94 | 12 | 371.3 |
| 2021 | 30.48 | 12 | 365.8 |
| 2022 | 28.27 | 12 | 339.2 |
| 2023 | 29.45 | 12 | 353.3 |
| 2024 | 30.17 | 12 | 362.0 |
| 2025 | 29.66 | 12 | 355.9 |

Sample mean 2012+: **30.50 Mt/month** (`E1_means.txt`).  
Pipeline Table-3 conversion qbar (mt/β from `out/exporter_responses.csv`): **30.69 Mt/month**.  
`06_figures.jl` hard-code: 30.69.

**No ~17 Mt/month shortfall** in current `aus_exports.csv`.

### E2. Mirror (partner-reported imports from Australia)

Built from `bilateral_*.csv` with `partner==Australia` only (`E2_mirror_by_year.csv`).  
**Not a full world mirror** — only pulled importers (JP, KR, IN, DE, NL, TR, UK, PH, …). Will understate Aus exports. Status: `E2_mirror_status.txt`.

### E3. ABS / IEA published totals

**NOT AVAILABLE** in the repo. No ABS AHECC or IEA annual Australian export tonnage files under `data/raw/` (`E3_published_comparison.txt`). Comparison table cannot be built without adding those data.

### E4. Comtrade field

From `00_pull_comtrade.jl`: **`netWgt`** (kg) → tonnes via **`/ 1000.0`**. Filters: HS 2701, World partner, `customsCode=C00`, `motCode=0`. Collapse by date with **maximum** (not sum). Details: `E4_comtrade_field.txt`.

---

## Task F — Replacement without a ratio

**Script:** `diagnostics/task_F.jl`  
**Seed:** 20260906 · Coverage set = same exporters as `03_leakage.jl` gates.

### F1. World total exports (sum of included exporters), Mt levels, binary \(z\)

| h | β (Mt) | SE OLS | MBB 90% |
|--:|-------:|-------:|---------|
| 0 | −0.8930 | 1.844 | [−4.796, +3.086] |
| 1 | +4.451 | 2.028 | [−0.2048, +6.843] |

Full \(h=0..18\): `F1_world_total_mt.csv`.  
Under full replacement, expect ≈0; point at \(h=0\) is negative but band covers zero and large positives.

### F2. Log-then-rescale vs levels in Mt, \(h=0\)

| exporter | Mt from log×qbar | Mt level LP |
|----------|-----------------:|------------:|
| Australia | −2.084 | −1.910 |
| Indonesia | +2.021 | +1.044 |
| USA | +0.2588 | +0.1557 |
| Colombia | −0.3115 | −0.6399 |
| Canada | −0.03729 | −0.02467 |

Full paths: `F2_log_vs_level_mt.csv`. Level and log-rescaled paths differ materially for Indonesia.

### F3–F4. Cumulative leakage (ratio recomputed inside each MBB draw)

| H | CumAUS (Mt) | CumCOMP (Mt) | cum leakage | MBB 90% of ratio |
|--:|------------:|-------------:|------------:|------------------|
| 0 | −2.084 | +1.931 | 0.9265 | [−2.359, +5.091] |
| 3 | −3.432 | +9.879 | 2.878 | [−5.065, +15.45] |
| 6 | −6.191 | +14.45 | 2.334 | [−8.142, +14.60] |
| 12 | −11.25 | +15.71 | 1.397 | [−11.43, +17.83] |
| 18 | −12.22 | +15.59 | 1.275 | [−13.37, +19.87] |

Source: `F3F4_cumulative_leakage.csv`.

### F5. How current Table 5 bands were built

**Ratio recomputed per bootstrap draw** in `03_leakage.jl` (loop assigns `BS[b,:] = lk` from `leakage_path`, then 5%/95% quantiles).  
**Not** formed by dividing marginal confidence bands.  
Evidence file: `F5_current_leakage_bands.txt`.

---

## Stop point

Tasks **A–F** complete. Tasks G–L not run yet (per instructions).

Reproduce:
```bash
cd coalEl
julia --project=. diagnostics/task_A.jl
julia --project=. diagnostics/task_B.jl   # longest
julia --project=. diagnostics/task_C.jl
julia --project=. diagnostics/task_D.jl
julia --project=. diagnostics/task_E.jl
julia --project=. diagnostics/task_F.jl
# or: julia --project=. diagnostics/run_AF.jl
```
