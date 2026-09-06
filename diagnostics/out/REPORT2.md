# Diagnostics report 2 (Tasks M, N, O, G, H, I, J, L)

**Ground rules:** pipeline scripts and `out/` untouched. All work under `diagnostics/`.  
**Harmonised LP:** \(z_t^w\), DROP 2022-02/03, 6 own lags, REA + 6 lags, month FE (`diagnostics/harm_lp.jl`).  
**MBB seed:** `20260906` · **Task M/H draws:** 10 000 · **block:** 12 months.

Inventory reminder: current Fig 1–2 / leakage use **binary \(z\)**, **keep 2022**, **no REA** (`03_leakage.jl`).

---

## Task M — Harmonise specification [BLOCKING]

**Script:** `diagnostics/task_M.jl`  
**Outputs:** `M1_exporters_current_vs_harm.csv`, `M2_leakage_current_vs_harm.csv`, `M3_japan_current_vs_harm.csv`, `M4_first_stage_check.txt`

### M4. Does harmonised Aus log \(h=0\) equal Table 2 (−0.1196)?

| estimator | β |
|-----------|--:|
| Harmonised (own lags + REA only) | **−0.1163** |
| `02_estimate` first_stage_F (q lags + **p lags** + REA) | **−0.1196** |
| Table 2 target | −0.1196 |

**Answer: NO for the Task-M control set; YES for the Table-2 control set.**  
Why: Table 2 includes **6 lags of log price** in addition to own quantity lags and REA. Task M’s paper-stated leakage/importer controls omit price lags. FS-style match is confirmed in `M4_first_stage_check.txt`.

### M1. Exporter Mt at \(h=0\) (current vs harmonised)

| exporter | current Mt [MBB] | harm log×qbar [MBB 10k] | harm levels |
|----------|-----------------:|-------------------------|------------:|
| Australia | −2.084 [−3.19, −0.62] | **−3.548 [−4.60, −0.91]** | −3.283 |
| Indonesia | +2.021 [−2.18, +6.63] | **−0.021 [−4.32, +4.32]** | −0.594 |
| USA | +0.259 | −0.005 [−3.19, +2.50] | −0.111 |
| Colombia | −0.311 | +0.181 [−3.05, +2.01] | −0.454 |
| Canada | −0.037 | −0.195 [−0.87, +0.21] | −0.233 |

Full \(h=0..18\): `M1_*.csv`.  

**Substance:** under the paper’s stated controls, the **Indonesian offset disappears** and the Australian shortfall **rises** (~3.5 Mt). The current Fig 1 replacement story is **specification-dependent**.

### M2. Leakage (ratio recomputed inside each of 10 000 draws)

| estimator | measure | H | current | harm | harm MBB 90% |
|-----------|---------|--:|--------:|-----:|--------------|
| log_qbar | impact | 0 | 0.9265 | **−0.0116** | [−3.49, +1.98] |
| levels | impact | 0 | — | −0.424 | [−3.96, +1.41] |
| log_qbar | cumulative | 6 | — | 0.649 | [−14.4, +16.9] |
| log_qbar | cumulative | 12 | — | 0.508 | [−21.5, +22.8] |
| world_agg | world Mt resp. | 0 | — | −3.674 | [−9.34, +3.32] |
| world_agg | implied leak | 0 | — | 2.036 | [−0.41, +5.00] |

Harmonised **impact leakage ≈ 0** (band covers 0 and 1). Current Table 5’s 0.93 does **not** survive harmonisation.

### M3. Japan \(h=0\) (harmonised)

| series | current β (log) | harm log | harm Mt |
|--------|----------------:|---------:|--------:|
| Total | −0.0160 | **−0.0713** | −0.990 |
| From Australia | −0.0390 | **−0.0593** | −0.529 |
| Non-Australia | — | **−0.1354** | −0.757 |

Under harmonisation, Japan’s **total** falls more, and **non-Australian** imports also fall — weaker “pure substitution” reading than the paper’s current binary-\(z\) figures.

---

## Task N — Gas control from Pink Sheet [BLOCKING]

**Script:** `diagnostics/task_N.jl`  
Extracted `Natural gas, Europe` and `Liquefied natural gas, Japan`, CPI-deflated like coal (`N_gas_on_panel.csv`).

### Price responses, \(z_w\) with **2022 RETAINED**

| h | (1) drop 2022 | (2) keep 2022 | (3) + step | (4a) gas EU | (4b) LNG JP | (4c) both |
|--:|--------------:|--------------:|-----------:|------------:|------------:|----------:|
| 0 | −0.0026 | **+0.0640** | +0.0642 | **+0.0452** | +0.0652 | **+0.0478** |
| 3 | −0.0593 | +0.1342 | +0.1357 | +0.0531 | +0.1219 | +0.0746 |
| 6 | +0.0177 | +0.2443 | +0.2475 | +0.1238 | +0.2040 | +0.1412 |
| 12 | +0.0811 | +0.1337 | +0.1360 | −0.0366 | −0.0089 | +0.0345 |

**Verdict (`N_verdict.txt`):** the positive \(h=0\) price response **survives** Europe gas and both-gas controls (+0.045 to +0.048).  
→ **2022 exclusion is not defensible on gas-confounding alone.** If the paper keeps 2022 in \(z_t\), it should report a **positive** Newcastle price response; if it drops 2022, it needs a **different** justification (e.g. Ukraine energy crisis as a distinct shock, not “gas prices in the regression”).

---

## Task O — Data gaps

**Script:** `diagnostics/task_O.jl`

### O1. Material gaps (estimation span)
- **Aus exports:** `2019-03`, `2019-07`, plus 2026-01…07 (panel longer than trade)
- **Price:** `2025-10`, `2026-07`
- **Competitors:** many early-sample gaps (USA/Colombia/Indonesia start later)
- Full list: `O1_missing_months.csv`

### O2. 2019 and event/horizon loss
- Aus missing **2019-03** and **2019-07**
- Event **2019-02** (NQ floods) **enters** \(h=0\) FS regression: **YES**
- Lost for that event at **\(h=1\)** and **\(h=5\)** (needs missing neighbours); other listed horizons still enter (`O2_event_horizon_impact.txt`)

### O3. Comtrade one-row-per-date
Assertion on current `aus_exports` / `exp_*` / `imports_*`: **does not fire** — all unique by date (`O3_comtrade_unique_date.txt`). Pipeline uses `maximum` by date after pull.

---

## Task G — Instrument unpredictability (harmonised sample)

**Script:** `diagnostics/task_GHIJL.jl` → `G_unpredictability.*`

| instrument | joint F (18 lags) | p-value | \(R^2\) | N |
|------------|------------------:|--------:|--------:|--:|
| \(z\) | 1.557 | 0.0830 | 0.2756 | 148 |
| \(z_w\) | 1.515 | 0.0964 | 0.2582 | 148 |

Neither rejects unpredictability at 5%; both are marginally significant at 10%. Not a clean pass, not a clear fail.

---

## Task H — Randomisation inference [PRIORITY]

**10 000** placebo event sets, Nov–Apr only, excluding true events; severity pattern preserved. Seed `20260906`.

| outcome (harmonised \(h=0\)) | β | randomisation \(p\) |
|------------------------------|--:|--------------------:|
| Aus exports | −0.1163 | **0.0031** |
| Price | −0.0049 | **0.9303** |

Aus export response is **highly unlikely** under seasonal placebos. Price response is **indistinguishable** from placebo (as expected under the drop-2022 null).

---

## Task I — Japan decomposition (harmonised)

**Outputs:** `I1_japan_non_aus.csv` (MBB \(N=1000\) for bands), `I2_accounting_gap.txt`

| series | \(h=0\) β | MBB 90% |
|--------|----------:|---------|
| non-Aus log | −0.1354 | (see CSV) |
| non-Aus Mt | −0.757 | (see CSV) |
| total log | −0.0713 | |
| from-Aus log | −0.0593 | |

**Accounting gap** (share \(=0.663\)):  
mechanical \(= 0.663 \times (-0.0593) = -0.0393\);  
observed total \(= -0.0713\);  
**gap \(= -0.0319\)** [MBB 90% −0.125, +0.045].  
Gap band covers 0 — cannot reject mechanical share accounting.

---

## Task J — Promised robustness (harmonised LP)

`J_robustness.csv`:

| spec | \(h\) | β q_aus | β price |
|------|------:|--------:|--------:|
| baseline \(z_w\) drop 2022 real | 0 | −0.1163 | −0.0049 |
| binary \(z\) drop 2022 | 0 | −0.0673 | −0.0127 |
| major only drop 2022 | 0 | −0.1422 | +0.0340 |
| \(z_w\) **keep** 2022 | 0 | −0.1157 | **+0.0649** |
| drop mines | 0 | −0.1173 | +0.0187 |
| nominal price | 0 | −0.1163 | −0.0087 |

**Proxy SVAR:** **EXISTS** — `out/irf_svar.csv` from `02_estimate.jl` (`J_proxy_svar.txt`).  
HIGH/LOW group: existing `08_estimate_by_group.jl` output noted; full bilateral re-run under harm LP not duplicated here (`J_high_low_note.txt`).

---

## Task L — Text / figure consistency

From `L_consistency.txt`:
1. **True price β = −0.0026.** Cite Table 2’s −0.003; “+0.00” is a bad gloss.
2. **Fig 1 bands:** YES for Australia & Indonesia only.
3. **Fig 3/4:** plot **Japan, Korea, India** — “Japan:” captions are wrong/incomplete.
4. **Fig 2:** ylim/clamp **[−2, 3]**; leakage.csv reaches **−6.03 … +4.39** → **clipped**.

---

## Bottom line for the paper

1. **Harmonising to the stated spec kills the Indonesia replacement / 0.93 leakage story** at impact; Aus shortfall looks larger; Japan looks more like a shortfall than pure switch.
2. **Keeping 2022 in \(z_t\) produces a positive price response that gas controls do not remove** — exclusion needs a non-gas justification, or the paper should report the positive price.
3. **Randomisation:** Aus export effect is real (\(p\approx 0.003\)); price null under drop-2022 is confirmed (\(p\approx 0.93\)).
4. Fix captions (Fig 3/4), Fig 2 clipping note, and “+0.00” price wording.

Reproduce:
```bash
cd coalEl
julia --project=. diagnostics/task_M.jl      # ~2 min, 10k MBB
julia --project=. diagnostics/task_N.jl
julia --project=. diagnostics/task_O.jl
julia --project=. diagnostics/task_GHIJL.jl  # includes H 10k placebos
```
