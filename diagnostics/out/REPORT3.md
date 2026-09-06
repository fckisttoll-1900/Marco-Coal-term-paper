# Diagnostics report 3 (Tasks P, Q, R, S, T)

**Ground rules:** pipeline scripts and `out/` untouched. All work under `diagnostics/`.  
**Harmonised LP:** \(z_t^w\), DROP 2022-02/03, 6 own lags, REA + 6 lags, month FE (`diagnostics/harm_lp.jl`).  
**MBB:** seed `20260906`, **10 000** draws, 12-month blocks (Tasks P–Q, R). Task T bands use 1 000 draws.

---

## Task P — Fix world-aggregate leakage cell [BLOCKING]

**Scripts:** `diagnostics/exporters_harm.jl`, `diagnostics/task_PQ.jl`  
**Outputs:** `P_definition.txt`, `P_world_leakage_fixed.csv`

### Exact definitions (as coded)

| object | definition |
|--------|------------|
| **Coverage set** | Australia, Indonesia, USA, Colombia, Canada (`03_leakage` density gates; Russia/SA/Mongolia fail) |
| **\(W^{\mathrm{inc}}_t\)** | sum of observed Mt over that set in month \(t\), **requiring Australia observed** — **includes Australia** |
| **\(W^{\mathrm{exc}}_t\)** | same sum **excluding Australia** |
| **Australia** | LP of Aus export Mt (levels) |
| **Leakage\(_h\)** | \(-\,\beta_{W^{\mathrm{exc}},h}\,/\,\beta_{\mathrm{Aus},h}\) when \(\beta_{\mathrm{Aus},h}< -0.01\); recomputed **inside each of 10 000 MBB draws** |

### Why M2’s world cell was inconsistent

M2 built \(Q^{\mathrm{world}}\) as **\(W^{\mathrm{inc}}\)** (includes Australia, response ≈ −3.67 Mt) then coded  
`implied_leak = 1 − ΔW/(−ΔAUS)`.  
When \(W\) includes Aus the correct identity is  
\(\mathrm{leak}=(\Delta W-\Delta\mathrm{AUS})/(-\Delta\mathrm{AUS})=1+\Delta W/(-\Delta\mathrm{AUS})\).  
M2’s minus sign on the \(\Delta W\) term produces a spurious ≈ +2 whenever \(\Delta W\approx\Delta\mathrm{AUS}<0\).

### Corrected paths, \(h=0..18\) (Mt levels; MBB 90%)

| h | world inc Aus | world ex Aus | Australia | leakage \(-\beta_{W^{\mathrm{exc}}}/\beta_{\mathrm{Aus}}\) |
|--:|--------------:|-------------:|----------:|----------------------------------------------------------:|
| 0 | −3.674 [−9.34, +3.32] | −0.526 [−6.14, +5.78] | −3.283 [−4.37, −1.00] | **−0.160** [−2.64, +2.62] |
| 1 | +1.454 | +3.674 | −2.171 | +1.693 |
| 2 | −1.112 | +0.682 | −1.609 | +0.424 |
| 3 | +1.553 | +0.710 | +0.326 | NaN (Aus ≥ 0) |
| 6 | +6.990 | +7.879 | −1.288 | +6.119 |
| 12 | −1.735 | −1.757 | −0.170 | −10.366 |
| 18 | +5.051 | +4.663 | +0.204 | NaN |

Full table: `P_world_leakage_fixed.csv`.

**Consistency check at \(h=0\):** Australia alone (−3.28) is the same order as world-including-Aus (−3.67); world-ex-Aus is near zero (−0.53). Impact leakage ≈ **0**, not +2. The M2 “implied leakage +2.036” was a **formula bug**, not an economic result.

---

## Task Q — Horizon paths [BLOCKING]

**Script:** `diagnostics/task_PQ.jl` (shared MBB with P)  
**Outputs:** `Q1_*.csv/txt`, `Q2_competitors_mt.csv`, `Q3_*.csv/txt`, `Q4_cumulative_leakage.csv`

### Q1. Australia Mt — catch-up at long \(h\)?

| h | β Mt [90%] |
|--:|-----------:|
| 0 | **−3.283** [−4.37, −1.00] |
| 6 | −1.288 [−3.70, +2.62] |
| 12 | −0.170 [−3.74, +3.55] |
| 15 | +1.061 [−3.65, +3.61] |
| 18 | +0.204 [−3.81, +3.59] |

Point estimates turn **mildly positive** at \(h\in\{9,10,15,16,18\}\), but **no long-horizon cell has a positive lower 90% band**.  
**Verdict: no credible catch-up** — recovery to ~0 within noise, not a statistically clear rebound.

### Q2. Competitors Mt and their sum

At \(h=0\) (levels, harmonised):

| exporter | β Mt [90%] |
|----------|-----------:|
| Indonesia | −0.594 [−4.67, +2.16] |
| USA | −0.111 [−2.96, +2.43] |
| Colombia | −0.454 [−3.59, +1.64] |
| Canada | −0.233 [−0.82, +0.14] |
| **SUM ex-Aus** | **−1.392** [−8.68, +3.24] |

(Note: sum of country LPs ≠ LP of \(W^{\mathrm{exc}}\) (−0.53) because missingness/samples differ.)

At \(h=6\), SUM ex-Aus flips to **+2.69** (Indonesia +1.67, Colombia +1.92) with wide bands covering 0. At \(h=12\), SUM ≈ −1.08. Full paths: `Q2_competitors_mt.csv`.

### Q3. Japan total / from-Aus / non-Aus (logs + Mt), \(h=0..12\)

Impact (\(h=0\)):

| series | log [90%] | Mt [90%] |
|--------|----------:|---------:|
| Total | −0.071 [−0.20, +0.02] | −0.990 [−2.80, +0.22] |
| From Australia | −0.059 [−0.21, −0.00] | −0.529 [−1.87, −0.04] |
| Non-Australia | −0.135 [−0.23, +0.06] | −0.757 [−1.32, +0.35] |

**Does non-Australian recover after \(h=0\)?**  
Yes in the **point path**: non-Aus log moves from −0.135 at \(h=0\) to ≈0 by \(h=1\) and briefly positive at \(h=3\) (+0.055); Mt from −0.76 to −0.09 then +0.26 at \(h=3\).  
Bands always cover 0 after impact. So: **short-run comovement down, then noisy mean-reversion — not a clean substitution boom**.

### Q4. Cumulative leakage every \(H=0..18\)

Coded as \(\mathrm{cum\text{-}leak}_H = -(\sum_{h=0}^{H}\beta^{\mathrm{comp}}_h)\,/\,(\sum_{h=0}^{H}\beta^{\mathrm{Aus}}_h)\) where \(\beta^{\mathrm{comp}}\) is the **sum of individual competitor LPs** (MBB ratio inside each draw).

| H | Cum Aus | Cum EX | cum leak [90%] |
|--:|--------:|-------:|----------------|
| 0 | −3.283 | −1.392 | **−0.424** [−4.02, +1.54] |
| 3 | −6.737 | −0.686 | −0.102 |
| 6 | −7.550 | +4.230 | **+0.560** [−15.2, +17.5] |
| 12 | −10.977 | +6.102 | **+0.556** [−19.6, +21.3] |
| 18 | −10.659 | +10.437 | **+0.979** [−18.5, +22.3] |

Point cumulative leakage rises toward ~1 by \(H=18\), but **intervals are huge** (cover 0 and values ≫1). Do not treat the path as precise.

---

## Task R — Randomisation inference, extended [BLOCKING]

**Script:** `diagnostics/task_R.jl` · **Output:** `R_randomisation.csv`  
10 000 placebos; event-month pool as in Task H.

| test | h | β | rand \(p\) |
|------|--:|--:|----------:|
| R1 price, **keep 2022** (\(z_w\)) | 0 | +0.0649 | **0.235** |
| R2 price, drop 2022 | 0 | −0.0049 | 0.930 |
| R2 price, drop 2022 | 3 | −0.0709 | 0.627 |
| R2 price, drop 2022 | 6 | −0.0164 | 0.943 |
| R3 Aus exports (log), drop 2022 | 0 | −0.1163 | **0.003** |
| R3 Aus exports | 3 | +0.0116 | 0.808 |
| R3 Aus exports | 6 | −0.0403 | 0.313 |
| R3 Aus exports | 12 | +0.0003 | 0.994 |
| R4 world ex-Aus Mt | 0 | −0.526 | 0.901 |
| R4 world ex-Aus Mt | 6 | +7.879 | 0.159 |
| R4 world ex-Aus Mt | 12 | −1.757 | 0.821 |

**Takeaway:** only the **Australian export impact** is randomised-significant. The keep-2022 price bump (+0.065) is **not** (\(p\approx0.23\)). Competitor/world-ex responses are indistinguishable from placebos.

---

## Task S — Unpredictability with fewer lags

**Script:** `diagnostics/task_S.jl` · **Outputs:** `S_unpredictability_lags.csv`, `S_interpretation.txt`  
Regress \(z\) / \(z_w\) on own lags + REA lags + month FE (same idea as Task G).

| instrument | lags | # regressors | F | p | \(R^2\) | N |
|------------|-----:|-------------:|--:|--:|--------:|--:|
| z | 6 | 18 | 1.557 | **0.083** | 0.276 | 148 |
| \(z_w\) | 6 | 18 | 1.515 | 0.096 | 0.258 | 148 |
| z | 3 | 9 | 1.965 | **0.048** | 0.202 | 154 |
| \(z_w\) | 3 | 9 | 1.607 | 0.119 | 0.173 | 154 |
| z | 1 | 3 | 2.177 | 0.093 | 0.129 | 161 |
| \(z_w\) | 1 | 3 | 1.414 | 0.241 | 0.103 | 161 |

**Verdict:** the marginal \(p\approx0.083\) at 6 lags for binary \(z\) **does not look like pure overfitting** — with **3 lags** the \(p\)-value **falls to 0.048** (stronger, not weaker). \(z_w\) stays weaker (\(p\approx0.10\)–0.24). Mild predictability of the narrative instrument remains a caveat; preferred \(z_w\) is less predictable than binary \(z\).

---

## Task T — HIGH/LOW groups under harmonised LP

**Script:** `diagnostics/task_T.jl` · **Output:** `T_high_low_harmonised.csv`  
HIGH = Japan; LOW = Germany, Netherlands, Turkey, UK, Philippines. Horizon \(h=0..12\).

### Individual LOW countries: inestimable

Harmonised LP requires continuous own lags + REA lags. Sparse Comtrade bilateral series leave almost no usable rows / event mass:

| country | finite World months | usable LP rows \(h=0\) | \(\lvert z\rvert\) in rows |
|---------|--------------------:|-----------------------:|---------------------------:|
| Japan | 161 | 131 | 3.0 |
| Germany | 108 | 7 | 0 |
| Netherlands | 89 | 0 | 0 |
| Turkey | 104 | 16 | 0.5 |
| UK | 112 | 22 | 0.5 |
| Philippines | 86 | 1 | 0 |

Individual LOW β cells are **NaN** (gate: ≥40 rows and \(\sum\lvert z\rvert\ge1.5\)). Report **HIGH (Japan)** and **LOW pooled Mt sum** only.

### HIGH — Japan (harmonised)

| h | total log | from-Aus log | non-Aus log | total Mt | from-Aus Mt | non-Aus Mt |
|--:|----------:|-------------:|------------:|---------:|------------:|-----------:|
| 0 | −0.071 | −0.059 | −0.135 | −0.990 | −0.529 | −0.757 |
| 6 | +0.050 | −0.032 | −0.004 | +0.502 | −0.321 | −0.095 |
| 12 | −0.021 | −0.005 | −0.077 | −0.227 | +0.028 | −0.401 |

(MBB 90% for Japan in CSV; from-Aus impact band excludes 0.)

### LOW — pooled sum of five reporters (Mt)

| h | total | from-Aus | non-Aus |
|--:|------:|---------:|--------:|
| 0 | **+0.302** [−4.64, +2.34] | +0.121 [−0.10, +0.50] | +0.659 [−2.97, +2.84] |
| 6 | −1.193 [−4.08, +2.32] | −0.133 [−0.36, +0.35] | −1.133 [−4.13, +2.06] |
| 12 | +1.717 [−3.44, +3.75] | +0.270 [−0.45, +0.45] | +1.523 [−3.20, +3.82] |

**Contrast:** Japan contracts on impact (total and both legs); the LOW pool’s point estimate is a small **positive** total/non-Aus response with bands covering 0 — no evidence of a shared LOW substitution pattern under the harmonised spec, and country-level LPs are not identified.

---

## Bottom line for the paper

1. **M2 world leakage +2 was a coding bug**; corrected impact leakage ≈ **0** (ex-Aus world flat; Aus ≈ −3.3 Mt).  
2. **No credible Australian catch-up**; competitor offsets noisy; Japan non-Aus **recovers toward zero after \(h=0\)** in points only.  
3. **Randomisation:** Aus exports \(h=0\) survive; price (keep or drop 2022) and world-ex do not.  
4. **Unpredictability:** fewer lags make binary-\(z\) predictability **stronger**, not weaker.  
5. **HIGH/LOW:** only Japan is estimable country-by-country; LOW needs pooling and still shows null bands.
