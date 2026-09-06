# Diagnostics report 4 (Tasks U, V, W, X)

**Ground rules:** pipeline scripts and `out/` untouched. All work under `diagnostics/`.  
**New baseline LP:** \(z_t^w\), **2022 RETAINED**, 6 own lags, REA + 6 lags, month FE (`z_baseline` in `harm_lp.jl`).  
**MBB / placebos:** seed `20260906`, **10 000** draws, 12-month blocks.

---

## Task U — Leakage on \(W^{\mathrm{exc}}\) [BLOCKING]

**Script:** `diagnostics/task_UVX.jl`  
**Outputs:** `U1_*.csv`, `U2_*.csv`, `U3_*.csv`, `U4_H3_band_check.txt`

**Definition (as coded):**  
\(\mathrm{leak}_h = -\beta_{W^{\mathrm{exc}},h}/\beta_{\mathrm{Aus},h}\) when \(\beta_{\mathrm{Aus}}<-0.01\);  
\(\mathrm{cum\text{-}leak}_H = -(\sum_{h=0}^{H}\beta_{W^{\mathrm{exc}},h})\,/\,(\sum_{h=0}^{H}\beta_{\mathrm{Aus},h})\).  
Numerator is the **LP of the \(W^{\mathrm{exc}}\) series**, not the sum of country LPs. Ratio recomputed inside each MBB draw.

### U1. Per-horizon leakage, keep-2022, \(h=0..6\)

| h | \(W^{\mathrm{exc}}\) Mt [90%] | Aus Mt [90%] | leakage [90%] |
|--:|------------------------------:|-------------:|--------------:|
| 0 | +1.486 [−4.33, +8.99] | −3.226 [−4.26, −1.25] | **0.461** [−1.70, +3.99] |
| 1 | **+6.359** [−2.53, +16.0] | −2.367 [−3.80, +1.07] | **2.687** [−0.86, +18.5] |
| 2 | +3.494 [−6.29, +13.0] | −1.537 [−4.48, +1.01] | 2.273 [−4.03, +27.2] |
| 3 | +3.949 [−6.57, +15.2] | +0.546 [−2.09, +2.68] | NaN (Aus ≥ 0) |
| 6 | +9.168 [−9.93, +16.4] | −1.516 [−3.64, +2.15] | 6.048 [−8.86, +44.7] |

Impact leakage ≈ **0.46** (band covers 0 and 1). The large \(h=1\) \(W^{\mathrm{exc}}\) point (+6.4 Mt) drives cumulative ratios upward; its band still covers 0.

### U2. Cumulative leakage, keep-2022, \(H=0..3\)

| H | Cum Aus | Cum \(W^{\mathrm{exc}}\) | cum leak [90%] |
|--:|--------:|-------------------------:|----------------|
| 0 | −3.226 | +1.486 | 0.461 [−1.80, +4.06] |
| 1 | −5.593 | +7.845 | 1.403 [−1.62, +7.92] |
| 2 | −7.130 | +11.339 | 1.590 [−1.94, +8.84] |
| 3 | −6.585 | +15.287 | **2.322** [**−5.36**, **+16.30**] |

### U3. Same cumulative, drop-2022 (robustness)

| H | Cum Aus | Cum \(W^{\mathrm{exc}}\) | cum leak [90%] |
|--:|--------:|-------------------------:|----------------|
| 0 | −3.283 | −0.526 | −0.160 [−2.82, +2.77] |
| 1 | −5.454 | +3.149 | 0.577 [−4.07, +6.39] |
| 2 | −7.063 | +3.830 | 0.542 [−3.76, +6.70] |
| 3 | −6.737 | +4.540 | **0.674** [−8.35, +12.12] |

Keeping 2022 raises the H=3 point from ~0.67 to ~2.32; bands remain enormous under both specs.

### U4. Does the H=3 cumulative band exclude 0? Exclude 0.73?

**No and no.** Keep-2022 H=3 band = [−5.36, +16.30] contains both **0** and **0.73**.

---

## Task V — Randomisation for horizons that matter [BLOCKING]

**Output:** `V_randomisation.csv` · 10 000 placebos · cyclone-season pool · keep-2022 \(z_w\) (10 events, pool 80).

| test | h | β | rand \(p\) |
|------|--:|--:|----------:|
| V1 \(W^{\mathrm{exc}}\) Mt | 1 | +6.359 | **0.196** |
| V1 \(W^{\mathrm{exc}}\) Mt | 2 | +3.494 | 0.410 |
| V1 \(W^{\mathrm{exc}}\) Mt | 3 | +3.949 | 0.312 |
| V2 Cum \(W^{\mathrm{exc}}\) through H=3 | 3 | +15.287 | **0.288** |
| V3 Australia Mt | 1 | −2.367 | **0.022** |
| V3 Australia Mt | 2 | −1.537 | 0.183 |
| V4 Japan from-Aus log | 0 | −0.053 | 0.209 |
| V5 Price (keep-2022) | 0 | +0.0649 | **0.235** |

**V5 confirms** Task R’s 0.235 under the new baseline control set (same estimate and \(p\)).  
**V1–V2:** the \(h=1\) / cumulative \(W^{\mathrm{exc}}\) cells that inflate U2 are **not** randomised-significant.  
**V3:** Aus Mt remains significant at \(h=1\) (\(p\approx0.02\)), not at \(h=2\).  
**V4:** Japan from-Aus impact is not randomised-significant under keep-2022.

---

## Task W — Which variables make \(z\) predictable?

**Script:** `diagnostics/task_W.jl` · **Outputs:** `W_blockF.csv`, `W_blockF.txt`  
3 lags each of: (a) lag \(z\), (b) \(\Delta\log q_{\mathrm{Aus}}\), (c) \(\Delta\log p\), (d) REA; month FE always retained. Block-F = unrestricted vs drop that block.

| instrument | block | F | p |
|------------|-------|--:|--:|
| z | joint a+b+c+d | 1.79 | 0.056 |
| z | (a) lag z | 1.23 | 0.300 |
| z | (b) Δlog q | 0.64 | 0.590 |
| z | **(c) Δlog p** | **4.40** | **0.006** |
| z | (d) REA | 0.66 | 0.581 |
| \(z_w\) | joint a+b+c+d | 2.48 | 0.006 |
| \(z_w\) | **(a) lag \(z_w\)** | **4.71** | **0.004** |
| \(z_w\) | (b) Δlog q | 1.01 | 0.390 |
| \(z_w\) | (c) Δlog p | 3.56 | 0.016 |
| \(z_w\) | (d) REA | 0.41 | 0.749 |

**Verdict**

- **Binary \(z\):** dominated by **(c) lagged price** → **exogeneity threat** (instrument predicted by recent coal-price moves).  
- **Preferred \(z_w\):** dominated by **(a) own lags** → primarily **event clustering** within cyclone seasons. Price block is still significant (\(p\approx0.016\)), so a **secondary** price channel remains; REA is irrelevant.  
- Overall for the paper’s \(z_w\): clustering is the main predictability story; do not ignore the mild price-lag association.

---

## Task X — Baseline table rebuild (keep-2022)

**Outputs:** `X_baseline_keep2022.csv`, `X_baseline_keep2022.md`  
This is the paper-facing table under the new baseline.

| object | h=0 | h=1 | h=2 | h=3 | h=6 | h=12 |
|--------|----:|----:|----:|----:|----:|-----:|
| Aus exports log | −0.116 [−0.15, −0.04] | −0.088 [−0.14, 0.04] | −0.051 [−0.15, 0.03] | 0.020 [−0.07, 0.09] | −0.048 [−0.12, 0.07] | −0.012 [−0.12, 0.11] |
| Aus exports Mt | −3.226 [−4.26, −1.25] | −2.367 [−3.80, 1.07] | −1.537 [−4.48, 1.01] | 0.546 [−2.09, 2.68] | −1.516 [−3.64, 2.15] | −0.484 [−3.44, 3.31] |
| Price (log real) | **0.065** [−0.14, 0.24] | 0.120 [−0.17, 0.44] | 0.114 [−0.21, 0.54] | 0.119 [−0.26, 0.58] | 0.212 [−0.34, 0.62] | 0.089 [−0.45, 0.49] |
| \(W^{\mathrm{exc}}\) Mt | 1.486 [−4.33, 8.99] | 6.359 [−2.53, 16.01] | 3.494 [−6.29, 12.96] | 3.949 [−6.57, 15.16] | 9.168 [−9.93, 16.40] | 2.110 [−16.27, 13.19] |
| Indonesia Mt | 1.560 [−3.75, 7.45] | 2.943 [−2.08, 9.58] | 1.319 [−4.12, 7.76] | 2.840 [−3.83, 9.77] | 2.658 [−5.01, 8.61] | −2.440 [−8.73, 7.08] |
| USA Mt | 0.282 [−2.32, 2.31] | 0.817 [−1.53, 3.69] | 0.586 [−1.98, 4.31] | 0.252 [−2.77, 3.78] | −0.262 [−3.73, 3.68] | 1.939 [−4.28, 4.33] |
| Colombia Mt | −0.566 [−3.02, 1.52] | **2.397** [0.05, 4.59] | 0.849 [−2.78, 2.57] | 0.308 [−2.11, 2.94] | 1.971 [−3.05, 3.55] | 1.318 [−3.29, 3.45] |
| Canada Mt | −0.095 [−0.55, 0.32] | 0.182 [−0.25, 0.70] | −0.110 [−0.57, 0.46] | −0.325 [−0.76, 0.27] | −0.081 [−0.61, 0.62] | −0.190 [−0.69, 0.68] |
| Japan total log | −0.023 [−0.15, 0.07] | 0.046 [−0.09, 0.13] | −0.011 [−0.14, 0.15] | 0.065 [−0.13, 0.19] | 0.062 [−0.13, 0.13] | −0.007 [−0.17, 0.15] |
| Japan from-Aus log | −0.053 [−0.15, 0.01] | −0.047 [−0.12, 0.05] | 0.051 [−0.11, 0.14] | 0.067 [−0.11, 0.16] | −0.002 [−0.11, 0.13] | −0.010 [−0.15, 0.14] |
| Japan non-Aus log | −0.036 [−0.19, 0.16] | 0.070 [−0.08, 0.25] | −0.070 [−0.18, 0.18] | 0.015 [−0.15, 0.26] | −0.024 [−0.18, 0.18] | 0.005 [−0.20, 0.20] |
| Japan total Mt | −0.275 [−2.18, 1.09] | 0.663 [−1.29, 1.92] | −0.355 [−2.11, 2.09] | 0.800 [−1.97, 2.67] | 0.807 [−1.96, 1.93] | −0.043 [−2.39, 2.22] |
| Japan from-Aus Mt | −0.479 [−1.36, 0.11] | −0.423 [−1.09, 0.49] | 0.425 [−1.04, 1.21] | 0.571 [−1.03, 1.45] | −0.039 [−1.04, 1.17] | −0.030 [−1.33, 1.28] |
| Japan non-Aus Mt | −0.209 [−1.07, 0.87] | 0.345 [−0.46, 1.32] | −0.432 [−1.02, 0.93] | 0.085 [−0.85, 1.41] | −0.106 [−0.99, 1.01] | 0.034 [−1.10, 1.12] |

(`Australia_mt` in the CSV duplicates Aus exports Mt.)

**Paper-facing contrasts vs drop-2022:** Aus impact still ≈ −3.2 Mt / −0.12 log with bands excluding 0; price impact is **positive** (+0.065) with band covering 0; \(W^{\mathrm{exc}}\) and Indonesia point toward offset but bands (and randomisation) do not pin them down; Colombia \(h=1\) is the only competitor cell whose 90% band excludes 0.

---

## Bottom line

1. **U:** Proper \(W^{\mathrm{exc}}\) leakage at impact ≈ **0.46**; H=3 cumulative ≈ **2.32** under keep-2022 — but the **90% band includes both 0 and 0.73**. Drop-2022 H=3 ≈ 0.67, same conclusion.  
2. **V:** Aus Mt \(h=1\) survives randomisation; \(W^{\mathrm{exc}}\) \(h=1\) / cum-H=3 and price \(h=0\) do **not**.  
3. **W:** For \(z_w\), predictability is mainly **own lags (clustering)**; binary \(z\) is mainly **price lags (exogeneity threat)**.  
4. **X:** Use `X_baseline_keep2022.csv` / `.md` as the keep-2022 paper numbers.
