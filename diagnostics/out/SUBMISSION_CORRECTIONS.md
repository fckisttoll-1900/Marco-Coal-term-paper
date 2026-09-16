# Minimum submission corrections (paste-ready)

Generated from clean rebuild + headline keep-2022 randomisation.
Do not edit Overleaf from this file automatically — paste manually.

## Pipeline status

- Modified: `01_build_panel.jl` (restores `z_w`: major=1, minor=0.5 from `events.csv`)
- Added: `diagnostics/task_headline_rand.jl`
- Rebuilt: `data/panel.csv` (now includes `z_w`)
- `load_panel()` succeeds: T=175 from 2012-01, sum(z_w)=6.0, N(z>0)=10 (8×0.5 + 2×1)
- Headline results: `diagnostics/out/headline_rand_keep2022.csv`, `.txt`
- Baseline point estimates on clean panel still round to the verified display values
  (Aus log −0.116; Aus Mt −3.23; Wexc +1.49; impact leakage 0.46).
- Aus Mt is a **direct levels** local projection on tonnes/1e6, not a log×mean conversion.

## Headline keep-2022 randomisation (10,000 draws; valid-row pool)

| Object | Units | Shock | Sample | h | β (clean) | rand p | pool |
|--------|-------|-------|--------|---|-----------|--------|------|
| Australian exports | log | z_w keep | keep-2022 | 0 | −0.116 | **0.0003** | 67 |
| Real Newcastle price | log | z_w keep | keep-2022 | 0 | +0.065 | **0.204** | 70 |
| Wexc | Mt | z_w keep | keep-2022 | 0 | +1.49 | **0.700** | 67 |

Procedure (`task_headline_rand.jl`, final):
Placebo event months are drawn from the November–April pool **restricted to valid
regression rows for that outcome and horizon**, while preserving the number of
**effective** treated months and their severity-weight multiset (here all ten
events enter h=0). Exact month-of-year composition and multi-month event
clustering are not preserved.

## Randomisation p-value audit (retain vs em dash)

### Table 1 — Australian exports (`tab:firststage`)

| Displayed | Verdict | Reason |
|-----------|---------|--------|
| h=0 p=0.003 | **REPLACE → 0.0003** | Was drop-2022 (8 events). Keep-2022 valid-row pool: 0.0003 |
| h=1 p=0.022 | **RETAIN** | Matches keep-2022 Aus **Mt** h=1 (V3: 0.0215, β=−2.37). Not a log test |
| h=2 p=0.183 | **RETAIN** | Matches keep-2022 Aus **Mt** h=2 (V3: 0.1827, β=−1.54) |
| h=6 p=0.313 | **REPLACE → —** | From drop-2022 Aus **log** (R3); keep log h=6 is −0.048 ≠ −0.040 |
| h=12 p=0.994 | **REPLACE → —** | From drop-2022 Aus log (R3); keep log h=12 is −0.012 ≠ 0.000 |

### Table 2 — Price (`tab:price`)

| Displayed | Verdict | Reason |
|-----------|---------|--------|
| h=0 p=0.235 | **REPLACE → 0.204** | Keep-2022 valid-row pool recomputed: 0.2038 (β=+0.065) |
| h=3 p=0.627 | **REPLACE → —** | From drop-2022 price h=3 (β=−0.071); keep is +0.119 |
| h=6 p=0.943 | **REPLACE → —** | From drop-2022 price h=6; keep is +0.212 |

### Table 3 — Replacement (`tab:leakage`)

| Displayed | Verdict | Reason |
|-----------|---------|--------|
| Wexc h=0 p=0.901 | **REPLACE → 0.700** | Was drop-2022 Wexc β=−0.53. Keep +1.49 valid-row pool → p=0.700 |
| Wexc h=1 p=0.196 | **RETAIN** | Matches keep-2022 Wexc h=1 (V1: 0.1956, β=+6.36) |
| Cum. leakage H=3 p=0.288 | **REPLACE → —** | p is for cum **Wexc** sum (15.29), not leakage rate 2.32 |

### Table 4 — Japan (`tab:japan`)

| Displayed | Verdict | Reason |
|-----------|---------|--------|
| from Aus p=0.209 | **RETAIN** | Matches keep-2022 Japan_from_aus_log h=0 (V4: 0.2094, β=−0.053) |

### Prose p-values to update

- “p=0.003” (Aus impact) → **0.0003** (keep-2022, valid-row pool)
- “p=0.235” (price impact) → **0.204**
- “p=0.901 on impact” (Wexc) → **0.700**
- Leave “p=0.196” (Wexc h=1) unchanged (not rerun)

## Figure 4 / REA claim — REMOVE

Do **not** rebuild Figure 4. Current “With REA” vs “No REA” contrast also changes
the sample/shock (drop-2022 vs the keep-2022 baseline). Recommend deleting
Figure 4 and the paragraph claiming omission of REA raises leakage to 0.93.

## Figure 2 note

`fig2_randomisation.pdf` was built from drop-2022 placebo draws (`H_placebo_draws.csv`,
n_events=8, old p≈0.003). After updating the keep-2022 impact p to 0.0003, redraw
Figure 2 from keep-2022 valid-row-pool placebos or drop the figure until redrawn.

## Bootstrap wording (verified)

`block_idx` draws overlapping 12-month blocks until the calendar length T is filled
(T=175 in the clean panel ≈ ceil(T/12)=15 blocks), then the LP keeps only rows with
finite outcomes/controls.

Recommended sentence:
“The bootstrap resamples the approximately 175-month calendar using about 15 draws
of 12-month moving blocks per replication, after which valid estimation rows are
selected.”
(If you prefer the prior “174” wording, it remains approximately correct.)

## Leave-one-event (no rerun)

C3 uses first-stage controls **with coal-price lags** (and full z_w), baseline
β≈−0.129; dropping April 2017 → −0.070. That is **not** the paper’s harmonised
baseline (−0.116, no price lags).

---

# Paste-ready Overleaf replacements

## A. Randomisation procedure (replace the paragraph that claims seasonal composition / out-of-season split)

```latex
We therefore report a randomisation test alongside every headline estimate.
Placebo event months are drawn from the November--April pool restricted to the
valid local-projection rows for that outcome and horizon, while preserving the
number of effective treated months and their severity-weight multiset; exact
month-of-year composition and multi-month event clustering are not preserved.
The test reports the share of 10{,}000 placebo draws producing a coefficient at
least as large in absolute value as the estimate. Where the bootstrap and the
randomisation test disagree, we say so.
```

## B. Table notes (Australian exports) — fix conversion claim + randomisation description

```latex
\item \textit{Notes:} Local projections as in equation~(\ref{eq:lp}). Brackets
are 90 percent moving-block bootstrap intervals (10{,}000 draws, 12-month
blocks). The Mt column is a direct levels local projection (tonnes), not a
rescaled log coefficient. Randomisation $p$-values use 10{,}000 placebo draws
from the November--April pool restricted to valid regression rows for that
outcome and horizon, preserving the multiset of effective treated-month severity
weights (not exact month-of-year composition or multi-month clustering).
```

## C. Table 1 cells (rand.\ $p$ column)

```latex
0  & $-0.116$                      & $-3.23$                       & $0.0003$ & $-3.28$ \\
   & \footnotesize$[-0.15,-0.04]$  & \footnotesize$[-4.26,-1.25]$  &         &         \\
1  & $-0.088$                      & $-2.37$                       & $0.022$ & $-2.17$ \\
   & \footnotesize$[-0.14,+0.04]$  & \footnotesize$[-3.80,+1.07]$  &         &         \\
2  & $-0.051$                      & $-1.54$                       & $0.183$ & $-1.61$ \\
   & \footnotesize$[-0.15,+0.03]$  & \footnotesize$[-4.48,+1.01]$  &         &         \\
3  & $+0.020$                      & $+0.55$                       &         & $+0.33$ \\
   & \footnotesize$[-0.07,+0.09]$  & \footnotesize$[-2.09,+2.68]$  &         &         \\
6  & $-0.048$                      & $-1.52$                       & ---     & $-1.29$ \\
12 & $-0.012$                      & $-0.48$                       & ---     & $-0.17$ \\
```

Optional clarity in notes if retaining h=1 and h=2 p-values:
``Randomisation $p$-values at $h=1,2$ are for the Mt specification.''

## D. First-stage prose (impact p + delete “survives every specification”)

```latex
Disruption months reduce Australian shipments by 11.6 log points on impact, or
about 3.2~Mt, with a 90 percent band of $[-0.15,-0.04]$ in logs and
$[-4.26,-1.25]$ in tonnes. The randomisation $p$-value is 0.0003: three in ten
thousand seasonal placebo assignments produce a contraction this large
(Figure~\ref{fig:randomisation}). The narrative measure is a relevant shifter of
Australian seaborne supply. The impact estimate is the one object in this design
estimated with enough precision to defend, but it is sensitive to the Cyclone
Debbie months: in a separate first-stage specification that also includes coal-price
lags, the coefficient is $-0.129$, and excluding April 2017 alone moves it to
$-0.070$.
```

## E. Price table — em dashes for misattributed horizons

```latex
0  & $+0.065$ \footnotesize$[-0.14,+0.24]$ & $0.204$ & $-0.003$ \\
1  & $+0.120$ \footnotesize$[-0.17,+0.44]$ &         &          \\
3  & $+0.119$ \footnotesize$[-0.26,+0.58]$ & ---     & $-0.059$ \\
6  & $+0.212$ \footnotesize$[-0.34,+0.62]$ & ---     & $+0.018$ \\
12 & $+0.089$ \footnotesize$[-0.45,+0.49]$ &         & $+0.081$ \\
```

## F. Leakage table

```latex
$W^{\mathrm{exc}}$, $h=0$ (Mt) & $+1.49$ \footnotesize$[-4.33,+8.99]$  & $0.700$ & $-0.53$ \\
$W^{\mathrm{exc}}$, $h=1$ (Mt) & $+6.36$ \footnotesize$[-2.53,+16.01]$ & $0.196$ & $+3.67$ \\
Impact leakage rate            & $0.46$  \footnotesize$[-1.70,+3.99]$  &         & $-0.16$ \\
Cumulative leakage, $H=3$      & $2.32$  \footnotesize$[-5.36,+16.30]$ & ---     & $0.67$  \\
```

And in the prose that cites competitor p-values:

```latex
Randomisation confirms the reading: the aggregate competitor response has
$p=0.700$ on impact and $p=0.196$ at $h=1$, so no horizon of the replacement
response is distinguishable from a seasonal placebo.
```

## G. Delete Figure 4 paragraph (REA / leakage 0.93)

Delete the paragraph beginning ``We note one specification result that bears on
how such estimates should be read. Omitting the global activity control...'' and
delete Figure~4 / any cross-reference claiming that omitting REA raises impact
leakage to 0.93. That contrast also changes the sample and shock definition
(drop-2022), so it is not an REA-only sensitivity.

Also delete parallel discussion-section sentences that repeat the 0.93 /
Indonesian $+2.0$~Mt REA-omission claim.

## H. Japan — delete “substitution gap”

```latex
Japanese imports from Australia fall by 5.3 log points on impact, imports from
other origins by 3.6, and total imports by 2.3. An Australia-only mechanical
benchmark using Japan's mean Australian import share of 66.3 percent would map
the Australian-origin response into a total-import response of
$0.663\times(-0.053)=-0.035$. The estimated total-import response is $-0.023$,
so the descriptive difference is $+0.012$ log points. This is not an accounting
identity: the total and origin-specific responses come from separate local
projections, and the bootstrap interval on the difference covers zero by a wide
margin. Excluding early 2022, the same descriptive comparison is $-0.032$ with
interval $[-0.125,+0.045]$, again covering zero, and in that specification
imports from other origins fall by more than imports from Australia.
```

## I. Discussion — impact p + leave-one wording

```latex
The impact estimate survives a randomisation test against November--April
placebo event sets at $p=0.0003$. Narrative identification of physical supply
interruptions works in this market, and the instrument measures what it claims
to measure.
```

```latex
Ten disruption months enter a sample of roughly 150 usable observations, and
only two carry full severity weight: March and April 2017, the two months of
Cyclone Debbie. In a separate first-stage sensitivity that includes coal-price
lags (unlike the baseline local projection), the coefficient is $-0.129$, and
dropping April 2017 alone moves it to $-0.070$. The first stage is, to a first
approximation, Cyclone Debbie.
```

Delete: ``survives every specification and inference procedure we ran.''

## J. Bootstrap (methods / discussion)

```latex
The bootstrap resamples the approximately 175-month calendar using about 15
draws of 12-month moving blocks per replication, after which valid estimation
rows are selected.
```

Replace claims of ``roughly fourteen blocks'' with ``about fifteen blocks''
where they refer to the moving-block scheme on the full post-2011 calendar.
