# One-shot: h=0 bootstrap percentiles for price and Aus export LP
# (same controls as 02_estimate first stage: month FE, REA, z_w, drop 2022-02/03)

using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random, Printf

Random.seed!(20260905)

const NLAGS = 6
const NBOOT = 2000
const BLOCK = 12
const SAMPLE_START = "2012-01"
const DROP = Set(["2022-02", "2022-03"])
const ROOT = @__DIR__

olsb(y, X) = hcat(ones(length(y)), X) \ y
ym(x) = first(string(x), 7)

function monthdummies(dates)
    T = length(dates); D = zeros(T, 11)
    for (i, d) in enumerate(dates)
        m = month(d isa Date ? d : Date(ym(d) * "-01"))
        m > 1 && (D[i, m - 1] = 1.0)
    end
    return D
end

function block_idx(T, block, rng)
    idx = Int[]
    while length(idx) < T
        s = rand(rng, 1:max(1, T - block + 1))
        append!(idx, s:min(s + block - 1, T))
    end
    return idx[1:T]
end

"""LP of Δy on z with lags of y, of q, of p, REA, month FE — match first_stage_F controls."""
function lp_h0(y, z, dates, q, p, rea)
    T = length(y); D = monthdummies(dates)
    rows = Int[]
    for t in (NLAGS + 1):T
        ok = isfinite(y[t]) && isfinite(y[t - 1]) && isfinite(z[t])
        ok = ok && all(isfinite, @view q[(t - NLAGS):(t - 1)]) &&
                   all(isfinite, @view p[(t - NLAGS):(t - 1)])
        ok = ok && isfinite(rea[t]) && all(isfinite, @view rea[(t - NLAGS):(t - 1)])
        ok && push!(rows, t)
    end
    length(rows) < 40 && return (NaN, NaN, NaN, Int[])
    n = length(rows)
    k = 2 * NLAGS + (NLAGS + 1) + size(D, 2)
    X = zeros(n, 1 + k)
    X[:, 1] = z[rows]
    c = 2
    for l in 1:NLAGS
        X[:, c] = q[rows .- l]; c += 1
        X[:, c] = p[rows .- l]; c += 1
    end
    X[:, c] = rea[rows]; c += 1
    for l in 1:NLAGS
        X[:, c] = rea[rows .- l]; c += 1
    end
    X[:, c:end] = D[rows, :]
    yy = [y[t] - y[t - 1] for t in rows]
    try
        b = olsb(yy, X)
        e = yy - hcat(ones(n), X) * b
        s2 = sum(e.^2) / (n - size(X, 2) - 1)
        XtX = hcat(ones(n), X)' * hcat(ones(n), X)
        se = sqrt(s2 * inv(XtX)[2, 2])
        return b[2], b[2] / se, (b[2] / se)^2, rows
    catch
        return (NaN, NaN, NaN, rows)
    end
end

panel = CSV.read(joinpath(ROOT, "data", "panel.csv"), DataFrame)
panel.ds = ym.(panel.date)
panel = panel[panel.ds .>= SAMPLE_START, :]
z = Float64.(coalesce.(panel.z_w, 0.0))
for (i, s) in enumerate(panel.ds)
    s in DROP && (z[i] = 0.0)
end
dates = [Date(s * "-01") for s in panel.ds]
lq = log.(Float64.(coalesce.(panel.q_aus, NaN)))
# real price if present
pcol = :p_real in propertynames(panel) ? :p_real : :p_usd
lp = log.(Float64.(coalesce.(panel[!, pcol], NaN)))
rea = Float64.(coalesce.(panel.rea, NaN))

bq, tq, Fq, _ = lp_h0(lq, z, dates, lq, lp, rea)
bp, tp, Fp, _ = lp_h0(lp, z, dates, lq, lp, rea)
@printf("point  q_aus β=%+.4f t=%+.2f F=%.2f\n", bq, tq, Fq)
@printf("point  price β=%+.4f t=%+.2f F=%.2f\n", bp, tp, Fp)
@printf("OLS SE(price)=%.4f   90%% normal approx=[%+.3f,%+.3f]\n",
        abs(bp / tp), bp - 1.645 * abs(bp / tp), bp + 1.645 * abs(bp / tp))
@printf("MDE (10%%,80%%)=%.3f\n", (1.645 + 0.842) * abs(bp / tp))

T = length(z)
Bq = fill(NaN, NBOOT); Bp = fill(NaN, NBOOT)
rng = MersenneTwister(20260905)

# --- moving-block bootstrap (may drop draws with singular FE) ---
print("block bootstrap $NBOOT ... ")
n_ok = 0
for b in 1:NBOOT
    idx = block_idx(T, BLOCK, rng)
    zb = z[idx]; lqb = lq[idx]; lpb = lp[idx]; reab = rea[idx]; db = dates[idx]
    sum(abs, zb) < 1.5 && continue
    Bq[b], _, _, _ = lp_h0(lqb, zb, db, lqb, lpb, reab)
    Bp[b], _, _, _ = lp_h0(lpb, zb, db, lqb, lpb, reab)
    if isfinite(Bq[b]) && isfinite(Bp[b])
        global n_ok += 1
    end
end
println("done ($n_ok joint finite)")

cq = filter(isfinite, Bq); cp = filter(isfinite, Bp)
@printf("boot n finite: q=%d price=%d\n", length(cq), length(cp))
@printf("BLOCK q_aus  h=0: point=%+.4f  90%% [%+.4f, %+.4f]  p50=%+.4f\n",
        bq, quantile(cq, 0.05), quantile(cq, 0.95), quantile(cq, 0.50))
@printf("BLOCK price  h=0: point=%+.4f  90%% [%+.4f, %+.4f]  p50=%+.4f\n",
        bp, quantile(cp, 0.05), quantile(cp, 0.95), quantile(cp, 0.50))

# --- wild bootstrap of residuals (design fixed; closer to 02_estimate) ---
function wild_beta(y, z, dates, q, p, rea, nboot, rng)
    beta0, _, _, rows = lp_h0(y, z, dates, q, p, rea)
    !isfinite(beta0) && return fill(NaN, nboot), NaN
    Tloc = length(y); D = monthdummies(dates)
    n = length(rows)
    k = 2 * NLAGS + (NLAGS + 1) + size(D, 2)
    X = zeros(n, 1 + k)
    X[:, 1] = z[rows]
    c = 2
    for l in 1:NLAGS
        X[:, c] = q[rows .- l]; c += 1
        X[:, c] = p[rows .- l]; c += 1
    end
    X[:, c] = rea[rows]; c += 1
    for l in 1:NLAGS
        X[:, c] = rea[rows .- l]; c += 1
    end
    X[:, c:end] = D[rows, :]
    Z = hcat(ones(n), X)
    yy = [y[t] - y[t - 1] for t in rows]
    bhat = Z \ yy
    e = yy - Z * bhat
    out = fill(NaN, nboot)
    for b in 1:nboot
        yb = Z * bhat + e .* rand(rng, [-1.0, 1.0], n)
        out[b] = (Z \ yb)[2]
    end
    return out, beta0
end

rng2 = MersenneTwister(20260905)
wq, _ = wild_beta(lq, z, dates, lq, lp, rea, NBOOT, rng2)
wp, _ = wild_beta(lp, z, dates, lq, lp, rea, NBOOT, rng2)
wq = filter(isfinite, wq); wp = filter(isfinite, wp)
@printf("WILD  q_aus  h=0: point=%+.4f  90%% [%+.4f, %+.4f]\n",
        bq, quantile(wq, 0.05), quantile(wq, 0.95))
@printf("WILD  price  h=0: point=%+.4f  90%% [%+.4f, %+.4f]\n",
        bp, quantile(wp, 0.05), quantile(wp, 0.95))

# exporter pipeline Australia Mt for substitution into text
exp = CSV.read(joinpath(ROOT, "out", "exporter_responses.csv"), DataFrame)
a0 = exp[(exp.exporter .== "Australia") .& (exp.horizon .== 0), :]
@printf("\nexporter pipeline Australia Mt h=0: %.3f [%.3f, %.3f]\n",
        a0.mt[1], a0.lo90_mt[1], a0.hi90_mt[1])
@printf("exporter pipeline Australia log β: %.4f\n", a0.beta[1])

# IEA thermal share
aus_th = 209.0; world_th = 1180.0
@printf("\nIEA Mid-Year 2025: Aus thermal exports %.0f Mt / world thermal ~%.0f Mt => s=%.3f\n",
        aus_th, world_th, aus_th / world_th)
