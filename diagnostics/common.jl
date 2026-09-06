# diagnostics/common.jl
# Shared helpers for diagnostics ONLY. Does not modify pipeline scripts.
# Conventions mirror 02_estimate.jl / 03_leakage.jl where noted.

using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random, Printf

const DIAG_ROOT = @__DIR__
const REPO      = dirname(DIAG_ROOT)
const DIAG_OUT  = joinpath(DIAG_ROOT, "out")
const RAW       = joinpath(REPO, "data", "raw")
const PANEL_PATH = joinpath(REPO, "data", "panel.csv")
const EVENTS_PATH = joinpath(REPO, "data", "events.csv")

const SAMPLE_START = "2012-01"
const NLAGS = 6
const BLOCK = 12
const BOOT_SEED = 20260906   # diagnostics bootstrap seed (reported in REPORT)
const NBOOT_B = 1000
const DROP_2022 = Set(["2022-02", "2022-03"])

mkpath(DIAG_OUT)

ym(x) = first(string(x), 7)

function monthdummies(dates)::Matrix{Float64}
    T = length(dates)
    D = zeros(T, 11)
    for (i, d) in enumerate(dates)
        m = month(d isa Date ? d : Date(ym(d) * "-01"))
        m > 1 && (D[i, m - 1] = 1.0)
    end
    return D
end

olsb(y, X) = hcat(ones(length(y)), X) \ y

function block_idx(T::Int, block::Int, rng)
    idx = Int[]
    while length(idx) < T
        s = rand(rng, 1:max(1, T - block + 1))
        append!(idx, s:min(s + block - 1, T))
    end
    return idx[1:T]
end

"Load panel from 2012-01 with log q, log real p, REA, three z variants."
function load_panel()
    panel = CSV.read(PANEL_PATH, DataFrame)
    panel.ds = ym.(panel.date)
    panel = panel[panel.ds .>= SAMPLE_START, :]
    sort!(panel, :date)
    dates = [d isa Date ? d : Date(ym(d) * "-01") for d in panel.date]
    lq = [ismissing(x) || x <= 0 ? NaN : log(Float64(x)) for x in panel.q_aus]
    lp = [ismissing(x) || x <= 0 ? NaN : log(Float64(x)) for x in panel.p_real]
    rea = [ismissing(x) ? NaN : Float64(x) for x in panel.rea]
    z_bin = Float64.(coalesce.(panel.z, 0.0))
    z_maj = Float64.(coalesce.(panel.z_major, 0.0))
    z_w   = Float64.(coalesce.(panel.z_w, 0.0))
    q_mt  = [ismissing(x) ? NaN : Float64(x) / 1e6 for x in panel.q_aus]  # Mt
    p_nom = [ismissing(x) || x <= 0 ? NaN : log(Float64(x)) for x in panel.p_usd]
    return (; panel, dates, ds = panel.ds, lq, lp, rea, z_bin, z_maj, z_w, q_mt, p_nom)
end

"""
Row set for first-stage-style LP at horizon h (mirrors 02_estimate.jl first_stage_F /
lpiv row filter): need y[t+h], y[t-1], z[t], lags of q and p, REA level+lags.
"""
function rows_fs(y, z, dates, q, p, rea, h::Int; nlags::Int = NLAGS)
    T = length(y)
    rows = Int[]
    for t in (nlags + 1):(T - h)
        ok = isfinite(y[t + h]) && isfinite(y[t - 1]) && isfinite(z[t])
        ok = ok && all(isfinite, @view q[(t - nlags):(t - 1)]) &&
                   all(isfinite, @view p[(t - nlags):(t - 1)])
        ok = ok && isfinite(rea[t]) && all(isfinite, @view rea[(t - nlags):(t - 1)])
        ok && push!(rows, t)
    end
    return rows
end

"Design matrix for FS-style LP: [z | q lags | p lags | rea | rea lags | month FE]."
function design_fs(rows, z, q, p, rea, dates; nlags::Int = NLAGS)
    D = monthdummies(dates)
    n = length(rows)
    k = 2 * nlags + (nlags + 1) + size(D, 2)
    X = zeros(n, 1 + k)
    X[:, 1] = z[rows]
    c = 2
    for l in 1:nlags
        X[:, c] = q[rows .- l]; c += 1
        X[:, c] = p[rows .- l]; c += 1
    end
    X[:, c] = rea[rows]; c += 1
    for l in 1:nlags
        X[:, c] = rea[rows .- l]; c += 1
    end
    X[:, c:end] = D[rows, :]
    return X
end

"""
Estimate FS-style LP coefficient on z at horizon h.
Returns (beta, se_ols, n, df) or NaNs if singular / too few obs.
"""
function lp_fs(y, z, dates, q, p, rea, h::Int; nlags::Int = NLAGS)
    rows = rows_fs(y, z, dates, q, p, rea, h; nlags = nlags)
    (length(rows) < 40 || sum(abs, z[rows]) < 1.5) && return (NaN, NaN, length(rows), NaN)
    X = design_fs(rows, z, q, p, rea, dates; nlags = nlags)
    yy = [y[t + h] - y[t - 1] for t in rows]
    try
        Z = hcat(ones(length(rows)), X)
        b = Z \ yy
        e = yy - Z * b
        k = size(Z, 2)
        df = length(rows) - k
        df <= 0 && return (NaN, NaN, length(rows), Float64(df))
        s2 = sum(e.^2) / df
        se = sqrt(s2 * inv(Z' * Z)[2, 2])
        return (b[2], se, length(rows), Float64(df))
    catch
        return (NaN, NaN, length(rows), NaN)
    end
end

"Newey-West SE for the shock coefficient (HAC on score of β_z)."
function nw_se_fs(y, z, dates, q, p, rea, h::Int; nlags::Int = NLAGS, bandwidth::Int = 1)
    rows = rows_fs(y, z, dates, q, p, rea, h; nlags = nlags)
    (length(rows) < 40 || sum(abs, z[rows]) < 1.5) && return NaN
    X = design_fs(rows, z, q, p, rea, dates; nlags = nlags)
    yy = [y[t + h] - y[t - 1] for t in rows]
    Z = hcat(ones(length(rows)), X)
    try
        b = Z \ yy
        e = yy - Z * b
        n = length(rows)
        # bread
        Qinv = inv(Z' * Z / n)
        # meat: NW on Z_i * e_i
        S = zeros(size(Z, 2), size(Z, 2))
        scores = [Z[i, :] * e[i] for i in 1:n]
        for i in 1:n
            S .+= scores[i] * scores[i]'
        end
        L = bandwidth
        for ell in 1:L
            w = 1.0 - ell / (L + 1)
            Gam = zeros(size(Z, 2), size(Z, 2))
            for i in (ell + 1):n
                Gam .+= scores[i] * scores[i - ell]'
            end
            S .+= w * (Gam + Gam')
        end
        S ./= n
        V = Qinv * S * Qinv / n
        return sqrt(V[2, 2])
    catch
        return NaN
    end
end

"""
Simple LP as in 03_leakage: y on z + own lags + month FE (no REA, no p lags).
"""
function rows_simple(y, z, h::Int; nlags::Int = NLAGS)
    T = length(y)
    rows = Int[]
    for t in (nlags + 1):(T - h)
        ok = isfinite(y[t + h]) && isfinite(y[t - 1]) && isfinite(z[t])
        ok && (ok = all(isfinite, @view y[(t - nlags):(t - 1)]))
        ok && push!(rows, t)
    end
    return rows
end

function lp_simple(y, z, dates, h::Int; nlags::Int = NLAGS)
    rows = rows_simple(y, z, h; nlags = nlags)
    (length(rows) < 40 || sum(z[rows]) < 2) && return (NaN, NaN, length(rows), NaN)
    D = monthdummies(dates)
    L = zeros(length(rows), nlags)
    for (i, t) in enumerate(rows), l in 1:nlags
        L[i, l] = y[t - l]
    end
    X = hcat(z[rows], L, D[rows, :])
    yy = [y[t + h] - y[t - 1] for t in rows]
    try
        Z = hcat(ones(length(rows)), X)
        b = Z \ yy
        e = yy - Z * b
        df = length(rows) - size(Z, 2)
        df <= 0 && return (NaN, NaN, length(rows), Float64(df))
        s2 = sum(e.^2) / df
        se = sqrt(s2 * inv(Z' * Z)[2, 2])
        return (b[2], se, length(rows), Float64(df))
    catch
        return (NaN, NaN, length(rows), NaN)
    end
end

function apply_drop!(z::Vector{Float64}, ds, dropset)
    zout = copy(z)
    for i in eachindex(zout)
        string(ds[i]) in dropset && (zout[i] = 0.0)
    end
    return zout
end

println("diagnostics/common.jl loaded; BOOT_SEED=$BOOT_SEED")
