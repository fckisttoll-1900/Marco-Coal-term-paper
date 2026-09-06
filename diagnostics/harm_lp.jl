# diagnostics/harm_lp.jl
# Harmonised LP: z_w, drop 2022-02/03, 6 own lags, REA+6 lags, month FE.
# Used by Tasks M, G, H, I, J, L. Does not modify pipeline scripts.

include(joinpath(@__DIR__, "common.jl"))

const NBOOT_M = 10_000
const NBOOT_H = 10_000

"""
Rows for harmonised LP at horizon h:
  y[t+h], y[t-1], z[t], y lags 1..nlags, rea[t] and rea lags all finite.
"""
function rows_harm(y, z, rea, h::Int; nlags::Int = NLAGS)
    T = length(y)
    rows = Int[]
    for t in (nlags + 1):(T - h)
        ok = isfinite(y[t + h]) && isfinite(y[t - 1]) && isfinite(z[t])
        ok = ok && all(isfinite, @view y[(t - nlags):(t - 1)])
        ok = ok && isfinite(rea[t]) && all(isfinite, @view rea[(t - nlags):(t - 1)])
        ok && push!(rows, t)
    end
    return rows
end

function design_harm(rows, z, y, rea, dates; nlags::Int = NLAGS)
    D = monthdummies(dates)
    n = length(rows)
    k = nlags + (nlags + 1) + size(D, 2)   # own lags + rea + rea lags + FE
    X = zeros(n, 1 + k)
    X[:, 1] = z[rows]
    c = 2
    for l in 1:nlags
        X[:, c] = y[rows .- l]; c += 1
    end
    X[:, c] = rea[rows]; c += 1
    for l in 1:nlags
        X[:, c] = rea[rows .- l]; c += 1
    end
    X[:, c:end] = D[rows, :]
    return X
end

"""
Harmonised LP coefficient path for horizons 0..H.
Returns Vector of length H+1 (NaN if inestimable).
"""
function lp_harm_path(y, z, rea, dates, H::Int; nlags::Int = NLAGS)
    β = fill(NaN, H + 1)
    for h in 0:H
        rows = rows_harm(y, z, rea, h; nlags = nlags)
        (length(rows) < 40 || sum(abs, z[rows]) < 1.5) && continue
        X = design_harm(rows, z, y, rea, dates; nlags = nlags)
        yy = [y[t + h] - y[t - 1] for t in rows]
        try
            β[h + 1] = olsb(yy, X)[2]
        catch
        end
    end
    return β
end

function lp_harm_h(y, z, rea, dates, h::Int; nlags::Int = NLAGS)
    rows = rows_harm(y, z, rea, h; nlags = nlags)
    (length(rows) < 40 || sum(abs, z[rows]) < 1.5) && return (NaN, NaN, 0)
    X = design_harm(rows, z, y, rea, dates; nlags = nlags)
    yy = [y[t + h] - y[t - 1] for t in rows]
    try
        Z = hcat(ones(length(rows)), X)
        b = Z \ yy
        e = yy - Z * b
        df = length(rows) - size(Z, 2)
        s2 = sum(e.^2) / max(df, 1)
        se = sqrt(s2 * inv(Z' * Z)[2, 2])
        return (b[2], se, length(rows))
    catch
        return (NaN, NaN, length(rows))
    end
end

"Legacy preferred instrument: z_w with 2022-02/03 zeroed (pre-REPORT4)."
function z_preferred(P)
    apply_drop!(P.z_w, P.ds, DROP_2022)
end

"New baseline (REPORT4+): z_w with 2022 RETAINED."
function z_baseline(P)
    copy(P.z_w)
end

"""
Align a (ds => tonnes) dict onto panel calendar; return tonnes vector (NaN gaps).
"""
function align_to_panel(tonnes_by_ds::Dict, ds_panel)
    [haskey(tonnes_by_ds, string(s)) ? Float64(tonnes_by_ds[string(s)]) : NaN
     for s in ds_panel]
end

println("harm_lp.jl loaded; NBOOT_M=$NBOOT_M NBOOT_H=$NBOOT_H")
