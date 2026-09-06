# diagnostics/task_W.jl — Which regressor blocks make z predictable?
# 3 lags; block-F for (a) lagged z, (b) Δlog Aus, (c) Δlog price, (d) REA.
include(joinpath(@__DIR__, "harm_lp.jl"))

try
    using Distributions
catch
    @warn "Distributions.jl unavailable; p-values will be NaN"
end

P = load_panel()
T = length(P.dates)
Dfe = monthdummies(P.dates)
const NL = 3

function fdist_cdf(F, d1, d2)
    try
        return cdf(FDist(d1, d2), F)
    catch
        return NaN
    end
end

"""
Build unrestricted model with all four blocks + month FE.
Returns rows, y, Xu (with intercept), and column index ranges for each block
(1-indexed into Xu after intercept column 1).
"""
function build_design(z)
    rows = Int[]
    for t in (NL + 2):T
        ok = isfinite(z[t])
        for l in 1:NL
            ok = ok && isfinite(z[t - l])
            ok = ok && isfinite(P.lq[t - l]) && isfinite(P.lq[t - l - 1])
            ok = ok && isfinite(P.lp[t - l]) && isfinite(P.lp[t - l - 1])
            ok = ok && isfinite(P.rea[t - l])
        end
        ok && push!(rows, t)
    end
    n = length(rows)
    # blocks: a=NL, b=NL, c=NL, d=NL, FE=11
    n_a = NL; n_b = NL; n_c = NL; n_d = NL; n_fe = size(Dfe, 2)
    X = zeros(n, n_a + n_b + n_c + n_d + n_fe)
    for (i, t) in enumerate(rows)
        c = 1
        for l in 1:NL
            X[i, c] = z[t - l]; c += 1
        end
        for l in 1:NL
            X[i, c] = P.lq[t - l] - P.lq[t - l - 1]; c += 1
        end
        for l in 1:NL
            X[i, c] = P.lp[t - l] - P.lp[t - l - 1]; c += 1
        end
        for l in 1:NL
            X[i, c] = P.rea[t - l]; c += 1
        end
        X[i, c:end] = Dfe[t, :]
    end
    y = z[rows]
    Xu = hcat(ones(n), X)
    # column ranges in Xu (excluding intercept at 1)
    ia = 2:(1 + n_a)
    ib = (2 + n_a):(1 + n_a + n_b)
    ic = (2 + n_a + n_b):(1 + n_a + n_b + n_c)
    id = (2 + n_a + n_b + n_c):(1 + n_a + n_b + n_c + n_d)
    # FE = remaining
    return y, Xu, n, Dict("a_lag_z" => collect(ia),
                          "b_dlog_q" => collect(ib),
                          "c_dlog_p" => collect(ic),
                          "d_rea" => collect(id))
end

function block_F(y, Xu, drop_cols)
    n, k_u = size(Xu)
    bu = Xu \ y
    eu = y - Xu * bu
    SSR_u = sum(eu.^2)
    keep = setdiff(1:k_u, drop_cols)
    Xr = Xu[:, keep]
    br = Xr \ y
    er = y - Xr * br
    SSR_r = sum(er.^2)
    q = length(drop_cols)
    df_u = n - k_u
    F = ((SSR_r - SSR_u) / q) / (SSR_u / df_u)
    p = 1 - fdist_cdf(F, q, df_u)
    return F, p, q, n, df_u, SSR_u, SSR_r
end

out = DataFrame(instrument = String[], block = String[],
                F = Float64[], p = Float64[], q = Int[], N = Int[], df = Int[])

open(joinpath(DIAG_OUT, "W_blockF.txt"), "w") do io
    println(io, "Task W: block-F tests at 3 lags + month FE always retained.")
    println(io, "Unrestricted: lag z + Δlog q + Δlog p + REA lags (3 each) + FE.")
    println(io, "")
    for (name, z) in (("z", P.z_bin), ("z_w", P.z_w))
        y, Xu, n, blocks = build_design(z)
        # joint on all economic + lag-z blocks (exclude FE+intercept)
        all_drop = vcat(blocks["a_lag_z"], blocks["b_dlog_q"],
                        blocks["c_dlog_p"], blocks["d_rea"])
        Fj, pj, qj, _, dfj, _, _ = block_F(y, Xu, all_drop)
        @printf(io, "=== %s  N=%d ===\n", name, n)
        @printf(io, "Joint (a+b+c+d vs FE only): F=%.4f p=%.4f q=%d\n", Fj, pj, qj)
        push!(out, (name, "joint_abcd", Fj, pj, qj, n, dfj))
        for (bname, cols) in (("a_lag_z", blocks["a_lag_z"]),
                              ("b_dlog_q", blocks["b_dlog_q"]),
                              ("c_dlog_p", blocks["c_dlog_p"]),
                              ("d_rea", blocks["d_rea"]))
            F, p, q, _, df, _, _ = block_F(y, Xu, cols)
            @printf(io, "  %-10s  F=%.4f  p=%.4f  q=%d\n", bname, F, p, q)
            push!(out, (name, bname, F, p, q, n, df))
        end
        println(io, "")
    end
    # verdict
    println(io, "Verdict rule: if (a) dominates → cyclone-season event clustering;")
    println(io, "if (c) or (d) dominate → exogeneity threat via price/REA.")
    for name in ("z", "z_w")
        sub = out[out.instrument .== name, :]
        blocks = ["a_lag_z", "b_dlog_q", "c_dlog_p", "d_rea"]
        ps = Dict(b => sub[sub.block .== b, :].p[1] for b in blocks)
        Fs = Dict(b => sub[sub.block .== b, :].F[1] for b in blocks)
        best = first(sort(collect(ps); by = x -> x[2]))[1]  # lowest p
        @printf(io, "%s: strongest block = %s (F=%.3f p=%.4f)\n",
                name, best, Fs[best], ps[best])
        if best == "a_lag_z" && ps["a_lag_z"] <= min(ps["c_dlog_p"], ps["d_rea"])
            println(io, "  → dominated by lagged z → event clustering, not price/REA threat.")
        elseif ps["c_dlog_p"] < 0.10 || ps["d_rea"] < 0.10
            println(io, "  → price and/or REA matter → exogeneity concern.")
        else
            println(io, "  → no clear exogeneity threat from (c)/(d); check (a)/(b).")
        end
    end
end

CSV.write(joinpath(DIAG_OUT, "W_blockF.csv"), out)
println("wrote W_blockF.*")
println("script: diagnostics/task_W.jl")
