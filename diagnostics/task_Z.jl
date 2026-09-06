# diagnostics/task_Z.jl — Price-block predictability without mine incidents
# Same design as Task W (3 lags), comparing full z_w vs weather-only (mines zeroed).
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
const MINE = Set(["2020-05", "2024-06"])

function fdist_cdf(F, d1, d2)
    try
        return cdf(FDist(d1, d2), F)
    catch
        return NaN
    end
end

function drop_mines(z)
    zout = copy(z)
    for i in eachindex(zout)
        string(P.ds[i]) in MINE && (zout[i] = 0.0)
    end
    return zout
end

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
    ia = 2:(1 + n_a)
    ib = (2 + n_a):(1 + n_a + n_b)
    ic = (2 + n_a + n_b):(1 + n_a + n_b + n_c)
    id = (2 + n_a + n_b + n_c):(1 + n_a + n_b + n_c + n_d)
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
    return F, p, q, n, df_u
end

out = DataFrame(instrument = String[], variant = String[], block = String[],
                F = Float64[], p = Float64[], q = Int[], N = Int[], df = Int[],
                n_events = Int[])

open(joinpath(DIAG_OUT, "Z_mine_drop_blockF.txt"), "w") do io
    println(io, "Task Z: Task-W block-F with/without Grosvenor mine months (2020-05, 2024-06).")
    println(io, "Focus: does lagged-price block (c) lose significance for z_w?")
    println(io, "")
    for (name, z0) in (("z", P.z_bin), ("z_w", P.z_w))
        for (vname, z) in (("full", z0), ("drop_mines", drop_mines(z0)))
            n_ev = count(>(0), z)
            y, Xu, n, blocks = build_design(z)
            all_drop = vcat(blocks["a_lag_z"], blocks["b_dlog_q"],
                            blocks["c_dlog_p"], blocks["d_rea"])
            @printf(io, "=== %s / %s  N=%d  events_gt0=%d ===\n", name, vname, n, n_ev)
            for (bname, cols) in (("joint_abcd", all_drop),
                                  ("a_lag_z", blocks["a_lag_z"]),
                                  ("b_dlog_q", blocks["b_dlog_q"]),
                                  ("c_dlog_p", blocks["c_dlog_p"]),
                                  ("d_rea", blocks["d_rea"]))
                F, p, q, _, df = block_F(y, Xu, cols)
                @printf(io, "  %-12s  F=%.4f  p=%.4f\n", bname, F, p)
                push!(out, (name, vname, bname, F, p, q, n, df, n_ev))
            end
            println(io, "")
        end
    end
    # Decision on z_w price block
    p_full = out[(out.instrument .== "z_w") .& (out.variant .== "full") .&
                 (out.block .== "c_dlog_p"), :].p[1]
    p_drop = out[(out.instrument .== "z_w") .& (out.variant .== "drop_mines") .&
                 (out.block .== "c_dlog_p"), :].p[1]
    F_full = out[(out.instrument .== "z_w") .& (out.variant .== "full") .&
                 (out.block .== "c_dlog_p"), :].F[1]
    F_drop = out[(out.instrument .== "z_w") .& (out.variant .== "drop_mines") .&
                 (out.block .== "c_dlog_p"), :].F[1]
    @printf(io, "z_w price block: full F=%.4f p=%.4f; drop_mines F=%.4f p=%.4f\n",
            F_full, p_full, F_drop, p_drop)
    # "loses significance" = was <0.05 (or <0.10) and becomes >=0.10, or clearly weaker past 5%
    # Paper Task W reported p=0.016; use conventional 5% and also report 10%.
    survives_5 = p_drop < 0.05
    survives_10 = p_drop < 0.10
    if !survives_10
        println(io, "VERDICT: price block LOSES significance without mines (p>=0.10).")
        println(io, "KEEP SENTENCE (i).")
    else
        println(io, "VERDICT: price block SURVIVES mine removal (still p<0.10",
                survives_5 ? " and p<0.05" : " but p>=0.05", ").")
        println(io, "KEEP SENTENCE (ii).")
    end
end

CSV.write(joinpath(DIAG_OUT, "Z_mine_drop_blockF.csv"), out)
println("wrote Z_mine_drop_blockF.*")
println("script: diagnostics/task_Z.jl")
