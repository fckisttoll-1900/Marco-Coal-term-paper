# diagnostics/task_S.jl — Unpredictability with fewer lags
include(joinpath(@__DIR__, "harm_lp.jl"))

P = load_panel()
T = length(P.dates)
Dfe = monthdummies(P.dates)

has_dist = try
    @eval Main using Distributions
    true
catch
    false
end

function ftest_z(z, nlags::Int)
    # need t-nlags-1 for Δ at lag nlags
    rows2 = Int[]
    for t in (nlags + 2):T
        ok = isfinite(z[t])
        for l in 1:nlags
            ok = ok && isfinite(P.lq[t - l]) && isfinite(P.lq[t - l - 1])
            ok = ok && isfinite(P.lp[t - l]) && isfinite(P.lp[t - l - 1])
            ok = ok && isfinite(P.rea[t - l])
        end
        ok && push!(rows2, t)
    end
    n = length(rows2)
    nreg = 3 * nlags   # Δq, Δp, REA lags
    X = zeros(n, nreg + size(Dfe, 2))
    for (i, t) in enumerate(rows2)
        c = 1
        for l in 1:nlags
            X[i, c] = P.lq[t - l] - P.lq[t - l - 1]; c += 1
        end
        for l in 1:nlags
            X[i, c] = P.lp[t - l] - P.lp[t - l - 1]; c += 1
        end
        for l in 1:nlags
            X[i, c] = P.rea[t - l]; c += 1
        end
        X[i, c:end] = Dfe[t, :]
    end
    y = z[rows2]
    Z = hcat(ones(n), X)
    e = y - Z * (Z \ y)
    k = size(Z, 2)
    Xr = hcat(ones(n), Dfe[rows2, :])
    er = y - Xr * (Xr \ y)
    SSR_r = sum(er.^2); SSR_u = sum(e.^2); q = nreg
    F = ((SSR_r - SSR_u) / q) / (SSR_u / (n - k))
    pval = has_dist ? (1 - cdf(Main.FDist(q, n - k), F)) : NaN
    R2 = 1 - SSR_u / sum((y .- mean(y)).^2)
    return (F, pval, R2, n, q)
end

out = DataFrame(instrument = String[], nlags = Int[], n_regressors = Int[],
                F = Float64[], p = Float64[], R2 = Float64[], N = Int[])
for nl in (6, 3, 1)
    for (name, z) in (("z", P.z_bin), ("z_w", P.z_w))
        F, p, R2, n, q = ftest_z(z, nl)
        push!(out, (name, nl, q, F, p, R2, n))
    end
end
CSV.write(joinpath(DIAG_OUT, "S_unpredictability_lags.csv"), out)

open(joinpath(DIAG_OUT, "S_interpretation.txt"), "w") do io
    println(io, out)
    println(io, "")
    println(io, "Interpretation:")
    r6 = out[(out.instrument .== "z") .& (out.nlags .== 6), :]
    r3 = out[(out.instrument .== "z") .& (out.nlags .== 3), :]
    r1 = out[(out.instrument .== "z") .& (out.nlags .== 1), :]
    @printf(io, "z: 6 lags p=%.4f; 3 lags p=%.4f; 1 lag p=%.4f\n",
            r6.p[1], r3.p[1], r1.p[1])
    # If fewer lags make p larger (less significant), marginal 6-lag result looks like overfitting.
    # If fewer lags make p smaller/stay small, more like genuine predictability.
    if r3.p[1] > 0.10 && r6.p[1] < 0.10
        println(io, "VERDICT: marginal p≈0.08 at 6 lags looks like OVERFITTING —")
        println(io, "significance disappears (or weakens past 10%) with 3 lags.")
    elseif r1.p[1] < 0.05 || r3.p[1] < 0.05
        println(io, "VERDICT: predictability SURVIVES fewer lags — more consistent with")
        println(io, "genuine forecastability than pure overfitting.")
    else
        println(io, "VERDICT: across lag lengths, evidence of predictability is WEAK/MARGINAL;")
        println(io, "the 6-lag p≈0.083 should not be read as strong rejection of exogeneity.")
    end
end
println("wrote S_*")
println("script: diagnostics/task_S.jl")
