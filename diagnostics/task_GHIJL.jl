# diagnostics/task_GHIJL.jl — G,H,I,J,L under HARMONISED specification
include(joinpath(@__DIR__, "harm_lp.jl"))

P = load_panel()
z_pref = z_preferred(P)
T = length(P.dates)
Dfe = monthdummies(P.dates)

has_dist = try
    @eval Main using Distributions
    true
catch
    false
end

# ============================ Task G ============================
println("G: instrument unpredictability...")
function ftest_z(z)
    rows2 = Int[]
    for t in (NLAGS + 2):T
        ok = isfinite(z[t])
        for l in 1:NLAGS
            ok = ok && isfinite(P.lq[t - l]) && isfinite(P.lq[t - l - 1])
            ok = ok && isfinite(P.lp[t - l]) && isfinite(P.lp[t - l - 1])
            ok = ok && isfinite(P.rea[t - l])
        end
        ok && push!(rows2, t)
    end
    n = length(rows2)
    X = zeros(n, 18 + size(Dfe, 2))
    for (i, t) in enumerate(rows2)
        c = 1
        for l in 1:NLAGS
            X[i, c] = P.lq[t - l] - P.lq[t - l - 1]; c += 1
        end
        for l in 1:NLAGS
            X[i, c] = P.lp[t - l] - P.lp[t - l - 1]; c += 1
        end
        for l in 1:NLAGS
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
    SSR_r = sum(er.^2); SSR_u = sum(e.^2); q = 18
    F = ((SSR_r - SSR_u) / q) / (SSR_u / (n - k))
    pval = has_dist ? (1 - cdf(Main.FDist(q, n - k), F)) : NaN
    R2 = 1 - SSR_u / sum((y .- mean(y)).^2)
    return (F, pval, R2, n)
end
Fg, pg, R2g, ng = ftest_z(P.z_bin)
Fw, pw, R2w, nw = ftest_z(P.z_w)
open(joinpath(DIAG_OUT, "G_unpredictability.txt"), "w") do io
    @printf(io, "z binary: F=%.4f p=%s R2=%.4f N=%d\n", Fg, string(pg), R2g, ng)
    @printf(io, "z_w:      F=%.4f p=%s R2=%.4f N=%d\n", Fw, string(pw), R2w, nw)
    println(io, "Controls: 6 lags Δlog q_aus, 6 lags Δlog p_real, 6 lags REA, month FE.")
    println(io, "Joint F on the 18 lag coefficients.")
end
CSV.write(joinpath(DIAG_OUT, "G_unpredictability.csv"),
          DataFrame(instrument = ["z", "z_w"], F = [Fg, Fw], p = [pg, pw],
                    R2 = [R2g, R2w], N = [ng, nw]))

# ============================ Task H ============================
println("H: randomisation inference ($NBOOT_H placebos)...")
true_idx = findall(>(0), z_pref)
n_events = length(true_idx)
season = Int[]
for i in 1:T
    m = month(P.dates[i])
    (m >= 11 || m <= 4) && push!(season, i)
end
true_set = Set(true_idx)
pool = [i for i in season if !(i in true_set)]
length(pool) < n_events && error("placebo pool too small: $(length(pool)) < $n_events")

βq_true, _, _ = lp_harm_h(P.lq, z_pref, P.rea, P.dates, 0)
βp_true, _, _ = lp_harm_h(P.lp, z_pref, P.rea, P.dates, 0)
true_vals = z_pref[true_idx]

rng = MersenneTwister(BOOT_SEED)
Bq = fill(NaN, NBOOT_H); Bp = fill(NaN, NBOOT_H)
for b in 1:NBOOT_H
    pick = pool[randperm(rng, length(pool))[1:n_events]]
    zpl = zeros(T)
    for (k, i) in enumerate(pick)
        zpl[i] = true_vals[k]
    end
    Bq[b], _, _ = lp_harm_h(P.lq, zpl, P.rea, P.dates, 0)
    Bp[b], _, _ = lp_harm_h(P.lp, zpl, P.rea, P.dates, 0)
    b % 2000 == 0 && @printf("  placebo %d/%d\n", b, NBOOT_H)
end
cq = filter(isfinite, Bq); cp = filter(isfinite, Bp)
pq = mean(abs.(cq) .>= abs(βq_true) - 0*βq_true)
pp = mean(abs.(cp) .>= abs(βp_true))
open(joinpath(DIAG_OUT, "H_randomisation.txt"), "w") do io
    @printf(io, "harmonised h=0 Aus exports: beta=%+.6f  rand_p=%.4f  (N_finite=%d)\n",
            βq_true, pq, length(cq))
    @printf(io, "harmonised h=0 price:       beta=%+.6f  rand_p=%.4f  (N_finite=%d)\n",
            βp_true, pp, length(cp))
    @printf(io, "n_events=%d  pool_size=%d  season=Nov-Apr  draws=%d  seed=%d\n",
            n_events, length(pool), NBOOT_H, BOOT_SEED)
end
CSV.write(joinpath(DIAG_OUT, "H_placebo_draws.csv"),
          DataFrame(draw = 1:NBOOT_H, beta_q = Bq, beta_p = Bp))

# ============================ Task I ============================
println("I: Japan decomposition...")
jp = CSV.read(joinpath(RAW, "bilateral_jpn.csv"), DataFrame)
jp.ds = first.(string.(jp.date), 7)
function jpull(part)
    sub = jp[jp.partner .== part, :]
    m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(sub))
    [haskey(m, string(s)) ? m[string(s)] : NaN for s in P.ds]
end
jw = jpull("World"); ja = jpull("Australia")
jn = [isfinite(jw[i]) && isfinite(ja[i]) ? jw[i] - ja[i] : NaN for i in 1:T]
jwl = [isfinite(x) && x > 0 ? log(x) : NaN for x in jw]
jal = [isfinite(x) && x > 0 ? log(x) : NaN for x in ja]
jnl = [isfinite(x) && x > 0 ? log(x) : NaN for x in jn]
jwt = jw ./ 1e6; jat = ja ./ 1e6; jnt = jn ./ 1e6

# Use 1000 MBB for I paths (10k would be very slow for 6 series × 13 h); note in file
const NBOOT_I = 1000
i1 = DataFrame(series = String[], h = Int[], beta = Float64[],
               mbb_lo = Float64[], mbb_hi = Float64[])
for (name, y) in (("non_aus_log", jnl), ("non_aus_mt", jnt),
                  ("total_log", jwl), ("from_aus_log", jal),
                  ("total_mt", jwt), ("from_aus_mt", jat))
    path = lp_harm_path(y, z_pref, P.rea, P.dates, 12)
    BS = fill(NaN, NBOOT_I, 13)
    rngi = MersenneTwister(BOOT_SEED)
    for b in 1:NBOOT_I
        ii = block_idx(T, BLOCK, rngi)
        BS[b, :] = lp_harm_path(y[ii], z_pref[ii], P.rea[ii], P.dates[ii], 12)
    end
    for h in 0:12
        c = filter(isfinite, BS[:, h + 1])
        push!(i1, (name, h, path[h + 1],
                   isempty(c) ? NaN : quantile(c, 0.05),
                   isempty(c) ? NaN : quantile(c, 0.95)))
    end
end
CSV.write(joinpath(DIAG_OUT, "I1_japan_non_aus.csv"), i1)

share = 0.663
β_aus0 = i1[(i1.series .== "from_aus_log") .& (i1.h .== 0), :].beta[1]
β_tot0 = i1[(i1.series .== "total_log") .& (i1.h .== 0), :].beta[1]
mech = share * β_aus0
gap = β_tot0 - mech
rngg = MersenneTwister(BOOT_SEED)
gaps = fill(NaN, NBOOT_I)
for b in 1:NBOOT_I
    ii = block_idx(T, BLOCK, rngg)
    ba = lp_harm_path(jal[ii], z_pref[ii], P.rea[ii], P.dates[ii], 0)[1]
    bt = lp_harm_path(jwl[ii], z_pref[ii], P.rea[ii], P.dates[ii], 0)[1]
    gaps[b] = bt - share * ba
end
cg = filter(isfinite, gaps)
open(joinpath(DIAG_OUT, "I2_accounting_gap.txt"), "w") do io
    @printf(io, "Aus share of Japan imports (paper): %.3f\n", share)
    @printf(io, "h=0 from-Aus log response: %+.6f\n", β_aus0)
    @printf(io, "mechanical total = share * fromAus: %+.6f\n", mech)
    @printf(io, "observed total log response: %+.6f\n", β_tot0)
    @printf(io, "gap (observed - mechanical): %+.6f\n", gap)
    @printf(io, "gap MBB 90%% [%+.6f, %+.6f]  (NBOOT_I=%d, recomputed per draw)\n",
            quantile(cg, 0.05), quantile(cg, 0.95), NBOOT_I)
end

# ============================ Task J ============================
println("J: specification robustness...")
mine = Set(["2020-05", "2024-06"])
function z_custom(; weighted=true, drop2022=true, drop_mines=false, major_only=false)
    z = major_only ? copy(P.z_maj) : (weighted ? copy(P.z_w) : copy(P.z_bin))
    drop2022 && (z = apply_drop!(z, P.ds, DROP_2022))
    if drop_mines
        for i in 1:T
            string(P.ds[i]) in mine && (z[i] = 0.0)
        end
    end
    return z
end
jrows = DataFrame(spec = String[], h = Int[],
                  beta_q_aus = Float64[], beta_price = Float64[])
specs = [
    ("baseline_zw_drop2022_real", z_custom(), P.lp),
    ("z_binary_drop2022_real", z_custom(weighted=false), P.lp),
    ("z_major_drop2022_real", z_custom(major_only=true), P.lp),
    ("zw_keep2022_real", z_custom(drop2022=false), P.lp),
    ("zw_drop2022_drop_mines_real", z_custom(drop_mines=true), P.lp),
    ("zw_drop2022_nominal_price", z_custom(), P.p_nom),
]
for (name, z, pseries) in specs
    for h in (0, 3)
        bq, _, _ = lp_harm_h(P.lq, z, P.rea, P.dates, h)
        bp, _, _ = lp_harm_h(pseries, z, P.rea, P.dates, h)
        push!(jrows, (name, h, bq, bp))
    end
end
CSV.write(joinpath(DIAG_OUT, "J_robustness.csv"), jrows)

group_file = joinpath(REPO, "out", "group_first_stage.txt")
open(joinpath(DIAG_OUT, "J_high_low_note.txt"), "w") do io
    if isfile(group_file)
        println(io, "08_estimate_by_group.jl / group_first_stage.txt:")
        print(io, read(group_file, String))
        println(io, "Full HIGH/LOW bilateral re-estimation under harm LP: not re-run here.")
    else
        println(io, "NOT AVAILABLE: group_first_stage.txt missing")
    end
end
svar = joinpath(REPO, "out", "irf_svar.csv")
open(joinpath(DIAG_OUT, "J_proxy_svar.txt"), "w") do io
    if isfile(svar)
        println(io, "EXISTS: out/irf_svar.csv from 02_estimate.jl")
        println(io, read(svar, String))
    else
        println(io, "NOT AVAILABLE: proxy SVAR output missing")
    end
end

# ============================ Task L ============================
println("L: text/figure consistency...")
fs = read(joinpath(REPO, "out", "first_stage.txt"), String)
lk = CSV.read(joinpath(REPO, "out", "leakage.csv"), DataFrame)
open(joinpath(DIAG_OUT, "L_consistency.txt"), "w") do io
    println(io, "L1: first_stage.txt:")
    println(io, fs)
    m = match(r"price\s+beta=([+\-0-9.eE]+)", fs)
    println(io, "True pipeline price beta: ", m === nothing ? "PARSE FAIL" : m.captures[1])
    println(io, "Section 3.1 '+0.00' rounds the near-zero estimate; Table 2 -0.003 is the signed value to cite.")
    println(io, "")
    println(io, "L2: YES — 06_figures.jl draws shaded bands only for Australia and Indonesia.")
    println(io, "")
    println(io, "L3: fig3 and fig4 plot Japan, Korea, and India. Captions that say only 'Japan:' are incomplete.")
    println(io, "")
    println(io, "L4: fig2 uses ylims!(ax, -2, 3) and clamps bands to [-2, 3].")
    finite = filter(isfinite, lk.leakage)
    @printf(io, "leakage.csv finite min/max: %.6f / %.6f\n", minimum(finite), maximum(finite))
    println(io, "Series is CLIPPED in the figure when |leakage| > 3.")
end

println("wrote G H I J L outputs")
println("script: diagnostics/task_GHIJL.jl")
