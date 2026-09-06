# diagnostics/task_UVX.jl
# Tasks U (W^exc leakage), V (randomisation), X (baseline table).
# New baseline: z_w, 2022 RETAINED, 6 own lags, REA+6, month FE.
include(joinpath(@__DIR__, "exporters_harm.jl"))

const NBOOT = 10_000
const H_X = 12
const H_U = 6

P = load_panel()
z_keep = z_baseline(P)
z_drop = z_preferred(P)
T = length(P.dates)
E = load_exporters_harm(P)

qband(v) = (c = filter(isfinite, v); isempty(c) ? (NaN, NaN) : (quantile(c, 0.05), quantile(c, 0.95)))
leak_ratio(wexc, aus) = (isfinite(aus) && aus < -0.01 && isfinite(wexc)) ? -(wexc) / aus : NaN
cum_leak(cw, ca) = (isfinite(ca) && abs(ca) > 1e-8 && isfinite(cw)) ? -cw / ca : NaN

# Japan series
jp = CSV.read(joinpath(RAW, "bilateral_jpn.csv"), DataFrame)
jp.ds = first.(string.(jp.date), 7)
function jpull(part)
    sub = jp[jp.partner .== part, :]
    m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(sub))
    [haskey(m, string(s)) ? m[string(s)] : NaN for s in P.ds]
end
jw = jpull("World"); ja = jpull("Australia")
jn = [isfinite(jw[i]) && isfinite(ja[i]) ? jw[i] - ja[i] : NaN for i in 1:T]
jp_series = Dict(
    "Japan_total_log" => [isfinite(x) && x > 0 ? log(x) : NaN for x in jw],
    "Japan_from_aus_log" => [isfinite(x) && x > 0 ? log(x) : NaN for x in ja],
    "Japan_non_aus_log" => [isfinite(x) && x > 0 ? log(x) : NaN for x in jn],
    "Japan_total_mt" => jw ./ 1e6,
    "Japan_from_aus_mt" => ja ./ 1e6,
    "Japan_non_aus_mt" => jn ./ 1e6,
)

# ---- shared point paths keep-2022 ----
println("UVX: point paths (keep-2022)...")
B_aus = lp_harm_path(E.Q_mt[:, E.aidx], z_keep, P.rea, P.dates, H_X)
B_wexc = lp_harm_path(E.W_exc, z_keep, P.rea, P.dates, H_X)
B_lq = lp_harm_path(P.lq, z_keep, P.rea, P.dates, H_X)
B_lp = lp_harm_path(P.lp, z_keep, P.rea, P.dates, H_X)
B_comp = fill(NaN, H_X + 1, E.J)
for j in 1:E.J
    B_comp[:, j] = lp_harm_path(E.Q_mt[:, j], z_keep, P.rea, P.dates, H_X)
end
B_jp = Dict(k => lp_harm_path(v, z_keep, P.rea, P.dates, H_X) for (k, v) in jp_series)

# drop-2022 points for U3
B_aus_d = lp_harm_path(E.Q_mt[:, E.aidx], z_drop, P.rea, P.dates, 3)
B_wexc_d = lp_harm_path(E.W_exc, z_drop, P.rea, P.dates, 3)

println("UVX: MBB keep-2022 N=$NBOOT ...")
rng = MersenneTwister(BOOT_SEED)
BS_aus = fill(NaN, NBOOT, H_X + 1)
BS_wexc = fill(NaN, NBOOT, H_X + 1)
BS_lq = fill(NaN, NBOOT, H_X + 1)
BS_lp = fill(NaN, NBOOT, H_X + 1)
BS_comp = fill(NaN, NBOOT, H_X + 1, E.J)
BS_jp = Dict(k => fill(NaN, NBOOT, H_X + 1) for k in keys(jp_series))
for b in 1:NBOOT
    ii = block_idx(T, BLOCK, rng)
    zb = z_keep[ii]; rb = P.rea[ii]; db = P.dates[ii]
    BS_aus[b, :] = lp_harm_path(E.Q_mt[ii, E.aidx], zb, rb, db, H_X)
    BS_wexc[b, :] = lp_harm_path(E.W_exc[ii], zb, rb, db, H_X)
    BS_lq[b, :] = lp_harm_path(P.lq[ii], zb, rb, db, H_X)
    BS_lp[b, :] = lp_harm_path(P.lp[ii], zb, rb, db, H_X)
    for j in 1:E.J
        BS_comp[b, :, j] = lp_harm_path(E.Q_mt[ii, j], zb, rb, db, H_X)
    end
    for (k, v) in jp_series
        BS_jp[k][b, :] = lp_harm_path(v[ii], zb, rb, db, H_X)
    end
    b % 2000 == 0 && @printf("  keep draw %d/%d\n", b, NBOOT)
end

println("UVX: MBB drop-2022 (U3) N=$NBOOT ...")
rngd = MersenneTwister(BOOT_SEED)
BS_aus_d = fill(NaN, NBOOT, 4)
BS_wexc_d = fill(NaN, NBOOT, 4)
for b in 1:NBOOT
    ii = block_idx(T, BLOCK, rngd)
    zb = z_drop[ii]; rb = P.rea[ii]; db = P.dates[ii]
    BS_aus_d[b, :] = lp_harm_path(E.Q_mt[ii, E.aidx], zb, rb, db, 3)
    BS_wexc_d[b, :] = lp_harm_path(E.W_exc[ii], zb, rb, db, 3)
    b % 2000 == 0 && @printf("  drop draw %d/%d\n", b, NBOOT)
end

# ===================== U =====================
println("UVX: writing U...")
u1 = DataFrame(h = Int[], wexc = Float64[], wexc_lo = Float64[], wexc_hi = Float64[],
               aus = Float64[], aus_lo = Float64[], aus_hi = Float64[],
               leakage = Float64[], leak_lo = Float64[], leak_hi = Float64[])
for h in 0:H_U
    lo_w, hi_w = qband(BS_wexc[:, h + 1])
    lo_a, hi_a = qband(BS_aus[:, h + 1])
    lk = leak_ratio(B_wexc[h + 1], B_aus[h + 1])
    boots = [leak_ratio(BS_wexc[b, h + 1], BS_aus[b, h + 1]) for b in 1:NBOOT]
    llo, lhi = qband(boots)
    push!(u1, (h, B_wexc[h + 1], lo_w, hi_w, B_aus[h + 1], lo_a, hi_a, lk, llo, lhi))
end
CSV.write(joinpath(DIAG_OUT, "U1_per_horizon_leakage_keep2022.csv"), u1)

function write_cum(path, Baus, Bwexc, BSaus, BSwexc; Hs = 0:3)
    out = DataFrame(H = Int[], CumAUS = Float64[], CumWexc = Float64[],
                    cum_leak = Float64[], lo90 = Float64[], hi90 = Float64[])
    for H in Hs
        cA = sum(Baus[h + 1] for h in 0:H)
        cW = sum(Bwexc[h + 1] for h in 0:H)
        lk = cum_leak(cW, cA)
        boots = fill(NaN, size(BSaus, 1))
        for b in 1:size(BSaus, 1)
            cAb = sum(BSaus[b, h + 1] for h in 0:H)
            cWb = sum(BSwexc[b, h + 1] for h in 0:H)
            boots[b] = cum_leak(cWb, cAb)
        end
        lo, hi = qband(boots)
        push!(out, (H, cA, cW, lk, lo, hi))
    end
    CSV.write(path, out)
    return out
end

u2 = write_cum(joinpath(DIAG_OUT, "U2_cum_leakage_keep2022.csv"),
               B_aus, B_wexc, BS_aus, BS_wexc)
u3 = write_cum(joinpath(DIAG_OUT, "U3_cum_leakage_drop2022.csv"),
               B_aus_d, B_wexc_d, BS_aus_d, BS_wexc_d)

open(joinpath(DIAG_OUT, "U4_H3_band_check.txt"), "w") do io
    r = u2[u2.H .== 3, :]
    lo, hi, pt = r.lo90[1], r.hi90[1], r.cum_leak[1]
    excl0 = isfinite(lo) && isfinite(hi) && (lo > 0 || hi < 0)
    excl073 = isfinite(lo) && isfinite(hi) && (0.73 < lo || 0.73 > hi)
    @printf(io, "KEEP-2022 cumulative leakage H=3:\n")
    @printf(io, "  point=%.4f  90%% band=[%.4f, %.4f]\n", pt, lo, hi)
    @printf(io, "  excludes zero? %s\n", excl0)
    @printf(io, "  excludes 0.73? %s\n", excl073)
    rd = u3[u3.H .== 3, :]
    @printf(io, "DROP-2022 cumulative leakage H=3 (robustness):\n")
    @printf(io, "  point=%.4f  90%% band=[%.4f, %.4f]\n",
            rd.cum_leak[1], rd.lo90[1], rd.hi90[1])
end

# ===================== X =====================
println("UVX: writing X baseline table...")
horizons_x = [0, 1, 2, 3, 6, 12]
xrows = DataFrame(object = String[], h = Int[], beta = Float64[],
                  lo90 = Float64[], hi90 = Float64[])

function push_x!(name, pathβ, BS)
    for h in horizons_x
        lo, hi = qband(BS[:, h + 1])
        push!(xrows, (name, h, pathβ[h + 1], lo, hi))
    end
end

push_x!("Aus_exports_log", B_lq, BS_lq)
push_x!("Aus_exports_mt", B_aus, BS_aus)
push_x!("Price_log_real", B_lp, BS_lp)
push_x!("W_exc_mt", B_wexc, BS_wexc)
for j in 1:E.J
    push_x!("$(E.labels[j])_mt", B_comp[:, j], BS_comp[:, :, j])
end
for (k, pathβ) in B_jp
    push_x!(k, pathβ, BS_jp[k])
end
CSV.write(joinpath(DIAG_OUT, "X_baseline_keep2022.csv"), xrows)

open(joinpath(DIAG_OUT, "X_baseline_keep2022.md"), "w") do io
    println(io, "# Baseline table (keep-2022 harmonised LP)")
    println(io, "")
    println(io, "Spec: ``z_w``, 2022 RETAINED, 6 own lags, REA+6 lags, month FE.")
    println(io, "MBB 90% bands, N=$(NBOOT), seed=$(BOOT_SEED), block=$(BLOCK).")
    println(io, "")
    println(io, "| object | h=0 | h=1 | h=2 | h=3 | h=6 | h=12 |")
    println(io, "|--------|----:|----:|----:|----:|----:|-----:|")
    objs = unique(xrows.object)
    for obj in objs
        cells = String[]
        for h in horizons_x
            r = xrows[(xrows.object .== obj) .& (xrows.h .== h), :]
            push!(cells, @sprintf("%.3f [%.2f, %.2f]", r.beta[1], r.lo90[1], r.hi90[1]))
        end
        println(io, "| ", obj, " | ", join(cells, " | "), " |")
    end
end

# ===================== V =====================
println("UVX: randomisation V (10000 placebos)...")

function placebo_pool(z_ref)
    true_idx = findall(>(0), z_ref)
    season = Int[i for i in 1:T if (m = month(P.dates[i]); m >= 11 || m <= 4)]
    pool = [i for i in season if !(i in Set(true_idx))]
    return true_idx, z_ref[true_idx], pool
end

function rand_p_h(y, z_ref, h; nboot = NBOOT, seed = BOOT_SEED)
    true_idx, true_vals, pool = placebo_pool(z_ref)
    n_events = length(true_idx)
    length(pool) < n_events && return (NaN, NaN, 0, n_events, length(pool))
    β0, _, _ = lp_harm_h(y, z_ref, P.rea, P.dates, h)
    !isfinite(β0) && return (β0, NaN, 0, n_events, length(pool))
    rng = MersenneTwister(seed)
    B = fill(NaN, nboot)
    for b in 1:nboot
        pick = pool[randperm(rng, length(pool))[1:n_events]]
        zpl = zeros(T)
        for (k, i) in enumerate(pick)
            zpl[i] = true_vals[k]
        end
        B[b], _, _ = lp_harm_h(y, zpl, P.rea, P.dates, h)
    end
    c = filter(isfinite, B)
    p = isempty(c) ? NaN : mean(abs.(c) .>= abs(β0))
    return (β0, p, length(c), n_events, length(pool))
end

function rand_p_cum_wexc(H; nboot = NBOOT, seed = BOOT_SEED)
    true_idx, true_vals, pool = placebo_pool(z_keep)
    n_events = length(true_idx)
    length(pool) < n_events && return (NaN, NaN, 0, n_events, length(pool))
    path0 = lp_harm_path(E.W_exc, z_keep, P.rea, P.dates, H)
    β0 = sum(path0[h + 1] for h in 0:H)
    !isfinite(β0) && return (β0, NaN, 0, n_events, length(pool))
    rng = MersenneTwister(seed)
    B = fill(NaN, nboot)
    for b in 1:nboot
        pick = pool[randperm(rng, length(pool))[1:n_events]]
        zpl = zeros(T)
        for (k, i) in enumerate(pick)
            zpl[i] = true_vals[k]
        end
        path = lp_harm_path(E.W_exc, zpl, P.rea, P.dates, H)
        B[b] = sum(path[h + 1] for h in 0:H)
    end
    c = filter(isfinite, B)
    p = isempty(c) ? NaN : mean(abs.(c) .>= abs(β0))
    return (β0, p, length(c), n_events, length(pool))
end

vrows = DataFrame(test = String[], h = Int[], beta = Float64[],
                  rand_p = Float64[], n_finite = Int[], n_events = Int[], pool = Int[])

v_tests = [
    ("V1_Wexc_mt", E.W_exc, 1),
    ("V1_Wexc_mt", E.W_exc, 2),
    ("V1_Wexc_mt", E.W_exc, 3),
    ("V3_Aus_mt", E.Q_mt[:, E.aidx], 1),
    ("V3_Aus_mt", E.Q_mt[:, E.aidx], 2),
    ("V4_Japan_from_aus_log", jp_series["Japan_from_aus_log"], 0),
    ("V5_price_keep2022", P.lp, 0),
]

for (name, y, h) in v_tests
    @printf("  %s h=%d ...\n", name, h)
    β, p, nf, ne, pool = rand_p_h(y, z_keep, h)
    push!(vrows, (name, h, β, p, nf, ne, pool))
end
@printf("  V2_CumWexc_H3 ...\n")
β, p, nf, ne, pool = rand_p_cum_wexc(3)
push!(vrows, ("V2_CumWexc", 3, β, p, nf, ne, pool))

CSV.write(joinpath(DIAG_OUT, "V_randomisation.csv"), vrows)

println("wrote U/V/X outputs")
println("script: diagnostics/task_UVX.jl")
