# diagnostics/task_PQ.jl — Fix world leakage (P) + horizon paths (Q)
# One 10000-draw MBB shared across P and Q for efficiency.
include(joinpath(@__DIR__, "exporters_harm.jl"))

const NBOOT = 10_000
const Hmax = 18

P = load_panel()
z_pref = z_preferred(P)
T = length(P.dates)
E = load_exporters_harm(P)

open(joinpath(DIAG_OUT, "P_definition.txt"), "w") do io
    println(io, "=== BUG IN M2 world_aggregate cell ===")
    println(io, "M2 built Qworld as SUM of ALL included exporters INCLUDING Australia,")
    println(io, "then coded implied_leak = 1 - ΔW / (-ΔAUS).")
    println(io, "Correct identity when W includes Aus: leakage = (ΔW - ΔAUS)/(-ΔAUS) = 1 + ΔW/(-ΔAUS).")
    println(io, "M2 used 1 - ΔW/(-ΔAUS), which flips the sign of the ΔW term → spurious ~+2.")
    println(io, "")
    println(io, "=== CORRECT DEFINITIONS (this task) ===")
    println(io, "Coverage set (03_leakage gates): ", join(E.labels, ", "))
    println(io, "Australia index: ", E.aidx, " (", E.labels[E.aidx], ")")
    println(io, "W_inc[t] = sum_j Q_mt[t,j] over finite j, requiring Aus observed")
    println(io, "W_exc[t] = sum_{j≠Aus} Q_mt[t,j] over finite j, requiring Aus observed")
    println(io, "Australia response: LP of Aus Mt levels (harmonised)")
    println(io, "Leakage_h = - (β_W_exc,h) / (β_Aus,h)   when β_Aus < 0")
    println(io, "Ratio recomputed inside each MBB draw (seed=$BOOT_SEED, N=$NBOOT, block=$BLOCK).")
    println(io, "Estimator: direct Mt levels under harmonised LP.")
end

println("PQ: point paths...")
B_aus = lp_harm_path(E.Q_mt[:, E.aidx], z_pref, P.rea, P.dates, Hmax)
B_winc = lp_harm_path(E.W_inc, z_pref, P.rea, P.dates, Hmax)
B_wexc = lp_harm_path(E.W_exc, z_pref, P.rea, P.dates, Hmax)
B_comp = fill(NaN, Hmax + 1, E.J)
for j in 1:E.J
    B_comp[:, j] = lp_harm_path(E.Q_mt[:, j], z_pref, P.rea, P.dates, Hmax)
end
B_sum_ex = [sum(B_comp[h + 1, j] for j in 1:E.J if j != E.aidx && isfinite(B_comp[h + 1, j]); init = 0.0)
            for h in 0:Hmax]

println("PQ: MBB $NBOOT ...")
rng = MersenneTwister(BOOT_SEED)
BS_aus = fill(NaN, NBOOT, Hmax + 1)
BS_winc = fill(NaN, NBOOT, Hmax + 1)
BS_wexc = fill(NaN, NBOOT, Hmax + 1)
BS_comp = fill(NaN, NBOOT, Hmax + 1, E.J)
for b in 1:NBOOT
    ii = block_idx(T, BLOCK, rng)
    zb = z_pref[ii]; rb = P.rea[ii]; db = P.dates[ii]
    BS_aus[b, :] = lp_harm_path(E.Q_mt[ii, E.aidx], zb, rb, db, Hmax)
    BS_winc[b, :] = lp_harm_path(E.W_inc[ii], zb, rb, db, Hmax)
    BS_wexc[b, :] = lp_harm_path(E.W_exc[ii], zb, rb, db, Hmax)
    for j in 1:E.J
        BS_comp[b, :, j] = lp_harm_path(E.Q_mt[ii, j], zb, rb, db, Hmax)
    end
    b % 2000 == 0 && @printf("  draw %d/%d\n", b, NBOOT)
end

qband(v) = (c = filter(isfinite, v); isempty(c) ? (NaN, NaN) : (quantile(c, 0.05), quantile(c, 0.95)))

function leak_pt(wexc, aus)
    (isfinite(aus) && aus < -0.01 && isfinite(wexc)) ? -(wexc) / aus : NaN
end

# ---- P table ----
pout = DataFrame(h = Int[],
                 world_inc_Aus = Float64[], winc_lo = Float64[], winc_hi = Float64[],
                 world_ex_Aus = Float64[], wexc_lo = Float64[], wexc_hi = Float64[],
                 Australia = Float64[], aus_lo = Float64[], aus_hi = Float64[],
                 leakage = Float64[], leak_lo = Float64[], leak_hi = Float64[])
for h in 0:Hmax
    lo_a, hi_a = qband(BS_aus[:, h + 1])
    lo_i, hi_i = qband(BS_winc[:, h + 1])
    lo_e, hi_e = qband(BS_wexc[:, h + 1])
    lk = leak_pt(B_wexc[h + 1], B_aus[h + 1])
    boots = [leak_pt(BS_wexc[b, h + 1], BS_aus[b, h + 1]) for b in 1:NBOOT]
    llo, lhi = qband(boots)
    push!(pout, (h,
                 B_winc[h + 1], lo_i, hi_i,
                 B_wexc[h + 1], lo_e, hi_e,
                 B_aus[h + 1], lo_a, hi_a,
                 lk, llo, lhi))
end
CSV.write(joinpath(DIAG_OUT, "P_world_leakage_fixed.csv"), pout)

# ---- Q1 Australia ----
q1 = DataFrame(h = Int[], beta_mt = Float64[], lo90 = Float64[], hi90 = Float64[])
for h in 0:Hmax
    lo, hi = qband(BS_aus[:, h + 1])
    push!(q1, (h, B_aus[h + 1], lo, hi))
end
CSV.write(joinpath(DIAG_OUT, "Q1_australia_mt.csv"), q1)
open(joinpath(DIAG_OUT, "Q1_catchup.txt"), "w") do io
    pos_long = [(h, B_aus[h + 1]) for h in 6:Hmax if isfinite(B_aus[h + 1]) && B_aus[h + 1] > 0]
    println(io, "Australia Mt path (levels, harmonised). Positive long-horizon points:")
    if isempty(pos_long)
        println(io, "NONE — no catch-up (no positive β at h=6..18).")
    else
        for (h, b) in pos_long
            @printf(io, "  h=%d β=%+.4f\n", h, b)
        end
    end
    @printf(io, "h=0: %+.4f  h=12: %+.4f  h=18: %+.4f\n",
            B_aus[1], B_aus[13], B_aus[19])
end

# ---- Q2 competitors ----
q2 = DataFrame(exporter = String[], h = Int[], beta_mt = Float64[],
               lo90 = Float64[], hi90 = Float64[])
for j in 1:E.J, h in 0:Hmax
    lo, hi = qband(BS_comp[:, h + 1, j])
    push!(q2, (E.labels[j], h, B_comp[h + 1, j], lo, hi))
end
# sum of competitors
for h in 0:Hmax
    boots = fill(NaN, NBOOT)
    for b in 1:NBOOT
        boots[b] = sum(BS_comp[b, h + 1, j] for j in 1:E.J if j != E.aidx && isfinite(BS_comp[b, h + 1, j]); init = 0.0)
    end
    lo, hi = qband(boots)
    push!(q2, ("SUM_ex_Australia", h, B_sum_ex[h + 1], lo, hi))
end
CSV.write(joinpath(DIAG_OUT, "Q2_competitors_mt.csv"), q2)

# ---- Q3 Japan ----
println("PQ: Japan paths...")
jp = CSV.read(joinpath(RAW, "bilateral_jpn.csv"), DataFrame)
jp.ds = first.(string.(jp.date), 7)
function jpull(part)
    sub = jp[jp.partner .== part, :]
    m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(sub))
    [haskey(m, string(s)) ? m[string(s)] : NaN for s in P.ds]
end
jw = jpull("World"); ja = jpull("Australia")
jn = [isfinite(jw[i]) && isfinite(ja[i]) ? jw[i] - ja[i] : NaN for i in 1:T]
series_j = Dict(
    "total_log" => [isfinite(x) && x > 0 ? log(x) : NaN for x in jw],
    "from_aus_log" => [isfinite(x) && x > 0 ? log(x) : NaN for x in ja],
    "non_aus_log" => [isfinite(x) && x > 0 ? log(x) : NaN for x in jn],
    "total_mt" => jw ./ 1e6,
    "from_aus_mt" => ja ./ 1e6,
    "non_aus_mt" => jn ./ 1e6,
)
q3 = DataFrame(series = String[], h = Int[], beta = Float64[], lo90 = Float64[], hi90 = Float64[])
for (name, y) in series_j
    path = lp_harm_path(y, z_pref, P.rea, P.dates, 12)
    BS = fill(NaN, NBOOT, 13)
    rngj = MersenneTwister(BOOT_SEED)
    for b in 1:NBOOT
        ii = block_idx(T, BLOCK, rngj)
        BS[b, :] = lp_harm_path(y[ii], z_pref[ii], P.rea[ii], P.dates[ii], 12)
    end
    for h in 0:12
        lo, hi = qband(BS[:, h + 1])
        push!(q3, (name, h, path[h + 1], lo, hi))
    end
end
CSV.write(joinpath(DIAG_OUT, "Q3_japan_paths.csv"), q3)
open(joinpath(DIAG_OUT, "Q3_non_aus_recovery.txt"), "w") do io
    println(io, "Japan non-Australian imports (harmonised), logs and Mt:")
    for unit in ("non_aus_log", "non_aus_mt")
        println(io, unit, ":")
        for h in 0:12
            r = q3[(q3.series .== unit) .& (q3.h .== h), :]
            @printf(io, "  h=%2d  β=%+.4f  [%+.4f, %+.4f]\n",
                    h, r.beta[1], r.lo90[1], r.hi90[1])
        end
        # recovery = later horizons less negative / positive vs h=0
        b0 = q3[(q3.series .== unit) .& (q3.h .== 0), :].beta[1]
        later = [q3[(q3.series .== unit) .& (q3.h .== h), :].beta[1] for h in 1:12]
        recovered = any(isfinite(b0) && isfinite(b) && b > b0 + 0.01 for b in later)  # improved vs impact
        pos_later = any(isfinite(b) && b > 0 for b in later)
        println(io, "  vs h=0: any later β > h=0+0.01? ", recovered)
        println(io, "  any later β > 0? ", pos_later)
    end
end

# ---- Q4 cumulative leakage every H ----
q4 = DataFrame(H = Int[], CumAUS = Float64[], CumEX = Float64[],
               cum_leak = Float64[], lo90 = Float64[], hi90 = Float64[])
for H in 0:Hmax
    cA = sum(B_aus[h + 1] for h in 0:H)
    cE = sum(B_wexc[h + 1] for h in 0:H)  # use world-ex path cumulative
    # Prefer summing competitor level responses for consistency with W_exc LP
    cE2 = sum(B_sum_ex[h + 1] for h in 0:H)
    lk = (isfinite(cA) && abs(cA) > 1e-8) ? -cE2 / cA : NaN
    boots = fill(NaN, NBOOT)
    for b in 1:NBOOT
        cAb = sum(BS_aus[b, h + 1] for h in 0:H)
        cEb = 0.0
        for h in 0:H, j in 1:E.J
            j == E.aidx && continue
            isfinite(BS_comp[b, h + 1, j]) && (cEb += BS_comp[b, h + 1, j])
        end
        boots[b] = (isfinite(cAb) && abs(cAb) > 1e-8) ? -cEb / cAb : NaN
    end
    lo, hi = qband(boots)
    push!(q4, (H, cA, cE2, lk, lo, hi))
end
CSV.write(joinpath(DIAG_OUT, "Q4_cumulative_leakage.csv"), q4)

println("wrote P and Q outputs")
println("script: diagnostics/task_PQ.jl")
