# diagnostics/task_M.jl — Harmonise 03/05 to paper specification
include(joinpath(@__DIR__, "harm_lp.jl"))

const EXPORTERS = [
    ("aus_exports.csv", "Australia", true),
    ("exp_idn.csv", "Indonesia", false),
    ("exp_usa.csv", "USA", false),
    ("exp_col.csv", "Colombia", false),
    ("exp_can.csv", "Canada", false),
    ("exp_rus.csv", "Russia", false),
    ("exp_zaf.csv", "South Africa", false),
    ("exp_mng.csv", "Mongolia", false),
]

function readexp(fname)
    path = joinpath(RAW, fname)
    !isfile(path) && return nothing
    d = CSV.read(path, DataFrame)
    d.ds = first.(string.(d.date), 7)
    d = d[d.ds .>= SAMPLE_START, :]
    isempty(d) && return nothing
    return d
end
cal_span(a, b) = [Dates.format(d, "yyyy-mm") for d in Date(a * "-01"):Month(1):Date(b * "-01")]

P = load_panel()
z_pref = z_preferred(P)
T = length(P.dates)

# --- select exporters with 03 coverage gates, align to PANEL calendar ---
series = Tuple{String,Bool,DataFrame}[]
for (f, lab, isa_) in EXPORTERS
    d = readexp(f)
    d === nothing && continue
    dens = nrow(d) / length(cal_span(d.ds[1], d.ds[end]))
    keep = isa_ || (nrow(d) >= 110 && d.ds[end] >= "2024-01" && dens >= 0.85)
    keep && push!(series, (lab, isa_, d))
end
labels = [s[1] for s in series]
isaus = [s[2] for s in series]
aidx = findfirst(isaus)

Q_t = fill(NaN, T, length(series))  # tonnes on panel calendar
for (j, (lab, _, d)) in enumerate(series)
    m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(d))
    for i in 1:T
        haskey(m, string(P.ds[i])) && (Q_t[i, j] = m[string(P.ds[i])])
    end
end
qbar = [mean(filter(isfinite, Q_t[:, j])) for j in 1:length(series)]
Q_mt = Q_t ./ 1e6
ylog = [log.(max.(Q_t[:, j], NaN)) for j in 1:length(series)]  # wrong - need elementwise
Ylog = map(j -> [isfinite(Q_t[i, j]) && Q_t[i, j] > 0 ? log(Q_t[i, j]) : NaN for i in 1:T],
           1:length(series))

# ---------- CURRENT (from pipeline out/) ----------
cur_exp = CSV.read(joinpath(REPO, "out", "exporter_responses.csv"), DataFrame)
cur_lk  = CSV.read(joinpath(REPO, "out", "leakage.csv"), DataFrame)
cur_imp = CSV.read(joinpath(REPO, "out", "importer_responses.csv"), DataFrame)

# ---------- HARMONISED point estimates ----------
println("M: harmonised point paths...")
Hmax = 18
B_log = fill(NaN, Hmax + 1, length(series))
B_lvl = fill(NaN, Hmax + 1, length(series))
for j in 1:length(series)
    B_log[:, j] = lp_harm_path(Ylog[j], z_pref, P.rea, P.dates, Hmax)
    B_lvl[:, j] = lp_harm_path(Q_mt[:, j], z_pref, P.rea, P.dates, Hmax)
end

# M4 check vs Table 2 / first_stage
β_harm_aus0 = B_log[1, aidx]
β_fs, _, _, _ = lp_fs(P.lq, z_pref, P.dates, P.lq, P.lp, P.rea, 0)
open(joinpath(DIAG_OUT, "M4_first_stage_check.txt"), "w") do io
    @printf(io, "harmonised Aus log export h=0 (own lags + REA, no price lags): %.8f\n", β_harm_aus0)
    @printf(io, "Table 2 / 02_estimate first_stage_F (q lags + p lags + REA): %.8f\n", β_fs)
    @printf(io, "target Table 2: -0.1196\n")
    @printf(io, "match_FS_style: %s\n", abs(β_fs - (-0.1196)) < 5e-4 ? "YES" : "NO")
    @printf(io, "match_harmonised_ownlags: %s\n", abs(β_harm_aus0 - (-0.1196)) < 5e-4 ? "YES" : "NO")
    if abs(β_harm_aus0 - (-0.1196)) >= 5e-4
        println(io, "WHY: Task M harmonised LP uses 6 OWN lags of the outcome + REA+6 lags + month FE.")
        println(io, "Table 2 first stage (02_estimate.jl first_stage_F) additionally includes 6 lags of")
        println(io, "log price (and uses log Aus exports as q). Different control set => different beta.")
        @printf(io, "FS-style beta equals Table 2: %.8f (should match).\n", β_fs)
    end
end

# ---------- MBB 10000 ----------
println("M: MBB $NBOOT_M draws (this takes a while)...")
rng = MersenneTwister(BOOT_SEED)
BS_log = fill(NaN, NBOOT_M, Hmax + 1, length(series))
BS_lvl = fill(NaN, NBOOT_M, Hmax + 1, length(series))
for b in 1:NBOOT_M
    ii = block_idx(T, BLOCK, rng)
    zb = z_pref[ii]; rb = P.rea[ii]; db = P.dates[ii]
    for j in 1:length(series)
        BS_log[b, :, j] = lp_harm_path(Ylog[j][ii], zb, rb, db, Hmax)
        BS_lvl[b, :, j] = lp_harm_path(Q_mt[ii, j], zb, rb, db, Hmax)
    end
    b % 1000 == 0 && @printf("  draw %d/%d\n", b, NBOOT_M)
end

qband(v) = (c = filter(isfinite, v); isempty(c) ? (NaN, NaN) : (quantile(c, 0.05), quantile(c, 0.95)))

# M1 table
m1 = DataFrame(exporter = String[], h = Int[],
               current_mt = Float64[], current_lo = Float64[], current_hi = Float64[],
               harm_log_mt = Float64[], harm_log_lo = Float64[], harm_log_hi = Float64[],
               harm_lvl_mt = Float64[], harm_lvl_lo = Float64[], harm_lvl_hi = Float64[])
for j in 1:length(series), h in 0:Hmax
    cur = cur_exp[(cur_exp.exporter .== labels[j]) .& (cur_exp.horizon .== h), :]
    cmt = isempty(cur) ? NaN : cur.mt[1]
    clo = isempty(cur) ? NaN : cur.lo90_mt[1]
    chi = isempty(cur) ? NaN : cur.hi90_mt[1]
    hmt = B_log[h + 1, j] * qbar[j] / 1e6
    hlo, hhi = qband(BS_log[:, h + 1, j] .* (qbar[j] / 1e6))
    lmt = B_lvl[h + 1, j]
    llo, lhi = qband(BS_lvl[:, h + 1, j])
    push!(m1, (labels[j], h, cmt, clo, chi, hmt, hlo, hhi, lmt, llo, lhi))
end
CSV.write(joinpath(DIAG_OUT, "M1_exporters_current_vs_harm.csv"), m1)

# M2 leakage helpers
function leak_from_B(Bmt_path, aidx, h)  # Bmt_path is vector over exporters at horizon h
    den = -Bmt_path[aidx]
    num = sum(Bmt_path[j] for j in 1:length(Bmt_path) if j != aidx && isfinite(Bmt_path[j]); init = 0.0)
    (isfinite(den) && den > 0.01) ? num / den : NaN
end
function cum_leak(Bmt_full, aidx, H)  # Bmt_full[h+1, j]
    cA = sum(Bmt_full[h + 1, aidx] for h in 0:H)
    cC = 0.0
    for j in 1:size(Bmt_full, 2), h in 0:H
        j == aidx && continue
        isfinite(Bmt_full[h + 1, j]) && (cC += Bmt_full[h + 1, j])
    end
    (isfinite(cA) && abs(cA) > 1e-8) ? -cC / cA : NaN
end

# point Mt paths
Mt_log = [B_log[h + 1, j] * qbar[j] / 1e6 for h in 0:Hmax, j in 1:length(series)]
Mt_lvl = copy(B_lvl)

m2 = DataFrame(kind = String[], H = Int[],
               current_leak = Float64[], current_lo = Float64[], current_hi = Float64[],
               harm_leak = Float64[], harm_lo = Float64[], harm_hi = Float64[])

# impact + cumulative for log-qbar and levels
for (kind, Mt0) in (("impact_log_qbar", Mt_log), ("impact_levels", Mt_lvl))
    for H in (0, 3, 6, 12, 18)
        if kind == "impact_log_qbar" || kind == "impact_levels"
            # for impact use H as horizon; for cumulative use sum — separate loops cleaner
        end
    end
end

# rewrite M2 cleanly
m2 = DataFrame(estimator = String[], measure = String[], H = Int[],
               current = Float64[], current_lo = Float64[], current_hi = Float64[],
               harm = Float64[], harm_lo = Float64[], harm_hi = Float64[])

for H in (0, 3, 6, 12, 18)
    # current impact from leakage.csv
    crow = cur_lk[cur_lk.horizon .== H, :]
    cpt = isempty(crow) ? NaN : crow.leakage[1]
    clo = isempty(crow) ? NaN : crow.lo90[1]
    chi = isempty(crow) ? NaN : crow.hi90[1]

    # harm impact log-qbar
    pt = leak_from_B([Mt_log[H + 1, j] for j in 1:length(series)], aidx, H)
    boots = fill(NaN, NBOOT_M)
    for b in 1:NBOOT_M
        path = [BS_log[b, H + 1, j] * qbar[j] / 1e6 for j in 1:length(series)]
        boots[b] = leak_from_B(path, aidx, H)
    end
    lo, hi = qband(boots)
    push!(m2, ("log_qbar", "impact", H, cpt, clo, chi, pt, lo, hi))

    # harm impact levels
    ptL = leak_from_B([Mt_lvl[H + 1, j] for j in 1:length(series)], aidx, H)
    for b in 1:NBOOT_M
        path = [BS_lvl[b, H + 1, j] for j in 1:length(series)]
        boots[b] = leak_from_B(path, aidx, H)
    end
    lo, hi = qband(boots)
    push!(m2, ("levels", "impact", H, NaN, NaN, NaN, ptL, lo, hi))

    # cumulative log-qbar
    ptC = cum_leak(Mt_log, aidx, H)
    for b in 1:NBOOT_M
        Bb = [BS_log[b, h + 1, j] * qbar[j] / 1e6 for h in 0:Hmax, j in 1:length(series)]
        boots[b] = cum_leak(Bb, aidx, H)
    end
    lo, hi = qband(boots)
    push!(m2, ("log_qbar", "cumulative", H, NaN, NaN, NaN, ptC, lo, hi))

    # cumulative levels
    ptCL = cum_leak(Mt_lvl, aidx, H)
    for b in 1:NBOOT_M
        Bb = BS_lvl[b, :, :]
        boots[b] = cum_leak(Bb, aidx, H)
    end
    lo, hi = qband(boots)
    push!(m2, ("levels", "cumulative", H, NaN, NaN, NaN, ptCL, lo, hi))
end

# world-aggregate F1-style under harmonised
Qworld = fill(NaN, T)
for i in 1:T
    isfinite(Q_mt[i, aidx]) || continue
    vals = [Q_mt[i, j] for j in 1:length(series) if isfinite(Q_mt[i, j])]
    Qworld[i] = sum(vals)
end
Bw = lp_harm_path(Qworld, z_pref, P.rea, P.dates, Hmax)
BSw = fill(NaN, NBOOT_M, Hmax + 1)
rng2 = MersenneTwister(BOOT_SEED)
for b in 1:NBOOT_M
    ii = block_idx(T, BLOCK, rng2)
    BSw[b, :] = lp_harm_path(Qworld[ii], z_pref[ii], P.rea[ii], P.dates[ii], Hmax)
end
for H in (0, 3, 6, 12, 18)
    # "leakage" analogue: 1 + world_response / (-aus_response) if aus falls
    # Under full replacement world≈0; report world response and implied replacement
    aus = Mt_log[H + 1, aidx]
    w = Bw[H + 1]
    # implied leakage from world: (world - aus) / (-aus) = 1 - world/(-aus) if measuring replacement of Aus shortfall
    # If ΔW = ΔAUS + ΔCOMP, leakage = ΔCOMP/(-ΔAUS) = (ΔW - ΔAUS)/(-ΔAUS) = 1 - ΔW/(-ΔAUS) when ΔAUS<0
    impl = (isfinite(aus) && aus < -0.01 && isfinite(w)) ? 1.0 - w / (-aus) : NaN
    boots = fill(NaN, NBOOT_M)
    for b in 1:NBOOT_M
        a = BS_log[b, H + 1, aidx] * qbar[aidx] / 1e6
        ww = BSw[b, H + 1]
        boots[b] = (isfinite(a) && a < -0.01 && isfinite(ww)) ? 1.0 - ww / (-a) : NaN
    end
    lo, hi = qband(boots)
    push!(m2, ("world_aggregate", "implied_leak_impact", H, NaN, NaN, NaN, impl, lo, hi))
    push!(m2, ("world_aggregate", "world_mt_response", H, NaN, NaN, NaN, w,
               qband(BSw[:, H + 1])[1], qband(BSw[:, H + 1])[2]))
end
CSV.write(joinpath(DIAG_OUT, "M2_leakage_current_vs_harm.csv"), m2)

# M3 Japan
println("M: Japan series...")
jp = CSV.read(joinpath(RAW, "bilateral_jpn.csv"), DataFrame)
jp.ds = first.(string.(jp.date), 7)
function jp_series(partner)
    sub = jp[jp.partner .== partner, :]
    m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(sub))
    [haskey(m, string(s)) ? m[string(s)] : NaN for s in P.ds]
end
j_world = jp_series("World")
j_aus = jp_series("Australia")
j_non = [isfinite(j_world[i]) && isfinite(j_aus[i]) ? j_world[i] - j_aus[i] : NaN for i in 1:T]
j_world_l = [isfinite(x) && x > 0 ? log(x) : NaN for x in j_world]
j_aus_l = [isfinite(x) && x > 0 ? log(x) : NaN for x in j_aus]
j_non_l = [isfinite(x) && x > 0 ? log(x) : NaN for x in j_non]
j_world_mt = j_world ./ 1e6
j_aus_mt = j_aus ./ 1e6
j_non_mt = j_non ./ 1e6

m3 = DataFrame(series = String[], h = Int[],
               current_beta = Float64[], current_lo = Float64[], current_hi = Float64[],
               harm_log = Float64[], harm_log_lo = Float64[], harm_log_hi = Float64[],
               harm_mt = Float64[], harm_mt_lo = Float64[], harm_mt_hi = Float64[])

function cur_jp(src, h)
    r = cur_imp[(cur_imp.importer .== "Japan") .& (cur_imp.source .== src) .& (cur_imp.horizon .== h), :]
    isempty(r) ? (NaN, NaN, NaN) : (r.beta[1], r.lo90[1], r.hi90[1])
end

jp_defs = [
    ("Japan_total", j_world_l, j_world_mt, "World"),
    ("Japan_from_Australia", j_aus_l, j_aus_mt, "Australia"),
    ("Japan_non_Australia", j_non_l, j_non_mt, nothing),
]
for (name, yl, ym, src) in jp_defs
    pl = lp_harm_path(yl, z_pref, P.rea, P.dates, 12)
    pm = lp_harm_path(ym, z_pref, P.rea, P.dates, 12)
    BSl = fill(NaN, NBOOT_M, 13)
    BSm = fill(NaN, NBOOT_M, 13)
    rngj = MersenneTwister(BOOT_SEED)
    for b in 1:NBOOT_M
        ii = block_idx(T, BLOCK, rngj)
        BSl[b, :] = lp_harm_path(yl[ii], z_pref[ii], P.rea[ii], P.dates[ii], 12)
        BSm[b, :] = lp_harm_path(ym[ii], z_pref[ii], P.rea[ii], P.dates[ii], 12)
    end
    for h in 0:12
        cb, clo, chi = src === nothing ? (NaN, NaN, NaN) : cur_jp(src, h)
        llo, lhi = qband(BSl[:, h + 1])
        mlo, mhi = qband(BSm[:, h + 1])
        push!(m3, (name, h, cb, clo, chi, pl[h + 1], llo, lhi, pm[h + 1], mlo, mhi))
    end
end
CSV.write(joinpath(DIAG_OUT, "M3_japan_current_vs_harm.csv"), m3)

println("wrote M1 M2 M3 M4")
println("script: diagnostics/task_M.jl")
