# diagnostics/task_F.jl — Replacement without a ratio
include(joinpath(@__DIR__, "common.jl"))

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
calendar(a, b) = [Dates.format(d, "yyyy-mm") for d in Date(a * "-01"):Month(1):Date(b * "-01")]

series = Tuple{String,Bool,DataFrame}[]
for (f, lab, isa_) in EXPORTERS
    d = readexp(f)
    d === nothing && continue
    dens = nrow(d) / length(calendar(d.ds[1], d.ds[end]))
    keep = isa_ || (nrow(d) >= 110 && d.ds[end] >= "2024-01" && dens >= 0.85)
    keep && push!(series, (lab, isa_, d))
end

start = maximum(d.ds[1] for (_, _, d) in series)
stop  = minimum(d.ds[end] for (_, _, d) in series)
dates_e = calendar(max(start, SAMPLE_START), stop)
J = length(series)
Q = fill(NaN, length(dates_e), J)   # tonnes
labels = String[]; isaus = Bool[]
for (j, (lab, isa_, d)) in enumerate(series)
    m = Dict(d.ds .=> Float64.(d.tonnes))
    for (i, s) in enumerate(dates_e)
        haskey(m, s) && (Q[i, j] = m[s])
    end
    push!(labels, lab); push!(isaus, isa_)
end
aidx = findfirst(isaus)
evset = Set(string.(ym.(CSV.read(EVENTS_PATH, DataFrame).date)))
z = Float64[s in evset ? 1.0 : 0.0 for s in dates_e]
dates_ed = [Date(s * "-01") for s in dates_e]
qbar = [mean(filter(isfinite, Q[:, j])) for j in 1:J]

# World total = sum across included exporters (same coverage set as leakage)
# Missing months: sum of available; if Australia missing, world NaN for that month
Qworld = fill(NaN, length(dates_e))
for i in 1:length(dates_e)
    vals = [Q[i, j] for j in 1:J if isfinite(Q[i, j])]
    isempty(vals) && continue
    # require Australia observed to avoid composition bias
    isfinite(Q[i, aidx]) || continue
    Qworld[i] = sum(vals)
end
Qworld_mt = Qworld ./ 1e6
Q_mt = Q ./ 1e6

# F1: LP of world total Mt
f1 = DataFrame(h = Int[], beta_mt = Float64[], se_ols = Float64[],
               mbb_lo90 = Float64[], mbb_hi90 = Float64[])
rng = MersenneTwister(BOOT_SEED)
for h in 0:18
    β, se, _, _ = lp_simple(Qworld_mt, z, dates_ed, h)
    mbb = fill(NaN, NBOOT_B)
    T = length(z)
    for b in 1:NBOOT_B
        ii = block_idx(T, BLOCK, rng)
        bb, _, _, _ = lp_simple(Qworld_mt[ii], z[ii], dates_ed[ii], h)
        mbb[b] = bb
    end
    c = filter(isfinite, mbb)
    push!(f1, (h, β, se,
               isempty(c) ? NaN : quantile(c, 0.05),
               isempty(c) ? NaN : quantile(c, 0.95)))
end
CSV.write(joinpath(DIAG_OUT, "F1_world_total_mt.csv"), f1)

# F2: level Mt LPs vs log-then-rescale
f2 = DataFrame(exporter = String[], h = Int[],
               beta_log = Float64[], mt_from_log = Float64[],
               beta_level_mt = Float64[])
for j in 1:J
    ylog = log.(Q[:, j])
    ymt = Q_mt[:, j]
    for h in 0:18
        bl, _, _, _ = lp_simple(ylog, z, dates_ed, h)
        bm, _, _, _ = lp_simple(ymt, z, dates_ed, h)
        push!(f2, (labels[j], h, bl, bl * qbar[j] / 1e6, bm))
    end
end
CSV.write(joinpath(DIAG_OUT, "F2_log_vs_level_mt.csv"), f2)

# F3/F4: cumulative responses and cumulative leakage with ratio-in-bootstrap
function path_mt_log(Q, z, dates, qbar)
    J = size(Q, 2)
    B = fill(NaN, 19, J)
    for j in 1:J, h in 0:18
        β, _, _, _ = lp_simple(log.(Q[:, j]), z, dates, h)
        B[h + 1, j] = β * qbar[j] / 1e6
    end
    return B
end

Bmt = path_mt_log(Q, z, dates_ed, qbar)
f3 = DataFrame(H = Int[], CumAUS = Float64[], CumCOMP = Float64[],
               cum_leakage = Float64[],
               cum_leak_mbb_lo90 = Float64[], cum_leak_mbb_hi90 = Float64[])

function cumuls(Bmt, aidx, H)
    cumA = sum(Bmt[h + 1, aidx] for h in 0:H)
    cumC = 0.0
    for j in 1:size(Bmt, 2)
        j == aidx && continue
        for h in 0:H
            isfinite(Bmt[h + 1, j]) && (cumC += Bmt[h + 1, j])
        end
    end
    leak = (isfinite(cumA) && abs(cumA) > 1e-8) ? -cumC / cumA : NaN
    return cumA, cumC, leak
end

for H in (0, 3, 6, 12, 18)
    cA, cC, lk = cumuls(Bmt, aidx, H)
    # bootstrap: recompute paths and ratio inside each draw
    rng = MersenneTwister(BOOT_SEED + H)
    T = length(z)
    boots = fill(NaN, NBOOT_B)
    for b in 1:NBOOT_B
        ii = block_idx(T, BLOCK, rng)
        Bb = path_mt_log(Q[ii, :], z[ii], dates_ed[ii], qbar)
        _, _, boots[b] = cumuls(Bb, aidx, H)
    end
    c = filter(isfinite, boots)
    push!(f3, (H, cA, cC, lk,
               isempty(c) ? NaN : quantile(c, 0.05),
               isempty(c) ? NaN : quantile(c, 0.95)))
end
CSV.write(joinpath(DIAG_OUT, "F3F4_cumulative_leakage.csv"), f3)

# F5: document current Table 5 construction
open(joinpath(DIAG_OUT, "F5_current_leakage_bands.txt"), "w") do io
    println(io, "CURRENT Table 5 / out/leakage.csv bands:")
    println(io, "Constructed in 03_leakage.jl by recomputing the leakage RATIO inside")
    println(io, "each moving-block bootstrap draw (BS[b,:] = lk from leakage_path),")
    println(io, "then taking 5th/95th percentiles of those draws.")
    println(io, "NOT constructed by dividing marginal confidence bands.")
    println(io, "Evidence: 03_leakage.jl lines ~221-230 (bootstrap loop) and ~297-300 (quantiles).")
end

println("wrote F1–F5 outputs")
println("script: diagnostics/task_F.jl")
