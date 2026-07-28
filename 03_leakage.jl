#=
03_leakage.jl
-------------
Estimate the LEAKAGE RATE: when an exogenous disruption cuts Australian coal
exports, what share of the lost volume do competing exporters replace?

Empirical counterpart to the rebound effect in Richter, Mendelevitch & Jotzo
(2018), who obtain 73.3% by simulation for a unilateral Australian export tax.
That number is a model output and has, to our knowledge, never been estimated.

METHOD
------
For each exporter j, local projections of log exports on the narrative shock:

    log Q_{j,t+h} - log Q_{j,t-1} = a + b_{j,h} z_t + lags + month dummies + e

Leakage at horizon h, converting log responses to tonnes with mean volumes:

                    sum_{j != AUS} b_{j,h} * Qbar_j
    leakage_h  =  ----------------------------------
                      | b_{AUS,h} * Qbar_AUS |

    0.00  no replacement (world supply falls one for one)
    0.73  the Richter et al. prediction
    1.00  complete replacement (world supply unchanged; price should not move)

Series are placed on a full monthly CALENDAR with NaN for missing months, so a
reporting gap costs observations without breaking the lag structure. Exporters
are selected on coverage rather than requiring a balanced panel across all of
them, which would collapse the sample to the worst-reported country.

Month dummies matter: cyclones cluster in the Australian wet season.

Run:  julia --project=. 03_leakage.jl
=#

using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random, Printf

Random.seed!(20260728)

const H            = 18
const NLAGS        = 6
const NBOOT        = 1000
const BLOCK        = 12
const SAMPLE_START = "2012-01"    # 2010-11 Comtrade data are unreliable
const MIN_MONTHS   = 110          # coverage needed to enter the main panel
const MIN_END      = "2024-01"    # must still be reporting recently
const MIN_DENSITY  = 0.85         # observed months / span (South Africa is 70%)
const RAW          = joinpath(@__DIR__, "data", "raw")
const OUT          = joinpath(@__DIR__, "out")

const EXPORTERS = [
    ("aus_exports.csv", "Australia",    true),
    ("exp_idn.csv",     "Indonesia",    false),
    ("exp_rus.csv",     "Russia",       false),
    ("exp_zaf.csv",     "South Africa", false),
    ("exp_col.csv",     "Colombia",     false),
    ("exp_usa.csv",     "USA",          false),
    ("exp_can.csv",     "Canada",       false),
    ("exp_mng.csv",     "Mongolia",     false),
]

olsb(y::Vector{Float64}, X::Matrix{Float64}) = hcat(ones(length(y)), X) \ y

function readexp(fname::String)
    path = joinpath(RAW, fname)
    isfile(path) || return nothing
    d = CSV.read(path, DataFrame)
    ncol(d) < 2 && return nothing
    ds = first.(string.(d[!, 1]), 7)
    tn = Float64.(coalesce.(d[!, 2], NaN))
    keep = .!isnan.(tn) .& (tn .> 0) .& (ds .>= SAMPLE_START)
    isempty(ds[keep]) && return nothing
    out = DataFrame(ds = ds[keep], tonnes = tn[keep])
    sort!(out, :ds)
    return unique(out, :ds)
end

"Full monthly calendar of \"YYYY-MM\" strings between two such strings."
function calendar(a::String, b::String)
    d0 = Date(a * "-01"); d1 = Date(b * "-01")
    return [Dates.format(d, "yyyy-mm") for d in d0:Month(1):d1]
end

function monthdummies(ds::Vector{String})
    D = zeros(length(ds), 11)
    for (i, s) in enumerate(ds)
        m = parse(Int, s[6:7])
        m > 1 && (D[i, m - 1] = 1.0)
    end
    return D
end

"""
    lp_betas(y, z, D, H, nlags)

LP coefficients on `z` for horizons 0..H. `y` may contain NaN: a row enters the
regression only if the outcome, the base period, every lag and the shock are
all observed.
"""
function lp_betas(y::Vector{Float64}, z::Vector{Float64}, D::Matrix{Float64},
                  H::Int, nlags::Int)
    T = length(y)
    β = fill(NaN, H + 1)
    for h in 0:H
        rows = Int[]
        for t in (nlags + 1):(T - h)
            ok = isfinite(y[t + h]) && isfinite(y[t - 1]) && isfinite(z[t])
            ok && (ok = all(isfinite, @view y[(t - nlags):(t - 1)]))
            ok && push!(rows, t)
        end
        length(rows) < 40 && continue
        sum(z[rows]) < 2 && continue          # need events in the estimation rows
        L = zeros(length(rows), nlags)
        for (i, t) in enumerate(rows), l in 1:nlags
            L[i, l] = y[t - l]
        end
        X  = hcat(z[rows], L, D[rows, :])
        yy = [y[t + h] - y[t - 1] for t in rows]
        β[h + 1] = olsb(yy, X)[2]
    end
    return β
end

function leakage_path(Q::Matrix{Float64}, isaus::Vector{Bool}, z::Vector{Float64},
                      D::Matrix{Float64}, H::Int, nlags::Int)
    J = size(Q, 2)
    B = fill(NaN, H + 1, J)
    for j in 1:J
        B[:, j] = lp_betas(log.(Q[:, j]), z, D, H, nlags)
    end
    qbar = [mean(filter(isfinite, Q[:, j])) for j in 1:J]

    leak = fill(NaN, H + 1); NUM = fill(NaN, H + 1); DEN = fill(NaN, H + 1)
    SHR = zeros(H + 1)
    a = findfirst(isaus)
    cap = sum(qbar[j] for j in 1:J if j != a)
    for h in 1:(H + 1)
        den = -B[h, a] * qbar[a]          # tonnes Australia loses (>0 if it falls)
        num = 0.0; contrib = 0.0
        for j in 1:J
            j == a && continue
            isnan(B[h, j]) && continue     # skip, do not abort the horizon
            num += B[h, j] * qbar[j]
            contrib += qbar[j]
        end
        NUM[h] = num; DEN[h] = den; SHR[h] = contrib / cap
        # ratio only meaningful when Australia measurably contracts
        (isfinite(den) && den > 0.2e6) && (leak[h] = num / den)
    end
    return leak, B, qbar, NUM, DEN, SHR
end

function block_idx(T::Int, block::Int, rng)
    idx = Int[]
    while length(idx) < T
        s = rand(rng, 1:max(1, T - block + 1))
        append!(idx, s:min(s + block - 1, T))
    end
    return idx[1:T]
end

# ------------------------------------------------------------------- data
println("="^72)
println("Leakage: do competing exporters replace disrupted Australian volume?")
println("="^72)

all_series = Tuple{String,Bool,DataFrame}[]
for (f, lab, isa_) in EXPORTERS
    d = readexp(f)
    d === nothing ? (@warn "missing, skipping" file = f) : push!(all_series, (lab, isa_, d))
end
isempty(all_series) && error("no exporter files -- run 00_pull_comtrade.jl first")

println("\nCOVERAGE")
@printf("%-14s %8s %10s %10s %9s %8s  %s\n", "", "months", "from", "to", "mean Mt", "dense", "status")
sel = Tuple{String,Bool,DataFrame}[]; excl = Tuple{String,Bool,DataFrame}[]
for (lab, isa_, d) in all_series
    dens = nrow(d) / length(calendar(d.ds[1], d.ds[end]))
    keep = isa_ || (nrow(d) >= MIN_MONTHS && d.ds[end] >= MIN_END && dens >= MIN_DENSITY)
    why  = keep ? "included" :
           d.ds[end] < MIN_END ? "EXCLUDED (stale)" :
           dens < MIN_DENSITY  ? "EXCLUDED (gappy)" : "EXCLUDED (short)"
    @printf("%-14s %8d %10s %10s %9.2f %7.0f%%  %s\n", lab, nrow(d), d.ds[1], d.ds[end],
            mean(d.tonnes) / 1e6, 100 * dens, why)
    push!(keep ? sel : excl, (lab, isa_, d))
end
length(sel) < 2 && error("need Australia plus at least one competitor with adequate coverage")

start = maximum([d.ds[1]   for (_, _, d) in sel])
stop  = minimum([d.ds[end] for (_, _, d) in sel])
dates = calendar(max(start, SAMPLE_START), stop)
length(dates) < 60 && error("common span is only $(length(dates)) months")

J = length(sel)
Q = fill(NaN, length(dates), J)
labels = String[]; isaus = Bool[]
for (j, (lab, isa_, d)) in enumerate(sel)
    m = Dict(d.ds .=> d.tonnes)
    for (i, s) in enumerate(dates)
        haskey(m, s) && (Q[i, j] = m[s])
    end
    push!(labels, lab); push!(isaus, isa_)
end

ev = CSV.read(joinpath(@__DIR__, "data", "events.csv"), DataFrame)
evset = Set(first.(string.(ev.date), 7))
z = Float64[s in evset ? 1.0 : 0.0 for s in dates]
D = monthdummies(dates)

@printf("\nspan           : %s to %s (%d months)\n", dates[1], dates[end], length(dates))
@printf("exporters      : %s\n", join(labels, ", "))
@printf("events in span : %d\n", Int(sum(z)))
obs = [count(isfinite, Q[:, j]) for j in 1:J]
@printf("observed months: %s\n\n", join(["$(labels[j])=$(obs[j])" for j in 1:J], ", "))
sum(z) < 3 && error("fewer than 3 events in span -- nothing to identify")

# ------------------------------------------------------------- estimation
leak, B, qbar, NUM, DEN, SHR = leakage_path(Q, isaus, z, D, H, NLAGS)

println("-"^72)
println("EXPORTER RESPONSES to a disruption (log points, + = expands)")
@printf("%-14s %9s %9s %9s %9s %10s\n", "", "h=0", "h=3", "h=6", "h=12", "mean Mt")
for j in 1:J
    @printf("%-14s %9.4f %9.4f %9.4f %9.4f %10.2f\n", labels[j],
            B[1, j], B[4, j], B[7, j], B[13, j], qbar[j] / 1e6)
end

print("\nbootstrapping ($NBOOT draws, block $BLOCK) ... ")
rng = MersenneTwister(20260728); T = length(dates)
BS = fill(NaN, NBOOT, H + 1)
for b in 1:NBOOT
    ii = block_idx(T, BLOCK, rng)
    lk, _, _, _, _, _ = leakage_path(Q[ii, :], isaus, z[ii], D[ii, :], H, NLAGS)
    BS[b, :] = lk
end
println("done")

println("\n" * "-"^72)
println("LEAKAGE RATE (share of lost Australian volume replaced by competitors)")
println("  0.00 = none    0.73 = Richter et al. (2018)    1.00 = complete")
@printf("\n%7s %10s %10s %9s %22s %5s\n",
        "horizon", "AUS (Mt)", "rest (Mt)", "leakage", "90% band", "cov")
for h in [0, 1, 3, 6, 12, 18]
    h > H && continue
    col = filter(isfinite, BS[:, h + 1])
    note = !isfinite(leak[h + 1]) ? "  <- AUS response too small/positive" : ""
    lo = isempty(col) ? NaN : quantile(col, 0.05)
    hi = isempty(col) ? NaN : quantile(col, 0.95)
    @printf("%7d %10.3f %10.3f %9s %22s %4.0f%%%s\n", h,
            -DEN[h + 1] / 1e6, NUM[h + 1] / 1e6,
            isfinite(leak[h+1]) ? @sprintf("%.3f", leak[h+1]) : "n/a",
            isfinite(lo) ? @sprintf("[%+7.3f, %+7.3f]", lo, hi) : "-",
            100 * SHR[h + 1], note)
end
println("\nAUS (Mt) is the tonnage Australia loses; rest (Mt) is what competitors add.")
println("These are the primary objects: the ratio is unstable when AUS is near zero.")

# excluded exporters, estimated individually on their own span
if !isempty(excl)
    println("\n" * "-"^72)
    println("EXCLUDED EXPORTERS (own sample; not in the leakage total above)")
    @printf("%-14s %9s %9s %9s %10s  %s\n", "", "h=0", "h=3", "h=6", "mean Mt", "span")
    for (lab, _, d) in excl
        cal = calendar(d.ds[1], d.ds[end])
        m = Dict(d.ds .=> d.tonnes)
        y = [haskey(m, s) ? m[s] : NaN for s in cal]
        zz = Float64[s in evset ? 1.0 : 0.0 for s in cal]
        bb = lp_betas(log.(y), zz, monthdummies(cal), H, NLAGS)
        @printf("%-14s %9.4f %9.4f %9.4f %10.2f  %s..%s\n", lab,
                bb[1], bb[4], bb[7], mean(d.tonnes) / 1e6, d.ds[1], d.ds[end])
    end
    tot_excl = sum(mean(d.tonnes) for (_, _, d) in excl) / 1e6
    @printf("\nExcluded capacity: %.1f Mt/month vs %.1f Mt/month included.\n",
            tot_excl, sum(qbar[.!isaus]) / 1e6)
    println("If excluded exporters also expand, true leakage is HIGHER than reported.")
end

mkpath(OUT)
lo90 = [(c = filter(isfinite, BS[:, h+1]); isempty(c) ? NaN : quantile(c, 0.05)) for h in 0:H]
hi90 = [(c = filter(isfinite, BS[:, h+1]); isempty(c) ? NaN : quantile(c, 0.95)) for h in 0:H]
CSV.write(joinpath(OUT, "leakage.csv"),
          DataFrame(horizon = 0:H, leakage = leak, lo90 = lo90, hi90 = hi90))
bt = DataFrame(horizon = 0:H)
for j in 1:J; bt[!, Symbol(replace(labels[j], " " => "_"))] = B[:, j]; end
CSV.write(joinpath(OUT, "exporter_responses.csv"), bt)
println("\nwrote out/leakage.csv and out/exporter_responses.csv")

println("""
Reading the results
-------------------
Richter et al. argue competitors CANNOT expand immediately (most run at
capacity), so leakage should be LOW at h=0 and RISE toward 73%. That shape
would validate the mechanism dynamically, not just on average.

The estimate is a LOWER BOUND on policy leakage, for three reasons: a cyclone
is transitory while a tax is permanent, so competitors respond less here;
Russia, South Africa and Mongolia drop out on reporting coverage; and buyers
may run down inventories rather than switch supplier at short horizons. Total
exports also cannot separate new production from diversion away from a
producer's own domestic market.
""")
