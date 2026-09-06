#=
08_estimate_by_group.jl
-----------------------
LP-IV demand elasticities separately for HIGH vs LOW Australian-import-share
importers (advisor add-on).

Quantity = log(World imports) from data/raw/bilateral_<stub>.csv
Price / shock / REA from data/panel.csv (same z_w, real price, DROP_EVENTS)

Groups (clean-margin only by default):
  HIGH : Japan, Taiwan (≥ ~50% Aus)
  LOW  : Germany, Netherlands, Turkey, UK, Philippines, Morocco (< 25%)

Also reports country-by-country and a pooled HIGH / LOW estimate
(country FE within the pool).

Run:  julia --project=. 08_estimate_by_group.jl
=#

using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random, Printf

Random.seed!(20260905)

const H            = 12
const NLAGS        = 6
const NBOOT        = 1000
const SAMPLE_START = "2012-01"
const F_PRICE_MIN  = 10.0
const DROP_EVENTS  = Set(["2022-02", "2022-03"])
const RAW          = joinpath(@__DIR__, "data", "raw")
const OUTDIR       = joinpath(@__DIR__, "out")

# clean-margin contrast set
const HIGH = ["Japan", "Taiwan"]
const LOW  = ["Germany", "Netherlands", "Turkey", "United Kingdom",
              "Philippines", "Morocco"]

const STUB = Dict(
    "Japan" => "jpn", "Taiwan" => "twn", "Korea" => "kor",
    "Germany" => "deu", "Netherlands" => "nld", "Turkey" => "tur",
    "United Kingdom" => "gbr", "Philippines" => "phl", "Morocco" => "mar",
    "Vietnam" => "vnm", "Brazil" => "bra", "India" => "ind",
)

# ---- helpers (calendar-safe LP-IV; same design as 02_estimate.jl) -----
olsb(y, X) = hcat(ones(length(y)), X) \ y
ym(x) = first(string(x), 7)

function monthdummies(dates)::Matrix{Float64}
    T = length(dates); D = zeros(T, 11)
    for (i, d) in enumerate(dates)
        m = month(d isa Date ? d : Date(ym(d) * "-01"))
        m > 1 && (D[i, m - 1] = 1.0)
    end
    return D
end

function lpiv(q, p, z, dates, H, nlags; rea = nothing, month_fe = true,
              cfe::Union{Nothing,Matrix{Float64}} = nothing,
              q_lags::Int = 1)
    # q_lags: how many lags of quantity to include. Sparse Comtrade reporters
    # (EU) often miss months; use q_lags=1 so only consecutive (t, t-1) needed.
    # Price / REA stay at `nlags` (those series are continuous on the panel).
    T = length(q)
    D = month_fe ? monthdummies(dates) : zeros(T, 0)
    eps = fill(NaN, H + 1); bq = fill(NaN, H + 1); bp = fill(NaN, H + 1)
    nobs = fill(0, H + 1)
    qL = max(q_lags, 1)
    need = max(nlags, qL)
    for h in 0:H
        rows = Int[]
        for t in (need + 1):(T - h)
            ok = isfinite(q[t + h]) && isfinite(q[t - 1]) &&
                 isfinite(p[t + h]) && isfinite(p[t - 1]) && isfinite(z[t])
            ok = ok && all(isfinite, @view q[(t - qL):(t - 1)])
            ok = ok && all(isfinite, @view p[(t - nlags):(t - 1)])
            if rea !== nothing
                ok = ok && isfinite(rea[t]) &&
                     all(isfinite, @view rea[(t - nlags):(t - 1)])
            end
            ok && push!(rows, t)
        end
        (length(rows) < 30 || sum(abs, z[rows]) < 1.0) && continue
        n = length(rows)
        kextra = (cfe === nothing ? 0 : size(cfe, 2))
        k = qL + nlags + (rea === nothing ? 0 : nlags + 1) + size(D, 2) + kextra
        X = zeros(n, 1 + k)
        X[:, 1] = z[rows]; c = 2
        for l in 1:qL
            X[:, c] = q[rows .- l]; c += 1
        end
        for l in 1:nlags
            X[:, c] = p[rows .- l]; c += 1
        end
        if rea !== nothing
            X[:, c] = rea[rows]; c += 1
            for l in 1:nlags
                X[:, c] = rea[rows .- l]; c += 1
            end
        end
        size(D, 2) > 0 && (X[:, c:c + size(D, 2) - 1] = D[rows, :]; c += size(D, 2))
        cfe !== nothing && (X[:, c:end] = cfe[rows, :])
        yq = [q[t + h] - q[t - 1] for t in rows]
        yp = [p[t + h] - p[t - 1] for t in rows]
        bq[h + 1] = olsb(yq, X)[2]
        bp[h + 1] = olsb(yp, X)[2]
        eps[h + 1] = bq[h + 1] / bp[h + 1]
        nobs[h + 1] = n
    end
    return eps, bq, bp, nobs
end

function first_stage_F(y, z, dates; q, p, rea = nothing, nlags = NLAGS, q_lags = 1)
    T = length(y); D = monthdummies(dates)
    qL = max(q_lags, 1); need = max(nlags, qL)
    rows = Int[]
    for t in (need + 1):T
        ok = isfinite(y[t]) && isfinite(y[t - 1]) && isfinite(z[t])
        ok = ok && all(isfinite, @view q[(t - qL):(t - 1)])
        ok = ok && all(isfinite, @view p[(t - nlags):(t - 1)])
        if rea !== nothing
            ok = ok && isfinite(rea[t]) &&
                 all(isfinite, @view rea[(t - nlags):(t - 1)])
        end
        ok && push!(rows, t)
    end
    length(rows) < 30 && return (NaN, NaN, NaN)
    n = length(rows)
    k = qL + nlags + (rea === nothing ? 0 : nlags + 1) + size(D, 2)
    X = zeros(n, 1 + k)
    X[:, 1] = z[rows]; c = 2
    for l in 1:qL
        X[:, c] = q[rows .- l]; c += 1
    end
    for l in 1:nlags
        X[:, c] = p[rows .- l]; c += 1
    end
    if rea !== nothing
        X[:, c] = rea[rows]; c += 1
        for l in 1:nlags
            X[:, c] = rea[rows .- l]; c += 1
        end
    end
    X[:, c:end] = D[rows, :]
    yy = [y[t] - y[t - 1] for t in rows]
    b = olsb(yy, X)
    e = yy - hcat(ones(n), X) * b
    s2 = sum(e.^2) / (n - size(X, 2) - 1)
    se = sqrt(s2 * inv(hcat(ones(n), X)' * hcat(ones(n), X))[2, 2])
    return b[2], b[2] / se, (b[2] / se)^2
end

function boot_lpiv(q, p, z, dates, H, nlags, nboot; rea = nothing, cfe = nothing,
                   q_lags = 1)
    T = length(q)
    out = fill(NaN, nboot, H + 1)
    for b in 1:nboot
        w = rand([-1.0, 1.0], T)
        e, _, _, _ = lpiv(q, p, z .* w, dates, H, nlags; rea = rea, cfe = cfe,
                          q_lags = q_lags)
        out[b, :] = e
    end
    return out
end

"World-import series on the panel calendar; NaN where missing."
function world_imports(stub::String, cal_ym::Vector{String})
    path = joinpath(RAW, "bilateral_$(stub).csv")
    !isfile(path) && return fill(NaN, length(cal_ym)), 0.0
    df = CSV.read(path, DataFrame)
    df.ds = ym.(df.date)
    w = df[string.(df.partner) .== "World", :]
    m = Dict(string.(w.ds) .=> Float64.(w.tonnes))
    a = df[string.(df.partner) .== "Australia", :]
    aus_share = isempty(a) || isempty(w) ? NaN :
                sum(Float64.(a.tonnes)) / sum(Float64.(w.tonnes))
    y = [haskey(m, s) ? m[s] : NaN for s in cal_ym]
    return y, aus_share
end

# ---- load panel -------------------------------------------------------
mkpath(OUTDIR)
panel = CSV.read(joinpath(@__DIR__, "data", "panel.csv"), DataFrame)
sort!(panel, :date)
panel = panel[ym.(panel.date) .>= SAMPLE_START, :]
panel.date = [d isa Date ? d : Date(ym(d) * "-01") for d in panel.date]
cal_ym = ym.(panel.date)

lp = [ismissing(x) || x <= 0 ? NaN : log(Float64(x)) for x in panel.p_real]
rea = [ismissing(x) ? NaN : Float64(x) for x in panel.rea]
z = Float64.(coalesce.(panel.z_w, 0.0))
for i in 1:length(z)
    cal_ym[i] in DROP_EVENTS && (z[i] = 0.0)
end
dates = panel.date
lq_aus = [ismissing(x) || x <= 0 ? NaN : log(Float64(x)) for x in panel.q_aus]

println("="^70)
println("LP-IV by Aus-import-share group (HIGH vs LOW)")
println("="^70)
@printf("sample %s..%s  shock=z_w  drop=%s\n\n",
        cal_ym[1], cal_ym[end], join(sort(collect(DROP_EVENTS)), ","))

# common price first stage (same for all — Wald denominator)
# use Aus exports as the q in the control set for the price equation
bp0, tp, Fp = first_stage_F(lp, z, dates; q = lq_aus, p = lp, rea = rea, q_lags = 6)
bq0, tq, Fq = first_stage_F(lq_aus, z, dates; q = lq_aus, p = lp, rea = rea, q_lags = 6)
price_ok = isfinite(Fp) && Fp >= F_PRICE_MIN && bp0 > 0
@printf("COMMON PRICE FIRST STAGE: β=%+.4f  t=%+.2f  F=%.2f  gate=%s\n",
        bp0, tp, Fp, price_ok)
@printf("COMMON q_aus FIRST STAGE: β=%+.4f  t=%+.2f  F=%.2f\n\n", bq0, tq, Fq)

results = DataFrame(group = String[], country = String[], horizon = Int[],
                    elasticity = Float64[], lo90 = Float64[], hi90 = Float64[],
                    bq = Float64[], bp = Float64[], nobs = Int[],
                    aus_share = Float64[], reliable = Bool[], q_lags = Int[])

function run_one(label, group, q_tons, aus_share; q_lags = 1)
    n_ok = count(isfinite, q_tons)
    n_ok < 40 && (@warn "skip (few obs)" country = label n = n_ok; return)
    lx = [isfinite(x) && x > 0 ? log(x) : NaN for x in q_tons]
    e, bqq, bpp, nn = lpiv(lx, lp, z, dates, H, NLAGS; rea = rea, q_lags = q_lags)
    bs = boot_lpiv(lx, lp, z, dates, H, NLAGS, NBOOT; rea = rea, q_lags = q_lags)
    bqi, tqi, Fqi = first_stage_F(lx, z, dates; q = lx, p = lp, rea = rea, q_lags = q_lags)
    @printf("%-18s [%s] Aus=%.0f%%  N=%d  qL=%d  dq/dz β=%+.3f F=%.1f  reliable=%s\n",
            label, group, 100 * aus_share, n_ok, q_lags, bqi, Fqi, price_ok)
    for h in [0, 3, 6, 12]
        colb = filter(!isnan, bs[:, h + 1])
        lo = isempty(colb) ? NaN : quantile(colb, 0.05)
        hi = isempty(colb) ? NaN : quantile(colb, 0.95)
        @printf("    h=%2d  eps=%+.3f  [%+.3f, %+.3f]  βq=%+.3f βp=%+.3f n=%d\n",
                h, e[h + 1], lo, hi, bqq[h + 1], bpp[h + 1], nn[h + 1])
        push!(results, (group, label, h, e[h + 1], lo, hi,
                        bqq[h + 1], bpp[h + 1], nn[h + 1], aus_share, price_ok, q_lags))
    end
end

println("-"^70)
println("COUNTRY-BY-COUNTRY")
println("-"^70)
# Japan has dense Comtrade → q_lags=6; sparse EU/PH reporters → q_lags=1
const QLAGS = Dict("Japan" => 6, "Taiwan" => 6, "Korea" => 6)
for (group, names_) in (("HIGH", HIGH), ("LOW", LOW))
    println("\n### $group")
    for name in names_
        stub = get(STUB, name, nothing)
        stub === nothing && continue
        y, share = world_imports(stub, cal_ym)
        if all(!isfinite, y)
            println("  $name  MISSING bilateral_$(stub).csv — pull with 04c")
            continue
        end
        isnan(share) || @printf("  (file Aus share %.1f%%)\n", 100 * share)
        ql = get(QLAGS, name, 1)
        run_one(name, group, y, coalesce(share, NaN); q_lags = ql)
    end
end

# ---- pooled HIGH / LOW with country FE --------------------------------
println("\n" * "-"^70)
println("POOLED (country FE within group)")
println("-"^70)

function pooled(group::String, names_::Vector{String}; q_lags = 1)
    qs = Float64[]; ps = Float64[]; zs = Float64[]; rs = Float64[]
    ds = Date[]; cids = Int[]
    id = 0
    shares = Float64[]
    used = String[]
    for name in names_
        stub = get(STUB, name, nothing)
        stub === nothing && continue
        y, share = world_imports(stub, cal_ym)
        all(!isfinite, y) && continue
        id += 1
        push!(used, name)
        push!(shares, share)
        for t in 1:length(cal_ym)
            push!(qs, isfinite(y[t]) && y[t] > 0 ? log(y[t]) : NaN)
            push!(ps, lp[t]); push!(zs, z[t]); push!(rs, rea[t])
            push!(ds, dates[t]); push!(cids, id)
        end
    end
    id < 1 && (@warn "no countries for pooled $group"; return)
    T = length(qs)
    G = id
    cfe = G > 1 ? zeros(T, G - 1) : nothing
    if cfe !== nothing
        for t in 1:T
            cids[t] > 1 && (cfe[t, cids[t] - 1] = 1.0)
        end
    end
    e, bqq, bpp, nn = lpiv(qs, ps, zs, ds, H, NLAGS; rea = rs, cfe = cfe, q_lags = q_lags)
    bs = boot_lpiv(qs, ps, zs, ds, H, NLAGS, NBOOT; rea = rs, cfe = cfe, q_lags = q_lags)
    @printf("\nPOOLED %s  countries=%s  mean Aus share=%.0f%%  qL=%d\n",
            group, join(used, ", "), 100 * mean(filter(!isnan, shares)), q_lags)
    for h in [0, 3, 6, 12]
        colb = filter(!isnan, bs[:, h + 1])
        lo = isempty(colb) ? NaN : quantile(colb, 0.05)
        hi = isempty(colb) ? NaN : quantile(colb, 0.95)
        @printf("  h=%2d  eps=%+.3f  [%+.3f, %+.3f]  βq=%+.3f βp=%+.3f n=%d\n",
                h, e[h + 1], lo, hi, bqq[h + 1], bpp[h + 1], nn[h + 1])
        push!(results, (group, "POOLED", h, e[h + 1], lo, hi,
                        bqq[h + 1], bpp[h + 1], nn[h + 1],
                        mean(filter(!isnan, shares)), price_ok, q_lags))
    end
end

pooled("HIGH", HIGH; q_lags = 6)
pooled("LOW", LOW; q_lags = 1)

CSV.write(joinpath(OUTDIR, "elasticities_by_group.csv"), results)
open(joinpath(OUTDIR, "group_first_stage.txt"), "w") do io
    @printf(io, "price β=%+.4f t=%+.2f F=%.2f gate=%s\n", bp0, tp, Fp, price_ok)
    @printf(io, "q_aus β=%+.4f t=%+.2f F=%.2f\n", bq0, tq, Fq)
end
println("\nwrote out/elasticities_by_group.csv")
println("wrote out/group_first_stage.txt")
if !price_ok
    println("\nVERDICT: price first stage still fails — group contrast is descriptive only.")
end
