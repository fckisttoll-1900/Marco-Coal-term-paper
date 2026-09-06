# diagnostics/task_T.jl — HIGH/LOW importer groups under harmonised LP
include(joinpath(@__DIR__, "harm_lp.jl"))

const NBOOT = 1000   # bands; full 10k would be very slow across many series
const Hmax = 12

P = load_panel()
z_pref = z_preferred(P)
T = length(P.dates)

# HIGH: Japan; LOW: DE, NL, TR, UK, PH
const GROUPS = [
    ("HIGH", "Japan", "bilateral_jpn.csv"),
    ("LOW", "Germany", "bilateral_deu.csv"),
    ("LOW", "Netherlands", "bilateral_nld.csv"),
    ("LOW", "Turkey", "bilateral_tur.csv"),
    ("LOW", "United Kingdom", "bilateral_gbr.csv"),
    ("LOW", "Philippines", "bilateral_phl.csv"),
]

function pull_partner(path, partner)
    !isfile(path) && return nothing
    d = CSV.read(path, DataFrame)
    !("partner" in names(d)) && return nothing
    d.ds = first.(string.(d.date), 7)
    sub = d[d.partner .== partner, :]
    isempty(sub) && return nothing
    m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(sub))
    return [haskey(m, string(s)) ? m[string(s)] : NaN for s in P.ds]
end

qband(v) = (c = filter(isfinite, v); isempty(c) ? (NaN, NaN) : (quantile(c, 0.05), quantile(c, 0.95)))

out = DataFrame(group = String[], country = String[], series = String[],
                h = Int[], beta = Float64[], lo90 = Float64[], hi90 = Float64[],
                n_obs_hint = Int[])

println("T: HIGH/LOW under harmonised LP...")
for (group, country, fname) in GROUPS
    path = joinpath(RAW, fname)
    w = pull_partner(path, "World")
    a = pull_partner(path, "Australia")
    if w === nothing || a === nothing
        @warn "missing bilateral" country fname
        continue
    end
    non = [isfinite(w[i]) && isfinite(a[i]) ? w[i] - a[i] : NaN for i in 1:T]
    defs = [
        ("total_log", [isfinite(x) && x > 0 ? log(x) : NaN for x in w]),
        ("from_aus_log", [isfinite(x) && x > 0 ? log(x) : NaN for x in a]),
        ("non_aus_log", [isfinite(x) && x > 0 ? log(x) : NaN for x in non]),
        ("total_mt", w ./ 1e6),
        ("from_aus_mt", a ./ 1e6),
        ("non_aus_mt", non ./ 1e6),
    ]
    for (sname, y) in defs
        n_ok = count(isfinite, y)
        pathβ = lp_harm_path(y, z_pref, P.rea, P.dates, Hmax)
        BS = fill(NaN, NBOOT, Hmax + 1)
        rng = MersenneTwister(BOOT_SEED)
        for b in 1:NBOOT
            ii = block_idx(T, BLOCK, rng)
            BS[b, :] = lp_harm_path(y[ii], z_pref[ii], P.rea[ii], P.dates[ii], Hmax)
        end
        for h in 0:Hmax
            lo, hi = qband(BS[:, h + 1])
            push!(out, (group, country, sname, h, pathβ[h + 1], lo, hi, n_ok))
        end
    end
    @printf("  done %s/%s\n", group, country)
end

# LOW group pooled: average of country betas (descriptive) at each h for key series
# Also estimate pooled LOW total by summing Mt across LOW reporters when all present
low_files = [(c, f) for (g, c, f) in GROUPS if g == "LOW"]
low_world = []
low_aus = []
for (c, f) in low_files
    w = pull_partner(joinpath(RAW, f), "World")
    a = pull_partner(joinpath(RAW, f), "Australia")
    (w === nothing || a === nothing) && continue
    push!(low_world, w); push!(low_aus, a)
end
if !isempty(low_world)
    # sum across LOW countries (Mt), monthwise when any present
    Wsum = fill(NaN, T); Asum = fill(NaN, T)
    for i in 1:T
        ws = [w[i] for w in low_world if isfinite(w[i])]
        as_ = [a[i] for a in low_aus if isfinite(a[i])]
        isempty(ws) || (Wsum[i] = sum(ws))
        isempty(as_) || (Asum[i] = sum(as_))
    end
    Nsum = [isfinite(Wsum[i]) && isfinite(Asum[i]) ? Wsum[i] - Asum[i] : NaN for i in 1:T]
    for (sname, y) in (("LOW_pool_total_mt", Wsum ./ 1e6),
                       ("LOW_pool_from_aus_mt", Asum ./ 1e6),
                       ("LOW_pool_non_aus_mt", Nsum ./ 1e6))
        pathβ = lp_harm_path(y, z_pref, P.rea, P.dates, Hmax)
        BS = fill(NaN, NBOOT, Hmax + 1)
        rng = MersenneTwister(BOOT_SEED)
        for b in 1:NBOOT
            ii = block_idx(T, BLOCK, rng)
            BS[b, :] = lp_harm_path(y[ii], z_pref[ii], P.rea[ii], P.dates[ii], Hmax)
        end
        for h in 0:Hmax
            lo, hi = qband(BS[:, h + 1])
            push!(out, ("LOW", "POOLED", sname, h, pathβ[h + 1], lo, hi, count(isfinite, y)))
        end
    end
end

CSV.write(joinpath(DIAG_OUT, "T_high_low_harmonised.csv"), out)
println("wrote T_high_low_harmonised.csv")
println("script: diagnostics/task_T.jl")
