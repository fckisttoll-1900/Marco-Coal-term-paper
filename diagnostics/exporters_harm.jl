# diagnostics/exporters_harm.jl
# Shared exporter panel for Tasks P–T (harmonised coverage gates).
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

"""Load exporters on panel calendar; return named tuple."""
function load_exporters_harm(P)
    series = Tuple{String,Bool,DataFrame}[]
    for (f, lab, isa_) in EXPORTERS
        d = readexp(f)
        d === nothing && continue
        dens = nrow(d) / length(cal_span(d.ds[1], d.ds[end]))
        keep = isa_ || (nrow(d) >= 110 && d.ds[end] >= "2024-01" && dens >= 0.85)
        keep && push!(series, (lab, isa_, d))
    end
    T = length(P.dates)
    J = length(series)
    labels = [s[1] for s in series]
    isaus = [s[2] for s in series]
    aidx = findfirst(isaus)
    Q_t = fill(NaN, T, J)
    for (j, (_, _, d)) in enumerate(series)
        m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(d))
        for i in 1:T
            haskey(m, string(P.ds[i])) && (Q_t[i, j] = m[string(P.ds[i])])
        end
    end
    qbar = [mean(filter(isfinite, Q_t[:, j])) for j in 1:J]
    Q_mt = Q_t ./ 1e6
    # world including Aus: sum of observed; require Aus observed
    W_inc = fill(NaN, T)
    W_exc = fill(NaN, T)
    for i in 1:T
        isfinite(Q_mt[i, aidx]) || continue
        vals_all = [Q_mt[i, j] for j in 1:J if isfinite(Q_mt[i, j])]
        vals_exc = [Q_mt[i, j] for j in 1:J if j != aidx && isfinite(Q_mt[i, j])]
        W_inc[i] = sum(vals_all)
        W_exc[i] = isempty(vals_exc) ? NaN : sum(vals_exc)
    end
    return (; labels, isaus, aidx, Q_mt, qbar, W_inc, W_exc, J)
end

println("exporters_harm.jl loaded")
