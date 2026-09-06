#=
07_import_shares.jl
-------------------
Answer the referee / advisor question: how dependent is each importer on
Australian coal, and can we split the sample into HIGH vs LOW dependence?

Uses the bilateral Comtrade pulls from 04_pull_bilateral.jl
(data/raw/bilateral_*.csv with columns date, partner, tonnes).

For each importer and year:
    aus_share = tonnes from Australia / tonnes from World

Also writes a country-level summary with a HIGH / MID / LOW label based on
average Aus share over the overlapping sample.

Output:
    out/import_shares_annual.csv
    out/import_shares_summary.csv

Run:  julia --project=. 07_import_shares.jl
=#

using CSV, DataFrames, Statistics, Printf, Dates

const RAW = joinpath(@__DIR__, "data", "raw")
const OUT = joinpath(@__DIR__, "out")

# Thresholds for the advisor contrast (tunable; print them in the summary).
const HIGH_CUT = 0.50   # ≥ 50% of imports from Australia
const LOW_CUT  = 0.25   # <  25% → LOW; else MID

function classify(share::Float64)
    share >= HIGH_CUT && return "HIGH"
    share <  LOW_CUT  && return "LOW"
    return "MID"
end

"Load one bilateral file; return nothing if missing/empty."
function load_bilateral(path::String)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    nrow(df) == 0 && return nothing
    df.date    = string.(df.date)
    df.partner = string.(df.partner)
    df.tonnes  = Float64.(df.tonnes)
    return df
end

function annual_shares(df::DataFrame, importer::String)
    df = copy(df)
    df.year = first.(df.date, 4)
    rows = NamedTuple[]
    for g in groupby(df, :year)
        w = g[g.partner .== "World", :tonnes]
        isempty(w) && continue
        world = sum(w)
        world > 0 || continue
        aus = sum(g[g.partner .== "Australia", :tonnes]; init = 0.0)
        # residual = world − sum of named origins (not used for share, just QA)
        push!(rows, (importer = importer,
                     year     = g.year[1],
                     world_mt = world / 1e6,
                     aus_mt   = aus / 1e6,
                     aus_share = aus / world))
    end
    return DataFrame(rows)
end

function main()
    mkpath(OUT)

    files = filter(f -> startswith(f, "bilateral_") && endswith(f, ".csv"),
                   readdir(RAW))
    isempty(files) && error("no bilateral_*.csv in $RAW -- run 04_pull_bilateral.jl")

    annual = DataFrame()
    println("="^66)
    println("Australian share of coal imports (HS 2701), by importer")
    println("="^66)
    @printf("HIGH ≥ %.0f%%   LOW < %.0f%%   else MID\n\n",
            100 * HIGH_CUT, 100 * LOW_CUT)

    summaries = NamedTuple[]

    for fname in sort(files)
        # bilateral_jpn.csv → Japan (title-case from ISO-ish stub)
        stub = replace(fname, r"^bilateral_" => "", r"\.csv$" => "")
        label = Dict("jpn" => "Japan", "kor" => "Korea", "twn" => "Taiwan",
                     "ind" => "India", "chn" => "China", "tur" => "Turkey",
                     "vnm" => "Vietnam", "deu" => "Germany", "nld" => "Netherlands",
                     "bra" => "Brazil", "tha" => "Thailand", "phl" => "Philippines",
                     "mys" => "Malaysia", "pak" => "Pakistan", "esp" => "Spain",
                     "gbr" => "United Kingdom", "mar" => "Morocco")
        importer = get(label, stub, uppercasefirst(stub))

        df = load_bilateral(joinpath(RAW, fname))
        df === nothing && continue
        ann = annual_shares(df, importer)
        isempty(ann) && continue
        append!(annual, ann)

        # overall share = total Aus / total World over available years
        share = sum(ann.aus_mt) / sum(ann.world_mt)
        grp   = classify(share)
        push!(summaries, (importer = importer,
                          file = fname,
                          years = nrow(ann),
                          span = "$(minimum(ann.year))-$(maximum(ann.year))",
                          world_mt_yr = mean(ann.world_mt),
                          aus_share = share,
                          group = grp))

        @printf("%-16s  Aus share %5.1f%%  [%s]  %s  (%.1f Mt/yr imports)\n",
                importer, 100 * share, ann.year[1] * "…" * ann.year[end],
                grp, mean(ann.world_mt))
        for r in eachrow(ann)
            @printf("    %s  %5.1f%%\n", r.year, 100 * r.aus_share)
        end
        println()
    end

    summary = DataFrame(summaries)
    sort!(summary, :aus_share, rev = true)
    sort!(annual, [:importer, :year])

    CSV.write(joinpath(OUT, "import_shares_annual.csv"), annual)
    CSV.write(joinpath(OUT, "import_shares_summary.csv"), summary)

    println("-"^66)
    println("GROUPING FOR THE PAPER")
    high = summary[summary.group .== "HIGH", :importer]
    mid  = summary[summary.group .== "MID",  :importer]
    low  = summary[summary.group .== "LOW",  :importer]
    println("  HIGH (Aus-dependent):     ", isempty(high) ? "(none yet)" : join(high, ", "))
    println("  MID:                      ", isempty(mid)  ? "(none yet)" : join(mid, ", "))
    println("  LOW  (not Aus-dependent): ", isempty(low)  ? "(none yet)" : join(low, ", "))
    println()
    println("wrote out/import_shares_annual.csv")
    println("wrote out/import_shares_summary.csv")

    if isempty(low)
        println("""

        No LOW-dependence importers on disk yet. Add them in 04_pull_bilateral.jl
        (Turkey, Vietnam, Germany, …) and re-pull, then re-run this script.
        """)
    end
end

main()
