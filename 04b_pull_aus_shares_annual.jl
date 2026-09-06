#=
04b_pull_aus_shares_annual.jl
-----------------------------
Fast path for the advisor question: annual Comtrade (C/A/HS) for
World + Australia only, enough to classify HIGH vs LOW dependence.

Writes data/raw/aus_share_annual.csv with columns:
    importer, year, world_t, aus_t, aus_share

Also refreshes out/import_shares_summary.csv style grouping.

Uses curl under the hood (Julia HTTP was TLS-flaky against Comtrade from
this environment). Preview endpoint: one year per call, no key needed.

Run:  julia --project=. 04b_pull_aus_shares_annual.jl
      julia --project=. 04b_pull_aus_shares_annual.jl 2015 2023
=#

using CSV, DataFrames, JSON3, Printf

const RAW = joinpath(@__DIR__, "data", "raw")
const OUT = joinpath(@__DIR__, "out")
const CMD = "2701"

# (stub, M49 code, label, clean consumption margin?)
const IMPORTERS = [
    ("jpn", "392", "Japan",       true),
    ("kor", "410", "Korea",       true),
    ("twn", "490", "Taiwan",      true),
    ("chn", "156", "China",       false),
    ("ind", "699", "India",       false),
    ("tur", "792", "Turkey",      true),
    ("vnm", "704", "Vietnam",     true),
    ("deu", "276", "Germany",     true),
    ("nld", "528", "Netherlands", true),
    ("bra", "76",  "Brazil",      true),
    ("tha", "764", "Thailand",    true),
    ("phl", "608", "Philippines", true),
    ("esp", "724", "Spain",       true),
    ("gbr", "826", "United Kingdom", true),
    ("mar", "504", "Morocco",     true),
]

const HIGH_CUT = 0.50
const LOW_CUT  = 0.25

classify(s) = s >= HIGH_CUT ? "HIGH" : s < LOW_CUT ? "LOW" : "MID"

function curl_json(url::String)
    tmp = tempname()
    try
        run(`curl -sS -m 90 -o $tmp $url`)
        return JSON3.read(read(tmp, String))
    finally
        isfile(tmp) && rm(tmp)
    end
end

function fetch_year(reporter::String, year::Int)
    url = "https://comtradeapi.un.org/public/v1/preview/C/A/HS?" *
          "reporterCode=$(reporter)&partnerCode=0,36&partner2Code=0&" *
          "period=$(year)&flowCode=M&cmdCode=$(CMD)&customsCode=C00&motCode=0"
    j = curl_json(url)
    data = haskey(j, :data) && j.data !== nothing ? collect(j.data) : Any[]
    world = 0.0; aus = 0.0
    for rec in data
        cc = string(get(rec, :customsCode, "C00"))
        (cc == "C00" || cc == "") || continue
        w = get(rec, :netWgt, nothing)
        (w === nothing || w === missing) && continue
        wv = try Float64(w) catch; continue end
        wv > 0 || continue
        pc = string(get(rec, :partnerCode, ""))
        pc == "0"  && (world = max(world, wv))
        pc == "36" && (aus   = max(aus, wv))
    end
    return world, aus
end

function main()
    y0 = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2015
    y1 = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2023
    mkpath(RAW); mkpath(OUT)

    rows = NamedTuple[]
    @printf("Annual Comtrade Aus shares, %d–%d (kg → share)\n\n", y0, y1)

    for (stub, code, label, clean) in IMPORTERS
        @printf("%s%s\n", label, clean ? "" : "  (domestic production)")
        for y in y0:y1
            world, aus = try
                fetch_year(code, y)
            catch err
                @warn "failed" importer = label year = y error = err
                sleep(1.0); continue
            end
            sleep(0.6)   # be polite to preview endpoint
            world == 0 && (@printf("  %d  -- no World total\n", y); continue)
            share = aus / world
            @printf("  %d  Aus share %5.1f%%  (World %.1f Mt)\n",
                    y, 100 * share, world / 1e9)   # kg → Mt: /1e9
            push!(rows, (importer = label, stub = stub, year = y,
                         world_mt = world / 1e9, aus_mt = aus / 1e9,
                         aus_share = share, clean_margin = clean))
        end
        println()
    end

    isempty(rows) && error("nothing retrieved from Comtrade")
    df = DataFrame(rows)
    CSV.write(joinpath(RAW, "aus_share_annual.csv"), df)

    # country summary
    sums = combine(groupby(df, [:importer, :stub, :clean_margin])) do g
        DataFrame(years = nrow(g),
                  span = "$(minimum(g.year))-$(maximum(g.year))",
                  world_mt_yr = mean(g.world_mt),
                  aus_share = sum(g.aus_mt) / sum(g.world_mt))
    end
    sums.group = classify.(sums.aus_share)
    sort!(sums, :aus_share, rev = true)
    CSV.write(joinpath(OUT, "import_shares_summary.csv"), sums)
    CSV.write(joinpath(OUT, "import_shares_annual.csv"), df)

    println("="^66)
    println("GROUPING")
    for g in ("HIGH", "MID", "LOW")
        names_ = sums[sums.group .== g, :importer]
        println("  $g: ", isempty(names_) ? "(none)" : join(names_, ", "))
    end
    println("\nwrote data/raw/aus_share_annual.csv")
    println("wrote out/import_shares_summary.csv")
    println("wrote out/import_shares_annual.csv")
end

main()
